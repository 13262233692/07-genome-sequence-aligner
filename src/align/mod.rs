use crate::fasta::Base;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TraceDirection {
    Diagonal,
    Up,
    Left,
}

#[derive(Clone)]
pub struct AlignmentConfig {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_open_penalty: i32,
    pub gap_extend_penalty: i32,
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        AlignmentConfig {
            match_score: 2,
            mismatch_penalty: -3,
            gap_open_penalty: -5,
            gap_extend_penalty: -2,
        }
    }
}

impl AlignmentConfig {
    pub fn simple(match_score: i32, mismatch_penalty: i32, gap_penalty: i32) -> Self {
        AlignmentConfig {
            match_score,
            mismatch_penalty,
            gap_open_penalty: gap_penalty,
            gap_extend_penalty: gap_penalty,
        }
    }
}

pub fn score_pair(a: Base, b: Base, config: &AlignmentConfig) -> i32 {
    if a == Base::N || b == Base::N {
        return config.mismatch_penalty / 2;
    }
    if a == b {
        config.match_score
    } else {
        config.mismatch_penalty
    }
}

pub struct DPMatrix {
    pub scores: Vec<i32>,
    pub traceback: Vec<TraceDirection>,
    pub rows: usize,
    pub cols: usize,
}

impl DPMatrix {
    pub fn new(rows: usize, cols: usize) -> Self {
        let size = rows * cols;
        DPMatrix {
            scores: vec![0i32; size],
            traceback: vec![TraceDirection::Diagonal; size],
            rows,
            cols,
        }
    }

    #[inline]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        i * self.cols + j
    }

    #[inline]
    pub fn get_score(&self, i: usize, j: usize) -> i32 {
        self.scores[self.idx(i, j)]
    }

    #[inline]
    pub fn set_score(&mut self, i: usize, j: usize, val: i32) {
        let idx = self.idx(i, j);
        self.scores[idx] = val;
    }

    #[inline]
    pub fn get_trace(&self, i: usize, j: usize) -> TraceDirection {
        self.traceback[self.idx(i, j)]
    }

    #[inline]
    pub fn set_trace(&mut self, i: usize, j: usize, dir: TraceDirection) {
        let idx = self.idx(i, j);
        self.traceback[idx] = dir;
    }
}

pub struct AlignmentResult {
    pub score: i32,
    pub query_aligned: Vec<u8>,
    pub target_aligned: Vec<u8>,
    pub cigar: String,
    pub mismatches: u32,
    pub insertions: u32,
    pub deletions: u32,
}

pub fn needleman_wunsch(
    query: &[Base],
    target: &[Base],
    config: &AlignmentConfig,
) -> AlignmentResult {
    let m = query.len() + 1;
    let n = target.len() + 1;
    let mut matrix = DPMatrix::new(m, n);

    for i in 1..m {
        matrix.set_score(i, 0, config.gap_open_penalty + (i as i32 - 1) * config.gap_extend_penalty);
        matrix.set_trace(i, 0, TraceDirection::Up);
    }
    for j in 1..n {
        matrix.set_score(0, j, config.gap_open_penalty + (j as i32 - 1) * config.gap_extend_penalty);
        matrix.set_trace(0, j, TraceDirection::Left);
    }

    let use_simd = std::is_x86_feature_detected!("avx2") && n > 32;

    for i in 1..m {
        if use_simd {
            #[cfg(target_arch = "x86_64")]
            unsafe {
                crate::simd::avx2::fill_row_simd(&mut matrix, i, query, target, config);
            }
            #[cfg(not(target_arch = "x86_64"))]
            {
                fill_row_scalar(&mut matrix, i, query, target, config);
            }
        } else {
            fill_row_scalar(&mut matrix, i, query, target, config);
        }
    }

    crate::traceback::traceback(&matrix, query, target)
}

fn fill_row_scalar(
    matrix: &mut DPMatrix,
    i: usize,
    query: &[Base],
    target: &[Base],
    config: &AlignmentConfig,
) {
    let q_base = query[i - 1];
    let cols = matrix.cols;

    for j in 1..cols {
        let diag = matrix.get_score(i - 1, j - 1) + score_pair(q_base, target[j - 1], config);
        let up = matrix.get_score(i - 1, j) + config.gap_extend_penalty;
        let left = matrix.get_score(i, j - 1) + config.gap_extend_penalty;

        let (best_score, best_dir) = if diag >= up && diag >= left {
            (diag, TraceDirection::Diagonal)
        } else if up >= left {
            (up, TraceDirection::Up)
        } else {
            (left, TraceDirection::Left)
        };

        matrix.set_score(i, j, best_score);
        matrix.set_trace(i, j, best_dir);
    }
}
