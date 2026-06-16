use crate::fasta::Base;
use std::fmt::Debug;
use std::ops::{Add, Sub};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum AffineState {
    Match,
    Insert,
    Delete,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct TraceDirection {
    pub state: AffineState,
}

impl TraceDirection {
    pub fn new_match() -> Self {
        TraceDirection { state: AffineState::Match }
    }
    pub fn new_insert() -> Self {
        TraceDirection { state: AffineState::Insert }
    }
    pub fn new_delete() -> Self {
        TraceDirection { state: AffineState::Delete }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ScoreKind {
    I16,
    I32,
    I64,
}

pub trait Score:
    Copy
    + Clone
    + Debug
    + PartialOrd
    + Ord
    + Eq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Into<i64>
    + Send
    + Sync
    + 'static
{
    const MIN: Self;
    const MAX: Self;
    const KIND: ScoreKind;
    fn checked_add(self, rhs: Self) -> Option<Self>;
    fn saturating_add(self, rhs: Self) -> Self;
    fn from_i32(v: i32) -> Self;
}

impl Score for i16 {
    const MIN: Self = i16::MIN;
    const MAX: Self = i16::MAX;
    const KIND: ScoreKind = ScoreKind::I16;

    #[inline]
    fn checked_add(self, rhs: Self) -> Option<Self> {
        i16::checked_add(self, rhs)
    }

    #[inline]
    fn saturating_add(self, rhs: Self) -> Self {
        i16::saturating_add(self, rhs)
    }

    #[inline]
    fn from_i32(v: i32) -> Self {
        v as i16
    }
}

impl Score for i32 {
    const MIN: Self = i32::MIN;
    const MAX: Self = i32::MAX;
    const KIND: ScoreKind = ScoreKind::I32;

    #[inline]
    fn checked_add(self, rhs: Self) -> Option<Self> {
        i32::checked_add(self, rhs)
    }

    #[inline]
    fn saturating_add(self, rhs: Self) -> Self {
        i32::saturating_add(self, rhs)
    }

    #[inline]
    fn from_i32(v: i32) -> Self {
        v
    }
}

impl Score for i64 {
    const MIN: Self = i64::MIN;
    const MAX: Self = i64::MAX;
    const KIND: ScoreKind = ScoreKind::I64;

    #[inline]
    fn checked_add(self, rhs: Self) -> Option<Self> {
        i64::checked_add(self, rhs)
    }

    #[inline]
    fn saturating_add(self, rhs: Self) -> Self {
        i64::saturating_add(self, rhs)
    }

    #[inline]
    fn from_i32(v: i32) -> Self {
        v as i64
    }
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

    pub fn required_score_kind(&self, query_len: usize, target_len: usize) -> ScoreKind {
        let max_len = query_len.max(target_len) as i64;
        let match_s = self.match_score as i64;
        let ext = self.gap_extend_penalty as i64;
        let open = self.gap_open_penalty as i64;
        let min_per_step = if ext < open { ext } else { open };
        let max_possible = match_s * max_len;
        let min_possible = min_per_step * max_len;

        if max_possible <= i16::MAX as i64 && min_possible >= i16::MIN as i64 {
            ScoreKind::I16
        } else if max_possible <= i32::MAX as i64 && min_possible >= i32::MIN as i64 {
            ScoreKind::I32
        } else {
            ScoreKind::I64
        }
    }
}

#[inline]
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

pub struct GenericDPMatrix<S: Score> {
    pub mat_m: Vec<S>,
    pub mat_x: Vec<S>,
    pub mat_y: Vec<S>,
    pub traceback: Vec<TraceDirection>,
    pub rows: usize,
    pub cols: usize,
}

impl<S: Score> GenericDPMatrix<S> {
    pub fn new(rows: usize, cols: usize) -> Self {
        let size = rows * cols;
        GenericDPMatrix {
            mat_m: vec![S::from_i32(0); size],
            mat_x: vec![S::from_i32(0); size],
            mat_y: vec![S::from_i32(0); size],
            traceback: vec![TraceDirection::new_match(); size],
            rows,
            cols,
        }
    }

    #[inline]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        debug_assert!(i < self.rows, "row out of bounds: i={}, rows={}", i, self.rows);
        debug_assert!(j < self.cols, "col out of bounds: j={}, cols={}", j, self.cols);
        i * self.cols + j
    }

    #[inline]
    pub fn get_m(&self, i: usize, j: usize) -> S {
        self.mat_m[self.idx(i, j)]
    }
    #[inline]
    pub fn set_m(&mut self, i: usize, j: usize, val: S) {
        let idx = self.idx(i, j);
        self.mat_m[idx] = val;
    }
    #[inline]
    pub fn get_x(&self, i: usize, j: usize) -> S {
        self.mat_x[self.idx(i, j)]
    }
    #[inline]
    pub fn set_x(&mut self, i: usize, j: usize, val: S) {
        let idx = self.idx(i, j);
        self.mat_x[idx] = val;
    }
    #[inline]
    pub fn get_y(&self, i: usize, j: usize) -> S {
        self.mat_y[self.idx(i, j)]
    }
    #[inline]
    pub fn set_y(&mut self, i: usize, j: usize, val: S) {
        let idx = self.idx(i, j);
        self.mat_y[idx] = val;
    }

    #[inline]
    pub fn get_score(&self, i: usize, j: usize) -> S {
        let k = self.idx(i, j);
        let m = self.mat_m[k];
        let x = self.mat_x[k];
        let y = self.mat_y[k];
        if m >= x && m >= y {
            m
        } else if x >= y {
            x
        } else {
            y
        }
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

pub type DPMatrix = GenericDPMatrix<i32>;

pub struct AlignmentResult {
    pub score: i64,
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
    let m = query.len();
    let n = target.len();
    let total = m.checked_mul(n).unwrap_or(usize::MAX);
    let too_large = total > 1_000_000 * 2000;

    if too_large {
        return self::hirschberg::hirschberg_align(query, target, config);
    }

    match config.required_score_kind(m, n) {
        ScoreKind::I16 => run_generic_nw::<i16>(query, target, config),
        ScoreKind::I32 => run_generic_nw::<i32>(query, target, config),
        ScoreKind::I64 => run_generic_nw::<i64>(query, target, config),
    }
}

fn run_generic_nw<S: Score>(
    query: &[Base],
    target: &[Base],
    config: &AlignmentConfig,
) -> AlignmentResult {
    let m = query.len() + 1;
    let n = target.len() + 1;

    let mut matrix = GenericDPMatrix::<S>::new(m, n);
    let gap_open = S::from_i32(config.gap_open_penalty);
    let inf = S::MIN.saturating_add(S::from_i32(1000));

    matrix.set_m(0, 0, S::from_i32(0));
    matrix.set_x(0, 0, inf);
    matrix.set_y(0, 0, inf);

    for i in 1..m {
        matrix.set_m(i, 0, inf);
        matrix.set_x(i, 0, gap_open.saturating_add(S::from_i32(i as i32 * config.gap_extend_penalty)));
        matrix.set_y(i, 0, inf);
        matrix.set_trace(i, 0, TraceDirection::new_insert());
    }
    for j in 1..n {
        matrix.set_m(0, j, inf);
        matrix.set_x(0, j, inf);
        matrix.set_y(0, j, gap_open.saturating_add(S::from_i32(j as i32 * config.gap_extend_penalty)));
        matrix.set_trace(0, j, TraceDirection::new_delete());
    }

    let use_simd = std::is_x86_feature_detected!("avx2") && n > 32 && S::KIND == ScoreKind::I32;

    for i in 1..m {
        if use_simd && S::KIND == ScoreKind::I32 {
            #[cfg(target_arch = "x86_64")]
            unsafe {
                let m32: *mut GenericDPMatrix<i32> =
                    &mut matrix as *mut GenericDPMatrix<S> as *mut GenericDPMatrix<i32>;
                crate::simd::avx2::fill_row_simd(&mut *m32, i, query, target, config);
            }
            #[cfg(not(target_arch = "x86_64"))]
            {
                fill_row_scalar::<S>(&mut matrix, i, query, target, config);
            }
        } else {
            fill_row_scalar::<S>(&mut matrix, i, query, target, config);
        }
    }

    let score_i64: i64 = matrix.get_score(m - 1, n - 1).into();
    let alignment = crate::traceback::traceback_generic(&matrix, query, target);
    AlignmentResult {
        score: score_i64,
        ..alignment
    }
}

fn fill_row_scalar<S: Score>(
    matrix: &mut GenericDPMatrix<S>,
    i: usize,
    query: &[Base],
    target: &[Base],
    config: &AlignmentConfig,
) {
    let q_base = query[i - 1];
    let cols = matrix.cols;
    let gap_open = S::from_i32(config.gap_open_penalty);
    let gap_ext = S::from_i32(config.gap_extend_penalty);
    let open_plus_ext = gap_open.saturating_add(gap_ext);

    for j in 1..cols {
        let sp = score_pair(q_base, target[j - 1], config);
        let s = S::from_i32(sp);

        let m_prev = matrix.get_m(i - 1, j - 1);
        let x_prev = matrix.get_x(i - 1, j - 1);
        let y_prev = matrix.get_y(i - 1, j - 1);
        let match_from = if m_prev >= x_prev && m_prev >= y_prev {
            m_prev
        } else if x_prev >= y_prev {
            x_prev
        } else {
            y_prev
        };
        let m_new = match_from.saturating_add(s);

        let m_up = matrix.get_m(i - 1, j).saturating_add(open_plus_ext);
        let x_up = matrix.get_x(i - 1, j).saturating_add(gap_ext);
        let x_new = if m_up >= x_up { m_up } else { x_up };

        let m_left = matrix.get_m(i, j - 1).saturating_add(open_plus_ext);
        let y_left = matrix.get_y(i, j - 1).saturating_add(gap_ext);
        let y_new = if m_left >= y_left { m_left } else { y_left };

        let best_state = if m_new >= x_new && m_new >= y_new {
            AffineState::Match
        } else if x_new >= y_new {
            AffineState::Insert
        } else {
            AffineState::Delete
        };

        matrix.set_m(i, j, m_new);
        matrix.set_x(i, j, x_new);
        matrix.set_y(i, j, y_new);
        matrix.set_trace(i, j, TraceDirection { state: best_state });
    }
}

pub mod hirschberg;
