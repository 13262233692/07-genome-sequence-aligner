use crate::align::{AffineState, AlignmentResult, GenericDPMatrix, Score};
use crate::fasta::Base;

pub fn traceback_generic<S: Score>(
    matrix: &GenericDPMatrix<S>,
    query: &[Base],
    target: &[Base],
) -> AlignmentResult {
    let m = matrix.rows - 1;
    let n = matrix.cols - 1;

    let end_m = matrix.get_m(m, n);
    let end_x = matrix.get_x(m, n);
    let end_y = matrix.get_y(m, n);
    let mut state = if end_m >= end_x && end_m >= end_y {
        AffineState::Match
    } else if end_x >= end_y {
        AffineState::Insert
    } else {
        AffineState::Delete
    };

    let mut i = m;
    let mut j = n;

    let mut query_aligned = Vec::new();
    let mut target_aligned = Vec::new();

    let mut mismatches = 0u32;
    let mut insertions = 0u32;
    let mut deletions = 0u32;

    while i > 0 || j > 0 {
        if i > 0 && j == 0 {
            state = AffineState::Insert;
        } else if i == 0 && j > 0 {
            state = AffineState::Delete;
        }
        match state {
            AffineState::Match => {
                let q_base = query[i - 1];
                let t_base = target[j - 1];
                query_aligned.push(q_base.to_byte());
                target_aligned.push(t_base.to_byte());
                if q_base != t_base {
                    mismatches += 1;
                }
                if i > 1 && j > 1 {
                    let prev_m = matrix.get_m(i - 1, j - 1);
                    let prev_x = matrix.get_x(i - 1, j - 1);
                    let prev_y = matrix.get_y(i - 1, j - 1);
                    state = if prev_m >= prev_x && prev_m >= prev_y {
                        AffineState::Match
                    } else if prev_x >= prev_y {
                        AffineState::Insert
                    } else {
                        AffineState::Delete
                    };
                } else if i == 1 && j == 1 {
                } else if i > 0 && j == 0 {
                    state = AffineState::Insert;
                } else {
                    state = AffineState::Delete;
                }
                i -= 1;
                j -= 1;
            }
            AffineState::Insert => {
                query_aligned.push(query[i - 1].to_byte());
                target_aligned.push(b'-');
                insertions += 1;
                if i > 1 {
                    let m_up = matrix.get_m(i - 1, j);
                    let x_up = matrix.get_x(i - 1, j);
                    if x_up >= m_up {
                        state = AffineState::Insert;
                    } else {
                        state = AffineState::Match;
                    }
                } else if j > 0 {
                    state = AffineState::Delete;
                }
                i -= 1;
            }
            AffineState::Delete => {
                query_aligned.push(b'-');
                target_aligned.push(target[j - 1].to_byte());
                deletions += 1;
                if j > 1 {
                    let m_left = matrix.get_m(i, j - 1);
                    let y_left = matrix.get_y(i, j - 1);
                    if y_left >= m_left {
                        state = AffineState::Delete;
                    } else {
                        state = AffineState::Match;
                    }
                } else if i > 0 {
                    state = AffineState::Insert;
                }
                j -= 1;
            }
        }
    }

    query_aligned.reverse();
    target_aligned.reverse();

    let cigar = build_cigar(&query_aligned, &target_aligned);
    let score_i64: i64 = matrix.get_score(matrix.rows - 1, matrix.cols - 1).into();

    AlignmentResult {
        score: score_i64,
        query_aligned,
        target_aligned,
        cigar,
        mismatches,
        insertions,
        deletions,
    }
}

pub fn traceback(
    matrix: &crate::align::DPMatrix,
    query: &[Base],
    target: &[Base],
) -> AlignmentResult {
    traceback_generic::<i32>(matrix, query, target)
}

fn build_cigar(query_aligned: &[u8], target_aligned: &[u8]) -> String {
    let mut cigar = String::new();
    let mut count: usize = 0;
    let mut last_op: Option<char> = None;

    for (q, t) in query_aligned.iter().zip(target_aligned.iter()) {
        let op = match (*q, *t) {
            (b'-', _) => 'D',
            (_, b'-') => 'I',
            _ => 'M',
        };

        if Some(op) == last_op {
            count += 1;
        } else {
            if let Some(prev) = last_op {
                cigar.push_str(&format!("{}{}", count, prev));
            }
            last_op = Some(op);
            count = 1;
        }
    }

    if let Some(prev) = last_op {
        cigar.push_str(&format!("{}{}", count, prev));
    }

    cigar
}

pub fn build_result_from_edits(
    edits: &[EditOp],
    query: &[Base],
    target: &[Base],
) -> AlignmentResult {
    let mut query_aligned = Vec::new();
    let mut target_aligned = Vec::new();
    let mut mismatches = 0u32;
    let mut insertions = 0u32;
    let mut deletions = 0u32;

    let mut qi = 0;
    let mut ti = 0;

    for op in edits {
        match op {
            EditOp::Match(b) => {
                query_aligned.push(query[qi].to_byte());
                target_aligned.push(target[ti].to_byte());
                if !*b {
                    mismatches += 1;
                }
                qi += 1;
                ti += 1;
            }
            EditOp::Insertion => {
                query_aligned.push(query[qi].to_byte());
                target_aligned.push(b'-');
                insertions += 1;
                qi += 1;
            }
            EditOp::Deletion => {
                query_aligned.push(b'-');
                target_aligned.push(target[ti].to_byte());
                deletions += 1;
                ti += 1;
            }
        }
    }

    let cigar = build_cigar(&query_aligned, &target_aligned);

    AlignmentResult {
        score: 0,
        query_aligned,
        target_aligned,
        cigar,
        mismatches,
        insertions,
        deletions,
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EditOp {
    Match(bool),
    Insertion,
    Deletion,
}
