use crate::align::{AlignmentConfig, AlignmentResult, Score, ScoreKind, score_pair};
use crate::fasta::Base;
use crate::traceback::{EditOp, build_result_from_edits};

fn last_col<S: Score>(query: &[Base], target: &[Base], config: &AlignmentConfig) -> Vec<S> {
    let m = query.len();
    let n = target.len();
    let gap_ext = config.gap_extend_penalty;
    let gap_open = config.gap_open_penalty;
    let zero = S::from_i32(0);

    let mut prev = vec![zero; n + 1];
    let mut curr = vec![zero; n + 1];

    for j in 1..=n {
        prev[j] = S::from_i32(gap_open + (j as i32 - 1) * gap_ext);
    }

    for i in 1..=m {
        curr[0] = S::from_i32(gap_open + (i as i32 - 1) * gap_ext);
        let q_base = query[i - 1];

        for j in 1..=n {
            let sp = score_pair(q_base, target[j - 1], config);
            let diag = prev[j - 1].saturating_add(S::from_i32(sp));
            let up = prev[j].saturating_add(S::from_i32(gap_ext));
            let left = curr[j - 1].saturating_add(S::from_i32(gap_ext));

            curr[j] = if diag >= up && diag >= left {
                diag
            } else if up >= left {
                up
            } else {
                left
            };
        }

        std::mem::swap(&mut prev, &mut curr);
    }

    prev
}

fn hirschberg_recurse<S: Score>(
    query: &[Base],
    target: &[Base],
    config: &AlignmentConfig,
) -> Vec<EditOp> {
    let m = query.len();
    let n = target.len();

    if m == 0 {
        let mut ops = Vec::new();
        for _ in 0..n {
            ops.push(EditOp::Deletion);
        }
        return ops;
    }
    if n == 0 {
        let mut ops = Vec::new();
        for _ in 0..m {
            ops.push(EditOp::Insertion);
        }
        return ops;
    }

    const THRESHOLD: usize = 2000;

    if m <= THRESHOLD || n <= THRESHOLD {
        return nw_small::<S>(query, target, config);
    }

    let mid = m / 2;

    let score_l = last_col::<S>(&query[..mid], target, config);

    let q_rev: Vec<Base> = query[mid..].iter().rev().cloned().collect();
    let t_rev: Vec<Base> = target.iter().rev().cloned().collect();
    let s = last_col::<S>(&q_rev, &t_rev, config);
    let score_r_rev: Vec<S> = s.into_iter().rev().collect();

    let mut best_j = 0;
    let mut best_sum = S::MIN;

    for j in 0..=n {
        let sum = score_l[j].saturating_add(score_r_rev[j]);
        if j == 0 || sum > best_sum {
            best_sum = sum;
            best_j = j;
        }
    }

    let mut left = hirschberg_recurse::<S>(&query[..mid], &target[..best_j], config);
    let mut right = hirschberg_recurse::<S>(&query[mid..], &target[best_j..], config);
    left.append(&mut right);
    left
}

fn nw_small<S: Score>(query: &[Base], target: &[Base], config: &AlignmentConfig) -> Vec<EditOp> {
    use crate::align::{GenericDPMatrix, TraceDirection};

    let m = query.len() + 1;
    let n = target.len() + 1;
    let gap_ext = config.gap_extend_penalty;
    let gap_open = config.gap_open_penalty;

    let mut matrix = GenericDPMatrix::<S>::new(m, n);

    for i in 1..m {
        matrix.set_score(
            i,
            0,
            S::from_i32(gap_open + (i as i32 - 1) * gap_ext),
        );
        matrix.set_trace(i, 0, TraceDirection::Up);
    }
    for j in 1..n {
        matrix.set_score(
            0,
            j,
            S::from_i32(gap_open + (j as i32 - 1) * gap_ext),
        );
        matrix.set_trace(0, j, TraceDirection::Left);
    }

    for i in 1..m {
        let q_base = query[i - 1];
        for j in 1..n {
            let sp = score_pair(q_base, target[j - 1], config);
            let diag = matrix.get_score(i - 1, j - 1).saturating_add(S::from_i32(sp));
            let up = matrix.get_score(i - 1, j).saturating_add(S::from_i32(gap_ext));
            let left = matrix.get_score(i, j - 1).saturating_add(S::from_i32(gap_ext));

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

    let mut ops = Vec::new();
    let mut i = m - 1;
    let mut j = n - 1;

    while i > 0 || j > 0 {
        let dir = if i > 0 && j > 0 {
            matrix.get_trace(i, j)
        } else if i > 0 {
            TraceDirection::Up
        } else {
            TraceDirection::Left
        };

        match dir {
            TraceDirection::Diagonal => {
                let qb = query[i - 1];
                let tb = target[j - 1];
                ops.push(EditOp::Match(qb == tb));
                i -= 1;
                j -= 1;
            }
            TraceDirection::Up => {
                ops.push(EditOp::Insertion);
                i -= 1;
            }
            TraceDirection::Left => {
                ops.push(EditOp::Deletion);
                j -= 1;
            }
        }
    }

    ops.reverse();
    ops
}

pub fn hirschberg_align(
    query: &[Base],
    target: &[Base],
    config: &AlignmentConfig,
) -> AlignmentResult {
    let kind = config.required_score_kind(query.len(), target.len());

    let edits = match kind {
        ScoreKind::I16 => hirschberg_recurse::<i16>(query, target, config),
        ScoreKind::I32 => hirschberg_recurse::<i32>(query, target, config),
        ScoreKind::I64 => hirschberg_recurse::<i64>(query, target, config),
    };

    let score = compute_score::<i64>(query, target, config);
    let mut result = build_result_from_edits(&edits, query, target);
    result.score = score;
    result
}

fn compute_score<S: Score>(
    query: &[Base],
    target: &[Base],
    config: &AlignmentConfig,
) -> i64 {
    let m = query.len();
    let n = target.len();
    let gap_ext = config.gap_extend_penalty;
    let gap_open = config.gap_open_penalty;
    let zero = S::from_i32(0);

    let mut prev = vec![zero; n + 1];
    let mut curr = vec![zero; n + 1];

    for j in 1..=n {
        prev[j] = S::from_i32(gap_open + (j as i32 - 1) * gap_ext);
    }

    for i in 1..=m {
        curr[0] = S::from_i32(gap_open + (i as i32 - 1) * gap_ext);
        let q_base = query[i - 1];

        for j in 1..=n {
            let sp = score_pair(q_base, target[j - 1], config);
            let diag = prev[j - 1].saturating_add(S::from_i32(sp));
            let up = prev[j].saturating_add(S::from_i32(gap_ext));
            let left = curr[j - 1].saturating_add(S::from_i32(gap_ext));

            curr[j] = if diag >= up && diag >= left {
                diag
            } else if up >= left {
                up
            } else {
                left
            };
        }

        std::mem::swap(&mut prev, &mut curr);
    }

    let final_score: i64 = prev[n].into();
    final_score
}
