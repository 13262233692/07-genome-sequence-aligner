use crate::align::{AffineState, AlignmentConfig, AlignmentResult, Score, ScoreKind, score_pair, TraceDirection};
use crate::fasta::Base;
use crate::traceback::{EditOp, build_result_from_edits};

fn last_col_affine<S: Score>(query: &[Base], target: &[Base], config: &AlignmentConfig) -> (Vec<S>, Vec<S>, Vec<S>) {
    let m = query.len();
    let n = target.len();
    let gap_ext = S::from_i32(config.gap_extend_penalty);
    let gap_open = S::from_i32(config.gap_open_penalty);
    let open_plus_ext = gap_open.saturating_add(gap_ext);
    let zero = S::from_i32(0);
    let inf = S::MIN.saturating_add(S::from_i32(1000));

    let mut prev_m = vec![zero; n + 1];
    let mut prev_x = vec![zero; n + 1];
    let mut prev_y = vec![zero; n + 1];
    let mut curr_m = vec![zero; n + 1];
    let mut curr_x = vec![zero; n + 1];
    let mut curr_y = vec![zero; n + 1];

    prev_m[0] = zero;
    prev_x[0] = inf;
    prev_y[0] = inf;
    for j in 1..=n {
        prev_m[j] = inf;
        prev_x[j] = inf;
        prev_y[j] = gap_open.saturating_add(S::from_i32(j as i32 * config.gap_extend_penalty));
    }

    for i in 1..=m {
        curr_m[0] = inf;
        curr_x[0] = gap_open.saturating_add(S::from_i32(i as i32 * config.gap_extend_penalty));
        curr_y[0] = inf;
        let q_base = query[i - 1];

        for j in 1..=n {
            let sp = score_pair(q_base, target[j - 1], config);
            let s = S::from_i32(sp);

            let m_prev = prev_m[j - 1];
            let x_prev = prev_x[j - 1];
            let y_prev = prev_y[j - 1];
            let match_from = if m_prev >= x_prev && m_prev >= y_prev {
                m_prev
            } else if x_prev >= y_prev {
                x_prev
            } else {
                y_prev
            };
            curr_m[j] = match_from.saturating_add(s);

            let m_up = prev_m[j].saturating_add(open_plus_ext);
            let x_up = prev_x[j].saturating_add(gap_ext);
            curr_x[j] = if m_up >= x_up { m_up } else { x_up };

            let m_left = curr_m[j - 1].saturating_add(open_plus_ext);
            let y_left = curr_y[j - 1].saturating_add(gap_ext);
            curr_y[j] = if m_left >= y_left { m_left } else { y_left };
        }

        std::mem::swap(&mut prev_m, &mut curr_m);
        std::mem::swap(&mut prev_x, &mut curr_x);
        std::mem::swap(&mut prev_y, &mut curr_y);
    }

    (prev_m, prev_x, prev_y)
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

    let (l_m, l_x, l_y) = last_col_affine::<S>(&query[..mid], target, config);

    let q_rev: Vec<Base> = query[mid..].iter().rev().cloned().collect();
    let t_rev: Vec<Base> = target.iter().rev().cloned().collect();
    let (r_m, r_x, r_y) = last_col_affine::<S>(&q_rev, &t_rev, config);
    let r_m_rev: Vec<S> = r_m.into_iter().rev().collect();
    let r_x_rev: Vec<S> = r_x.into_iter().rev().collect();
    let r_y_rev: Vec<S> = r_y.into_iter().rev().collect();

    let mut best_j = 0;
    let mut best_sum = S::MIN;

    for j in 0..=n {
        let m_sum = l_m[j].saturating_add(r_m_rev[j]);
        let x_sum = l_x[j].saturating_add(r_x_rev[j]);
        let y_sum = l_y[j].saturating_add(r_y_rev[j]);
        let local_max = if m_sum >= x_sum && m_sum >= y_sum {
            m_sum
        } else if x_sum >= y_sum {
            x_sum
        } else {
            y_sum
        };
        if j == 0 || local_max > best_sum {
            best_sum = local_max;
            best_j = j;
        }
    }

    let mut left = hirschberg_recurse::<S>(&query[..mid], &target[..best_j], config);
    let mut right = hirschberg_recurse::<S>(&query[mid..], &target[best_j..], config);
    left.append(&mut right);
    left
}

fn nw_small<S: Score>(query: &[Base], target: &[Base], config: &AlignmentConfig) -> Vec<EditOp> {
    use crate::align::GenericDPMatrix;

    let m = query.len() + 1;
    let n = target.len() + 1;
    let gap_ext = config.gap_extend_penalty;
    let gap_open = config.gap_open_penalty;
    let open_plus_ext = gap_open + gap_ext;

    let mut matrix = GenericDPMatrix::<S>::new(m, n);
    let inf = S::MIN.saturating_add(S::from_i32(1000));

    matrix.set_m(0, 0, S::from_i32(0));
    matrix.set_x(0, 0, inf);
    matrix.set_y(0, 0, inf);

    for i in 1..m {
        matrix.set_m(i, 0, inf);
        matrix.set_x(
            i,
            0,
            S::from_i32(gap_open + i as i32 * gap_ext),
        );
        matrix.set_y(i, 0, inf);
        matrix.set_trace(i, 0, TraceDirection::new_insert());
    }
    for j in 1..n {
        matrix.set_m(0, j, inf);
        matrix.set_x(0, j, inf);
        matrix.set_y(
            0,
            j,
            S::from_i32(gap_open + j as i32 * gap_ext),
        );
        matrix.set_trace(0, j, TraceDirection::new_delete());
    }

    let gap_ext_s = S::from_i32(gap_ext);
    let open_plus_ext_s = S::from_i32(open_plus_ext);

    for i in 1..m {
        let q_base = query[i - 1];
        for j in 1..n {
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

            let m_up = matrix.get_m(i - 1, j).saturating_add(open_plus_ext_s);
            let x_up = matrix.get_x(i - 1, j).saturating_add(gap_ext_s);
            let x_new = if m_up >= x_up { m_up } else { x_up };

            let m_left = matrix.get_m(i, j - 1).saturating_add(open_plus_ext_s);
            let y_left = matrix.get_y(i, j - 1).saturating_add(gap_ext_s);
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

    let mut ops = Vec::new();
    let end_m = matrix.get_m(m - 1, n - 1);
    let end_x = matrix.get_x(m - 1, n - 1);
    let end_y = matrix.get_y(m - 1, n - 1);
    let mut state = if end_m >= end_x && end_m >= end_y {
        AffineState::Match
    } else if end_x >= end_y {
        AffineState::Insert
    } else {
        AffineState::Delete
    };

    let mut i = m - 1;
    let mut j = n - 1;

    while i > 0 || j > 0 {
        if i > 0 && j == 0 {
            state = AffineState::Insert;
        } else if i == 0 && j > 0 {
            state = AffineState::Delete;
        }
        match state {
            AffineState::Match => {
                let qb = query[i - 1];
                let tb = target[j - 1];
                ops.push(EditOp::Match(qb == tb));
                if i > 1 && j > 1 {
                    let pm = matrix.get_m(i - 1, j - 1);
                    let px = matrix.get_x(i - 1, j - 1);
                    let py = matrix.get_y(i - 1, j - 1);
                    state = if pm >= px && pm >= py {
                        AffineState::Match
                    } else if px >= py {
                        AffineState::Insert
                    } else {
                        AffineState::Delete
                    };
                } else if i == 1 && j == 0 {
                    state = AffineState::Insert;
                } else if i == 0 && j == 1 {
                    state = AffineState::Delete;
                }
                i -= 1;
                j -= 1;
            }
            AffineState::Insert => {
                ops.push(EditOp::Insertion);
                if i > 1 {
                    let pm = matrix.get_m(i - 1, j);
                    let px = matrix.get_x(i - 1, j);
                    state = if px >= pm { AffineState::Insert } else { AffineState::Match };
                } else if j > 0 {
                    state = AffineState::Delete;
                }
                i -= 1;
            }
            AffineState::Delete => {
                ops.push(EditOp::Deletion);
                if j > 1 {
                    let pm = matrix.get_m(i, j - 1);
                    let py = matrix.get_y(i, j - 1);
                    state = if py >= pm { AffineState::Delete } else { AffineState::Match };
                } else if i > 0 {
                    state = AffineState::Insert;
                }
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

    let score = compute_score_affine::<i64>(query, target, config);
    let mut result = build_result_from_edits(&edits, query, target);
    result.score = score;
    result
}

fn compute_score_affine<S: Score>(
    query: &[Base],
    target: &[Base],
    config: &AlignmentConfig,
) -> i64 {
    let (m_last, x_last, y_last) = last_col_affine::<S>(query, target, config);
    let n = target.len();
    let m_v = m_last[n];
    let x_v = x_last[n];
    let y_v = y_last[n];
    let final_score = if m_v >= x_v && m_v >= y_v {
        m_v
    } else if x_v >= y_v {
        x_v
    } else {
        y_v
    };
    let out: i64 = final_score.into();
    out
}
