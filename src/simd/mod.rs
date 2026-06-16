#[cfg(target_arch = "x86_64")]
pub mod avx2 {
    use crate::align::{score_pair, AffineState, AlignmentConfig, GenericDPMatrix, TraceDirection};
    use crate::fasta::Base;
    use std::arch::x86_64::*;

    #[inline]
    #[target_feature(enable = "avx2")]
    pub unsafe fn fill_row_simd(
        matrix: &mut GenericDPMatrix<i32>,
        i: usize,
        query: &[Base],
        target: &[Base],
        config: &AlignmentConfig,
    ) {
        let q_base = query[i - 1];
        let gap_open = config.gap_open_penalty;
        let gap_ext = config.gap_extend_penalty;
        let open_plus_ext = gap_open + gap_ext;
        let cols = matrix.cols;

        let gap_ext_vec = _mm256_set1_epi32(gap_ext);
        let open_plus_ext_vec = _mm256_set1_epi32(open_plus_ext);

        let simd_width = 8usize;
        let j_end = cols;
        let j_simd_end = 1 + ((j_end - 1) / simd_width) * simd_width;

        let row_m_ptr = matrix.mat_m.as_mut_ptr().add(i * cols);
        let row_x_ptr = matrix.mat_x.as_mut_ptr().add(i * cols);
        let row_y_ptr = matrix.mat_y.as_mut_ptr().add(i * cols);
        let row_trace_ptr = matrix.traceback.as_mut_ptr().add(i * cols);

        let prev_m_ptr = matrix.mat_m.as_ptr().add((i - 1) * cols);
        let prev_x_ptr = matrix.mat_x.as_ptr().add((i - 1) * cols);

        let mut j = 1usize;
        while j < j_simd_end {
            let prev_m_diag = _mm256_loadu_si256(prev_m_ptr.add(j - 1) as *const __m256i);
            let prev_x_diag = _mm256_loadu_si256(prev_x_ptr.add(j - 1) as *const __m256i);
            let prev_y_diag = _mm256_loadu_si256(matrix.mat_y.as_ptr().add((i - 1) * cols).add(j - 1) as *const __m256i);

            let max_mxy = {
                let m_vs_x = _mm256_cmpgt_epi32(prev_m_diag, prev_x_diag);
                let best_mx = _mm256_blendv_epi8(prev_x_diag, prev_m_diag, m_vs_x);
                let mx_vs_y = _mm256_cmpgt_epi32(best_mx, prev_y_diag);
                _mm256_blendv_epi8(prev_y_diag, best_mx, mx_vs_y)
            };

            let mut sub_scores = [0i32; 8];
            for k in 0..8 {
                sub_scores[k] = score_pair(q_base, target[j + k - 1], config);
            }
            let sub_vec = _mm256_loadu_si256(sub_scores.as_ptr() as *const __m256i);

            let m_new_vec = _mm256_add_epi32(max_mxy, sub_vec);

            let prev_m_up = _mm256_loadu_si256(prev_m_ptr.add(j) as *const __m256i);
            let prev_x_up = _mm256_loadu_si256(prev_x_ptr.add(j) as *const __m256i);
            let m_plus_open_ext_up = _mm256_add_epi32(prev_m_up, open_plus_ext_vec);
            let x_plus_ext_up = _mm256_add_epi32(prev_x_up, gap_ext_vec);
            let x_up_better = _mm256_cmpgt_epi32(m_plus_open_ext_up, x_plus_ext_up);
            let x_new_vec = _mm256_blendv_epi8(x_plus_ext_up, m_plus_open_ext_up, x_up_better);

            _mm256_storeu_si256(row_m_ptr.add(j) as *mut __m256i, m_new_vec);
            _mm256_storeu_si256(row_x_ptr.add(j) as *mut __m256i, x_new_vec);

            let mut y_left_arr = [0i32; 8];
            for k in 0..8 {
                let jk = j + k;
                let cand_m_y = (*row_m_ptr.add(jk - 1)).saturating_add(open_plus_ext);
                let cand_y_y = (*row_y_ptr.add(jk - 1)).saturating_add(gap_ext);
                let y_val = if cand_m_y >= cand_y_y { cand_m_y } else { cand_y_y };
                *row_y_ptr.add(jk) = y_val;
                y_left_arr[k] = y_val;
            }

            let m_new_arr: &mut [i32; 8] = &mut [0i32; 8];
            let x_new_arr: &mut [i32; 8] = &mut [0i32; 8];
            _mm256_storeu_si256(m_new_arr.as_mut_ptr() as *mut __m256i, m_new_vec);
            _mm256_storeu_si256(x_new_arr.as_mut_ptr() as *mut __m256i, x_new_vec);

            for k in 0..8 {
                let jk = j + k;
                let m_v = m_new_arr[k];
                let x_v = x_new_arr[k];
                let y_v = y_left_arr[k];
                let best_state = if m_v >= x_v && m_v >= y_v {
                    AffineState::Match
                } else if x_v >= y_v {
                    AffineState::Insert
                } else {
                    AffineState::Delete
                };
                *row_trace_ptr.add(jk) = TraceDirection { state: best_state };
            }

            j += simd_width;
        }

        while j < j_end {
            let sp = score_pair(q_base, target[j - 1], config);

            let m_prev = *prev_m_ptr.add(j - 1);
            let x_prev = *prev_x_ptr.add(j - 1);
            let y_prev = *matrix.mat_y.as_ptr().add((i - 1) * cols).add(j - 1);
            let match_from = if m_prev >= x_prev && m_prev >= y_prev {
                m_prev
            } else if x_prev >= y_prev {
                x_prev
            } else {
                y_prev
            };
            let m_new = match_from.saturating_add(sp);

            let m_up = (*prev_m_ptr.add(j)).saturating_add(open_plus_ext);
            let x_up = (*prev_x_ptr.add(j)).saturating_add(gap_ext);
            let x_new = if m_up >= x_up { m_up } else { x_up };

            let m_left = (*row_m_ptr.add(j - 1)).saturating_add(open_plus_ext);
            let y_left = (*row_y_ptr.add(j - 1)).saturating_add(gap_ext);
            let y_new = if m_left >= y_left { m_left } else { y_left };

            let best_state = if m_new >= x_new && m_new >= y_new {
                AffineState::Match
            } else if x_new >= y_new {
                AffineState::Insert
            } else {
                AffineState::Delete
            };

            *row_m_ptr.add(j) = m_new;
            *row_x_ptr.add(j) = x_new;
            *row_y_ptr.add(j) = y_new;
            *row_trace_ptr.add(j) = TraceDirection { state: best_state };

            j += 1;
        }
    }
}

pub fn has_avx2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}
