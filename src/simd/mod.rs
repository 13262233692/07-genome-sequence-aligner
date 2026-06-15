#[cfg(target_arch = "x86_64")]
pub mod avx2 {
    use crate::align::{score_pair, AlignmentConfig, DPMatrix, TraceDirection};
    use crate::fasta::Base;
    use std::arch::x86_64::*;

    #[inline]
    #[target_feature(enable = "avx2")]
    pub unsafe fn fill_row_simd(
        matrix: &mut DPMatrix,
        i: usize,
        query: &[Base],
        target: &[Base],
        config: &AlignmentConfig,
    ) {
        let q_base = query[i - 1];
        let gap_ext = config.gap_extend_penalty;
        let gap_ext_vec = _mm256_set1_epi32(gap_ext);
        let cols = matrix.cols;

        let simd_width = 8usize;
        let j_end = cols;
        let j_simd_end = 1 + ((j_end - 1) / simd_width) * simd_width;

        let row_scores_ptr = matrix.scores.as_mut_ptr().add(i * cols);
        let prev_row_scores_ptr = matrix.scores.as_ptr().add((i - 1) * cols);
        let row_trace_ptr = matrix.traceback.as_mut_ptr().add(i * cols);

        let mut j = 1usize;
        while j < j_simd_end {
            let score_diag = _mm256_loadu_si256(prev_row_scores_ptr.add(j - 1) as *const __m256i);
            let score_up = _mm256_loadu_si256(prev_row_scores_ptr.add(j) as *const __m256i);

            let mut sub_scores = [0i32; 8];
            for k in 0..8 {
                sub_scores[k] = score_pair(q_base, target[j + k - 1], config);
            }
            let sub_vec = _mm256_loadu_si256(sub_scores.as_ptr() as *const __m256i);

            let diag_vec = _mm256_add_epi32(score_diag, sub_vec);
            let up_vec = _mm256_add_epi32(score_up, gap_ext_vec);

            let left_base = *row_scores_ptr.add(j - 1);
            let mut left_arr = [left_base; 8];
            for k in 0..7 {
                left_arr[k + 1] = left_arr[k] + gap_ext;
            }
            let left_vec = _mm256_loadu_si256(left_arr.as_ptr() as *const __m256i);

            let diag_vs_up = _mm256_cmpgt_epi32(diag_vec, up_vec);
            let best_du = _mm256_blendv_epi8(up_vec, diag_vec, diag_vs_up);

            let best_du_vs_left = _mm256_cmpgt_epi32(best_du, left_vec);
            let best = _mm256_blendv_epi8(left_vec, best_du, best_du_vs_left);

            _mm256_storeu_si256(row_scores_ptr.add(j) as *mut __m256i, best);

            let diag_vs_left = _mm256_cmpgt_epi32(diag_vec, left_vec);
            let diag_best = _mm256_and_si256(diag_vs_up, diag_vs_left);

            let not_diag_vs_up = _mm256_andnot_si256(diag_vs_up, _mm256_set1_epi32(-1));
            let up_best = _mm256_and_si256(not_diag_vs_up, _mm256_cmpgt_epi32(up_vec, left_vec));

            let mut diag_best_arr = [0i32; 8];
            let mut up_best_arr = [0i32; 8];
            _mm256_storeu_si256(diag_best_arr.as_mut_ptr() as *mut __m256i, diag_best);
            _mm256_storeu_si256(up_best_arr.as_mut_ptr() as *mut __m256i, up_best);

            for k in 0..8 {
                let d_mask = (diag_best_arr[k] >> 31) & 1;
                let u_mask = (up_best_arr[k] >> 31) & 1;
                let dir = if d_mask != 0 {
                    TraceDirection::Diagonal
                } else if u_mask != 0 {
                    TraceDirection::Up
                } else {
                    TraceDirection::Left
                };
                *row_trace_ptr.add(j + k) = dir;
            }

            j += simd_width;
        }

        while j < j_end {
            let diag = *prev_row_scores_ptr.add(j - 1) + score_pair(q_base, target[j - 1], config);
            let up = *prev_row_scores_ptr.add(j) + gap_ext;
            let left = *row_scores_ptr.add(j - 1) + gap_ext;

            let (best_score, best_dir) = if diag >= up && diag >= left {
                (diag, TraceDirection::Diagonal)
            } else if up >= left {
                (up, TraceDirection::Up)
            } else {
                (left, TraceDirection::Left)
            };

            *row_scores_ptr.add(j) = best_score;
            *row_trace_ptr.add(j) = best_dir;

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
