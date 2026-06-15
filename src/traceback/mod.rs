use crate::align::{AlignmentResult, DPMatrix, TraceDirection};
use crate::fasta::Base;

pub fn traceback(matrix: &DPMatrix, query: &[Base], target: &[Base]) -> AlignmentResult {
    let mut i = matrix.rows - 1;
    let mut j = matrix.cols - 1;

    let mut query_aligned = Vec::new();
    let mut target_aligned = Vec::new();

    let mut mismatches = 0u32;
    let mut insertions = 0u32;
    let mut deletions = 0u32;

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
                let q_base = query[i - 1];
                let t_base = target[j - 1];
                query_aligned.push(q_base.to_byte());
                target_aligned.push(t_base.to_byte());
                if q_base != t_base {
                    mismatches += 1;
                }
                i -= 1;
                j -= 1;
            }
            TraceDirection::Up => {
                query_aligned.push(query[i - 1].to_byte());
                target_aligned.push(b'-');
                insertions += 1;
                i -= 1;
            }
            TraceDirection::Left => {
                query_aligned.push(b'-');
                target_aligned.push(target[j - 1].to_byte());
                deletions += 1;
                j -= 1;
            }
        }
    }

    query_aligned.reverse();
    target_aligned.reverse();

    let cigar = build_cigar(&query_aligned, &target_aligned);
    let score = matrix.get_score(matrix.rows - 1, matrix.cols - 1);

    AlignmentResult {
        score,
        query_aligned,
        target_aligned,
        cigar,
        mismatches,
        insertions,
        deletions,
    }
}

fn build_cigar(query_aligned: &[u8], target_aligned: &[u8]) -> String {
    let mut cigar = String::new();
    let mut count: usize = 0;
    let mut last_op: Option<char> = None;

    for (q, t) in query_aligned.iter().zip(target_aligned.iter()) {
        let op = match (*q, *t) {
            (b'-', _) => 'I',
            (_, b'-') => 'D',
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
