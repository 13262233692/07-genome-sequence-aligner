use genome_aligner::fasta::{Base, PackedSequence};
use genome_aligner::align::{AlignmentConfig, needleman_wunsch, score_pair, ScoreKind};

#[test]
fn test_base_from_byte() {
    assert_eq!(Base::from_byte(b'A'), Some(Base::A));
    assert_eq!(Base::from_byte(b'a'), Some(Base::A));
    assert_eq!(Base::from_byte(b'T'), Some(Base::T));
    assert_eq!(Base::from_byte(b'C'), Some(Base::C));
    assert_eq!(Base::from_byte(b'G'), Some(Base::G));
    assert_eq!(Base::from_byte(b'N'), Some(Base::N));
    assert_eq!(Base::from_byte(b'R'), Some(Base::N));
    assert_eq!(Base::from_byte(b'X'), None);
}

#[test]
fn test_packed_sequence_roundtrip() {
    let mut seq = PackedSequence::new("test".to_string());
    let bases = [Base::A, Base::T, Base::C, Base::G, Base::N, Base::A, Base::T];
    for &b in &bases {
        seq.push(b);
    }
    assert_eq!(seq.len(), 7);
    for (i, &expected) in bases.iter().enumerate() {
        assert_eq!(seq.get(i), Some(expected), "Mismatch at index {}", i);
    }
    assert_eq!(seq.get(7), None);
}

#[test]
fn test_packed_sequence_large() {
    let mut seq = PackedSequence::new("large".to_string());
    let pattern = [Base::A, Base::T, Base::C, Base::G, Base::N];
    for _ in 0..1000 {
        for &b in &pattern {
            seq.push(b);
        }
    }
    assert_eq!(seq.len(), 5000);
    for i in 0..5000 {
        assert_eq!(seq.get(i), Some(pattern[i % 5]), "Mismatch at index {}", i);
    }
}

#[test]
fn test_score_pair() {
    let config = AlignmentConfig::default();
    assert_eq!(score_pair(Base::A, Base::A, &config), 2);
    assert_eq!(score_pair(Base::A, Base::T, &config), -3);
    assert_eq!(score_pair(Base::N, Base::A, &config), -1);
    assert_eq!(score_pair(Base::A, Base::N, &config), -1);
}

#[test]
fn test_perfect_alignment() {
    let query = vec![Base::A, Base::T, Base::C, Base::G];
    let target = vec![Base::A, Base::T, Base::C, Base::G];
    let config = AlignmentConfig::default();
    let result = needleman_wunsch(&query, &target, &config);
    assert_eq!(result.score, 8);
    assert_eq!(result.mismatches, 0);
    assert_eq!(result.insertions, 0);
    assert_eq!(result.deletions, 0);
    assert_eq!(result.cigar, "4M");
}

#[test]
fn test_single_mismatch() {
    let query = vec![Base::A, Base::T, Base::C, Base::G];
    let target = vec![Base::A, Base::T, Base::G, Base::G];
    let config = AlignmentConfig::default();
    let result = needleman_wunsch(&query, &target, &config);
    assert_eq!(result.mismatches, 1);
    assert_eq!(result.cigar, "4M");
}

#[test]
fn test_insertion() {
    let query = vec![Base::A, Base::A, Base::T, Base::C, Base::G];
    let target = vec![Base::A, Base::T, Base::C, Base::G];
    let config = AlignmentConfig::default();
    let result = needleman_wunsch(&query, &target, &config);
    assert!(result.insertions >= 1);
}

#[test]
fn test_deletion() {
    let query = vec![Base::A, Base::T, Base::C, Base::G];
    let target = vec![Base::A, Base::A, Base::T, Base::C, Base::G];
    let config = AlignmentConfig::default();
    let result = needleman_wunsch(&query, &target, &config);
    assert!(result.deletions >= 1);
}

#[test]
fn test_n_base_alignment() {
    let query = vec![Base::A, Base::N, Base::C, Base::G];
    let target = vec![Base::A, Base::T, Base::C, Base::G];
    let config = AlignmentConfig::default();
    let result = needleman_wunsch(&query, &target, &config);
    assert_eq!(result.mismatches, 1);
}

#[test]
fn test_alignment_result_aligned_sequences() {
    let query = vec![Base::A, Base::T, Base::C, Base::G];
    let target = vec![Base::A, Base::T, Base::C, Base::G];
    let config = AlignmentConfig::default();
    let result = needleman_wunsch(&query, &target, &config);
    assert_eq!(result.query_aligned, vec![b'A', b'T', b'C', b'G']);
    assert_eq!(result.target_aligned, vec![b'A', b'T', b'C', b'G']);
}

#[test]
fn test_score_kind_selection_i16() {
    let config = AlignmentConfig::default();
    let kind = config.required_score_kind(5000, 5000);
    assert_eq!(kind, ScoreKind::I16);
}

#[test]
fn test_score_kind_selection_i32() {
    let config = AlignmentConfig::default();
    let kind = config.required_score_kind(2_000_000, 2_000_000);
    assert_eq!(kind, ScoreKind::I32);
}

#[test]
fn test_score_kind_selection_i64() {
    let mut config = AlignmentConfig::default();
    config.match_score = 100;
    let kind = config.required_score_kind(3_000_000_000, 3_000_000_000);
    assert_eq!(kind, ScoreKind::I64);
}

#[test]
fn test_medium_sequence_alignment() {
    let seq_len = 500;
    let query: Vec<Base> = (0..seq_len)
        .map(|i| match i % 4 {
            0 => Base::A,
            1 => Base::T,
            2 => Base::C,
            _ => Base::G,
        })
        .collect();
    let mut target = query.clone();
    target[100] = Base::C;
    target[401] = Base::A;

    let config = AlignmentConfig::default();
    let result = needleman_wunsch(&query, &target, &config);

    assert_eq!(result.mismatches, 2);
    assert_eq!(result.insertions, 0);
    assert_eq!(result.deletions, 0);
    assert_eq!(result.cigar, "500M");
}

#[test]
fn test_saturating_add_no_panic_i32() {
    let big_seq: Vec<Base> = (0..2000).map(|_| Base::A).collect();
    let target = big_seq.clone();
    let mut config = AlignmentConfig::default();
    config.match_score = i16::MAX as i32;
    config.gap_extend_penalty = 0;
    config.gap_open_penalty = 0;

    let result = needleman_wunsch(&big_seq, &target, &config);
    assert!(result.mismatches == 0);
    assert!(result.score > 0);
}

#[test]
fn test_hirschberg_small_sequences() {
    let query = vec![Base::A, Base::T, Base::C, Base::G, Base::A, Base::T, Base::C, Base::G];
    let mut target = query.clone();
    target[3] = Base::A;

    let config = AlignmentConfig::default();
    let result = genome_aligner::align::hirschberg::hirschberg_align(&query, &target, &config);

    assert_eq!(result.mismatches, 1);
    assert_eq!(result.insertions, 0);
    assert_eq!(result.deletions, 0);
    assert_eq!(result.cigar, "8M");
}

#[test]
fn test_hirschberg_with_indels() {
    let query = vec![Base::A, Base::A, Base::T, Base::C, Base::G];
    let target = vec![Base::A, Base::T, Base::C, Base::G];
    let config = AlignmentConfig::default();
    let result = genome_aligner::align::hirschberg::hirschberg_align(&query, &target, &config);
    assert!(result.insertions >= 1);
}

#[test]
fn test_hirschberg_perfect_match() {
    let seq: Vec<Base> = (0..500)
        .map(|i| match i % 4 {
            0 => Base::A,
            1 => Base::T,
            2 => Base::C,
            _ => Base::G,
        })
        .collect();
    let target = seq.clone();
    let config = AlignmentConfig::default();
    let result = genome_aligner::align::hirschberg::hirschberg_align(&seq, &target, &config);

    assert_eq!(result.mismatches, 0);
    assert_eq!(result.insertions, 0);
    assert_eq!(result.deletions, 0);
}

#[test]
fn test_hirschberg_vs_standard_consistency() {
    let query: Vec<Base> = (0..100)
        .map(|i| match i % 4 {
            0 => Base::A,
            1 => Base::T,
            2 => Base::C,
            _ => Base::G,
        })
        .collect();
    let mut target = query.clone();
    for i in 0..10 {
        let idx = i * 9;
        target[idx] = match target[idx] {
            Base::A => Base::T,
            Base::T => Base::C,
            Base::C => Base::G,
            Base::G => Base::A,
            Base::N => Base::N,
        };
    }

    let config = AlignmentConfig::default();

    let standard = needleman_wunsch(&query, &target, &config);
    let hirsch = genome_aligner::align::hirschberg::hirschberg_align(&query, &target, &config);

    assert_eq!(standard.mismatches, hirsch.mismatches);
    assert_eq!(standard.insertions, hirsch.insertions);
    assert_eq!(standard.deletions, hirsch.deletions);
    assert_eq!(standard.cigar, hirsch.cigar);
    assert_eq!(standard.score, hirsch.score);
}

#[test]
fn test_boundary_checked_mul_no_overflow() {
    let a: usize = usize::MAX / 2;
    let b: usize = 4;
    let c = a.checked_mul(b);
    assert_eq!(c, None);

    let m = 1000usize;
    let n = 1000usize;
    let total = m.checked_mul(n).unwrap_or(usize::MAX);
    assert_eq!(total, 1_000_000);
}

#[test]
fn test_large_alignment_saturating_mismatch_tolerance() {
    let n = 3000usize;
    let query: Vec<Base> = (0..n).map(|i| {
        if i % 3 == 0 { Base::A }
        else if i % 3 == 1 { Base::T }
        else { Base::C }
    }).collect();
    let target: Vec<Base> = (0..n).map(|i| {
        if i % 3 == 0 { Base::G }
        else if i % 3 == 1 { Base::C }
        else { Base::T }
    }).collect();

    let mut config = AlignmentConfig::default();
    config.match_score = 1;
    config.mismatch_penalty = -1;
    config.gap_extend_penalty = -2;
    config.gap_open_penalty = -2;

    let result = needleman_wunsch(&query, &target, &config);
    assert!(result.mismatches > 0);
    assert!(result.cigar.len() > 0);
}

#[test]
fn test_affine_gap_open_vs_extend() {
    let query = vec![Base::A, Base::A, Base::A, Base::A, Base::T, Base::C, Base::G];
    let target = vec![Base::A, Base::T, Base::C, Base::G];

    let mut config = AlignmentConfig::default();
    config.match_score = 2;
    config.mismatch_penalty = -3;
    config.gap_open_penalty = -10;
    config.gap_extend_penalty = -1;

    let result = needleman_wunsch(&query, &target, &config);
    assert_eq!(result.insertions, 3);
    assert_eq!(result.cigar, "3I4M");

    let expected_score: i64 = 4 * 2 + (-10) + 3 * (-1);
    assert_eq!(result.score, expected_score);
}

#[test]
fn test_affine_multiple_gaps_heavier_open_penalty() {
    let query = vec![Base::A, Base::T, Base::A, Base::A, Base::C, Base::G];
    let target = vec![Base::A, Base::T, Base::C, Base::G];

    let mut config = AlignmentConfig::default();
    config.match_score = 2;
    config.mismatch_penalty = -3;
    config.gap_open_penalty = -100;
    config.gap_extend_penalty = -1;

    let result = needleman_wunsch(&query, &target, &config);
    assert_eq!(result.insertions, 2);
    assert_eq!(result.deletions, 0);

    let gap_cost: i64 = -100 + 2 * (-1);
    let expected_score: i64 = 4 * 2 + gap_cost;
    assert_eq!(result.score, expected_score);
}

#[test]
fn test_affine_vs_linear_different_result() {
    let query = vec![Base::A, Base::A, Base::A, Base::T, Base::C, Base::G];
    let target = vec![Base::A, Base::T, Base::C, Base::G];

    let mut affine_heavy = AlignmentConfig::default();
    affine_heavy.match_score = 2;
    affine_heavy.mismatch_penalty = -3;
    affine_heavy.gap_open_penalty = -20;
    affine_heavy.gap_extend_penalty = -1;

    let mut affine_light = AlignmentConfig::default();
    affine_light.match_score = 2;
    affine_light.mismatch_penalty = -3;
    affine_light.gap_open_penalty = -1;
    affine_light.gap_extend_penalty = -1;

    let r_heavy = needleman_wunsch(&query, &target, &affine_heavy);
    let r_light = needleman_wunsch(&query, &target, &affine_light);

    assert!(r_heavy.insertions >= 2);
    assert!(r_light.insertions >= 2);
    assert!(r_heavy.score < r_light.score);
}
