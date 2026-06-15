use genome_aligner::fasta::{Base, PackedSequence};
use genome_aligner::align::{AlignmentConfig, needleman_wunsch, score_pair};

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
