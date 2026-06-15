pub mod fasta;
pub mod align;
pub mod simd;
pub mod traceback;
pub mod sam;

pub use align::{AlignmentConfig, AlignmentResult, needleman_wunsch};
pub use fasta::{Base, PackedSequence, SequenceRecord, FastaStream, read_fasta_file};
pub use sam::{SamHeader, SamRecord, write_sam_output};
