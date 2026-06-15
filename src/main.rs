use clap::Parser;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

use genome_aligner::fasta::{Base, FastaStream, SequenceRecord};
use genome_aligner::align::{AlignmentConfig, needleman_wunsch};
use genome_aligner::sam::write_sam_output;

#[derive(Parser, Debug)]
#[command(name = "genome-aligner", version, about = "High-performance genome sequence alignment tool")]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(clap::Subcommand, Debug)]
enum Command {
    Align {
        #[arg(short, long, help = "Path to query FASTA file")]
        query: PathBuf,
        #[arg(short, long, help = "Path to target/reference FASTA file")]
        target: PathBuf,
        #[arg(short = 'M', long, default_value = "2", help = "Match score")]
        match_score: i32,
        #[arg(short = 'x', long, default_value = "-3", help = "Mismatch penalty")]
        mismatch_penalty: i32,
        #[arg(short = 'G', long, default_value = "-5", help = "Gap open penalty")]
        gap_open: i32,
        #[arg(short = 'E', long, default_value = "-2", help = "Gap extend penalty")]
        gap_extend: i32,
        #[arg(short = 'O', long, help = "Output SAM file path (stdout if omitted)")]
        output: Option<PathBuf>,
        #[arg(long, default_value = "1048576", help = "Chunk size for streaming reads")]
        chunk_size: usize,
    },
    Info {
        #[arg(short, long, help = "Path to FASTA file")]
        input: PathBuf,
    },
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    match args.command {
        Command::Align {
            query,
            target,
            match_score,
            mismatch_penalty,
            gap_open,
            gap_extend,
            output,
            chunk_size,
        } => {
            let config = AlignmentConfig {
                match_score,
                mismatch_penalty,
                gap_open_penalty: gap_open,
                gap_extend_penalty: gap_extend,
            };

            eprintln!("Loading query sequences from {:?}", query);
            let query_records = load_all_records(&query, chunk_size)?;
            eprintln!("Loaded {} query sequence(s)", query_records.len());

            eprintln!("Loading target sequences from {:?}", target);
            let target_records = load_all_records(&target, chunk_size)?;
            eprintln!("Loaded {} target sequence(s)", target_records.len());

            let writer: Box<dyn Write> = match output {
                Some(path) => Box::new(BufWriter::new(std::fs::File::create(path)?)),
                None => Box::new(BufWriter::new(io::stdout())),
            };
            let mut writer = writer;

            let simd_available = genome_aligner::simd::has_avx2();
            if simd_available {
                eprintln!("AVX2 SIMD acceleration: ENABLED");
            } else {
                eprintln!("AVX2 SIMD acceleration: not available, using scalar fallback");
            }

            for qrec in &query_records {
                for trec in &target_records {
                    eprintln!(
                        "Aligning {} ({}bp) vs {} ({}bp)...",
                        qrec.name,
                        qrec.sequence.len(),
                        trec.name,
                        trec.sequence.len()
                    );

                    let query_bases: Vec<Base> = qrec.sequence.iter().collect();
                    let target_bases: Vec<Base> = trec.sequence.iter().collect();

                    let result = needleman_wunsch(&query_bases, &target_bases, &config);

                    eprintln!(
                        "  Score: {} | Mismatches: {} | Insertions: {} | Deletions: {} | CIGAR: {}",
                        result.score, result.mismatches, result.insertions, result.deletions, result.cigar
                    );

                    write_sam_output(&qrec.sequence, &trec.sequence, &result, &mut writer)?;
                }
            }

            writer.flush()?;
            eprintln!("Alignment complete.");
        }
        Command::Info { input } => {
            let records = load_all_records(&input, 1048576)?;
            println!("FASTA file: {:?}", input);
            println!("Number of sequences: {}", records.len());
            for rec in &records {
                println!(
                    "  {} : {} bp",
                    rec.name,
                    rec.sequence.len()
                );
            }
        }
    }

    Ok(())
}

fn load_all_records(path: &PathBuf, chunk_size: usize) -> io::Result<Vec<SequenceRecord>> {
    let stream = FastaStream::from_file(path, chunk_size)?;
    let mut records = Vec::new();
    for result in stream {
        let rec = result?;
        records.push(rec);
    }
    Ok(records)
}
