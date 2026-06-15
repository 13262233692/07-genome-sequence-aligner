use crate::align::AlignmentResult;
use crate::fasta::PackedSequence;
use std::io::{self, Write};

pub struct SamRecord {
    pub qname: String,
    pub flag: u16,
    pub rname: String,
    pub pos: u32,
    pub mapq: u8,
    pub cigar: String,
    pub rnext: String,
    pub pnext: u32,
    pub tlen: i32,
    pub seq: String,
    pub qual: String,
}

impl SamRecord {
    pub fn from_alignment(
        query_name: &str,
        target_name: &str,
        result: &AlignmentResult,
    ) -> Self {
        let seq: String = result
            .query_aligned
            .iter()
            .map(|&b| if b == b'-' { '.' as char } else { b as char })
            .collect();

        let qual = "*".to_string();

        SamRecord {
            qname: query_name.to_string(),
            flag: 0,
            rname: target_name.to_string(),
            pos: 1,
            mapq: 60,
            cigar: result.cigar.clone(),
            rnext: "*".to_string(),
            pnext: 0,
            tlen: 0,
            seq,
            qual,
        }
    }

    pub fn write(&self, writer: &mut dyn Write) -> io::Result<()> {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.qname,
            self.flag,
            self.rname,
            self.pos,
            self.mapq,
            self.cigar,
            self.rnext,
            self.pnext,
            self.tlen,
            self.seq,
            self.qual
        )
    }
}

pub struct SamHeader {
    pub version: String,
    pub references: Vec<(String, u64)>,
}

impl SamHeader {
    pub fn new() -> Self {
        SamHeader {
            version: "1.6".to_string(),
            references: Vec::new(),
        }
    }

    pub fn add_reference(&mut self, name: String, length: u64) {
        self.references.push((name, length));
    }

    pub fn write(&self, writer: &mut dyn Write) -> io::Result<()> {
        writeln!(writer, "@HD\tVN:{}\tSO:unknown", self.version)?;
        for (name, length) in &self.references {
            writeln!(writer, "@SQ\tSN:{}\tLN:{}", name, length)?;
        }
        Ok(())
    }
}

pub fn write_sam_output(
    query: &PackedSequence,
    target: &PackedSequence,
    result: &AlignmentResult,
    writer: &mut dyn Write,
) -> io::Result<()> {
    let mut header = SamHeader::new();
    header.add_reference(target.name().to_string(), target.len() as u64);
    header.write(writer)?;

    let record = SamRecord::from_alignment(query.name(), target.name(), result);
    record.write(writer)?;

    Ok(())
}
