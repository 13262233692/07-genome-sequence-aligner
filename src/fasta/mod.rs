use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[repr(u8)]
pub enum Base {
    A = 0,
    T = 1,
    C = 2,
    G = 3,
    N = 4,
}

impl Base {
    pub fn from_byte(b: u8) -> Option<Self> {
        match b.to_ascii_uppercase() {
            b'A' => Some(Base::A),
            b'T' => Some(Base::T),
            b'C' => Some(Base::C),
            b'G' => Some(Base::G),
            b'N' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' => {
                Some(Base::N)
            }
            _ => None,
        }
    }

    pub fn to_byte(self) -> u8 {
        match self {
            Base::A => b'A',
            Base::T => b'T',
            Base::C => b'C',
            Base::G => b'G',
            Base::N => b'N',
        }
    }

    pub fn as_u8(self) -> u8 {
        self as u8
    }
}

impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_byte() as char)
    }
}

const BASES_PER_BYTE: usize = 4;

#[derive(Clone)]
pub struct PackedSequence {
    data: Vec<u8>,
    n_bitmap: Vec<u8>,
    len: usize,
    name: String,
}

impl PackedSequence {
    pub fn new(name: String) -> Self {
        PackedSequence {
            data: Vec::new(),
            n_bitmap: Vec::new(),
            len: 0,
            name,
        }
    }

    pub fn with_capacity(name: String, capacity: usize) -> Self {
        let byte_cap = (capacity + BASES_PER_BYTE - 1) / BASES_PER_BYTE;
        PackedSequence {
            data: Vec::with_capacity(byte_cap),
            n_bitmap: Vec::with_capacity(byte_cap),
            len: 0,
            name,
        }
    }

    pub fn push(&mut self, base: Base) {
        let byte_idx = self.len / BASES_PER_BYTE;
        let bit_offset = (self.len % BASES_PER_BYTE) * 2;

        if bit_offset == 0 {
            self.data.push(0);
            self.n_bitmap.push(0);
        }

        let val = base.as_u8() & 0b11;
        self.data[byte_idx] |= val << bit_offset;

        if base == Base::N {
            self.n_bitmap[byte_idx] |= 1u8 << bit_offset;
        }

        self.len += 1;
    }

    pub fn get(&self, idx: usize) -> Option<Base> {
        if idx >= self.len {
            return None;
        }
        let byte_idx = idx / BASES_PER_BYTE;
        let bit_offset = (idx % BASES_PER_BYTE) * 2;

        let is_n = (self.n_bitmap[byte_idx] >> bit_offset) & 0b11 != 0;
        if is_n {
            return Some(Base::N);
        }

        let val = (self.data[byte_idx] >> bit_offset) & 0b11;
        match val {
            0 => Some(Base::A),
            1 => Some(Base::T),
            2 => Some(Base::C),
            3 => Some(Base::G),
            _ => unreachable!(),
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn iter(&self) -> PackedSequenceIter<'_> {
        PackedSequenceIter { seq: self, pos: 0 }
    }

    pub fn to_string_lossy(&self) -> String {
        self.iter().map(|b| b.to_byte() as char).collect()
    }
}

pub struct PackedSequenceIter<'a> {
    seq: &'a PackedSequence,
    pos: usize,
}

impl<'a> Iterator for PackedSequenceIter<'a> {
    type Item = Base;

    fn next(&mut self) -> Option<Self::Item> {
        let base = self.seq.get(self.pos)?;
        self.pos += 1;
        Some(base)
    }
}

pub struct SequenceRecord {
    pub name: String,
    pub sequence: PackedSequence,
}

pub struct FastaReader {
    reader: BufReader<File>,
    chunk_size: usize,
    current_name: Option<String>,
    exhausted: bool,
}

impl FastaReader {
    pub fn from_file<P: AsRef<Path>>(path: P, chunk_size: usize) -> io::Result<Self> {
        let file = File::open(path)?;
        Ok(FastaReader {
            reader: BufReader::with_capacity(chunk_size * 2, file),
            chunk_size,
            current_name: None,
            exhausted: false,
        })
    }

    fn read_next_record(&mut self) -> io::Result<Option<SequenceRecord>> {
        if self.exhausted {
            return Ok(None);
        }

        let mut line = Vec::new();

        if self.current_name.is_none() {
            loop {
                line.clear();
                let bytes_read = self.reader.read_until(b'\n', &mut line)?;
                if bytes_read == 0 {
                    self.exhausted = true;
                    return Ok(None);
                }
                let line_str = String::from_utf8_lossy(&line);
                let trimmed = line_str.trim();
                if trimmed.is_empty() {
                    continue;
                }
                if trimmed.starts_with('>') {
                    let name = trimmed[1..].split_whitespace().next().unwrap_or("").to_string();
                    self.current_name = Some(name);
                    break;
                } else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "FASTA file does not start with header line",
                    ));
                }
            }
        }

        let name = self.current_name.clone().unwrap();
        let mut seq = PackedSequence::with_capacity(name.clone(), self.chunk_size);
        let mut total_bases = 0usize;

        loop {
            line.clear();
            let bytes_read = self.reader.read_until(b'\n', &mut line)?;
            if bytes_read == 0 {
                self.exhausted = true;
                self.current_name = None;
                break;
            }

            let line_str = String::from_utf8_lossy(&line);
            let trimmed = line_str.trim();

            if trimmed.is_empty() {
                continue;
            }

            if trimmed.starts_with('>') {
                let new_name = trimmed[1..].split_whitespace().next().unwrap_or("").to_string();
                self.current_name = Some(new_name);
                break;
            }

            for b in trimmed.bytes() {
                if let Some(base) = Base::from_byte(b) {
                    seq.push(base);
                    total_bases += 1;
                }
            }

            if total_bases >= self.chunk_size {
                break;
            }
        }

        if seq.is_empty() {
            return self.read_next_record();
        }

        Ok(Some(SequenceRecord { name, sequence: seq }))
    }
}

pub struct FastaStream {
    reader: FastaReader,
    done: bool,
}

impl FastaStream {
    pub fn from_file<P: AsRef<Path>>(path: P, chunk_size: usize) -> io::Result<Self> {
        let reader = FastaReader::from_file(path, chunk_size)?;
        Ok(FastaStream {
            reader,
            done: false,
        })
    }
}

impl Iterator for FastaStream {
    type Item = io::Result<SequenceRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        match self.reader.read_next_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => {
                self.done = true;
                None
            }
            Err(e) => {
                self.done = true;
                Some(Err(e))
            }
        }
    }
}

pub fn read_fasta_file<P: AsRef<Path>>(path: P) -> io::Result<Vec<SequenceRecord>> {
    let chunk_size = 1024 * 1024;
    let stream = FastaStream::from_file(path, chunk_size)?;
    let mut records = Vec::new();
    for result in stream {
        records.push(result?);
    }
    Ok(records)
}
