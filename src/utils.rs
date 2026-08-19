use crate::types::*;

use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};


/**
 * Get all neighbors (kmers with edit distance 1) of a given kmer value
 * Returns a hashmap of neighbor kmer value to NeighborInfo
 */
pub fn _get_neighbors(value: u64, value_size: u8, bidirectional: bool) -> HashMap<u64, NeighborInfo> {
    // get all the values with edit distance 1 from the input value

    let mut neighbors: HashMap<u64, NeighborInfo> = HashMap::new();
    let bases = [0, 1, 2, 3]; // A, C, G, T

    for i in 0..value_size {
        let shift = i * 2;
        let previous_base: u8 = if i == value_size - 1 {
            4 // N (unknown)
        } else {
            ((value >> (shift + 2)) & 0b11) as u8
        };
        let next_base: u8 = if i == 0 {
            4 // N (unknown)
        } else {
            ((value >> (shift - 2)) & 0b11) as u8
        };

        // Substitutions
        for &b in &bases {
            let current_base = (value >> shift) & 0b11;

            if b != current_base {
                let neighbor = (value & !(0b11 << shift)) | (b << shift);
                neighbors.insert(neighbor, NeighborInfo {
                    op: BASES_TO_SUBSTITUTION[current_base as usize][b as usize].unwrap(),
                    prev_base: previous_base,
                    next_base,
                    position: i,
                });
            }
        }

        // Indels
        for &b in &bases {
            if shift == 0 && b == (value >> shift) & 0b11 {
                continue; // skip the original base for the first position
            }

            let left_part = (value >> (shift + 2)) << ((shift + 2));
            let right_part = value & ((1 << (shift + 2)) - 1);
            let neighbor_insert = left_part | (b << shift) | (right_part >> 2);
            neighbors.entry(neighbor_insert)
                .and_modify(|info| {
                    if info.op != BASES_TO_INSERTION[b as usize].unwrap() {
                        info.op = EditOperation::AMBIGUOUS
                    }
                })
                .or_insert(NeighborInfo {
                    op: BASES_TO_INSERTION[b as usize].unwrap(),
                    prev_base: previous_base,
                    next_base,
                    position: i,
                });



            let right_part = value & ((1 << shift) - 1);
            let neighbor_delete = left_part | (right_part << 2) | b;
            let original_base = (value >> shift) & 0b11;
            neighbors.entry(neighbor_delete)
                .and_modify(|info|
                    if info.op != BASES_TO_DELETION[original_base as usize].unwrap() {
                        info.op = EditOperation::AMBIGUOUS
                    }
                )
                .or_insert(NeighborInfo {
                    op: BASES_TO_DELETION[original_base as usize].unwrap(),
                    prev_base: previous_base,
                    next_base,
                    position: i,
                });
        }
    }

    neighbors
}

pub fn _kmer_to_string(kmer: u64, k: u8) -> String {
    // for debugging: convert a kmer to a string

    let mut s = Vec::with_capacity(k as usize);
    for i in (0..k).rev() {
        let shift = i * 2;
        let base = ((kmer >> shift) & 0b11) as usize;
        s.push(crate::types::SEQ_TO_BYTE[base]);
    }
    String::from_utf8(s).unwrap()
}

pub fn _show_neighbors(kmer: u64, k: u8, bidirectional: bool) {
    // for debugging: print all the neighbors of a value

    let neighbors = _get_neighbors(kmer, k, bidirectional);
    for (neighbor, info) in neighbors {
        println!("Neighbor: {}, Operation: {}", _kmer_to_string(neighbor, k), sbs96_str(&(info.op, info.prev_base, info.next_base)));
    }
}

pub fn is_fastx_file(file_path: &str) -> bool {
    // Check if a file is in FASTA or FASTQ format based on its extension
    let lower_path = file_path.to_lowercase();
    let fastx_extensions = [".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz",
                            ".fq", ".fnq", ".fastq", ".fq.gz", ".fnq.gz", ".fastq.gz", ".bam"];
    fastx_extensions.iter().any(|ext| lower_path.ends_with(ext))
}

fn fastx_size_multiplier(file_path: &str) -> u64 {
    const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

    let mut file = match File::open(file_path) {
        Ok(file) => file,
        Err(_) => return 1,
    };
    let mut first_two_bytes = [0u8; 2];
    if file.read_exact(&mut first_two_bytes).is_err() {
        return 1;
    }

    let is_gzip = first_two_bytes == GZIP_MAGIC;
    let first_fastx_byte = if is_gzip {
        if file.seek(SeekFrom::Start(0)).is_err() {
            return 1;
        }
        let mut decoder = MultiGzDecoder::new(file);
        let mut first_byte = [0u8; 1];
        if decoder.read_exact(&mut first_byte).is_err() {
            return 1;
        }
        first_byte[0]
    } else {
        first_two_bytes[0]
    };

    match (first_fastx_byte, is_gzip) {
        (b'@', false) => 1, // FASTQ
        (b'@', true) => 4,  // FASTQ.GZ
        (b'>', false) => 2, // FASTA
        (b'>', true) => 8,  // FASTA.GZ
        _ => 1,
    }
}

/**
 * Estimate a suitable subsampling rate `-c` from raw sequencing input files.
 * File types are detected from their contents rather than their extensions.
 * Disk sizes are weighted by 1x for FASTQ, 4x for gzipped FASTQ, 2x for
 * FASTA, and 8x for gzipped FASTA.
 * Returns ceiling(total_weighted_size / 10 GiB * 1000).
 *
 * This is chosen so that the number of sketched (k,v)-mers is around 10M, for
 * efficient loading and in-memory processing.
 *
 * Returns (used_c, total_weighted_size).
 */
pub fn estimate_c_from_raw_files(files: &[&str]) -> (usize, u64) {
    const TEN_GB: u64 = 10 * 1024 * 1024 * 1024;

    let total_size = files.iter().fold(0u64, |total, file_path| {
        let size = std::fs::metadata(file_path).map(|metadata| metadata.len()).unwrap_or(0);
        total.saturating_add(size.saturating_mul(fastx_size_multiplier(file_path)))
    });

    if total_size == 0 {
        return (1000, 0);
    }

    // let chunks = total_size.div_ceil(SIXTEEN_GB);
    // ((chunks as usize) * 1000, total_size)

    (((total_size as f64 / TEN_GB as f64) * 1000.).ceil() as usize, total_size)
}
