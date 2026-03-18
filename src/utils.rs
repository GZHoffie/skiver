use crate::types::*;

use std::collections::HashMap;


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

pub fn is_sketch_file(file_path: &str) -> bool {
    // Check if a file is a kv-mer sketch file based on its extension
    let lower_path = file_path.to_lowercase();
    lower_path.ends_with(".kvmer")
}

/**
 * Estimate a suitable subsampling rate `-c` from raw sequencing input files.
 * For .gz files, the decompressed size is estimated as 4x the compressed size.
 * Returns ceiling(total_estimated_size / 16G) * 1000, with a minimum of 1000.
 * 
 * This is chosen so that the size of the memory usage is roughly under 2GB, for 
 * efficient loading and in-memory processing.
 * 
 * Returns (used_c, total_estimated_size).
 */
pub fn estimate_c_from_raw_files(files: &[&str]) -> (usize, u64) {
    const SIXTEEN_GB: u64 = 16 * 1024 * 1024 * 1024;
    const GZ_FACTOR: u64 = 4; // Estimated decompressed size is 4x compressed size for .gz files

    let total_size: u64 = files.iter()
        .filter(|f| is_fastx_file(f))
        .map(|f| {
            let size = std::fs::metadata(f).map(|m| m.len()).unwrap_or(0);
            if f.to_lowercase().ends_with(".gz") { size * GZ_FACTOR } else { size }
        })
        .sum();

    if total_size == 0 {
        return (1000, 0);
    }

    let chunks = total_size.div_ceil(SIXTEEN_GB);
    ((chunks as usize) * 1000, total_size)
}
