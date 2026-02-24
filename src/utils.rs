use crate::types::*;

use std::collections::HashMap;


/**
 * Get all neighbors (kmers with edit distance 1) of a given kmer value
 * Returns a hashmap of neighbor kmer value to (EditOperation, previous_base, next_base)
 */
pub fn _get_neighbors(value: u64, value_size: u8, bidirectional: bool) -> HashMap<u64, (EditOperation, u8, u8)> {
    // get all the values with edit distance 1 from the input value

    let mut neighbors: HashMap<u64, (EditOperation, u8, u8)> = HashMap::new();
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
                neighbors.insert(neighbor, (BASES_TO_SUBSTITUTION[current_base as usize][b as usize].unwrap(), previous_base, next_base));
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
                .and_modify(|(op, _prev, _next)| {
                    if *op != BASES_TO_INSERTION[b as usize].unwrap() {
                        *op = EditOperation::AMBIGUOUS
                    }
                })
                .or_insert(
                    (BASES_TO_INSERTION[b as usize].unwrap(), previous_base, next_base)
                );
            
            
            
            let right_part = value & ((1 << shift) - 1);
            let neighbor_delete = left_part | (right_part << 2) | b;
            let original_base = (value >> shift) & 0b11;
            neighbors.entry(neighbor_delete)
                .and_modify(|(op, _prev, _next)| 
                    if *op != BASES_TO_DELETION[original_base as usize].unwrap() {
                        *op = EditOperation::AMBIGUOUS
                    }
                )
                .or_insert(
                    (BASES_TO_DELETION[original_base as usize].unwrap(), previous_base, next_base)
                );
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
    for (neighbor, op) in neighbors {
        println!("Neighbor: {}, Operation: {}", _kmer_to_string(neighbor, k), sbs96_str(&op));
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
