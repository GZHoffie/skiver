use log::{info, warn, error};
use needletail::parse_fastx_file;
use rust_htslib::{bam, bam::Read as BamRead}; // Added rust-htslib
use serde::{Serialize, Deserialize};
//use rayon::prelude::*;

use std::fs::File;
use std::io::BufWriter;
use std::io::{prelude::*, BufReader};

use std::collections::HashMap;
use std::collections::HashSet;

use crate::{seeding::*, types::*, utils::*, constants::*};
use crate::summary::{ErrorSummary, ErrorSpectrumSummary, PhredScoreSummary, ReadPositionSummary};

/// kv-mer statistics for downstream analysis.
pub struct KVmerStats {
    pub k: u8,
    pub v: u8,
    pub keys: Vec<u64>,
    pub consensus_values: Vec<u64>,
    pub error_summary: ErrorSummary,
    pub error_spectrum: ErrorSpectrumSummary,
    pub phred_summary: PhredScoreSummary,
    pub read_position_summary: ReadPositionSummary,
}


#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct KVmerSet {
    pub key_size: u8,
    pub value_size: u8,
    pub kv_size: u8,
    pub num_kvmers: u32,

    /// key -> value -> list of per-observation metadata.
    /// The count of a (key, value) pair is `info_list.len()`.
    pub key_value_qual_map: HashMap<u64, HashMap<u64, Vec<ValueInfo>>>,

    // utilities to extract key and value from a kmer hash
    key_mask: u64,
    value_mask: u64,

    // whether both forward and reverse complement of the reads are included
    bidirectional: bool,
}


impl KVmerSet {
    pub fn new(key_size: u8, value_size: u8, bidirectional: bool) -> Self {
        assert!(key_size <= 32 && value_size <= 32, "Currently, we only support k, v <= 32.");

        let v_mask = (1 << (value_size * 2)) - 1;
        let k_mask = ((1 << (key_size * 2)) - 1) << (value_size * 2);

        KVmerSet {
            key_size,
            value_size,
            kv_size: key_size + value_size,
            num_kvmers: 0,
            key_value_qual_map: HashMap::new(),
            key_mask: k_mask,
            value_mask: v_mask,
            bidirectional,
        }
    }


    pub fn to_value_string(&self, kmer: u64) -> String {
        // for debugging: convert a kmer to a string

        let mut s = Vec::with_capacity(self.value_size as usize);
        for i in (0..self.value_size).rev() {
            let shift = i * 2;
            let base = ((kmer >> shift) & 0b11) as usize;
            s.push(crate::types::SEQ_TO_BYTE[base]);
        }
        String::from_utf8(s).unwrap()
    }

    pub fn to_key_string(&self, kmer: u64) -> String {
        // for debugging: convert a kmer to a string

        let mut s = Vec::with_capacity(self.key_size as usize);
        for i in (0..self.key_size).rev() {
            let shift = i * 2;
            let base = ((kmer >> shift) & 0b11) as usize;
            s.push(crate::types::SEQ_TO_BYTE[base]);
        }
        String::from_utf8(s).unwrap()
    }

    pub fn homopolymer_length(&self, key: u64, value: u64) -> u32 {
        let mut longest_homopolymer: u32 = 1;
        let mut current_homopolymer: u32 = 1;

        // Find the longest homopolymer at the end of the key
        let mut last_base = key & 0b11;
        for i in 1..self.key_size {
            let shift = i * 2;
            let base = (key >> shift) & 0b11;
            if base == last_base {
                current_homopolymer += 1;
            } else {
                break;
            }
        }
        // Extend the homopolymer into the value
        for i in (0..self.value_size).rev() {
            let shift = i * 2;
            let base = (value >> shift) & 0b11;
            if base == last_base {
                current_homopolymer += 1;
            } else {
                if current_homopolymer > longest_homopolymer {
                    longest_homopolymer = current_homopolymer;
                }
                current_homopolymer = 1;
                last_base = base;
            }
        }

        if current_homopolymer > longest_homopolymer {
            longest_homopolymer = current_homopolymer;
        }

        longest_homopolymer
    }


    /// Record a batch of (key, value, ValueInfo) triples.
    pub fn add_kv_qual_vector(&mut self, key_vec: &[u64], value_vec: &[u64], info_vec: &[ValueInfo]) {
        assert!(key_vec.len() == value_vec.len() && key_vec.len() == info_vec.len(),
                "Key, value, and info vectors must have the same length.");
        for ((&key, &value), info) in key_vec.iter().zip(value_vec.iter()).zip(info_vec.iter()) {
            self.key_value_qual_map
                .entry(key).or_insert_with(HashMap::new)
                .entry(value).or_insert_with(Vec::new)
                .push(info.clone());
        }
        self.num_kvmers += key_vec.len() as u32;
    }


    fn extract_markers_masked(&self, string: &[u8], key_vec: &mut Vec<u64>, value_vec: &mut Vec<u64>, c: usize, trim_front: usize, trim_back: usize, value_info_vec: &mut Vec<ValueInfo>) {
        let start = std::cmp::min(trim_front, string.len());
        let end = string.len().saturating_sub(trim_back);
        let string_trimmed = &string[start..end];
        // extract sketched kv-mers from the given sequence string
        #[cfg(any(target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                use crate::avx2_seeding::*;
                unsafe {
                    extract_markers_avx2_masked(string_trimmed, key_vec, value_vec, value_info_vec, c, self.key_size as usize, self.value_size as usize, self.bidirectional);
                }
            } else {
                fmh_seeds_masked(string_trimmed, key_vec, value_vec, value_info_vec, c, self.key_size as usize, self.value_size as usize, self.bidirectional);
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            fmh_seeds_masked(string_trimmed, key_vec, value_vec, value_info_vec, c, self.key_size as usize, self.value_size as usize, self.bidirectional);
        }
    }

    /// Like `extract_markers_masked`, but also extracts quality scores and builds `ValueInfo`.
    fn extract_markers_masked_with_qual(&self, string: &[u8], qual: &[u8], key_vec: &mut Vec<u64>, value_vec: &mut Vec<u64>, info_vec: &mut Vec<ValueInfo>, c: usize, trim_front: usize, trim_back: usize) {
        let start = std::cmp::min(trim_front, string.len());
        let end = string.len().saturating_sub(trim_back);
        let string_trimmed = &string[start..end];
        let qual_trimmed = &qual[start..end];
        #[cfg(any(target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                use crate::avx2_seeding::*;
                unsafe {
                    extract_markers_avx2_masked_with_qual(string_trimmed, qual_trimmed, key_vec, value_vec, info_vec, c, self.key_size as usize, self.value_size as usize, self.bidirectional);
                }
            } else {
                fmh_seeds_masked_with_qual(string_trimmed, qual_trimmed, key_vec, value_vec, info_vec, c, self.key_size as usize, self.value_size as usize, self.bidirectional);
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            fmh_seeds_masked_with_qual(string_trimmed, qual_trimmed, key_vec, value_vec, info_vec, c, self.key_size as usize, self.value_size as usize, self.bidirectional);
        }
    }

    // MODIFIED: Added BAM/SAM support
    pub fn add_file_to_kvmer_set(
        &mut self,
        seq_file: &str,
        c: usize,
        trim_front: usize,
        trim_back: usize,
    ) {
        let seq_file_clone = seq_file.to_string();

        if seq_file_clone.ends_with(".bam") || seq_file_clone.ends_with(".sam") {
            match bam::Reader::from_path(&seq_file_clone) {
                Ok(mut reader) => {
                    if !self.bidirectional {
                        // [FIXME] Correct the coverage estimation when using forward strand only with BAM/SAM input files
                        warn!("Using --forward-only with BAM/SAM input files may make the estimation of true coverage inaccurate.")
                    }
                    for record_result in reader.records() {
                        match record_result {
                            Ok(record) => {
                                let seq = record.seq().as_bytes();
                                let qual = record.qual().to_vec();
                                let mut key_vec: Vec<u64> = Vec::new();
                                let mut value_vec: Vec<u64> = Vec::new();
                                let mut info_vec: Vec<ValueInfo> = Vec::new();
                                self.extract_markers_masked_with_qual(&seq, &qual, &mut key_vec, &mut value_vec, &mut info_vec, c, trim_front, trim_back);
                                self.add_kv_qual_vector(&key_vec, &value_vec, &info_vec);
                            }
                            Err(e) => warn!("Error reading BAM/SAM record: {}", e),
                        }
                    }
                }
                Err(e) => error!("{} is not a valid BAM/SAM file (Error: {}); skipping.", seq_file_clone, e),
            }
        } else {
            let reader = parse_fastx_file(&seq_file_clone);
            if !reader.is_ok() {
                error!("{} is not a valid fasta/fastq file; skipping.", seq_file_clone);
                return;
            }
            let mut reader = reader.unwrap();
            while let Some(record) = reader.next() {
                match record {
                    Ok(record) => {
                        let mut key_vec: Vec<u64> = Vec::new();
                        let mut value_vec: Vec<u64> = Vec::new();
                        if let Some(qual) = record.qual() {
                            // FASTQ: record quality scores alongside k,v-mers.
                            let mut info_vec: Vec<ValueInfo> = Vec::new();
                            self.extract_markers_masked_with_qual(&record.seq(), qual, &mut key_vec, &mut value_vec, &mut info_vec, c, trim_front, trim_back);
                            self.add_kv_qual_vector(&key_vec, &value_vec, &info_vec);
                        } else {
                            // FASTA: no quality scores; record position/strand but empty qual.
                            let mut info_vec: Vec<ValueInfo> = Vec::new();
                            self.extract_markers_masked(&record.seq(), &mut key_vec, &mut value_vec, c, trim_front, trim_back, &mut info_vec);
                            self.add_kv_qual_vector(&key_vec, &value_vec, &info_vec);
                        }
                    }
                    Err(e) => warn!("Error reading record: {}", e),
                }
            }
        }
    }

    pub fn containment_index(&self, other: &KVmerSet) -> (f64, f64) {
        // check the key containment index and key-value pair containment index
        // each key/ key-value pair is counted once
        let mut shared_keys = 0;
        let mut shared_key_values = 0;
        let mut total_key_values = 0;

        for (key, value_map) in &self.key_value_qual_map {
            if let Some(other_value_map) = other.key_value_qual_map.get(key) {
                shared_keys += 1;

                for (value, _qual_list) in value_map {
                    if let Some(_other_qual_list) = other_value_map.get(value) {
                        shared_key_values += 1;
                    }
                }
            }
            total_key_values += value_map.len();
        }

        let key_containment = if self.key_value_qual_map.is_empty() {
            0.0
        } else {
            shared_keys as f64 / self.key_value_qual_map.len() as f64
        };

        let key_value_containment = if total_key_values == 0 {
            0.0
        } else {
            shared_key_values as f64 / total_key_values as f64
        };

        (key_containment, key_value_containment)
    }

    /**
     * Find the number of one-edit neighbors of the consensus value[0:v].
     * [FIXME] Optimize this function.
     */
    fn _num_consensus_up_to_v(&self, consensus: u64, v: u8, _bidirectional: bool, value_map: &HashMap<u64, Vec<ValueInfo>>) -> u32 {
        let consensus_up_to_v = consensus >> ((self.value_size - v) * 2);

        let mut num_consensus: u32 = 0;
        for (neighbor, info_list) in value_map {
            let _neighbors_up_to_v = neighbor >> ((self.value_size - v) * 2);
            if _neighbors_up_to_v == consensus_up_to_v {
                num_consensus += info_list.len() as u32;
            }
        }
        num_consensus
    }

    /// Walk every observation for a given key base-by-base against `consensus`.
    /// For each observation (qual string):
    ///   - empty qual string (FASTA source) → skipped entirely
    ///   - at each position p (0 = first/leftmost base = MSB pair):
    ///       match   → increment qscore_correct[phred] and continue
    ///       mismatch → increment qscore_error[phred] and stop this observation
    fn accumulate_qscore_calibration(
        consensus: u64,
        value_size: u8,
        value_map: &HashMap<u64, Vec<ValueInfo>>,
    ) -> (HashMap<u8, u64>, HashMap<u8, u64>) {
        let mut qscore_correct: HashMap<u8, u64> = HashMap::new();
        let mut qscore_error: HashMap<u8, u64> = HashMap::new();
        for (value, info_list) in value_map {
            for info in info_list {
                if info.qual.is_empty() {
                    continue; // no quality data (FASTA source)
                }
                for p in 0..value_size as usize {
                    let bit_shift = 2 * (value_size as usize - 1 - p);
                    let value_base     = (value     >> bit_shift) & 0b11;
                    let consensus_base = (consensus >> bit_shift) & 0b11;
                    let phred = info.qual[p].saturating_sub(33);
                    if value_base == consensus_base {
                        *qscore_correct.entry(phred).or_insert(0) += 1;
                    } else {
                        *qscore_error.entry(phred).or_insert(0) += 1;
                        break;
                    }
                }
            }
        }
        (qscore_correct, qscore_error)
    }

    /// Walk every observation for a given key base-by-base against `consensus`.
    /// For each position p in the value:
    ///   - correct → record pos_from_start and pos_from_end in correct maps
    ///   - mismatch → record in error maps and stop (same break logic as qscore calibration)
    /// Empty qual observations (FASTA) are skipped.
    fn accumulate_read_position_calibration(
        consensus: u64,
        value_size: u8,
        value_map: &HashMap<u64, Vec<ValueInfo>>,
    ) -> (HashMap<u32, u64>, HashMap<u32, u64>, HashMap<u32, u64>, HashMap<u32, u64>) {
        let mut correct_from_start: HashMap<u32, u64> = HashMap::new();
        let mut correct_from_end: HashMap<u32, u64> = HashMap::new();
        let mut error_from_start: HashMap<u32, u64> = HashMap::new();
        let mut error_from_end: HashMap<u32, u64> = HashMap::new();
        for (value, info_list) in value_map {
            for info in info_list {
                if info.qual.is_empty() {
                    continue;
                }
                for p in 0..value_size as usize {
                    let bit_shift = 2 * (value_size as usize - 1 - p);
                    let value_base     = (value     >> bit_shift) & 0b11;
                    let consensus_base = (consensus >> bit_shift) & 0b11;
                    let (pos_from_start, pos_from_end) = if info.is_forward {
                        (info.start_index + p as u32,
                         info.dist_to_read_end.saturating_sub(1 + p as u32))
                    } else {
                        (info.start_index.saturating_sub(p as u32),
                         info.dist_to_read_end + p as u32)
                    };
                    if value_base == consensus_base {
                        *correct_from_start.entry(pos_from_start).or_insert(0) += 1;
                        *correct_from_end.entry(pos_from_end).or_insert(0) += 1;
                    } else {
                        *error_from_start.entry(pos_from_start).or_insert(0) += 1;
                        *error_from_end.entry(pos_from_end).or_insert(0) += 1;
                        break;
                    }
                }
            }
        }
        (correct_from_start, correct_from_end, error_from_start, error_from_end)
    }

    pub fn get_stats(&self, threshold: u32) -> KVmerStats {
        let mut keys: Vec<u64> = Vec::new();
        let mut consensus_values: Vec<u64> = Vec::new();
        let mut error_summary = ErrorSummary::new(self.value_size as usize);
        let mut error_spectrum = ErrorSpectrumSummary::new();
        let mut phred_summary = PhredScoreSummary::new();
        let mut read_position_summary = ReadPositionSummary::new();

        for (key, value_map) in &self.key_value_qual_map {

            let mut max_count = 0;
            let mut sum_count = 0;
            let mut max_value: u64 = 0;

            // find the consensus value
            for (value, info_list) in value_map {
                let count = info_list.len() as u32;
                sum_count += count;
                if count > max_count {
                    max_count = count;
                    max_value = *value;
                }
            }

            // skip low coverage keys
            if sum_count <= threshold {
                continue;
            }

            let mut error_count_map: HashMap<(EditOperation, u8, u8), u32> = HashMap::new();
            let neighbors = _get_neighbors(max_value, self.value_size, self.bidirectional);
            if neighbors.contains_key(&max_value) {
                // This would confound the X=0 case
                continue;
            }

            let per_v_consensus: Vec<u32> = (1..=self.value_size)
                .map(|v| self._num_consensus_up_to_v(max_value, v, self.bidirectional, value_map))
                .collect();

            let mut num_neighbors = 0;
            for (value, info_list) in value_map {
                let count = info_list.len() as u32;
                if *value != max_value && neighbors.contains_key(value) {
                    let info = neighbors.get(value).unwrap();
                    let entry = error_count_map.entry((info.op, info.prev_base, info.next_base)).or_insert(0);
                    *entry += count;
                    num_neighbors += count;
                }
            }

            let (key_correct, key_error) = Self::accumulate_qscore_calibration(max_value, self.value_size, value_map);
            let (pos_correct_start, pos_correct_end, pos_error_start, pos_error_end) =
                Self::accumulate_read_position_calibration(max_value, self.value_size, value_map);

            keys.push(*key);
            consensus_values.push(max_value);
            error_summary.update(max_count, sum_count, num_neighbors, &per_v_consensus);
            error_spectrum.update(error_count_map);
            phred_summary.update(key_correct, key_error);
            read_position_summary.update(pos_correct_start, pos_correct_end, pos_error_start, pos_error_end);
        }

        KVmerStats {
            k: self.key_size,
            v: self.value_size,
            keys,
            consensus_values,
            error_summary,
            error_spectrum,
            phred_summary,
            read_position_summary,
        }
    }

    #[allow(unused)]
    pub fn get_stats_with_reference(&self, threshold: u32, reference: &KVmerSet) -> KVmerStats {
        let mut keys: Vec<u64> = Vec::new();
        let mut consensus_values: Vec<u64> = Vec::new();
        let mut error_summary = ErrorSummary::new(self.value_size as usize);
        let mut error_spectrum = ErrorSpectrumSummary::new();
        let mut phred_summary = PhredScoreSummary::new();
        let mut read_position_summary = ReadPositionSummary::new();

        // for debugging: the number of k-mers that the read set shares with the reference
        let mut shared_kmer_count: u32 = 0;

        for (key, ref_value_map) in &reference.key_value_qual_map {

            if !self.key_value_qual_map.contains_key(&key) {
                continue;
            }

            let consensus_value = *ref_value_map.keys().next().unwrap();
            let value_map = self.key_value_qual_map.get(&key).unwrap();

            let mut max_count = 0;
            let mut sum_count = 0;

            // find the consensus value
            for (value, info_list) in value_map {
                let count = info_list.len() as u32;
                sum_count += count;
                if count > max_count {
                    max_count = count;
                }
            }
            let consensus_count = value_map.get(&consensus_value).map_or(0, |q| q.len() as u32);
            shared_kmer_count += sum_count;

            if ref_value_map.len() > 1 {
                // skip non-unique reference kv-mers
                continue;
            }

            // [FIXME] skip if max_value != consensus_value?

            // skip low coverage keys
            if sum_count <= threshold {
                continue;
            }

            let mut error_count_map: HashMap<(EditOperation, u8, u8), u32> = HashMap::new();
            let neighbors = _get_neighbors(consensus_value, self.value_size, self.bidirectional);
            if neighbors.contains_key(&consensus_value) {
                // This would confound the X=0 case
                continue;
            }

            let per_v_consensus: Vec<u32> = (1..=self.value_size)
                .map(|v| self._num_consensus_up_to_v(consensus_value, v, self.bidirectional, value_map))
                .collect();

            let mut num_neighbors = 0;
            for (value, info_list) in value_map {
                let count = info_list.len() as u32;
                if *value != consensus_value && neighbors.contains_key(value) {
                    let info = neighbors.get(value).unwrap();
                    let entry = error_count_map.entry((info.op, info.prev_base, info.next_base)).or_insert(0);
                    *entry += count;
                    num_neighbors += count;
                }
            }

            let (key_correct, key_error) = Self::accumulate_qscore_calibration(consensus_value, self.value_size, value_map);
            let (pos_correct_start, pos_correct_end, pos_error_start, pos_error_end) =
                Self::accumulate_read_position_calibration(consensus_value, self.value_size, value_map);

            keys.push(*key);
            consensus_values.push(consensus_value);
            error_summary.update(consensus_count, sum_count, num_neighbors, &per_v_consensus);
            error_spectrum.update(error_count_map);
            phred_summary.update(key_correct, key_error);
            read_position_summary.update(pos_correct_start, pos_correct_end, pos_error_start, pos_error_end);
        }

        //println!("Total count of kvmers that match reference: {}", shared_kmer_count);
        //println!("Number of kvmers in read set: {}", self.num_kvmers);
        //println!("Proportion of kvmers that match reference: {:.4}%", shared_kmer_count as f64 / self.num_kvmers as f64 * 100.);

        KVmerStats {
            k: self.key_size,
            v: self.value_size,
            keys,
            consensus_values,
            error_summary,
            error_spectrum,
            phred_summary,
            read_position_summary,
        }
    }

    pub fn output_stats(&self, output_path: &String, stats: &KVmerStats, show_error_types: bool, show_error_vs_v: bool) {
        // create file for output
        let mut writer = File::create(&output_path).unwrap();
        // general info
        write!(writer, "key,consensus_value,homopolymer_length,consensus_count,neighbor_count,total_count").unwrap();
        // errors
        if show_error_types {
            for op in ALL_OPERATIONS {
                write!(writer, ",{:?}", op).unwrap();
            }
        }
        // for p vs. v regression
        if show_error_vs_v {
            for v in 1..=self.value_size {
                write!(writer, ",consensus_count_up_to_v{}", v).unwrap();
            }
        }

        writeln!(writer).unwrap();


        for i in 0..stats.keys.len() {
            write!(writer,
                "{},{},{},{},{},{}",
                self.to_key_string(stats.keys[i]),
                self.to_value_string(stats.consensus_values[i]),
                self.homopolymer_length(stats.keys[i], stats.consensus_values[i]),
                stats.error_summary.consensus_counts[i],
                stats.error_summary.neighbor_counts[i],
                stats.error_summary.total_counts[i],
            ).unwrap();
            if show_error_types {
                for op in ALL_OPERATIONS.iter() {
                    let mut total_count: u32 = 0;
                    for prev_base in 0..5 {
                        for next_base in 0..5 {
                            let count = stats.error_spectrum.error_counts[i].get(&(*op, prev_base, next_base)).unwrap_or(&0);
                            total_count += *count;
                        }
                    }
                    write!(writer, ",{}", total_count).unwrap();
                }
            }
            if show_error_vs_v {
                for v in 1..=self.value_size {
                    let consensus_count_up_to_v = stats.error_summary.consensus_up_to_v_counts[(v - 1) as usize][i];
                    write!(writer, ",{}", consensus_count_up_to_v).unwrap();
                }
            }
            writeln!(writer).unwrap();
        }
    }

    pub fn dump(&self, output_dir: &str) {

        //let mut file = &mut File::create_new(output_dir).unwrap();
        let mut writer = BufWriter::new(
            File::create(&output_dir)
                .expect(&format!("{} path not valid; exiting ", output_dir)),
        );
        //let config = bincode::config::standard().with_big_endian().with_fixed_int_encoding();

        bincode::serialize_into(&mut writer, &self).unwrap();
        info!("Sketching complete.");
    }

    pub fn load(&mut self, input_file: &str) {
        let file = File::open(input_file).expect(&format!("The sketch `{}` could not be opened. Exiting", input_file));
        let reader = BufReader::with_capacity(10_000_000, file);
        //let reader = BufReader::new(file);
        let that: KVmerSet = bincode::deserialize_from(reader)
            .expect(&format!(
                "The sketch `{}` is not a valid sketch.",
                &input_file
            ));

        // load the data into self
        if self.key_size != that.key_size || self.value_size != that.value_size {
            warn!("Key size or value size does not match when loading KVmerSet from file. Skipping input file {}.", input_file);
        } else {
            for (kmer, value_map) in that.key_value_qual_map {
                let entry = self.key_value_qual_map.entry(kmer).or_insert_with(HashMap::new);
                for (value, info_list) in value_map {
                    entry.entry(value).or_insert_with(Vec::new).extend(info_list);
                }
            }
            self.num_kvmers += that.num_kvmers;
        }
    }

}
