use log::{info, warn, error};
use needletail::parse_fastx_file;
#[cfg(feature = "parallel-gzip")]
use needletail::parse_fastx_reader;
#[cfg(all(feature = "bam", not(target_os = "windows")))]
use rust_htslib::{bam, bam::Read as BamRead};
use serde::{Serialize, Deserialize};
//use rayon::prelude::*;

use std::fs::File;
use std::io::{BufWriter, BufReader};

use std::collections::HashMap;

use crate::{seeding::*, types::*};
use crate::summary::{KVmerSummary, ErrorSpectrumSummary, GCContentSummary, PhredScoreSummary, ReadPositionSummary};

/// kv-mer statistics for downstream analysis.
pub struct KVmerStats {
    pub k: u8,
    pub v: u8,
    pub keys: Vec<u64>,
    pub consensus_values: Vec<u64>,
    pub kvmer_summary: KVmerSummary,
    pub error_spectrum: ErrorSpectrumSummary,
    pub phred_summary: PhredScoreSummary,
    pub gc_content_summary: GCContentSummary,
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

    fn add_kv_qual_buffers(&mut self, key_vec: &mut Vec<u64>, value_vec: &mut Vec<u64>, info_vec: &mut Vec<ValueInfo>) {
        debug_assert!(key_vec.len() == value_vec.len() && key_vec.len() == info_vec.len());
        self.num_kvmers += key_vec.len() as u32;
        for ((key, value), info) in key_vec.drain(..).zip(value_vec.drain(..)).zip(info_vec.drain(..)) {
            self.key_value_qual_map
                .entry(key).or_insert_with(HashMap::new)
                .entry(value).or_insert_with(Vec::new)
                .push(info);
        }
    }

    fn merge(&mut self, other: KVmerSet) {
        if self.key_value_qual_map.is_empty() && self.num_kvmers == 0 {
            *self = other;
            return;
        }
        self.num_kvmers += other.num_kvmers;
        for (key, value_map) in other.key_value_qual_map {
            let target = self.key_value_qual_map.entry(key).or_insert_with(HashMap::new);
            for (value, info_list) in value_map {
                target.entry(value).or_insert_with(Vec::new).extend(info_list);
            }
        }
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
        let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
        self.add_file_to_kvmer_set_with_threads(seq_file, c, trim_front, trim_back, threads);
    }

    pub fn add_file_to_kvmer_set_with_threads(
        &mut self,
        seq_file: &str,
        c: usize,
        trim_front: usize,
        trim_back: usize,
        threads: usize,
    ) {
        let seq_file_clone = seq_file.to_string();

        #[cfg(all(feature = "bam", not(target_os = "windows")))]
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
            return;
        }

        #[cfg(any(not(feature = "bam"), target_os = "windows"))]
        if seq_file_clone.ends_with(".bam") || seq_file_clone.ends_with(".sam") {
            error!("BAM/SAM support is not enabled in this build; skipping {}.", seq_file_clone);
            return;
        }

        {
            #[cfg(feature = "parallel-gzip")]
            let lower_path = seq_file_clone.to_lowercase();
            #[cfg(feature = "parallel-gzip")]
            let mut worker_count = threads.max(1);
            #[cfg(not(feature = "parallel-gzip"))]
            let worker_count = threads.max(1);
            #[cfg(feature = "parallel-gzip")]
            let (reader, gzip_threads) = if lower_path.ends_with(".gz") && threads > 1 {
                let gzip_threads = (threads * 3 / 4).max(1).min(threads - 1);
                worker_count = (threads - gzip_threads).max(1);
                match rapidgzip_core::Decoder::builder().decoder_threads(gzip_threads).build() {
                    Ok(decoder) => match decoder.open(&seq_file_clone) {
                        Ok(reader) => (parse_fastx_reader(reader), Some(gzip_threads)),
                        Err(e) => {
                            warn!("Could not open {} with rapidgzip-core ({}); falling back to needletail's gzip decoder.", seq_file_clone, e);
                            (parse_fastx_file(&seq_file_clone), None)
                        }
                    },
                    Err(e) => {
                        warn!("Could not configure rapidgzip-core ({}); falling back to needletail's gzip decoder.", e);
                        (parse_fastx_file(&seq_file_clone), None)
                    }
                }
            } else {
                (parse_fastx_file(&seq_file_clone), None)
            };
            #[cfg(not(feature = "parallel-gzip"))]
            let reader = parse_fastx_file(&seq_file_clone);
            if !reader.is_ok() {
                error!("{} is not a valid fasta/fastq file; skipping.", seq_file_clone);
                return;
            }
            let mut reader = reader.unwrap();
            if worker_count == 1 {
                let mut key_vec = Vec::new();
                let mut value_vec = Vec::new();
                let mut info_vec = Vec::new();
                while let Some(record) = reader.next() {
                    match record {
                        Ok(record) => {
                            let seq = record.seq();
                            if let Some(qual) = record.qual() {
                                self.extract_markers_masked_with_qual(&seq, qual, &mut key_vec, &mut value_vec, &mut info_vec, c, trim_front, trim_back);
                            } else {
                                self.extract_markers_masked(&seq, &mut key_vec, &mut value_vec, c, trim_front, trim_back, &mut info_vec);
                            }
                            self.add_kv_qual_buffers(&mut key_vec, &mut value_vec, &mut info_vec);
                        }
                        Err(e) => warn!("Error reading record: {}", e),
                    }
                }
                return;
            }
            let (sender, receiver) = crossbeam_channel::bounded::<Vec<(Vec<u8>, Option<Vec<u8>>)>>(worker_count * 2);
            let key_size = self.key_size;
            let value_size = self.value_size;
            let bidirectional = self.bidirectional;
            let (worker_sets, parse_failed) = std::thread::scope(|scope| {
                let mut handles = Vec::with_capacity(worker_count);
                for _ in 0..worker_count {
                    let receiver = receiver.clone();
                    handles.push(scope.spawn(move || {
                        let mut local = KVmerSet::new(key_size, value_size, bidirectional);
                        let mut key_vec = Vec::new();
                        let mut value_vec = Vec::new();
                        let mut info_vec = Vec::new();
                        for batch in receiver {
                            for (seq, qual) in batch {
                                key_vec.clear();
                                value_vec.clear();
                                info_vec.clear();
                                if let Some(qual) = qual {
                                    local.extract_markers_masked_with_qual(&seq, &qual, &mut key_vec, &mut value_vec, &mut info_vec, c, trim_front, trim_back);
                                } else {
                                    local.extract_markers_masked(&seq, &mut key_vec, &mut value_vec, c, trim_front, trim_back, &mut info_vec);
                                }
                                local.add_kv_qual_buffers(&mut key_vec, &mut value_vec, &mut info_vec);
                            }
                        }
                        local
                    }));
                }
                drop(receiver);

                let mut batch = Vec::with_capacity(1024);
                let mut parse_failed = false;
                while let Some(record) = reader.next() {
                    match record {
                        Ok(record) => {
                            batch.push((record.seq().into_owned(), record.qual().map(|qual| qual.to_vec())));
                            if batch.len() == 1024 {
                                sender.send(std::mem::replace(&mut batch, Vec::with_capacity(1024))).unwrap();
                            }
                        }
                        Err(e) => {
                            warn!("Error reading record: {}", e);
                            parse_failed = true;
                        }
                    }
                }
                if !batch.is_empty() {
                    sender.send(batch).unwrap();
                }
                drop(sender);
                (handles.into_iter().map(|handle| handle.join().unwrap()).collect::<Vec<_>>(), parse_failed)
            });
            if parse_failed {
                error!("Parsing failed for {}; discarding partial results.", seq_file_clone);
                return;
            }
            #[cfg(feature = "parallel-gzip")]
            if let Some(gzip_threads) = gzip_threads {
                info!("Parallel gzip decoding completed using up to {} thread(s).", gzip_threads);
            }
            for worker_set in worker_sets {
                self.merge(worker_set);
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

    pub fn get_stats(&self, threshold: u32) -> KVmerStats {
        let mut keys: Vec<u64> = Vec::new();
        let mut consensus_values: Vec<u64> = Vec::new();
        let mut kvmer_summary = KVmerSummary::new(self.value_size as usize);
        let mut error_spectrum = ErrorSpectrumSummary::new(self.value_size as usize);
        let mut phred_summary = PhredScoreSummary::new();
        let mut gc_content_summary = GCContentSummary::new();
        let mut read_position_summary = ReadPositionSummary::new();

        let mut sorted_keys: Vec<u64> = self.key_value_qual_map.keys().copied().collect();
        sorted_keys.sort_unstable();
        for key in sorted_keys {
            let value_map = &self.key_value_qual_map[&key];
            // find the consensus (most frequent) value
            let mut max_count = 0;
            let mut sum_count = 0;
            let mut max_value: u64 = 0;
            for (value, info_list) in value_map {
                let count = info_list.len() as u32;
                sum_count += count;
                if count > max_count || (count == max_count && *value < max_value) {
                    max_count = count;
                    max_value = *value;
                }
            }

            // skip low coverage keys
            if sum_count <= threshold {
                continue;
            }

            if kvmer_summary.update(key, max_value, self.key_size, self.value_size, self.bidirectional, value_map) {
                keys.push(key);
                consensus_values.push(max_value);
                error_spectrum.update(kvmer_summary.error_counts_per_key.last().unwrap().clone(), kvmer_summary.forward_error_counts_per_key.last().unwrap().clone());
                phred_summary.update(max_value, self.value_size, value_map);
                gc_content_summary.update(max_value, self.value_size, value_map);
                read_position_summary.update(max_value, self.value_size, value_map);
            }
        }

        KVmerStats {
            k: self.key_size,
            v: self.value_size,
            keys,
            consensus_values,
            kvmer_summary,
            error_spectrum,
            phred_summary,
            gc_content_summary,
            read_position_summary,
        }
    }

    #[allow(unused)]
    pub fn get_stats_with_reference(&self, threshold: u32, reference: &KVmerSet) -> KVmerStats {
        let mut keys: Vec<u64> = Vec::new();
        let mut consensus_values: Vec<u64> = Vec::new();
        let mut kvmer_summary = KVmerSummary::new(self.value_size as usize);
        let mut error_spectrum = ErrorSpectrumSummary::new(self.value_size as usize);
        let mut phred_summary = PhredScoreSummary::new();
        let mut gc_content_summary = GCContentSummary::new();
        let mut read_position_summary = ReadPositionSummary::new();

        // for debugging: the number of k-mers that the read set shares with the reference
        let mut shared_kmer_count: u32 = 0;

        let mut sorted_keys: Vec<u64> = reference.key_value_qual_map.keys().copied().collect();
        sorted_keys.sort_unstable();
        for key in sorted_keys {
            let ref_value_map = &reference.key_value_qual_map[&key];

            if !self.key_value_qual_map.contains_key(&key) {
                continue;
            }

            let consensus_value = *ref_value_map.keys().next().unwrap();
            let value_map = self.key_value_qual_map.get(&key).unwrap();

            let sum_count: u32 = value_map.values().map(|v| v.len() as u32).sum();
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

            if kvmer_summary.update(key, consensus_value, self.key_size, self.value_size, self.bidirectional, value_map) {
                keys.push(key);
                consensus_values.push(consensus_value);
                error_spectrum.update(kvmer_summary.error_counts_per_key.last().unwrap().clone(), kvmer_summary.forward_error_counts_per_key.last().unwrap().clone());
                phred_summary.update(consensus_value, self.value_size, value_map);
                gc_content_summary.update(consensus_value, self.value_size, value_map);
                read_position_summary.update(consensus_value, self.value_size, value_map);
            }
        }

        //println!("Total count of kvmers that match reference: {}", shared_kmer_count);
        //println!("Number of kvmers in read set: {}", self.num_kvmers);
        //println!("Proportion of kvmers that match reference: {:.4}%", shared_kmer_count as f64 / self.num_kvmers as f64 * 100.);

        KVmerStats {
            k: self.key_size,
            v: self.value_size,
            keys,
            consensus_values,
            kvmer_summary,
            error_spectrum,
            phred_summary,
            gc_content_summary,
            read_position_summary,
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

    pub fn load(&mut self, input_file: &str) -> Result<(), String> {
        let file = File::open(input_file)
            .map_err(|e| format!("Unable to open file `{}`: {}", input_file, e))?;
        let reader = BufReader::with_capacity(10_000_000, file);
        //let reader = BufReader::new(file);
        let that: KVmerSet = bincode::deserialize_from(reader)
            .map_err(|_| format!(
                "The input `{}` is not recognized as either FASTA/FASTQ or a valid kv-mer sketch. The sketch may have been generated by an older version of skiver",
                &input_file
            ))?;

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
        Ok(())
    }

}
