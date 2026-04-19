use std::collections::HashMap;
use crate::types::{EditOperation, ALL_OPERATIONS, SEQ_TO_BYTE, SEQ_TO_CHAR, ValueInfo, NeighborInfo};
use crate::utils::_get_neighbors;
use log::info;

fn logit(p: f64) -> f64 { (p / (1.0 - p)).ln() }
fn sigmoid(x: f64) -> f64 { 1.0 / (1.0 + (-x).exp()) }

/// Find δ such that the weighted average of `sigmoid(logit(raw_rate_i) + δ)` equals `target`.
/// `bins` is a slice of `(raw_rate, weight)` pairs; bins with `raw_rate == 0` contribute 0
/// regardless of δ and are excluded from the search (they don't affect the result).
/// Returns δ = 0.0 if all bins have zero rate or total weight is zero.
fn find_logit_shift(bins: &[(f64, f64)], target: f64) -> f64 {
    let total_weight: f64 = bins.iter().map(|(_, w)| w).sum();
    if total_weight == 0.0 { return 0.0; }

    let weighted_avg = |delta: f64| -> f64 {
        bins.iter().map(|(r, w)| {
            let norm = if *r <= 0.0 { 0.0 } else if *r >= 1.0 { 1.0 } else {
                sigmoid(logit(*r) + delta)
            };
            norm * w
        }).sum::<f64>() / total_weight
    };

    let mut lo = -50.0_f64;
    let mut hi =  50.0_f64;
    for _ in 0..64 {
        let mid = (lo + hi) * 0.5;
        if weighted_avg(mid) < target { lo = mid; } else { hi = mid; }
    }
    (lo + hi) * 0.5
}

/// Per-key error-rate statistics.
/// Corresponds to `KVmerStats` fields: `consensus_counts`, `total_counts`,
/// `neighbor_counts`, `consensus_up_to_v_counts`.
pub struct ErrorSummary {
    pub consensus_counts: Vec<u32>,
    pub total_counts: Vec<u32>,
    pub neighbor_counts: Vec<u32>,
    /// Per-key consensus counts for each value prefix length, indexed `[v-1][key_idx]`.
    pub consensus_up_to_v_counts: Vec<Vec<u32>>,
    pub key_strings: Vec<String>,
    pub value_strings: Vec<String>,
    pub second_value_strings: Vec<String>,
    pub second_counts: Vec<u32>,
    pub homopolymer_lengths: Vec<u32>,
    pub error_counts_per_key: Vec<HashMap<NeighborInfo, u32>>,
    pub forward_error_counts_per_key: Vec<HashMap<NeighborInfo, u32>>,
    v: usize,
}

impl ErrorSummary {
    pub fn new(v: usize) -> Self {
        ErrorSummary {
            consensus_counts: Vec::new(),
            total_counts: Vec::new(),
            neighbor_counts: Vec::new(),
            consensus_up_to_v_counts: vec![Vec::new(); v],
            key_strings: Vec::new(),
            value_strings: Vec::new(),
            second_value_strings: Vec::new(),
            second_counts: Vec::new(),
            homopolymer_lengths: Vec::new(),
            error_counts_per_key: Vec::new(),
            forward_error_counts_per_key: Vec::new(),
            v,
        }
    }

    fn to_kmer_string(kmer: u64, size: u8) -> String {
        let mut s = Vec::with_capacity(size as usize);
        for i in (0..size).rev() {
            s.push(SEQ_TO_BYTE[((kmer >> (i * 2)) & 0b11) as usize]);
        }
        String::from_utf8(s).unwrap()
    }

    fn homopolymer_length(key: u64, key_size: u8, value: u64, value_size: u8) -> u32 {
        let mut longest: u32 = 1;
        let mut current: u32 = 1;
        let mut last_base = key & 0b11;
        for i in 1..key_size {
            let base = (key >> (i * 2)) & 0b11;
            if base == last_base { current += 1; } else { break; }
        }
        for i in (0..value_size).rev() {
            let base = (value >> (i * 2)) & 0b11;
            if base == last_base {
                current += 1;
            } else {
                if current > longest { longest = current; }
                current = 1;
                last_base = base;
            }
        }
        if current > longest { longest = current; }
        longest
    }

    fn num_consensus_up_to_v(consensus: u64, v: u8, value_size: u8, value_map: &HashMap<u64, Vec<ValueInfo>>) -> u32 {
        let prefix = consensus >> ((value_size - v) * 2);
        value_map.iter().map(|(neighbor, info_list)| {
            if (neighbor >> ((value_size - v) * 2)) == prefix { info_list.len() as u32 } else { 0 }
        }).sum()
    }

    /// Accumulate one key's error statistics, computing all derived values from the raw inputs.
    /// Returns `false` (and skips insertion) if `consensus` is its own one-edit neighbor,
    /// which would confound the X=0 error case.
    pub fn update(
        &mut self,
        key: u64,
        consensus: u64,
        key_size: u8,
        value_size: u8,
        bidirectional: bool,
        value_map: &HashMap<u64, Vec<ValueInfo>>,
    ) -> bool {
        let consensus_count = value_map.get(&consensus).map_or(0, |q| q.len() as u32);
        let sum_count: u32 = value_map.values().map(|v| v.len() as u32).sum();
        let key_string = Self::to_kmer_string(key, key_size);
        let value_string = Self::to_kmer_string(consensus, value_size);
        let homopolymer_length = Self::homopolymer_length(key, key_size, consensus, value_size);

        // neighbors filter: skip keys whose consensus is its own neighbor
        let neighbors = _get_neighbors(consensus, value_size, bidirectional);
        if neighbors.contains_key(&consensus) {
            return false;
        }

        let per_v_consensus: Vec<u32> = (1..=value_size)
            .map(|v| Self::num_consensus_up_to_v(consensus, v, value_size, value_map))
            .collect();

        let mut error_count_map: HashMap<NeighborInfo, u32> = HashMap::new();
        let mut forward_error_count_map: HashMap<NeighborInfo, u32> = HashMap::new();
        let mut num_neighbors: u32 = 0;
        for (value, info_list) in value_map {
            let count = info_list.len() as u32;
            if *value != consensus {
                if let Some(info) = neighbors.get(value) {
                    *error_count_map.entry(*info).or_insert(0) += count;
                    num_neighbors += count;
                    let forward_count = info_list.iter().filter(|i| i.is_forward).count() as u32;
                    *forward_error_count_map.entry(*info).or_insert(0) += forward_count;
                }
            }
        }

        // find second most common value (highest-count non-consensus value)
        let second = value_map.iter()
            .filter(|&(&v, _)| v != consensus)
            .max_by_key(|(_, info_list)| info_list.len());
        let (second_value_string, second_count) = match second {
            Some((&v, info_list)) => (Self::to_kmer_string(v, value_size), info_list.len() as u32),
            None => (String::new(), 0),
        };

        // store
        self.consensus_counts.push(consensus_count);
        self.total_counts.push(sum_count);
        self.neighbor_counts.push(num_neighbors);
        for (j, &c) in per_v_consensus.iter().enumerate() {
            if j < self.consensus_up_to_v_counts.len() {
                self.consensus_up_to_v_counts[j].push(c);
            }
        }
        self.key_strings.push(key_string);
        self.value_strings.push(value_string);
        self.second_value_strings.push(second_value_string);
        self.second_counts.push(second_count);
        self.homopolymer_lengths.push(homopolymer_length);
        self.error_counts_per_key.push(error_count_map);
        self.forward_error_counts_per_key.push(forward_error_count_map);

        true
    }
}

impl ErrorSummary {
    pub fn to_csv(&self, indices: Option<&[usize]>) -> String {
        use std::fmt::Write;
        use std::collections::HashSet;
        let n = self.consensus_counts.len();
        let index_set: HashSet<usize> = match indices {
            Some(idx) => idx.iter().copied().collect(),
            None => (0..n).collect(),
        };
        let mut out = String::new();
        write!(out, "key,consensus_value,passes_filter,homopolymer_length,consensus_count,neighbor_count,total_count").unwrap();
        for op in ALL_OPERATIONS {
            write!(out, ",{:?}", op).unwrap();
        }
        for v in 1..=self.v {
            write!(out, ",consensus_count_up_to_v{}", v).unwrap();
        }
        writeln!(out).unwrap();
        for i in 0..n {
            write!(out,
                "{},{},{},{},{},{},{}",
                self.key_strings[i],
                self.value_strings[i],
                index_set.contains(&i),
                self.homopolymer_lengths[i],
                self.consensus_counts[i],
                self.neighbor_counts[i],
                self.total_counts[i],
            ).unwrap();
            for op in ALL_OPERATIONS.iter() {
                let total_count: u32 = self.error_counts_per_key[i].iter()
                    .filter(|(ni, _)| ni.op == *op)
                    .map(|(_, &c)| c)
                    .sum();
                write!(out, ",{}", total_count).unwrap();
            }
            for v in 1..=self.v {
                let consensus_count_up_to_v = self.consensus_up_to_v_counts[v - 1][i];
                write!(out, ",{}", consensus_count_up_to_v).unwrap();
            }
            writeln!(out).unwrap();
        }
        out
    }
}

/// Per-key error-type spectrum statistics.
/// Corresponds to `KVmerStats` field: `error_counts`.
pub struct ErrorSpectrumSummary {
    pub error_counts: Vec<HashMap<NeighborInfo, u32>>,
    pub forward_error_counts: Vec<HashMap<NeighborInfo, u32>>,
    v: usize,
}

impl ErrorSpectrumSummary {
    pub fn new(v: usize) -> Self {
        ErrorSpectrumSummary {
            error_counts: Vec::new(),
            forward_error_counts: Vec::new(),
            v,
        }
    }

    /// Accumulate one key's per-operation error counts.
    pub fn update(&mut self, error_map: HashMap<NeighborInfo, u32>, forward_error_map: HashMap<NeighborInfo, u32>) {
        self.error_counts.push(error_map);
        self.forward_error_counts.push(forward_error_map);
    }
}

impl ErrorSpectrumSummary {
    pub fn to_dependence_on_t_csv(&self, indices: Option<&[usize]>, k: usize, ignore_smallest_t: usize, ignore_largest_t: usize) -> String {
        use std::fmt::Write;
        let all: Vec<usize>;
        let indices = match indices {
            Some(idx) => idx,
            None => { all = (0..self.error_counts.len()).collect(); &all }
        };
        // Aggregate counts by (op, prev_base, next_base, position) for the given indices.
        let mut totals: HashMap<(EditOperation, u8, u8, u8), u64> = HashMap::new();
        for &i in indices {
            for (ni, &count) in &self.error_counts[i] {
                *totals.entry((ni.op, ni.prev_base, ni.next_base, ni.position)).or_insert(0) += count as u64;
            }
        }

        let v_min = 1 + ignore_smallest_t;
        let v_max = self.v.saturating_sub(ignore_largest_t);
        let mut out = String::new();
        write!(out, "operation,prev_base,next_base,total").unwrap();
        for pos in v_min..=v_max {
            write!(out, ",freq_at_t{}", k + pos).unwrap();
        }
        writeln!(out).unwrap();

        for &op in ALL_OPERATIONS.iter() {
            for prev_base in 0u8..4 {
                for next_base in 0u8..4 {
                    let counts: Vec<u64> = (v_min..=v_max)
                        .map(|pos| totals.get(&(op, prev_base, next_base, pos as u8)).copied().unwrap_or(0))
                        .collect();
                    let total: u64 = counts.iter().sum();
                    if total > 0 {
                        write!(out, "{},{},{},{}",
                            op,
                            SEQ_TO_CHAR[prev_base as usize],
                            SEQ_TO_CHAR[next_base as usize],
                            total,
                        ).unwrap();
                        for c in &counts {
                            write!(out, ",{}", c).unwrap();
                        }
                        writeln!(out).unwrap();
                    }
                }
            }
        }
        out
    }

    pub fn to_csv(&self, indices: Option<&[usize]>) -> String {
        use std::fmt::Write;
        let all: Vec<usize>;
        let indices = match indices {
            Some(idx) => idx,
            None => { all = (0..self.error_counts.len()).collect(); &all }
        };
        // Aggregate total and forward-strand counts by (op, prev_base, next_base).
        let mut totals: HashMap<(EditOperation, u8, u8), u64> = HashMap::new();
        let mut forward_totals: HashMap<(EditOperation, u8, u8), u64> = HashMap::new();
        for &i in indices {
            for (ni, &count) in &self.error_counts[i] {
                *totals.entry((ni.op, ni.prev_base, ni.next_base)).or_insert(0) += count as u64;
            }
            for (ni, &count) in &self.forward_error_counts[i] {
                *forward_totals.entry((ni.op, ni.prev_base, ni.next_base)).or_insert(0) += count as u64;
            }
        }

        let mut out = String::new();
        writeln!(out, "operation,prev_base,next_base,total,forward").unwrap();

        for &op in ALL_OPERATIONS.iter() {
            for prev_base in 0u8..4 {
                for next_base in 0u8..4 {
                    let key = (op, prev_base, next_base);
                    let total = totals.get(&key).copied().unwrap_or(0);
                    if total > 0 {
                        let forward = forward_totals.get(&key).copied().unwrap_or(0);
                        writeln!(out, "{},{},{},{},{}",
                            op,
                            SEQ_TO_CHAR[prev_base as usize],
                            SEQ_TO_CHAR[next_base as usize],
                            total,
                            forward,
                        ).unwrap();
                    }
                }
            }
        }
        out
    }
}

/// Phred quality-score calibration statistics broken down by error type.
/// `observed` counts every quality-score occurrence across all values (consensus + neighbors).
/// `substitution`, `insertion`, `deletion` count only the neighbor observations by type.
pub struct PhredScoreSummary {
    pub observed: HashMap<u8, u64>,
    pub substitution: HashMap<u8, u64>,
    pub insertion: HashMap<u8, u64>,
    pub deletion: HashMap<u8, u64>,
    pub observed_per_key: Vec<HashMap<u8, u64>>,
    pub substitution_per_key: Vec<HashMap<u8, u64>>,
    pub insertion_per_key: Vec<HashMap<u8, u64>>,
    pub deletion_per_key: Vec<HashMap<u8, u64>>,
}

impl PhredScoreSummary {
    pub fn new() -> Self {
        PhredScoreSummary {
            observed: HashMap::new(),
            substitution: HashMap::new(),
            insertion: HashMap::new(),
            deletion: HashMap::new(),
            observed_per_key: Vec::new(),
            substitution_per_key: Vec::new(),
            insertion_per_key: Vec::new(),
            deletion_per_key: Vec::new(),
        }
    }

    /// Accumulate one key's Phred calibration data using edit-distance-1 neighbors.
    /// Every value (consensus and neighbors alike) contributes the quality score of every base
    /// position to `observed`. Neighbor values additionally increment the appropriate error-type
    /// bucket at the quality score of the specific error position (ni.position, right-indexed).
    pub fn update(
        &mut self,
        consensus: u64,
        value_size: u8,
        bidirectional: bool,
        value_map: &HashMap<u64, Vec<ValueInfo>>,
    ) {
        let neighbors = _get_neighbors(consensus, value_size, bidirectional);

        let mut key_observed: HashMap<u8, u64> = HashMap::new();
        let mut key_substitution: HashMap<u8, u64> = HashMap::new();
        let mut key_insertion: HashMap<u8, u64> = HashMap::new();
        let mut key_deletion: HashMap<u8, u64> = HashMap::new();

        for (value, info_list) in value_map {
            if *value == consensus {
                for info in info_list {
                    if info.qual.is_empty() { continue; }
                    for p in 0..value_size as usize {
                        let phred = info.qual[p].saturating_sub(33);
                        *key_observed.entry(phred).or_insert(0) += 1;
                    }
                }
            } else if let Some(ni) = neighbors.get(value) {
                // ni.position is right-indexed (0 = LSB); convert to left-indexed qual offset
                let qual_idx = (value_size as usize).saturating_sub(1 + ni.position as usize);
                for info in info_list {
                    if info.qual.is_empty() { continue; }
                    // All base positions contribute to observed.
                    for p in 0..value_size as usize {
                        let phred = info.qual[p].saturating_sub(33);
                        *key_observed.entry(phred).or_insert(0) += 1;
                    }
                    // Error type is recorded at the specific error position.
                    if qual_idx < info.qual.len() {
                        let phred = info.qual[qual_idx].saturating_sub(33);
                        match ni.op {
                            EditOperation::AC | EditOperation::AG | EditOperation::AT |
                            EditOperation::CA | EditOperation::CG | EditOperation::CT |
                            EditOperation::GA | EditOperation::GC | EditOperation::GT |
                            EditOperation::TA | EditOperation::TC | EditOperation::TG => {
                                *key_substitution.entry(phred).or_insert(0) += 1;
                            }
                            EditOperation::_A | EditOperation::_C |
                            EditOperation::_G | EditOperation::_T => {
                                *key_insertion.entry(phred).or_insert(0) += 1;
                            }
                            EditOperation::A_ | EditOperation::C_ |
                            EditOperation::G_ | EditOperation::T_ => {
                                *key_deletion.entry(phred).or_insert(0) += 1;
                            }
                            EditOperation::AMBIGUOUS => {}
                        }
                    }
                }
            }
        }

        for (&q, &c) in &key_observed     { *self.observed.entry(q).or_insert(0)      += c; }
        for (&q, &c) in &key_substitution { *self.substitution.entry(q).or_insert(0)  += c; }
        for (&q, &c) in &key_insertion    { *self.insertion.entry(q).or_insert(0)      += c; }
        for (&q, &c) in &key_deletion     { *self.deletion.entry(q).or_insert(0)       += c; }
        self.observed_per_key.push(key_observed);
        self.substitution_per_key.push(key_substitution);
        self.insertion_per_key.push(key_insertion);
        self.deletion_per_key.push(key_deletion);
    }
}

/// Read-position error calibration statistics broken down by error type.
/// `observed` counts every read-position occurrence across all values (consensus + neighbors).
/// `substitution`, `insertion`, `deletion` count only the neighbor observations by type.
pub struct ReadPositionSummary {
    pub observed_from_start_per_key: Vec<HashMap<u32, u64>>,
    pub observed_from_end_per_key: Vec<HashMap<u32, u64>>,
    pub substitution_from_start_per_key: Vec<HashMap<u32, u64>>,
    pub substitution_from_end_per_key: Vec<HashMap<u32, u64>>,
    pub insertion_from_start_per_key: Vec<HashMap<u32, u64>>,
    pub insertion_from_end_per_key: Vec<HashMap<u32, u64>>,
    pub deletion_from_start_per_key: Vec<HashMap<u32, u64>>,
    pub deletion_from_end_per_key: Vec<HashMap<u32, u64>>,
}

impl ReadPositionSummary {
    pub fn new() -> Self {
        ReadPositionSummary {
            observed_from_start_per_key: Vec::new(),
            observed_from_end_per_key: Vec::new(),
            substitution_from_start_per_key: Vec::new(),
            substitution_from_end_per_key: Vec::new(),
            insertion_from_start_per_key: Vec::new(),
            insertion_from_end_per_key: Vec::new(),
            deletion_from_start_per_key: Vec::new(),
            deletion_from_end_per_key: Vec::new(),
        }
    }

    /// Accumulate one key's read-position calibration data using edit-distance-1 neighbors.
    /// Every value (consensus and neighbors alike) contributes the read position of every base
    /// to `observed`. Neighbor values additionally increment the appropriate error-type bucket
    /// at the read position of the specific error position (ni.position, right-indexed).
    /// Deletions are attributed to both the previous-base and next-base flanking positions.
    pub fn update(
        &mut self,
        consensus: u64,
        value_size: u8,
        bidirectional: bool,
        value_map: &HashMap<u64, Vec<ValueInfo>>,
    ) {
        let neighbors = _get_neighbors(consensus, value_size, bidirectional);

        let mut observed_from_start: HashMap<u32, u64> = HashMap::new();
        let mut observed_from_end: HashMap<u32, u64> = HashMap::new();
        let mut substitution_from_start: HashMap<u32, u64> = HashMap::new();
        let mut substitution_from_end: HashMap<u32, u64> = HashMap::new();
        let mut insertion_from_start: HashMap<u32, u64> = HashMap::new();
        let mut insertion_from_end: HashMap<u32, u64> = HashMap::new();
        let mut deletion_from_start: HashMap<u32, u64> = HashMap::new();
        let mut deletion_from_end: HashMap<u32, u64> = HashMap::new();

        // Helper: convert left-indexed position p and ValueInfo to read positions.
        let read_pos = |info: &ValueInfo, p: usize| -> (u32, u32) {
            if info.is_forward {
                (info.start_index + p as u32,
                 info.dist_to_read_end.saturating_sub(1 + p as u32))
            } else {
                (info.start_index.saturating_sub(p as u32),
                 info.dist_to_read_end + p as u32)
            }
        };

        for (value, info_list) in value_map {
            if *value == consensus {
                for info in info_list {
                    if info.qual.is_empty() { continue; }
                    for p in 0..value_size as usize {
                        let (s, e) = read_pos(info, p);
                        *observed_from_start.entry(s).or_insert(0) += 1;
                        *observed_from_end.entry(e).or_insert(0) += 1;
                    }
                }
            } else if let Some(ni) = neighbors.get(value) {
                // ni.position is right-indexed; convert to left-indexed position p
                let p = (value_size as usize).saturating_sub(1 + ni.position as usize);
                for info in info_list {
                    if info.qual.is_empty() { continue; }
                    // All base positions contribute to observed.
                    for q in 0..value_size as usize {
                        let (s, e) = read_pos(info, q);
                        *observed_from_start.entry(s).or_insert(0) += 1;
                        *observed_from_end.entry(e).or_insert(0) += 1;
                    }
                    // Error type is recorded at the specific error position.
                    let (ps, pe) = read_pos(info, p);
                    match ni.op {
                        EditOperation::AC | EditOperation::AG | EditOperation::AT |
                        EditOperation::CA | EditOperation::CG | EditOperation::CT |
                        EditOperation::GA | EditOperation::GC | EditOperation::GT |
                        EditOperation::TA | EditOperation::TC | EditOperation::TG => {
                            *substitution_from_start.entry(ps).or_insert(0) += 1;
                            *substitution_from_end.entry(pe).or_insert(0) += 1;
                        }
                        EditOperation::_A | EditOperation::_C |
                        EditOperation::_G | EditOperation::_T => {
                            *insertion_from_start.entry(ps).or_insert(0) += 1;
                            *insertion_from_end.entry(pe).or_insert(0) += 1;
                        }
                        EditOperation::A_ | EditOperation::C_ |
                        EditOperation::G_ | EditOperation::T_ => {
                            // Attribute deletion to both flanking positions.
                            let p_prev = (value_size as usize).saturating_sub(1 + ni.prev_base as usize);
                            let p_next = (value_size as usize).saturating_sub(1 + ni.next_base as usize);
                            for &dp in &[p_prev, p_next] {
                                let (dps, dpe) = read_pos(info, dp);
                                *deletion_from_start.entry(dps).or_insert(0) += 1;
                                *deletion_from_end.entry(dpe).or_insert(0) += 1;
                            }
                        }
                        EditOperation::AMBIGUOUS => {}
                    }
                }
            }
        }

        self.observed_from_start_per_key.push(observed_from_start);
        self.observed_from_end_per_key.push(observed_from_end);
        self.substitution_from_start_per_key.push(substitution_from_start);
        self.substitution_from_end_per_key.push(substitution_from_end);
        self.insertion_from_start_per_key.push(insertion_from_start);
        self.insertion_from_end_per_key.push(insertion_from_end);
        self.deletion_from_start_per_key.push(deletion_from_start);
        self.deletion_from_end_per_key.push(deletion_from_end);
    }
}

impl ReadPositionSummary {
    /// `per_base_error_rate` is used for the same normalisation as `PhredScoreSummary::to_csv`:
    /// start and end positions are normalised independently via a logit-shift so that each
    /// direction's weighted average of `normalized_error_rate` equals `per_base_error_rate`.
    /// `normalized_error_rate` is always in `[0, 1]`.
    pub fn to_csv(&self, indices: Option<&[usize]>, per_base_error_rate: f64) -> String {
        use std::fmt::Write;
        let all: Vec<usize>;
        let indices = match indices {
            Some(idx) => idx,
            None => { all = (0..self.observed_from_start_per_key.len()).collect(); &all }
        };
        let mut observed_from_start: HashMap<u32, u64> = HashMap::new();
        let mut observed_from_end: HashMap<u32, u64> = HashMap::new();
        let mut substitution_from_start: HashMap<u32, u64> = HashMap::new();
        let mut substitution_from_end: HashMap<u32, u64> = HashMap::new();
        let mut insertion_from_start: HashMap<u32, u64> = HashMap::new();
        let mut insertion_from_end: HashMap<u32, u64> = HashMap::new();
        let mut deletion_from_start: HashMap<u32, u64> = HashMap::new();
        let mut deletion_from_end: HashMap<u32, u64> = HashMap::new();

        for &i in indices {
            for (&pos, &c) in &self.observed_from_start_per_key[i]     { *observed_from_start.entry(pos).or_insert(0)     += c; }
            for (&pos, &c) in &self.observed_from_end_per_key[i]       { *observed_from_end.entry(pos).or_insert(0)       += c; }
            for (&pos, &c) in &self.substitution_from_start_per_key[i] { *substitution_from_start.entry(pos).or_insert(0) += c; }
            for (&pos, &c) in &self.substitution_from_end_per_key[i]   { *substitution_from_end.entry(pos).or_insert(0)   += c; }
            for (&pos, &c) in &self.insertion_from_start_per_key[i]    { *insertion_from_start.entry(pos).or_insert(0)    += c; }
            for (&pos, &c) in &self.insertion_from_end_per_key[i]      { *insertion_from_end.entry(pos).or_insert(0)      += c; }
            for (&pos, &c) in &self.deletion_from_start_per_key[i]     { *deletion_from_start.entry(pos).or_insert(0)     += c; }
            for (&pos, &c) in &self.deletion_from_end_per_key[i]       { *deletion_from_end.entry(pos).or_insert(0)       += c; }
        }

        // Compute logit-shift δ for each direction independently.
        let bins_for = |obs_map: &HashMap<u32, u64>,
                        sub_map: &HashMap<u32, u64>,
                        ins_map: &HashMap<u32, u64>,
                        del_map: &HashMap<u32, u64>| -> f64 {
            let bins: Vec<(f64, f64)> = obs_map.iter().map(|(pos, &no)| {
                let ne = sub_map.get(pos).copied().unwrap_or(0)
                    + ins_map.get(pos).copied().unwrap_or(0)
                    + del_map.get(pos).copied().unwrap_or(0);
                (ne as f64 / no as f64, no as f64)
            }).collect();
            find_logit_shift(&bins, per_base_error_rate)
        };
        let delta_start = bins_for(&observed_from_start, &substitution_from_start, &insertion_from_start, &deletion_from_start);
        let delta_end   = bins_for(&observed_from_end,   &substitution_from_end,   &insertion_from_end,   &deletion_from_end);

        let normalize_bin = |no: u64, ns: u64, ni: u64, nd: u64, delta: f64| -> (f64, f64, f64, f64) {
            if no == 0 { return (0.0, 0.0, 0.0, 0.0); }
            let ne = ns + ni + nd;
            let raw = ne as f64 / no as f64;
            let norm_err = if raw <= 0.0 { 0.0 } else if raw >= 1.0 { 1.0 } else { sigmoid(logit(raw) + delta) };
            let (norm_sub, norm_ins, norm_del) = if ne > 0 {
                let s = norm_err / ne as f64;
                (ns as f64 * s, ni as f64 * s, nd as f64 * s)
            } else { (0.0, 0.0, 0.0) };
            (norm_sub, norm_ins, norm_del, norm_err)
        };

        let mut out = String::new();
        writeln!(out, "index,from_start,num_observed,num_substitution,num_insertion,num_deletion,normalized_substitution_rate,normalized_insertion_rate,normalized_deletion_rate,normalized_error_rate").unwrap();

        let mut start_positions: Vec<u32> = observed_from_start.keys().copied().collect();
        start_positions.sort();
        start_positions.dedup();
        for pos in start_positions {
            let no = observed_from_start.get(&pos).copied().unwrap_or(0);
            let ns = substitution_from_start.get(&pos).copied().unwrap_or(0);
            let ni = insertion_from_start.get(&pos).copied().unwrap_or(0);
            let nd = deletion_from_start.get(&pos).copied().unwrap_or(0);
            let (norm_sub, norm_ins, norm_del, norm_err) = normalize_bin(no, ns, ni, nd, delta_start);
            writeln!(out, "{},true,{},{},{},{},{:.6},{:.6},{:.6},{:.6}",
                pos, no, ns, ni, nd, norm_sub, norm_ins, norm_del, norm_err).unwrap();
        }

        let mut end_positions: Vec<u32> = observed_from_end.keys().copied().collect();
        end_positions.sort();
        end_positions.dedup();
        for pos in end_positions {
            let no = observed_from_end.get(&pos).copied().unwrap_or(0);
            let ns = substitution_from_end.get(&pos).copied().unwrap_or(0);
            let ni = insertion_from_end.get(&pos).copied().unwrap_or(0);
            let nd = deletion_from_end.get(&pos).copied().unwrap_or(0);
            let (norm_sub, norm_ins, norm_del, norm_err) = normalize_bin(no, ns, ni, nd, delta_end);
            writeln!(out, "{},false,{},{},{},{},{:.6},{:.6},{:.6},{:.6}",
                pos, no, ns, ni, nd, norm_sub, norm_ins, norm_del, norm_err).unwrap();
        }

        out
    }
}

impl PhredScoreSummary {
    /// `per_base_error_rate` is the global error rate from the hazard-ratio model
    /// (i.e. `1 - exp(-lambda)`) used to normalise the per-qscore rates via a logit-shift so
    /// that `sum(normalized_error_rate * num_observed) / sum(num_observed) == per_base_error_rate`.
    /// `normalized_error_rate` is always in `[0, 1]`.
    pub fn to_csv(&self, indices: Option<&[usize]>, per_base_error_rate: f64) -> String {
        use std::fmt::Write;
        let all: Vec<usize>;
        let indices = match indices {
            Some(idx) => idx,
            None => { all = (0..self.observed_per_key.len()).collect(); &all }
        };
        let mut observed: HashMap<u8, u64> = HashMap::new();
        let mut substitution: HashMap<u8, u64> = HashMap::new();
        let mut insertion: HashMap<u8, u64> = HashMap::new();
        let mut deletion: HashMap<u8, u64> = HashMap::new();
        for &i in indices {
            for (&q, &c) in &self.observed_per_key[i]     { *observed.entry(q).or_insert(0)     += c; }
            for (&q, &c) in &self.substitution_per_key[i] { *substitution.entry(q).or_insert(0) += c; }
            for (&q, &c) in &self.insertion_per_key[i]    { *insertion.entry(q).or_insert(0)    += c; }
            for (&q, &c) in &self.deletion_per_key[i]     { *deletion.entry(q).or_insert(0)     += c; }
        }

        // Build (raw_rate, weight) pairs for the logit-shift normalisation.
        let bins: Vec<(f64, f64)> = observed.iter().map(|(q, &no)| {
            let ne = substitution.get(q).copied().unwrap_or(0)
                + insertion.get(q).copied().unwrap_or(0)
                + deletion.get(q).copied().unwrap_or(0);
            (ne as f64 / no as f64, no as f64)
        }).collect();
        let delta = find_logit_shift(&bins, per_base_error_rate);

        let mut scores: Vec<u8> = observed.keys().copied().collect();
        scores.sort();
        scores.dedup();
        let mut out = String::new();
        writeln!(out, "qscore,empirical_qscore,num_observed,num_substitution,num_insertion,num_deletion,normalized_substitution_rate,normalized_insertion_rate,normalized_deletion_rate,normalized_error_rate").unwrap();
        for q in scores {
            let num_observed     = observed.get(&q).copied().unwrap_or(0);
            let num_substitution = substitution.get(&q).copied().unwrap_or(0);
            let num_insertion    = insertion.get(&q).copied().unwrap_or(0);
            let num_deletion     = deletion.get(&q).copied().unwrap_or(0);
            let num_error = num_substitution + num_insertion + num_deletion;
            let norm_err = if num_observed > 0 {
                let raw = num_error as f64 / num_observed as f64;
                if raw <= 0.0 { 0.0 } else if raw >= 1.0 { 1.0 } else { sigmoid(logit(raw) + delta) }
            } else { 0.0 };
            // Distribute norm_err proportionally among sub/ins/del.
            let (norm_sub, norm_ins, norm_del) = if num_error > 0 {
                let scale = norm_err / num_error as f64;
                (num_substitution as f64 * scale,
                 num_insertion    as f64 * scale,
                 num_deletion     as f64 * scale)
            } else { (0.0, 0.0, 0.0) };
            let empirical_q = if norm_err > 0.0 { -10.0 * norm_err.log10() } else { f64::INFINITY };
            writeln!(out, "{},{:.4},{},{},{},{},{:.6},{:.6},{:.6},{:.6}",
                q, empirical_q,
                num_observed, num_substitution, num_insertion, num_deletion,
                norm_sub, norm_ins, norm_del, norm_err,
            ).unwrap();
        }
        out
    }
}
