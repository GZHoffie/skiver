use std::collections::HashMap;
use crate::types::{EditOperation, ALL_OPERATIONS, SEQ_TO_BYTE, SEQ_TO_CHAR, ValueInfo, NeighborInfo};
use crate::utils::_get_neighbors;
use log::info;


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

/// Read-position error calibration statistics.
/// Stores per-key counts of correct/erroneous bases indexed by position from the
/// start or end of the read.
pub struct ReadPositionSummary {
    pub correct_from_start_per_key: Vec<HashMap<u32, u64>>,
    pub correct_from_end_per_key: Vec<HashMap<u32, u64>>,
    pub error_from_start_per_key: Vec<HashMap<u32, u64>>,
    pub error_from_end_per_key: Vec<HashMap<u32, u64>>,
}

impl ReadPositionSummary {
    pub fn new() -> Self {
        ReadPositionSummary {
            correct_from_start_per_key: Vec::new(),
            correct_from_end_per_key: Vec::new(),
            error_from_start_per_key: Vec::new(),
            error_from_end_per_key: Vec::new(),
        }
    }

    /// If `last_base_only` is true, only the last base of each value is considered.
    pub fn update(&mut self, consensus: u64, value_size: u8, value_map: &HashMap<u64, Vec<ValueInfo>>, last_base_only: bool) {
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

                    let record = !last_base_only || (p == (value_size as usize - 1) && last_base_only);
                    if value_base == consensus_base {
                        if record {
                            *correct_from_start.entry(pos_from_start).or_insert(0) += 1;
                            *correct_from_end.entry(pos_from_end).or_insert(0) += 1;
                        }
                    } else {
                        if record {
                            *error_from_start.entry(pos_from_start).or_insert(0) += 1;
                            *error_from_end.entry(pos_from_end).or_insert(0) += 1;
                        }
                        break;
                    }
                    
                }
            }
        }
        self.correct_from_start_per_key.push(correct_from_start);
        self.correct_from_end_per_key.push(correct_from_end);
        self.error_from_start_per_key.push(error_from_start);
        self.error_from_end_per_key.push(error_from_end);
    }
}

impl ReadPositionSummary {
    pub fn to_csv(&self, indices: Option<&[usize]>) -> String {
        use std::fmt::Write;
        let all: Vec<usize>;
        let indices = match indices {
            Some(idx) => idx,
            None => { all = (0..self.correct_from_start_per_key.len()).collect(); &all }
        };
        let mut correct_from_start: HashMap<u32, u64> = HashMap::new();
        let mut correct_from_end: HashMap<u32, u64> = HashMap::new();
        let mut error_from_start: HashMap<u32, u64> = HashMap::new();
        let mut error_from_end: HashMap<u32, u64> = HashMap::new();

        for &i in indices {
            for (&pos, &c) in &self.correct_from_start_per_key[i] { *correct_from_start.entry(pos).or_insert(0) += c; }
            for (&pos, &c) in &self.correct_from_end_per_key[i]   { *correct_from_end.entry(pos).or_insert(0) += c; }
            for (&pos, &e) in &self.error_from_start_per_key[i]   { *error_from_start.entry(pos).or_insert(0) += e; }
            for (&pos, &e) in &self.error_from_end_per_key[i]     { *error_from_end.entry(pos).or_insert(0) += e; }
        }

        let mut out = String::new();
        writeln!(out, "index,from_start,num_correct,num_error,hazard_rate").unwrap();

        let mut start_positions: Vec<u32> = correct_from_start.keys().chain(error_from_start.keys()).copied().collect();
        start_positions.sort();
        start_positions.dedup();
        for pos in start_positions {
            let nc = correct_from_start.get(&pos).copied().unwrap_or(0);
            let ne = error_from_start.get(&pos).copied().unwrap_or(0);
            let error_rate = if nc + ne > 0 { ne as f64 / (nc + ne) as f64 } else { 0.0 };
            writeln!(out, "{},true,{},{},{:.6}", pos, nc, ne, error_rate).unwrap();
        }

        let mut end_positions: Vec<u32> = correct_from_end.keys().chain(error_from_end.keys()).copied().collect();
        end_positions.sort();
        end_positions.dedup();
        for pos in end_positions {
            let nc = correct_from_end.get(&pos).copied().unwrap_or(0);
            let ne = error_from_end.get(&pos).copied().unwrap_or(0);
            let error_rate = if nc + ne > 0 { ne as f64 / (nc + ne) as f64 } else { 0.0 };
            writeln!(out, "{},false,{},{},{:.6}", pos, nc, ne, error_rate).unwrap();
        }

        out
    }
}

impl PhredScoreSummary {
    /// Returns `(prop_sub, prop_ins, prop_del)` — the fraction of each error type among all
    /// errors observed for the selected indices.  Falls back to equal thirds if no errors seen.
    pub fn error_proportions(&self, indices: Option<&[usize]>) -> (f64, f64, f64) {
        let all: Vec<usize>;
        let indices = match indices {
            Some(idx) => idx,
            None => { all = (0..self.substitution_per_key.len()).collect(); &all }
        };
        let total_sub: u64 = indices.iter().flat_map(|&i| self.substitution_per_key[i].values().copied()).sum();
        let total_ins: u64 = indices.iter().flat_map(|&i| self.insertion_per_key[i].values().copied()).sum();
        let total_del: u64 = indices.iter().flat_map(|&i| self.deletion_per_key[i].values().copied()).sum();
        let total = total_sub + total_ins + total_del;
        if total == 0 {
            return (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0);
        }
        (total_sub as f64 / total as f64,
         total_ins as f64 / total as f64,
         total_del as f64 / total as f64)
    }

    /// Estimates per-quality-score error rates using a Bayesian approach:
    ///   Pr[error_type | Q] = Pr[Q | error_type] * Pr[error_type] / Pr[Q]
    /// where Pr[Q | error_type] is estimated from the accumulated counts, Pr[Q] from `observed`,
    /// and Pr[error_type] = per_base_error_rate * proportion_of_that_error_type.
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

        let total_sub: u64 = substitution.values().sum();
        let total_ins: u64 = insertion.values().sum();
        let total_del: u64 = deletion.values().sum();
        let total_errors = total_sub + total_ins + total_del;
        let total_observed: u64 = observed.values().sum();

        // Pr[error_type] = per_base_error_rate * proportion_t
        let (prop_sub, prop_ins, prop_del) = if total_errors > 0 {
            (total_sub as f64 / total_errors as f64,
             total_ins as f64 / total_errors as f64,
             total_del as f64 / total_errors as f64)
        } else {
            (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
        };
        let prior_sub = per_base_error_rate * prop_sub;
        let prior_ins = per_base_error_rate * prop_ins;
        let prior_del = per_base_error_rate * prop_del;

        let mut scores: Vec<u8> = observed.keys().copied().collect();
        scores.sort();
        scores.dedup();
        let mut out = String::new();
        writeln!(out, "qscore,empirical_qscore,num_observed,num_substitution,num_insertion,num_deletion,bayesian_substitution_rate,bayesian_insertion_rate,bayesian_deletion_rate,bayesian_error_rate").unwrap();
        for q in scores {
            let num_observed = observed.get(&q).copied().unwrap_or(0);
            let num_sub      = substitution.get(&q).copied().unwrap_or(0);
            let num_ins      = insertion.get(&q).copied().unwrap_or(0);
            let num_del      = deletion.get(&q).copied().unwrap_or(0);

            // Pr[Q] and Pr[Q | error_type]
            let pr_q           = if total_observed > 0 { num_observed as f64 / total_observed as f64 } else { 0.0 };
            let pr_q_given_sub = if total_sub > 0 { num_sub as f64 / total_sub as f64 } else { 0.0 };
            let pr_q_given_ins = if total_ins > 0 { num_ins as f64 / total_ins as f64 } else { 0.0 };
            let pr_q_given_del = if total_del > 0 { num_del as f64 / total_del as f64 } else { 0.0 };

            // Pr[error_type | Q] = Pr[Q | error_type] * Pr[error_type] / Pr[Q]
            let (norm_sub, norm_ins, norm_del, norm_err) = if pr_q > 0.0 {
                let ns = (pr_q_given_sub * prior_sub / pr_q).min(1.0);
                let ni = (pr_q_given_ins * prior_ins / pr_q).min(1.0);
                let nd = (pr_q_given_del * prior_del / pr_q).min(1.0);
                let ne = (ns + ni + nd).min(1.0);
                (ns, ni, nd, ne)
            } else {
                (0.0, 0.0, 0.0, 0.0)
            };

            let empirical_q = if norm_err > 0.0 { -10.0 * norm_err.log10() } else { f64::INFINITY };
            writeln!(out, "{},{:.4},{},{},{},{},{:.6},{:.6},{:.6},{:.6}",
                q, empirical_q,
                num_observed, num_sub, num_ins, num_del,
                norm_sub, norm_ins, norm_del, norm_err,
            ).unwrap();
        }
        out
    }
}
