use std::cell::RefCell;
use std::collections::HashMap;
use crate::types::{EditOperation, ALL_OPERATIONS, SEQ_TO_BYTE, SEQ_TO_CHAR, ValueInfo, NeighborInfo};
use crate::utils::_get_neighbors;
use crate::constants::EPSILON;

/// Estimate lambda from hazard ratios with beta held fixed.
/// Uses the cloglog linearisation: log(-log(1-h(t))) ≈ (beta-1)*log(t) + log(lambda*beta).
/// Positions where the denominator was zero (None) are skipped.
fn fit_lambda_given_beta(hazard_ratios: &[Option<f32>], beta: f32, k: usize) -> f32 {
    if beta <= 0.0 {
        return 0.0;
    }
    let intercepts: Vec<f32> = hazard_ratios.iter().enumerate()
        .filter_map(|(i, &hr)| hr.map(|h| {
            let y = (-(-(h.clamp(EPSILON, 1.0 - EPSILON))).ln_1p()).ln();
            let x = (i as f32 + k as f32).ln();
            y - (beta - 1.0) * x
        }))
        .collect();
    if intercepts.is_empty() {
        return 0.0;
    }
    let mean_intercept = intercepts.iter().sum::<f32>() / intercepts.len() as f32;
    (mean_intercept.exp() / beta).max(0.0)
}

fn confidence_interval(values: Option<&Vec<f32>>) -> (f32, f32) {
    let Some(values) = values else { return (0.0, 0.0); };
    if values.is_empty() {
        return (0.0, 0.0);
    }
    let mut values = values.clone();
    values.sort_by(f32::total_cmp);
    let n = values.len();
    (values[(n as f32 * 0.05) as usize], values[(n as f32 * 0.95) as usize])
}

/// Per-key error-rate statistics.
/// Corresponds to `KVmerStats` fields: `consensus_counts`, `total_counts`,
/// `neighbor_counts`, `consensus_up_to_v_counts`.
pub struct KVmerSummary {
    pub consensus_counts: Vec<u32>,
    pub total_counts: Vec<u32>,
    pub neighbor_counts: Vec<u32>,
    /// Per-key consensus counts for each value prefix length, indexed `[v-1][key_idx]`.
    pub consensus_up_to_v_counts: Vec<Vec<u32>>,
    pub key_strings: Vec<String>,
    pub value_strings: Vec<String>,
    pub homopolymer_lengths: Vec<u32>,
    pub error_counts_per_key: Vec<HashMap<NeighborInfo, u32>>,
    pub forward_error_counts_per_key: Vec<HashMap<NeighborInfo, u32>>,
    v: usize,
}

impl KVmerSummary {
    pub fn new(v: usize) -> Self {
        KVmerSummary {
            consensus_counts: Vec::new(),
            total_counts: Vec::new(),
            neighbor_counts: Vec::new(),
            consensus_up_to_v_counts: vec![Vec::new(); v],
            key_strings: Vec::new(),
            value_strings: Vec::new(),
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
        self.homopolymer_lengths.push(homopolymer_length);
        self.error_counts_per_key.push(error_count_map);
        self.forward_error_counts_per_key.push(forward_error_count_map);

        true
    }
}

impl KVmerSummary {
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

/// Phred quality-score calibration statistics.
/// Stores per-key correct/error counts indexed by (qscore, 0-based position in value).
/// The sequential scan stops at the first mismatch, so position p is only recorded
/// when all positions 0..p-1 matched consensus — implementing the hazard model.
pub struct PhredScoreSummary {
    /// Per-key correct counts indexed by (qscore, 0-based position in value).
    pub correct_pos_per_key: Vec<HashMap<(u8, u8), u64>>,
    /// Per-key error counts indexed by (qscore, 0-based position in value).
    pub error_pos_per_key: Vec<HashMap<(u8, u8), u64>>,
    pub bootstrap_lambdas: RefCell<HashMap<u8, Vec<f32>>>,
}

impl PhredScoreSummary {
    pub fn new() -> Self {
        PhredScoreSummary {
            correct_pos_per_key: Vec::new(),
            error_pos_per_key: Vec::new(),
            bootstrap_lambdas: RefCell::new(HashMap::new()),
        }
    }

    /// Accumulate one key's Phred calibration data across all value positions.
    pub fn update(&mut self, consensus: u64, value_size: u8, value_map: &HashMap<u64, Vec<ValueInfo>>) {
        let mut key_correct_pos: HashMap<(u8, u8), u64> = HashMap::new();
        let mut key_error_pos: HashMap<(u8, u8), u64> = HashMap::new();
        for (value, info_list) in value_map {
            for info in info_list {
                if info.qual.is_empty() {
                    continue;
                }
                for p in 0..value_size as usize {
                    let bit_shift = 2 * (value_size as usize - 1 - p);
                    let value_base     = (value     >> bit_shift) & 0b11;
                    let consensus_base = (consensus >> bit_shift) & 0b11;
                    let phred = info.qual[p].saturating_sub(33);
                    if value_base == consensus_base {
                        *key_correct_pos.entry((phred, p as u8)).or_insert(0) += 1;
                    } else {
                        *key_error_pos.entry((phred, p as u8)).or_insert(0) += 1;
                        break;
                    }
                }
            }
        }
        self.correct_pos_per_key.push(key_correct_pos);
        self.error_pos_per_key.push(key_error_pos);
    }

    pub fn clear_bootstrap_results(&self) {
        self.bootstrap_lambdas.borrow_mut().clear();
    }

    pub fn bootstrap_with_indices(&self, indices_sample: &[usize], estimated_lambda: f32, estimated_beta: f32, v_min: usize, v_max: usize, k: usize) {
        for (q, lambda) in self._lambda_with_indices(indices_sample, estimated_lambda, estimated_beta, v_min, v_max, k) {
            self.bootstrap_lambdas.borrow_mut().entry(q).or_default().push(lambda);
        }
    }

    fn _lambda_with_indices(&self, indices_sample: &[usize], estimated_lambda: f32, estimated_beta: f32, v_min: usize, v_max: usize, k: usize) -> HashMap<u8, f32> {
        let mut correct_pos: HashMap<(u8, u8), u64> = HashMap::new();
        let mut error_pos: HashMap<(u8, u8), u64> = HashMap::new();
        for &i in indices_sample {
            for (&key, &c) in &self.correct_pos_per_key[i] {
                *correct_pos.entry(key).or_insert(0) += c;
            }
            for (&key, &e) in &self.error_pos_per_key[i] {
                *error_pos.entry(key).or_insert(0) += e;
            }
        }

        let mut qscores: Vec<u8> = correct_pos.keys().chain(error_pos.keys())
            .map(|&(q, _)| q)
            .collect();
        qscores.sort_unstable();
        qscores.dedup();

        let mut lambdas = HashMap::new();
        for q in qscores {
            let mut hazard_ratios: Vec<Option<f32>> = Vec::new();
            for p in (v_min - 1)..v_max {
                let c = *correct_pos.get(&(q, p as u8)).unwrap_or(&0);
                let e = *error_pos.get(&(q, p as u8)).unwrap_or(&0);
                let total = c + e;
                hazard_ratios.push(if total > 0 { Some(e as f32 / total as f32) } else { None });
            }
            let lambda = fit_lambda_given_beta(&hazard_ratios, estimated_beta, k);
            lambdas.insert(q, if lambda.is_finite() { lambda } else { estimated_lambda });
        }
        lambdas
    }

    pub fn to_csv(&self, indices: &[usize], v_min: usize, v_max: usize, global_beta: f32, k: usize) -> String {
        let lambdas = self._lambda_with_indices(indices, 0.0, global_beta, v_min, v_max, k);
        let bootstrap_lambdas = self.bootstrap_lambdas.borrow();
        let mut csv = String::from("qscore,empirical_qscore,per_base_error_rate,per_base_error_rate_5-95th_percentile,num_correct,num_error\n");
        let mut qscores: Vec<u8> = lambdas.keys().copied().collect();
        qscores.sort_unstable();
        for q in qscores {
            let total_correct: u64 = (v_min - 1..v_max).map(|p| self.sum_pos(indices, &(q, p as u8), true)).sum();
            let total_error: u64 = (v_min - 1..v_max).map(|p| self.sum_pos(indices, &(q, p as u8), false)).sum();
            if total_correct + total_error == 0 {
                continue;
            }
            let lambda = lambdas[&q];
            let per_base_error_rate = 1.0 - (-(lambda as f64)).exp();
            let empirical_qscore = -10.0 * per_base_error_rate.log10();
            let ci = confidence_interval(bootstrap_lambdas.get(&q));
            let per_base_error_rate_ci = (1.0 - (-(ci.0 as f64)).exp(), 1.0 - (-(ci.1 as f64)).exp());
            csv.push_str(&format!(
                "{},{:.3},{:.6},{:.6}~{:.6},{},{}\n",
                q, empirical_qscore, per_base_error_rate, per_base_error_rate_ci.0, per_base_error_rate_ci.1, total_correct, total_error,
            ));
        }
        csv
    }

    fn sum_pos(&self, indices: &[usize], key: &(u8, u8), correct: bool) -> u64 {
        indices.iter().map(|&i| {
            let map = if correct { &self.correct_pos_per_key[i] } else { &self.error_pos_per_key[i] };
            map.get(key).copied().unwrap_or(0)
        }).sum()
    }
}

/// GC-content error calibration statistics.
/// Stores per-key correct/error counts indexed by (gc_content_percent, 0-based position in value).
/// Uses the same hazard model as PhredScoreSummary: position p is only recorded when all
/// positions 0..p-1 matched consensus.
pub struct GCContentSummary {
    /// Per-key correct counts indexed by (gc_content %, 0-based position in value).
    pub correct_pos_per_key: Vec<HashMap<(u8, u8), u64>>,
    /// Per-key error counts indexed by (gc_content %, 0-based position in value).
    pub error_pos_per_key: Vec<HashMap<(u8, u8), u64>>,
    pub bootstrap_lambdas: RefCell<HashMap<(u8, u8), Vec<f32>>>,
}

impl GCContentSummary {
    pub fn new() -> Self {
        GCContentSummary {
            correct_pos_per_key: Vec::new(),
            error_pos_per_key: Vec::new(),
            bootstrap_lambdas: RefCell::new(HashMap::new()),
        }
    }

    pub fn update(&mut self, consensus: u64, value_size: u8, value_map: &HashMap<u64, Vec<ValueInfo>>) {
        let mut key_correct_pos: HashMap<(u8, u8), u64> = HashMap::new();
        let mut key_error_pos: HashMap<(u8, u8), u64> = HashMap::new();
        for (value, info_list) in value_map {
            for info in info_list {
                if info.qual.is_empty() {
                    continue;
                }
                for p in 0..value_size as usize {
                    let bit_shift = 2 * (value_size as usize - 1 - p);
                    let value_base     = (value     >> bit_shift) & 0b11;
                    let consensus_base = (consensus >> bit_shift) & 0b11;
                    if value_base == consensus_base {
                        *key_correct_pos.entry((info.gc_content, p as u8)).or_insert(0) += 1;
                    } else {
                        *key_error_pos.entry((info.gc_content, p as u8)).or_insert(0) += 1;
                        break;
                    }
                }
            }
        }
        self.correct_pos_per_key.push(key_correct_pos);
        self.error_pos_per_key.push(key_error_pos);
    }

    pub fn clear_bootstrap_results(&self) {
        self.bootstrap_lambdas.borrow_mut().clear();
    }

    pub fn bootstrap_with_indices(&self, indices_sample: &[usize], estimated_lambda: f32, estimated_beta: f32, v_min: usize, v_max: usize, k: usize, step: u8) {
        for (bin, lambda) in self._lambda_with_indices(indices_sample, estimated_lambda, estimated_beta, v_min, v_max, k, step) {
            self.bootstrap_lambdas.borrow_mut().entry(bin).or_default().push(lambda);
        }
    }

    fn _lambda_with_indices(&self, indices_sample: &[usize], estimated_lambda: f32, estimated_beta: f32, v_min: usize, v_max: usize, k: usize, step: u8) -> HashMap<(u8, u8), f32> {
        let mut correct_pos: HashMap<(u8, u8), u64> = HashMap::new();
        let mut error_pos: HashMap<(u8, u8), u64> = HashMap::new();
        for &i in indices_sample {
            for (&key, &c) in &self.correct_pos_per_key[i] {
                *correct_pos.entry(key).or_insert(0) += c;
            }
            for (&key, &e) in &self.error_pos_per_key[i] {
                *error_pos.entry(key).or_insert(0) += e;
            }
        }

        let mut lambdas = HashMap::new();
        let mut gc_min: u8 = 0;
        while gc_min <= 100 {
            let gc_max = gc_min.saturating_add(step).min(101);

            let mut pos_correct: HashMap<u8, u64> = HashMap::new();
            let mut pos_error: HashMap<u8, u64> = HashMap::new();
            for gc in gc_min..gc_max {
                for p in (v_min - 1)..v_max {
                    if let Some(&c) = correct_pos.get(&(gc, p as u8)) {
                        *pos_correct.entry(p as u8).or_insert(0) += c;
                    }
                    if let Some(&e) = error_pos.get(&(gc, p as u8)) {
                        *pos_error.entry(p as u8).or_insert(0) += e;
                    }
                }
            }

            let mut hazard_ratios: Vec<Option<f32>> = Vec::new();
            let mut total_correct = 0u64;
            let mut total_error = 0u64;
            for p in (v_min - 1)..v_max {
                let c = pos_correct.get(&(p as u8)).copied().unwrap_or(0);
                let e = pos_error.get(&(p as u8)).copied().unwrap_or(0);
                total_correct += c;
                total_error += e;
                let tot = c + e;
                hazard_ratios.push(if tot > 0 { Some(e as f32 / tot as f32) } else { None });
            }

            if total_correct + total_error > 0 {
                let lambda = fit_lambda_given_beta(&hazard_ratios, estimated_beta, k);
                lambdas.insert((gc_min, gc_max), if lambda.is_finite() { lambda } else { estimated_lambda });
            }

            if gc_max == 101 { break; }
            gc_min = gc_max;
        }
        lambdas
    }

    pub fn to_csv(&self, indices: &[usize], v_min: usize, v_max: usize, global_beta: f32, k: usize, step: u8) -> String {
        let lambdas = self._lambda_with_indices(indices, 0.0, global_beta, v_min, v_max, k, step);
        let bootstrap_lambdas = self.bootstrap_lambdas.borrow();
        let mut csv = String::from("gc_content_min,gc_content_max_exclusive,per_base_error_rate,per_base_error_rate_5-95th_percentile,num_correct,num_error\n");
        let mut bins: Vec<(u8, u8)> = lambdas.keys().copied().collect();
        bins.sort_unstable();
        for (gc_min, gc_max) in bins {
            let mut total_correct = 0u64;
            let mut total_error = 0u64;
            for gc in gc_min..gc_max {
                for p in (v_min - 1)..v_max {
                    total_correct += self.sum_pos(indices, &(gc, p as u8), true);
                    total_error += self.sum_pos(indices, &(gc, p as u8), false);
                }
            }
            let lambda = lambdas[&(gc_min, gc_max)];
            let per_base_error_rate = 1.0 - (-(lambda as f64)).exp();
            let lambda_ci = confidence_interval(bootstrap_lambdas.get(&(gc_min, gc_max)));
            let per_base_error_rate_ci = (1.0 - (-(lambda_ci.0 as f64)).exp(), 1.0 - (-(lambda_ci.1 as f64)).exp());
            csv.push_str(&format!(
                "{},{},{:.6},{:.6}~{:.6},{},{}\n",
                gc_min, gc_max, per_base_error_rate, per_base_error_rate_ci.0, per_base_error_rate_ci.1, total_correct, total_error,
            ));
        }
        csv
    }

    fn sum_pos(&self, indices: &[usize], key: &(u8, u8), correct: bool) -> u64 {
        indices.iter().map(|&i| {
            let map = if correct { &self.correct_pos_per_key[i] } else { &self.error_pos_per_key[i] };
            map.get(key).copied().unwrap_or(0)
        }).sum()
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

    pub fn update(&mut self, consensus: u64, value_size: u8, value_map: &HashMap<u64, Vec<ValueInfo>>) {
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
        writeln!(out, "index,from_start,num_correct,num_error,error_rate").unwrap();

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
