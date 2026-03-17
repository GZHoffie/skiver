use std::fmt;
use std::collections::HashMap;
use crate::types::{EditOperation, ALL_OPERATIONS, SEQ_TO_CHAR};

/// Per-key error-rate statistics.
/// Corresponds to `KVmerStats` fields: `consensus_counts`, `total_counts`,
/// `neighbor_counts`, `consensus_up_to_v_counts`.
pub struct ErrorSummary {
    pub consensus_counts: Vec<u32>,
    pub total_counts: Vec<u32>,
    pub neighbor_counts: Vec<u32>,
    /// Per-key consensus counts for each value prefix length, indexed `[v-1][key_idx]`.
    pub consensus_up_to_v_counts: Vec<Vec<u32>>,
    v: usize,
}

impl ErrorSummary {
    pub fn new(v: usize) -> Self {
        ErrorSummary {
            consensus_counts: Vec::new(),
            total_counts: Vec::new(),
            neighbor_counts: Vec::new(),
            consensus_up_to_v_counts: vec![Vec::new(); v],
            v,
        }
    }

    /// Accumulate one key's error statistics.
    /// `per_v_consensus[j]` is the consensus count for value prefix length `j + 1`.
    pub fn update(
        &mut self,
        consensus_count: u32,
        total_count: u32,
        neighbor_count: u32,
        per_v_consensus: &[u32],
    ) {
        self.consensus_counts.push(consensus_count);
        self.total_counts.push(total_count);
        self.neighbor_counts.push(neighbor_count);
        for (j, &c) in per_v_consensus.iter().enumerate() {
            if j < self.consensus_up_to_v_counts.len() {
                self.consensus_up_to_v_counts[j].push(c);
            }
        }
    }
}

impl fmt::Display for ErrorSummary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let total_consensus: u64 = self.consensus_counts.iter().map(|&c| c as u64).sum();
        let total_count: u64 = self.total_counts.iter().map(|&c| c as u64).sum();
        let total_neighbor: u64 = self.neighbor_counts.iter().map(|&c| c as u64).sum();
        writeln!(f, "num_keys,total_consensus,total_count,total_neighbor")?;
        writeln!(f, "{},{},{},{}", self.consensus_counts.len(), total_consensus, total_count, total_neighbor)?;
        writeln!(f, "v,consensus_up_to_v")?;
        for j in 0..self.v {
            let sum: u64 = self.consensus_up_to_v_counts[j].iter().map(|&c| c as u64).sum();
            writeln!(f, "{},{}", j + 1, sum)?;
        }
        Ok(())
    }
}

/// Per-key error-type spectrum statistics.
/// Corresponds to `KVmerStats` field: `error_counts`.
pub struct ErrorSpectrumSummary {
    pub error_counts: Vec<HashMap<(EditOperation, u8, u8), u32>>,
}

impl ErrorSpectrumSummary {
    pub fn new() -> Self {
        ErrorSpectrumSummary {
            error_counts: Vec::new(),
        }
    }

    /// Accumulate one key's per-operation error counts.
    pub fn update(&mut self, error_map: HashMap<(EditOperation, u8, u8), u32>) {
        self.error_counts.push(error_map);
    }
}

impl fmt::Display for ErrorSpectrumSummary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut totals: HashMap<(EditOperation, u8, u8), u64> = HashMap::new();
        for map in &self.error_counts {
            for (&key, &count) in map {
                *totals.entry(key).or_insert(0) += count as u64;
            }
        }
        writeln!(f, "operation,prev_base,next_base,count")?;
        for &op in ALL_OPERATIONS.iter() {
            for prev_base in 0u8..4 {
                for next_base in 0u8..4 {
                    let count = totals.get(&(op, prev_base, next_base)).copied().unwrap_or(0);
                    if count > 0 {
                        writeln!(
                            f,
                            "{},{},{},{}",
                            op,
                            SEQ_TO_CHAR[prev_base as usize],
                            SEQ_TO_CHAR[next_base as usize],
                            count
                        )?;
                    }
                }
            }
        }
        Ok(())
    }
}

/// Phred quality-score calibration statistics.
/// Corresponds to `KVmerStats` fields: `qscore_correct`, `qscore_error`,
/// `qscore_correct_per_key`, `qscore_error_per_key`.
pub struct PhredScoreSummary {
    pub correct: HashMap<u8, u64>,
    pub error: HashMap<u8, u64>,
    pub correct_per_key: Vec<HashMap<u8, u64>>,
    pub error_per_key: Vec<HashMap<u8, u64>>,
}

impl PhredScoreSummary {
    pub fn new() -> Self {
        PhredScoreSummary {
            correct: HashMap::new(),
            error: HashMap::new(),
            correct_per_key: Vec::new(),
            error_per_key: Vec::new(),
        }
    }

    /// Accumulate one key's Phred calibration data.
    pub fn update(&mut self, key_correct: HashMap<u8, u64>, key_error: HashMap<u8, u64>) {
        for (&q, &c) in &key_correct { *self.correct.entry(q).or_insert(0) += c; }
        for (&q, &e) in &key_error   { *self.error.entry(q).or_insert(0)   += e; }
        self.correct_per_key.push(key_correct);
        self.error_per_key.push(key_error);
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

    pub fn update(
        &mut self,
        correct_from_start: HashMap<u32, u64>,
        correct_from_end: HashMap<u32, u64>,
        error_from_start: HashMap<u32, u64>,
        error_from_end: HashMap<u32, u64>,
    ) {
        self.correct_from_start_per_key.push(correct_from_start);
        self.correct_from_end_per_key.push(correct_from_end);
        self.error_from_start_per_key.push(error_from_start);
        self.error_from_end_per_key.push(error_from_end);
    }
}

impl fmt::Display for ReadPositionSummary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut correct_from_start: HashMap<u32, u64> = HashMap::new();
        let mut correct_from_end: HashMap<u32, u64> = HashMap::new();
        let mut error_from_start: HashMap<u32, u64> = HashMap::new();
        let mut error_from_end: HashMap<u32, u64> = HashMap::new();

        for map in &self.correct_from_start_per_key {
            for (&pos, &c) in map { *correct_from_start.entry(pos).or_insert(0) += c; }
        }
        for map in &self.correct_from_end_per_key {
            for (&pos, &c) in map { *correct_from_end.entry(pos).or_insert(0) += c; }
        }
        for map in &self.error_from_start_per_key {
            for (&pos, &e) in map { *error_from_start.entry(pos).or_insert(0) += e; }
        }
        for map in &self.error_from_end_per_key {
            for (&pos, &e) in map { *error_from_end.entry(pos).or_insert(0) += e; }
        }

        writeln!(f, "index,from_start,num_correct,num_error")?;

        let mut start_positions: Vec<u32> = correct_from_start.keys().chain(error_from_start.keys()).copied().collect();
        start_positions.sort();
        start_positions.dedup();
        for pos in start_positions {
            writeln!(f, "{},true,{},{}", pos,
                correct_from_start.get(&pos).copied().unwrap_or(0),
                error_from_start.get(&pos).copied().unwrap_or(0))?;
        }

        let mut end_positions: Vec<u32> = correct_from_end.keys().chain(error_from_end.keys()).copied().collect();
        end_positions.sort();
        end_positions.dedup();
        for pos in end_positions {
            writeln!(f, "{},false,{},{}", pos,
                correct_from_end.get(&pos).copied().unwrap_or(0),
                error_from_end.get(&pos).copied().unwrap_or(0))?;
        }

        Ok(())
    }
}

impl fmt::Display for PhredScoreSummary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "phred,correct,error")?;
        let mut scores: Vec<u8> = self.correct.keys().chain(self.error.keys()).copied().collect();
        scores.sort();
        scores.dedup();
        for q in scores {
            writeln!(
                f,
                "{},{},{}",
                q,
                self.correct.get(&q).copied().unwrap_or(0),
                self.error.get(&q).copied().unwrap_or(0),
            )?;
        }
        Ok(())
    }
}
