// This file contains multiple implementations from sylph (https://github.com/bluenote-1577/sylph). Below is their license.

/*
MIT License

Copyright (c) 2023 Jim Shaw

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

use std::collections::HashMap;
use std::fmt;

pub type Kmer = u64;

/**
 * A lookup table to convert a byte to a 2-bit sequence.
 * Adopted from https://github.com/bluenote-1577/sylph/blob/main/src/types.rs
 * 
 * A/a -> 0
 * C/c -> 1
 * G/g -> 2
 * T/t -> 3, U/u -> 3
 */
pub const BYTE_TO_SEQ: [u8; 256] = {
    let mut arr = [0u8; 256];

    arr[b'A' as usize] = 0;
    arr[b'C' as usize] = 1;
    arr[b'G' as usize] = 2;
    arr[b'T' as usize] = 3;
    arr[b'U' as usize] = 3;

    arr[b'a' as usize] = 0;
    arr[b'c' as usize] = 1;
    arr[b'g' as usize] = 2;
    arr[b't' as usize] = 3;
    arr[b'u' as usize] = 3;

    arr
};

pub const SEQ_TO_BYTE: [u8; 4] = [b'A', b'C', b'G', b'T'];
pub const SEQ_TO_CHAR: [char; 5] = ['A', 'C', 'G', 'T', 'N'];
// A -> T (3), C -> G (2), G -> C (1), T -> A (0), N -> N (4)
pub const SEQ_TO_COMPLEMENT_BIN: [u8; 5] = [3, 2, 1, 0, 4];


#[derive(Hash, PartialEq, Eq, Debug, Clone, Copy)]
pub enum EditOperation {
    /* SUBSTITUTION */
    AC,
    AG,
    AT,

    CA,
    CG,
    CT,

    GA,
    GC,
    GT,

    TA,
    TC,
    TG,

    /* INSERTION */
    _A,
    _C,
    _G,
    _T,

    /* DELETION */
    A_,
    C_,
    G_,
    T_,

    AMBIGUOUS, // when multiple operations can lead to the same neighbor
}

impl fmt::Display for EditOperation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EditOperation::AC => write!(f, "A>C"),
            EditOperation::AG => write!(f, "A>G"),
            EditOperation::AT => write!(f, "A>T"),

            EditOperation::CA => write!(f, "C>A"),
            EditOperation::CG => write!(f, "C>G"),
            EditOperation::CT => write!(f, "C>T"),

            EditOperation::GA => write!(f, "G>A"),
            EditOperation::GC => write!(f, "G>C"),
            EditOperation::GT => write!(f, "G>T"),

            EditOperation::TA => write!(f, "T>A"),
            EditOperation::TC => write!(f, "T>C"),
            EditOperation::TG => write!(f, "T>G"),

            EditOperation::_A => write!(f, "->A"),
            EditOperation::_C => write!(f, "->C"),
            EditOperation::_G => write!(f, "->G"),
            EditOperation::_T => write!(f, "->T"),

            EditOperation::A_ => write!(f, "A>-"),
            EditOperation::C_ => write!(f, "C>-"),
            EditOperation::G_ => write!(f, "G>-"),
            EditOperation::T_ => write!(f, "T>-"),

            EditOperation::AMBIGUOUS => write!(f, "AMBIGUOUS"),
        }
    }
}

// 2-D array to map (from, to) -> EditOperation
pub const BASES_TO_SUBSTITUTION: [[Option<EditOperation>; 4]; 4] = {
    let mut arr = [[None; 4]; 4];

    arr[0][1] = Some(EditOperation::AC);
    arr[0][2] = Some(EditOperation::AG);
    arr[0][3] = Some(EditOperation::AT);

    arr[1][0] = Some(EditOperation::CA);
    arr[1][2] = Some(EditOperation::CG);
    arr[1][3] = Some(EditOperation::CT);

    arr[2][0] = Some(EditOperation::GA);
    arr[2][1] = Some(EditOperation::GC);
    arr[2][3] = Some(EditOperation::GT);

    arr[3][0] = Some(EditOperation::TA);
    arr[3][1] = Some(EditOperation::TC);
    arr[3][2] = Some(EditOperation::TG);

    arr
};

// Use in the case we include both forward and reverse complements of the reads
pub const BASES_TO_SUBSTITUTION_CANONICAL: [[Option<EditOperation>; 4]; 4] = {
    let mut arr = [[None; 4]; 4];

    arr[0][1] = Some(EditOperation::TG);
    arr[0][2] = Some(EditOperation::TC);
    arr[0][3] = Some(EditOperation::TA);

    arr[1][0] = Some(EditOperation::CA);
    arr[1][2] = Some(EditOperation::CG);
    arr[1][3] = Some(EditOperation::CT);

    arr[2][0] = Some(EditOperation::CT);
    arr[2][1] = Some(EditOperation::CG);
    arr[2][3] = Some(EditOperation::CA);

    arr[3][0] = Some(EditOperation::TA);
    arr[3][1] = Some(EditOperation::TC);
    arr[3][2] = Some(EditOperation::TG);

    arr
};


pub const BASES_TO_INSERTION: [Option<EditOperation>; 4] = [
    Some(EditOperation::_A),
    Some(EditOperation::_C),
    Some(EditOperation::_G),
    Some(EditOperation::_T),
];

pub const BASES_TO_INSERTION_CANONICAL: [Option<EditOperation>; 4] = [
    Some(EditOperation::_T),
    Some(EditOperation::_G),
    Some(EditOperation::_G),
    Some(EditOperation::_T),
];

pub const BASES_TO_DELETION: [Option<EditOperation>; 4] = [
    Some(EditOperation::A_),
    Some(EditOperation::C_),
    Some(EditOperation::G_),
    Some(EditOperation::T_),
];

pub const BASES_TO_DELETION_CANONICAL: [Option<EditOperation>; 4] = [
    Some(EditOperation::T_),
    Some(EditOperation::G_),
    Some(EditOperation::G_),
    Some(EditOperation::T_),
];

pub const ALL_OPERATIONS: [EditOperation; 20] = [
    EditOperation::AC,
    EditOperation::AG,
    EditOperation::AT,

    EditOperation::GA,
    EditOperation::GC,
    EditOperation::GT,

    EditOperation::CA,
    EditOperation::CG,
    EditOperation::CT,

    EditOperation::TA,
    EditOperation::TC,
    EditOperation::TG,

    EditOperation::_A,
    EditOperation::_C,
    EditOperation::_G,
    EditOperation::_T,

    EditOperation::A_,
    EditOperation::C_,
    EditOperation::G_,
    EditOperation::T_,
];


pub const ALL_OPERATIONS_CANONICAL: [EditOperation; 10] = [
    EditOperation::CA,
    EditOperation::CG,
    EditOperation::CT,
    EditOperation::TA,
    EditOperation::TC,
    EditOperation::TG,

    EditOperation::_G,
    EditOperation::_T,

    EditOperation::G_,
    EditOperation::T_,
];

pub fn sbs96_str(op: &(EditOperation, u8, u8)) -> String {
    format!("{}[{}]{}", SEQ_TO_CHAR[op.1 as usize], op.0, SEQ_TO_CHAR[op.2 as usize])
}

/**
 * kv-mer statistics for downstream analysis.
 */
pub struct KVmerStats {
    pub k: u8,
    pub v: u8,

    pub keys: Vec<u64>,
    pub consensus_values: Vec<u64>,

    pub consensus_counts: Vec<u32>,
    pub total_counts: Vec<u32>,
    pub neighbor_counts: Vec<u32>,
    pub error_counts: Vec<HashMap<(EditOperation, u8, u8), u32>>,

    pub consensus_up_to_v_counts: Vec<Vec<u32>>,

    /// Quality-score calibration: for each Phred score, how many bases agreed
    /// with the consensus value (walking left-to-right, stopping at first mismatch).
    pub qscore_correct: HashMap<u8, u64>,
    /// Quality-score calibration: for each Phred score, how many bases were the
    /// first mismatch against the consensus (one per value observation at most).
    pub qscore_error: HashMap<u8, u64>,

    /// Per-key qscore correct counts (parallel to `keys`), enabling index-based filtering.
    pub qscore_correct_per_key: Vec<HashMap<u8, u64>>,
    /// Per-key qscore error counts (parallel to `keys`), enabling index-based filtering.
    pub qscore_error_per_key: Vec<HashMap<u8, u64>>,
}


#[derive(Clone)]
pub struct SequenceInfo {
    pub seq: Vec<u8>,
}
