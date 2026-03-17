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

use crate::types::*;
use std::collections::HashSet;
use std::collections::HashMap;

#[inline]
pub fn mm_hash64_masked(kmer: u64, mask: Option<u64>) -> u64 {
    //TODO this is bugged. Fix after release
    let mut key = kmer;
    if let Some(mask) = mask {
        key = kmer & mask;
    }
    key = !key.wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    key
}


pub fn decode(byte: u64) -> u8 {
    if byte == 0 {
        return b'A';
    } else if byte == 1 {
        return b'C';
    } else if byte == 2 {
        return b'G';
    } else if byte == 3 {
        return b'T';
    } else {
        panic!("decoding failed")
    }
}
pub fn print_string(kmer: u64, k: usize) {
    let mut bytes = vec![];
    let mask = 3;
    for i in 0..k {
        let val = kmer >> 2 * i;
        let val = val & mask;
        bytes.push(decode(val));
    }
    dbg!(std::str::from_utf8(&bytes.into_iter().rev().collect::<Vec<u8>>()).unwrap());
}
#[inline]
fn _position_min<T: Ord>(slice: &[T]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value1.cmp(value0))
        .map(|(idx, _)| idx)
}

pub fn fmh_seeds_masked(
    string: &[u8],
    keys_vec: &mut Vec<u64>,
    values_vec: &mut Vec<u64>,
    value_info_vec: &mut Vec<ValueInfo>,
    c: usize,
    k: usize,
    v: usize,
    bidirectional: bool,
) {
    type MarkerBits = u64;
    if string.len() < k + v {
        return;
    }

    let mut rolling_key_f: MarkerBits = 0;
    let mut rolling_key_r: MarkerBits = 0;
    let mut rolling_value_f: MarkerBits = 0;
    let mut rolling_value_r: MarkerBits = 0;


    let key_mask = (1u64 << (2 * k)) - 1;
    let value_mask = (1u64 << (2 * v)) - 1;

    let len = string.len();

    let threshold_marker = u64::MAX / (c as u64);

    // Initialize keys
    for i in 0..k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        let _nuc_r = 3 - nuc_f;
        rolling_key_f <<= 2;
        rolling_key_f |= nuc_f;
    }

    // Initialize values
    for i in 0..v {
        let nuc_f = BYTE_TO_SEQ[string[i + k - 1] as usize] as u64;
        let _nuc_r = 3 - nuc_f;
        rolling_value_f <<= 2;
        rolling_value_f |= nuc_f;
    }

    // initialize for reverse complement
    if bidirectional {
        for i in 0..v - 1 {
            let nuc_r = 3 - BYTE_TO_SEQ[string[i] as usize] as u64;
            rolling_value_r >>= 2;
            rolling_value_r |= nuc_r << (2 * (v - 1));
        }

        for i in 0..k {
            let nuc_r = 3 - BYTE_TO_SEQ[string[i + v - 1] as usize] as u64;
            rolling_key_r >>= 2;
            rolling_key_r |= nuc_r << (2 * (k - 1));
        }
    }

    // Iterate through the string
    for i in k + v - 1..len {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        

        // append the first base of the value to the key
        let first_base_v = (rolling_value_f >> (2 * (v - 1))) & 0b11;

        rolling_key_f <<= 2;
        rolling_key_f |= first_base_v;
        rolling_key_f &= key_mask;

        rolling_value_f <<= 2;
        rolling_value_f |= nuc_f;
        rolling_value_f &= value_mask;

        let hash_f = mm_hash64_masked(rolling_key_f, None);

        if hash_f < threshold_marker {
            keys_vec.push(rolling_key_f as u64);
            values_vec.push(rolling_value_f as u64);
            value_info_vec.push(ValueInfo { qual: vec![], start_index: (i - v + 1) as u32, dist_to_read_end: len as u32 - (i - v + 1) as u32, is_forward: true });
        }

        if bidirectional {
            //let nuc_value_r = 3 - BYTE_TO_SEQ[string[i + v - 1] as usize] as u64;
            //let nuc_key_r = 3 - BYTE_TO_SEQ[string[i + k + v - 1] as usize] as u64;
            let nuc_r = 3 - nuc_f;

            let last_base_k = rolling_key_r & 0b11;

            rolling_value_r >>= 2;
            rolling_value_r |= last_base_k << (2 * (v - 1));
            rolling_value_r &= value_mask;

            rolling_key_r >>= 2;
            rolling_key_r |= nuc_r << (2 * (k - 1));
            rolling_key_r &= key_mask;



            let hash_r = mm_hash64_masked(rolling_key_r, None);
            if hash_r < threshold_marker {
                keys_vec.push(rolling_key_r as u64);
                values_vec.push(rolling_value_r as u64);
                value_info_vec.push(ValueInfo { qual: vec![], start_index: (i - k + 1) as u32, dist_to_read_end: len as u32 - (i - k + 1) as u32, is_forward: false });
            }
        }
    }
}

/// Like `fmh_seeds_masked`, but also records the Phred quality scores for each
/// selected value.
///
/// `quals_vec` receives one entry per selected k,v-mer:
///   - Forward strand: `qual[i-v+1 ..= i]`  (v bytes, in value-position order)
///   - RC strand:      `[qual[i-k], qual[i-k-1], ..., qual[i-k-v+1]]`
///                     (reversed, so index 0 = quality of RC value position 0)
///
/// Callers must ensure `qual.len() == string.len()`.
pub fn fmh_seeds_masked_with_qual(
    string: &[u8],
    qual: &[u8],
    keys_vec: &mut Vec<u64>,
    values_vec: &mut Vec<u64>,
    value_info_vec: &mut Vec<ValueInfo>,
    c: usize,
    k: usize,
    v: usize,
    bidirectional: bool,
) {
    type MarkerBits = u64;
    if string.len() < k + v {
        return;
    }
    debug_assert_eq!(string.len(), qual.len());

    let mut rolling_key_f: MarkerBits = 0;
    let mut rolling_key_r: MarkerBits = 0;
    let mut rolling_value_f: MarkerBits = 0;
    let mut rolling_value_r: MarkerBits = 0;

    let key_mask = (1u64 << (2 * k)) - 1;
    let value_mask = (1u64 << (2 * v)) - 1;
    let len = string.len();
    let threshold_marker = u64::MAX / (c as u64);

    // Initialize keys (k-1 bases)
    for i in 0..k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        rolling_key_f <<= 2;
        rolling_key_f |= nuc_f;
    }

    // Initialize values (v bases)
    for i in 0..v {
        let nuc_f = BYTE_TO_SEQ[string[i + k - 1] as usize] as u64;
        rolling_value_f <<= 2;
        rolling_value_f |= nuc_f;
    }

    // Initialize RC
    if bidirectional {
        for i in 0..v - 1 {
            let nuc_r = 3 - BYTE_TO_SEQ[string[i] as usize] as u64;
            rolling_value_r >>= 2;
            rolling_value_r |= nuc_r << (2 * (v - 1));
        }
        for i in 0..k {
            let nuc_r = 3 - BYTE_TO_SEQ[string[i + v - 1] as usize] as u64;
            rolling_key_r >>= 2;
            rolling_key_r |= nuc_r << (2 * (k - 1));
        }
    }

    for i in k + v - 1..len {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;

        // Advance forward key/value
        let first_base_v = (rolling_value_f >> (2 * (v - 1))) & 0b11;
        rolling_key_f <<= 2;
        rolling_key_f |= first_base_v;
        rolling_key_f &= key_mask;
        rolling_value_f <<= 2;
        rolling_value_f |= nuc_f;
        rolling_value_f &= value_mask;

        let hash_f = mm_hash64_masked(rolling_key_f, None);
        if hash_f < threshold_marker {
            keys_vec.push(rolling_key_f);
            values_vec.push(rolling_value_f);
            // Value covers read positions [i-v+1, i]; position 0 = i-v+1.
            value_info_vec.push(ValueInfo { qual: qual[i - v + 1..=i].to_vec(), start_index: (i - v + 1) as u32, dist_to_read_end: len as u32 - (i - v + 1) as u32, is_forward: true });
        }

        if bidirectional {
            let nuc_r = 3 - nuc_f;
            let last_base_k = rolling_key_r & 0b11;

            rolling_value_r >>= 2;
            rolling_value_r |= last_base_k << (2 * (v - 1));
            rolling_value_r &= value_mask;

            rolling_key_r >>= 2;
            rolling_key_r |= nuc_r << (2 * (k - 1));
            rolling_key_r &= key_mask;

            let hash_r = mm_hash64_masked(rolling_key_r, None);
            if hash_r < threshold_marker {
                keys_vec.push(rolling_key_r);
                values_vec.push(rolling_value_r);
                // RC value position p corresponds to forward read position (i-k-p).
                // Quality string is in RC-value-position order: p=0 → qual[i-k].
                let rc_qual: Vec<u8> = (0..v).map(|p| qual[i - k - p]).collect();
                value_info_vec.push(ValueInfo { qual: rc_qual, start_index: (i - k + 1) as u32, dist_to_read_end: len as u32 - (i - k + 1) as u32, is_forward: false });
            }
        }
    }
}

pub fn count_seeds_in_set(
    string: &[u8],
    k: usize,
    kmer_count: &mut HashMap<u64, u32>,
    kmer_set: &HashSet<u64>,
    bidirectional: bool,
) {
    type MarkerBits = u64;
    if string.len() < k {
        return;
    }

    let marker_k = k;
    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_r_marker: MarkerBits = 0;

    let marker_reverse_shift_dist = 2 * (marker_k - 1);
    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));
    let len = string.len();
    //    let threshold = i64::MIN + (u64::MAX / (c as u64)) as i64;
    //    let threshold_marker = i64::MIN + (u64::MAX / sketch_params.marker_c as u64) as i64;

    for i in 0..marker_k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        //        let nuc_f = KmerEnc::encode(string[i]
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        //        rolling_kmer_r = KmerEnc::rc(rolling_kmer_f, k);
        if bidirectional {
            rolling_kmer_r_marker >>= 2;
            rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
        }
    }
    for i in marker_k-1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_f_marker &= marker_mask;

        if kmer_set.contains(&rolling_kmer_f_marker) {
            *kmer_count.entry(rolling_kmer_f_marker).or_insert(0) += 1;
        }
        if bidirectional {
            rolling_kmer_r_marker >>= 2;
            rolling_kmer_r_marker &= marker_rev_mask;
            rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;    

            if kmer_set.contains(&rolling_kmer_r_marker) {
                *kmer_count.entry(rolling_kmer_r_marker).or_insert(0) += 1;
            }
        }
    }
}

        
