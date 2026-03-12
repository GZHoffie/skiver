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

use std::arch::x86_64::*;
use crate::types::*;

/**
 * A fast hash function for 64-bit kmers using AVX2 instructions.
 */
#[inline]
#[target_feature(enable = "avx2")]
pub unsafe fn mm_hash256_masked(kmer: __m256i, mask: Option<i64>) -> __m256i {
    // mask the kmer so that only the masked bits are used in the hash
    let mut key = kmer;
    if let Some(mask) = mask {
        let mask_vec = _mm256_set_epi64x(mask, mask, mask, mask);
        key = _mm256_and_si256(kmer, mask_vec);
    }

    let s1 = _mm256_slli_epi64(key, 21);
    key = _mm256_add_epi64(key, s1);
    
    key = _mm256_xor_si256(key, _mm256_cmpeq_epi64(key, key));

    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 24));
    let s2 = _mm256_slli_epi64(key, 3);
    let s3 = _mm256_slli_epi64(key, 8);

    key = _mm256_add_epi64(key, s2);
    key = _mm256_add_epi64(key, s3);
    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 14));
    let s4 = _mm256_slli_epi64(key, 2);
    let s5 = _mm256_slli_epi64(key, 4);
    key = _mm256_add_epi64(key, s4);
    key = _mm256_add_epi64(key, s5);
    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 28));

    let s6 = _mm256_slli_epi64(key, 31);
    key = _mm256_add_epi64(key, s6);

    return key;
}

/**
 * Shift each 64-bit lane in a __m256i by 2*k - 2 bits to the left.
 * This is used to update the reverse kmer in the rolling hash.
 */
#[target_feature(enable = "avx2")]
pub unsafe fn _shift_mm256_left_by_k(kmer: __m256i, k: usize) -> __m256i {
    // shift left by 2*k - 2
    let shifted = match k {
        1 => { _mm256_slli_epi64(kmer, 0) }
        2 => { _mm256_slli_epi64(kmer, 2) }
        3 => { _mm256_slli_epi64(kmer, 4) }
        4 => { _mm256_slli_epi64(kmer, 6) }
        5 => { _mm256_slli_epi64(kmer, 8) }
        6 => { _mm256_slli_epi64(kmer, 10) }
        7 => { _mm256_slli_epi64(kmer, 12) }
        8 => { _mm256_slli_epi64(kmer, 14) }
        9 => { _mm256_slli_epi64(kmer, 16) }
        10 => { _mm256_slli_epi64(kmer, 18) }
        11 => { _mm256_slli_epi64(kmer, 20) }
        12 => { _mm256_slli_epi64(kmer, 22) }
        13 => { _mm256_slli_epi64(kmer, 24) }
        14 => { _mm256_slli_epi64(kmer, 26) }
        15 => { _mm256_slli_epi64(kmer, 28) }
        16 => { _mm256_slli_epi64(kmer, 30) }
        17 => { _mm256_slli_epi64(kmer, 32) }
        18 => { _mm256_slli_epi64(kmer, 34) }
        19 => { _mm256_slli_epi64(kmer, 36) }
        20 => { _mm256_slli_epi64(kmer, 38) }
        21 => { _mm256_slli_epi64(kmer, 40) }
        22 => { _mm256_slli_epi64(kmer, 42) }
        23 => { _mm256_slli_epi64(kmer, 44) }
        24 => { _mm256_slli_epi64(kmer, 46) }
        25 => { _mm256_slli_epi64(kmer, 48) }
        26 => { _mm256_slli_epi64(kmer, 50) }
        27 => { _mm256_slli_epi64(kmer, 52) }
        28 => { _mm256_slli_epi64(kmer, 54) }
        29 => { _mm256_slli_epi64(kmer, 56) }
        30 => { _mm256_slli_epi64(kmer, 58) }
        31 => { _mm256_slli_epi64(kmer, 60) }
        32 => { _mm256_slli_epi64(kmer, 62) }
        
        //21 => { _mm256_slli_epi64(kmer, 40) }
        //31 => { _mm256_slli_epi64(kmer, 60) }
        _ => { panic!() }
    };
    return shifted;
}

#[target_feature(enable = "avx2")]
pub unsafe fn _shift_mm256_right_by_k(kmer: __m256i, k: usize) -> __m256i {
    // shift right by 2*k - 2
    let shifted = match k {
        1 => { _mm256_srli_epi64(kmer, 0) }
        2 => { _mm256_srli_epi64(kmer, 2) }
        3 => { _mm256_srli_epi64(kmer, 4) }
        4 => { _mm256_srli_epi64(kmer, 6) }
        5 => { _mm256_srli_epi64(kmer, 8) }
        6 => { _mm256_srli_epi64(kmer, 10) }
        7 => { _mm256_srli_epi64(kmer, 12) }
        8 => { _mm256_srli_epi64(kmer, 14) }
        9 => { _mm256_srli_epi64(kmer, 16) }
        10 => { _mm256_srli_epi64(kmer, 18) }
        11 => { _mm256_srli_epi64(kmer, 20) }
        12 => { _mm256_srli_epi64(kmer, 22) }
        13 => { _mm256_srli_epi64(kmer, 24) }
        14 => { _mm256_srli_epi64(kmer, 26) }
        15 => { _mm256_srli_epi64(kmer, 28) }
        16 => { _mm256_srli_epi64(kmer, 30) }
        17 => { _mm256_srli_epi64(kmer, 32) }
        18 => { _mm256_srli_epi64(kmer, 34) }
        19 => { _mm256_srli_epi64(kmer, 36) }
        20 => { _mm256_srli_epi64(kmer, 38) }
        21 => { _mm256_srli_epi64(kmer, 40) }
        22 => { _mm256_srli_epi64(kmer, 42) }
        23 => { _mm256_srli_epi64(kmer, 44) }
        24 => { _mm256_srli_epi64(kmer, 46) }
        25 => { _mm256_srli_epi64(kmer, 48) }
        26 => { _mm256_srli_epi64(kmer, 50) }
        27 => { _mm256_srli_epi64(kmer, 52) }
        28 => { _mm256_srli_epi64(kmer, 54) }
        29 => { _mm256_srli_epi64(kmer, 56) }
        30 => { _mm256_srli_epi64(kmer, 58) }
        31 => { _mm256_srli_epi64(kmer, 60) }
        32 => { _mm256_srli_epi64(kmer, 62) }
        _ => { panic!() }
    };
    return shifted;
}


/**
 * Extract kmers using FracMinHash from a DNA string using AVX2 instructions.
 * The k-mers are extracted from both the forward and reverse strands.
 */
#[target_feature(enable = "avx2")]
pub unsafe fn extract_markers_avx2_masked(string: &[u8], keys_vec: &mut Vec<u64>, values_vec: &mut Vec<u64>, c: usize, k: usize, v: usize, bidirectional: bool) { unsafe {
    let t = k + v;

    if string.len() < t {
        return;
    }

    // divide the string into 4 parts for parallel processing
    let len = (string.len() - t + 1) / 4;
    if len <= 0 {
        return;
    }
    let string1 = &string[0..len + t - 1];
    let string2 = &string[len..2 * len + t - 1];
    let string3 = &string[2 * len..3 * len + t - 1];
    let string4 = &string[3 * len..4 * len + t - 1];

    // storing keys and values
    let mut rolling_key_f = _mm256_set_epi64x(0, 0, 0, 0);
    let mut rolling_value_f = _mm256_set_epi64x(0, 0, 0, 0);
    let mut rolling_key_r = _mm256_set_epi64x(0, 0, 0, 0);
    let mut rolling_value_r = _mm256_set_epi64x(0, 0, 0, 0);

    let rev_sub = _mm256_set_epi64x(3, 3, 3, 3);

    // Initialize keys
    for i in 0..k - 1 {
        let nuc_f1 = BYTE_TO_SEQ[string1[i] as usize] as i64;
        let nuc_f2 = BYTE_TO_SEQ[string2[i] as usize] as i64;
        let nuc_f3 = BYTE_TO_SEQ[string3[i] as usize] as i64;
        let nuc_f4 = BYTE_TO_SEQ[string4[i] as usize] as i64;
        // f_nucs = [nuc_f1, nuc_f2, nuc_f3, nuc_f4]
        let f_nucs = _mm256_set_epi64x(nuc_f4, nuc_f3, nuc_f2, nuc_f1);
        // rolling_key_f = (rolling_key_f << 2)
        rolling_key_f = _mm256_slli_epi64(rolling_key_f, 2);
        // rolling_key_f = rolling_key_f | f_nucs
        rolling_key_f = _mm256_or_si256(rolling_key_f, f_nucs);
    }

    // Initialize values
    for i in 0..v {
        let nuc_f1 = BYTE_TO_SEQ[string1[i + k - 1] as usize] as i64;
        let nuc_f2 = BYTE_TO_SEQ[string2[i + k - 1] as usize] as i64;
        let nuc_f3 = BYTE_TO_SEQ[string3[i + k - 1] as usize] as i64;
        let nuc_f4 = BYTE_TO_SEQ[string4[i + k - 1] as usize] as i64;
        // f_nucs = [nuc_f1, nuc_f2, nuc_f3, nuc_f4]
        let f_nucs = _mm256_set_epi64x(nuc_f4, nuc_f3, nuc_f2, nuc_f1);
        // value_f = (value_f << 2)
        rolling_value_f = _mm256_slli_epi64(rolling_value_f, 2);
        // value_f = value_f | f_nucs
        rolling_value_f = _mm256_or_si256(rolling_value_f, f_nucs);
    }

    if bidirectional {
        // initialize key
        for i in 0..v - 1 {
            let nuc_r1 = 3 - BYTE_TO_SEQ[string1[i] as usize] as i64;
            let nuc_r2 = 3 - BYTE_TO_SEQ[string2[i] as usize] as i64;
            let nuc_r3 = 3 - BYTE_TO_SEQ[string3[i] as usize] as i64;
            let nuc_r4 = 3 - BYTE_TO_SEQ[string4[i] as usize] as i64;
            // r_nucs = [nuc_r1, nuc_r2, nuc_r3, nuc_r4]
            let r_nucs = _mm256_set_epi64x(nuc_r4, nuc_r3, nuc_r2, nuc_r1);
            // rolling_value_r = (rolling_value_r >> 2)
            rolling_value_r = _mm256_srli_epi64(rolling_value_r, 2);
            // rolling_value_r = rolling_value_r | (r_nucs << (2 * (v - 1)))
            let shift_value_r = _shift_mm256_left_by_k(r_nucs, v);
            rolling_value_r = _mm256_or_si256(rolling_value_r, shift_value_r);
        }

        for i in 0..k {
            let nuc_r1 = 3 - BYTE_TO_SEQ[string1[i + v - 1] as usize] as i64;
            let nuc_r2 = 3 - BYTE_TO_SEQ[string2[i + v - 1] as usize] as i64;
            let nuc_r3 = 3 - BYTE_TO_SEQ[string3[i + v - 1] as usize] as i64;
            let nuc_r4 = 3 - BYTE_TO_SEQ[string4[i + v - 1] as usize] as i64;
            // r_nucs = [nuc_r1, nuc_r2, nuc_r3, nuc_r4]
            let r_nucs = _mm256_set_epi64x(nuc_r4, nuc_r3, nuc_r2, nuc_r1);
            // rolling_key_r = (rolling_key_r >> 2)
            rolling_key_r = _mm256_srli_epi64(rolling_key_r, 2);
            // rolling_key_r = rolling_key_r | (r_nucs << (2 * (k - 1)))
            let shift_nuc_r = _shift_mm256_left_by_k(r_nucs, k);
            rolling_key_r = _mm256_or_si256(rolling_key_r, shift_nuc_r);
        }
    }

    let key_mask = ((1u64 << (2 * k)) - 1) as i64;
    let value_mask = ((1u64 << (2 * v)) - 1) as i64;
    let threshold_marker = u64::MAX / c as u64;

    let mm256_key_mask = _mm256_set_epi64x(key_mask, key_mask, key_mask, key_mask);
    let mm256_value_mask = _mm256_set_epi64x(value_mask, value_mask, value_mask, value_mask);

    // iterate over the string
    for i in k + v - 1..(len + t - 1) {
        let nuc_f1 = BYTE_TO_SEQ[string1[i] as usize] as i64;
        let nuc_f2 = BYTE_TO_SEQ[string2[i] as usize] as i64;
        let nuc_f3 = BYTE_TO_SEQ[string3[i] as usize] as i64;
        let nuc_f4 = BYTE_TO_SEQ[string4[i] as usize] as i64;


        let f_nucs = _mm256_set_epi64x(nuc_f4, nuc_f3, nuc_f2, nuc_f1);

        // find the first base of the value
        // first_base_v = (rolling_value_f >> (2 * (v - 1))) & 0b11
        let first_base_v = _mm256_and_si256(_shift_mm256_right_by_k(rolling_value_f, v), rev_sub);

        // f_marker = ((f_marker << 2) | f_nuc) & marker_mask
        rolling_key_f = _mm256_slli_epi64(rolling_key_f, 2);
        rolling_key_f = _mm256_or_si256(rolling_key_f, first_base_v);
        rolling_key_f = _mm256_and_si256(rolling_key_f, mm256_key_mask);

        rolling_value_f = _mm256_slli_epi64(rolling_value_f, 2);
        rolling_value_f = _mm256_or_si256(rolling_value_f, f_nucs);
        rolling_value_f = _mm256_and_si256(rolling_value_f, mm256_value_mask);

        if bidirectional {
            // find the last base of the key
            let last_base_k = _mm256_and_si256(rolling_key_r, rev_sub);
            let r_nucs = _mm256_sub_epi64(rev_sub, f_nucs);

            rolling_value_r = _mm256_srli_epi64(rolling_value_r, 2);
            let shift_value_r = _shift_mm256_left_by_k(last_base_k, v);
            rolling_value_r = _mm256_and_si256(rolling_value_r, mm256_value_mask);
            rolling_value_r = _mm256_or_si256(rolling_value_r, shift_value_r);

            rolling_key_r = _mm256_srli_epi64(rolling_key_r, 2);
            let shift_nuc_r = _shift_mm256_left_by_k(r_nucs, k);
            rolling_key_r = _mm256_and_si256(rolling_key_r, mm256_key_mask);
            rolling_key_r = _mm256_or_si256(rolling_key_r, shift_nuc_r);
        }


        let hash_key_f = mm_hash256_masked(rolling_key_f, None);
        let h1 = _mm256_extract_epi64(hash_key_f, 0) as u64;
        let h2 = _mm256_extract_epi64(hash_key_f, 1) as u64;
        let h3 = _mm256_extract_epi64(hash_key_f, 2) as u64;
        let h4 = _mm256_extract_epi64(hash_key_f, 3) as u64;

        let key1 = _mm256_extract_epi64(rolling_key_f, 0) as u64;
        let key2 = _mm256_extract_epi64(rolling_key_f, 1) as u64;
        let key3 = _mm256_extract_epi64(rolling_key_f, 2) as u64;
        let key4 = _mm256_extract_epi64(rolling_key_f, 3) as u64;

        let value1 = _mm256_extract_epi64(rolling_value_f, 0) as u64;
        let value2 = _mm256_extract_epi64(rolling_value_f, 1) as u64;
        let value3 = _mm256_extract_epi64(rolling_value_f, 2) as u64;
        let value4 = _mm256_extract_epi64(rolling_value_f, 3) as u64;

        if h1 < threshold_marker {
            keys_vec.push(key1 as u64);
            values_vec.push(value1 as u64);
        }
        if h2 < threshold_marker {
            keys_vec.push(key2 as u64);
            values_vec.push(value2 as u64);
        }
        if h3 < threshold_marker {
            keys_vec.push(key3 as u64);
            values_vec.push(value3 as u64);
        }
        if h4 < threshold_marker {
            keys_vec.push(key4 as u64);
            values_vec.push(value4 as u64);
        }

        if bidirectional {
            let hash_key_r = mm_hash256_masked(rolling_key_r, None);
            let hr1 = _mm256_extract_epi64(hash_key_r, 0) as u64;
            let hr2 = _mm256_extract_epi64(hash_key_r, 1) as u64;
            let hr3 = _mm256_extract_epi64(hash_key_r, 2) as u64;
            let hr4 = _mm256_extract_epi64(hash_key_r, 3) as u64;

            let rkey1 = _mm256_extract_epi64(rolling_key_r, 0) as u64;
            let rkey2 = _mm256_extract_epi64(rolling_key_r, 1) as u64;
            let rkey3 = _mm256_extract_epi64(rolling_key_r, 2) as u64;
            let rkey4 = _mm256_extract_epi64(rolling_key_r, 3) as u64;

            let rvalue1 = _mm256_extract_epi64(rolling_value_r, 0) as u64;
            let rvalue2 = _mm256_extract_epi64(rolling_value_r, 1) as u64;
            let rvalue3 = _mm256_extract_epi64(rolling_value_r, 2) as u64;
            let rvalue4 = _mm256_extract_epi64(rolling_value_r, 3) as u64;

            if hr1 < threshold_marker {
                keys_vec.push(rkey1 as u64);
                values_vec.push(rvalue1 as u64);
            }
            if hr2 < threshold_marker {
                keys_vec.push(rkey2 as u64);
                values_vec.push(rvalue2 as u64);
            }
            if hr3 < threshold_marker {
                keys_vec.push(rkey3 as u64);
                values_vec.push(rvalue3 as u64);
            }
            if hr4 < threshold_marker {
                keys_vec.push(rkey4 as u64);
                values_vec.push(rvalue4 as u64);
            }
        }
    }
}}

/**
 * Like `extract_markers_avx2_masked`, but also records the Phred quality scores
 * for each selected value.
 *
 * The string is split into 4 lanes.  Lane `j` (0-indexed) starts at offset
 * `j * lane_len` in the original string/qual slice.  At loop iteration `i`:
 *   - Forward value for lane j covers original positions [j*lane_len + i-v+1 .. j*lane_len + i].
 *   - RC value for lane j covers original positions [j*lane_len + i-k-v+1 .. j*lane_len + i-k]
 *     in *reversed* order (position 0 of RC value → original position j*lane_len + i-k).
 */
#[target_feature(enable = "avx2")]
pub unsafe fn extract_markers_avx2_masked_with_qual(string: &[u8], qual: &[u8], keys_vec: &mut Vec<u64>, values_vec: &mut Vec<u64>, quals_vec: &mut Vec<Vec<u8>>, c: usize, k: usize, v: usize, bidirectional: bool) { unsafe {
    let t = k + v;

    if string.len() < t {
        return;
    }

    let len = (string.len() - t + 1) / 4;
    if len <= 0 {
        return;
    }

    // Lane offsets into the original string/qual.
    let offsets = [0usize, len, 2 * len, 3 * len];

    let string1 = &string[offsets[0]..offsets[0] + len + t - 1];
    let string2 = &string[offsets[1]..offsets[1] + len + t - 1];
    let string3 = &string[offsets[2]..offsets[2] + len + t - 1];
    let string4 = &string[offsets[3]..offsets[3] + len + t - 1];

    let mut rolling_key_f = _mm256_set_epi64x(0, 0, 0, 0);
    let mut rolling_value_f = _mm256_set_epi64x(0, 0, 0, 0);
    let mut rolling_key_r = _mm256_set_epi64x(0, 0, 0, 0);
    let mut rolling_value_r = _mm256_set_epi64x(0, 0, 0, 0);

    let rev_sub = _mm256_set_epi64x(3, 3, 3, 3);

    // Initialize keys
    for i in 0..k - 1 {
        let nuc_f1 = BYTE_TO_SEQ[string1[i] as usize] as i64;
        let nuc_f2 = BYTE_TO_SEQ[string2[i] as usize] as i64;
        let nuc_f3 = BYTE_TO_SEQ[string3[i] as usize] as i64;
        let nuc_f4 = BYTE_TO_SEQ[string4[i] as usize] as i64;
        let f_nucs = _mm256_set_epi64x(nuc_f4, nuc_f3, nuc_f2, nuc_f1);
        rolling_key_f = _mm256_slli_epi64(rolling_key_f, 2);
        rolling_key_f = _mm256_or_si256(rolling_key_f, f_nucs);
    }

    // Initialize values
    for i in 0..v {
        let nuc_f1 = BYTE_TO_SEQ[string1[i + k - 1] as usize] as i64;
        let nuc_f2 = BYTE_TO_SEQ[string2[i + k - 1] as usize] as i64;
        let nuc_f3 = BYTE_TO_SEQ[string3[i + k - 1] as usize] as i64;
        let nuc_f4 = BYTE_TO_SEQ[string4[i + k - 1] as usize] as i64;
        let f_nucs = _mm256_set_epi64x(nuc_f4, nuc_f3, nuc_f2, nuc_f1);
        rolling_value_f = _mm256_slli_epi64(rolling_value_f, 2);
        rolling_value_f = _mm256_or_si256(rolling_value_f, f_nucs);
    }

    if bidirectional {
        for i in 0..v - 1 {
            let nuc_r1 = 3 - BYTE_TO_SEQ[string1[i] as usize] as i64;
            let nuc_r2 = 3 - BYTE_TO_SEQ[string2[i] as usize] as i64;
            let nuc_r3 = 3 - BYTE_TO_SEQ[string3[i] as usize] as i64;
            let nuc_r4 = 3 - BYTE_TO_SEQ[string4[i] as usize] as i64;
            let r_nucs = _mm256_set_epi64x(nuc_r4, nuc_r3, nuc_r2, nuc_r1);
            rolling_value_r = _mm256_srli_epi64(rolling_value_r, 2);
            let shift_value_r = _shift_mm256_left_by_k(r_nucs, v);
            rolling_value_r = _mm256_or_si256(rolling_value_r, shift_value_r);
        }
        for i in 0..k {
            let nuc_r1 = 3 - BYTE_TO_SEQ[string1[i + v - 1] as usize] as i64;
            let nuc_r2 = 3 - BYTE_TO_SEQ[string2[i + v - 1] as usize] as i64;
            let nuc_r3 = 3 - BYTE_TO_SEQ[string3[i + v - 1] as usize] as i64;
            let nuc_r4 = 3 - BYTE_TO_SEQ[string4[i + v - 1] as usize] as i64;
            let r_nucs = _mm256_set_epi64x(nuc_r4, nuc_r3, nuc_r2, nuc_r1);
            rolling_key_r = _mm256_srli_epi64(rolling_key_r, 2);
            let shift_nuc_r = _shift_mm256_left_by_k(r_nucs, k);
            rolling_key_r = _mm256_or_si256(rolling_key_r, shift_nuc_r);
        }
    }

    let key_mask = ((1u64 << (2 * k)) - 1) as i64;
    let value_mask = ((1u64 << (2 * v)) - 1) as i64;
    let threshold_marker = u64::MAX / c as u64;

    let mm256_key_mask = _mm256_set_epi64x(key_mask, key_mask, key_mask, key_mask);
    let mm256_value_mask = _mm256_set_epi64x(value_mask, value_mask, value_mask, value_mask);

    for i in k + v - 1..(len + t - 1) {
        let nuc_f1 = BYTE_TO_SEQ[string1[i] as usize] as i64;
        let nuc_f2 = BYTE_TO_SEQ[string2[i] as usize] as i64;
        let nuc_f3 = BYTE_TO_SEQ[string3[i] as usize] as i64;
        let nuc_f4 = BYTE_TO_SEQ[string4[i] as usize] as i64;

        let f_nucs = _mm256_set_epi64x(nuc_f4, nuc_f3, nuc_f2, nuc_f1);

        let first_base_v = _mm256_and_si256(_shift_mm256_right_by_k(rolling_value_f, v), rev_sub);

        rolling_key_f = _mm256_slli_epi64(rolling_key_f, 2);
        rolling_key_f = _mm256_or_si256(rolling_key_f, first_base_v);
        rolling_key_f = _mm256_and_si256(rolling_key_f, mm256_key_mask);

        rolling_value_f = _mm256_slli_epi64(rolling_value_f, 2);
        rolling_value_f = _mm256_or_si256(rolling_value_f, f_nucs);
        rolling_value_f = _mm256_and_si256(rolling_value_f, mm256_value_mask);

        if bidirectional {
            let last_base_k = _mm256_and_si256(rolling_key_r, rev_sub);
            let r_nucs = _mm256_sub_epi64(rev_sub, f_nucs);

            rolling_value_r = _mm256_srli_epi64(rolling_value_r, 2);
            let shift_value_r = _shift_mm256_left_by_k(last_base_k, v);
            rolling_value_r = _mm256_and_si256(rolling_value_r, mm256_value_mask);
            rolling_value_r = _mm256_or_si256(rolling_value_r, shift_value_r);

            rolling_key_r = _mm256_srli_epi64(rolling_key_r, 2);
            let shift_nuc_r = _shift_mm256_left_by_k(r_nucs, k);
            rolling_key_r = _mm256_and_si256(rolling_key_r, mm256_key_mask);
            rolling_key_r = _mm256_or_si256(rolling_key_r, shift_nuc_r);
        }

        let hash_key_f = mm_hash256_masked(rolling_key_f, None);
        let h1 = _mm256_extract_epi64(hash_key_f, 0) as u64;
        let h2 = _mm256_extract_epi64(hash_key_f, 1) as u64;
        let h3 = _mm256_extract_epi64(hash_key_f, 2) as u64;
        let h4 = _mm256_extract_epi64(hash_key_f, 3) as u64;

        let key1 = _mm256_extract_epi64(rolling_key_f, 0) as u64;
        let key2 = _mm256_extract_epi64(rolling_key_f, 1) as u64;
        let key3 = _mm256_extract_epi64(rolling_key_f, 2) as u64;
        let key4 = _mm256_extract_epi64(rolling_key_f, 3) as u64;

        let value1 = _mm256_extract_epi64(rolling_value_f, 0) as u64;
        let value2 = _mm256_extract_epi64(rolling_value_f, 1) as u64;
        let value3 = _mm256_extract_epi64(rolling_value_f, 2) as u64;
        let value4 = _mm256_extract_epi64(rolling_value_f, 3) as u64;

        // Forward qual: lane j covers original positions [offsets[j]+i-v+1 .. offsets[j]+i].
        if h1 < threshold_marker {
            keys_vec.push(key1);
            values_vec.push(value1);
            quals_vec.push(qual[offsets[0] + i - v + 1..=offsets[0] + i].to_vec());
        }
        if h2 < threshold_marker {
            keys_vec.push(key2);
            values_vec.push(value2);
            quals_vec.push(qual[offsets[1] + i - v + 1..=offsets[1] + i].to_vec());
        }
        if h3 < threshold_marker {
            keys_vec.push(key3);
            values_vec.push(value3);
            quals_vec.push(qual[offsets[2] + i - v + 1..=offsets[2] + i].to_vec());
        }
        if h4 < threshold_marker {
            keys_vec.push(key4);
            values_vec.push(value4);
            quals_vec.push(qual[offsets[3] + i - v + 1..=offsets[3] + i].to_vec());
        }

        if bidirectional {
            let hash_key_r = mm_hash256_masked(rolling_key_r, None);
            let hr1 = _mm256_extract_epi64(hash_key_r, 0) as u64;
            let hr2 = _mm256_extract_epi64(hash_key_r, 1) as u64;
            let hr3 = _mm256_extract_epi64(hash_key_r, 2) as u64;
            let hr4 = _mm256_extract_epi64(hash_key_r, 3) as u64;

            let rkey1 = _mm256_extract_epi64(rolling_key_r, 0) as u64;
            let rkey2 = _mm256_extract_epi64(rolling_key_r, 1) as u64;
            let rkey3 = _mm256_extract_epi64(rolling_key_r, 2) as u64;
            let rkey4 = _mm256_extract_epi64(rolling_key_r, 3) as u64;

            let rvalue1 = _mm256_extract_epi64(rolling_value_r, 0) as u64;
            let rvalue2 = _mm256_extract_epi64(rolling_value_r, 1) as u64;
            let rvalue3 = _mm256_extract_epi64(rolling_value_r, 2) as u64;
            let rvalue4 = _mm256_extract_epi64(rolling_value_r, 3) as u64;

            // RC qual: lane j, RC value position p → original position offsets[j]+i-k-p.
            if hr1 < threshold_marker {
                keys_vec.push(rkey1);
                values_vec.push(rvalue1);
                quals_vec.push((0..v).map(|p| qual[offsets[0] + i - k - p]).collect());
            }
            if hr2 < threshold_marker {
                keys_vec.push(rkey2);
                values_vec.push(rvalue2);
                quals_vec.push((0..v).map(|p| qual[offsets[1] + i - k - p]).collect());
            }
            if hr3 < threshold_marker {
                keys_vec.push(rkey3);
                values_vec.push(rvalue3);
                quals_vec.push((0..v).map(|p| qual[offsets[2] + i - k - p]).collect());
            }
            if hr4 < threshold_marker {
                keys_vec.push(rkey4);
                values_vec.push(rvalue4);
                quals_vec.push((0..v).map(|p| qual[offsets[3] + i - k - p]).collect());
            }
        }
    }
}}