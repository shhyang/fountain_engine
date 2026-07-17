// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

//! Binary matrix storage for inactive-variable coefficients in the master system.
//!
//! Each sparse equation has one row; column `seq` is the coefficient of inactive
//! variable `seq` (0..num_inactive). Rows are bit-packed (`u64` words) for fast
//! GF(2) row XOR during BP peeling and back-substitution.

const WORD_BITS: usize = 64;

/// Row-major GF(2) matrix: `rows[equ_id][seq]` is 0 or 1.
///
/// Rows are allocated with fixed width `max_inactive_num` (same as
/// [`SparseSystem`](super::sparse_equation::SparseSystem) inactive budget).
/// Only columns `0..num_inactive` are logically live during decoding.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BinaryMatrix {
    width: usize,
    words_per_row: usize,
    rows: Vec<u64>,
}

impl BinaryMatrix {
    /// Create an empty matrix with the given maximum inactive column count.
    #[must_use]
    pub fn new(width: usize) -> Self {
        Self {
            width,
            words_per_row: words_for_cols(width),
            rows: Vec::new(),
        }
    }

    /// Maximum number of inactive columns (allocated row width).
    #[must_use]
    pub fn width(&self) -> usize {
        self.width
    }

    /// Number of equation rows stored.
    #[must_use]
    pub fn num_rows(&self) -> usize {
        self.rows.len() / self.words_per_row
    }

    /// `u64` words per allocated row (from max inactive width).
    #[must_use]
    pub fn words_per_row(&self) -> usize {
        self.words_per_row
    }

    /// Packed row `equ_id` (`words_per_row` words; only `0..num_inactive` columns are live).
    #[must_use]
    pub fn row_words(&self, equ_id: usize) -> &[u64] {
        assert!(equ_id < self.num_rows(), "row index out of range");
        let start = equ_id * self.words_per_row;
        &self.rows[start..start + self.words_per_row]
    }

    /// Copy the live prefix of row `equ_id` into `dst` (length ≥ `words_for_cols(num_inactive)`).
    pub fn copy_row_words(&self, equ_id: usize, dst: &mut [u64], num_inactive: usize) {
        assert!(
            num_inactive <= self.width,
            "num_inactive exceeds matrix width"
        );
        let nwords = words_for_cols(num_inactive);
        assert!(dst.len() >= nwords, "dst shorter than row word count");
        dst[..nwords].copy_from_slice(&self.row_words(equ_id)[..nwords]);
    }

    /// Append a row from packed words; returns the new `equ_id`.
    pub fn append_row_from_words(&mut self, src: &[u64], num_inactive: usize) -> usize {
        assert!(
            num_inactive <= self.width,
            "num_inactive exceeds matrix width"
        );
        let nwords = words_for_cols(num_inactive);
        assert!(src.len() >= nwords, "src shorter than row word count");
        let equ_id = self.push_row();
        self.row_words_mut(equ_id)[..nwords].copy_from_slice(&src[..nwords]);
        equ_id
    }

    /// XOR row `src` into row `dst` over columns `0..num_inactive`.
    pub fn xor_row_words(&mut self, dst: usize, src: usize, num_inactive: usize) {
        self.xor_rows(dst, src, num_inactive);
    }

    /// XOR packed `src` into row `dst` over columns `0..num_inactive`.
    pub fn xor_row_from_words(&mut self, dst: usize, src: &[u64], num_inactive: usize) {
        assert!(
            num_inactive <= self.width,
            "num_inactive exceeds matrix width"
        );
        assert!(dst < self.num_rows(), "row index out of range");
        let nwords = words_for_cols(num_inactive);
        assert!(src.len() >= nwords, "src shorter than row word count");
        let full_words = num_inactive / WORD_BITS;
        let rem = num_inactive % WORD_BITS;
        let dst_start = dst * self.words_per_row;
        for (i, &word) in src[..full_words].iter().enumerate() {
            self.rows[dst_start + i] ^= word;
        }
        if rem > 0 {
            let mask = (1_u64 << rem) - 1;
            self.rows[dst_start + full_words] ^= src[full_words] & mask;
        }
    }

    /// Copy row `src_row` from `src` into row `dst_row` of `self` (live prefix `0..num_inactive`).
    pub fn copy_row_from(
        &mut self,
        dst_row: usize,
        src: &BinaryMatrix,
        src_row: usize,
        num_inactive: usize,
    ) {
        assert!(num_inactive <= self.width && num_inactive <= src.width);
        let nwords = words_for_cols(num_inactive);
        self.row_words_mut(dst_row)[..nwords].copy_from_slice(&src.row_words(src_row)[..nwords]);
    }

    /// Append every row of `src` (`0..src.num_rows()`).
    pub fn append_all_rows_from(&mut self, src: &BinaryMatrix, num_inactive: usize) {
        for src_row in 0..src.num_rows() {
            let mut words = vec![0_u64; self.words_per_row];
            src.copy_row_words(src_row, &mut words, num_inactive);
            self.append_row_from_words(&words, num_inactive);
        }
    }

    /// Append all rows of `src` by moving its packed storage (consumes `src`).
    ///
    /// Requires matching allocated row width (`width` / `words_per_row`). Prefer this over
    /// [`append_all_rows_from`](Self::append_all_rows_from) when `src` is no longer needed.
    pub fn append_matrix(&mut self, mut src: BinaryMatrix) {
        assert_eq!(
            self.words_per_row, src.words_per_row,
            "append_matrix: words_per_row mismatch"
        );
        assert_eq!(self.width, src.width, "append_matrix: width mismatch");
        self.rows.append(&mut src.rows);
    }

    /// Expand all rows to dense `Vec<Vec<u8>>` (validation / GF(256) fallback).
    #[must_use]
    pub fn to_dense_rows(&self, num_inactive: usize) -> Vec<Vec<u8>> {
        (0..self.num_rows())
            .map(|r| self.row_bytes(r, num_inactive))
            .collect()
    }

    /// Append a zero row; returns the new equation index (`equ_id`).
    pub fn push_row(&mut self) -> usize {
        let equ_id = self.num_rows();
        self.rows.resize(self.rows.len() + self.words_per_row, 0);
        equ_id
    }

    /// clear row
    pub fn clear_row(&mut self, row_idx: usize) {
        assert!(row_idx < self.num_rows(), "row index out of range");
        let words = self.row_words_mut(row_idx);
        for word in words {
            *word = 0;
        }
    }

    /// Coefficient at `(equ_id, seq)`.
    #[must_use]
    pub fn get(&self, equ_id: usize, seq: usize) -> bool {
        assert!(seq < self.width, "inactive seq out of range");
        let words = self.row_words(equ_id);
        let word = seq / WORD_BITS;
        let bit = seq % WORD_BITS;
        (words[word] >> bit) & 1 == 1
    }

    /// Set coefficient at `(equ_id, seq)` to 1.
    pub fn set(&mut self, equ_id: usize, seq: usize) {
        assert!(seq < self.width, "inactive seq out of range");
        let words = self.row_words_mut(equ_id);
        let word = seq / WORD_BITS;
        let bit = seq % WORD_BITS;
        words[word] |= 1_u64 << bit;
    }

    /// Toggle coefficient at `(equ_id, seq)` in GF(2).
    pub fn flip(&mut self, equ_id: usize, seq: usize) {
        assert!(seq < self.width, "inactive seq out of range");
        let words = self.row_words_mut(equ_id);
        let word = seq / WORD_BITS;
        let bit = seq % WORD_BITS;
        words[word] ^= 1_u64 << bit;
    }

    /// Clear coefficient at `(equ_id, seq)`.
    pub fn clear(&mut self, equ_id: usize, seq: usize) {
        assert!(seq < self.width, "inactive seq out of range");
        let words = self.row_words_mut(equ_id);
        let word = seq / WORD_BITS;
        let bit = seq % WORD_BITS;
        words[word] &= !(1_u64 << bit);
    }

    /// XOR columns `0..num_inactive` of `src` into `dst` (`dst ^= src`).
    pub fn xor_rows(&mut self, dst: usize, src: usize, num_inactive: usize) {
        assert!(
            num_inactive <= self.width,
            "num_inactive exceeds matrix width"
        );
        assert!(
            dst < self.num_rows() && src < self.num_rows(),
            "row index out of range"
        );
        if dst == src || num_inactive == 0 {
            return;
        }
        let words_per_row = self.words_per_row;
        let dst_start = dst * words_per_row;
        let src_start = src * words_per_row;
        let full_words = num_inactive / WORD_BITS;
        let rem = num_inactive % WORD_BITS;

        for i in 0..full_words {
            self.rows[dst_start + i] ^= self.rows[src_start + i];
        }
        if rem > 0 {
            let mask = (1_u64 << rem) - 1;
            self.rows[dst_start + full_words] ^= self.rows[src_start + full_words] & mask;
        }
    }

    /// XOR a dense prefix of `other` (length `num_inactive`) into row `equ_id`.
    pub fn xor_row_bytes(&mut self, equ_id: usize, other: &[u8], num_inactive: usize) {
        assert!(
            num_inactive <= self.width,
            "num_inactive exceeds matrix width"
        );
        assert!(
            other.len() >= num_inactive,
            "other row shorter than num_inactive"
        );
        for (seq, &val) in other[..num_inactive].iter().enumerate() {
            if val != 0 {
                self.flip(equ_id, seq);
            }
        }
    }

    /// Export row `equ_id` as `num_inactive` bytes (0/1), for GE / HDPC helpers.
    #[must_use]
    pub fn row_bytes(&self, equ_id: usize, num_inactive: usize) -> Vec<u8> {
        assert!(
            num_inactive <= self.width,
            "num_inactive exceeds matrix width"
        );
        let mut out = vec![0_u8; num_inactive];
        for seq in self.iter_set_bits(equ_id, num_inactive) {
            out[seq] = 1;
        }
        out
    }

    /// Iterate set column indices `seq` in `0..num_inactive` for row `equ_id`.
    #[must_use]
    pub fn iter_set_bits(&self, equ_id: usize, num_inactive: usize) -> RowSetBitsIter<'_> {
        assert!(
            num_inactive <= self.width,
            "num_inactive exceeds matrix width"
        );
        RowSetBitsIter {
            words: self.row_words(equ_id),
            num_inactive,
            word_idx: 0,
            sub_word_offset: 0,
        }
    }

    /// Packed row `equ_id` as mutable words (`words_per_row`; only `0..num_inactive` columns are live).
    pub fn row_words_mut(&mut self, row_id: usize) -> &mut [u64] {
        assert!(row_id < self.num_rows(), "row index out of range");
        let start = row_id * self.words_per_row;
        let end = start + self.words_per_row;
        &mut self.rows[start..end]
    }

    /// Coefficient at `(row, logical_col)` when columns are permuted by `q`.
    #[must_use]
    pub fn get_logical(&self, row: usize, q: &[usize], logical_col: usize) -> bool {
        self.get(row, q[logical_col])
    }

    /// Set coefficient at `(row, logical_col)` via column permutation `q`.
    pub fn set_logical(&mut self, row: usize, q: &[usize], logical_col: usize, val: bool) {
        if val {
            self.set(row, q[logical_col]);
        } else {
            self.clear(row, q[logical_col]);
        }
    }

    /// Swap two rows in place.
    pub fn swap_rows(&mut self, a: usize, b: usize) {
        if a == b {
            return;
        }
        assert!(a < self.num_rows() && b < self.num_rows());
        let wpr = self.words_per_row;
        let a_start = a * wpr;
        let b_start = b * wpr;
        for i in 0..wpr {
            self.rows.swap(a_start + i, b_start + i);
        }
    }

    /// Keep only the first `n` rows.
    pub fn truncate_rows(&mut self, n: usize) {
        assert!(n <= self.num_rows());
        self.rows.truncate(n * self.words_per_row);
    }

    /// XOR logical columns `start_logical_col..num_inactive` of `src` into `dst`.
    /// `q_inv[physical_col]` is the logical column index (must match current `q`).
    fn xor_rows_logical_upper(
        &mut self,
        dst: usize,
        src: usize,
        q_inv: &[usize],
        start_logical_col: usize,
        num_inactive: usize,
    ) {
        let full_words = num_inactive / WORD_BITS;
        let rem = num_inactive % WORD_BITS;

        for word_idx in 0..full_words {
            let base = word_idx * WORD_BITS;
            let mut src_w = self.row_words(src)[word_idx];
            if src_w == 0 {
                continue;
            }
            let mut apply = 0_u64;
            while src_w != 0 {
                let tz = src_w.trailing_zeros() as usize;
                let phys = base + tz;
                if q_inv[phys] >= start_logical_col {
                    apply |= 1_u64 << tz;
                }
                src_w &= src_w - 1;
            }
            if apply != 0 {
                self.row_words_mut(dst)[word_idx] ^= apply;
            }
        }

        if rem == 0 {
            return;
        }
        let word_idx = full_words;
        let base = full_words * WORD_BITS;
        let rem_mask = (1_u64 << rem) - 1;
        let mut src_w = self.row_words(src)[word_idx] & rem_mask;
        if src_w == 0 {
            return;
        }
        let mut apply = 0_u64;
        while src_w != 0 {
            let tz = src_w.trailing_zeros() as usize;
            let phys = base + tz;
            if q_inv[phys] >= start_logical_col {
                apply |= 1_u64 << tz;
            }
            src_w &= src_w - 1;
        }
        if apply != 0 {
            self.row_words_mut(dst)[word_idx] ^= apply & rem_mask;
        }
    }

    fn build_q_inv(q: &[usize], n: usize) -> Vec<usize> {
        let mut q_inv = vec![0; n];
        for (logical, &phys) in q.iter().enumerate().take(n) {
            q_inv[phys] = logical;
        }
        q_inv
    }

    fn q_inv_swap(q_inv: &mut [usize], pi: usize, pj: usize, i: usize, j: usize) {
        q_inv[pi] = j;
        q_inv[pj] = i;
    }

    /// Incremental LU over GF(2) on packed rows (same semantics as
    /// [`Matrix::lu_decomp_incr`](crate::algebra::linear_algebra::Matrix::lu_decomp_incr) with [`GF2`](crate::algebra::finite_field::GF2)).
    pub fn lu_decomp_incr_gf2(
        &mut self,
        q: &mut [usize],
        r: usize,
        num_inactive: usize,
    ) -> (Vec<usize>, usize) {
        let m = self.num_rows();
        let n = num_inactive;
        let mut p: Vec<usize> = (0..m).collect();
        let mut q_inv = Self::build_q_inv(q, n);

        for i in 0..r {
            for k in r..m {
                let l = self.get_logical(k, q, i);
                //self.set_logical(k, q, i, l); // this is not necessary, because l is got from the same position.
                if l {
                    self.xor_rows_logical_upper(k, i, &q_inv, i + 1, n);
                }
            }
        }

        let mut i = r;
        for j in r..n {
            let mut pivot_found = false;
            for k in i..m {
                if self.get_logical(k, q, j) {
                    p.swap(i, k);
                    self.swap_rows(i, k);
                    pivot_found = true;
                    break;
                }
            }
            if pivot_found {
                let pi = q[i];
                let pj = q[j];
                q.swap(i, j);
                Self::q_inv_swap(&mut q_inv, pi, pj, i, j);
                for k in i + 1..m {
                    let l = self.get_logical(k, q, i);
                    //self.set_logical(k, q, i, l); // this is not necessary, because l is got from the same position.
                    if l {
                        self.xor_rows_logical_upper(k, i, &q_inv, i + 1, n);
                    }
                }
                i += 1;
                if i == m {
                    break;
                }
            }
        }
        (p, i)
    }
}

/// Iterator over set column indices in a packed GF(2) row.
pub struct RowSetBitsIter<'a> {
    words: &'a [u64],
    num_inactive: usize,
    word_idx: usize,
    sub_word_offset: usize,
}

impl Iterator for RowSetBitsIter<'_> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        loop {
            let global_base = self.word_idx * WORD_BITS;
            if global_base >= self.num_inactive {
                return None;
            }
            let bits_in_word = (self.num_inactive - global_base).min(WORD_BITS);
            if self.sub_word_offset >= bits_in_word {
                self.word_idx += 1;
                self.sub_word_offset = 0;
                continue;
            }
            let mut w = self.words[self.word_idx];
            if bits_in_word < WORD_BITS {
                w &= (1_u64 << bits_in_word) - 1;
            }
            w >>= self.sub_word_offset;
            if w != 0 {
                let tz = w.trailing_zeros() as usize;
                let seq = global_base + self.sub_word_offset + tz;
                self.sub_word_offset += tz + 1;
                return Some(seq);
            }
            self.word_idx += 1;
            self.sub_word_offset = 0;
        }
    }
}

fn words_for_cols(cols: usize) -> usize {
    cols.div_ceil(WORD_BITS).max(1)
}

#[cfg(test)]
mod tests {
    use super::BinaryMatrix;

    #[test]
    fn zero_width_matrix_num_rows() {
        let m = BinaryMatrix::new(0);
        assert_eq!(m.num_rows(), 0);
    }

    #[test]
    fn push_row_and_get_set_flip() {
        let mut m = BinaryMatrix::new(10);
        let r0 = m.push_row();
        assert_eq!(r0, 0);
        assert!(!m.get(0, 3));
        m.set(0, 3);
        assert!(m.get(0, 3));
        m.flip(0, 3);
        assert!(!m.get(0, 3));
        m.flip(0, 3);
        assert!(m.get(0, 3));
    }

    #[test]
    fn xor_rows_matches_byte_xor() {
        let mut m = BinaryMatrix::new(100);
        let a = m.push_row();
        let b = m.push_row();
        m.set(a, 1);
        m.set(a, 63);
        m.set(a, 64);
        m.set(b, 1);
        m.set(b, 2);
        m.set(b, 99);

        m.xor_rows(b, a, 100);
        assert!(m.get(a, 1));
        assert!(!m.get(b, 1)); // 1 XOR 1
        assert!(m.get(b, 2));
        assert!(m.get(b, 63));
        assert!(m.get(b, 64));
        assert!(m.get(b, 99));
    }

    #[test]
    fn row_bytes_round_trip() {
        let mut m = BinaryMatrix::new(70);
        let r = m.push_row();
        m.set(r, 0);
        m.set(r, 63);
        m.set(r, 64);
        m.set(r, 69);
        let bytes = m.row_bytes(r, 70);
        assert_eq!(bytes.len(), 70);
        assert_eq!(bytes[0], 1);
        assert_eq!(bytes[63], 1);
        assert_eq!(bytes[64], 1);
        assert_eq!(bytes[69], 1);
        assert_eq!(bytes[1], 0);
    }

    #[test]
    fn iter_set_bits_matches_row_bytes() {
        let mut m = BinaryMatrix::new(100);
        let r = m.push_row();
        for seq in [0, 1, 63, 64, 99] {
            m.set(r, seq);
        }
        let bytes = m.row_bytes(r, 100);
        let bits: Vec<usize> = m.iter_set_bits(r, 100).collect();
        assert_eq!(bits, vec![0, 1, 63, 64, 99]);
        for (seq, &byte) in bytes.iter().enumerate() {
            assert_eq!(byte, if bits.contains(&seq) { 1 } else { 0 });
        }
    }

    #[test]
    fn xor_row_bytes() {
        let mut m = BinaryMatrix::new(8);
        let r = m.push_row();
        m.xor_row_bytes(r, &[1, 0, 1, 0, 0, 0, 0, 1], 8);
        assert!(m.get(r, 0));
        assert!(!m.get(r, 1));
        assert!(m.get(r, 2));
        assert!(m.get(r, 7));
    }

    #[test]
    fn multiple_rows_independent() {
        let mut m = BinaryMatrix::new(5);
        let r0 = m.push_row();
        let r1 = m.push_row();
        m.set(r0, 2);
        m.set(r1, 4);
        assert!(m.get(r0, 2));
        assert!(!m.get(r1, 2));
        assert!(m.get(r1, 4));
        m.xor_rows(r1, r0, 5);
        assert!(m.get(r1, 2));
        assert!(m.get(r1, 4));
    }

    #[test]
    fn row_words_matches_row_bytes() {
        let mut m = BinaryMatrix::new(130);
        let r = m.push_row();
        m.set(r, 0);
        m.set(r, 63);
        m.set(r, 64);
        m.set(r, 129);
        let num_inactive = 130;
        let bytes = m.row_bytes(r, num_inactive);
        let mut copied = vec![0_u64; m.words_per_row()];
        m.copy_row_words(r, &mut copied, num_inactive);
        for (seq, &byte) in bytes.iter().enumerate() {
            let word = seq / 64;
            let bit = seq % 64;
            let bit_set = (copied[word] >> bit) & 1 == 1;
            assert_eq!(bit_set, byte == 1, "seq {seq}");
        }
    }

    #[test]
    fn append_matrix_matches_append_all_rows_from() {
        let n = 20;
        let mut left = BinaryMatrix::new(n);
        let r0 = left.push_row();
        left.set(r0, 3);
        left.set(r0, 17);

        let mut tail = BinaryMatrix::new(n);
        let t0 = tail.push_row();
        tail.set(t0, 1);
        let t1 = tail.push_row();
        tail.set(t1, 19);

        let mut by_ref = left.clone();
        by_ref.append_all_rows_from(&tail, n);

        left.append_matrix(tail);
        assert_eq!(left.num_rows(), by_ref.num_rows());
        assert_eq!(left.row_bytes(0, n), by_ref.row_bytes(0, n));
        assert_eq!(left.row_bytes(1, n), by_ref.row_bytes(1, n));
        assert_eq!(left.row_bytes(2, n), by_ref.row_bytes(2, n));
    }

    #[test]
    fn append_row_from_words_round_trip() {
        let mut m = BinaryMatrix::new(20);
        let r0 = m.push_row();
        m.set(r0, 3);
        m.set(r0, 17);
        let mut words = vec![0_u64; m.words_per_row()];
        m.copy_row_words(r0, &mut words, 20);
        let r1 = m.append_row_from_words(&words, 20);
        assert_eq!(m.row_bytes(r0, 20), m.row_bytes(r1, 20));
    }

    #[test]
    fn xor_row_from_words_matches_xor_rows() {
        let mut m = BinaryMatrix::new(100);
        let a = m.push_row();
        let b = m.push_row();
        m.set(a, 1);
        m.set(a, 63);
        m.set(b, 1);
        m.set(b, 2);
        let mut src = vec![0_u64; m.words_per_row()];
        m.copy_row_words(a, &mut src, 100);
        m.xor_row_from_words(b, &src, 100);
        let mut m2 = BinaryMatrix::new(100);
        let a2 = m2.push_row();
        let b2 = m2.push_row();
        m2.set(a2, 1);
        m2.set(a2, 63);
        m2.set(b2, 1);
        m2.set(b2, 2);
        m2.xor_row_words(b2, a2, 100);
        assert_eq!(m.row_bytes(b, 100), m2.row_bytes(b2, 100));
    }

    #[test]
    fn clear_bit() {
        let mut m = BinaryMatrix::new(4);
        let r = m.push_row();
        m.set(r, 2);
        assert!(m.get(r, 2));
        m.clear(r, 2);
        assert!(!m.get(r, 2));
        assert_eq!(m.width(), 4);
    }

    #[test]
    fn swap_rows_preserves_other_rows() {
        let mut m = BinaryMatrix::new(8);
        let r0 = m.push_row();
        let r1 = m.push_row();
        let r2 = m.push_row();
        m.set(r0, 1);
        m.set(r1, 3);
        m.set(r2, 5);
        m.swap_rows(r0, r2);
        assert!(m.get(r0, 5));
        assert!(m.get(r1, 3));
        assert!(m.get(r2, 1));
    }

    #[test]
    fn lu_decomp_incr_gf2_matches_byte_lu() {
        use crate::algebra::finite_field::GF2;
        use crate::algebra::linear_algebra::Matrix;

        let n = 32;
        let m = 28;
        let mut bytes: Vec<Vec<u8>> = (0..m)
            .map(|i: usize| {
                (0..n)
                    .map(|j: usize| {
                        if i.wrapping_mul(17).wrapping_add(j.wrapping_mul(13)) % 5 != 0 {
                            0
                        } else {
                            1
                        }
                    })
                    .collect()
            })
            .collect();
        let mut packed = BinaryMatrix::new(n);
        for row in &bytes {
            let mut words = vec![0_u64; packed.words_per_row()];
            for (j, &v) in row.iter().enumerate() {
                if v != 0 {
                    words[j / 64] |= 1_u64 << (j % 64);
                }
            }
            packed.append_row_from_words(&words, n);
        }

        let mut q_byte = (0..n).collect::<Vec<_>>();
        let mut q_pack = q_byte.clone();
        let r0 = 0;

        let (p_byte, r_byte) = Matrix::lu_decomp_incr(&GF2::new(), &mut bytes, &mut q_byte, r0);
        let (p_pack, r_pack) = packed.lu_decomp_incr_gf2(&mut q_pack, r0, n);

        assert_eq!(r_byte, r_pack);
        assert_eq!(p_byte, p_pack);
        assert_eq!(q_byte, q_pack);

        bytes.truncate(r_byte);
        packed.truncate_rows(r_pack);
        for (i, byte_row) in bytes.iter().enumerate() {
            assert_eq!(
                *byte_row,
                packed.row_bytes(i, n),
                "row {i} mismatch after LU"
            );
        }
    }

    #[test]
    fn lu_decomp_incr_gf2_incremental_matches_byte() {
        use crate::algebra::finite_field::GF2;
        use crate::algebra::linear_algebra::Matrix;

        let n = 24;
        let mut bytes: Vec<Vec<u8>> = (0..20)
            .map(|i| (0..n).map(|j| u8::from((i + j) % 4 == 0)).collect())
            .collect();
        let mut packed = BinaryMatrix::new(n);
        for row in &bytes {
            let mut words = vec![0_u64; packed.words_per_row()];
            for (j, &v) in row.iter().enumerate() {
                if v != 0 {
                    words[j / 64] |= 1_u64 << (j % 64);
                }
            }
            packed.append_row_from_words(&words, n);
        }

        let mut q_byte = (0..n).collect::<Vec<_>>();
        let mut q_pack = q_byte.clone();

        let (p1, r1) = Matrix::lu_decomp_incr(&GF2::new(), &mut bytes, &mut q_byte, 0);
        let (p2, r2) = packed.lu_decomp_incr_gf2(&mut q_pack, 0, n);
        assert_eq!(r1, r2);
        assert_eq!(p1, p2);

        // append 5 more rows and decompose from r1
        for i in 0..5 {
            let row: Vec<u8> = (0..n).map(|j| u8::from((i + j + 7) % 3 == 0)).collect();
            bytes.push(row.clone());
            let mut words = vec![0_u64; packed.words_per_row()];
            for (j, &v) in row.iter().enumerate() {
                if v != 0 {
                    words[j / 64] |= 1_u64 << (j % 64);
                }
            }
            packed.append_row_from_words(&words, n);
        }

        let (p_byte, r_byte) = Matrix::lu_decomp_incr(&GF2::new(), &mut bytes, &mut q_byte, r1);
        let (p_pack, r_pack) = packed.lu_decomp_incr_gf2(&mut q_pack, r1, n);
        assert_eq!(r_byte, r_pack);
        assert_eq!(p_byte, p_pack);
        assert_eq!(q_byte, q_pack);
    }
}
