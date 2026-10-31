// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE for details.

use crate::algebra::finite_field::GF2;
use crate::algebra::finite_field::GF256;
use crate::algebra::linear_algebra::Matrix;
use crate::core::binary_matrix::BinaryMatrix;
use crate::data_manager::DataManager;
use crate::traits::LDPC;
use crate::types::CodeParams;

/// High-Density Parity-Check (HDPC) precode interface.
///
/// Implementors define how the HDPC matrix is multiplied against data,
/// supporting dense, binary, and sparse representations.
///
/// ## Finite field
///
/// Return [`gf_poly`](HDPC::gf_poly) so the engine can call
/// [`DataManager::config_finite_field`](crate::data_manager::DataManager::config_finite_field)
/// before precoding / inactivation LU:
///
/// | `gf_poly()` | HDPC arithmetic | LU / manager scalars |
/// |-------------|-----------------|----------------------|
/// | [`GF2_FIELD_POLY`](crate::types::GF2_FIELD_POLY) | XOR only (binary matrices) | `GF2`, pivots in `{0,1}` |
/// | `0x11D`, … | GF(256) (`mul_alpha`, table multiply) | Full GF(256) |
///
/// The `gf` argument on [`mul_binary`](HDPC::mul_binary) / [`mul_sparse`](HDPC::mul_sparse) is
/// `manager.gf256()` after configuration; binary implementors may ignore it.
pub trait HDPC {
    /// Multiplies the HDPC matrix by data vectors, writing results into `y` via the data manager.
    ///
    /// `x` and `y` are half-open contiguous data-ID spans: `x` indexes input columns in order,
    /// and `y.len()` must equal `params.h` (output HDPC rows).
    ///
    /// **2.1 breaking change vs 2.0.x:** arguments are [`Range`](std::ops::Range) spans instead of
    /// `&[usize]` slices (IDs must be contiguous).
    fn mul_data(
        &self,
        _manager: &mut DataManager,
        _params: &CodeParams,
        _x: std::ops::Range<usize>,
        _y: std::ops::Range<usize>,
    ) {
    }
    /// Multiplies the HDPC matrix by a binary matrix given as column vectors. Returns the product rows.
    ///
    /// `gf` is the manager's field after [`gf_poly`](HDPC::gf_poly) was applied (`None` for GF(2)-only HDPC).
    #[allow(deprecated)]
    fn mul_binary(
        &self,
        gf: Option<&GF256>,
        params: &CodeParams,
        n: usize,
        v: &dyn Fn(usize) -> Vec<u8>,
    ) -> Vec<Vec<u8>> {
        let s = |i: usize| {
            v(i).iter()
                .enumerate()
                .filter(|&(_, &v)| v == 1)
                .map(|(j, _)| j)
                .collect::<Vec<_>>()
        };
        self.mul_sparse(gf, params, n, &s)
    }

    /// Multiplies the HDPC matrix by a prebuilt \(\tilde{G}\) (`kl` rows of length `num_inactive`).
    ///
    /// Default collects rows via [`mul_binary`](HDPC::mul_binary); GF(256) precodes should override.
    fn mul_binary_from_rows(
        &self,
        gf: Option<&GF256>,
        params: &CodeParams,
        rows: &[Vec<u8>],
    ) -> Vec<Vec<u8>> {
        let n = rows.first().map_or(0, Vec::len);
        let v = |i: usize| rows[i].clone();
        self.mul_binary(gf, params, n, &v)
    }

    /// Multiplies the HDPC matrix by \(\tilde{G}\) rows `0..kl` filled on demand.
    ///
    /// `fill_row(i, dst)` must overwrite `dst[..num_inactive]` with row `i`
    /// (`kl` = [`CodeParams::num_message_ldpc`](crate::types::CodeParams::num_message_ldpc)).
    /// GF(256) precodes should override. The default copies every row and calls
    /// [`mul_binary_from_rows`](HDPC::mul_binary_from_rows); that fallback is not the decode path.
    fn mul_binary_streamed(
        &self,
        gf: Option<&GF256>,
        params: &CodeParams,
        num_inactive: usize,
        fill_row: &mut dyn FnMut(usize, &mut [u8]),
    ) -> Vec<Vec<u8>> {
        let kl = params.num_message_ldpc();
        let mut row = vec![0u8; num_inactive];
        let mut rows = Vec::with_capacity(kl);
        for i in 0..kl {
            fill_row(i, &mut row);
            rows.push(row.clone());
        }
        self.mul_binary_from_rows(gf, params, &rows)
    }

    /// Multiplies the HDPC matrix by packed GF(2) rows of \(\tilde{G}\).
    ///
    /// `fill_words(i, dst)` overwrites the live `u64` words of row `i` (length
    /// `num_inactive.div_ceil(64)`). The product zeroes `dst` before each call.
    /// This is the decode path: rows are weight tens of bits, and the running
    /// accumulator is dense. Default expands the words and calls
    /// [`mul_binary_supports`](Self::mul_binary_supports).
    fn mul_binary_words(
        &self,
        gf: Option<&GF256>,
        params: &CodeParams,
        num_inactive: usize,
        fill_words: &mut dyn FnMut(usize, &mut [u64]),
    ) -> Vec<Vec<u8>> {
        let nwords = num_inactive.div_ceil(64).max(1);
        let mut words = vec![0u64; nwords];
        self.mul_binary_supports(gf, params, num_inactive, &mut |i, cols| {
            words.fill(0);
            fill_words(i, &mut words);
            cols.clear();
            let rem = num_inactive % 64;
            for (wi, &word) in words.iter().enumerate() {
                let mut w = word;
                if wi + 1 == nwords && rem != 0 {
                    w &= (1_u64 << rem) - 1;
                }
                let base = wi * 64;
                while w != 0 {
                    let tz = w.trailing_zeros() as usize;
                    cols.push((base + tz) as u32);
                    w &= w - 1;
                }
            }
        })
    }

    /// Multiplies the HDPC matrix by the supports of \(\tilde{G}\) rows `0..kl`.
    ///
    /// `support_into(i, cols)` must leave `cols` holding the set column indices of row `i`
    /// (values in `0..num_inactive`). The product clears `cols` before each call.
    /// Default expands those columns to bytes and calls [`mul_binary_streamed`](Self::mul_binary_streamed).
    /// GF(256) precodes should override and XOR the supports straight into the accumulator.
    fn mul_binary_supports(
        &self,
        gf: Option<&GF256>,
        params: &CodeParams,
        num_inactive: usize,
        support_into: &mut dyn FnMut(usize, &mut Vec<u32>),
    ) -> Vec<Vec<u8>> {
        let mut cols = Vec::new();
        self.mul_binary_streamed(gf, params, num_inactive, &mut |i, dst| {
            dst[..num_inactive].fill(0);
            cols.clear();
            support_into(i, &mut cols);
            for &col in &cols {
                dst[col as usize] = 1;
            }
        })
    }

    /// Multiplies the HDPC matrix by a packed binary matrix (`tilde_g`: row `i` = column `i`).
    ///
    /// Default falls back to [`mul_binary`](HDPC::mul_binary) via byte rows. Binary precodes
    /// (e.g. Raptor-10) should override with a word-XOR Gray-code implementation.
    fn mul_binary_packed(
        &self,
        gf: Option<&GF256>,
        params: &CodeParams,
        num_inactive: usize,
        tilde_g: &BinaryMatrix,
    ) -> BinaryMatrix {
        let v = |i: usize| tilde_g.row_bytes(i, num_inactive);
        let dense = self.mul_binary(gf, params, num_inactive, &v);
        let mut out = BinaryMatrix::new(num_inactive);
        for row in dense {
            let nwords = out.words_per_row();
            let mut words = vec![0_u64; nwords];
            for (seq, &byte) in row.iter().enumerate().take(num_inactive) {
                if byte != 0 {
                    let word = seq / 64;
                    let bit = seq % 64;
                    words[word] |= 1_u64 << bit;
                }
            }
            out.append_row_from_words(&words, num_inactive);
        }
        out
    }

    /// Multiplies the HDPC matrix by a sparse matrix given as non-zero index lists per column.
    /// This interface is deprecated. Use [`mul_binary`](HDPC::mul_binary) instead.
    #[deprecated(since = "1.1.0", note = "use mul_binary instead")]
    fn mul_sparse(
        &self,
        _gf: Option<&GF256>,
        _params: &CodeParams,
        _n: usize,
        _s: &dyn Fn(usize) -> Vec<usize>,
    ) -> Vec<Vec<u8>> {
        vec![]
    }
    /// Variant of [`mul_sparse`](HDPC::mul_sparse) using the sub-matrix `S_h` dimensions.
    fn mul_sparse_sh(
        &self,
        _gf: Option<&GF256>,
        _params: &CodeParams,
        _s: &dyn Fn(usize) -> Vec<usize>,
    ) -> Vec<Vec<u8>> {
        vec![]
    }

    /// LU of \(I' + D_s S_h\) for ordinary precoding (doc `alg:precode`).
    ///
    /// Returns row permutation `p` (len `h`) and an `h × h` matrix in the same
    /// in-place L/U layout as [`Matrix::lu_decomp`].
    fn lu_idssh(
        &self,
        gf: Option<&GF256>,
        params: &CodeParams,
        ldpc: &dyn LDPC,
    ) -> (Vec<usize>, Vec<Vec<u8>>) {
        let sh_column = |row: usize| {
            let mut cols = Vec::new();
            ldpc.for_each_inactive_in_row(row, &mut |id| {
                if id >= params.b {
                    cols.push(id - params.b);
                }
            });
            cols
        };

        let mut m = self.mul_sparse_sh(gf, params, &sh_column);
        for (i, row) in m.iter_mut().enumerate().take(params.h) {
            row[i] ^= 1; // GF(256) addition is XOR in this codebase
        }

        let (p, r) = match gf {
            Some(gf) => Matrix::lu_decomp(gf, &mut m),
            None => Matrix::lu_decomp(&GF2::new(), &mut m),
        };
        assert_eq!(r, params.h, "I' + D_s S_h singular, rank {r}");

        (p, m)
    }

    /// Primitive polynomial for GF(256) sessions, or [`GF2_FIELD_POLY`](crate::types::GF2_FIELD_POLY) for XOR-only HDPC.
    fn gf_poly(&self) -> u16 {
        0x11D_u16
    }
}
