// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use crate::algebra::finite_field::GF2;
use crate::algebra::finite_field::GF256;
use crate::algebra::linear_algebra::Matrix;
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
    /// Multiplies the HDPC matrix by data vectors, writing results into `y_ids` via the data manager.
    fn mul_data(
        &self,
        _manager: &mut DataManager,
        _params: &CodeParams,
        _x_ids: &[usize],
        _y_ids: &[usize],
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
            ldpc.inactive_row(row)
                .into_iter()
                .filter(|&id| id >= params.b)
                .map(|id| id - params.b)
                .collect::<Vec<_>>()
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
