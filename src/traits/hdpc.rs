// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use crate::algebra::finite_field::GF256;
use crate::data_manager::DataManager;
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

    /// Primitive polynomial for GF(256) sessions, or [`GF2_FIELD_POLY`](crate::types::GF2_FIELD_POLY) for XOR-only HDPC.
    fn gf_poly(&self) -> u16 {
        0x11D_u16
    }
}
