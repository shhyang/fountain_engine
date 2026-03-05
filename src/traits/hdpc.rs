// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use crate::data_manager::DataManager;
//use crate::algebra::finite_field::GF256;
use crate::types::CodeParams;

/// High-Density Parity-Check (HDPC) precode interface.
///
/// Implementors define how the HDPC matrix is multiplied against data,
/// supporting dense, binary, and sparse representations.
pub trait HDPC {
    /// Multiplies the HDPC matrix by data vectors, writing results into `y_ids` via the data manager.
    fn mul_data(&self, _manager: &mut DataManager, _params: &CodeParams, _x_ids: &[usize], _y_ids: &[usize]) {
    }
    /// Multiplies the HDPC matrix by a binary matrix given as column vectors. Returns the product rows.
    fn mul_binary(&self, params: &CodeParams, n: usize, v: &dyn Fn(usize) -> Vec<u8>) -> Vec<Vec<u8>> {
        let s = |i: usize| v(i)
            .iter()
            .enumerate()
            .filter(|&(_, &v)| v == 1)
            .map(|(j, _)| j).collect::<Vec<_>>();
        self.mul_sparse(params, n, &s)
    }
    /// Multiplies the HDPC matrix by a sparse matrix given as non-zero index lists per column.
    fn mul_sparse(&self, _params: &CodeParams, _n: usize, _s: &dyn Fn(usize) -> Vec<usize>) -> Vec<Vec<u8>> {
        vec![]
    }
    /// Variant of [`mul_sparse`](HDPC::mul_sparse) using the sub-matrix `S_h` dimensions.
    fn mul_sparse_sh(&self, _params: &CodeParams, _s: &dyn Fn(usize) -> Vec<usize>) -> Vec<Vec<u8>> {
        vec![]
    }
}
