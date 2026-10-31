// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE for details.

use crate::data_manager::DataManager;
use crate::types::CodeParams;

/// LDPC trait for the parity check matrices S_a and S_i.
///
/// Prefer the enumerator APIs (`for_each_active` / `for_each_inactive` and the
/// oriented `*_in_column` / `*_in_row` helpers) so regular sparse patterns can
/// stream nonzeros without allocating adjacency lists. The `active_*` /
/// `inactive_*` methods that return [`Vec`] remain for 2.0.x compatibility;
/// default `for_each_*` implementations bridge through them so older
/// implementors keep working.
pub trait LDPC {
    /// Invoke `f(check_row, var_col)` for each nonzero of \(S_a\).
    ///
    /// Default: column-major walk via [`for_each_active_in_column`](Self::for_each_active_in_column).
    fn for_each_active(&self, params: &CodeParams, f: &mut dyn FnMut(usize, usize)) {
        for var in 0..params.a {
            self.for_each_active_in_column(var, &mut |check| f(check, var));
        }
    }

    /// Invoke `f(check_row, var_col)` for each nonzero of \(S_i\)
    /// (`var_col` in `0..params.num_inactive()`, i.e. message inactive plus HDPC).
    ///
    /// Default: column-major walk via [`for_each_inactive_in_column`](Self::for_each_inactive_in_column).
    fn for_each_inactive(&self, params: &CodeParams, f: &mut dyn FnMut(usize, usize)) {
        for var in 0..params.num_inactive() {
            self.for_each_inactive_in_column(var, &mut |check| f(check, var));
        }
    }

    /// Nonzeros of \(S_a\) in column `var_col`: call `f(check_row)`.
    ///
    /// Default bridges to [`active_column`](Self::active_column) for 2.0.x implementors.
    fn for_each_active_in_column(&self, var_col: usize, f: &mut dyn FnMut(usize)) {
        for check in self.active_column(var_col) {
            f(check);
        }
    }

    /// Nonzeros of \(S_i\) in column `var_col`: call `f(check_row)`.
    ///
    /// Default bridges to [`inactive_column`](Self::inactive_column) for 2.0.x implementors.
    fn for_each_inactive_in_column(&self, var_col: usize, f: &mut dyn FnMut(usize)) {
        for check in self.inactive_column(var_col) {
            f(check);
        }
    }

    /// Nonzeros of \(S_a\) in row `check_row`: call `f(var_col)`.
    ///
    /// Default bridges to [`active_row`](Self::active_row) for 2.0.x implementors.
    fn for_each_active_in_row(&self, check_row: usize, f: &mut dyn FnMut(usize)) {
        for var in self.active_row(check_row) {
            f(var);
        }
    }

    /// Nonzeros of \(S_i\) in row `check_row`: call `f(var_col)`.
    ///
    /// Default bridges to [`inactive_row`](Self::inactive_row) for 2.0.x implementors.
    fn for_each_inactive_in_row(&self, check_row: usize, f: &mut dyn FnMut(usize)) {
        for var in self.inactive_row(check_row) {
            f(var);
        }
    }

    /// Encode the message vectors with the LDPC matrix (default: per-edge XOR into LDPC vars).
    ///
    /// Only message variables participate (`S_a` over `0..a` and \(S_i\) over `0..b`);
    /// HDPC columns of \(S_i\) are applied later in ordinary precoding.
    fn encode(&self, manager: &mut DataManager, params: &CodeParams) {
        manager.ensure_zero_range(manager.data_id_range_of_ldpc_variable());
        self.for_each_active(params, &mut |check, var| {
            let src = manager.data_id_of_active_variable(var);
            let dst = manager.data_id_of_ldpc_variable(check);
            manager.add_one_to_vector(src, dst);
        });
        for var in 0..params.b {
            self.for_each_inactive_in_column(var, &mut |check| {
                let src = manager.data_id_of_inactive_variable(var);
                let dst = manager.data_id_of_ldpc_variable(check);
                manager.add_one_to_vector(src, dst);
            });
        }
    }

    /// Non-zero entries for a column of S_a.
    fn active_column(&self, _var_col: usize) -> Vec<usize> {
        vec![]
    }

    /// Non-zero entries for a row of S_a.
    fn active_row(&self, _check_row: usize) -> Vec<usize> {
        vec![]
    }

    /// Non-zero entries for a column of S_i.
    fn inactive_column(&self, _var_col: usize) -> Vec<usize> {
        vec![]
    }

    /// Non-zero entries for a row of S_i.
    fn inactive_row(&self, _check_row: usize) -> Vec<usize> {
        vec![]
    }
}
