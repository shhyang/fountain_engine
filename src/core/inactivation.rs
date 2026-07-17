// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use crate::algebra::linear_algebra::Matrix;
use crate::core::binary_matrix::BinaryMatrix;
use crate::data_manager::DataManager;
use crate::traits::HDPC;
use crate::types::CodeParams;

use crate::algebra::finite_field::GF2;

/// System for solving inactive vectors in the GE phase (doc eq:inact_dec_1, sec:solv_inactive).
/// Formed from unused coded vectors: each row is the inactive coefficient vector of an unused
/// sparse equation, with RHS given by the equation's data_id.
pub struct InactiveSystem {
    params: CodeParams,
    num_inactive: usize,
    /// GF(2) packed rows until first LU; `None` after materialization or for GF(256).
    matrix_f_packed: Option<BinaryMatrix>,
    /// Dense byte rows for LU (always used after materialize, or exclusively for GF(256)).
    matrix_f: Vec<Vec<u8>>,
    /// RHS: data_id of the coded vector for each row (X_u' in the doc).
    rhs_data_ids: Vec<usize>,
    /// Column permutation of the matrix (used by incremental LU).
    q: Vec<usize>,
    /// Rows `0..num_binary_rows` are master-exported `{0,1}` rows; remainder may be GF(256) HDPC.
    num_binary_rows: usize,
    /// Current rank of the matrix.
    r: usize,
}

impl InactiveSystem {
    /// Create an empty inactive system with a given number of inactive variables (columns).
    pub fn new(params: &CodeParams) -> Self {
        Self {
            params: params.clone(),
            num_inactive: 0,
            matrix_f_packed: None,
            matrix_f: Vec::new(),
            rhs_data_ids: Vec::new(),
            q: Vec::new(),
            num_binary_rows: 0,
            r: 0,
        }
    }

    /// Initialize the system from unused equations.
    ///
    /// Master inactive rows are always accumulated in packed GF(2) form (even on GF(256) precodes).
    pub fn init(&mut self, num_inactive: usize) {
        self.num_inactive = num_inactive;
        self.matrix_f = Vec::new();
        self.rhs_data_ids = Vec::new();
        self.q = (0..num_inactive).collect::<Vec<_>>();
        self.num_binary_rows = 0;
        self.r = 0;
        self.matrix_f_packed = Some(BinaryMatrix::new(num_inactive));
    }

    /// Append a dense byte row (GF(256) path or after packed materialization).
    /// Only called when matrix_f_packet is none.
    fn add_row(&mut self, row: Vec<u8>, rhs_id: usize) {
        //self.materialize_dense_matrix_f();
        self.matrix_f.push(row);
        self.rhs_data_ids.push(rhs_id);
    }

    /// Append a row by copying packed words from the master inactive matrix.
    pub fn add_row_from_master(&mut self, master: &BinaryMatrix, equ_id: usize, rhs_id: usize) {
        if let Some(packed) = &mut self.matrix_f_packed {
            let n = self.num_inactive;
            let mut words = vec![0_u64; packed.words_per_row()];
            master.copy_row_words(equ_id, &mut words, n);
            packed.append_row_from_words(&words, n);
            self.rhs_data_ids.push(rhs_id);
        } else {
            self.add_row(master.row_bytes(equ_id, self.num_inactive), rhs_id);
        }
    }

    fn materialize_packed_prefix(&mut self) {
        if let Some(packed) = self.matrix_f_packed.take() {
            let mut prefix = packed.to_dense_rows(self.num_inactive);
            prefix.append(&mut self.matrix_f);
            self.matrix_f = prefix;
        }
    }

    /*
    fn materialize_dense_matrix_f(&mut self) {
        if let Some(packed) = self.matrix_f_packed.take() {
            let mut dense = packed.to_dense_rows(self.num_inactive);
            self.matrix_f.append(&mut dense);
        }
    }*/

    fn matrix_row_count(&self) -> usize {
        self.matrix_f_packed
            .as_ref()
            .map_or(0, BinaryMatrix::num_rows)
            + self.matrix_f.len()
    }

    /// True when new packed master rows were added after an HDPC dense tail was appended.
    fn has_stale_hdpc_tail(&self) -> bool {
        !self.matrix_f.is_empty()
            && self
                .matrix_f_packed
                .as_ref()
                .is_some_and(|packed| self.r < packed.num_rows())
    }

    /// Drop GF(256) HDPC rows from a prior solve attempt (row set changed).
    fn clear_hdpc_tail(&mut self) {
        let tail = self.matrix_f.len();
        if tail == 0 {
            return;
        }
        let keep = self.rhs_data_ids.len().saturating_sub(tail);
        self.rhs_data_ids.truncate(keep);
        self.matrix_f.clear();
    }

    fn lu_get_coeff(&self, row: usize, logical_col: usize) -> u8 {
        if let Some(packed) = &self.matrix_f_packed {
            if row < packed.num_rows() {
                u8::from(packed.get_logical(row, &self.q, logical_col))
            } else {
                self.matrix_f[row - packed.num_rows()][self.q[logical_col]]
            }
        } else {
            self.matrix_f[row][self.q[logical_col]]
        }
    }

    pub fn append_hdpc_constraints(
        &mut self,
        hdpc: &dyn HDPC,
        tilde_g_rows: &[Vec<u8>],
        manager: &mut DataManager,
    ) {
        let mut hdpc_constraint_data_ids = Vec::with_capacity(self.params.h);
        for _ in 0..self.params.h {
            hdpc_constraint_data_ids.push(manager.temp_data_id());
        }
        let variable_ids_active: Vec<usize> = manager.data_id_range_of_active_variable();
        hdpc.mul_data(
            manager,
            &self.params,
            &variable_ids_active,
            &hdpc_constraint_data_ids,
        );
        self.rhs_data_ids.append(&mut hdpc_constraint_data_ids);

        let mut matrix_d = hdpc.mul_binary_from_rows(manager.gf256(), &self.params, tilde_g_rows);

        let b = self.params.b;
        for (i, row) in matrix_d.iter_mut().enumerate().take(self.params.h) {
            row[i + b] ^= 1;
        }
        self.num_binary_rows = self.matrix_row_count();
        self.matrix_f.append(&mut matrix_d);
    }

    /// GF(2) fast path: packed \(\tilde{G}\) and [`HDPC::mul_binary_packed`].
    pub fn append_hdpc_constraints_packed(
        &mut self,
        hdpc: &dyn HDPC,
        tilde_g: &BinaryMatrix,
        manager: &mut DataManager,
    ) {
        let mut hdpc_constraint_data_ids = Vec::with_capacity(self.params.h);
        for _ in 0..self.params.h {
            hdpc_constraint_data_ids.push(manager.temp_data_id());
        }
        let variable_ids_active: Vec<usize> = manager.data_id_range_of_active_variable();
        hdpc.mul_data(
            manager,
            &self.params,
            &variable_ids_active,
            &hdpc_constraint_data_ids,
        );
        self.rhs_data_ids.append(&mut hdpc_constraint_data_ids);

        let mut matrix_d =
            hdpc.mul_binary_packed(manager.gf256(), &self.params, self.num_inactive, tilde_g);

        let b = self.params.b;
        for i in 0..self.params.h {
            matrix_d.flip(i, i + b);
        }

        if let Some(packed) = &mut self.matrix_f_packed {
            packed.append_all_rows_from(&matrix_d, self.num_inactive);
        } else {
            //self.materialize_dense_matrix_f(); // this is not necessary, because matrix_f_packet == none.
            self.matrix_f
                .extend(matrix_d.to_dense_rows(self.num_inactive));
        }
        self.num_binary_rows = self.matrix_row_count();
    }

    /// Check that the current LU factorization has a non-zero diagonal (full rank).
    pub(crate) fn has_full_rank_factorization(&self) -> bool {
        if self.matrix_row_count() < self.num_inactive {
            return false;
        }
        (0..self.num_inactive).all(|j| self.lu_get_coeff(j, j) != 0)
    }

    /// Perform incremental LU decomposition on F (doc: sec:unused_coded_vecs, `lu_decomp_incr`).
    /// Removes redundant rows and their RHS vectors via `DataManager::remove` when
    /// `remove_redundant_rhs` is true.
    /// Returns the current rank.
    pub fn lu_decomposition(
        &mut self,
        manager: &mut DataManager,
        remove_redundant_rhs: bool,
    ) -> usize {
        if self.has_stale_hdpc_tail() {
            self.clear_hdpc_tail();
        }

        let (p, r) = if let Some(packed) = &mut self.matrix_f_packed {
            if self.matrix_f.is_empty() {
                packed.lu_decomp_incr_gf2(&mut self.q, self.r, self.num_inactive)
            } else {
                // HDPC dense tail: materialize packed master rows, then one mixed LU pass.
                // Keeping packed+dense split through pass-2 LU leaves the packed prefix out of
                // sync with column swaps / dense eliminations, so `lu_solve_dense` reads wrong
                // coefficients (see `test_hdpc_precode_1` in fountain_scheme).
                self.materialize_packed_prefix();
                Matrix::lu_decomp_incr_mixed(
                    manager.gf256().unwrap(),
                    &mut self.matrix_f,
                    &mut self.q,
                    self.r,
                    self.num_binary_rows,
                )
            }
        } else if manager.gf256().is_some() {
            Matrix::lu_decomp_incr_mixed(
                manager.gf256().unwrap(),
                &mut self.matrix_f,
                &mut self.q,
                self.r,
                self.num_binary_rows,
            )
        } else {
            // matrix_f_packet is none, gf256() is none
            //self.materialize_dense_matrix_f(); // this is not necessary, because matrix_f is already dense
            Matrix::lu_decomp_incr(&GF2::new(), &mut self.matrix_f, &mut self.q, self.r)
        };

        // remove redundancy and apply row permutation on rhs_data_ids
        let ids_new: Vec<_> = p.iter().take(r).map(|&pi| self.rhs_data_ids[pi]).collect();
        for &pi in p.iter().skip(r) {
            if remove_redundant_rhs {
                manager.remove(self.rhs_data_ids[pi]);
            }
        }
        self.rhs_data_ids = ids_new;
        if let Some(packed) = &mut self.matrix_f_packed {
            let packed_keep = r.min(packed.num_rows());
            packed.truncate_rows(packed_keep);
            self.matrix_f.truncate(r.saturating_sub(packed_keep));
        } else {
            self.matrix_f.truncate(r);
        }
        if self.num_binary_rows > r {
            self.num_binary_rows = r;
        }
        self.r = r;
        r
    }

    /// True when inactive rows are only packed GF(2) master exports (no GF(256) HDPC tail).
    pub(crate) fn is_binary_only(&self) -> bool {
        self.matrix_f.is_empty() && self.matrix_f_packed.is_some()
    }

    /// Only call when rank equals `num_inactive`. Forward + backward substitution to solve
    /// F B_i = X_u'. Returns data_ids for each inactive variable (index = inactive sequence).
    pub fn lu_solve(&mut self, manager: &mut DataManager) -> Vec<usize> {
        // G4.2: GF(256) precodes still use packed GF(2) LU/solve when HDPC rows were not needed.
        if self.is_binary_only() {
            self.lu_solve_gf2_packed(manager)
        } else {
            if self.matrix_f_packed.is_some() {
                self.materialize_packed_prefix();
            }
            self.lu_solve_dense(manager)
        }
    }

    fn lu_solve_gf2_packed(&mut self, manager: &mut DataManager) -> Vec<usize> {
        let q = self.q.clone();
        let n = self.rhs_data_ids.len();
        if self
            .matrix_f_packed
            .as_ref()
            .is_none_or(|p| p.num_rows() != n)
        {
            panic!("The number of rows in F must equal the number of RHS vectors");
        }

        for (j, &phys_col) in q.iter().enumerate().take(n - 1) {
            for i in j + 1..n {
                let add = self.matrix_f_packed.as_ref().unwrap().get(i, phys_col);
                if add {
                    let src = self.rhs_data_ids[j];
                    let dst = self.rhs_data_ids[i];
                    manager.add_one_to_vector(src, dst);
                }
            }
        }

        for j in (0..n).rev() {
            let phys_col = q[j];
            if !self.matrix_f_packed.as_ref().unwrap().get(j, phys_col) {
                panic!(
                    "Singular matrix in inactive system: diagonal element at position {} is zero",
                    j
                );
            }
            for i in (0..j).rev() {
                if self.matrix_f_packed.as_ref().unwrap().get(i, phys_col) {
                    let src = self.rhs_data_ids[j];
                    let dst = self.rhs_data_ids[i];
                    manager.add_one_to_vector(src, dst);
                }
            }
        }

        let mut ids_new = vec![0; self.r];
        for i in 0..self.r {
            ids_new[q[i]] = self.rhs_data_ids[i];
        }
        ids_new
    }

    fn lu_solve_dense(&mut self, manager: &mut DataManager) -> Vec<usize> {
        let n = self.matrix_row_count();
        if n != self.rhs_data_ids.len() {
            panic!("The number of rows in F must equal the number of RHS vectors");
        }

        // Forward substitution (L part): solve Ly = b
        for j in 0..n - 1 {
            for i in j + 1..n {
                let l = self.lu_get_coeff(i, j);
                if l == 0 {
                    continue;
                }
                if l == 1 {
                    manager.add_one_to_vector(self.rhs_data_ids[j], self.rhs_data_ids[i]);
                } else {
                    manager.mul_add(self.rhs_data_ids[j], l, self.rhs_data_ids[i]);
                }
            }
        }

        // Backward substitution (U part): solve Ux = y
        for j in (0..n).rev() {
            let l = self.lu_get_coeff(j, j);
            if l == 0 {
                panic!(
                    "Singular matrix in inactive system: diagonal element at position {} is zero",
                    j
                );
            }

            if l != 1 {
                manager.divide_scalar(l, self.rhs_data_ids[j]);
            }

            for i in (0..j).rev() {
                let l = self.lu_get_coeff(i, j);
                if l == 0 {
                    continue;
                }
                if l == 1 {
                    manager.add_one_to_vector(self.rhs_data_ids[j], self.rhs_data_ids[i]);
                } else {
                    manager.mul_add(self.rhs_data_ids[j], l, self.rhs_data_ids[i]);
                }
            }
        }

        // Apply column permutation on rhs_data_ids
        let mut ids_new = vec![0; self.r];
        for i in 0..self.r {
            ids_new[self.q[i]] = self.rhs_data_ids[i];
        }
        ids_new
    }
}
