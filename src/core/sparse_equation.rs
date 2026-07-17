// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

//! Master system implementation for fountain code encoding and decoding.
//!
//! Design follows docs/doc-erasurecodes.org, section "Representation of the System"
//! (anchor: `sec:representation_of_the_system`).
//!
//! The master system is the linear system (eq:sysdec) with variable vectors
//! (B_a, L, B_b, P) and three block rows: LT encoding, LDPC constraints, HDPC constraints.
//! This module implements the representation: variable IDs sorted by status (decoded,
//! active, inactive), sparse equations (data_id, decoded/active coefficients),
//! inactive coefficients in [`BinaryMatrix`], and BP phase operations.

#![allow(dead_code)]

use super::binary_matrix::BinaryMatrix;
use super::tanner_graph::TannerGraphDegree2;
use crate::data_manager::DataManager;
use crate::traits::LDPC;
use crate::types::{CodeParams, InactivationStrategy};
use smallvec::SmallVec;
use std::collections::VecDeque;

/// Inline capacity for LT/LDPC active variable lists (typical degree is small).
type ActiveVarIds = SmallVec<[usize; 8]>;
/// Per-variable adjacency to sparse equations (inline for typical low degree).
type VarAdjacentList = SmallVec<[usize; 4]>;

// ---------- Sparse equation (one row of the first two block rows) ----------

/// Status of a sparse equation in the master system.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum EquationStatus {
    /// Equation has been used to solve a variable vector.
    Decoded,
    /// Not yet used.
    Unused,
}

/// Status of a variable vector in the master system.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum VariableStatus {
    /// Not yet decoded; may be resolved by BP or inactivated.
    Active,
    /// Decoded by BP decoding.
    Decoded { equ_id: usize },
    /// Inactivated (pre-inactive or dynamically inactivated).
    Inactive { seq: usize },
}

/// One sparse equation: corresponds to an LDPC constraint or an LT-coded vector.
/// Stored as: data_id (RHS) and active variable IDs. Decoded/inactive coeffs live in
/// [`SparseSystem::inactive_coeff_matrix`] (row = `equ_id`).
#[derive(Clone, Debug)]
pub struct SparseEquation {
    /// Data vector ID of the coded vector.
    data_id: usize,
    /// Variable IDs that are still active in this equation.
    active_var_ids: ActiveVarIds,
}

impl SparseEquation {
    /// Degree = number of active variable vectors in this equation.
    fn degree(&self) -> usize {
        self.active_var_ids.len()
    }

    fn decode_variable(&mut self, var_id: usize) {
        if let Some(pos) = self.active_var_ids.iter().position(|&x| x == var_id) {
            self.active_var_ids.swap_remove(pos);
        } else {
            eprintln!(
                "Unexpected: variable {} is not active in equation {}",
                var_id, self.data_id
            );
        }
    }

    fn inactivate_variable(&mut self, var_id: usize) {
        if let Some(pos) = self.active_var_ids.iter().position(|&x| x == var_id) {
            self.active_var_ids.swap_remove(pos);
        } else {
            eprintln!(
                "Unexpected: variable {} is not active in equation {}",
                var_id, self.data_id
            );
        }
    }
}

struct IDPosition {
    ids: Vec<usize>,
    pos: Vec<usize>,
}

impl IDPosition {
    pub fn new(len: usize) -> Self {
        let ids: Vec<usize> = (0..len).collect();
        let pos: Vec<usize> = (0..len).collect();
        Self { ids, pos }
    }

    pub fn push_next(&mut self) {
        let id = self.ids.len();
        self.ids.push(id);
        self.pos.push(id);
    }

    pub fn swap(&mut self, i: usize, j: usize) {
        let vi = self.ids[i];
        let vj = self.ids[j];
        self.ids[i] = vj;
        self.ids[j] = vi;
        self.pos[vi] = j;
        self.pos[vj] = i;
    }
}

// ---------- Master system ----------

/// Master system state as in doc sec:representation_of_the_system.
///
/// - `var_ids`: variable vector IDs; layout is
///   [decoded (0..num_decoded)] [active] [inactive (num_total - num_inactive ..)].
/// - `var_adjacent[var_id]`: indices into `sparse_equations` for equations involving that variable.
/// - `equation_ids`: equation indices; [decoded (0..num_decoded)] [unused].
/// - `inactive_coeff_matrix`: GF(2) coefficients for inactive variables, one row per equation.
/// - `solvable_equations`: queue of equation indices with degree 1.
pub struct SparseSystem {
    params: CodeParams,
    /// Number of variables already decoded.
    num_decoded: usize,
    /// Number of variables originally inactivated.
    num_original_inactive: usize,
    /// Number of inactive variables (pre-inactive + dynamically inactivated).
    num_inactive: usize,
    /// Max allowed inactive count (cap for dynamic inactivation).
    max_inactive_num: usize,

    /// Inactivation strategy
    inac_strategy: InactivationStrategy,

    /// Variable IDs reordered by status: decoded, then active, then inactive.
    var_ids: IDPosition,
    /// For each variable id (0..num_total), list of sparse equation indices adjacent to it.
    var_adjacent: Vec<VarAdjacentList>,

    /// Inactive-variable coefficient matrix (row index = equation id).
    inactive_coeff_matrix: BinaryMatrix,
    /// All sparse equations (LDPC + LT coded).
    sparse_equations: Vec<SparseEquation>,
    /// Equation indices reordered by status: decoded first, then unused.
    equation_ids: IDPosition,
    /// Indices of equations with degree 1 (solvable).
    solvable_equations: VecDeque<usize>,
    /// Reused target data IDs for `bp_solve` / `reverse_bp` / `forward_bp`.
    scratch_targets: Vec<usize>,
    /// Precomputed decoded-equation peers for Original back-substitution (`build_bs_peers`).
    bs_peers: Option<Vec<VarAdjacentList>>,
}

// ---------- BP phase ----------
impl SparseSystem {
    /// Initialize variable parameters as in doc alg:var_list.
    pub fn new(
        params: &CodeParams,
        max_inactive_num: usize,
        inac_strategy: InactivationStrategy,
    ) -> Self {
        let num_total = params.num_total();
        let num_decoded = 0;
        let num_inactive = params.num_inactive(); //do not use params.num_pre_inactive(); 

        // Variable adjacent list for each variable (length num_total).
        let var_adjacent = (0..num_total).map(|_| VarAdjacentList::new()).collect();

        Self {
            params: params.clone(),
            num_decoded,
            num_original_inactive: num_inactive,
            num_inactive,
            max_inactive_num,
            inac_strategy,
            var_ids: IDPosition::new(num_total),
            var_adjacent,
            inactive_coeff_matrix: BinaryMatrix::new(max_inactive_num),
            sparse_equations: Vec::new(),
            equation_ids: IDPosition::new(0),
            solvable_equations: VecDeque::new(),
            scratch_targets: Vec::with_capacity(16),
            bs_peers: None,
        }
    }

    /// Column index in the inactive matrix.
    ///
    /// GE uses dynamic columns first (newest at 0), then pre-inactive — same layout as
    /// legacy `inactive_adj_matrix_row`, not the internal matrix column order.
    fn inactive_matrix_seq(&self, var_id_seq: usize) -> usize {
        let n = self.num_inactive;
        let m = n - self.num_original_inactive;
        if var_id_seq < m {
            n - var_id_seq - 1
        } else {
            var_id_seq - m
        }
    }

    /// Export an equation's inactive coefficients in GE / LU column order.
    fn inactive_coeff_row_ge(&self, equ_id: usize) -> Vec<u8> {
        self.inactive_coeff_matrix
            .row_bytes(equ_id, self.num_inactive)
        //let mut ge = vec![0u8; n];
        //for (internal, &c) in raw.iter().enumerate() {
        //    if c != 0 {
        //        ge[self.inactive_ge_col(internal)] = 1;
        //    }
        //}
        //ge
    }

    /// Get the status of a variable vector.
    fn var_status(&self, var_id: usize) -> VariableStatus {
        let pos = self.var_ids.pos[var_id];
        let total = self.params.num_total();
        if pos < self.num_decoded {
            VariableStatus::Decoded {
                equ_id: self.equation_ids.ids[pos],
            }
        } else if pos >= total - self.num_inactive {
            VariableStatus::Inactive {
                seq: self.inactive_matrix_seq(pos - (total - self.num_inactive)),
            }
        } else {
            VariableStatus::Active
        }
    }

    /// Add LDPC constraint to the system (doc alg:add_ldpc_constraint).
    fn add_ldpc_constraint_by_id(
        &mut self,
        manager: &mut DataManager,
        ldpc: &dyn LDPC,
        ldpc_id: usize,
    ) {
        let data_id = manager.coded_data_id(self.params.k + ldpc_id);
        let mut active_var_ids = ldpc.active_row(ldpc_id);
        let inactive_var_ids = ldpc.inactive_row(ldpc_id);

        active_var_ids.push(self.params.a + ldpc_id);

        let next_equation_id = self.sparse_equations.len();
        self.inactive_coeff_matrix.push_row();
        for &idx in &inactive_var_ids {
            self.inactive_coeff_matrix.set(next_equation_id, idx);
        }

        for &var_id in &active_var_ids {
            self.var_adjacent[var_id].push(next_equation_id);
        }

        self.sparse_equations.push(SparseEquation {
            data_id,
            active_var_ids: active_var_ids.into(),
        });
        self.equation_ids.push_next();

        if self.sparse_equations.last().unwrap().degree() == 1 {
            self.solvable_equations.push_back(next_equation_id);
        }
    }

    /// Add all LDPC constraints (convenience).
    pub fn add_ldpc_constraints(&mut self, manager: &mut DataManager, ldpc: &dyn LDPC) {
        for ldpc_id in 0..self.params.l {
            self.add_ldpc_constraint_by_id(manager, ldpc, ldpc_id);
        }
    }

    /// Add one LT-coded vector to the system (doc alg:add_coded_vector).
    /// `coded_id` is in [k+l+h, ...); `data_id` is the data vector for this coded vector.
    /// `active_ids` and `inactive_ids` are the degree set for this coded_id (variable indices).
    pub fn add_lt_coded_vector(
        &mut self,
        manager: &mut DataManager,
        data_id: usize,
        degree_set: &[usize],
    ) -> usize {
        let next_equation_id = self.sparse_equations.len();
        self.inactive_coeff_matrix.push_row();

        let mut active_var_ids: ActiveVarIds = SmallVec::new();
        let num_inactive = self.num_inactive;

        for &var_id in degree_set {
            let status = self.var_status(var_id);
            match status {
                VariableStatus::Decoded { equ_id } => {
                    let eq = &self.sparse_equations[equ_id];
                    manager.add_one_to_vector(eq.data_id, data_id);
                    self.inactive_coeff_matrix
                        .xor_rows(next_equation_id, equ_id, num_inactive);
                }
                VariableStatus::Active => {
                    active_var_ids.push(var_id);
                    self.var_adjacent[var_id].push(next_equation_id);
                }
                VariableStatus::Inactive { seq } => {
                    self.inactive_coeff_matrix.flip(next_equation_id, seq);
                }
            }
        }

        self.sparse_equations.push(SparseEquation {
            data_id,
            active_var_ids,
        });
        self.equation_ids.push_next();

        if self.sparse_equations.last().unwrap().degree() == 1 {
            self.solvable_equations.push_back(next_equation_id);
        }
        next_equation_id
    }

    /// Solve one degree-1 sparse equation (doc alg:bpdecoding BP_solve).
    /// `equ_id` is the index into `sparse_equations`. Call only when equation has degree 1 before. But it may become degree 0 after back substitution.
    fn bp_solve(&mut self, manager: &mut DataManager, equ_id: usize) {
        let Some(var_id) = self.sparse_equations[equ_id].active_var_ids.pop() else {
            return;
        };
        // Move variable into decoded region.
        let var_pos = self.var_ids.pos[var_id];
        self.var_ids.swap(var_pos, self.num_decoded);

        // Move equation into decoded region.
        let equ_pos = self.equation_ids.pos[equ_id];
        self.equation_ids.swap(equ_pos, self.num_decoded);

        self.num_decoded += 1;

        let var_data_id = manager.data_id_of_variable_vector(var_id);
        manager.move_to(self.sparse_equations[equ_id].data_id, var_data_id);
        self.sparse_equations[equ_id].data_id = var_data_id;

        let num_inactive = self.num_inactive;
        //let mut adj_data_ids = Vec::new();
        self.scratch_targets.clear();

        for &equ_id_adj in &self.var_adjacent[var_id] {
            if equ_id_adj == equ_id {
                continue;
            }
            self.inactive_coeff_matrix
                .xor_rows(equ_id_adj, equ_id, num_inactive);
            let eq_adj = &mut self.sparse_equations[equ_id_adj];
            //adj_data_ids.push(eq_adj.data_id);
            self.scratch_targets.push(eq_adj.data_id);
            eq_adj.decode_variable(var_id);
            if eq_adj.degree() == 1 {
                self.solvable_equations.push_back(equ_id_adj);
            }
        }
        //manager.broadcast_add_owned(var_data_id, adj_data_ids);
        self.broadcast_src_to_scratch_targets(manager, var_data_id);
    }

    /// Inactivate a variable (doc alg:inactivation).
    fn inactivate(&mut self, var_id: usize) {
        let pos = self.var_ids.pos[var_id];
        self.var_ids
            .swap(pos, self.params.num_total() - self.num_inactive - 1);
        self.num_inactive += 1;
        // pay attention to the sequence number of the inactive variable, and the new position of the variable in the var_ids vector.
        let new_seq = self.num_inactive - 1;

        for &equ_id in &self.var_adjacent[var_id] {
            let eq = &mut self.sparse_equations[equ_id];
            eq.inactivate_variable(var_id);
            self.inactive_coeff_matrix.flip(equ_id, new_seq);

            if eq.degree() == 1 {
                self.solvable_equations.push_back(equ_id);
            }
        }
    }

    /// Pop a solvable equation index, if any.
    fn next_solvable(&mut self) -> Option<usize> {
        self.solvable_equations.pop_front()
    }

    /// BP decoding with inactivation (doc alg:bpdecoding BP_inactivation).
    pub fn run_bp_inactivation(&mut self, manager: &mut DataManager) -> usize {
        loop {
            if let Some(equ_id) = self.next_solvable() {
                self.bp_solve(manager, equ_id);
            } else if let Some(var_id) = self.next_inactivation() {
                self.inactivate(var_id);
            } else {
                break;
            }
        }
        self.num_decoded
    }

    /// Whether all variable vectors are either decoded or inactive (BP phase complete).
    pub fn is_bp_complete(&self) -> bool {
        self.num_decoded + self.num_inactive == self.params.num_total()
    }

    /// Current number of inactive variables.
    pub fn num_inactive(&self) -> usize {
        self.num_inactive
    }
}

// ---------- Inactive system (GE phase) ----------
impl SparseSystem {
    /// Get the inactivation coefficient row of a sparse equation (GE column order).
    pub fn inactive_coeff_row(&self, equ_id: usize) -> Vec<u8> {
        self.inactive_coeff_row_ge(equ_id)
    }

    /// Master inactive coefficient storage (row = sparse `equ_id`).
    #[must_use]
    pub fn inactive_coeff_matrix(&self) -> &BinaryMatrix {
        &self.inactive_coeff_matrix
    }

    /// Generate the system of inactive variables from unused equations (doc: phase_change).
    /// Packed GF(2) path: copies `u64` row words from the master matrix (no `row_bytes`).
    pub fn inactive_system_from_unused_equations_packed(
        &self,
        add_row: &mut dyn FnMut(&BinaryMatrix, usize, usize),
    ) {
        for &equ_id in &self.equation_ids.ids[self.num_decoded..] {
            let data_id = self.sparse_equations[equ_id].data_id;
            add_row(&self.inactive_coeff_matrix, equ_id, data_id);
        }
    }

    /// Generate the system of inactive variables from unused equations (doc: phase_change,
    /// "forming equations of unused coded vectors"). Collects unused sparse equations and forms
    /// F B_i = X_u' with matrix_f and rhs_data_ids.
    /// Call when BP is complete (is_bp_complete()).
    pub fn inactive_system_from_unused_equations(&self, add_row: &mut dyn FnMut(Vec<u8>, usize)) {
        for &equ_id in &self.equation_ids.ids[self.num_decoded..] {
            let eq = &self.sparse_equations[equ_id];
            let row = self.inactive_coeff_row_ge(equ_id);
            add_row(row, eq.data_id);
        }
    }

    /// used by the new inactive_sover
    pub fn inactive_coefficients_matrix_from_unused_equations(&self) -> (BinaryMatrix, Vec<usize>) {
        let mut matrix = BinaryMatrix::new(self.num_inactive);
        let mut rhs_ids = Vec::new();
        for &equ_id in &self.equation_ids.ids[self.num_decoded..] {
            let data_id = self.sparse_equations[equ_id].data_id;
            matrix.append_row_from_words(
                self.inactive_coeff_matrix.row_words(equ_id),
                self.num_inactive,
            );
            rhs_ids.push(data_id);
        }
        (matrix, rhs_ids)
    }

    /// Referenced from doc eq:tildeD, sec:solv_inactive.
    /// Return a binary vector of length num_inactive for each variable vector from 0 to k+l-1,
    /// If the variable vector is decoded, return the inactive coefficient vector of the coded vector used to decode it.
    /// If the variable vector is inactive, return a vector with all zero entries except for the index of the inactive variable among all inactive variables.
    pub fn tilde_g_row(&self, var_id: usize) -> Vec<u8> {
        let status = self.var_status(var_id);
        match status {
            VariableStatus::Decoded { equ_id } => self.inactive_coeff_row_ge(equ_id),
            VariableStatus::Inactive { seq } => {
                let mut row = vec![0u8; self.num_inactive];
                row[seq] = 1;
                row
            }
            VariableStatus::Active => {
                unreachable!()
            }
        }
    }

    /// Materialize all `kl` rows of \(\tilde{G}\) in one pass (doc eq:tildeD).
    ///
    /// Row `i` matches [`tilde_g_row`](Self::tilde_g_row)(`i`). Used by the GF(256) HDPC path
    /// so `rq_mul_binary` can index rows without a per-index callback.
    pub fn tilde_g_rows(&self) -> Vec<Vec<u8>> {
        let kl = self.params.num_message_ldpc();
        let n = self.num_inactive;
        let mut rows = Vec::with_capacity(kl);
        for var_id in 0..kl {
            match self.var_status(var_id) {
                VariableStatus::Decoded { equ_id } => {
                    rows.push(self.inactive_coeff_row_ge(equ_id));
                }
                VariableStatus::Inactive { seq } => {
                    let mut row = vec![0u8; n];
                    row[seq] = 1;
                    rows.push(row);
                }
                VariableStatus::Active => {
                    unreachable!("tilde_g_rows requires BP complete (no active vars)")
                }
            }
        }
        rows
    }

    /// Materialize \(\tilde{G}\) as a packed `kl × num_inactive` matrix (doc eq:tildeD).
    ///
    /// Row `i` matches [`tilde_g_row`](Self::tilde_g_row)(`i`) without byte expansion.
    pub fn tilde_g_packed(&self) -> BinaryMatrix {
        let kl = self.params.num_message_ldpc();
        let n = self.num_inactive;
        let mut g = BinaryMatrix::new(n);
        for var_id in 0..kl {
            let row = g.push_row();
            match self.var_status(var_id) {
                VariableStatus::Decoded { equ_id } => {
                    g.copy_row_from(row, &self.inactive_coeff_matrix, equ_id, n);
                }
                VariableStatus::Inactive { seq } => {
                    g.set(row, seq);
                }
                VariableStatus::Active => {
                    unreachable!("tilde_g_packed requires BP complete (no active vars)")
                }
            }
        }
        g
    }
}
// ---------- Inactivation strategies ----------
impl SparseSystem {
    fn next_inactivation(&mut self) -> Option<usize> {
        if self.num_inactive < self.max_inactive_num
            && self.sparse_equations.len() >= self.params.num_message_ldpc()
            && self.num_decoded < self.params.num_total() - self.num_inactive
        {
            match self.inac_strategy {
                InactivationStrategy::ByIndex => self.inact_strategy_by_index(),
                InactivationStrategy::FirstActive => self.inact_strategy_first_active(),
                InactivationStrategy::LastActive => self.inact_strategy_last_active(),
                InactivationStrategy::MinDegree => self.inact_strategy_min_degree(),
                InactivationStrategy::TrailRun => self.inact_strategy_trail_run(),
            }
        } else {
            None
        }
    }

    fn inact_strategy_first_active(&mut self) -> Option<usize> {
        Some(self.var_ids.ids[self.num_decoded])
    }

    fn inact_strategy_by_index(&mut self) -> Option<usize> {
        let active_end = self.params.num_total() - self.num_inactive;
        for i in 0..self.params.num_active() {
            let pos = self.var_ids.pos[i];
            if pos < self.num_decoded || pos >= active_end {
                continue;
            }
            return Some(self.var_ids.ids[pos]);
        }
        None
    }

    fn inact_strategy_last_active(&mut self) -> Option<usize> {
        Some(self.var_ids.ids[self.params.num_total() - self.num_inactive - 1])
    }

    fn inact_strategy_min_degree(&mut self) -> Option<usize> {
        let mut min_degree = usize::MAX;
        let mut min_equ_id = usize::MAX;
        for pos in self.num_decoded..self.sparse_equations.len() {
            let equ_id = self.equation_ids.ids[pos];
            let equ = &self.sparse_equations[equ_id];
            if equ.degree() < min_degree {
                min_degree = equ.degree();
                min_equ_id = equ_id;
                if min_degree <= 2 {
                    break;
                }
            }
        }
        if min_degree >= 2 {
            Some(self.sparse_equations[min_equ_id].active_var_ids[0])
        } else {
            eprintln!("Unexpected: min_degree is {}, less than 2.", min_degree);
            None
        }
    }

    fn inact_strategy_trail_run(&mut self) -> Option<usize> {
        let mut graph = self.build_tanner_graph_for_trial_bp();
        if graph.size() == 0 {
            return self.inact_strategy_last_active();
        }
        graph.best_inactivation_by_trial_run()
    }

    /// Build a Tanner graph over active variables and unused equations for trial BP
    /// (doc: "generate a Tanner graph from the representation of the system").
    /// Only active variable nodes and check nodes are included; inactive coefficients
    /// are not needed for trial runs.
    pub fn build_tanner_graph_for_trial_bp(&self) -> TannerGraphDegree2 {
        let start = self.num_decoded;
        let end = self.params.num_total() - self.num_inactive;
        let active_var_ids: Vec<usize> = self.var_ids.ids[start..end].to_vec();
        let equ_deg_2: Vec<Vec<usize>> = self.equation_ids.ids[self.num_decoded..]
            .iter()
            .filter(|&equ_id| self.sparse_equations[*equ_id].degree() == 2)
            .map(|&equ_id| {
                self.sparse_equations[equ_id]
                    .active_var_ids
                    .iter()
                    .copied()
                    .collect()
            })
            .collect();
        TannerGraphDegree2::from_active_and_equations(&active_var_ids, &equ_deg_2)
    }
}
// ---------- Back substitution ----------
impl SparseSystem {
    /// Precompute decoded-equation peer lists for [`Self::reverse_bp`] / [`Self::forward_bp`].
    ///
    /// Call once when BP completes (e.g. at phase change). Clears any previous cache.
    pub fn build_bs_peers(&mut self) {
        let n = self.num_decoded;
        let mut peers: Vec<VarAdjacentList> = (0..n).map(|_| VarAdjacentList::new()).collect();
        for (pos, peer) in peers.iter_mut().enumerate().take(n.saturating_sub(1)) {
            let var_id = self.var_ids.ids[pos];
            for &equ_id_adj in &self.var_adjacent[var_id] {
                let ep = self.equation_ids.pos[equ_id_adj];
                if ep > pos && ep < n {
                    peer.push(equ_id_adj);
                }
            }
        }
        self.bs_peers = Some(peers);
    }

    /// Drop cached BS peer lists (optional, after back-substitution).
    pub fn clear_bs_peers(&mut self) {
        self.bs_peers = None;
    }

    fn broadcast_src_to_scratch_targets(&mut self, manager: &mut DataManager, src_id: usize) {
        match self.scratch_targets.len() {
            0 => {}
            1 => manager.add_one_to_vector(src_id, self.scratch_targets[0]),
            2 => {
                manager.add_one_to_vector(src_id, self.scratch_targets[0]);
                manager.add_one_to_vector(src_id, self.scratch_targets[1]);
            }
            3 => {
                manager.add_one_to_vector(src_id, self.scratch_targets[0]);
                manager.add_one_to_vector(src_id, self.scratch_targets[1]);
                manager.add_one_to_vector(src_id, self.scratch_targets[2]);
            }
            _ => manager.broadcast_add_owned(src_id, self.scratch_targets.split_off(0)),
        }
    }

    /// Iterate decoded variables as `(var_id, equ_id)` pairs.
    /// Use [`BinaryMatrix::iter_set_bits`] on [`Self::inactive_coeff_matrix`] for sparse coeffs.
    pub fn iter_decoded_var(&self) -> impl Iterator<Item = (usize, usize)> {
        self.var_ids.ids[..self.num_decoded]
            .iter()
            .enumerate()
            .map(move |(pos, &var_id)| {
                let equ_id = self.equation_ids.ids[pos];
                (var_id, equ_id)
            })
    }

    /// Iterate through all inactive variables and return the variable id and its
    /// GE / LU solution column index.
    pub fn iter_inactive_var(&self) -> impl Iterator<Item = (usize, usize)> {
        let inactive_start = self.params.num_total() - self.num_inactive;
        self.var_ids.ids[inactive_start..self.params.num_total()]
            .iter()
            .enumerate()
            .map(|(pos, &var_id)| (var_id, self.inactive_matrix_seq(pos)))
    }

    /// Reverse BP decoding of the decoded coded vectors before back substitution.
    pub fn reverse_bp(&mut self, manager: &mut DataManager) {
        let peers = self
            .bs_peers
            .take()
            .expect("build_bs_peers must be called before reverse_bp");
        let num_inactive = self.num_inactive;
        for pos in (0..self.num_decoded - 1).rev() {
            let var_id = self.var_ids.ids[pos];
            let equ_id = self.equation_ids.ids[pos];
            //let mut adj_data_ids = Vec::new();
            self.scratch_targets.clear();

            for &equ_id_adj in &peers[pos] {
                self.inactive_coeff_matrix
                    .xor_rows(equ_id_adj, equ_id, num_inactive);
                //adj_data_ids.push(eq_adj.data_id);
                self.scratch_targets
                    .push(self.sparse_equations[equ_id_adj].data_id);
            }

            let var_data_id = manager.data_id_of_variable_vector(var_id);
            //manager.broadcast_add_owned(var_data_id, adj_data_ids);
            self.broadcast_src_to_scratch_targets(manager, var_data_id);
        }
        self.bs_peers = Some(peers);
    }

    /// Forward BP decoding of the decoded coded vectors after back substitution.
    pub fn forward_bp(&mut self, manager: &mut DataManager) {
        let peers = self
            .bs_peers
            .take()
            .expect("build_bs_peers must be called before forward_bp");
        for (pos, peer) in peers.iter().enumerate().take(self.num_decoded - 1) {
            let var_id = self.var_ids.ids[pos];
            self.scratch_targets.clear();
            for &equ_id_adj in peer {
                self.scratch_targets
                    .push(self.sparse_equations[equ_id_adj].data_id);
            }

            let var_data_id = manager.data_id_of_variable_vector(var_id);
            self.broadcast_src_to_scratch_targets(manager, var_data_id);
        }
        self.bs_peers = Some(peers);
    }
}

#[cfg(test)]
mod inactive_seq_tests {
    use super::*;
    use crate::types::{CodeParams, InactivationStrategy};

    #[test]
    fn pre_inactive_vars_use_internal_matrix_columns() {
        let params = CodeParams::new(5, 2, 0, 0);
        let sys = SparseSystem::new(&params, 8, InactivationStrategy::ByIndex);
        for (var_id, expected_internal) in [(2, 0), (3, 1), (4, 2)] {
            match sys.var_status(var_id) {
                VariableStatus::Inactive { seq } => assert_eq!(seq, expected_internal),
                other => panic!("var {var_id} expected inactive, got {other:?}"),
            }
        }
    }

    #[test]
    fn dynamic_inactive_internal_seq_matches_inactivate_flip() {
        let params = CodeParams::new(5, 2, 0, 0);
        let mut sys = SparseSystem::new(&params, 8, InactivationStrategy::ByIndex);

        sys.inactivate(0);
        match sys.var_status(0) {
            VariableStatus::Inactive { seq } => assert_eq!(seq, params.num_inactive()),
            other => panic!("var 0 expected inactive after first inactivation, got {other:?}"),
        }
        assert_eq!(sys.num_inactive(), params.num_inactive() + 1);

        sys.inactivate(1);
        match sys.var_status(1) {
            VariableStatus::Inactive { seq } => {
                assert_eq!(seq, params.num_inactive() + 1);
            }
            other => panic!("var 1 expected inactive after second inactivation, got {other:?}"),
        }
    }

    #[test]
    fn iter_inactive_var_uses_ge_lu_column_indices() {
        let params = CodeParams::new(5, 2, 0, 0);
        let mut sys = SparseSystem::new(&params, 8, InactivationStrategy::ByIndex);
        sys.inactivate(0);
        sys.inactivate(1);

        let mut ge_cols: Vec<_> = sys.iter_inactive_var().map(|(_, col)| col).collect();
        ge_cols.sort_unstable();
        assert_eq!(ge_cols, (0..sys.num_inactive()).collect::<Vec<_>>());
    }
}
