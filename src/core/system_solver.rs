// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

//! System-based decoder that works on `MasterSystem` and `InactiveSystem`.
//!
//! Uses the sparse master-system representation (doc sec:representation_of_the_system)
//! and the inactive system F B_i = X_u' (doc sec:solv_inactive).

use super::inactivation::InactiveSystem;
use super::sparse_equation::SparseSystem;
use crate::data_manager::DataManager;
use crate::traits::{CodeScheme, HDPC};
use crate::types::{CodeParams, DecodeStatus, DegreeSetFn, SolverType, SubstitutionMethod};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolvingPhase {
    BP,
    GE,
    BS,
}

pub struct SystemSolver {
    params: CodeParams,
    hdpc: Option<Box<dyn HDPC>>,
    sparse_system: SparseSystem,
    gen_degree_set: DegreeSetFn,
    inactive_system: InactiveSystem,
    pub phase: SolvingPhase,
    pub status: DecodeStatus,
    bs_method: SubstitutionMethod,
    /// True after GE `add_row_from_master`; HDPC path rebuilds only when set.
    inactive_dirty: bool,
}

impl SystemSolver {
    pub fn new<T: CodeScheme>(
        custom: &T,
        manager: &mut DataManager,
        solver_type: SolverType,
    ) -> Self {
        let params = custom.get_params();
        let gen_degree_set = custom.create_degree_set_fn();
        let (hdpc, ldpc_opt) = custom.create_precode();
        let decoding_config = custom.decoding_config();
        let max_inactive_num = decoding_config.max_inactive_num;
        let inac_strategy = decoding_config.inac_strategy;
        let bs_method = decoding_config.subs_method;

        // degree_set generator is not stored in `MasterSystem` for now; we only keep it here.
        let mut master = SparseSystem::new(&params, max_inactive_num, inac_strategy);

        if let Some(ldpc) = ldpc_opt {
            master.add_ldpc_constraints(manager, ldpc.as_ref());
        }

        if let Some(hdpc_box) = &hdpc {
            // See `precode_encode`: `GF2_FIELD_POLY` for binary HDPC, else GF(256).
            manager.config_finite_field(hdpc_box.gf_poly());
        }

        let inactive_system = InactiveSystem::new(&params);

        let num_k = params.k;
        let num_a = params.a;

        let mut solver = Self {
            params,
            hdpc,
            sparse_system: master,
            gen_degree_set,
            inactive_system,
            phase: SolvingPhase::BP,
            status: DecodeStatus::NotDecoded,
            bs_method,
            inactive_dirty: false,
        };

        match solver_type {
            SolverType::OrdEnc => {
                unreachable!()
            }
            SolverType::OrdDec => {
                //add padding vectors as decoded variables
                for var_id in num_a - manager.num_padding()..num_a {
                    let new_data_id = manager.coded_data_id(var_id);
                    manager.ensure_zero_one(new_data_id);
                    solver
                        .sparse_system
                        .add_lt_coded_vector(manager, new_data_id, &[var_id]);
                }
                let _ = solver.sparse_system.run_bp_inactivation(manager);
            }
            SolverType::SysEnc => {
                for coded_id in 0..manager.num_source() {
                    let new_data_id = manager.coded_data_id(coded_id);
                    manager.copy_to(coded_id, new_data_id);
                    solver.add_coded_vector(manager, coded_id, new_data_id);
                }
                for coded_id in manager.num_source()..num_k {
                    let new_data_id = manager.coded_data_id(coded_id);
                    manager.ensure_zero_one(new_data_id);
                    solver.add_coded_vector(manager, coded_id, new_data_id);
                }
            }
            SolverType::SysDec => {
                //add padding vectors as all-zero received coded vectors
                for coded_id in manager.num_source()..num_k {
                    let new_data_id = manager.coded_data_id(coded_id);
                    manager.ensure_zero_one(new_data_id);
                    solver.add_coded_vector(manager, coded_id, new_data_id);
                }
            }
        }

        solver
    }

    /// Solve the inactive system.
    fn solve_inactive(&mut self, manager: &mut DataManager) -> DecodeStatus {
        //let remove_redundant = !matches!(self.phase, SolvingPhase::GE);
        let remove_redundant = false;
        let mut r = self
            .inactive_system
            .lu_decomposition(manager, remove_redundant);
        let num_inactive = self.sparse_system.num_inactive();
        if r < num_inactive && self.hdpc.is_some() && num_inactive - r <= self.params.h {
            if self.inactive_dirty {
                self.rebuild_inactive_from_unused();
                r = self
                    .inactive_system
                    .lu_decomposition(manager, remove_redundant);
            }
            // G4.2: skip HDPC append when packed GF(2) rows already have full rank.
            if r < num_inactive {
                if manager.gf256().is_none() {
                    let tilde_g = self.sparse_system.tilde_g_packed();
                    self.inactive_system.append_hdpc_constraints_packed(
                        self.hdpc.as_deref().unwrap(),
                        &tilde_g,
                        manager,
                    );
                } else {
                    let tilde_g_rows = self.sparse_system.tilde_g_rows();
                    self.inactive_system.append_hdpc_constraints(
                        self.hdpc.as_deref().unwrap(),
                        &tilde_g_rows,
                        manager,
                    );
                }
                r = self
                    .inactive_system
                    .lu_decomposition(manager, remove_redundant);
            }
            self.hdpc = None;
        }
        if r < num_inactive {
            return DecodeStatus::NotDecoded;
        }
        if !self.inactive_system.has_full_rank_factorization() {
            self.rebuild_inactive_from_unused();
            r = self
                .inactive_system
                .lu_decomposition(manager, remove_redundant);
            if r < num_inactive {
                return DecodeStatus::NotDecoded;
            }
        }
        let data_ids_sol_inactive = self.inactive_system.lu_solve(manager);

        // enter BS phase
        self.phase = SolvingPhase::BS;

        match self.bs_method {
            SubstitutionMethod::Direct => {
                self.direct_back_substitution(manager, data_ids_sol_inactive);
            }
            SubstitutionMethod::Original => {
                self.original_back_substitution(manager, data_ids_sol_inactive);
            }
        }

        DecodeStatus::Decoded
    }

    fn direct_back_substitution(
        &mut self,
        manager: &mut DataManager,
        data_ids_sol_inactive: Vec<usize>,
    ) {
        let num_inactive = self.sparse_system.num_inactive();
        let coeff_matrix = self.sparse_system.inactive_coeff_matrix();
        let mut inactive_data_ids = Vec::with_capacity(8);
        for (var_id, equ_id) in self.sparse_system.iter_decoded_var() {
            let target = manager.data_id_of_variable_vector(var_id);
            let mut iter = coeff_matrix.iter_set_bits(equ_id, num_inactive);
            match iter.next() {
                None => {}
                Some(seq0) => {
                    let id0 = data_ids_sol_inactive[seq0];
                    match iter.next() {
                        None => manager.add_one_to_vector(id0, target),
                        Some(seq1) => {
                            let id1 = data_ids_sol_inactive[seq1];
                            match iter.next() {
                                None => manager.add_two_to_vector(id0, id1, target),
                                Some(seq2) => {
                                    let id2 = data_ids_sol_inactive[seq2];
                                    match iter.next() {
                                        None => manager.add_three_to_vector(id0, id1, id2, target),
                                        Some(seq3) => {
                                            inactive_data_ids.clear();
                                            inactive_data_ids.push(id0);
                                            inactive_data_ids.push(id1);
                                            inactive_data_ids.push(id2);
                                            inactive_data_ids.push(data_ids_sol_inactive[seq3]);
                                            for seq in iter {
                                                inactive_data_ids.push(data_ids_sol_inactive[seq]);
                                            }
                                            manager.add_to_vector(&inactive_data_ids, target);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        for (var_id, pos) in self.sparse_system.iter_inactive_var() {
            //dbg!("inactive variable", var_id, pos);
            //dbg!(manager.get_data_vector(data_ids_sol_inactive[pos]));
            manager.move_to(
                data_ids_sol_inactive[pos],
                manager.data_id_of_variable_vector(var_id),
            );
        }
    }

    fn original_back_substitution(
        &mut self,
        manager: &mut DataManager,
        data_ids_sol_inactive: Vec<usize>,
    ) {
        self.sparse_system.reverse_bp(manager);
        self.direct_back_substitution(manager, data_ids_sol_inactive);
        self.sparse_system.forward_bp(manager);
        self.sparse_system.clear_bs_peers();
    }

    /// Rebuild the inactive system from all unused sparse equations (fresh matrix, `r = 0`).
    fn rebuild_inactive_from_unused(&mut self) {
        self.inactive_system.init(self.sparse_system.num_inactive());
        let mut add_row = |master: &crate::core::binary_matrix::BinaryMatrix, equ_id, data_id| {
            self.inactive_system
                .add_row_from_master(master, equ_id, data_id);
        };
        self.sparse_system
            .inactive_system_from_unused_equations_packed(&mut add_row);
        self.inactive_dirty = false;
    }

    /// Called when BP phase on `MasterSystem` is complete and there are inactive variables.
    fn phase_change(&mut self, _manager: &DataManager) {
        self.rebuild_inactive_from_unused();
        self.sparse_system.build_bs_peers();
        self.phase = SolvingPhase::GE;
    }

    /// Add one coded vector and advance the BP/GE phases.
    pub fn add_coded_vector(&mut self, manager: &mut DataManager, coded_id: usize, data_id: usize) {
        let degree_set = (self.gen_degree_set)(coded_id);

        let equ_id = self
            .sparse_system
            .add_lt_coded_vector(manager, data_id, &degree_set);

        self.status = match self.phase {
            SolvingPhase::BP => {
                let num_decoded = self.sparse_system.run_bp_inactivation(manager);
                if num_decoded < self.params.num_total() {
                    if self.sparse_system.is_bp_complete() {
                        self.phase_change(manager);
                        self.solve_inactive(manager)
                    } else {
                        DecodeStatus::NotDecoded
                    }
                } else {
                    DecodeStatus::Decoded
                }
            }
            SolvingPhase::GE => {
                self.inactive_system.add_row_from_master(
                    self.sparse_system.inactive_coeff_matrix(),
                    equ_id,
                    data_id,
                );
                self.inactive_dirty = true;
                self.solve_inactive(manager)
            }
            SolvingPhase::BS => DecodeStatus::Decoded,
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_manager::DataManager;
    use crate::traits::{CodeScheme, PrecodePair};
    use crate::types::{CodeParams, CodeType, DecodingConfig, DegreeSetFn, Operation, SolverType};

    #[derive(Clone)]
    struct DummySystematic {
        k: usize,
    }

    impl CodeScheme for DummySystematic {
        fn get_params(&self) -> CodeParams {
            CodeParams::new(self.k, self.k, 0, 0)
        }

        fn code_type(&self) -> CodeType {
            CodeType::Systematic
        }

        fn create_degree_set_fn(&self) -> DegreeSetFn {
            Box::new(|_| vec![0])
        }

        fn create_precode(&self) -> PrecodePair {
            (None, None)
        }

        fn decoding_config(&self) -> DecodingConfig {
            DecodingConfig::default()
        }
    }

    fn ensure_zero_ids(mgr: &DataManager) -> Vec<usize> {
        let mut ids = Vec::new();
        for op in mgr.get_operations() {
            match op {
                Operation::EnsureZero { list_id } => ids.extend(list_id.iter().copied()),
                Operation::EnsureZeroOne { id } => ids.push(*id),
                _ => {}
            }
        }
        ids
    }

    fn coded_vector_data_ids(mgr: &DataManager, coded_ids: &[usize]) -> Vec<usize> {
        mgr.get_operations()
            .iter()
            .filter_map(|op| match op {
                Operation::InfoCodedVector { coded_id, data_id }
                    if coded_ids.contains(coded_id) =>
                {
                    Some(*data_id)
                }
                _ => None,
            })
            .collect()
    }

    #[test]
    fn sysenc_zero_fills_padding_coded_slots() {
        let scheme = DummySystematic { k: 12 };
        let mut mgr = DataManager::new();
        mgr.config_from(scheme.get_params(), SolverType::SysEnc);
        mgr.set_num_source(10);
        let _ = SystemSolver::new(&scheme, &mut mgr, SolverType::SysEnc);
        let zeroed = ensure_zero_ids(&mgr);
        for data_id in coded_vector_data_ids(&mgr, &[10, 11]) {
            assert!(
                zeroed.contains(&data_id),
                "padding coded data id {data_id} should be zeroed"
            );
        }
        assert_eq!(
            mgr.get_operations()
                .iter()
                .filter(|op| matches!(op, Operation::CopyTo { .. }))
                .count(),
            10
        );
    }

    #[test]
    fn sysenc_without_padding_copies_all_message_vectors() {
        let scheme = DummySystematic { k: 12 };
        let mut mgr = DataManager::new();
        mgr.config_from(scheme.get_params(), SolverType::SysEnc);
        let _ = SystemSolver::new(&scheme, &mut mgr, SolverType::SysEnc);
        assert_eq!(
            mgr.get_operations()
                .iter()
                .filter(|op| matches!(op, Operation::CopyTo { .. }))
                .count(),
            12
        );
        let zeroed = ensure_zero_ids(&mgr);
        for data_id in coded_vector_data_ids(&mgr, &[10, 11]) {
            assert!(
                !zeroed.contains(&data_id),
                "payload coded data id {data_id} should come from copy_to, not ensure_zero"
            );
        }
    }
}
