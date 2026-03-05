// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

//use crate::fountain_code::degree_set::DegreeSet;
use super::bp::BPDecoder;
use crate::algebra::{finite_field::GF256, linear_algebra::Matrix};
use crate::data_manager::DataManager;
use crate::traits::CodeScheme;
use crate::traits::HDPC;
use crate::types::CodeParams;
use crate::types::DecodeStatus;

struct InactiveSolver {
    matrix_d: Vec<Vec<u8>>,
    ids_hdpc: Vec<usize>, // data vector IDs for hdpc constraints
    matrix_f: Vec<Vec<u8>>,
    ids_unused: Vec<usize>, // data vector IDs for unused coded vectors
    q: Vec<usize>,          // column permutation of the matrix
    r: usize,               // rank of the matrix
}

impl InactiveSolver {
    fn new(num_inactive: usize) -> Self {
        let matrix_d = vec![];
        let ids_hdpc = vec![];
        let matrix_f: Vec<Vec<u8>> = Vec::new();
        let ids_unused = vec![];
        let q = (0..num_inactive).collect::<Vec<_>>();
        Self {
            matrix_d,
            ids_hdpc,
            matrix_f,
            ids_unused,
            q,
            r: 0,
        }
    }

    fn add_row(&mut self, row: Vec<u8>, id: usize) {
        self.matrix_f.push(row);
        self.ids_unused.push(id);
    }

    fn lu_decomposition(&mut self, manager: &mut DataManager) -> usize {
        //dbg!("lu_decomposition", &self.matrix_f);
        //dbg!("ids_unused", &self.ids_unused);
        let (p, r) =
            Matrix::lu_decomp_incr(&GF256::default(), &mut self.matrix_f, &mut self.q, self.r);
        // remove redundancy and apply row permutation on ids_unused
        let mut ids_new = Vec::with_capacity(r);
        for &pi in &p[..r] {
            ids_new.push(self.ids_unused[pi]);
        }
        for &pi in &p[r..] {
            manager.remove(self.ids_unused[pi]);
        }
        self.ids_unused = ids_new;
        self.matrix_f.truncate(r);
        self.r = r;
        //dbg!("lu_decomposition", &r);
        r
    }

    fn has_hdpc_constraints(&self) -> bool {
        !self.matrix_d.is_empty()
    }

    fn append_hdpc_constraints(&mut self) {
        self.matrix_f.append(&mut self.matrix_d);
        self.ids_unused.append(&mut self.ids_hdpc);
    }

    // only called when matrix_a has rank equal to the number of inactive variables
    fn lu_solve_incr(&mut self, manager: &mut DataManager) -> Vec<usize> {
        // interchanging columns of the matrix_a according to the permutation vector q
        //let mut matrix_a_new = vec![vec![0u8; self.matrix_f[0].len()]; self.matrix_f.len()];
        //for i in 0..self.matrix_f.len() {
        //    for j in 0..self.matrix_f[0].len() {
        //        matrix_a_new[i][j] = self.matrix_f[i][self.q[j]];
        //    }
        //}
        //dbg!("lu_solve", &self.ids_unused);
        //lu_solve_incr!(manager, &self.matrix_f, &self.ids_unused, &self.q);

        if self.matrix_f.len() != self.ids_unused.len() {
            panic!("The number of rows in A must be equal to the number of target IDs");
        }

        let n = self.matrix_f.len();

        // Forward substitution (L part): solve Ly = b
        for j in 0..n - 1 {
            for i in j + 1..n {
                let l = self.matrix_f[i][self.q[j]];
                if l != 0 {
                    manager.mul_add(self.ids_unused[j], l, self.ids_unused[i]);
                }
            }
        }

        // Backward substitution (U part): solve Ux = y
        for j in (0..n).rev() {
            let l = self.matrix_f[j][self.q[j]];
            if l == 0 {
                panic!(
                    "Singular matrix: diagonal element at position {} is zero",
                    j
                );
            }

            // Scale the current row by the inverse of the diagonal element
            manager.divide_scalar(l, self.ids_unused[j]);

            // Subtract scaled row from previous rows
            for i in (0..j).rev() {
                let l = self.matrix_f[i][self.q[j]];
                if l != 0 {
                    manager.mul_add(self.ids_unused[j], l, self.ids_unused[i]);
                }
            }
        }

        //dbg!("lu_solve_incr", &self.ids_unused, &self.q);
        // apply column permutation on ids_unused
        let mut ids_new = vec![0; self.r];
        for i in 0..self.r {
            ids_new[self.q[i]] = self.ids_unused[i];
        }
        ids_new
    }
}

pub struct Solver {
    params: CodeParams,
    bp_decoder: BPDecoder,
    //ldpc_type: LDPCType,
    hdpc: Option<Box<dyn HDPC>>,
    inactive_solver: Option<InactiveSolver>,
    pub phase: DecodePhase,
    pub status: DecodeStatus,
}

#[derive(Debug, PartialEq)]
pub enum DecodePhase {
    BP,
    GE,
}

impl Solver {
    pub fn new<T: CodeScheme>(custom: &T, manager: &mut DataManager) -> Self {
        let params = custom.get_params();
        let gen_degree_set = custom.create_degree_set_fn();
        let (hdpc, ldpc) = custom.create_precode();
        let max_inactive_num = custom.max_inactive_num();
        let mut bp_decoder = BPDecoder::new(&params, gen_degree_set, max_inactive_num);

        if let Some(ldpc) = ldpc {
            //let ldpc = LDPCType::new(ldpc_type, params.clone());
            bp_decoder.add_ldpc_check_node(manager, ldpc);
        }
        Self {
            params,
            bp_decoder,
            hdpc,
            inactive_solver: None,
            phase: DecodePhase::BP,
            status: DecodeStatus::NotDecoded,
        }
    }

    // called when BP phase is complete and not decoded, i.e., num_inactive > 0
    fn phase_change(&mut self, manager: &mut DataManager) {
        let mut inactive_solver = InactiveSolver::new(self.bp_decoder.num_inactive());

        // append rows of F matrix
        //dbg!("iter_unused_check");
        for (row, data_id) in self.bp_decoder.iter_unused_check() {
            //dbg!("add_row", &data_id, &row);
            //dbg!("data", manager.get_data_vector(data_id));
            inactive_solver.add_row(row, data_id);
        }

        // calculate the hdpc constraints for solving inactive variables
        if let Some(hdpc) = self.hdpc.as_mut() {
            // calculate D_d X_d
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
            //for i in &hdpc_constraint_data_ids {
            //   dbg!("hdpc_constraint_data_ids", i, manager.get_data_vector(*i));
            //}
            // calculate \tilde D
            let indices_g_row = |i: usize| self.bp_decoder.indices_tilde_g_row(i);

            //for i in 0..self.params.num_message_ldpc() {
            //    dbg!(i, indices_g_row(i));
            //}

            let mut matrix_d = hdpc.mul_sparse(
                //&GF256::default(),
                &self.params,
                //self.params.num_message_ldpc(),
                self.bp_decoder.num_inactive(),
                &indices_g_row,
            );
            let end_index = self.bp_decoder.num_inactive() - self.params.h;
            for (i, row) in matrix_d.iter_mut().take(self.params.h).enumerate() {
                row[end_index + i] ^= 1;
            }
            //dbg!("matrix_d", &matrix_d);

            inactive_solver.matrix_d = matrix_d;
            inactive_solver.ids_hdpc = hdpc_constraint_data_ids;
        }
        // enter GE phase
        self.inactive_solver = Some(inactive_solver);
        self.phase = DecodePhase::GE;
    }

    fn solve_inactive(&mut self, manager: &mut DataManager) -> DecodeStatus {
        let num_inactive = self.bp_decoder.num_inactive();
        let inactive_solver = self.inactive_solver.as_mut().unwrap();
        let mut r = inactive_solver.lu_decomposition(manager);

        //dbg!("r", r, num_inactive);
        if r < num_inactive
            && num_inactive - r <= self.params.h
            && inactive_solver.has_hdpc_constraints()
        {
            //dbg!("append_hdpc_constraints");
            inactive_solver.append_hdpc_constraints();
            r = inactive_solver.lu_decomposition(manager);
        }

        if r < num_inactive {
            return DecodeStatus::NotDecoded;
        }

        let cod_ids_sol_inactive = inactive_solver.lu_solve_incr(manager);

        for (var_id, seq, idx) in self.bp_decoder.iter_inactive_var() {
            //dbg!("solved inactive", &var_id, &cod_ids_sol_inactive[seq], &idx);
            let data_ids = idx
                .iter()
                .map(|&id| manager.data_id_of_variable_vector(id))
                .collect::<Vec<_>>();
            manager.broadcast_add(cod_ids_sol_inactive[seq], &data_ids);
            //dbg!("broadcast_add", manager.get_data_vector(cod_ids_sol_inactive[seq]));
            //if var_id < self.params.k {
            manager.move_to(
                cod_ids_sol_inactive[seq],
                manager.data_id_of_variable_vector(var_id),
            );
            //}
        }

        DecodeStatus::Decoded
    }

    pub fn add_coded_vector(&mut self, manager: &mut DataManager, coded_id: usize, data_id: usize) {
        //let data_id = self.params.solver_coded_id_to_data_id(coded_id);
        let check_id = self.bp_decoder.add_coded_vector(manager, coded_id, data_id);

        self.status = match self.phase {
            DecodePhase::BP => {
                let num_decoded_vars = self.bp_decoder.run(manager);
                if num_decoded_vars < self.params.num_total() {
                    // not decoded,check for phase change
                    if self.bp_decoder.is_complete() {
                        //dbg!("phase_change");
                        self.phase_change(manager);
                        self.solve_inactive(manager)
                    } else {
                        DecodeStatus::NotDecoded
                    }
                } else {
                    DecodeStatus::Decoded
                }
            }
            DecodePhase::GE => {
                let row = self.bp_decoder.inactive_adj_matrix_row(check_id);
                self.inactive_solver.as_mut().unwrap().add_row(row, data_id);
                self.solve_inactive(manager)
            }
        };
    }
}
