// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use crate::data_manager::DataManager;
use crate::traits::LDPC;
use crate::types::{CodeParams, DegreeSetFn};
use std::collections::VecDeque;

const CHECK_NODE_MAX_DEGREE: usize = 30;

// Tanner graph implementation for fountain code decoding.

#[derive(PartialEq, Eq, Debug)]
/// Variable node state in the Tanner graph.
pub enum VarState {
    Active,
    Inactive { seq: usize },
    Decoded { check_id: usize },
}

#[derive(PartialEq, Eq)]
pub enum CheckState {
    Decoded,
    Unused,
}

/// Check node in the Tanner graph
pub struct CheckNode {
    //id: usize, // id of check node is consecutive from 0 to k+l+h-1
    //pub coded_id: usize, // the id of coded vector of a check node is not necessarily consecutive, and starts from k+l+h.
    data_id: usize, // the id of data vector of a check node, which is not necessarily the same as the coded_id.
    state: CheckState,
    adjacent: Vec<usize>,
    degree: usize, // Cached degree for performance
}

impl CheckNode {
    fn degree(&self) -> usize {
        self.degree
    }

    fn add_adjacent(&mut self, var_id: usize) {
        if !self.adjacent.contains(&var_id) {
            self.adjacent.push(var_id);
            self.degree += 1;
        }
    }

    fn remove_adjacent(&mut self, var_id: usize) {
        if let Some(pos) = self.adjacent.iter().position(|&x| x == var_id) {
            self.adjacent.swap_remove(pos);
            self.degree -= 1;
        }
    }

    // Increment the cached degree
    /*fn increment_degree(&mut self) {
        self.degree += 1;
    }*/

    // Decrement the cached degree
    /*fn decrement_degree(&mut self) {
        if self.degree > 0 {
            self.degree -= 1;
        } else {
            panic!("Degree cannot be decremented below 0");
        }
    }*/
}

/// Tanner graph implementation for fountain code decoding
pub struct BPDecoder {
    params: CodeParams,
    gen_degree_set: DegreeSetFn,
    //degree_set:  Box<dyn DegreeSet>,
    variable_states: Vec<VarState>,
    variable_adjacent: Vec<Vec<usize>>,
    inactive_adj_matrix: Vec<Vec<u8>>,
    checks: Vec<CheckNode>,
    decodable_checks: VecDeque<usize>,
    num_decoded_vars: usize, // number of decoded variable nodes
    next_inactive_seq: usize,
    max_inactive_num: usize,
}

impl BPDecoder {
    /// Create a new Tanner graph with k variable nodes
    pub fn new(params: &CodeParams, gen_degree_set: DegreeSetFn, max_inactive_num: usize) -> Self {
        let total_variables = params.num_total();
        let mut variable_states = Vec::with_capacity(total_variables);
        let mut variable_adjacent = Vec::with_capacity(total_variables);

        // Initialize variable nodes
        for _i in 0..params.a {
            // active message variables
            variable_states.push(VarState::Active);
            variable_adjacent.push(Vec::new());
        }
        for _i in params.k..params.num_message_ldpc() {
            // ldpc variables
            variable_states.push(VarState::Active);
            variable_adjacent.push(Vec::new());
        }
        for i in params.a..params.k {
            // inactive message variables
            variable_states.push(VarState::Inactive {
                seq: max_inactive_num - params.num_inactive() + (i - params.a),
            });
            variable_adjacent.push(Vec::new());
        }

        for i in params.num_message_ldpc()..params.num_total() {
            // hdpc variables
            variable_states.push(VarState::Inactive {
                seq: max_inactive_num - params.num_inactive() + (i - params.num_active()),
            });
            variable_adjacent.push(Vec::new());
        }

        Self {
            params: params.clone(),
            gen_degree_set,
            variable_states,
            variable_adjacent,
            inactive_adj_matrix: Vec::new(),
            checks: Vec::new(),
            decodable_checks: VecDeque::new(),
            num_decoded_vars: 0,
            next_inactive_seq: 0,
            max_inactive_num,
        }
    }

    pub fn add_ldpc_check_node(&mut self, manager: &mut DataManager, ldpc: Box<dyn LDPC>) {
        //, data_id_list: Vec<usize>) {
        //let coded_id_list = vec![0; self.params.l];
        for ldpc_id in 0..self.params.l {
            let data_id = manager.coded_data_id(ldpc_id + self.params.k); //self.params.ldpc_id_to_data_id(ldpc_id); //data_id_list[ldpc_id];
            let check_id = self.add_check_node(data_id);
            // add edges between the check node and active message nodes
            let active_neighbors = ldpc.active_row(ldpc_id);
            //dbg!("LDPC check_id: {}, active_neighbors: {:?}", &check_id, &active_neighbors);
            for &neighbor in &active_neighbors {
                self.add_edge_active(neighbor, check_id);
            }

            // add the edge between the coded vector and the LDPC variable node
            self.add_edge_active(ldpc_id + self.params.a, check_id);

            // add edges between the check node and inactive variable nodes
            if self.params.num_inactive() > 0 {
                let inactive_neighbors = ldpc.inactive_row(ldpc_id);
                for &neighbor in &inactive_neighbors {
                    self.add_edge_inactive(neighbor + self.params.num_active(), check_id);
                }
            }
        }
    }
    /// add a check node with coded_id and return the id of the check node
    fn add_check_node(&mut self, data_id: usize) -> usize {
        self.inactive_adj_matrix
            .push(vec![0u8; self.max_inactive_num]);
        self.checks.push(CheckNode {
            //id: self.checks.len()-1,
            //coded_id,
            data_id,
            state: CheckState::Unused,
            adjacent: Vec::with_capacity(CHECK_NODE_MAX_DEGREE),
            degree: 0, // Initial degree is 0
                       //inactive_coeff: Vec::new(),
        });
        self.checks.len() - 1
    }

    fn add_decodable(&mut self, check_id: usize) {
        self.decodable_checks.push_back(check_id);
    }

    fn next_decodable(&mut self) -> Option<usize> {
        self.decodable_checks.pop_front()
    }

    pub fn num_inactive(&self) -> usize {
        self.next_inactive_seq + self.params.num_inactive()
    }

    /// Add an edge between variable node and check node
    fn add_edge_inactive(&mut self, var_id: usize, check_id: usize) {
        let seq = match self.variable_states[var_id] {
            VarState::Inactive { seq } => seq,
            _ => {
                unreachable!()
            }
        };
        self.inactive_adj_matrix[check_id][seq] ^= 1;
    }

    fn add_edge_active(&mut self, var: usize, check: usize) {
        match self.variable_states[var] {
            VarState::Active => {
                self.variable_adjacent[var].push(check);
                //if self.checks[check].adjacent.insert(var) {
                // Only increment degree if the edge was actually added (not already present)
                //    self.checks[check].increment_degree();
                //}
                self.checks[check].add_adjacent(var);
            }
            _ => {
                unreachable!()
            }
        }
    }

    fn xor_inactive_vars(&mut self, check1: usize, check2: usize) {
        for i in 0..self.max_inactive_num {
            self.inactive_adj_matrix[check2][i] ^= self.inactive_adj_matrix[check1][i];
        }
    }

    /*
    fn remove_edge_active(&mut self, var_id: usize, check_id: usize) {
        // it is not necessary to remove the edge from the variable adjacent set
        //self.variable_adjacent[var_id].remove(&check_id);
        if self.checks[check_id].adjacent.remove(&var_id) {
            // Only decrement degree if the edge was actually removed (was present)
            self.checks[check_id].decrement_degree();
        }
    } */

    /// Get the coded id of from a check node
    pub fn get_data_id(&self, check_id: usize) -> usize {
        self.checks[check_id].data_id
    }

    /*fn variable_id_from_seq(&self, seq_in: usize) -> usize {
        let seq_out = if seq_in < self.next_inactive_seq {
            seq_in
        } else {
            seq_in - self.next_inactive_seq + self.params.max_inactive_num() - self.params.num_inactive()
        };
        for i in 0..self.variables.len() {
            if let VarState::Inactive {seq} = self.variables[i].state {
                if seq == seq_out {
                    return i;
                }
            }
        }
        unreachable!()
    }*/

    fn solve_check(&mut self, manager: &mut DataManager, check_id: usize) {
        if self.checks[check_id].degree() == 0 {
            // no active variable nodes, skip
            return;
        }
        let var_id = *self.checks[check_id].adjacent.first().unwrap();
        match self.variable_states[var_id] {
            VarState::Active => self.variable_states[var_id] = VarState::Decoded { check_id },
            VarState::Decoded { check_id: _ } => {
                unreachable!()
            }
            VarState::Inactive { seq: _ } => {
                unreachable!()
            }
        }

        //dbg!("solve check", &var_id, &check_id);
        //dbg!(&self.inactive_adj_matrix_row(check_id));

        self.checks[check_id].state = CheckState::Decoded;
        self.num_decoded_vars += 1;

        let var_data_id = manager.data_id_of_variable_vector(var_id);
        // dbg!(&var_id, &self.get_data_id(check_id), &var_data_id);
        // get value of decoded variable vector from the coded vector
        manager.move_to(self.get_data_id(check_id), var_data_id);
        // dbg!("value of decoded vector", manager.get_data_vector(var_data_id));

        let check_to_subs: Vec<usize> = self.variable_adjacent[var_id].to_vec();
        for sub_id in check_to_subs {
            if sub_id == check_id {
                continue;
            }
            // perform back substitution
            manager.add_to_vector(&[var_data_id], self.get_data_id(sub_id));
            // xor inactive variables
            self.xor_inactive_vars(check_id, sub_id);
            // dbg!("inactive coefficient", &sub_id, manager.get_data_vector(self.get_data_id(sub_id)), &self.inactive_adj_matrix[sub_id]);
            // remove edge between the variable node and the check node
            //self.remove_edge_active(var_id, sub_id);
            self.checks[sub_id].remove_adjacent(var_id);
            // add to decodable checks if degree is 1
            if self.checks[sub_id].degree() == 1 {
                self.add_decodable(sub_id);
            }
        }
    }

    /// Add a coded vector to the BP decoder
    pub fn add_coded_vector(
        &mut self,
        manager: &mut DataManager,
        coded_id: usize,
        data_id: usize,
    ) -> usize {
        let degree_set = (self.gen_degree_set)(coded_id);
        let check_id = self.add_check_node(data_id);

        for &var_id in &degree_set {
            match self.variable_states[var_id] {
                VarState::Active => self.add_edge_active(var_id, check_id),
                VarState::Inactive { seq: _ } => self.add_edge_inactive(var_id, check_id),
                VarState::Decoded {
                    check_id: decoded_by,
                } => {
                    let var_data_id = manager.data_id_of_variable_vector(var_id);
                    manager.add_to_vector(&[var_data_id], data_id);
                    self.xor_inactive_vars(decoded_by, check_id);
                }
            }
        }

        //dbg!("inactive coefficient", &self.inactive_adj_matrix[check_id]);
        //dbg!("value of coded vector after adding", manager.get_vector(data_id));

        if self.checks[check_id].degree() == 1 {
            self.add_decodable(check_id);
        }
        check_id
    }

    /// Add a coded vector that connects to a single variable (degree-1 bootstrap).
    pub fn add_degree_one_coded_vector(
        &mut self,
        manager: &mut DataManager,
        data_id: usize,
        var_id: usize,
    ) -> usize {
        let check_id = self.add_check_node(data_id);
        match self.variable_states[var_id] {
            VarState::Active => self.add_edge_active(var_id, check_id),
            VarState::Inactive { seq: _ } => self.add_edge_inactive(var_id, check_id),
            VarState::Decoded {
                check_id: decoded_by,
            } => {
                let var_data_id = manager.data_id_of_variable_vector(var_id);
                manager.add_to_vector(&[var_data_id], data_id);
                self.xor_inactive_vars(decoded_by, check_id);
            }
        }
        if self.checks[check_id].degree() == 1 {
            self.add_decodable(check_id);
        }
        check_id
    }

    /// Run BP decoding
    pub fn run(&mut self, manager: &mut DataManager) -> usize {
        loop {
            if let Some(check_id) = self.next_decodable() {
                self.solve_check(manager, check_id);
            } else if let Some(var_id) = self.next_inactivation() {
                //dbg!("inactivate variable", var_id);
                self.inactivate(var_id);
            } else {
                break;
            }
        }
        self.num_decoded_vars
    }

    /// Find a variable node to inactivate.
    /// This function determines when and which variable node to inactivate.
    fn next_inactivation(&mut self) -> Option<usize> {
        if self.num_inactive() < self.max_inactive_num
            && self.checks.len() >= self.params.num_message_ldpc()
        {
            // Find a variable node to inactivate
            for i in 0..self.variable_states.len() {
                if self.variable_states[i] == VarState::Active {
                    return Some(i);
                }
            }
        }
        None
    }

    /// Inactivate a variable node
    /// Make sure to call this function only when the variable node is active.
    fn inactivate(&mut self, var_id: usize) {
        // Assign a sequence number based on current inactive count
        //dbg!("inactivate", &var_id);
        self.variable_states[var_id] = VarState::Inactive {
            seq: self.next_inactive_seq,
        };
        self.next_inactive_seq += 1;
        // update the related check nodes
        let checks_to_subs: Vec<usize> = self.variable_adjacent[var_id].to_vec();
        for sub_id in checks_to_subs {
            //self.remove_edge_active(var_id, sub_id);
            self.checks[sub_id].remove_adjacent(var_id);
            self.add_edge_inactive(var_id, sub_id);
            if self.checks[sub_id].degree() == 1 {
                self.add_decodable(sub_id);
            }
        }
    }

    /// Whether BP decoding is complete
    pub fn is_complete(&self) -> bool {
        self.num_decoded_vars + self.num_inactive() == self.params.num_total()
    }

    /// convert inactive sequence number to be consecutive
    fn inactive_seq_to_consecutive(&self, seq: usize) -> usize {
        if seq < self.max_inactive_num - self.params.num_inactive() {
            seq
        } else {
            seq + self.next_inactive_seq - (self.max_inactive_num - self.params.num_inactive())
        }
    }
    /// The following functions should be called only when the BP decoding is complete.
    /// For each variable node from 0 to k+l-1, return the non-zero indices
    pub fn indices_tilde_g_row(&self, var_id: usize) -> Vec<usize> {
        let mut indices = Vec::with_capacity(self.num_inactive());
        //let var_id = self.params.message_ldpc_vector_idx(idx);
        //if var_id > self.params.num_message_ldpc() {
        //    panic!("idx > num_message_ldpc");
        //}
        match self.variable_states[var_id] {
            VarState::Inactive { seq } => {
                indices.push(self.inactive_seq_to_consecutive(seq));
            }
            VarState::Decoded { check_id } => {
                indices = self.inactive_adj_matrix[check_id]
                    .iter()
                    .enumerate()
                    .filter(|&(_, &v)| v == 1)
                    .map(|(i, _)| self.inactive_seq_to_consecutive(i))
                    .collect()
            }
            VarState::Active => unreachable!(),
        }
        indices
    }

    // Get the inactive sequence numbers for a check node (commented out)
    // fn inactive_seqs_for_check(&self, check_id: usize) -> Vec<usize> { ... }

    /// Get a row of the inactive matrix for a check node with only the used inactive sequence columns
    pub fn inactive_adj_matrix_row(&self, check_id: usize) -> Vec<u8> {
        let row = &self.inactive_adj_matrix[check_id];
        let mut filtered = Vec::with_capacity(self.num_inactive());
        filtered.extend_from_slice(&row[..self.next_inactive_seq]);
        filtered.extend_from_slice(&row[self.max_inactive_num - self.params.num_inactive()..]);
        filtered
    }

    /// Iterate over unused check nodes. Output the inactive matrix row and the vector id.
    pub fn iter_unused_check(&mut self) -> impl Iterator<Item = (Vec<u8>, usize)> {
        self.checks
            .iter()
            .enumerate()
            .filter(|&(_, check)| check.state == CheckState::Unused)
            .map(|(i, check)| (self.inactive_adj_matrix_row(i), check.data_id))
    }

    /// Iterate over inactive variable nodes. Output the inactive sequence number and the variable node id.
    pub fn iter_inactive_var(&self) -> impl Iterator<Item = (usize, usize, Vec<usize>)> + '_ {
        self.variable_states
            .iter()
            .enumerate()
            .filter_map(move |(i, state)| {
                if let VarState::Inactive { seq } = state {
                    let idx = self.var_indices_inactive_adj_column(*seq);
                    let seq_out = self.inactive_seq_to_consecutive(*seq);
                    Some((i, seq_out, idx))
                } else {
                    None
                }
            })
    }

    /// For each inactive sequence number j, return the variable node ids i which are decoded by the check nodes c,
    /// such that the (c,j) entry of the inactive matrix is 1.
    fn var_indices_inactive_adj_column(&self, seq: usize) -> Vec<usize> {
        let mut indices = Vec::with_capacity(self.checks.len());
        for i in 0..self.inactive_adj_matrix.len() {
            if self.inactive_adj_matrix[i][seq] == 1 && self.checks[i].state == CheckState::Decoded
            {
                if self.checks[i].degree() != 1 {
                    panic!("Check node degree is not 1");
                }
                let var_id = *self.checks[i].adjacent.first().unwrap();
                //if var_id < self.params.k { // if LDPC vectors are expected to be decoded, return the variable node id
                indices.push(var_id);
                //}
            }
        }
        indices
    }
}
