// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use crate::data_manager::DataManager;
//use crate::types::Operation;
use crate::core::{precode_encode, Solver};
use crate::traits::{CodeScheme, DataOperator};
use crate::types::{CodeParams, CodeType, DecodeStatus, DegreeSetFn, SolverType};

/// *Fountain Code Encoder*
/// The encoder is used to encode the message vectors with optional precoding.
/// To use the encoder, a data manager is needed, which implements the `DataManager` trait.
pub struct Encoder {
    params: CodeParams,
    /// Manages all data vectors (message, coded, precode) and their operations.
    pub manager: DataManager,
    gen_degree_set: DegreeSetFn,
    code_type: CodeType,
}

impl Encoder {
    /// Create a new encoder with the given code configuration.
    pub fn new<T: CodeScheme>(custom: T) -> Self {
        let manager = DataManager::new();
        Self::initialize(custom, manager)
    }

    /// Creates a new encoder with a data operator for immediate operation execution.
    pub fn new_with_operator<T: CodeScheme>(custom: T, operator: Box<dyn DataOperator>) -> Self {
        let manager = DataManager::new_with_operator(operator);
        Self::initialize(custom, manager)
    }

    fn initialize<T: CodeScheme>(custom: T, mut manager: DataManager) -> Self {
        let params = custom.get_params();
        let gen_degree_set = custom.create_degree_set_fn();
        let code_type = custom.code_type();
        let solver_type = match code_type {
            CodeType::Systematic => SolverType::SysEnc,
            CodeType::Ordinary => SolverType::OrdEnc,
        };
        manager.config_from(params.clone(), solver_type);

        match code_type {
            CodeType::Ordinary => {
                //dbg!("ordinary encode");
                if params.l + params.h > 0 {
                    //dbg!("precode encode");
                    // move the inactive message vectors to the inactive message variable data vectors
                    for i in (params.a..params.k).rev() {
                        manager.move_to(i, manager.data_id_of_inactive_variable(i - params.a));
                    }
                    precode_encode(&mut manager, &params, &custom);
                }
            }
            CodeType::Systematic => {
                //dbg!("systematic encoding");
                let mut solver = Solver::new(&custom, &mut manager);
                // solve active vectors
                //for coded_id in 0..params.a {
                //    let new_data_id = manager.data_id_of_active_variable(coded_id);
                //    manager.copy_to(coded_id, new_data_id);
                //    solver.add_coded_vector(manager, coded_id, new_data_id);
                //}
                //assert!(solver.phase == DecodePhase::GE);
                // solve inactive vectors
                // for coded_id in params.a..params.k {
                for coded_id in 0..params.k {
                    let new_data_id = manager.coded_data_id(coded_id);
                    manager.copy_to(coded_id, new_data_id);
                    solver.add_coded_vector(&mut manager, coded_id, new_data_id);
                }
                if solver.status == DecodeStatus::NotDecoded {
                    panic!("systematic encoding failed");
                }

                for coded_id in 0..params.k {
                    manager.assign_data_id(coded_id, coded_id);
                }
            }
        }

        Self {
            params,
            manager,
            gen_degree_set,
            code_type,
        }
    }

    // (Commented-out data manager functions; doc above referred to them.)
    /*
    pub fn get_operations(&self) -> &[Operation] {
        self.manager.get_operations()
    }

    pub fn clear_operations(&mut self) {
        self.manager.clear_operations()
    }

    pub fn move_new_operations(&mut self) -> Vec<Operation> {
        self.manager.move_new_operations()
    }

    pub fn move_operator(&mut self) -> Box<dyn DataOperator> {
        self.manager.move_operator()
    }
    */
    //pub fn set_operator(&mut self, operator: Box<dyn DataOperator>) {
    //    self.manager.set_operator(operator);
    //}

    /// Returns a reference to the data vector at the given ID via the data manager.
    pub fn get_data_vector(&self, data_id: usize) -> &[u8] {
        self.manager.get_data_vector(data_id)
    }

    /// Generate the next coded vector. This function is supposed to be called after the precoding is done.
    ///
    /// # Arguments
    ///
    /// * `manager` - The data manager to use for encoding.
    ///   Return the data id of the coded vector.
    pub fn encode_coded_vector(&mut self, coded_id: usize) -> Option<usize> {
        if coded_id < self.params.k {
            //if let Some(data_ids) = self.msg_vec_data_ids.as_ref() {
            //    return data_ids[coded_id];
            //} else { // generate a new message vector for systematic encoding
            if self.code_type == CodeType::Systematic {
                return Some(coded_id);
            } else {
                eprintln!("coded id {} is less than k for ordinary encoding", coded_id);
                return None;
            }
            //}
        } else if coded_id < self.params.num_total() {
            eprintln!("coded id {} is out of range for encoding", coded_id);
            return None;
        }

        if coded_id < self.params.num_message_ldpc() {
            return Some(self.manager.data_id_of_ldpc_variable(coded_id - self.params.k));
        } else if coded_id < self.params.num_total() {
            return Some(self.manager.data_id_of_hdpc_variable(coded_id - self.params.num_message_ldpc()));
        } 

        // coded_id >= self.params.num_total()

        let data_id = self.manager.coded_data_id(coded_id);
        //dbg!("encode_coded_vector", &coded_id, &data_id);
        let (active_indices, inactive_indices) = (self.gen_degree_set)(coded_id);
        let active_data_ids = active_indices
            .iter()
            .map(|&id| self.manager.data_id_of_active_variable(id))
            .collect::<Vec<_>>();
        let inactive_data_ids = inactive_indices
            .iter()
            .map(|&id| self.manager.data_id_of_inactive_variable(id))
            .collect::<Vec<_>>();
        //dbg!("encode_coded_vector", &coded_id, &active_indices, &inactive_indices);
        //manager.ensure_zero(&[data_id]);
        self.manager.add_to_vector(&active_data_ids, data_id);
        self.manager.add_to_vector(&inactive_data_ids, data_id);
        self.manager.encode_coded_vector(coded_id, data_id);
        //dbg!("encode_coded_vector", manager.get_vector(data_id));
        Some(data_id)
    }
}
