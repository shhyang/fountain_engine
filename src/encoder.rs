// Copyright (c) 2025 Shenghao Yang.
// All rights reserved.

use crate::data_manager::DataManager;
//use crate::types::Operation;
use crate::core::{precode_encode, solver::Solver};
use crate::traits::{CodeScheme, DataOperator};
use crate::types::{CodeParams, CodeType, DecodeStatus, DegreeSetFn, SolverType};

/// *Fountain Code Encoder*
/// The encoder is used to encode the message vectors with optional precoding.
/// To use the encoder, a data manager is needed, which implements the `DataManager` trait.
pub struct Encoder {
    params: CodeParams,
    pub manager: DataManager,
    gen_degree_set: DegreeSetFn,
    code_type: CodeType,
}

impl Encoder {
    /// Create a new encoder with the given code configuration.
    pub fn new<T: CodeScheme>(custom: &T) -> Self {
        let manager = DataManager::new();
        Self::initialize(custom, manager)
    }

    /// Configure the encoder and data manager only; does **not** run precoding.
    pub fn new_without_precoding<T: CodeScheme>(custom: &T) -> Self {
        let manager = DataManager::new();
        Self::initialize_without_precoding(custom, manager)
    }

    pub fn new_with_operator<T: CodeScheme>(custom: &T, operator: Box<dyn DataOperator>) -> Self {
        let manager = DataManager::new_with_operator(operator);
        Self::initialize(custom, manager)
    }

    /// Like [`Self::new_with_operator`], but does not record operations (execute-only).
    pub fn new_with_operator_execute_only<T: CodeScheme>(
        custom: &T,
        operator: Box<dyn DataOperator>,
    ) -> Self {
        let manager = DataManager::new_with_operator_execute_only(operator);
        Self::initialize(custom, manager)
    }

    /// Like [`Self::new_with_operator`], but leaves precoding to [`Self::precode_encode`].
    pub fn new_with_operator_without_precoding<T: CodeScheme>(
        custom: &T,
        operator: Box<dyn DataOperator>,
    ) -> Self {
        let manager = DataManager::new_with_operator(operator);
        Self::initialize_without_precoding(custom, manager)
    }

    fn initialize<T: CodeScheme>(custom: &T, manager: DataManager) -> Self {
        let mut enc = Self::initialize_without_precoding(custom, manager);
        enc.precode_encode(custom);
        enc
    }

    fn initialize_without_precoding<T: CodeScheme>(custom: &T, mut manager: DataManager) -> Self {
        let params = custom.get_params();
        let gen_degree_set = custom.create_degree_set_fn();
        let code_type = custom.code_type();
        let solver_type = match code_type {
            CodeType::Systematic => SolverType::SysEnc,
            CodeType::Ordinary => SolverType::OrdEnc,
        };
        manager.config_from(params.clone(), solver_type);

        Self {
            params,
            manager,
            gen_degree_set,
            code_type,
        }
    }

    /// Run LDPC/HDPC precoding (ordinary) or systematic encoding solve.
    ///
    /// Must be called once after [`Self::new_without_precoding`] /
    /// [`Self::new_with_operator_without_precoding`] and before LT encoding.
    pub fn precode_encode<T: CodeScheme>(&mut self, custom: &T) {
        match self.code_type {
            CodeType::Ordinary => {
                if self.params.l + self.params.h > 0 {
                    for i in (self.params.a..self.params.k).rev() {
                        self.manager.move_to(
                            i,
                            self.manager.data_id_of_inactive_variable(i - self.params.a),
                        );
                    }
                    precode_encode(&mut self.manager, &self.params, custom);
                }
            }
            CodeType::Systematic => {
                let mut solver = Solver::new(custom, &mut self.manager);
                for coded_id in 0..self.params.k {
                    let new_data_id = self.manager.coded_data_id(coded_id);
                    self.manager.copy_to(coded_id, new_data_id);
                    solver.add_coded_vector(&mut self.manager, coded_id, new_data_id);
                }
                if solver.status == DecodeStatus::NotDecoded {
                    panic!("systematic encoding failed");
                }

                for coded_id in 0..self.params.k {
                    self.manager.assign_data_id(coded_id, coded_id);
                }
            }
        }
    }

    pub fn get_data_vector(&self, data_id: usize) -> &[u8] {
        self.manager.get_data_vector(data_id)
    }

    /// Generate the next coded vector after precoding completes.
    ///
    /// Returns the coded vector's data id, or `None` if `coded_id` is invalid for ordinary encoding.
    pub fn encode_coded_vector(&mut self, coded_id: usize) -> Option<usize> {
        if coded_id < self.params.k {
            if self.code_type == CodeType::Systematic {
                return Some(coded_id);
            } else {
                dbg!("coded id {} is less than k for ordinary encoding", coded_id);
                return None;
            }
        } else if coded_id < self.params.num_total() {
            dbg!("coded id {} is out of range for encoding", coded_id);
            return None;
        }

        if coded_id < self.params.num_message_ldpc() {
            return Some(
                self.manager
                    .data_id_of_ldpc_variable(coded_id - self.params.k),
            );
        } else if coded_id < self.params.num_total() {
            return Some(
                self.manager
                    .data_id_of_hdpc_variable(coded_id - self.params.num_message_ldpc()),
            );
        }

        let data_id = self.manager.coded_data_id(coded_id);
        let degree_set = (self.gen_degree_set)(coded_id);
        let data_ids = degree_set
            .iter()
            .map(|&id| self.manager.data_id_of_variable_vector(id))
            .collect::<Vec<_>>();
        self.manager.add_to_vector(&data_ids, data_id);
        self.manager.encode_coded_vector(coded_id, data_id);
        Some(data_id)
    }
}
