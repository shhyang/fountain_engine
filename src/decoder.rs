// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use crate::core::system_solver::SystemSolver;
use crate::data_manager::DataManager;
use crate::traits::{CodeScheme, DataOperator};
use crate::types::{CodeParams, CodeType, DecodeStatus, DegreeSetFn, SolverType};

/// Fountain code decoder using BP + inactivation + Gaussian elimination.
///
/// Accepts coded vectors one at a time and attempts to recover the original
/// source vectors. Supports both systematic and ordinary decoding modes.
pub struct Decoder {
    params: CodeParams,
    /// Public access to the underlying data manager for retrieving results.
    pub manager: DataManager,
    solver: SystemSolver,
    received_msg_vectors: Option<Vec<bool>>,
    gen_degree_set: Option<DegreeSetFn>,
}

impl Decoder {
    /// Creates a new decoder from a code scheme, without a data operator.
    pub fn new<T: CodeScheme>(custom: &T) -> Self {
        let num_source = custom.get_params().k;
        Self::new_with_num_source(custom, num_source)
    }

    /// Creates a decoder with application source count K (payload symbols).
    ///
    /// Internal block size remains `scheme.get_params().k` (K′). When `num_source < K′`,
    /// implicit padding is installed during solver initialization.
    pub fn new_with_num_source<T: CodeScheme>(custom: &T, num_source: usize) -> Self {
        let manager = DataManager::new();
        Self::initialize(custom, manager, num_source)
    }

    /// Creates a new decoder with a data operator for immediate operation execution.
    pub fn new_with_operator<T: CodeScheme>(custom: &T, operator: Box<dyn DataOperator>) -> Self {
        let num_source = custom.get_params().k;
        Self::new_with_operator_and_num_source(custom, operator, num_source)
    }

    /// Like [`Self::new_with_num_source`], with a data operator attached.
    pub fn new_with_operator_and_num_source<T: CodeScheme>(
        custom: &T,
        operator: Box<dyn DataOperator>,
        num_source: usize,
    ) -> Self {
        let manager = DataManager::new_with_operator(operator);
        Self::initialize(custom, manager, num_source)
    }

    /// Like [`Self::new_with_operator`], but does not record operations (execute-only).
    pub fn new_with_operator_execute_only<T: CodeScheme>(
        custom: &T,
        operator: Box<dyn DataOperator>,
    ) -> Self {
        let num_source = custom.get_params().k;
        Self::new_with_operator_execute_only_and_num_source(custom, operator, num_source)
    }

    /// Like [`Self::new_with_operator_and_num_source`], but skips recording operations.
    pub fn new_with_operator_execute_only_and_num_source<T: CodeScheme>(
        custom: &T,
        operator: Box<dyn DataOperator>,
        num_source: usize,
    ) -> Self {
        let manager = DataManager::new_with_operator_execute_only(operator);
        Self::initialize(custom, manager, num_source)
    }

    fn initialize<T: CodeScheme>(custom: &T, mut manager: DataManager, num_source: usize) -> Self {
        let params = custom.get_params();
        let mut data_ids = Option::None;
        let mut gen_degree_set = Option::None;
        let code_type = custom.code_type();
        let solver_type = match code_type {
            CodeType::Systematic => SolverType::SysDec,
            CodeType::Ordinary => SolverType::OrdDec,
        };
        manager.config_from(params.clone(), solver_type);
        manager.set_num_source(num_source);

        if code_type == CodeType::Systematic {
            //dbg!("systematic decoding");
            data_ids = Some(vec![false; num_source]);
            gen_degree_set = Some(custom.create_degree_set_fn());
        } else {
            //dbg!("ordinary decoding");
        }
        let solver = SystemSolver::new(custom, &mut manager, solver_type);
        Self {
            params,
            manager,
            solver,
            received_msg_vectors: data_ids,
            gen_degree_set,
        }
    }

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

    /// Returns the current decode status, accounting for systematic direct reception.
    pub fn decode_status(&self) -> DecodeStatus {
        if let Some(received_msg_vectors) = self.received_msg_vectors.as_ref() {
            let num_received = received_msg_vectors.iter().filter(|&id| *id).count();
            if num_received == self.params.num_message() {
                return DecodeStatus::Decoded;
            }
        }
        self.solver.status
    }

    /// Inserts a coded vector's data and feeds it to the solver. Returns the updated decode status.
    pub fn add_coded_vector(&mut self, coded_id: usize, vector: &[u8]) -> DecodeStatus {
        self.manager.insert_coded_vector(coded_id, vector);
        self.add_coded_id(coded_id)
    }

    /// Feeds a previously-inserted coded vector ID to the solver. Returns the updated decode status.
    pub fn add_coded_id(&mut self, coded_id: usize) -> DecodeStatus {
        let data_id = self.manager.insert_coded_id(coded_id);
        self.manager.add_coded_vector(coded_id, data_id);

        //let num_msg = self.params.num_message();
        let num_msg = self.manager.num_source();
        if coded_id < num_msg {
            // systematic decoding
            if let Some(received_msg_vectors) = self.received_msg_vectors.as_mut() {
                // copy the coded vector to the message vector
                //let msg_data_id = manager.data_id_of_message_vector(coded_id);
                self.manager.copy_to(data_id, coded_id);
                // check the number of received message vectors
                received_msg_vectors[coded_id] = true;
                let num_received = received_msg_vectors.iter().filter(|&id| *id).count();
                if num_received == num_msg {
                    return DecodeStatus::Decoded;
                }
            } else {
                eprintln!(
                    "coded_id: {} is only supported for systematic decoding",
                    coded_id
                );
                return self.solver.status;
            }
        } else if coded_id < self.params.num_total() {
            eprintln!("coded_id: {} is for parity-check constraints", coded_id);
            return self.solver.status;
        }

        self.solver
            .add_coded_vector(&mut self.manager, coded_id, data_id);

        if self.solver.status == DecodeStatus::Decoded {
            if let Some(received_msg_vectors) = self.received_msg_vectors.as_mut() {
                let gen_degree_set = self.gen_degree_set.as_mut().unwrap();
                // generate the missing message vectors
                for (msg_id, received) in received_msg_vectors
                    .iter_mut()
                    .enumerate()
                    .take(self.params.num_message())
                {
                    if *received {
                        continue;
                    }
                    let degree_set = (gen_degree_set)(msg_id);
                    let data_ids = degree_set
                        .iter()
                        .map(|&id| self.manager.data_id_of_variable_vector(id))
                        .collect::<Vec<_>>();
                    self.manager.ensure_zero_one(msg_id);
                    self.manager.add_to_vector_owned(data_ids, msg_id);
                    *received = true;
                }
            } else {
                // ordinary decoding
                self.manager.restore_for_ordinary();
            }
        }

        self.solver.status
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::traits::{CodeScheme, PrecodePair};
    use crate::types::{DecodingConfig, DegreeSetFn};

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

    #[test]
    fn new_with_num_source_registers_systematic_padding_coded_ids() {
        let scheme = DummySystematic { k: 12 };
        let decoder = Decoder::new_with_num_source(&scheme, 10);
        assert_eq!(decoder.manager.num_source(), 10);
        assert!(decoder.manager.has_padding());
        assert!(decoder.manager.data_id_of_coded_vector(10).is_some());
        assert!(decoder.manager.data_id_of_coded_vector(11).is_some());
        assert!(decoder.manager.data_id_of_coded_vector(12).is_none());
    }

    #[test]
    fn new_without_num_source_has_no_padding_bootstrap() {
        let scheme = DummySystematic { k: 12 };
        let decoder = Decoder::new(&scheme);
        assert!(!decoder.manager.has_padding());
        assert!(decoder.manager.data_id_of_coded_vector(10).is_none());
    }
}
