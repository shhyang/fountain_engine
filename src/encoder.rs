// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use crate::data_manager::DataManager;
//use crate::types::Operation;
use crate::core::ordinary_precode_encode;
use crate::core::system_solver::SystemSolver;
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
    /// Create a new encoder and run precoding immediately (same as before).
    pub fn new<T: CodeScheme>(custom: &T) -> Self {
        let num_source = custom.get_params().k;
        Self::new_with_num_source(custom, num_source)
    }

    /// Creates an encoder with application source count K (payload symbols).
    ///
    /// Internal block size remains `scheme.get_params().k` (K′). When `num_source < K′`,
    /// padding is installed during [`Self::precode_encode`].
    pub fn new_with_num_source<T: CodeScheme>(custom: &T, num_source: usize) -> Self {
        let manager = DataManager::new();
        Self::initialize(custom, manager, num_source)
    }

    /// Configure the encoder and data manager only; does **not** run precoding.
    pub fn new_without_precoding<T: CodeScheme>(custom: &T) -> Self {
        let num_source = custom.get_params().k;
        Self::new_without_precoding_with_num_source(custom, num_source)
    }

    /// Like [`Self::new_with_num_source`], but leaves precoding to [`Self::precode_encode`].
    pub fn new_without_precoding_with_num_source<T: CodeScheme>(
        custom: &T,
        num_source: usize,
    ) -> Self {
        let manager = DataManager::new();
        Self::initialize_without_precoding(custom, manager, num_source)
    }

    /// Creates a new encoder with a data operator and runs precoding immediately.
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

    /// Like [`Self::new_with_operator`], but leaves precoding to [`Self::precode_encode`].
    pub fn new_with_operator_without_precoding<T: CodeScheme>(
        custom: &T,
        operator: Box<dyn DataOperator>,
    ) -> Self {
        let num_source = custom.get_params().k;
        Self::new_with_operator_without_precoding_with_num_source(custom, operator, num_source)
    }

    /// Like [`Self::new_without_precoding_with_num_source`], with a data operator attached.
    pub fn new_with_operator_without_precoding_with_num_source<T: CodeScheme>(
        custom: &T,
        operator: Box<dyn DataOperator>,
        num_source: usize,
    ) -> Self {
        let manager = DataManager::new_with_operator(operator);
        Self::initialize_without_precoding(custom, manager, num_source)
    }

    fn initialize<T: CodeScheme>(custom: &T, manager: DataManager, num_source: usize) -> Self {
        let mut enc = Self::initialize_without_precoding(custom, manager, num_source);
        enc.precode_encode(custom);
        enc
    }

    fn initialize_without_precoding<T: CodeScheme>(
        custom: &T,
        mut manager: DataManager,
        num_source: usize,
    ) -> Self {
        let params = custom.get_params();
        let gen_degree_set = custom.create_degree_set_fn();
        let code_type = custom.code_type();
        let solver_type = match code_type {
            CodeType::Systematic => SolverType::SysEnc,
            CodeType::Ordinary => SolverType::OrdEnc,
        };
        manager.config_from(params.clone(), solver_type);
        manager.set_num_source(num_source);

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
                //dbg!("ordinary encode");
                if self.params.has_precode() {
                    self.manager.prepare_for_ordinary();
                    //dbg!("precode encode");
                    ordinary_precode_encode(&mut self.manager, custom);
                }
            }
            CodeType::Systematic => {
                //dbg!("systematic encoding");
                let solver = SystemSolver::new(custom, &mut self.manager, SolverType::SysEnc);
                // solve active vectors
                //for coded_id in 0..params.a {
                //    let new_data_id = manager.data_id_of_active_variable(coded_id);
                //    manager.copy_to(coded_id, new_data_id);
                //    solver.add_coded_vector(manager, coded_id, new_data_id);
                //}
                //assert!(solver.phase == DecodePhase::GE);
                // solve inactive vectors

                if solver.status == DecodeStatus::NotDecoded {
                    panic!("systematic encoding failed");
                }

                for coded_id in 0..self.params.k {
                    self.manager.assign_data_id(coded_id, coded_id);
                }
            }
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

        //todo: fix the following logic for returning the parity-check matrix part.
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

        // coded_id >= self.params.num_total()

        let data_id = self.manager.coded_data_id(coded_id);
        //dbg!("encode_coded_vector", &coded_id, &data_id);
        let degree_set = (self.gen_degree_set)(coded_id);
        let data_ids = degree_set
            .iter()
            .map(|&id| self.manager.data_id_of_variable_vector(id))
            .collect::<Vec<_>>();
        //dbg!("encode_coded_vector", &coded_id, &active_indices, &inactive_indices);
        //manager.ensure_zero(&[data_id]);
        self.manager.add_to_vector_owned(data_ids, data_id);
        self.manager.encode_coded_vector(coded_id, data_id);
        //dbg!("encode_coded_vector", manager.get_vector(data_id));
        Some(data_id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::traits::{CodeScheme, PrecodePair};
    use crate::types::{DecodingConfig, DegreeSetFn, Operation};

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

    #[derive(Clone)]
    struct DummyOrdinary {
        k: usize,
        a: usize,
        b: usize,
    }

    impl CodeScheme for DummyOrdinary {
        fn get_params(&self) -> CodeParams {
            CodeParams::new(self.k, self.a, 1, self.b)
        }

        fn code_type(&self) -> CodeType {
            CodeType::Ordinary
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

    fn ensure_zero_ids(enc: &Encoder) -> Option<Vec<usize>> {
        enc.manager.get_operations().iter().find_map(|op| match op {
            Operation::EnsureZero { list_id } => Some(list_id.clone()),
            _ => None,
        })
    }

    #[test]
    fn new_with_num_source_sets_padding_count() {
        let scheme = DummySystematic { k: 12 };
        let enc = Encoder::new_without_precoding_with_num_source(&scheme, 10);
        assert_eq!(enc.manager.num_source(), 10);
        assert!(enc.manager.has_padding());
    }

    #[test]
    fn new_without_num_source_has_no_padding() {
        let scheme = DummySystematic { k: 12 };
        let enc = Encoder::new_without_precoding(&scheme);
        assert_eq!(enc.manager.num_source(), 12);
        assert!(!enc.manager.has_padding());
    }

    #[test]
    fn ordinary_precode_prepares_padding_slots() {
        let scheme = DummyOrdinary { k: 12, a: 11, b: 1 };
        let mut enc = Encoder::new_without_precoding_with_num_source(&scheme, 10);
        enc.precode_encode(&scheme);
        assert_eq!(ensure_zero_ids(&enc), Some(vec![9, 10]));
    }
}
