// Copyright (c) 2025 Shenghao Yang.
// All rights reserved.

use crate::core::solver::Solver;
use crate::data_manager::DataManager;
use crate::traits::{CodeScheme, DataOperator};
use crate::types::{CodeParams, CodeType, DecodeStatus, DegreeSetFn, SolverType};

pub struct Decoder {
    params: CodeParams,
    pub manager: DataManager,
    solver: Solver,
    received_msg_vectors: Option<Vec<bool>>,
    gen_degree_set: Option<DegreeSetFn>,
}

impl Decoder {
    pub fn new<T: CodeScheme>(custom: T) -> Self {
        let manager = DataManager::new();
        Self::initialize(custom, manager)
    }

    pub fn new_with_operator<T: CodeScheme>(custom: T, operator: Box<dyn DataOperator>) -> Self {
        let manager = DataManager::new_with_operator(operator);
        Self::initialize(custom, manager)
    }

    fn initialize<T: CodeScheme>(custom: T, mut manager: DataManager) -> Self {
        let params = custom.get_params();
        let mut data_ids = Option::None;
        let mut gen_degree_set = Option::None;
        let code_type = custom.code_type();
        let solver_type = match code_type {
            CodeType::Systematic => SolverType::SysDec,
            CodeType::Ordinary => SolverType::OrdDec,
        };
        manager.config_from(params.clone(), solver_type);
        if code_type == CodeType::Systematic {
            //dbg!("systematic decoding");
            data_ids = Some(vec![false; params.k]);
            gen_degree_set = Some(custom.create_degree_set_fn());
        } else {
            //dbg!("ordinary decoding");
        }
        let solver = Solver::new(&custom, &mut manager);
        Self {
            params,
            manager,
            solver,
            received_msg_vectors: data_ids,
            gen_degree_set,
        }
    }

    pub fn decode_status(&self) -> DecodeStatus {
        if let Some(received_msg_vectors) = self.received_msg_vectors.as_ref() {
            let num_received = received_msg_vectors.iter().filter(|&id| *id).count();
            if num_received == self.params.k {
                return DecodeStatus::Decoded;
            }
        }
        self.solver.status
    }

    pub fn add_coded_vector(&mut self, coded_id: usize, vector: &[u8]) -> DecodeStatus {
        self.manager.insert_coded_vector(coded_id, vector);
        self.add_coded_id(coded_id)
    }

    pub fn add_coded_id(&mut self, coded_id: usize) -> DecodeStatus {
        let data_id = self.manager.insert_coded_id(coded_id);
        self.manager.add_coded_vector(coded_id, data_id);

        if coded_id < self.params.k {
            // systematic decoding
            if let Some(received_msg_vectors) = self.received_msg_vectors.as_mut() {
                // copy the coded vector to the message vector
                //let msg_data_id = manager.data_id_of_message_vector(coded_id);
                self.manager.copy_to(data_id, coded_id);
                // check the number of received message vectors
                received_msg_vectors[coded_id] = true;
                let num_received = received_msg_vectors.iter().filter(|&id| *id).count();
                if num_received == self.params.k {
                    return DecodeStatus::Decoded;
                }
            } else {
                dbg!(
                    "coded_id: {} is only supported for systematic decoding",
                    coded_id
                );
                return self.solver.status;
            }
        } else if coded_id < self.params.num_total() {
            dbg!(
                "coded_id: {} is only for parity-check constraints",
                coded_id
            );
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
                    .take(self.params.k)
                {
                    if *received {
                        continue;
                    }
                    let degree_set = (gen_degree_set)(msg_id);
                    let data_ids = degree_set
                        .iter()
                        .map(|&id| self.manager.data_id_of_variable_vector(id))
                        .collect::<Vec<_>>();
                    self.manager.ensure_zero(&[msg_id]);
                    self.manager.add_to_vector(&data_ids, msg_id);
                    *received = true;
                }
            } else {
                // ordinary decoding
                for i in 0..self.params.b {
                    self.manager.move_to(
                        self.manager.data_id_of_inactive_variable(i),
                        i + self.params.a,
                    );
                }
            }
        }

        self.solver.status
    }
}
