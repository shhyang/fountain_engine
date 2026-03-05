// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use crate::traits::DataOperator;
//se crate::data_operators::VecDataOperater;
use crate::algebra::finite_field::GF256;
use crate::types::{CodeParams, SolverType, Operation};
use std::collections::HashMap;

/// Starting point for dynamically-allocated data vector IDs.
const MAX_DATA_VECTOR_ID: usize = 10000000;


/// Interface between the encoder/decoder and the underlying data storage.
///
/// `DataManager` assigns and tracks unique data vector IDs, records
/// [`Operation`]s for delayed or immediate execution, and optionally
/// delegates to a [`DataOperator`] for actual vector storage and arithmetic.
pub struct DataManager {
    /// Code parameters defining the vector layout.
    params: CodeParams,
    /// Base offset for variable (intermediate) data vector IDs.
    variable_data_id_0: usize,
    /// Chronological log of all recorded operations.
    operations: Vec<Operation>,
    /// Index into `operations` up to which operations have been retrieved.
    last_retrieved_index: usize,
    //temp_data_id: usize,
    /// Next available data vector ID for dynamic allocation.
    next_data_id: usize,
    /// Maps external coded vector IDs to internal data vector IDs.
    coded_id_to_data_id: HashMap<usize, usize>,
    /// Optional operator that executes operations on actual data.
    operator: Option<Box<dyn DataOperator>>,
    /// Count of coded vectors inserted so far.
    pub coded_vector_inserted: usize,
}

impl Default for DataManager {
    fn default() -> Self {
        Self::new()
    }
}

impl DataManager {
    /// Create a new DataManager
    pub fn new() -> Self {
        Self {
            params: CodeParams::new(0, 0, 0, 0),
            variable_data_id_0: 0,
            operations: Vec::new(),
            //temp_data_id: 0,
            next_data_id: MAX_DATA_VECTOR_ID,
            coded_id_to_data_id: HashMap::new(),
            operator: None,
            last_retrieved_index: 0,
            coded_vector_inserted: 0,
        }
    }

    /// Creates a new `DataManager` with the given data operator for immediate execution.
    pub fn new_with_operator(operator: Box<dyn DataOperator>) -> Self {
        Self {
            params: CodeParams::new(0, 0, 0, 0),
            variable_data_id_0: 0,
            operations: Vec::new(),
            //temp_data_id: 0,
            next_data_id: MAX_DATA_VECTOR_ID,
            coded_id_to_data_id: HashMap::new(),
            operator: Some(operator),
            last_retrieved_index: 0,
            coded_vector_inserted: 0,
        }
    }

    /// Configures the manager with code parameters and solver type, setting up ID ranges.
    pub fn config_from(&mut self, params: CodeParams, solver_type: SolverType) {
        self.params = params;
        match solver_type {
            SolverType::OrdEnc => {
                self.variable_data_id_0 = 0;
            }
            SolverType::OrdDec => {
                self.variable_data_id_0 = 0;
            }
            SolverType::SysEnc => {
                self.variable_data_id_0 = self.params.k;
            }
            SolverType::SysDec => {
                self.variable_data_id_0 = self.params.k;
            }
        }
        self.next_data_id = self.variable_data_id_0 + self.params.num_total();
        //self.temp_data_id = self.next_data_id;
        //self.next_data_id += 1;
    }

    /// Registers a coded vector ID, allocating a new data ID if not already mapped. Returns the data ID.
    pub fn insert_coded_id(&mut self, coded_id: usize) -> usize {
        self.coded_vector_inserted += 1;
        if let Some(data_id) = self.coded_id_to_data_id.get(&coded_id) {
            *data_id
        } else {
            let data_id = self.next_data_id;
            self.next_data_id += 1;
            self.coded_id_to_data_id.insert(coded_id, data_id);
            data_id
        }
    }

    /// Explicitly maps a coded vector ID to an existing data vector ID.
    pub fn assign_data_id(&mut self, coded_id: usize, data_id: usize) {
        self.coded_id_to_data_id.insert(coded_id, data_id);
    }

    /// Allocates a fresh zero-initialized data vector for a coded vector and returns its data ID.
    pub fn coded_data_id(&mut self, coded_id: usize) -> usize {
        let data_id = self.next_data_id;
        self.next_data_id += 1;
        self.coded_id_to_data_id.insert(coded_id, data_id);
        self.ensure_zero(&[data_id]);
        data_id
    }

    /// Allocates and returns a fresh zero-initialized temporary data vector ID.
    pub fn temp_data_id(&mut self) -> usize {
        let data_id = self.next_data_id;
        self.next_data_id += 1;
        self.ensure_zero(&[data_id]);
        data_id
    }

    /// Looks up the data vector ID for a previously inserted coded vector, if any.
    pub fn data_id_of_coded_vector(&self, coded_id: usize) -> Option<usize> {
        self.coded_id_to_data_id.get(&coded_id).cloned()
    }

    /*
    pub fn get_data_id(&self, coded_id: usize) -> Option<usize> {
        self.coded_id_to_data_id.get(&coded_id).cloned()
    }

    /// Check if a coded vector has been allocated
    pub fn has_coded_vector(&self, coded_id: usize) -> bool {
        self.coded_id_to_data_id.contains_key(&coded_id)
    }

    /// Get all coded ID to data ID mappings
    pub fn get_coded_vector_mappings(&self) -> &HashMap<usize, usize> {
        &self.coded_id_to_data_id
    }
    */

    /// Returns the data vector ID for a variable vector at the given index.
    pub fn data_id_of_variable_vector(&self, var_id: usize) -> usize {
        var_id + self.variable_data_id_0
    }

    /// Returns data IDs for all active variable vectors (source active + LDPC).
    pub fn data_id_range_of_active_variable(&self) -> Vec<usize> {
        (self.variable_data_id_0..self.variable_data_id_0 + self.params.num_active()).collect()
    }

    /// Returns data IDs for all message and LDPC variable vectors.
    pub fn data_id_range_of_msg_ldpc_variable(&self) -> Vec<usize> {
        (self.variable_data_id_0..self.variable_data_id_0 + self.params.num_message_ldpc())
            .collect()
    }

    /// Returns data IDs for all LDPC variable vectors.
    pub fn data_id_range_of_ldpc_variable(&self) -> Vec<usize> {
        (self.variable_data_id_0 + self.params.a
            ..self.variable_data_id_0 + self.params.a + self.params.l)
            .collect()
    }

    /// Returns data IDs for all HDPC variable vectors.
    pub fn data_id_range_of_hdpc_variable(&self) -> Vec<usize> {
        (self.variable_data_id_0 + self.params.num_message_ldpc()
            ..self.variable_data_id_0 + self.params.num_total())
            .collect()
    }

    /// Returns the data vector ID for the active variable at the given index.
    pub fn data_id_of_active_variable(&self, idx: usize) -> usize {
        self.variable_data_id_0 + idx
    }

    /// Returns the data vector ID for the inactive variable at the given index.
    pub fn data_id_of_inactive_variable(&self, idx: usize) -> usize {
        self.variable_data_id_0 + self.params.num_active() + idx
    }

    /// Get the data vector ID for the LDPC variable vector at index idx
    pub fn data_id_of_ldpc_variable(&mut self, idx: usize) -> usize {
        self.variable_data_id_0 + self.params.a + idx
    }

    
    /// Returns the data vector ID for the HDPC variable at the given index.
    pub fn data_id_of_hdpc_variable(&self, idx: usize) -> usize {
        self.variable_data_id_0 + idx + self.params.num_message_ldpc()
    }
    

    /// Get all stored operations
    pub fn get_operations(&self) -> &[Operation] {
        &self.operations
    }

    /// Clear all stored operations
    pub fn clear_operations(&mut self) {
        self.operations.clear();
        self.last_retrieved_index = 0;
    }

    /// Get new operations
    pub fn move_new_operations(&mut self) -> Vec<Operation> {
        let start = self.last_retrieved_index;
        let end = self.operations.len();
        self.last_retrieved_index = end;
        self.operations[start..end].to_vec()
    }

    /// Takes the data operator out of this manager, returning ownership to the caller.
    pub fn move_operator(&mut self) -> Box<dyn DataOperator> {
        self.operator.take().unwrap()
    }

    /// Installs a data operator into this manager for immediate operation execution.
    pub fn set_operator(&mut self, operator: Box<dyn DataOperator>) {
        self.operator = Some(operator);
    }

    fn save_operation(&mut self, operation: Operation) {
        self.operations.push(operation.clone());
        if let Some(ref mut operator) = self.operator {
            operator.execute(&operation);
            //self.last_retrieved_index = self.operations.len();
        }
    }

    /// Returns a copy of the variable vector at the given variable index. Requires an operator.
    pub fn get_variable_vector(&self, var_id: usize) -> Vec<u8> {
        if let Some(operator) = self.operator.as_ref() {
            operator
                .get_vector(self.data_id_of_variable_vector(var_id))
                .to_vec()
        } else {
            panic!("Operator is not set");
        }
    }

    /// Returns a reference to the data vector with the given ID. Requires an operator.
    pub fn get_data_vector(&self, data_id: usize) -> &[u8] {
        if let Some(operator) = self.operator.as_ref() {
            operator.get_vector(data_id)
        } else {
            panic!("Operator is not set");
        }
    }

    /// Returns a copy of the coded vector with the given coded ID. Requires an operator.
    pub fn get_coded_vector(&self, coded_id: usize) -> Vec<u8> {
        if let Some(data_id) = self.coded_id_to_data_id.get(&coded_id) {
            if let Some(operator) = self.operator.as_ref() {
                operator.get_vector(*data_id).to_vec()
            } else {
                panic!("Operator is not set");
            }
        } else {
            panic!("Coded vector with ID {} does not exist", coded_id);
        }
    }

    /// Inserts raw data into the operator at the given data ID. Used for testing.
    pub fn insert_data_vector(&mut self, data_id: usize, vector: &[u8]) {
        if let Some(operator) = self.operator.as_mut() {
            operator.insert_vector(vector, data_id);
        } else {
            panic!("Operator is not set");
        }
    }

    /// Inserts a coded vector by coded ID, allocating a data ID and storing via the operator. Used for testing.
    pub fn insert_coded_vector(&mut self, coded_id: usize, vector: &[u8]) {
        let data_id = self.insert_coded_id(coded_id);
        if let Some(operator) = self.operator.as_mut() {
            operator.insert_vector(vector, data_id);
        } else {
            panic!("Operator is not set");
        }
    }

    //pub fn permute_vectors(&mut self, ids: &[usize], perm: &[usize]) {
    //    if ids.len() != perm.len() {
    //        panic!("The number of IDs must be equal to the number of permutations");
    //    }

    //    self.operations.push(Operation::PermuteVectors {
    //        ids: ids.to_vec(),
    //        perm: perm.to_vec(),
    //    });
    //}

    /// Records an [`Operation::EnsureZero`] to zero the given vectors.
    pub fn ensure_zero(&mut self, list_id: &[usize]) {
        self.save_operation(Operation::EnsureZero {
            list_id: list_id.to_vec(),
        });
    }

    /// Records an [`Operation::MultiplyAlpha`] to multiply a vector by the primitive element.
    pub fn multiply_alpha(&mut self, id: usize) {
        self.save_operation(Operation::MultiplyAlpha { id });
    }

    /// Multiplies a vector by a GF(256) scalar, with special-case optimizations for 0, 1, and alpha.
    pub fn multiply_scalar(&mut self, scalar: u8, id: usize) {
        if scalar == 0 {
            self.save_operation(Operation::EnsureZero { list_id: vec![id] });
        } else if scalar == 1 {
        } else if scalar == GF256::default().alpha {
            self.save_operation(Operation::MultiplyAlpha { id });
        } else {
            self.save_operation(Operation::MultiplyScalar { scalar, id });
        }
    }

    /// Divides a vector by a GF(256) scalar (multiplies by its inverse).
    pub fn divide_scalar(&mut self, scalar: u8, id: usize) {
        if scalar != 1 {
            let gf = GF256::default();
            let inverse = gf.inverse(scalar);
            if inverse == gf.alpha {
                self.save_operation(Operation::MultiplyAlpha { id });
            } else {
                self.save_operation(Operation::MultiplyScalar {
                    scalar: inverse,
                    id,
                });
            }
        }
    }

    /// XORs (adds) multiple source vectors into a single target vector.
    pub fn add_to_vector(&mut self, list_id: &[usize], target_id: usize) {
        self.save_operation(Operation::AddToVector {
            list_id: list_id.to_vec(),
            target_id,
        });
    }

    /// XORs a single source vector into each of the given target vectors.
    pub fn broadcast_add(&mut self, src_id: usize, target_ids: &[usize]) {
        self.save_operation(Operation::BroadcastAdd {
            src_id,
            target_ids: target_ids.to_vec(),
        });
    }

    /// Computes `target += scalar * src` over GF(256), with fast-paths for scalar 0 and 1.
    pub fn mul_add(&mut self, src_id: usize, scalar: u8, target_id: usize) {
        if scalar == 0 {
        } else if scalar == 1 {
            self.save_operation(Operation::AddToVector {
                list_id: vec![src_id],
                target_id,
            });
        } else {
            self.save_operation(Operation::MulAdd {
                src_id,
                scalar,
                target_id,
            });
        }
    }

    /// Moves a vector from `src_id` to `target_id`. No-op if they are equal.
    pub fn move_to(&mut self, src_id: usize, target_id: usize) {
        if src_id == target_id {
            return;
        }
        self.save_operation(Operation::MoveTo { src_id, target_id });
    }

    /// Copies a vector from `src_id` to `target_id`. No-op if they are equal.
    pub fn copy_to(&mut self, src_id: usize, target_id: usize) {
        if src_id == target_id {
            return;
        }
        self.save_operation(Operation::CopyTo { src_id, target_id });
    }

    /// Records a [`Operation::Remove`] to deallocate a vector.
    pub fn remove(&mut self, id: usize) {
        self.save_operation(Operation::Remove { id });
    }

    /// Records an informational marker associating a coded vector with its data vector (for decoding).
    pub fn add_coded_vector(&mut self, coded_id: usize, data_id: usize) {
        self.save_operation(Operation::InfoCodedVector { coded_id, data_id });
    }

    /// Records an informational marker associating a coded vector with its data vector (for encoding).
    pub fn encode_coded_vector(&mut self, coded_id: usize, data_id: usize) {
        self.save_operation(Operation::InfoCodedVector { coded_id, data_id });
    }
}

/// Macro for LU solve operations that works with any DataManager implementation
/// This macro implements forward and backward substitution for solving Ax = b in-place
/// where A is an LU-decomposed matrix and b are the target vectors.
#[macro_export]
macro_rules! lu_solve {
    ($manager:expr, $matrix_a:expr, $target_ids:expr) => {{
        let matrix_a = $matrix_a;
        let target_ids = $target_ids;

        if matrix_a.len() != target_ids.len() {
            panic!("The number of rows in A must be equal to the number of target IDs");
        }

        let n = matrix_a.len();

        // Forward substitution (L part): solve Ly = b
        for j in 0..n - 1 {
            for i in j + 1..n {
                if matrix_a[i][j] != 0 {
                    $manager.mul_add(target_ids[j], matrix_a[i][j], target_ids[i]);
                }
            }
        }

        // Backward substitution (U part): solve Ux = y
        for j in (0..n).rev() {
            if matrix_a[j][j] == 0 {
                panic!(
                    "Singular matrix: diagonal element at position {} is zero",
                    j
                );
            }

            // Scale the current row by the inverse of the diagonal element
            $manager.divide_scalar(matrix_a[j][j], target_ids[j]);

            // Subtract scaled row from previous rows
            for i in (0..j).rev() {
                if matrix_a[i][j] != 0 {
                    $manager.mul_add(target_ids[j], matrix_a[i][j], target_ids[i]);
                }
            }
        }
    }};
}

/// Incremental LU-solve macro using a column permutation vector `q`.
///
/// Like [`lu_solve!`] but accesses columns through `q[j]` to support
/// the incremental LU decomposition produced by [`Matrix::lu_decomp_incr`](crate::algebra::linear_algebra::Matrix::lu_decomp_incr).
#[macro_export]
macro_rules! lu_solve_incr {
    ($manager:expr, $matrix_a:expr, $target_ids:expr, $q:expr) => {{
        let matrix_a = $matrix_a;
        let target_ids = $target_ids;
        let q = $q;

        if matrix_a.len() != target_ids.len() {
            panic!("The number of rows in A must be equal to the number of target IDs");
        }

        let n = matrix_a.len();

        // Forward substitution (L part): solve Ly = b
        for j in 0..n - 1 {
            for i in j + 1..n {
                let l = matrix_a[i][q[j]];
                if l != 0 {
                    $manager.mul_add(target_ids[j], l, target_ids[i]);
                }
            }
        }

        // Backward substitution (U part): solve Ux = y
        for j in (0..n).rev() {
            let l = matrix_a[j][q[j]];
            if l == 0 {
                panic!(
                    "Singular matrix: diagonal element at position {} is zero",
                    j
                );
            }

            // Scale the current row by the inverse of the diagonal element
            $manager.divide_scalar(l, target_ids[j]);

            // Subtract scaled row from previous rows
            for i in (0..j).rev() {
                let l = matrix_a[i][q[j]];
                if l != 0 {
                    $manager.mul_add(target_ids[j], l, target_ids[i]);
                }
            }
        }
    }};
}
