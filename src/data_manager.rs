// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use crate::traits::DataOperator;
//se crate::data_operators::VecDataOperater;
use crate::algebra::finite_field::{Field, GF256};
use crate::types::{CodeParams, GF2_FIELD_POLY, Operation, SolverType};
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
    /// When `false`, operations are executed on the operator but not appended to `operations`.
    /// Use for production decode with an attached operator when trace/replay is not needed.
    record_operations: bool,
    /// Count of coded vectors inserted so far.
    pub num_coded_vector_inserted: usize,
    /// GF(256) instance for this session; unset when configured with [`GF2_FIELD_POLY`].
    gf256: Option<GF256>,
    /// Application source count K (payload symbols). Defaults to [`CodeParams::k`] in [`Self::config_from`].
    num_source: usize,
    /// Session direction from the last [`Self::config_from`] call (encode vs decode).
    solver_type: SolverType,
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
            num_coded_vector_inserted: 0,
            gf256: None,
            record_operations: true,
            num_source: 0,
            solver_type: SolverType::OrdEnc,
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
            num_coded_vector_inserted: 0,
            gf256: None,
            record_operations: true,
            num_source: 0,
            solver_type: SolverType::OrdEnc,
        }
    }

    /// Like [`Self::new_with_operator`], but skips recording operations (execute-only).
    ///
    /// Suitable when an operator is attached and operation traces are not needed.
    pub fn new_with_operator_execute_only(operator: Box<dyn DataOperator>) -> Self {
        let mut manager = Self::new_with_operator(operator);
        manager.record_operations = false;
        manager
    }

    /// When `false`, `save_operation` still runs the operator but
    /// does not push to the operation log.
    pub fn set_record_operations(&mut self, record: bool) {
        self.record_operations = record;
    }

    /// Returns whether operations are appended to the log after execution.
    #[must_use]
    pub fn records_operations(&self) -> bool {
        self.record_operations
    }

    /// Disables operation logging; operations are still executed on the attached operator.
    pub fn execute_only(&mut self) {
        self.record_operations = false;
    }

    /// Configures the manager with code parameters and solver type, setting up ID ranges.
    ///
    pub fn config_from(&mut self, params: CodeParams, solver_type: SolverType) {
        self.params = params;
        self.solver_type = solver_type;
        self.num_source = self.params.k;
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

    /// Sets the application source block size K (payload symbol count).
    ///
    /// Must satisfy `k <= params.k` (internal block size K′). When `k < params.k`, implicit
    /// padding applies; encode/decode installation is handled by the encoder, decoder, and
    /// solver once those paths call this method.
    ///
    /// Default after [`Self::config_from`] is `num_source == params.k` (no padding).
    pub fn set_num_source(&mut self, k: usize) {
        assert!(
            k <= self.params.k,
            "num_source ({k}) must be <= block_k ({})",
            self.params.k
        );
        assert!(
            k >= self.params.b,
            "num_source ({k}) must be >= inactive message count ({})",
            self.params.b
        );
        self.num_source = k;
    }

    /// Application source count K (payload symbols).
    #[must_use]
    pub fn num_source(&self) -> usize {
        self.num_source
    }

    /// Implicit zero padding count: `params.k − num_source`.
    #[must_use]
    pub fn num_padding(&self) -> usize {
        self.params.k.saturating_sub(self.num_source)
    }

    /// Whether this session uses implicit padding (`num_source < params.k`).
    #[must_use]
    pub fn has_padding(&self) -> bool {
        self.num_padding() > 0
    }

    /// Solver type from the last [`Self::config_from`] call.
    #[must_use]
    pub fn solver_type(&self) -> SolverType {
        self.solver_type
    }

    /// Configures the finite field for LU solves and GF(256) scalar ops on this session.
    ///
    /// - `pp == [`GF2_FIELD_POLY`](crate::types::GF2_FIELD_POLY)`: GF(2) mode — `gf256()` is `None`;
    ///   matrix elimination uses XOR; binary HDPC must not request non-trivial scalars.
    /// - Any other allowed primitive polynomial (e.g. `0x11D`): builds GF(256) tables for
    ///   `multiply_scalar`, `divide_scalar`, and LU over full `u8` coefficients.
    ///
    /// Propagates `pp` to the [`DataOperator`](crate::traits::DataOperator) when attached.
    pub fn config_finite_field(&mut self, pp: u16) {
        if pp == GF2_FIELD_POLY {
            self.gf256 = None;
        } else {
            self.gf256 = Some(GF256::new_with_primitive_polynomial(pp));
        }
        if let Some(operator) = self.operator.as_mut() {
            operator.config_finite_field(pp);
        }
    }

    /// Same session field as `gf` (clones tables); propagates to the data operator.
    pub fn config_finite_field_from(&mut self, gf: &GF256) {
        self.gf256 = Some(gf.clone());
        if let Some(operator) = self.operator.as_mut() {
            operator.config_finite_field_from(gf);
        }
    }

    /// GF(256) used for scalar multiply, divide, and consistency with matrix solves on this manager.
    #[inline]
    pub fn gf256(&self) -> Option<&GF256> {
        self.gf256.as_ref()
    }

    /// Registers a coded vector ID, allocating a new data ID if not already mapped. Returns the data ID.
    pub fn insert_coded_id(&mut self, coded_id: usize) -> usize {
        self.num_coded_vector_inserted += 1;
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
        if let Some(ref mut operator) = self.operator {
            operator.execute(&operation);
        }
        if self.record_operations {
            self.operations.push(operation);
        }
    }

    fn save_add_to_vector(&mut self, list_id: Vec<usize>, target_id: usize) {
        self.save_operation(Operation::AddToVector { list_id, target_id });
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
        } else if let Some(gf) = self.gf256() {
            if scalar == gf.primitive_element() {
                self.save_operation(Operation::MultiplyAlpha { id });
            } else {
                self.save_operation(Operation::MultiplyScalar { scalar, id });
            }
        } else {
            panic!("GF(256) is not set");
        }
    }

    /// Divides a vector by a GF(256) scalar (multiplies by its inverse).
    pub fn divide_scalar(&mut self, scalar: u8, id: usize) {
        if scalar == 1 {
        } else if let Some(gf) = self.gf256() {
            let inverse = gf.inverse(scalar);
            if inverse == gf.primitive_element() {
                self.save_operation(Operation::MultiplyAlpha { id });
            } else {
                self.save_operation(Operation::MultiplyScalar {
                    scalar: inverse,
                    id,
                });
            }
        } else {
            panic!("GF(256) is not set");
        }
    }

    /// XORs (adds) multiple source vectors into a single target vector.
    pub fn add_to_vector(&mut self, list_id: &[usize], target_id: usize) {
        match list_id.len() {
            0 => {}
            1 => self.add_one_to_vector(list_id[0], target_id),
            2 => self.add_two_to_vector(list_id[0], list_id[1], target_id),
            3 => self.add_three_to_vector(list_id[0], list_id[1], list_id[2], target_id),
            _ => self.save_add_to_vector(list_id.to_vec(), target_id),
        }
    }

    /// XOR one source vector into a target (hot path for GF(2) LU / single-source BS).
    pub fn add_one_to_vector(&mut self, src_id: usize, target_id: usize) {
        self.save_operation(Operation::AddOneToVector { src_id, target_id });
    }

    /// XOR two source vectors into a target (hot path for sparse back-substitution).
    pub fn add_two_to_vector(&mut self, s0: usize, s1: usize, target_id: usize) {
        if s0 == s1 {
            return;
        }
        self.save_operation(Operation::AddTwoToVector { s0, s1, target_id });
    }

    /// XOR three source vectors into a target (hot path for sparse back-substitution).
    pub fn add_three_to_vector(&mut self, s0: usize, s1: usize, s2: usize, target_id: usize) {
        self.save_operation(Operation::AddThreeToVector {
            s0,
            s1,
            s2,
            target_id,
        });
    }

    /// XOR sources into a target, taking ownership of `list_id` (no slice copy).
    pub fn add_to_vector_owned(&mut self, list_id: Vec<usize>, target_id: usize) {
        match list_id.len() {
            0 => {}
            1 => self.add_one_to_vector(list_id[0], target_id),
            2 => self.add_two_to_vector(list_id[0], list_id[1], target_id),
            3 => self.add_three_to_vector(list_id[0], list_id[1], list_id[2], target_id),
            _ => self.save_add_to_vector(list_id, target_id),
        }
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
            self.add_one_to_vector(src_id, target_id);
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

    /*
    /// Is a padding vector?
    pub fn is_padding(&self, var_id: usize) -> bool {
        if self.params.p == 0 {
            return false;
        }
        match self.padding_config {
            PaddingConfig::AtEnd => {
                if self.params.p <= self.params.b {
                    var_id >= self.params.num_message_ldpc() - self.params.p && var_id < self.params.num_message_ldpc()
                } else {
                    (var_id >= self.params.num_message_ldpc() - self.params.b && var_id < self.params.num_message_ldpc()) || (var_id >= self.params.k - self.params.p && var_id < self.params.a)
                }
            }
            PaddingConfig::AsActive => {
                var_id >= self.params.a - self.params.p && var_id < self.params.a
            }
        }
    }*/

    /// For ordinary coding with padding, prepare the data vector IDs so that padding vectors
    /// occupy the tail of the active message range before precoding, then move inactive
    /// message vectors into inactive variable slots.
    pub fn prepare_for_ordinary(&mut self) {
        for i in 1..=self.params.b {
            self.move_to(self.num_source - i, self.data_id_of_inactive_variable(self.params.b - i));
        }
        if self.has_padding() {
            self.ensure_zero(&(self.num_source - self.params.b..self.params.a).collect::<Vec<_>>());
        }
            // move the inactive message vectors to the inactive message variable data vectors
            //for i in (self.params.a..self.params.k).rev() {
            //    self.move_to(i, self.data_id_of_inactive_variable(i - self.params.a));
            //}
       // } else {
       //     let p_start = self.params.a - self.params.p;
       //     let num_msg = self.params.num_message();
       //     for i in (p_start..num_msg).rev() {
       //         self.move_to(i, self.data_id_of_inactive_variable(i - p_start));
       //     }
       //     // append the padding vectors to the active message range
       //     self.ensure_zero(&(p_start..self.params.a).collect::<Vec<_>>());
       // }
    }

    /// Restore the data vector IDs after ordinary decoding completes.
    pub fn restore_for_ordinary(&mut self) {
        //if self.params.p == 0 {
        let num_a = self.params.a - self.num_padding();
        for i in 0..self.params.b {
            self.move_to(self.data_id_of_inactive_variable(i), i + num_a);
        }
        //} else if self.padding_config == PaddingConfig::AtEnd {
       //     if self.params.p < self.params.b {
       //         let p_end = self.params.b - self.params.p;
       //         for i in 0..p_end {
       //             self.move_to(self.data_id_of_inactive_variable(i), i + self.params.a);
       //         }
       //     }
       // } else {
       //     let p_start = self.params.a - self.params.p;
       //     for i in 0..self.params.b {
       //         self.move_to(self.data_id_of_inactive_variable(i), i + p_start);
       //     }
       // }
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

#[cfg(test)]
mod num_source_tests {
    use super::*;
    use crate::types::SolverType;

    #[test]
    fn config_from_defaults_num_source_to_block_k() {
        let mut mgr = DataManager::new();
        mgr.config_from(CodeParams::new(12, 12, 0, 0), SolverType::SysEnc);
        assert_eq!(mgr.num_source(), 12);
        assert_eq!(mgr.num_padding(), 0);
        assert!(!mgr.has_padding());
        assert_eq!(mgr.solver_type(), SolverType::SysEnc);
    }

    #[test]
    fn set_num_source_lowers_k_and_sets_padding_count() {
        let mut mgr = DataManager::new();
        mgr.config_from(CodeParams::new(12, 11, 0, 0), SolverType::OrdEnc);
        mgr.set_num_source(10);
        assert_eq!(mgr.num_source(), 10);
        assert_eq!(mgr.num_padding(), 2);
        assert!(mgr.has_padding());
    }

    #[test]
    #[should_panic(expected = "num_source (13) must be <= block_k (12)")]
    fn set_num_source_panics_when_k_exceeds_block_k() {
        let mut mgr = DataManager::new();
        mgr.config_from(CodeParams::new(12, 12, 0, 0), SolverType::SysDec);
        mgr.set_num_source(13);
    }
}
