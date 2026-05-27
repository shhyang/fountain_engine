// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use crate::types::Operation;

/// Traits for data operations used by encoder and decoder.
/// They hide the implementation details of the vector operations,  
/// including the data storage and caculations.
/// They also provide data security by hiding the data operations.
/// The encoder and decoder use indices to identify the vectors.
/// The message vectors are indexed from 0 to k-1, which are supposed consecutive.
/// The data vectors are indexed by a unique vector ID, which is used by both the encoder and decoder.
/// In our implementation using Rust, usize is naturally used as the type of vector IDs.
/* pub trait DataOperations {
    /// Sets multiple packets to zero.
    fn ensure_zero(&mut self, list_id: &[usize]);

    /// Multiplies a vector by the alpha value.
    fn multiply_alpha(&mut self, id: usize);

    /// Multiplies a vector by a scalar value.
    fn multiply_scalar(&mut self, scalar: u8, id: usize);

    /// Divides a vector by a scalar value.
    fn divide_scalar(&mut self, scalar: u8, id: usize);

    /// Adds multiple vectors (from `list_id`) to a target vector (`target_id`).
    /// Equivalent to `add_many_to_one(list_id, id)` in the spec.
    fn add_to_vector(&mut self, list_id: &[usize], target_id: usize);

    /// Adds a source vector (`source_id`) to multiple target vectors (`target_ids`).
    /// Equivalent to `add_one_to_many(id, list_id)` in the spec.
    fn broadcast_add(&mut self, source_id: usize, target_ids: &[usize]);

    /// Adds a scalar multiple of a source vector (`source_id`) to multiple target vectors (`target_ids`).
    fn mul_add(&mut self, source_id: usize, scalar: u8, target_id: usize);

    /// Move a coded vector to a variable vector.
    fn move_to(&mut self, cod_id: usize, var_id: usize);

    /// Copy a vector to a target vector.
    fn copy_to(&mut self, src_id: usize, target_id: usize);

    /// Remove a vector.
    fn remove(&mut self, id: usize);
}
*/
/// This trait defines the interface for data operators used to store and retrieve the vectors.
/// It is used by the data manager to store and retrieve the vectors for the purpose of testing with data operator only.
/// This interface is not used by the encoder and decoder and is not required by the data manager in the future.
/// The data operator should be implemented for specific applications. Vector data management, vector operation acceleration, multi-threading, etc. should be considered in the data operator implementation, but not in the coding library.
pub trait DataOperator {
    /// Stores a byte vector under the given data ID.
    fn insert_vector(&mut self, _vector: &[u8], _data_id: usize) {}

    /// Retrieves a reference to the byte vector stored at the given data ID.
    fn get_vector(&self, _data_id: usize) -> &[u8] {
        &[]
    }

    /*
    fn insert_vector_owned(&mut self, vector: Vec<u8>, data_id: usize) {
        // do nothing by default
    }
    fn get_vector_owned(&mut self, data_id: usize) -> Vec<u8> {
        // return an empty vector by default
        vec![]
    }
    */

    /// Executes a single [`Operation`] on the stored vectors.
    fn execute(&mut self, operation: &Operation);

    /// Default: no-op. In-memory operators should rebuild GF(256) tables when `pp` is a real
    /// primitive polynomial; ignore [`crate::types::GF2_FIELD_POLY`].
    fn config_finite_field(&mut self, pp: u16) {
        let _ = pp;
    }

    /// Use the same GF(256) tables as `gf` (default: [`config_finite_field`](Self::config_finite_field) with `gf.primitive_polynomial()`).
    fn config_finite_field_from(&mut self, gf: &crate::algebra::finite_field::GF256) {
        self.config_finite_field(gf.primitive_polynomial());
    }
}
