// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under the MIT License. See LICENSE for details.

use crate::types::Operation;

/// Trait for storing vectors and applying fountain [`Operation`]s.
///
/// Production (execute-only) paths should call the **direct** methods
/// (`mul_add`, `add_one_to_vector`, …). Those have default implementations that
/// build an [`Operation`] and call [`Self::execute`], so existing operators keep
/// working. Fast backends (slab/simd) override the direct methods to skip the
/// enum build + match.
///
/// Delayed / logged sessions keep using [`Self::execute`] with a recorded
/// [`Operation`] list.
pub trait DataOperator {
    /// Stores a byte vector under the given data ID.
    fn insert_vector(&mut self, _vector: &[u8], _data_id: usize) {}

    /// Retrieves a reference to the byte vector stored at the given data ID.
    fn get_vector(&self, _data_id: usize) -> &[u8] {
        &[]
    }

    /// Executes a single [`Operation`] on the stored vectors.
    fn execute(&mut self, operation: &Operation);

    /// Align stored-field arithmetic with [`crate::DataManager::config_finite_field`].
    /// Default: no-op. In-memory operators should rebuild GF(256) tables when `pp` is a real
    /// primitive polynomial; ignore [`crate::types::GF2_FIELD_POLY`].
    fn config_finite_field(&mut self, pp: u16) {
        let _ = pp;
    }

    /// Use the same GF(256) tables as `gf` (default: [`config_finite_field`](Self::config_finite_field) with `gf.primitive_polynomial()`).
    fn config_finite_field_from(&mut self, gf: &crate::algebra::finite_field::GF256) {
        self.config_finite_field(gf.primitive_polynomial());
    }

    // --- Direct ops: defaults route through `execute` ---

    #[deprecated(note = "use ensure_zero_range or ensure_zero_one")]
    fn ensure_zero(&mut self, list_id: &[usize]) {
        self.execute(&Operation::EnsureZero {
            list_id: list_id.to_vec(),
        });
    }

    fn ensure_zero_one(&mut self, id: usize) {
        self.execute(&Operation::EnsureZeroOne { id });
    }

    /// Zero consecutive vectors in the half-open range `[range.start, range.end)`.
    fn ensure_zero_range(&mut self, range: std::ops::Range<usize>) {
        self.execute(&Operation::EnsureZeroRange {
            start: range.start,
            end: range.end,
        });
    }

    fn multiply_alpha(&mut self, id: usize) {
        self.execute(&Operation::MultiplyAlpha { id });
    }

    fn multiply_scalar(&mut self, scalar: u8, id: usize) {
        self.execute(&Operation::MultiplyScalar { scalar, id });
    }

    fn add_one_to_vector(&mut self, src_id: usize, target_id: usize) {
        self.execute(&Operation::AddOneToVector { src_id, target_id });
    }

    fn add_two_to_vector(&mut self, s0: usize, s1: usize, target_id: usize) {
        self.execute(&Operation::AddTwoToVector { s0, s1, target_id });
    }

    fn add_three_to_vector(&mut self, s0: usize, s1: usize, s2: usize, target_id: usize) {
        self.execute(&Operation::AddThreeToVector {
            s0,
            s1,
            s2,
            target_id,
        });
    }

    fn add_to_vector(&mut self, list_id: &[usize], target_id: usize) {
        self.execute(&Operation::AddToVector {
            list_id: list_id.to_vec(),
            target_id,
        });
    }

    fn broadcast_add(&mut self, src_id: usize, target_ids: &[usize]) {
        self.execute(&Operation::BroadcastAdd {
            src_id,
            target_ids: target_ids.to_vec(),
        });
    }

    fn mul_add(&mut self, src_id: usize, scalar: u8, target_id: usize) {
        self.execute(&Operation::MulAdd {
            src_id,
            scalar,
            target_id,
        });
    }

    fn move_to(&mut self, src_id: usize, target_id: usize) {
        self.execute(&Operation::MoveTo { src_id, target_id });
    }

    /// Move `count` consecutive vectors starting at `src_start` to `dst_start`.
    fn move_to_n(&mut self, src_start: usize, dst_start: usize, count: usize) {
        self.execute(&Operation::MoveToN {
            src_start,
            dst_start,
            count,
        });
    }

    fn copy_to(&mut self, src_id: usize, target_id: usize) {
        self.execute(&Operation::CopyTo { src_id, target_id });
    }

    /// Copy `count` consecutive vectors starting at `src_start` to `dst_start`.
    fn copy_to_n(&mut self, src_start: usize, dst_start: usize, count: usize) {
        self.execute(&Operation::CopyToN {
            src_start,
            dst_start,
            count,
        });
    }

    fn remove(&mut self, id: usize) {
        self.execute(&Operation::Remove { id });
    }

    /// Informational coded↔data association (no byte work for most backends).
    fn info_coded_vector(&mut self, coded_id: usize, data_id: usize) {
        self.execute(&Operation::InfoCodedVector { coded_id, data_id });
    }
}
