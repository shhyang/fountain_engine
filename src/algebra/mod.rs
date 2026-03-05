// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

//! Algebraic primitives for fountain code arithmetic.
//!
//! Provides GF(2^8) finite field operations, matrix/linear-system solvers,
//! and binary (GF(2)) vector utilities.

/// GF(2^8) finite field with lookup-table-based arithmetic.
pub mod finite_field;
/// Matrix operations and LU-based linear system solvers over GF(256).
pub mod linear_algebra;
/// Binary (GF(2)) vector XOR operations.
pub mod binary_vector;

pub use finite_field::*;
pub use linear_algebra::*;
pub use binary_vector::*;