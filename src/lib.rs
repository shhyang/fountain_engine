// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

//! Fountain Code Core Library
//!
//! This library provides the core algorithms, traits, and data structures
//! for fountain code encoding and decoding, built on the
//! [`core::system_solver::SystemSolver`] (BP + inactivation + Gaussian elimination) stack.

/// Finite field arithmetic, linear algebra, and binary vector utilities.
pub mod algebra;
mod core;
mod data_manager;
/// Fountain code decoder (BP + inactivation + Gaussian elimination).
pub mod decoder;
/// Fountain code encoder with optional precoding.
pub mod encoder;
/// Trait definitions for code schemes, data operators, LDPC, and HDPC precodes.
pub mod traits;
/// Core type definitions: code parameters, operation variants, and status enums.
pub mod types;

pub use algebra::*;
pub use core::binary_matrix::BinaryMatrix;
pub use data_manager::DataManager;
pub use decoder::*;
pub use encoder::*;
pub use traits::*;
pub use types::*;
