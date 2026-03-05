// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

//! Trait abstractions for fountain code components.
//!
//! Defines the pluggable interfaces for code schemes, data operators,
//! and LDPC/HDPC precodes.

/// Code scheme configuration trait.
pub mod code_scheme;
/// Data operator trait for vector storage and execution.
pub mod data_operator;
/// High-Density Parity-Check precode trait.
pub mod hdpc;
/// Low-Density Parity-Check precode trait.
pub mod ldpc;

pub use code_scheme::*;
pub use data_operator::*;
pub use hdpc::*;
pub use ldpc::*;
