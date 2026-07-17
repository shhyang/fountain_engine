// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

/// Maximum number of inactive vectors allowed in the system.
pub const MAX_INACTIVE_NUM: usize = 100; // Adjust this value as needed

/// Sentinel for [`HDPC::gf_poly`](crate::traits::HDPC::gf_poly) and
/// [`DataManager::config_finite_field`](crate::DataManager::config_finite_field).
///
/// Not a real degree-8 polynomial. It selects a **GF(2)-only** session:
/// - `DataManager::gf256()` is `None`; LU decomposition uses [`GF2`](crate::algebra::finite_field::GF2).
/// - HDPC matrix multiply uses XOR only (`mul_data` / `mul_binary`); no GF(256) scalars.
/// - `divide_scalar` / `multiply_scalar` on the manager only allow trivial scalars (e.g. pivot `1`).
///
/// Used by binary precodes (e.g. Raptor-10 `R10HDPC`). GF(256) HDPC (Raptor-Q style) returns a real
/// polynomial such as `0x11D`.
pub const GF2_FIELD_POLY: u16 = 0x001;

/// Type alias for degree set generator function (symbol index -> active and inactive indices).
pub type DegreeSetFn = Box<dyn FnMut(usize) -> Vec<usize>>;

/// A struct for decoding configuration.
#[derive(Clone, Debug)]
pub struct DecodingConfig {
    pub max_inactive_num: usize,
    //move this value to CodeParams
    //pub num_padding: usize,
    pub inac_strategy: InactivationStrategy,
    pub subs_method: SubstitutionMethod,
}

impl Default for DecodingConfig {
    fn default() -> Self {
        Self {
            max_inactive_num: 0,
            //move this value to CodeParams
            //num_padding: 0,
            inac_strategy: InactivationStrategy::ByIndex,
            subs_method: SubstitutionMethod::Direct,
        }
    }
}

impl DecodingConfig {
    pub fn with_max_inact_num(mut self, max_inactive_num: usize) -> Self {
        self.max_inactive_num = max_inactive_num;
        self
    }
}

/// Parameters defining the structure of a fountain code.
///
/// Specifies the number of message, active, inactive, LDPC, and HDPC vectors
/// that together determine the code's parity-check matrix layout.
#[derive(Clone, Debug)]
pub struct CodeParams {
    /// Total number of message vectors.
    pub k: usize,
    /// Number of active message vectors (decoded via BP).
    pub a: usize,
    /// Number of inactive message vectors (decoded via GE), equal to `k - a`.
    pub b: usize,
    /// Number of LDPC (Low-Density Parity-Check) constraint vectors.
    pub l: usize,
    /// Number of HDPC (High-Density Parity-Check) constraint vectors.
    pub h: usize,
    // Number of pre-inactive variable vectors.
    pub i: usize,
    // Number of padding vectors.
    //pub p: usize,
}

impl CodeParams {
    /// Creates new code parameters with pre-inactivation.
    /// - `b` is computed as `k - a`.
    /// - `i` is computed as `b + h`.
    pub fn new(k: usize, a: usize, l: usize, h: usize) -> Self {
        let b = k - a;
        let i = b + h;
        Self { k, a, b, l, h, i }
    }

    /// Creates new code parameters without pre-inactivation.
    pub fn new_without_pre_inact(k: usize, l: usize, h: usize) -> Self {
        Self {
            k,
            a: k,
            b: 0,
            l,
            h,
            i: 0,
        }
    }
    /*
    /// Create new code parameters with padding vectors.
    pub fn with_padding(mut self, p: usize) -> Self {
        if p < self.a {
            self.p = p;
        } else {
            eprintln!(
                "with_padding: padding count {p} must be less than active message vectors (a={}); ignoring",
                self.a
            );
        }
        self
    }*/

    /// Returns the number of active variable vectors (`a + l`).
    pub fn num_active(&self) -> usize {
        self.a + self.l
    }

    /// Returns the number of pre-inactive variable vectors (`i`).
    pub fn num_pre_inactive(&self) -> usize {
        self.i
    }

    /// Returns the number of inactive variable vectors (`b + h`).
    pub fn num_inactive(&self) -> usize {
        self.b + self.h
    }

    /// Returns the combined count of message and LDPC vectors (`k + l`).
    pub fn num_message_ldpc(&self) -> usize {
        self.k + self.l
    }

    /// Returns the total number of variable vectors (`k + l + h`).
    pub fn num_total(&self) -> usize {
        self.k + self.l + self.h
    }

    pub fn has_precode(&self) -> bool {
        self.l > 0 || self.h > 0
    }

    pub fn num_message(&self) -> usize {
        self.k
    }
}

/// Distinguishes between systematic and ordinary fountain codes.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum CodeType {
    /// Ordinary (non-systematic) encoding and decoding.
    Ordinary,
    /// Systematic encoding and decoding, where source symbols appear unmodified in the output.
    Systematic,
}

/// Identifies which solver configuration to use based on code type and direction.
#[derive(Copy, Clone, PartialEq, Debug)]
pub enum SolverType {
    /// Ordinary encoding solver.
    OrdEnc,
    /// Ordinary decoding solver.
    OrdDec,
    /// Systematic encoding solver.
    SysEnc,
    /// Systematic decoding solver.
    SysDec,
}

/// Strategy for selecting which variable to inactivate during BP decoding.
#[derive(Copy, Clone, PartialEq, Debug)]
pub enum InactivationStrategy {
    /// Inactivate by original variable index order.
    ByIndex,
    /// Inactivate the first remaining active variable.
    FirstActive,
    /// Inactivate the last remaining active variable.
    LastActive,
    /// Inactivate the variable with minimum remaining degree.
    MinDegree,
    /// Use a trial-run heuristic to choose the best inactivation.
    TrailRun,
}

/// Method used for back-substitution after Gaussian elimination.
#[derive(Copy, Clone, PartialEq, Debug)]
pub enum SubstitutionMethod {
    /// Solve directly using the GE-reduced matrix.
    Direct,
    /// Solve using the original constraint matrix.
    Original,
}

/// Result of the decoding process.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum DecodeStatus {
    /// All source vectors have been successfully recovered.
    Decoded,
    /// Decoding is not yet complete.
    NotDecoded,
}

/// Data vector operations that can be recorded for delayed execution.
///
/// Each variant represents a single atomic operation on data vectors
/// identified by their unique IDs. The [`DataManager`](crate::DataManager)
/// records these operations and optionally forwards them to a
/// [`DataOperator`](crate::traits::DataOperator) for immediate execution.
#[derive(Debug, Clone)]
pub enum Operation {
    /// Zero out the vectors with the given IDs.
    EnsureZero {
        /// IDs of vectors to zero.
        list_id: Vec<usize>,
    },
    EnsureZeroOne {
        /// ID of the vector to zero.
        id: usize,
    },
    /// Multiply a vector element-wise by the primitive element alpha in GF(256).
    MultiplyAlpha {
        /// ID of the vector to multiply.
        id: usize,
    },
    /// Multiply a vector element-wise by a GF(256) scalar.
    MultiplyScalar {
        /// GF(256) scalar multiplier.
        scalar: u8,
        /// ID of the vector to multiply.
        id: usize,
    },
    /// XOR one source vector into a target (GF(2) hot path; no `Vec` allocation).
    AddOneToVector {
        /// ID of the source vector.
        src_id: usize,
        /// ID of the destination vector.
        target_id: usize,
    },
    /// XOR two source vectors into a target (GF(2) hot path; no `Vec` allocation).
    AddTwoToVector {
        /// IDs of source vectors.
        s0: usize,
        s1: usize,
        /// ID of the destination vector.
        target_id: usize,
    },
    /// XOR three source vectors into a target (GF(2) hot path; no `Vec` allocation).
    AddThreeToVector {
        /// IDs of source vectors.
        s0: usize,
        s1: usize,
        s2: usize,
        /// ID of the destination vector.
        target_id: usize,
    },
    /// XOR (add in GF(2)/GF(256)) multiple source vectors into a target vector.
    AddToVector {
        /// IDs of source vectors to XOR into the target.
        list_id: Vec<usize>,
        /// ID of the destination vector.
        target_id: usize,
    },
    /// XOR a single source vector into each of several target vectors.
    BroadcastAdd {
        /// ID of the source vector.
        src_id: usize,
        /// IDs of destination vectors.
        target_ids: Vec<usize>,
    },
    /// Multiply source vector by a scalar and XOR the result into the target: `target += scalar * src`.
    MulAdd {
        /// ID of the source vector.
        src_id: usize,
        /// GF(256) scalar multiplier.
        scalar: u8,
        /// ID of the destination vector.
        target_id: usize,
    },
    /// Move a vector from `src_id` to `target_id`, invalidating the source.
    MoveTo {
        /// ID of the source vector.
        src_id: usize,
        /// ID of the destination vector.
        target_id: usize,
    },
    /// Copy a vector from `src_id` to `target_id`, keeping the source intact.
    CopyTo {
        /// ID of the source vector.
        src_id: usize,
        /// ID of the destination vector.
        target_id: usize,
    },
    /// Remove (deallocate) a vector.
    Remove {
        /// ID of the vector to remove.
        id: usize,
    },
    /// Informational marker associating a coded vector ID with its data vector ID.
    InfoCodedVector {
        /// External coded vector identifier.
        coded_id: usize,
        /// Internal data vector identifier.
        data_id: usize,
    },
}
