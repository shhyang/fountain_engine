// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use super::{HDPC, LDPC};
use crate::types::{
    BackSubstitutionMethod, CodeParams, CodeType, DegreeSetFn, InactivationStrategy,
};

/// Trait for code configuration
pub trait CodeScheme {
    /// Get code parameters
    fn get_params(&self) -> CodeParams;

    /// Get code type (systematic or ordinary)
    fn code_type(&self) -> CodeType;

    /// Create degree set generator
    fn create_degree_set_fn(&self) -> DegreeSetFn;

    /// Create precode components
    #[allow(clippy::type_complexity)]
    fn create_precode(&self) -> (Option<Box<dyn HDPC>>, Option<Box<dyn LDPC>>);

    /// Maximum number of inactivations
    fn max_inactive_num(&self) -> usize{
        self.get_params().h + self.get_params().b 
    }

    /// Inactivation strategy
    fn inac_strategy(&self) -> InactivationStrategy {
        InactivationStrategy::ByIndex
    }

    /// Back substitution method
    fn back_substitution_method(&self) -> BackSubstitutionMethod {
        BackSubstitutionMethod::Direct
    }
}
