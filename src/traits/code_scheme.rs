// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use super::{HDPC, LDPC};
use crate::types::{
    SubstitutionMethod, CodeParams, CodeType, DegreeSetFn, InactivationStrategy, DecodingConfig, PreInactivation,
};

pub type PrecodePair = (Option<Box<dyn HDPC>>, Option<Box<dyn LDPC>>);

/// Trait for code configuration
pub trait CodeScheme {
    /// Get code parameters
    fn get_params(&self) -> CodeParams;

    /// Get code type (systematic or ordinary)
    fn code_type(&self) -> CodeType;

    /// Create degree set generator
    fn create_degree_set_fn(&self) -> DegreeSetFn;

    /// Create precode components
    fn create_precode(&self) -> PrecodePair;

    /// Maximum number of inactivations
    fn max_inactive_num(&self) -> usize{
        self.get_params().h + self.get_params().b 
    }

    /// Decoding configuration
    fn decoding_config(&self) -> DecodingConfig {
        DecodingConfig {
            max_inactive_num: self.max_inactive_num(),
            num_padding: 0,
            pre_inactivation: PreInactivation::YesPre,
            inac_strategy: InactivationStrategy::ByIndex,
            subs_method: SubstitutionMethod::Direct,
        }
    }
}
