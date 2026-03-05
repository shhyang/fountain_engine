// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

/// Utilities for element-wise binary (GF(2)) vector operations.
pub struct BinaryVector;

impl BinaryVector {
    /// XORs vector `b` into `a` element-wise in place. Panics if lengths differ.
    pub fn add_inplace(a: &mut [u8], b: &[u8]) {
        assert_eq!(a.len(), b.len());
        for i in 0..a.len() {
            a[i] ^= b[i];
        }
    }
}
