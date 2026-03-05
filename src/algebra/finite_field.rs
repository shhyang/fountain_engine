// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

//! Finite field operations for erasure coding
//!
//! This module provides basic finite field arithmetic operations
//! over GF(2^m) fields commonly used in erasure coding.

/// GF256 implementation using lookup tables for efficient multiplication and division
pub struct GF256 {
    /// Multiplication table for GF(256)
    mul_table: [[u8; 256]; 256],
    /// Inverse table for GF(256)
    inv_table: [u8; 256],
    /// Primitive element alpha
    pub alpha: u8,
    /// Primitive polynomial without the first
    pub primitive_polynomial: u8,
}

impl GF256 {
    /// Helper function to multiply two field elements
    fn gf_multiply(a: u8, b: u8, primitive_polynomial: u8) -> u8 {
        let mut result = 0u8;
        let mut a = a;
        let mut b = b;

        while b != 0 {
            if b & 1 != 0 {
                result ^= a;
            }
            let carry = a & 0x80;
            a <<= 1;
            if carry != 0 {
                a ^= primitive_polynomial;
            }
            b >>= 1;
        }
        result
    }

    /// Create GF256 with the standard primitive polynomial x^8 + x^4 + x^3 + x + 1
    /// with the hex value of 0x11B
    pub fn new(primitive_polynomial: u16) -> Self {
        let mut mul_table = [[0u8; 256]; 256];
        let mut inv_table = [0u8; 256];

        // Generate multiplication table
        for (i, row) in mul_table.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = Self::gf_multiply(i as u8, j as u8, primitive_polynomial as u8);
            }
        }
        // Generate inverse table
        for i in 1..256 {
            for j in 1..256 {
                if mul_table[i][j] == 1 {
                    inv_table[i] = j as u8;
                    break;
                }
            }
        }
        Self {
            mul_table,
            inv_table,
            alpha: 0x02_u8,
            primitive_polynomial: primitive_polynomial as u8,
        }
    }

    /// Add two field elements (XOR operation)
    #[inline]
    pub fn add(&self, a: u8, b: u8) -> u8 {
        a ^ b
    }

    /// Multiply two field elements using lookup table
    #[inline]
    pub fn mul_lookup(&self, a: u8, b: u8) -> u8 {
        self.mul_table[a as usize][b as usize]
    }

    /// Multiplies two field elements with fast-paths for 0 and 1.
    pub fn multiply(&self, a: u8, b: u8) -> u8 {
        match a {
            0 => 0,
            1 => b,
            _ => self.mul_table[a as usize][b as usize],
        }
    }

    /// Multiplies a field element by the primitive element alpha, using bit-shift optimization.
    #[inline]
    pub fn mul_alpha(&self, a: u8) -> u8 {
        if a == 0 {
            return 0;
        }
        if a == 1 {
            return self.alpha;
        }
        if a < 128 {
            return a << 1;
        }
        self.add(a << 1, self.primitive_polynomial)
    }

    /// Get the multiplicative inverse of a field element
    #[inline]
    pub fn inverse(&self, a: u8) -> u8 {
        if a == 0 {
            panic!("Inverse of zero does not exist in finite field");
        }
        self.inv_table[a as usize]
    }

    /// Divide two field elements: a / b
    #[inline]
    pub fn divide(&self, a: u8, b: u8) -> u8 {
        if b == 0 {
            panic!("Division by zero in finite field");
        }
        self.multiply(a, self.inverse(b))
    }

    /// Vector addition in-place: `result[i] = result[i] + vec[i]`
    pub fn vector_addition_inplace(&self, result: &mut [u8], vec: &[u8]) {
        for (i, &vec_elem) in vec.iter().enumerate() {
            result[i] = self.add(result[i], vec_elem);
        }
    }

    /// Multiply a vector by alpha in-place: `result[i] = result[i] * alpha`
    pub fn multiply_alpha_inplace(&self, result: &mut [u8]) {
        for elem in result.iter_mut() {
            *elem = self.mul_alpha(*elem);
        }
    }

    /// Multiply a vector by a scalar in-place: `result[i] = scalar * result[i]`
    pub fn scalar_vector_multiply_inplace(&self, scalar: u8, result: &mut [u8]) {
        for elem in result.iter_mut() {
            *elem = self.multiply(scalar, *elem);
        }
    }
}

impl Default for GF256 {
    fn default() -> Self {
        Self::new(0x11B)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gf256_basic_operations() {
        let field = GF256::default();

        // Test addition (XOR)
        assert_eq!(field.add(5, 3), 6);
        assert_eq!(field.add(0, 7), 7);
        assert_eq!(field.add(7, 7), 0);

        // Test multiplication
        assert_eq!(field.multiply(0, 5), 0);
        assert_eq!(field.multiply(5, 0), 0);
        assert_eq!(field.multiply(1, 5), 5);
        assert_eq!(field.multiply(5, 1), 5);

        // Test inverse
        assert_eq!(field.multiply(5, field.inverse(5)), 1);
        assert_eq!(field.multiply(7, field.inverse(7)), 1);
    }

    #[test]
    #[should_panic(expected = "Inverse of zero")]
    fn test_inverse_of_zero() {
        let field = GF256::default();
        field.inverse(0);
    }

    #[test]
    fn test_gf256_multiplication_table() {
        let field = GF256::default();

        // Test that multiplication table is correct
        assert_eq!(field.multiply(2, 0), 0);
        assert_eq!(field.multiply(2, 3), 6);
        assert_eq!(field.multiply(7, 11), 49); // Corrected expected value
        assert_eq!(field.multiply(255, 1), 255);
        assert_eq!(field.multiply(255, 255), 19); // Corrected expected value
    }

    #[test]
    fn test_gf256_inverse_table() {
        let field = GF256::default();

        // Test that inverse table is correct
        for i in 1u8..=255u8 {
            let inv = field.inverse(i);
            assert_eq!(field.multiply(i, inv), 1);
        }
    }

    #[test]
    fn test_gf256_distributive_law() {
        let field = GF256::default();

        // Test distributive law: a * (b + c) = a * b + a * c
        let a = 5;
        let b = 3;
        let c = 7;

        let left = field.multiply(a, field.add(b, c));
        let right = field.add(field.multiply(a, b), field.multiply(a, c));

        assert_eq!(left, right);
    }

    #[test]
    fn test_gf256_associative_law() {
        let field = GF256::default();

        // Test associative law: (a * b) * c = a * (b * c)
        let a = 5;
        let b = 3;
        let c = 7;

        let left = field.multiply(field.multiply(a, b), c);
        let right = field.multiply(a, field.multiply(b, c));

        assert_eq!(left, right);
    }

    #[test]
    fn test_gf256_commutative_law() {
        let field = GF256::default();

        // Test commutative law: a * b = b * a
        let a = 5;
        let b = 3;

        assert_eq!(field.multiply(a, b), field.multiply(b, a));
    }

    #[test]
    fn test_gf256_zero_and_identity() {
        let field = GF256::default();

        // Test zero element properties
        for i in 0u8..=255u8 {
            assert_eq!(field.add(0, i), i);
            assert_eq!(field.add(i, 0), i);
            assert_eq!(field.multiply(0, i), 0);
            assert_eq!(field.multiply(i, 0), 0);
        }

        // Test identity element properties
        for i in 0u8..=255u8 {
            assert_eq!(field.multiply(1, i), i);
            assert_eq!(field.multiply(i, 1), i);
        }
    }

    #[test]
    fn test_gf256_division_consistency() {
        let field = GF256::default();

        // Test that division is consistent with multiplication
        for i in 1u8..=255u8 {
            for j in 1u8..=255u8 {
                let product = field.multiply(i, j);
                let quotient = field.multiply(product, field.inverse(i));
                assert_eq!(quotient, j);
            }
        }
    }

    #[test]
    fn test_mul_alpha_consistency() {
        let field = GF256::default();

        // Test that mul_alpha(x) equals multiply(x, alpha) for all x
        for x in 0u8..=255u8 {
            let result_mul_alpha = field.mul_alpha(x);
            let result_multiply = field.multiply(x, field.alpha);
            assert_eq!(
                result_mul_alpha, result_multiply,
                "mul_alpha({}) = {} but multiply({}, alpha) = {}",
                x, result_mul_alpha, x, result_multiply
            );
        }
    }
}
