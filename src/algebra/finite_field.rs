// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

/// Trait for finite field operations over u8 elements.
pub trait Field {
    fn add(&self, a: u8, b: u8) -> u8 {
        a ^ b
    }
    /// Default is bitwise AND — only meaningful when elements are in `{0, 1}`.
    /// [`GF2`] and [`GF256`] override this.
    fn mul(&self, a: u8, b: u8) -> u8 {
        a & b
    }
    fn inverse(&self, a: u8) -> u8 {
        a
    }
    fn divide(&self, a: u8, _b: u8) -> u8 {
        a
    }
}

/// Monic primitive degree-8 polynomials over GF(2) for generator `0x02` (not AES `0x11B`).
pub const PRIMITIVE_POLYNOMIALS: [u16; 16] = [
    0x11D, 0x12B, 0x12D, 0x14D, 0x15F, 0x163, 0x165, 0x169, 0x171, 0x187, 0x18D, 0x1A9, 0x1C3,
    0x1CF, 0x1E7, 0x1F5,
];

/// Returns true if `pp` is one of [`PRIMITIVE_POLYNOMIALS`].
#[must_use]
pub fn is_primitive_polynomial(pp: u16) -> bool {
    PRIMITIVE_POLYNOMIALS.contains(&pp)
}

/// Error from [`GF256::try_new_with_primitive_polynomial`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InvalidPrimitivePolynomial(pub u16);

impl std::fmt::Display for InvalidPrimitivePolynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "invalid primitive polynomial: 0x{:03X}", self.0)
    }
}

impl std::error::Error for InvalidPrimitivePolynomial {}

/// GF(2) for LU when [`crate::types::GF2_FIELD_POLY`] is configured (binary HDPC / 0–1 matrices).
pub struct GF2 {}

impl Field for GF2 {
    fn mul(&self, a: u8, b: u8) -> u8 {
        if a == 0 || b == 0 { 0 } else { 1 }
    }

    fn inverse(&self, a: u8) -> u8 {
        if a == 0 {
            panic!("Inverse of zero does not exist in GF(2)");
        }
        1
    }

    fn divide(&self, a: u8, b: u8) -> u8 {
        if b == 0 {
            panic!("Division by zero in GF(2)");
        }
        a
    }
}

impl GF2 {
    pub fn new() -> Self {
        Self {}
    }
}

impl Default for GF2 {
    fn default() -> Self {
        Self::new()
    }
}

/// GF256 implementation using lookup tables for efficient multiplication and division.
#[derive(Clone)]
pub struct GF256 {
    /// Multiplication table for GF(256)
    mul_table: [[u8; 256]; 256],
    /// Inverse table for GF(256)
    inv_table: [u8; 256],
    /// Primitive element alpha
    /// Primitive polynomial without the first bit
    primitive_polynomial: u8,
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

    /// Builds GF(256) for a primitive polynomial in [`PRIMITIVE_POLYNOMIALS`].
    ///
    /// # Errors
    ///
    /// Returns [`InvalidPrimitivePolynomial`] if `pp` is not whitelisted (e.g. AES `0x11B`).
    pub fn try_new_with_primitive_polynomial(pp: u16) -> Result<Self, InvalidPrimitivePolynomial> {
        if !is_primitive_polynomial(pp) {
            return Err(InvalidPrimitivePolynomial(pp));
        }
        Ok(Self::from_primitive_polynomial_unchecked(pp))
    }

    /// Like [`try_new_with_primitive_polynomial`](Self::try_new_with_primitive_polynomial) but panics on invalid `pp`.
    pub fn new_with_primitive_polynomial(pp: u16) -> Self {
        Self::try_new_with_primitive_polynomial(pp).unwrap_or_else(|e| panic!("{e}"))
    }

    fn from_primitive_polynomial_unchecked(pp: u16) -> Self {
        let poly_byte = (pp & 0xFF) as u8;
        let mut mul_table = [[0u8; 256]; 256];
        let mut inv_table = [0u8; 256];

        for (i, row) in mul_table.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = Self::gf_multiply(i as u8, j as u8, poly_byte);
            }
        }
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
            primitive_polynomial: poly_byte,
        }
    }

    pub fn primitive_element(&self) -> u8 {
        0x02_u8
    }

    pub fn primitive_polynomial(&self) -> u16 {
        0x100 + self.primitive_polynomial as u16
    }

    /// Multiply two field elements using lookup table
    #[inline]
    pub fn mul_lookup(&self, a: u8, b: u8) -> u8 {
        self.mul_table[a as usize][b as usize]
    }

    /// Multiplies a field element by the primitive element alpha, using bit-shift optimization.
    #[inline]
    pub fn mul_alpha(&self, a: u8) -> u8 {
        if a == 0 {
            return 0;
        }
        if a == 1 {
            return self.primitive_element();
        }
        if a < 128 {
            return a << 1;
        }
        self.add(a << 1, self.primitive_polynomial)
    }
}

impl Field for GF256 {
    fn mul(&self, a: u8, b: u8) -> u8 {
        match a {
            0 => 0,
            1 => b,
            _ => self.mul_table[a as usize][b as usize],
        }
    }
    fn inverse(&self, a: u8) -> u8 {
        if a == 0 {
            panic!("Inverse of zero does not exist in finite field");
        }
        self.inv_table[a as usize]
    }
    fn divide(&self, a: u8, b: u8) -> u8 {
        if b == 0 {
            panic!("Division by zero in finite field");
        }
        self.mul(a, self.inverse(b))
    }
}

impl Default for GF256 {
    fn default() -> Self {
        Self::new_with_primitive_polynomial(0x11D)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primitive_polynomial_whitelist() {
        assert_eq!(PRIMITIVE_POLYNOMIALS.len(), 16);
        assert!(is_primitive_polynomial(0x11D));
        assert!(!is_primitive_polynomial(0x11B));
    }

    #[test]
    fn try_new_rejects_aes_irreducible() {
        assert!(matches!(
            GF256::try_new_with_primitive_polynomial(0x11B),
            Err(InvalidPrimitivePolynomial(0x11B))
        ));
    }

    #[test]
    fn test_gf256_basic_operations() {
        let field = GF256::default();

        // Test addition (XOR)
        assert_eq!(field.add(5, 3), 6);
        assert_eq!(field.add(0, 7), 7);
        assert_eq!(field.add(7, 7), 0);

        // Test multiplication
        assert_eq!(field.mul(0, 5), 0);
        assert_eq!(field.mul(5, 0), 0);
        assert_eq!(field.mul(1, 5), 5);
        assert_eq!(field.mul(5, 1), 5);

        // Test inverse
        assert_eq!(field.mul(5, field.inverse(5)), 1);
        assert_eq!(field.mul(7, field.inverse(7)), 1);
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
        assert_eq!(field.mul(2, 0), 0);
        assert_eq!(field.mul(2, 3), 6);
        assert_eq!(field.mul(7, 11), 49); // Corrected expected value
        assert_eq!(field.mul(255, 1), 255);
        assert_eq!(field.mul(255, 255), 226); // Corrected expected value
    }

    #[test]
    fn test_gf256_inverse_table() {
        let field = GF256::default();

        // Test that inverse table is correct
        for i in 1u8..=255u8 {
            let inv = field.inverse(i);
            assert_eq!(field.mul(i, inv), 1);
        }
    }

    #[test]
    fn test_gf256_distributive_law() {
        let field = GF256::default();

        // Test distributive law: a * (b + c) = a * b + a * c
        let a = 5;
        let b = 3;
        let c = 7;

        let left = field.mul(a, field.add(b, c));
        let right = field.add(field.mul(a, b), field.mul(a, c));

        assert_eq!(left, right);
    }

    #[test]
    fn test_gf256_associative_law() {
        let field = GF256::default();

        // Test associative law: (a * b) * c = a * (b * c)
        let a = 5;
        let b = 3;
        let c = 7;

        let left = field.mul(field.mul(a, b), c);
        let right = field.mul(a, field.mul(b, c));

        assert_eq!(left, right);
    }

    #[test]
    fn test_gf256_commutative_law() {
        let field = GF256::default();

        // Test commutative law: a * b = b * a
        let a = 5;
        let b = 3;

        assert_eq!(field.mul(a, b), field.mul(b, a));
    }

    #[test]
    fn test_gf256_zero_and_identity() {
        let field = GF256::default();

        // Test zero element properties
        for i in 0u8..=255u8 {
            assert_eq!(field.add(0, i), i);
            assert_eq!(field.add(i, 0), i);
            assert_eq!(field.mul(0, i), 0);
            assert_eq!(field.mul(i, 0), 0);
        }

        // Test identity element properties
        for i in 0u8..=255u8 {
            assert_eq!(field.mul(1, i), i);
            assert_eq!(field.mul(i, 1), i);
        }
    }

    #[test]
    fn test_gf256_division_consistency() {
        let field = GF256::default();

        // Test that division is consistent with multiplication
        for i in 1u8..=255u8 {
            for j in 1u8..=255u8 {
                let product = field.mul(i, j);
                let quotient = field.mul(product, field.inverse(i));
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
            let result_multiply = field.mul(x, field.primitive_element());
            assert_eq!(
                result_mul_alpha, result_multiply,
                "mul_alpha({}) = {} but multiply({}, alpha) = {}",
                x, result_mul_alpha, x, result_multiply
            );
        }
    }
}
