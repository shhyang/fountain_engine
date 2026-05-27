// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

#![allow(clippy::needless_range_loop)]

use crate::algebra::finite_field::Field;
use crate::algebra::finite_field::GF256;

pub struct Vector;

impl Vector {
    pub fn add_inplace(a: &mut [u8], b: &[u8]) {
        assert_eq!(a.len(), b.len());
        for i in 0..a.len() {
            a[i] ^= b[i];
        }
    }

    /// Multiply a vector by alpha in-place: `result[i] = result[i] * alpha`
    pub fn multiply_alpha_inplace(field: &GF256, result: &mut [u8]) {
        for elem in result.iter_mut() {
            *elem = field.mul_alpha(*elem);
        }
    }

    /// Multiply a vector by a scalar in-place: `result[i] = scalar * result[i]`
    pub fn scalar_vector_multiply_inplace<F: Field>(field: &F, scalar: u8, result: &mut [u8]) {
        for elem in result.iter_mut() {
            *elem = field.mul(scalar, *elem);
        }
    }
}

/// Matrix operations over GF(256), including multiplication, permutation, and LU decomposition.
pub struct Matrix;

impl Matrix {
    /// Multiplies two matrices `a` (m×n) and `b` (n×p) over GF(256), returning the m×p product.
    pub fn multiply<F: Field>(field: &F, a: &[Vec<u8>], b: &[Vec<u8>]) -> Vec<Vec<u8>> {
        let m = a.len();
        let n = if m > 0 { a[0].len() } else { 0 };
        if n != b.len() {
            panic!("The number of columns in A must be equal to the number of rows in B");
        }
        let p = if n > 0 { b[0].len() } else { 0 };
        let mut result = vec![vec![0u8; p]; m];
        for (i, row) in result.iter_mut().enumerate() {
            for (j, &a_ij) in a[i].iter().take(n).enumerate() {
                for (k, &val) in b[j].iter().enumerate() {
                    row[k] ^= field.mul(a_ij, val);
                }
            }
        }
        result
    }

    /// Permuate the rows of a matrix inplace using swap
    pub fn permute_rows_inplace(a: &mut [Vec<u8>], p: &[usize]) {
        let n = a.len();
        let mut visited = vec![false; n];
        for i in 0..n {
            if visited[i] || p[i] == i {
                continue;
            }
            let mut j = i;
            while !visited[j] {
                visited[j] = true;
                let k = p[j];
                if k != i {
                    a.swap(j, k);
                }
                j = k;
            }
        }
    }

    /// Perform LU decomposition of A. Return the permutation vector p and the rank r.
    /// This modifies A in-place to store the LU decomposition.
    pub fn lu_decomp<F: Field>(field: &F, a: &mut [Vec<u8>]) -> (Vec<usize>, usize) {
        let m = a.len();
        let n = if m > 0 { a[0].len() } else { 0 };
        let mut p: Vec<usize> = (0..m).collect();
        let mut i = 0;

        for j in 0..n {
            let mut pivot_found = false;
            for k in i..a.len() {
                if a[k][j] != 0 {
                    p.swap(i, k);
                    a.swap(i, k);
                    pivot_found = true;
                    break;
                }
            }
            if pivot_found {
                // (i,j) entry is non-zero
                for k in i + 1..m {
                    let l = field.divide(a[k][j], a[i][j]);
                    a[k][j] = 0;
                    a[k][i] = l;

                    // Update the rest of the row
                    for col in (j + 1)..n {
                        a[k][col] = field.add(a[k][col], field.mul(l, a[i][col]));
                    }
                }
                i += 1;
                if i == m {
                    break;
                }
            }
        }

        (p, i)
    }

    /// Perform LU decomposition of A incrementally. Return the permutation vector p and the rank r.
    /// This modifies A in-place to store the LU decomposition, and update q in-place, so that
    /// UQ is the upper triangular matrix of the LU decomposition.
    pub fn lu_decomp_incr<F: Field>(
        field: &F,
        a: &mut [Vec<u8>],
        q: &mut [usize],
        r: usize,
    ) -> (Vec<usize>, usize) {
        let m = a.len();
        let n = if m > 0 { a[0].len() } else { 0 };
        let mut p = (0..m).collect::<Vec<_>>();

        for i in 0..r {
            for k in r..m {
                let l = field.divide(a[k][q[i]], a[i][q[i]]);
                a[k][q[i]] = l;
                for col in i + 1..n {
                    a[k][q[col]] = field.add(a[k][q[col]], field.mul(l, a[i][q[col]]));
                }
            }
        }

        let mut i = r;

        for j in r..n {
            let mut pivot_found = false;
            for k in i..a.len() {
                if a[k][q[j]] != 0 {
                    p.swap(i, k);
                    a.swap(i, k);
                    pivot_found = true;
                    break;
                }
            }

            if pivot_found {
                // at (i,j)
                q.swap(i, j);
                for k in i + 1..m {
                    let l = field.divide(a[k][q[i]], a[i][q[i]]);
                    a[k][q[i]] = l;

                    // Update the rest of the row
                    for col in i + 1..n {
                        a[k][q[col]] = field.add(a[k][q[col]], field.mul(l, a[i][q[col]]));
                    }
                }
                i += 1;
                if i == m {
                    break;
                }
            }
        }
        (p, i)
    }
}

/// Linear solver for systems over GF(256)
/// This is not used in the code, but it is a useful function for testing.
pub struct LinearSys;

impl LinearSys {
    /// Solve the linear system Ax = b.
    /// A is a matrix, and b is a vector.
    /// Returns Some(x) if the system has a solution, None otherwise.
    pub fn lin_solve<F: Field>(
        field: &F,
        a: &mut [Vec<u8>],
        b: &mut [Vec<u8>],
    ) -> Result<(), String> {
        let m = a.len();
        let n = if m > 0 { a[0].len() } else { 0 };

        let (p, r) = Matrix::lu_decomp(field, a);

        if r == n {
            // Take only the first n rows for solving
            let a_n = &mut a[..n];
            Matrix::permute_rows_inplace(b, &p);
            let b_n = &mut b[..n];
            Self::lu_solve(field, a_n, b_n)?;
            Ok(())
        } else {
            Err(format!("The matrix is not invertible, rank = {}", r))
        }
    }

    /// Solve the linear system Ax = b given the LU decomposition of A.
    /// This performs forward and backward substitution.
    pub fn lu_solve<F: Field>(field: &F, a: &[Vec<u8>], b: &mut [Vec<u8>]) -> Result<(), String> {
        let n = a.len();
        if n != b.len() {
            return Err("The number of rows in A and b must be the same".to_string());
        }

        /// Perform combined vector operation in-place: result[i] = result[i] + scalar * vec[i]
        #[inline]
        fn combined_vector_operation_inplace<F: Field>(
            field: &F,
            b: &mut [Vec<u8>],
            i: usize,
            scalar: u8,
            j: usize,
        ) {
            for k in 0..b[i].len() {
                b[i][k] = field.add(b[i][k], field.mul(scalar, b[j][k]));
            }
        }

        // Forward substitution (L part)
        for j in 0..n - 1 {
            for i in j + 1..n {
                combined_vector_operation_inplace(field, b, i, a[i][j], j);
            }
        }

        // Backward substitution (U part)
        for j in (0..n).rev() {
            Vector::scalar_vector_multiply_inplace(field, field.inverse(a[j][j]), &mut b[j]);
            for i in (0..j).rev() {
                combined_vector_operation_inplace(field, b, i, a[i][j], j);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_linear_system() {
        let field = GF256::default();

        let mut a = vec![vec![1, 1], vec![2, 1]];

        let x = vec![vec![3, 7, 19], vec![5, 6, 20]];

        let mut b = Matrix::multiply(&field, &a, &x);
        let result = LinearSys::lin_solve(&field, &mut a, &mut b);
        assert!(result.is_ok());

        assert_eq!(b, x);
    }

    #[test]
    fn test_singular_matrix() {
        let field = GF256::default();

        // Test a singular matrix (rank < n)
        let mut a = vec![
            vec![1, 1],
            vec![1, 1], // Same as first row
        ];
        let mut b = vec![vec![3], vec![5]];

        let result = LinearSys::lin_solve(&field, &mut a, &mut b);
        assert!(result.is_err());
    }
}
