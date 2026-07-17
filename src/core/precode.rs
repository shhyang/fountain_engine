// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use crate::data_manager::DataManager;
use crate::traits::CodeScheme;
use crate::traits::HDPC;
use crate::traits::LDPC;
use crate::types::CodeParams;

/// *Precode Encoder*
/// The precode encoder is used to encode the message vectors with a precode similar to that of RQ precoding.
/// To use the encoder, a data manager is needed, which implements the `DataManager` trait.
pub fn ordinary_precode_encode<T: CodeScheme>(manager: &mut DataManager, custom: &T) {
    let (hdpc, ldpc) = custom.create_precode();
    let params = custom.get_params();

    if let Some(ldpc) = ldpc.as_ref() {
        ldpc.encode(manager, &params);

        // Calculate D * [B_a; L; B_b] where D is the HDPC matrix, B is message packets, L is LDPC packets
        if let Some(hdpc) = hdpc {
            // GF(256) HDPC: real `pp`; binary HDPC (e.g. R10): `GF2_FIELD_POLY` → LU over GF(2).
            manager.config_finite_field(hdpc.gf_poly());
            let variable_ids_kl: Vec<usize> = manager.data_id_range_of_msg_ldpc_variable();
            let mut hdpc_constraint_data_ids = Vec::with_capacity(params.h);
            for _ in 0..params.h {
                hdpc_constraint_data_ids.push(manager.temp_data_id());
            }
            hdpc.mul_data(
                manager,
                &params,
                &variable_ids_kl,
                &hdpc_constraint_data_ids,
            );

            hdpc_solve(
                manager,
                &params,
                hdpc.as_ref(),
                ldpc.as_ref(),
                &hdpc_constraint_data_ids,
            );
        }
    }
}

/// Solve `(I' + D_s S_h) Y = X` for HDPC variable vectors.
fn hdpc_solve(
    manager: &mut DataManager,
    params: &CodeParams,
    hdpc: &dyn HDPC,
    ldpc: &dyn LDPC,
    hdpc_ids: &[usize],
) {
    // Performance (not implemented): cache `(p, lu)` from [`HDPC::lu_idssh`] keyed by
    // `(CodeParams, HDPC/LDPC scheme identity, gf_poly)` — static per K; per-encode work
    // remains RHS permutation, `lu_solve`, and LDPC substitution below.
    let (p, mut lu) = hdpc.lu_idssh(manager.gf256(), params, ldpc);

    let variable_ids_hdpc: Vec<usize> = manager.data_id_range_of_hdpc_variable();
    for i in 0..params.h {
        manager.move_to(hdpc_ids[p[i]], variable_ids_hdpc[i]);
    }

    lu_solve(manager, &mut lu, &variable_ids_hdpc);

    // Substitute P into L = S_h P + L'
    for (var_col, &var_id) in variable_ids_hdpc.iter().take(params.h).enumerate() {
        let ids = ldpc.inactive_column(var_col + params.b);
        let ids = ids
            .iter()
            .map(|&id| manager.data_id_of_ldpc_variable(id))
            .collect::<Vec<_>>();
        manager.broadcast_add_owned(var_id, ids);
    }
}

/// Forward and backward substitution for an LU-decomposed matrix on data vectors.
fn lu_solve(manager: &mut DataManager, matrix_a: &mut [Vec<u8>], target_ids: &[usize]) {
    if matrix_a.len() != target_ids.len() {
        panic!("The number of rows in A must be equal to the number of target IDs");
    }

    let n = matrix_a.len();

    for j in 0..n - 1 {
        for i in j + 1..n {
            if matrix_a[i][j] != 0 {
                manager.mul_add(target_ids[j], matrix_a[i][j], target_ids[i]);
            }
        }
    }

    for j in (0..n).rev() {
        if matrix_a[j][j] == 0 {
            panic!(
                "Singular matrix: diagonal element at position {} is zero",
                j
            );
        }

        manager.divide_scalar(matrix_a[j][j], target_ids[j]);

        for i in (0..j).rev() {
            if matrix_a[i][j] != 0 {
                manager.mul_add(target_ids[j], matrix_a[i][j], target_ids[i]);
            }
        }
    }
}
