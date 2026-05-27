// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use crate::algebra::linear_algebra::Matrix;
use crate::data_manager::DataManager;
use crate::traits::CodeScheme;
use crate::traits::HDPC;
use crate::traits::LDPC;
use crate::types::CodeParams;

use crate::algebra::finite_field::GF2;
/// *Precode Encoder*
/// The precode encoder is used to encode the message vectors with a precode similar to that of RQ precoding.
/// To use the encoder, a data manager is needed, which implements the `DataManager` trait.
pub fn precode_encode<T: CodeScheme>(manager: &mut DataManager, params: &CodeParams, custom: &T) {
    let (hdpc, ldpc) = custom.create_precode();
    //let ldpc = LDPCType::new(ldpc_type, params.clone());

    //ldpc_encode(manager, params.a, params.b, ldpc.as_ref());
    if let Some(ldpc) = ldpc.as_ref() {
        ldpc.encode(manager, params);
    }

    // Calculate D * [B_a; L; B_b] where D is the HDPC matrix, B is message packets, L is LDPC packets
    if let Some(hdpc) = hdpc {
        //dbg!("hdpc encode");
        // GF(256) HDPC: real `pp`; binary HDPC (e.g. R10): `GF2_FIELD_POLY` → LU over GF(2).
        manager.config_finite_field(hdpc.gf_poly());
        let variable_ids_kl: Vec<usize> = manager.data_id_range_of_msg_ldpc_variable();
        let mut hdpc_constraint_data_ids = Vec::with_capacity(params.h);
        for _ in 0..params.h {
            hdpc_constraint_data_ids.push(manager.temp_data_id());
        }
        hdpc.mul_data(manager, params, &variable_ids_kl, &hdpc_constraint_data_ids);

        //for i in params.hdpc_range() {
        //     dbg!("mul_vector", manager.get_vector(i));
        //}

        hdpc_solve(
            manager,
            params,
            &hdpc,
            ldpc.as_ref().unwrap(),
            &hdpc_constraint_data_ids,
        );
    }
}

/*fn ldpc_encode(manager: &mut DataManager, a: usize, b: usize, ldpc: &dyn LDPC) {
    // Calculate S_aB_a with the result added to L
    manager.ensure_zero(&manager.data_id_range_of_ldpc_variable());
    for msg_id in 0..a {
        let adj_checks = ldpc.active_column(msg_id);
        let ldpc_ids = adj_checks.iter().map(|&id| manager.data_id_of_ldpc_variable(id)).collect::<Vec<_>>();
        //dbg!("LDPC encoder", &msg_id, &ldpc_ids);
        manager.broadcast_add(manager.data_id_of_active_variable(msg_id), &ldpc_ids);
    }
    // Calculate S_bB_b with the result added to L
    for msg_id in 0..b {
        let adj_checks = ldpc.inactive_column(msg_id);
        //dbg!("LDPC encoder", &msg_id, &adj_checks);
        let ldpc_ids = adj_checks.iter().map(|&id| manager.data_id_of_ldpc_variable(id)).collect::<Vec<_>>();
        //dbg!("LDPC encoder", &msg_id, &ldpc_ids);
        manager.broadcast_add(manager.data_id_of_inactive_variable(msg_id), &ldpc_ids);
    }
    //dbg!("LDPC encoder", manager.get_vector(21));
}*/

/// Solve (I'+D_sS_h)Y = X
#[allow(clippy::borrowed_box)]
fn hdpc_solve(
    manager: &mut DataManager,
    params: &CodeParams,
    hdpc: &Box<dyn HDPC>,
    ldpc: &Box<dyn LDPC>,
    hdpc_ids: &[usize],
) {
    // todo: test the calculation of D_s S_h as [D_s D_b] [S_h // 0]
    let ldpc_adj_check_inactive = |row: usize| {
        //if row < params.l {
        // keep only entries less than b
        ldpc.inactive_row(row)
            .iter()
            .filter(|&id| *id >= params.b)
            .map(|&id| id - params.b)
            .collect::<Vec<_>>()
        //} else {
        //    vec![]
        //}
    };

    let mut idssh = hdpc.mul_sparse_sh(manager.gf256(), params, &ldpc_adj_check_inactive);
    //dbg!("idssh", &idssh);

    for (i, row) in idssh.iter_mut().take(params.h).enumerate() {
        row[i] ^= 1;
    }
    let (p, r) = match manager.gf256() {
        Some(gf) => Matrix::lu_decomp(gf, &mut idssh),
        None => Matrix::lu_decomp(&GF2::new(), &mut idssh),
    };
    //dbg!("p", &p);
    //dbg!("idssh", &idssh);
    if r < params.h {
        panic!("The matrix I'+D_sS_h is not invertible, rank = {}", r);
    }

    let variable_ids_hdpc: Vec<usize> = manager.data_id_range_of_hdpc_variable();
    for i in 0..params.h {
        manager.move_to(hdpc_ids[p[i]], variable_ids_hdpc[i]);
    }

    lu_solve(manager, &mut idssh, &variable_ids_hdpc);
    // Substitute P into L = S_h P + L'
    for (var_col, &var_id) in variable_ids_hdpc.iter().take(params.h).enumerate() {
        let ids = ldpc.inactive_column(var_col + params.b);
        let ids = ids
            .iter()
            .map(|&id| manager.data_id_of_ldpc_variable(id))
            .collect::<Vec<_>>();
        manager.broadcast_add(var_id, &ids);
    }
}

/// Macro for LU solve operations that works with the data manager.
/// This macro implements forward and backward substitution for solving Ax = b in-place
/// where A is an LU-decomposed matrix and b are the target vectors.
fn lu_solve(manager: &mut DataManager, matrix_a: &mut [Vec<u8>], target_ids: &[usize]) {
    if matrix_a.len() != target_ids.len() {
        panic!("The number of rows in A must be equal to the number of target IDs");
    }

    let n = matrix_a.len();

    // Forward substitution (L part): solve Ly = b
    for j in 0..n - 1 {
        for i in j + 1..n {
            if matrix_a[i][j] != 0 {
                manager.mul_add(target_ids[j], matrix_a[i][j], target_ids[i]);
            }
        }
    }

    // Backward substitution (U part): solve Ux = y
    for j in (0..n).rev() {
        if matrix_a[j][j] == 0 {
            panic!(
                "Singular matrix: diagonal element at position {} is zero",
                j
            );
        }

        // Scale the current row by the inverse of the diagonal element
        manager.divide_scalar(matrix_a[j][j], target_ids[j]);

        // Subtract scaled row from previous rows
        for i in (0..j).rev() {
            if matrix_a[i][j] != 0 {
                manager.mul_add(target_ids[j], matrix_a[i][j], target_ids[i]);
            }
        }
    }
}
