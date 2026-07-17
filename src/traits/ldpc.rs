// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

use crate::data_manager::DataManager;
use crate::types::CodeParams;

/// LDPC trait for the parity check matrices S_a and S_i.
pub trait LDPC {
    /// Non-zero entries for a column of S_a
    fn active_column(&self, _var_col: usize) -> Vec<usize> {
        vec![]
    }
    /// Non-zero entries for a row of S_a
    fn active_row(&self, _check_row: usize) -> Vec<usize> {
        vec![]
    }
    /// Non-zero entries for a column of S_i
    fn inactive_column(&self, _var_col: usize) -> Vec<usize> {
        vec![]
    }
    /// Non-zero entries for a row of S_i
    fn inactive_row(&self, _check_row: usize) -> Vec<usize> {
        vec![]
    }

    /// Encode the message vectors with the LDPC matrix with default encoding method
    fn encode(&self, manager: &mut DataManager, params: &CodeParams) {
        // Calculate S_aB_a with the result added to L
        manager.ensure_zero(&manager.data_id_range_of_ldpc_variable());
        for msg_id in 0..params.a {
            let adj_checks = self.active_column(msg_id);
            let ldpc_ids = adj_checks
                .iter()
                .map(|&id| manager.data_id_of_ldpc_variable(id))
                .collect::<Vec<_>>();
            //dbg!("LDPC encoder", &msg_id, &ldpc_ids);
            manager.broadcast_add_owned(manager.data_id_of_active_variable(msg_id), ldpc_ids);
        }
        // Calculate S_bB_b with the result added to L
        for msg_id in 0..params.b {
            let adj_checks = self.inactive_column(msg_id);
            //dbg!("LDPC encoder", &msg_id, &adj_checks);
            let ldpc_ids = adj_checks
                .iter()
                .map(|&id| manager.data_id_of_ldpc_variable(id))
                .collect::<Vec<_>>();
            //dbg!("LDPC encoder", &msg_id, &ldpc_ids);
            manager.broadcast_add_owned(manager.data_id_of_inactive_variable(msg_id), ldpc_ids);
        }
        //dbg!("LDPC encoder", manager.get_vector(21));
    }
}
