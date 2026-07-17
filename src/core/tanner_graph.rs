// Copyright (c) 2025 Shenghao Yang. All rights reserved.
// Licensed under AGPL-3.0 or commercial license. See LICENSE for details.

//! Lightweight Tanner graph for trial BP decoding (doc sec:inac_strategies).
//!
//! Contains only active variable nodes and check nodes (unused equations).
//! Used to run trial BP without computing inactive coefficients, e.g. to evaluate
//! |N(v)| = number of variables that can be solved by inactivating v.

use std::collections::HashMap;

/// Tanner graph over active variables and unused equations only.
/// Variable indices inside the graph are *local* (0..num_active); use
/// `active_var_ids()` to map to global MasterSystem var IDs.
#[derive(Clone, Debug)]
pub struct TannerGraphDegree2 {
    /// Global variable IDs; index = local index.
    active_var_ids: Vec<usize>,
    /// For each check, list of local var indices (adjacent variables).
    check_adjacent: Vec<Vec<usize>>,
    /// For each variable (local index), list of check indices it belongs to.
    var_adjacent: Vec<Vec<usize>>,
}

impl TannerGraphDegree2 {
    /// Build from list of active variable IDs (global) and unused equations.
    /// Each equation is a list of global var IDs that appear in that equation.
    pub fn from_active_and_equations(active_var_ids: &[usize], equations: &[Vec<usize>]) -> Self {
        let n_var = active_var_ids.len();
        let global_to_local: HashMap<usize, usize> = active_var_ids
            .iter()
            .enumerate()
            .map(|(i, &v)| (v, i))
            .collect();

        let mut var_adjacent = vec![Vec::new(); n_var];
        let check_adjacent: Vec<Vec<usize>> = equations
            .iter()
            .map(|eq_vars| {
                eq_vars
                    .iter()
                    .filter_map(|&v| global_to_local.get(&v).copied())
                    .collect()
            })
            .collect();

        for (c, adj) in check_adjacent.iter().enumerate() {
            for &v_local in adj {
                var_adjacent[v_local].push(c);
            }
        }

        Self {
            active_var_ids: active_var_ids.to_vec(),
            check_adjacent,
            var_adjacent,
        }
    }

    fn remove_variable(&mut self, var_local: usize) -> usize {
        let mut solved = 1;
        let var_adj = self.var_adjacent[var_local].clone();
        self.var_adjacent[var_local].clear();
        for &c in &var_adj {
            if let Some(pos) = self.check_adjacent[c].iter().position(|&x| x == var_local) {
                self.check_adjacent[c].swap_remove(pos);
            }
            if self.check_adjacent[c].len() == 1 {
                let other_var = self.check_adjacent[c][0];
                solved += self.remove_variable(other_var);
            }
        }
        solved
    }

    /// Run BP with inactivation of variable `var_local` until no degree-1 check remains. Returns the number of
    /// variables that would be solved (peeled). Mutates the graph.
    pub fn run_bp_with_inactivation(&mut self, var_local: usize) -> usize {
        self.remove_variable(var_local)
    }

    /// Best inactivation by trial run: for each variable with degree > 0,
    /// run BP with that variable inactivated; return the global var_id
    /// that maximizes |N(v)|.
    pub fn best_inactivation_by_trial_run(&mut self) -> Option<usize> {
        let mut best_var_id = None;
        let mut best_count = 0_usize;
        let num_vars = self.var_adjacent.len();
        for v_local in 0..num_vars {
            if self.var_adjacent[v_local].is_empty() {
                continue;
            }
            let count = self.run_bp_with_inactivation(v_local);
            if count > best_count {
                best_count = count;
                best_var_id = Some(self.active_var_ids[v_local]);
            }
        }
        best_var_id
    }

    /// Size
    pub fn size(&self) -> usize {
        self.check_adjacent.len()
    }
}
