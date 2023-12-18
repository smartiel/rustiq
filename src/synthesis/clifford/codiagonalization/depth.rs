use std::cmp::{max, min};

use super::common::{make_full_rank, permute_circuit};
use crate::routines::f2_linalg::{colop, diagonalize, rowop};
use crate::structures::{CliffordCircuit, CliffordGate, GraphState, Metric, PauliSet};
use crate::synthesis::clifford::graph_state::synthesize_graph_state;
use petgraph::algo::maximum_matching;
use petgraph::prelude::*;

fn score_matrix(
    table: &Vec<Vec<bool>>,
    qbits_used: &Vec<bool>,
    row: bool,
    rank: usize,
) -> Vec<Vec<i32>> {
    if row {
        let mut score_mat = vec![vec![0; table.len()]; table.len()];
        for i in 0..table.len() {
            for j in 0..table.len() {
                if i != j && !qbits_used[i] && !qbits_used[j] {
                    let score: i32 = table[i]
                        .iter()
                        .zip(table[j].iter())
                        .map(|(x_i, x_j)| (if *x_i { 1 } else { 0 }) * (if *x_j { 1 } else { -1 }))
                        .sum();
                    score_mat[i][j] = score;
                } else {
                    score_mat[i][j] = -1;
                }
            }
        }
        return score_mat;
    }
    let mut score_mat = vec![vec![0; rank]; rank];
    let k = rank;
    for i in 0..k {
        for j in 0..k {
            if i != j && !qbits_used[i] && !qbits_used[j] {
                let mut sum = 0;
                for m in 0..table.len() {
                    if table[m][i] {
                        sum += if table[m][j] { 1 } else { -1 };
                    }
                }
                score_mat[i][j] = sum;
            } else {
                score_mat[i][j] = -1;
            }
        }
    }
    return score_mat;
}

fn pick_best_operation(
    scores_rows: &Vec<Vec<i32>>,
    scores_cols: &Vec<Vec<i32>>,
) -> (bool, i32, (usize, usize)) {
    let mut best_score = -1;
    let mut best_control = 0;
    let mut best_target = 0;
    let mut is_row = true;
    for i in 0..scores_rows.len() {
        for j in 0..scores_rows.len() {
            if scores_rows[i][j] > best_score {
                best_control = i;
                best_target = j;
                best_score = scores_rows[i][j];
            }
        }
    }
    for i in 0..scores_cols.len() {
        for j in 0..scores_cols.len() {
            if scores_cols[i][j] > best_score {
                best_control = i;
                best_target = j;
                best_score = scores_cols[i][j];
                is_row = false;
            }
        }
    }
    return (is_row, best_score, (best_control, best_target));
}

fn has_ones(matrix: &Vec<Vec<bool>>) -> bool {
    for row in matrix.iter() {
        for elem in row.iter() {
            if *elem {
                return true;
            }
        }
    }
    return false;
}

fn run_matching(
    table: &Vec<Vec<bool>>,
    row_used: &Vec<bool>,
    col_used: &Vec<bool>,
) -> CliffordCircuit {
    let ncols = col_used.len();
    let nrows = row_used.len();
    let mut mgraph: UnGraph<(), i32> = UnGraph::new_undirected();
    for _ in 0..ncols + nrows {
        mgraph.add_node(());
    }
    for i in 0..ncols {
        for j in 0..nrows {
            if !col_used[i] && !row_used[j] && table[j][i] {
                mgraph.add_edge(NodeIndex::new(i), NodeIndex::new(ncols + j), 1);
            }
        }
    }
    let matching = maximum_matching(&mgraph);
    let mut circuit_piece = CliffordCircuit::new(nrows + ncols);
    for (qbit1, qbit2) in matching.edges() {
        let minq = min(qbit1.index(), qbit2.index());
        let maxq = max(qbit1.index(), qbit2.index());
        circuit_piece.gates.push(CliffordGate::CNOT(minq, maxq));
    }
    circuit_piece
}

fn reduce_x_part(pauli_set: &PauliSet) -> (CliffordCircuit, Vec<usize>, GraphState) {
    let (mut circuit, row_perm, rank, (mut z_table, mut x_table)) = make_full_rank(&pauli_set);
    let mut table: Vec<Vec<bool>> = x_table.iter().skip(rank).map(|r| r.clone()).collect();
    let mut cnot_circuit = CliffordCircuit::new(circuit.nqbits);
    while has_ones(&table) {
        let mut used_rows = vec![false; table.len()];
        let mut used_cols = vec![false; rank];
        loop {
            let score_rows = score_matrix(&table, &used_rows, true, rank);
            let score_cols = score_matrix(&table, &used_cols, false, rank);
            let (is_row, score, (i, j)) = pick_best_operation(&score_rows, &score_cols);
            if score <= 0 {
                break;
            }
            if is_row {
                cnot_circuit
                    .gates
                    .push(CliffordGate::CNOT(rank + i, rank + j));
                rowop(&mut table, i, j);
                used_rows[i] = true;
                used_rows[j] = true;
            } else {
                cnot_circuit.gates.push(CliffordGate::CNOT(i, j));
                colop(&mut table, j, i);
                used_cols[i] = true;
                used_cols[j] = true;
            }
        }
        let piece = run_matching(&table, &used_rows, &used_cols);
        cnot_circuit.extend_with(&piece);
        for gate in piece.gates.iter() {
            if let CliffordGate::CNOT(i, j) = gate {
                table[j - rank][*i] ^= true;
            }
        }
    }
    // update X & Z table in order to extract a graph state
    for gate in cnot_circuit.gates.iter() {
        if let CliffordGate::CNOT(i, j) = gate {
            rowop(&mut z_table, *j, *i);
            rowop(&mut x_table, *i, *j);
        }
    }
    diagonalize(&mut x_table, &mut z_table, rank);
    let mut graph_state = GraphState::new(rank);
    for col in 0..rank {
        for row in 0..rank {
            graph_state.adj[row][col] = z_table[row][col];
        }
    }
    for col in 0..rank {
        for row in 0..rank {
            assert_eq!(graph_state.adj[row][col], graph_state.adj[col][row]);
        }
    }

    let permuted_cnots = permute_circuit(&cnot_circuit, &row_perm);
    circuit.extend_with(&permuted_cnots);
    return (circuit, row_perm, graph_state);
}

pub fn codiagonalize_depth(pauli_set: &PauliSet) -> CliffordCircuit {
    let (mut circuit, perm, graph) = reduce_x_part(pauli_set);
    let gs_synth = synthesize_graph_state(&graph, &Metric::DEPTH, 0);
    let gs_synth = permute_circuit(&gs_synth, &perm);
    circuit.extend_with(&gs_synth.dagger());
    for i in 0..graph.n {
        circuit.gates.push(CliffordGate::H(perm[i]));
    }
    return circuit;
}

#[cfg(test)]
mod codiag_depth_tests {

    use super::*;
    use crate::structures::pauli_like::PauliLike;

    use rand::Rng;
    fn random_instance(n: usize, m: usize) -> PauliSet {
        let mut rng = rand::thread_rng();
        let mut pset = PauliSet::new(n);
        for _ in 0..m {
            let mut vec: Vec<bool> = Vec::new();
            for _ in 0..n {
                vec.push(rng.gen::<bool>());
            }
            for _ in 0..n {
                vec.push(false);
            }
            pset.insert_vec_bool(&vec, false);
        }
        for _ in 0..n * n {
            let i = rng.gen::<usize>() % n;
            loop {
                let j = rng.gen::<usize>() % n;
                if j != i {
                    pset.cnot(i, j);
                    let g2 = rng.gen::<bool>();
                    if g2 {
                        pset.h(j);
                    } else {
                        pset.s(j);
                    }
                    break;
                }
            }
            let g1 = rng.gen::<bool>();
            if g1 {
                pset.h(i);
            } else {
                pset.s(i);
            }
        }
        pset
    }

    #[test]
    fn test_thin() {
        for _ in 0..10 {
            let instance = random_instance(50, 20);
            let mut copy_instance = instance.clone();
            let circuit = codiagonalize_depth(&instance);
            copy_instance.conjugate_with_circuit(&circuit);
            for i in 0..instance.len() {
                let (_, vec) = copy_instance.get_as_vec_bool(i);
                for j in 0..instance.n {
                    assert!(!vec[j]);
                }
            }
        }
    }
    #[test]
    fn test_thick() {
        for _ in 0..10 {
            let instance = random_instance(20, 50);
            let mut copy_instance = instance.clone();
            let circuit = codiagonalize_depth(&instance);
            copy_instance.conjugate_with_circuit(&circuit);
            for i in 0..instance.len() {
                let (_, vec) = copy_instance.get_as_vec_bool(i);
                for j in 0..instance.n {
                    assert!(!vec[j]);
                }
            }
        }
    }

    #[test]
    fn test_small() {
        for _ in 0..10 {
            let instance = random_instance(10, 6);
            let mut copy_instance = instance.clone();
            let circuit = codiagonalize_depth(&instance);
            copy_instance.conjugate_with_circuit(&circuit);
            for i in 0..instance.len() {
                let (_, vec) = copy_instance.get_as_vec_bool(i);
                for j in 0..instance.n {
                    assert!(!vec[j]);
                }
            }
        }
    }
}
