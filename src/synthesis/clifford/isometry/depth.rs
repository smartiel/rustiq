use super::common::decompose;
use crate::routines::f2_linalg::{count_ones_except_diag, lu_facto, rowop, transpose, Matrix};
use crate::structures::{CliffordCircuit, CliffordGate, GraphState, IsometryTableau, PauliLike};
use crate::synthesis::clifford::graph_state::depth::get_czs;
/// Generates all the possible row operations that can be applied to matrix B in order to preserve its lower triangular structure
/// Here n is the number of columns and n+k is the number of rows
fn allowed_row_ops(n: usize, k: usize) -> Vec<(usize, usize)> {
    let mut moves = Vec::new();
    for i in 0..n + k {
        for j in 0..n + k {
            if i != j {
                if i > j {
                    moves.push((j, i));
                } else {
                    if i >= n {
                        moves.push((j, i));
                    }
                }
            }
        }
    }
    moves
}

fn score_matrix(
    graph: &mut GraphState,
    b_matrix: &mut Matrix,
    qubits_used: &Vec<bool>,
    moves: &Vec<(usize, usize)>,
) -> Vec<Vec<i32>> {
    let base_value: i32 = (graph.count_ones() + count_ones_except_diag(b_matrix)) as i32;
    let mut scores = vec![vec![-1; b_matrix.len()]; b_matrix.len()];
    for (j, i) in moves.iter() {
        if !qubits_used[*i] && !qubits_used[*j] {
            graph.cnot(*i, *j);
            rowop(b_matrix, *j, *i);
            scores[*i][*j] =
                base_value - (graph.count_ones() + count_ones_except_diag(b_matrix)) as i32;
            rowop(b_matrix, *j, *i);
            graph.cnot(*i, *j);
        }
    }
    return scores;
}

fn pick_best_operation(
    scores: &Vec<Vec<i32>>,
    moves: &Vec<(usize, usize)>,
) -> (i32, (usize, usize)) {
    let mut best_score = 0;
    let mut best_qubits: (usize, usize) = (0, 0);
    for (i, j) in moves.iter() {
        if scores[*j][*i] > best_score {
            best_score = scores[*j][*i];
            best_qubits = (*j, *i);
        }
    }
    return (best_score, best_qubits);
}
fn graph_state_and_b_synthesis(graph: &mut GraphState, b_matrix: &mut Matrix) -> CliffordCircuit {
    let row_ops = allowed_row_ops(
        b_matrix.first().unwrap().len(),
        b_matrix.len() - b_matrix.first().unwrap().len(),
    );
    let mut output = CliffordCircuit::new(b_matrix.len());
    for i in 0..graph.n {
        if graph.adj[i][i] {
            graph.s(i);
            output.gates.push(CliffordGate::S(i));
        }
    }
    while (graph.count_ones() + count_ones_except_diag(b_matrix)) > 0 {
        let mut qubits_used: Vec<bool> = vec![false; graph.n];
        loop {
            let scores = score_matrix(graph, b_matrix, &qubits_used, &row_ops);
            let (best_move, (control, target)) = pick_best_operation(&scores, &row_ops);
            if best_move == 0 {
                break;
            }
            graph.cnot(control, target);
            rowop(b_matrix, target, control);
            output.gates.push(CliffordGate::CNOT(control, target));
            qubits_used[control] = true;
            qubits_used[target] = true;
        }
        let cz_circuit = get_czs(&graph, &qubits_used);
        output.extend_with(&cz_circuit);
        graph.conjugate_with_circuit(&cz_circuit);
    }

    output
}

pub fn isometry_depth_synthesis(isometry: &IsometryTableau) -> CliffordCircuit {
    let (g_k, g_n, b, h_circuit) = decompose(&isometry);
    let (mut l, u, _, ops) = lu_facto(&transpose(&b));
    let mut output = CliffordCircuit::new(isometry.n + isometry.k);
    let mut gn_as_gs = GraphState::from_adj(g_n);
    gn_as_gs.conjugate_with_circuit(&ops);
    let mut gk_as_gs = GraphState::from_adj(g_k);

    // Synthesize gn + l
    let gn_circuit = graph_state_and_b_synthesis(&mut gn_as_gs, &mut l);
    // Synthesize gk + u^T
    let gk_circuit = graph_state_and_b_synthesis(&mut gk_as_gs, &mut transpose(&u));

    output.extend_with(&gk_circuit);
    for qbit in 0..isometry.n + isometry.k {
        output.gates.push(CliffordGate::H(qbit));
    }
    output.extend_with(&gn_circuit.dagger());
    output.extend_with(&ops.dagger());
    output.extend_with(&h_circuit);
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::routines::f2_linalg::Matrix;
    use crate::synthesis::clifford::isometry::common::extract_abcd;
    use rand::Rng;

    fn print_matrix(matrix: &Matrix) {
        for row in matrix.iter() {
            for elem in row.iter() {
                if *elem {
                    print!("1");
                } else {
                    print!("0");
                }
            }
            println!("");
        }
    }
    #[test]
    fn test_graph_state_and_b_synthesis_rectangular() {
        let mut rng = rand::thread_rng();
        let n = 5;
        let k = 5;
        let mut graph = GraphState::random(n + k);
        let mut b_matrix: Matrix = vec![vec![false; n]; n + k];
        for i in 0..n {
            b_matrix[i][i] = true;
        }

        for _ in 0..(n + k) * (n + k) {
            let i = rng.gen::<usize>() % (n + k - 1);
            let j = rng.gen::<usize>() % (n + k - i - 1) + i + 1;
            assert!(j > i);
            rowop(&mut b_matrix, i, j);
        }
        let circuit = graph_state_and_b_synthesis(&mut graph.clone(), &mut b_matrix.clone());
        graph.conjugate_with_circuit(&circuit);
        for gate in circuit.gates.iter() {
            match gate {
                CliffordGate::CNOT(i, j) => rowop(&mut b_matrix, *j, *i),
                _ => (),
            }
        }
        println!("=== After de-synthesis ===");
        println!("Graph:");
        print_matrix(&graph.adj);
        println!("B:");
        print_matrix(&b_matrix);
        for i in 0..n {
            assert_eq!(b_matrix[i][i], true);
            for l in 0..n {
                if l != i {
                    assert_eq!(b_matrix[i][l], false);
                }
            }
        }
        for i in n..n + k {
            for l in 0..n {
                assert_eq!(b_matrix[i][l], false);
            }
        }
        for i in 0..n + k {
            for l in 0..n + k {
                assert_eq!(graph.adj[i][l], false);
            }
        }
    }

    #[test]
    fn test_clifford_depth_synthesis() {
        for _ in 0..20 {
            let n = 10;
            let k = 0;
            let isometry = IsometryTableau::random(n, k);
            let (ref_a, ref_b, ref_c, ref_d) = extract_abcd(&isometry);
            let circuit = isometry_depth_synthesis(&isometry);
            let mut simulated = IsometryTableau::new(n, k);
            simulated.conjugate_with_circuit(&circuit);
            let (a, b, c, d) = extract_abcd(&simulated);

            assert_eq!(ref_a, a);
            assert_eq!(ref_b, b);
            assert_eq!(ref_c, c);
            assert_eq!(ref_d, d);
        }
    }

    #[test]
    fn test_isometry_depth_synthesis() {
        for _ in 0..20 {
            let n = 20;
            let k = 10;
            let mut isometry = IsometryTableau::random(n, k);
            isometry.normalize_inplace();
            let (ref_a, ref_b, ref_c, ref_d) = extract_abcd(&isometry);
            let circuit = isometry_depth_synthesis(&isometry);
            let mut simulated = IsometryTableau::new(n, k);
            simulated.conjugate_with_circuit(&circuit);
            simulated.normalize_inplace();
            let (a, b, c, d) = extract_abcd(&simulated);

            assert_eq!(ref_a, a);
            assert_eq!(ref_b, b);
            assert_eq!(ref_c, c);
            assert_eq!(ref_d, d);
        }
    }
}
