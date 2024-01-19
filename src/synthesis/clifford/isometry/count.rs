use super::common::decompose;
use crate::routines::decoding::information_set_decoding;
use crate::routines::f2_linalg::{lu_facto, rowop, transpose, Matrix};
use crate::structures::{CliffordCircuit, CliffordGate, GraphState, IsometryTableau, PauliLike};

#[derive(Clone, Copy)]
enum Type {
    Cz(usize, usize),
    Cnot(usize, usize),
    SCnotS(usize, usize),
}

fn gather_parities(circuit: &CliffordCircuit, n: usize, k: usize) -> (Vec<Vec<bool>>, Vec<Type>) {
    let mut graph_state = GraphState::new(n);
    let mut b = vec![vec![false; k]; n];
    for i in 0..std::cmp::min(n, k) {
        b[i][i] = true;
    }
    let mut parities = Vec::new();
    let mut moves = Vec::new();
    let mut ei = vec![false; n + k];
    for i in 0..n {
        ei[i] = true;
        parities.push(ei.clone());
        ei[i] = false;
        moves.push(Type::Cz(0, i));
    }
    for i in 0..k {
        ei[i + n] = true;
        parities.push(ei.clone());
        ei[i + n] = false;
        moves.push(Type::Cnot(0, i));
    }
    for (index, gate) in circuit.gates.iter().enumerate() {
        graph_state.conjugate_with_gate(&gate);
        match gate {
            CliffordGate::CNOT(i, j) => {
                rowop(&mut b, *j, *i);
                for parity in parities.iter_mut() {
                    parity[*i] ^= parity[*j];
                }
                ei[*j] = true;
                parities.push(ei.clone());
                ei[*j] = false;
                moves.push(Type::Cz(index + 1, *j));

                let mut par = graph_state.adj[*i].clone();
                par.extend_from_slice(&b[*i]);
                parities.push(par);
                moves.push(Type::Cnot(index + 1, *i));

                graph_state.s(*i);
                let mut par = graph_state.adj[*i].clone();
                par.extend(&b[*i]);
                parities.push(par);
                graph_state.s(*i);
                moves.push(Type::SCnotS(index + 1, *i));
            }

            CliffordGate::CZ(i, j) => {
                let mut par = graph_state.adj[*i].clone();
                par.extend_from_slice(&b[*i]);
                parities.push(par);
                moves.push(Type::Cnot(index + 1, *i));

                graph_state.s(*i);
                let mut par = graph_state.adj[*i].clone();
                par.extend_from_slice(&b[*i]);
                parities.push(par);
                graph_state.s(*i);
                moves.push(Type::SCnotS(index + 1, *i));

                let mut par = graph_state.adj[*j].clone();
                for b in b[*j].iter() {
                    par.push(*b);
                }
                parities.push(par);
                moves.push(Type::Cnot(index + 1, *j));

                graph_state.s(*j);
                let mut par = graph_state.adj[*j].clone();
                for b in b[*j].iter() {
                    par.push(*b);
                }
                parities.push(par);
                graph_state.s(*j);
                moves.push(Type::SCnotS(index + 1, *j));
            }
            _ => {}
        }
    }

    return (parities, moves);
}

fn graph_state_and_b_synthesis(
    graph_adj: &Matrix,
    b_matrix: &Matrix,
    niter: usize,
) -> CliffordCircuit {
    let graph = GraphState::from_adj(graph_adj.clone());
    let mut circuit = CliffordCircuit::new(graph.n);
    let k = b_matrix.first().unwrap().len();
    for i in 0..graph.n {
        if i > 0 {
            let parity_len = i + std::cmp::min(i, k);
            let (mut parities, moves) = gather_parities(&circuit, i, std::cmp::min(i, k));
            let mut target = vec![false; parity_len];
            for j in 0..i {
                target[j] = graph.adj[i][j];
                if j < k {
                    target[j + i] = b_matrix[i][j];
                }
            }
            let solution = information_set_decoding(&mut parities, &mut target, niter);
            let mut new_circuit = CliffordCircuit::new(graph.n);
            let moves: Vec<Type> = solution
                .iter()
                .enumerate()
                .filter(|(_a, b)| **b)
                .map(|(a, _)| moves[a])
                .collect();
            for mov in moves.iter() {
                match mov {
                    Type::Cnot(gindex, qbit) => {
                        if *gindex == 0 {
                            new_circuit.gates.push(CliffordGate::CNOT(i, *qbit));
                        }
                    }
                    Type::Cz(gindex, qbit) => {
                        if *gindex == 0 {
                            new_circuit.gates.push(CliffordGate::CZ(i, *qbit));
                        }
                    }
                    Type::SCnotS(gindex, qbit) => {
                        if *gindex == 0 {
                            new_circuit.gates.push(CliffordGate::S(*qbit));
                            new_circuit.gates.push(CliffordGate::CNOT(i, *qbit));
                            new_circuit.gates.push(CliffordGate::S(*qbit));
                        }
                    }
                }
            }
            for k in 0..circuit.gates.len() {
                new_circuit.gates.push(circuit.gates[k].clone());
                for mov in moves.iter() {
                    match mov {
                        Type::Cnot(gindex, qbit) => {
                            if *gindex == k + 1 {
                                new_circuit.gates.push(CliffordGate::CNOT(i, *qbit));
                            }
                        }
                        Type::Cz(gindex, qbit) => {
                            if *gindex == k + 1 {
                                new_circuit.gates.push(CliffordGate::CZ(i, *qbit));
                            }
                        }
                        Type::SCnotS(gindex, qbit) => {
                            if *gindex == k + 1 {
                                new_circuit.gates.push(CliffordGate::S(*qbit));
                                new_circuit.gates.push(CliffordGate::CNOT(i, *qbit));
                                new_circuit.gates.push(CliffordGate::S(*qbit));
                            }
                        }
                    }
                }
            }
            circuit = new_circuit;
        }
    }

    let mut simulated = GraphState::new(graph.n);
    simulated.conjugate_with_circuit(&circuit);
    for i in 0..graph.n {
        if simulated.adj[i][i] != graph.adj[i][i] {
            circuit.gates.push(CliffordGate::S(i));
        }
    }
    circuit.dagger()
}

pub fn isometry_count_synthesis(isometry: &IsometryTableau, niter: usize) -> CliffordCircuit {
    let (g_k, g_n, b, h_circuit) = decompose(&isometry);
    let (l, u, _, ops) = lu_facto(&transpose(&b));
    let mut output = CliffordCircuit::new(isometry.n + isometry.k);
    let mut gn_as_gs = GraphState::from_adj(g_n);
    gn_as_gs.conjugate_with_circuit(&ops);
    let gn_circuit = graph_state_and_b_synthesis(&gn_as_gs.adj, &l, niter);

    let gk_circuit = graph_state_and_b_synthesis(&g_k, &transpose(&u), niter);
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
    #[test]
    fn test_graph_state_and_b_synthesis() {
        let mut rng = rand::thread_rng();
        let n = 20;
        let mut graph_adj = GraphState::random(n);
        let mut b_matrix: Matrix = vec![vec![false; n]; n];
        for i in 0..n {
            b_matrix[i][i] = true;
        }

        for _ in 0..n * n {
            let i = rng.gen::<usize>() % (n - 1);
            let j = rng.gen::<usize>() % (n - i - 1) + i + 1;
            assert!(j > i);
            rowop(&mut b_matrix, i, j);
        }
        let circuit = graph_state_and_b_synthesis(&graph_adj.adj, &b_matrix, 1);
        graph_adj.conjugate_with_circuit(&circuit);
        for gate in circuit.gates.iter() {
            match gate {
                CliffordGate::CNOT(i, j) => rowop(&mut b_matrix, *j, *i),
                _ => (),
            }
        }
        for i in 0..n {
            assert_eq!(b_matrix[i][i], true);
            for l in 0..n {
                if l != i {
                    assert_eq!(b_matrix[i][l], false);
                }
            }
        }
        for i in 0..n {
            for l in 0..n {
                assert_eq!(graph_adj.adj[i][l], false);
            }
        }
    }
    #[test]
    fn test_graph_state_and_b_synthesis_rectangular() {
        let mut rng = rand::thread_rng();
        let n = 20;
        let k = 5;
        let mut graph_adj = GraphState::random(n + k);
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
        let circuit = graph_state_and_b_synthesis(&graph_adj.adj, &b_matrix, 1);
        graph_adj.conjugate_with_circuit(&circuit);
        for gate in circuit.gates.iter() {
            match gate {
                CliffordGate::CNOT(i, j) => rowop(&mut b_matrix, *j, *i),
                _ => (),
            }
        }
        println!("=== After de-synthesis ===");
        println!("Graph:");
        print_matrix(&graph_adj.adj);
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
                assert_eq!(graph_adj.adj[i][l], false);
            }
        }
    }

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
    fn test_clifford_count_synthesis() {
        for _ in 0..20 {
            let n = 20;
            let k = 0;
            let isometry = IsometryTableau::random(n, k);
            let (ref_a, ref_b, ref_c, ref_d) = extract_abcd(&isometry);
            let circuit = isometry_count_synthesis(&isometry, 1);
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
    fn test_isometry_count_synthesis() {
        for _ in 0..20 {
            let n = 20;
            let k = 10;
            let mut isometry = IsometryTableau::random(n, k);
            let circuit = isometry_count_synthesis(&isometry, 1);
            isometry.normalize_inplace();
            let (ref_a, ref_b, ref_c, ref_d) = extract_abcd(&isometry);
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
