use crate::structures::clifford_circuit::{CliffordCircuit, CliffordGate};
use crate::structures::graph_state::GraphState;
use crate::structures::pauli_like::PauliLike;
use petgraph::algo::maximum_matching;
use petgraph::prelude::*;

fn score_matrix(graph: &mut GraphState, qubits_used: &Vec<bool>) -> Vec<Vec<i32>> {
    let base_value: i32 = graph.count_ones() as i32;
    let mut scores = vec![vec![-1; graph.n]; graph.n];
    for i in 0..graph.n {
        for j in i + 1..graph.n {
            if !qubits_used[i] && !qubits_used[j] {
                graph.cnot(i, j);
                scores[i][j] = base_value - graph.count_ones() as i32;
                graph.cnot(i, j);
            }
        }
    }
    scores
}

fn pick_best_operation(scores: &Vec<Vec<i32>>) -> (i32, (usize, usize)) {
    let mut best_score = 0;
    let mut best_qubits: (usize, usize) = (0, 0);
    for i in 0..scores.len() {
        for j in 0..scores.len() {
            if scores[i][j] > best_score {
                best_score = scores[i][j];
                best_qubits = (i, j);
            }
        }
    }
    return (best_score, best_qubits);
}

pub fn get_czs(graph: &GraphState, qubits_used: &Vec<bool>) -> CliffordCircuit {
    let mut mgraph: UnGraph<(), i32> = UnGraph::new_undirected();
    for _ in 0..graph.n {
        mgraph.add_node(());
    }
    for i in 0..graph.n {
        for j in 0..graph.n {
            if !qubits_used[i] && !qubits_used[j] && graph.adj[i][j] {
                mgraph.add_edge(NodeIndex::new(i), NodeIndex::new(j), 1);
            }
        }
    }
    let matching = maximum_matching(&mgraph);
    let mut circuit_piece = CliffordCircuit::new(graph.n);
    for (qbit1, qbit2) in matching.edges() {
        circuit_piece
            .gates
            .push(CliffordGate::CZ(qbit1.index(), qbit2.index()));
    }
    circuit_piece
}

pub fn synthesize_graph_state_depth(input_graph: &GraphState) -> CliffordCircuit {
    let mut graph = input_graph.clone();
    let mut circuit = CliffordCircuit::new(graph.n);
    for i in 0..graph.n {
        if graph.adj[i][i] {
            graph.s(i);
            circuit.gates.push(CliffordGate::S(i));
        }
    }
    while graph.count_ones() > 0 {
        let mut qubits_used: Vec<bool> = vec![false; graph.n];
        loop {
            let scores = score_matrix(&mut graph, &qubits_used);
            let (best_move, (control, target)) = pick_best_operation(&scores);
            if best_move == 0 {
                break;
            }
            graph.cnot(control, target);
            circuit.gates.push(CliffordGate::CNOT(control, target));
            qubits_used[control] = true;
            qubits_used[target] = true;
        }
        let cz_circuit = get_czs(&graph, &qubits_used);
        circuit.extend_with(&cz_circuit);
        graph.conjugate_with_circuit(&cz_circuit);
    }

    return circuit.dagger();
}
