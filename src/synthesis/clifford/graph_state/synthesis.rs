use super::count::synthesize_graph_state_count;
use super::depth::synthesize_graph_state_depth;
/// This module contains the necessary methods to synthesize graph states
use crate::structures::{CliffordCircuit, CliffordGate, GraphState, Metric, PauliSet};

use super::utils::extract_graph_state;

pub fn synthesize_graph_state(
    graph: &GraphState,
    metric: &Metric,
    niter: usize,
) -> CliffordCircuit {
    match metric {
        Metric::DEPTH => synthesize_graph_state_depth(graph),
        Metric::COUNT => synthesize_graph_state_count(graph, niter),
    }
}

pub fn synthesize_stabilizer_state(
    stabilizers: &PauliSet,
    metric: &Metric,
    niter: usize,
) -> CliffordCircuit {
    // First, extracting the graph out of the stabilizers by making sure that the X part is full rank.
    let (graph, h_circuit) = extract_graph_state(stabilizers);
    let circuit = match metric {
        Metric::DEPTH => synthesize_graph_state_depth(&graph),
        Metric::COUNT => synthesize_graph_state_count(&graph, niter),
    };
    let mut output = CliffordCircuit::new(stabilizers.n);
    for i in 0..stabilizers.n {
        output.gates.push(CliffordGate::H(i));
    }
    output.extend_with(&circuit);
    output.extend_with(&h_circuit);

    output
}

#[cfg(test)]
mod stab_state_tests {

    use super::*;
    use crate::structures::pauli_like::PauliLike;
    use crate::structures::IsometryTableau;

    #[test]
    fn test_gs_synthesis_count() {
        let n = 10;
        let gs = GraphState::random(n);
        let circuit = synthesize_graph_state(&gs, &Metric::COUNT, 100);
        let mut graph = GraphState::new(n);
        graph.conjugate_with_circuit(&circuit);
        assert_eq!(gs.adj, graph.adj);
    }

    #[test]
    fn test_gs_problematic_example() {
        let n = 3;
        let gs = GraphState {
            adj: vec![
                vec![true, true, true],
                vec![true, true, true],
                vec![true, true, false],
            ],
            n: 3,
        };
        let circuit = synthesize_graph_state(&gs, &Metric::COUNT, 100);
        let mut graph = GraphState::new(n);
        graph.conjugate_with_circuit(&circuit);
        println!("{:?}", circuit);
        assert_eq!(gs.adj, graph.adj);
    }

    #[test]
    fn test_gs_synthesis_depth() {
        for _ in 0..10 {
            let n = 10;
            let gs = GraphState::random(n);
            let circuit_depth = synthesize_graph_state(&gs, &Metric::DEPTH, 100);
            let mut graph = GraphState::new(n);
            graph.conjugate_with_circuit(&circuit_depth);
            assert_eq!(gs.adj, graph.adj);
        }
    }
    #[test]
    fn test_stab_synthesis_depth() {
        for _ in 0..10 {
            let n = 20;
            let mut iso = IsometryTableau::random(0, n);
            iso.normalize_inplace();
            let circuit_depth = synthesize_stabilizer_state(&iso.stabilizers, &Metric::DEPTH, 100);
            let mut check_iso = IsometryTableau::new(0, n);
            check_iso.conjugate_with_circuit(&circuit_depth);
            check_iso.normalize_inplace();
            assert_eq!(check_iso, iso);
        }
    }
}
