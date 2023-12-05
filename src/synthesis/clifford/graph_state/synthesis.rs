use super::count::synthesize_graph_state_count;
use super::depth::synthesize_graph_state_depth;
/// This module contains the necessary methods to synthesize graph states
use crate::structures::clifford_circuit::CliffordCircuit;
use crate::structures::graph_state::GraphState;
use crate::structures::metric::Metric;

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

#[cfg(test)]
mod graph_state_tests {

    use super::*;
    use crate::structures::pauli_like::PauliLike;
    #[test]
    fn test_synthesis() {
        let n = 10;
        let gs = GraphState::random(n);
        let circuit = synthesize_graph_state(&gs, &Metric::COUNT, 100);
        let mut graph = GraphState::new(n);
        graph.conjugate_with_circuit(&circuit);
        assert_eq!(gs.adj, graph.adj);
    }
    #[test]
    fn test_synthesis_2() {
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
    fn test_synthesis_depth() {
        let n = 10;
        let gs = GraphState::random(n);
        let circuit_depth = synthesize_graph_state(&gs, &Metric::DEPTH, 100);
        let mut graph = GraphState::new(n);
        graph.conjugate_with_circuit(&circuit_depth);
        assert_eq!(gs.adj, graph.adj);
    }
}
