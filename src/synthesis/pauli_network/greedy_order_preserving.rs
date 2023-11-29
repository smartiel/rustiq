use super::greedy_pauli_network::{single_synthesis_step, Metric};
/// This module implements a data structure encoding the commutation DAG of a sequence of Pauli operators.
use crate::structures::clifford_circuit::CliffordCircuit;
use crate::structures::pauli_like::PauliLike;
use crate::structures::pauli_set::PauliSet;

use petgraph::prelude::*;
use pyo3::prelude::*;

pub type Dag = DiGraph<usize, ()>;

/// Constructs an anti-commutation Dag from a set of operators
fn build_dag_from_pauli_set(pauli_set: &PauliSet) -> Dag {
    let mut dag = Dag::new();
    let node_indices: Vec<NodeIndex> = (0..pauli_set.len()).map(|i| dag.add_node(i)).collect();

    for i in 0..pauli_set.len() {
        for j in 0..i {
            if !pauli_set.commute(i, j) {
                dag.add_edge(node_indices[i], node_indices[j], ());
            }
        }
    }
    return dag;
}

/// Computes the list of operators that can be synthesized
fn _get_front_layer(dag: &Dag) -> Vec<NodeIndex> {
    return dag
        .node_indices()
        .filter(|node| dag.neighbors(*node).collect::<Vec<_>>().len() == 0)
        .collect();
}

pub struct PauliDag {
    /// A global set containing all the operators
    pauli_set: PauliSet,
    /// The dag structure
    dag: Dag,
    /// A Pauli set containing only the front layer of the set
    /// (operators pairwise commute in this set)
    front_layer: PauliSet,
}

impl PauliDag {
    /// Constructs a PauliDag from a PauliSet
    pub fn from_pauli_set(pauli_set: PauliSet) -> Self {
        let nqbits = pauli_set.n;
        let dag = build_dag_from_pauli_set(&pauli_set);
        let indices_front_layer = _get_front_layer(&dag);

        let mut front_layer = PauliSet::new(nqbits);
        for index in indices_front_layer {
            let (phase, pstring) = pauli_set.get(*dag.node_weight(index).unwrap());
            front_layer.insert(&pstring, phase);
        }
        Self {
            pauli_set,
            dag,
            front_layer,
        }
    }

    /// Constructs a PauliDag from a slice of axes
    pub fn from_slice(axes: &[String]) -> Self {
        Self::from_pauli_set(PauliSet::from_slice(axes))
    }

    fn update_front_layer(&mut self) {
        // Popping the trivial operators from the front_layer
        self.front_layer.clear();
        // Removing nodes that were synthesized at this round
        loop {
            let node_count = self.dag.node_count();
            self.dag.retain_nodes(|graph, node_index| {
                if graph.first_edge(node_index, Direction::Outgoing) == None {
                    return self
                        .pauli_set
                        .support_size(*graph.node_weight(node_index).unwrap())
                        > 1;
                }
                return true;
            });
            if self.dag.node_count() == node_count {
                break;
            }
        }
        // Updating the front layer bucket
        for index in _get_front_layer(&self.dag) {
            let (phase, pstring) = self.pauli_set.get(*self.dag.node_weight(index).unwrap());
            self.front_layer.insert(&pstring, phase);
        }
    }

    /// Performs a single synthesis step
    fn single_step_synthesis(&mut self, metric: &Metric) -> CliffordCircuit {
        self.front_layer.support_size_sort();
        let circuit = single_synthesis_step(&mut self.front_layer, metric);
        // Updating the global set of operators
        self.pauli_set.conjugate_with_circuit(&circuit);
        // Updating the front layer
        self.update_front_layer();
        return circuit;
    }
}

/// A greedy synthesis method that preserves the sequence of operator up to swap of commuting operators.
#[pyfunction]
pub fn pauli_network_synthesis_no_permutation(
    axes: Vec<String>,
    metric: Metric,
) -> CliffordCircuit {
    let mut dag = PauliDag::from_slice(&axes);
    let mut circuit = CliffordCircuit::new(dag.pauli_set.n);
    while dag.front_layer.len() > 0 {
        circuit.extend_with(&dag.single_step_synthesis(&metric));
    }
    return circuit;
}

#[cfg(test)]
mod greedy_synthesis_tests {
    use super::*;
    use crate::structures::clifford_circuit::CliffordGate;
    use std::collections::HashSet;
    fn check_circuit(input: &[String], circuit: &CliffordCircuit) -> bool {
        let mut hit_map: HashSet<usize> = HashSet::new();
        let mut bucket = PauliSet::from_slice(input);
        for i in 0..bucket.len() {
            if bucket.support_size(i) == 1 {
                hit_map.insert(i);
            }
        }
        for gate in circuit.gates.iter() {
            match gate {
                CliffordGate::SqrtX(i) => {
                    bucket.sqrt_x(*i);
                }
                CliffordGate::SqrtXd(i) => {
                    bucket.sqrt_xd(*i);
                }
                CliffordGate::S(i) => {
                    bucket.s(*i);
                }
                CliffordGate::Sd(i) => {
                    bucket.sd(*i);
                }
                CliffordGate::H(i) => {
                    bucket.h(*i);
                }
                CliffordGate::CNOT(i, j) => {
                    bucket.cnot(*i, *j);
                }
            }
            for i in 0..bucket.len() {
                if bucket.support_size(i) == 1 {
                    hit_map.insert(i);
                }
            }
        }
        println!("Synthesized {} operators", hit_map.len());
        println!("{:?}", bucket);
        return hit_map.len() == input.len();
    }
    #[test]
    fn count_synthesis() {
        let input = vec!["XX".to_owned(), "ZZ".to_owned(), "YY".to_owned()];
        let circuit = pauli_network_synthesis_no_permutation(input.clone(), Metric::COUNT);
        println!("{circuit:?}");
        assert!(check_circuit(&input, &circuit))
    }
    #[test]
    fn depth_synthesis() {
        let input = vec!["XX".to_owned(), "ZZ".to_owned(), "YY".to_owned()];
        let circuit = pauli_network_synthesis_no_permutation(input.clone(), Metric::DEPTH);
        println!("{circuit:?}");
        assert!(check_circuit(&input, &circuit))
    }
}
