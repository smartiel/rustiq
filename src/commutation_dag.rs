/// This module implements a data structure encoding the commutation DAG of a sequence of Pauli operators.
use super::circuit::Circuit;
use super::greedy_pauli_network::{single_synthesis_step, Metric};
use super::pauli_set::PauliSet;
use petgraph::prelude::*;

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

fn _get_front_layer(dag: &Dag) -> Vec<usize> {
    return dag
        .node_indices()
        .filter(|node| dag.neighbors(*node).collect::<Vec<_>>().len() == 0)
        .map(|node| node.index())
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
            let (phase, pstring) = pauli_set.get(index);
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
    /// Performs a single synthesis step
    fn single_step_synthesis(&mut self, metric: &Metric) -> Circuit {
        let circuit = single_synthesis_step(&mut self.front_layer, metric);
        // Updating the global set of operators
        self.pauli_set.conjugate_with_circuit(&circuit);
        // Popping the trivial operators from the front_layer
        let current_front_layer = _get_front_layer(&self.dag);
        while self.front_layer.support_size(0) == 1 {
            let wit = self.front_layer.get(0);
            self.front_layer.pop();
            for index in current_front_layer.iter() {
                if self.pauli_set.get(*index) == wit {
                    self.dag.remove_node(NodeIndex::new(*index));
                    break;
                }
            }
        }
        for index in _get_front_layer(&self.dag) {
            let (phase, pstring) = self.pauli_set.get(index);
            self.front_layer.insert(&pstring, phase);
        }
        return circuit;
    }
}

/// A greedy synthesis method that preserves the sequence of operator up to swap of commuting operators.
pub fn pauli_network_synthesis_no_permutation(axes: &[String], metric: &Metric) -> Circuit {
    let mut dag = PauliDag::from_slice(axes);
    let mut circuit = Circuit::new(dag.pauli_set.n);
    while dag.dag.node_count() > 0 && dag.front_layer.len() > 0 {
        let piece = dag.single_step_synthesis(metric);
        circuit.extend_with(&piece);
    }
    return circuit;
}

#[cfg(test)]
mod greedy_synthesis_tests {
    use super::super::circuit::Gate;
    use super::*;
    use std::collections::HashSet;
    fn check_circuit(input: &[String], circuit: &Circuit) -> bool {
        let mut hit_map: HashSet<usize> = HashSet::new();
        let mut bucket = PauliSet::from_slice(input);
        for i in 0..bucket.len() {
            if bucket.support_size(i) == 1 {
                hit_map.insert(i);
            }
        }
        for gate in circuit.gates.iter() {
            match gate {
                Gate::SqrtX(i) => {
                    bucket.sqrt_x(*i);
                }
                Gate::SqrtXd(i) => {
                    bucket.sqrt_xd(*i);
                }
                Gate::S(i) => {
                    bucket.s(*i);
                }
                Gate::Sd(i) => {
                    bucket.sd(*i);
                }
                Gate::H(i) => {
                    bucket.h(*i);
                }
                Gate::CNOT(i, j) => {
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
}
