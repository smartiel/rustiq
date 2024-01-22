use super::greedy_pauli_network::{chunk_to_circuit, conjugate_with_chunk, ALL_CHUNKS};
use crate::structures::{
    hardware::SteinerTree, CliffordCircuit, CliffordGate, HardwareGraph, PauliSet,
};
use petgraph::algo::maximum_matching;
use petgraph::prelude::*;
use std::collections::HashMap;
struct SteinerData {
    data: Vec<(SteinerTree, Vec<usize>, usize)>,
    bucket: PauliSet,
}

impl SteinerData {
    fn from_bucket(bucket: PauliSet, graph: &HardwareGraph) -> Self {
        let mut data = Vec::new();
        for i in 0..bucket.len() {
            let support = bucket.get_support(i);
            let tree = graph.get_steiner_tree(&support);
            data.push((tree, support, i));
        }
        Self { data, bucket }
    }
    fn len(&self) -> usize {
        return self.data.len();
    }

    fn get_cost(&self) -> f64 {
        let mut cost = 0.;
        for (tree, _, _) in self.data.iter() {
            cost -= (0.5 as f64).powi(tree.cnot_cost() as i32);
        }
        cost
    }
    fn get_chunk_cost(&mut self, chunk: &[Option<CliffordGate>; 3], qbits: &[usize; 2]) -> f64 {
        self.conjugate_with(chunk, qbits, false);
        let cost = self.get_cost();
        self.conjugate_with(chunk, qbits, true);
        return cost;
    }
    fn conjugate_with(
        &mut self,
        chunk: &[Option<CliffordGate>; 3],
        qbits: &[usize; 2],
        reverse: bool,
    ) {
        conjugate_with_chunk(&mut self.bucket, chunk, qbits[0], qbits[1], reverse);
        for i in 0..self.data.len() {
            let new_support = self.bucket.get_support(self.data[i].2);
            self.data[i].0.update_tree(&new_support, qbits);
            self.data[i].1 = new_support;
        }
    }
    fn pop_trivial(&mut self) -> usize {
        let init_size = self.data.len();
        self.data.retain(|(_, support, _)| support.len() > 1);
        return init_size - self.data.len();
    }
}

fn build_graph(
    data: &mut SteinerData,
    hardware_graph: &HardwareGraph,
) -> (
    UnGraph<(), f64>,
    HashMap<(usize, usize), [Option<CliffordGate>; 3]>,
) {
    let mut graph: UnGraph<(), f64> = UnGraph::new_undirected();
    let mut best_chunks: HashMap<(usize, usize), [Option<CliffordGate>; 3]> = HashMap::new();
    for _ in 0..hardware_graph.len() {
        graph.add_node(());
    }
    for edge in hardware_graph.graph.raw_edges() {
        let qbit1 = edge.source().index();
        let qbit2 = edge.target().index();
        let mut best_chunk: Option<[Option<CliffordGate>; 3]> = None;
        let mut best_delta: Option<f64> = None;
        for chunk in ALL_CHUNKS.iter() {
            let cost_delta = data.get_chunk_cost(chunk, &[qbit1, qbit2]);
            if best_chunk.is_none() || cost_delta < best_delta.unwrap() {
                best_delta = Some(cost_delta);
                best_chunk = Some(chunk.clone());
            }
        }
        if let Some(delta) = best_delta {
            if delta < 0. {
                best_chunks.insert((qbit1, qbit2), best_chunk.unwrap());
                graph.add_edge(NodeIndex::new(qbit1), NodeIndex::new(qbit2), -delta);
            }
        }
    }
    return (graph, best_chunks);
}

fn single_synthesis_step_depth(
    data: &mut SteinerData,
    hardware_graph: &HardwareGraph,
) -> CliffordCircuit {
    let (graph, best_chunks) = build_graph(data, hardware_graph);
    let matching = maximum_matching(&graph);
    let mut circuit_piece = CliffordCircuit::new(hardware_graph.len());
    for (qbit1, qbit2) in matching.edges() {
        let chunk = best_chunks[&(qbit1.index(), qbit2.index())];
        circuit_piece.extend_with(&chunk_to_circuit(
            &chunk,
            qbit1.index(),
            qbit2.index(),
            hardware_graph.len(),
        ));
        data.conjugate_with(&chunk, &[qbit1.index(), qbit2.index()], false);
    }

    return circuit_piece;
}

pub fn aa_pauli_network_synthesis(
    bucket: PauliSet,
    couplings: &Vec<(usize, usize)>,
) -> CliffordCircuit {
    let hardware_graph = HardwareGraph::from_couplings(couplings);
    let mut output = CliffordCircuit::new(hardware_graph.len());
    let mut data = SteinerData::from_bucket(bucket, &hardware_graph);
    while data.len() > 0 {
        data.pop_trivial();
        let piece = single_synthesis_step_depth(&mut data, &hardware_graph);
        output.extend_with(&piece);
    }
    output
}
