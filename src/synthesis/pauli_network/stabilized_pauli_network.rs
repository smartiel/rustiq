use petgraph::algo::maximum_matching;
use petgraph::prelude::*;
use std::collections::HashMap;

use super::greedy_pauli_network::{chunk_to_circuit, conjugate_with_chunk, Metric, ALL_CHUNKS};
use crate::structures::clifford_circuit::{CliffordCircuit, CliffordGate};
use crate::structures::pauli_set::PauliSet;

fn weight(vec: &Vec<bool>) -> usize {
    let n: usize = vec.len() / 2;
    let (z, x) = vec.split_at(n);
    return z
        .iter()
        .zip(x.iter())
        .map(|(a, b)| (a | b))
        .filter(|c| *c)
        .count();
}

fn xor_vec_bool(target: &mut Vec<bool>, control: &Vec<bool>) {
    for (a, b) in target.iter_mut().zip(control.iter()) {
        *a ^= b;
    }
}

fn weight_with_decoding(vec: &Vec<bool>, stabilizers: &Vec<Vec<bool>>) -> (usize, Vec<bool>) {
    let mut result = vec.clone();
    let mut witness = vec![false; stabilizers.len()];
    let mut w = weight(vec);
    loop {
        let old_w = w;
        for (index, stab) in stabilizers.iter().enumerate() {
            xor_vec_bool(&mut result, stab);
            let new_w = weight(&result);
            if new_w < w {
                witness[index] ^= true;
                w = new_w;
            } else {
                xor_vec_bool(&mut result, stab);
            }
        }
        if old_w == w {
            break;
        }
    }
    return (w, witness);
}

fn build_graph(
    bucket: &mut PauliSet,
    stabilizers: &mut PauliSet,
) -> (
    UnGraph<(), i32>,
    HashMap<(usize, usize), [Option<CliffordGate>; 3]>,
) {
    let mut graph: UnGraph<(), i32> = UnGraph::new_undirected();
    let mut best_chunks: HashMap<(usize, usize), [Option<CliffordGate>; 3]> = HashMap::new();
    for _ in 0..bucket.n {
        graph.add_node(());
    }
    for qbit1 in 0..bucket.n {
        for qbit2 in (qbit1 + 1)..bucket.n {
            // computing the initial identity count
            let init_id_count = bucket.count_id(qbit1) + bucket.count_id(qbit2);
            let mut max_score = 0;
            let mut best_chunk: [Option<CliffordGate>; 3] = [None; 3];
            for chunk in ALL_CHUNKS {
                // conjugating with the chunk
                conjugate_with_chunk(bucket, &chunk, qbit1, qbit2, false);
                let new_count = bucket.count_id(qbit1) + bucket.count_id(qbit2);
                let score: i32 = new_count as i32 - init_id_count as i32;
                if score > max_score {
                    max_score = score;
                    best_chunk = chunk.clone();
                }
                best_chunks.insert((qbit1, qbit2), best_chunk);
                // undoing the conjugation
                conjugate_with_chunk(bucket, &chunk, qbit1, qbit2, true);
            }
            // If there exists a chunk that improves the score, we add an edge labeled with the score in the graph
            if max_score > 0 {
                graph.add_edge(NodeIndex::new(qbit1), NodeIndex::new(qbit2), max_score);
            }
        }
    }
    return (graph, best_chunks);
}

fn single_synthesis_step_depth(
    bucket: &mut PauliSet,
    stabilizers: &mut PauliSet,
) -> CliffordCircuit {
    let (graph, best_chunks) = build_graph(bucket);
    let matching = maximum_matching(&graph);
    let mut circuit_piece = CliffordCircuit::new(bucket.n);
    for (qbit1, qbit2) in matching.edges() {
        let chunk = best_chunks[&(qbit1.index(), qbit2.index())];
        circuit_piece.extend_with(&chunk_to_circuit(
            &chunk,
            qbit1.index(),
            qbit2.index(),
            bucket.n,
        ));
        conjugate_with_chunk(bucket, &chunk, qbit1.index(), qbit2.index(), false);
    }

    return circuit_piece;
}
