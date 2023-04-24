use petgraph::algo::maximum_matching;
use petgraph::prelude::*;

use std::collections::HashMap;

use super::circuit::{Circuit, Gate};
use super::pauli_set::PauliSet;

#[derive(Debug)]
pub enum Metric {
    COUNT,
    DEPTH,
}

const ALL_CHUNKS: [[Option<Gate>; 3]; 18] = [
    [None, None, Some(Gate::CNOT(0, 1))],
    [None, None, Some(Gate::CNOT(1, 0))],
    [None, Some(Gate::S(1)), Some(Gate::CNOT(0, 1))],
    [None, Some(Gate::S(1)), Some(Gate::CNOT(1, 0))],
    [None, Some(Gate::H(1)), Some(Gate::CNOT(0, 1))],
    [None, Some(Gate::H(1)), Some(Gate::CNOT(1, 0))],
    [Some(Gate::S(0)), None, Some(Gate::CNOT(0, 1))],
    [Some(Gate::S(0)), None, Some(Gate::CNOT(1, 0))],
    [Some(Gate::S(0)), Some(Gate::S(1)), Some(Gate::CNOT(0, 1))],
    [Some(Gate::S(0)), Some(Gate::S(1)), Some(Gate::CNOT(1, 0))],
    [Some(Gate::S(0)), Some(Gate::H(1)), Some(Gate::CNOT(0, 1))],
    [Some(Gate::S(0)), Some(Gate::H(1)), Some(Gate::CNOT(1, 0))],
    [Some(Gate::H(0)), None, Some(Gate::CNOT(0, 1))],
    [Some(Gate::H(0)), None, Some(Gate::CNOT(1, 0))],
    [Some(Gate::H(0)), Some(Gate::S(1)), Some(Gate::CNOT(0, 1))],
    [Some(Gate::H(0)), Some(Gate::S(1)), Some(Gate::CNOT(1, 0))],
    [Some(Gate::H(0)), Some(Gate::H(1)), Some(Gate::CNOT(0, 1))],
    [Some(Gate::H(0)), Some(Gate::H(1)), Some(Gate::CNOT(1, 0))],
];

fn chunk_to_circuit(
    chunk: &[Option<Gate>; 3],
    qbit1: usize,
    qbit2: usize,
    nqbits: usize,
) -> Circuit {
    let mut circuit_piece = Circuit::new(nqbits);
    for gate in chunk {
        match gate {
            Some(Gate::S(i)) => {
                if *i == 0 {
                    circuit_piece.gates.push(Gate::S(qbit1));
                } else {
                    circuit_piece.gates.push(Gate::S(qbit2));
                }
            }
            Some(Gate::H(i)) => {
                if *i == 0 {
                    circuit_piece.gates.push(Gate::H(qbit1));
                } else {
                    circuit_piece.gates.push(Gate::H(qbit2));
                }
            }
            Some(Gate::CNOT(i, _)) => {
                if *i == 0 {
                    circuit_piece.gates.push(Gate::CNOT(qbit1, qbit2));
                } else {
                    circuit_piece.gates.push(Gate::CNOT(qbit2, qbit1));
                }
            }
            _ => {}
        }
    }
    return circuit_piece;
}

fn conjugate_with_chunk(
    bucket: &mut PauliSet,
    chunk: &[Option<Gate>; 3],
    qbit1: usize,
    qbit2: usize,
    reverse: bool,
) {
    if reverse {
        for gate in chunk.iter().rev() {
            match gate {
                Some(Gate::S(i)) => {
                    if *i == 0 {
                        bucket.s(qbit1);
                    } else {
                        bucket.s(qbit2);
                    }
                }
                Some(Gate::H(i)) => {
                    if *i == 0 {
                        bucket.h(qbit1);
                    } else {
                        bucket.h(qbit2);
                    }
                }
                Some(Gate::CNOT(i, _)) => {
                    if *i == 0 {
                        bucket.cnot(qbit1, qbit2);
                    } else {
                        bucket.cnot(qbit2, qbit1);
                    }
                }
                _ => {}
            }
        }
    } else {
        for gate in chunk.iter() {
            match gate {
                Some(Gate::S(i)) => {
                    if *i == 0 {
                        bucket.s(qbit1);
                    } else {
                        bucket.s(qbit2);
                    }
                }
                Some(Gate::H(i)) => {
                    if *i == 0 {
                        bucket.h(qbit1);
                    } else {
                        bucket.h(qbit2);
                    }
                }
                Some(Gate::CNOT(i, _)) => {
                    if *i == 0 {
                        bucket.cnot(qbit1, qbit2);
                    } else {
                        bucket.cnot(qbit2, qbit1);
                    }
                }
                _ => {}
            }
        }
    }
}

fn single_synthesis_step_count(bucket: &mut PauliSet) -> Circuit {
    let mut max_score = -1;
    let mut best_chunk: [Option<Gate>; 3] = [None; 3];
    let mut best_args: [usize; 2] = [0, 0];
    for qbit1 in 0..bucket.n {
        for qbit2 in 0..qbit1 {
            // computing the initial identity count
            let init_id_count = bucket.count_id(qbit1) + bucket.count_id(qbit2);
            for chunk in ALL_CHUNKS {
                // conjugating with the chunk
                conjugate_with_chunk(bucket, &chunk, qbit1, qbit2, false);
                let new_count = bucket.count_id(qbit1) + bucket.count_id(qbit2);
                let score: i32 = new_count as i32 - init_id_count as i32;
                if score > max_score {
                    max_score = score;
                    best_chunk = chunk.clone();
                    best_args = [qbit1, qbit2];
                }
                conjugate_with_chunk(bucket, &chunk, qbit1, qbit2, true);
            }
        }
    }
    conjugate_with_chunk(bucket, &best_chunk, best_args[0], best_args[1], false);

    return chunk_to_circuit(&best_chunk, best_args[0], best_args[1], bucket.n);
}

fn build_graph(
    bucket: &mut PauliSet,
) -> (UnGraph<(), i32>, HashMap<(usize, usize), [Option<Gate>; 3]>) {
    let mut graph: UnGraph<(), i32> = UnGraph::new_undirected();
    let mut best_chunks: HashMap<(usize, usize), [Option<Gate>; 3]> = HashMap::new();
    for _ in 0..bucket.n {
        graph.add_node(());
    }
    for qbit1 in 0..bucket.n {
        for qbit2 in (qbit1 + 1)..bucket.n {
            // computing the initial identity count
            let init_id_count = bucket.count_id(qbit1) + bucket.count_id(qbit2);
            let mut max_score = -1;
            let mut best_chunk: [Option<Gate>; 3] = [None; 3];
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
                conjugate_with_chunk(bucket, &chunk, qbit1, qbit2, true);
            }
            graph.add_edge(NodeIndex::new(qbit1), NodeIndex::new(qbit2), max_score);
        }
    }
    return (graph, best_chunks);
}

fn single_synthesis_step_depth(bucket: &mut PauliSet) -> Circuit {
    let (graph, best_chunks) = build_graph(bucket);
    let matching = maximum_matching(&graph);
    let mut circuit_piece = Circuit::new(bucket.n);

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

fn single_synthesis_step(bucket: &mut PauliSet, metric: &Metric) -> Circuit {
    return match metric {
        Metric::COUNT => single_synthesis_step_count(bucket),
        Metric::DEPTH => single_synthesis_step_depth(bucket),
    };
}

/// Builds a Pauli network for a collection of Pauli operators
pub fn pauli_network_synthesis(axes: &[String], metric: &Metric) -> Circuit {
    if axes.is_empty() {
        return Circuit::new(0);
    }
    let nqbits = axes[0].len();
    let mut bucket = PauliSet::new(nqbits);
    let mut output = Circuit::new(nqbits);
    for axis in axes {
        bucket.insert(axis);
    }
    loop {
        bucket.support_size_sort();
        while bucket.support_size(0) == 1 {
            bucket.pop();
        }
        if bucket.len() == 0 {
            break;
        }
        let circuit_piece = single_synthesis_step(&mut bucket, &metric);
        output.extend_with(&circuit_piece);
    }
    return output;
}

#[cfg(test)]
mod greedy_synthesis_tests {
    use super::*;
    use std::collections::HashSet;
    fn check_circuit(input: &[String], circuit: &Circuit) -> bool {
        let mut hit_map: HashSet<usize> = HashSet::new();
        let mut bucket = PauliSet::from_slice(circuit.nqbits, input);
        for i in 0..bucket.len() {
            if bucket.support_size(i) == 1 {
                hit_map.insert(i);
            }
        }
        for gate in circuit.gates.iter() {
            match gate {
                Gate::S(i) => {
                    bucket.s(*i);
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
        return hit_map.len() == input.len();
    }

    #[test]
    fn count_synthesis() {
        let input = ["XX".to_owned(), "ZZ".to_owned(), "YY".to_owned()];
        let circuit = pauli_network_synthesis(&input, &Metric::COUNT);
        assert!(check_circuit(&input, &circuit))
    }
    #[test]
    fn depth_synthesis() {
        let input = ["XX".to_owned(), "ZZ".to_owned(), "YY".to_owned()];
        let circuit = pauli_network_synthesis(&input, &Metric::DEPTH);
        assert!(check_circuit(&input, &circuit))
    }
}
