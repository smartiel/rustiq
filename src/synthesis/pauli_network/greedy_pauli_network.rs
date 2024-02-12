use petgraph::algo::maximum_matching;
use petgraph::prelude::*;
use std::collections::HashMap;

use crate::structures::clifford_circuit::{CliffordCircuit, CliffordGate};
use crate::structures::metric::Metric;
use crate::structures::pauli_like::PauliLike;
use crate::structures::pauli_set::PauliSet;

use super::chunks::{Chunk, ALL_CHUNKS};

pub fn chunk_to_circuit(
    chunk: &Chunk,
    qbit1: usize,
    qbit2: usize,
    nqbits: usize,
) -> CliffordCircuit {
    let mut circuit_piece = CliffordCircuit::new(nqbits);
    for gate in chunk {
        match gate {
            Some(CliffordGate::S(i)) => {
                if *i == 0 {
                    circuit_piece.gates.push(CliffordGate::S(qbit1));
                } else {
                    circuit_piece.gates.push(CliffordGate::S(qbit2));
                }
            }
            Some(CliffordGate::Sd(i)) => {
                if *i == 0 {
                    circuit_piece.gates.push(CliffordGate::Sd(qbit1));
                } else {
                    circuit_piece.gates.push(CliffordGate::Sd(qbit2));
                }
            }
            Some(CliffordGate::SqrtX(i)) => {
                if *i == 0 {
                    circuit_piece.gates.push(CliffordGate::SqrtX(qbit1));
                } else {
                    circuit_piece.gates.push(CliffordGate::SqrtX(qbit2));
                }
            }
            Some(CliffordGate::SqrtXd(i)) => {
                if *i == 0 {
                    circuit_piece.gates.push(CliffordGate::SqrtXd(qbit1));
                } else {
                    circuit_piece.gates.push(CliffordGate::SqrtXd(qbit2));
                }
            }
            Some(CliffordGate::H(i)) => {
                if *i == 0 {
                    circuit_piece.gates.push(CliffordGate::H(qbit1));
                } else {
                    circuit_piece.gates.push(CliffordGate::H(qbit2));
                }
            }
            Some(CliffordGate::CNOT(i, _)) => {
                if *i == 0 {
                    circuit_piece.gates.push(CliffordGate::CNOT(qbit1, qbit2));
                } else {
                    circuit_piece.gates.push(CliffordGate::CNOT(qbit2, qbit1));
                }
            }
            _ => {}
        }
    }
    return circuit_piece;
}

pub fn conjugate_with_chunk(
    bucket: &mut PauliSet,
    chunk: &Chunk,
    qbit1: usize,
    qbit2: usize,
    reverse: bool,
) {
    if reverse {
        for gate in chunk.iter().rev() {
            match gate {
                Some(CliffordGate::S(i)) => {
                    if *i == 0 {
                        bucket.s(qbit1);
                    } else {
                        bucket.s(qbit2);
                    }
                }
                Some(CliffordGate::H(i)) => {
                    if *i == 0 {
                        bucket.h(qbit1);
                    } else {
                        bucket.h(qbit2);
                    }
                }
                Some(CliffordGate::SqrtX(i)) => {
                    if *i == 0 {
                        bucket.sqrt_x(qbit1);
                    } else {
                        bucket.sqrt_x(qbit2);
                    }
                }
                Some(CliffordGate::CNOT(i, _)) => {
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
                Some(CliffordGate::S(i)) => {
                    if *i == 0 {
                        bucket.s(qbit1);
                    } else {
                        bucket.s(qbit2);
                    }
                }
                Some(CliffordGate::SqrtX(i)) => {
                    if *i == 0 {
                        bucket.sqrt_x(qbit1);
                    } else {
                        bucket.sqrt_x(qbit2);
                    }
                }
                Some(CliffordGate::H(i)) => {
                    if *i == 0 {
                        bucket.h(qbit1);
                    } else {
                        bucket.h(qbit2);
                    }
                }
                Some(CliffordGate::CNOT(i, _)) => {
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

fn single_synthesis_step_count(bucket: &mut PauliSet) -> CliffordCircuit {
    let mut max_score = -1;
    let mut best_chunk: Chunk = [None; 3];
    let mut best_args: [usize; 2] = [0, 0];
    let support = bucket.get_support(0);
    for i1 in 0..support.len() {
        for i2 in  0..i1 {
            let qbit1 = support[i1];
            let qbit2 = support[i2];
            let init_id_count_1 = bucket.count_id(qbit1) ;
            let init_id_count_2 = bucket.count_id(qbit2) ;
            for chunk in ALL_CHUNKS {
                // conjugating with the chunk
                conjugate_with_chunk(bucket, &chunk, qbit1, qbit2, false);
                let new_count_1 = bucket.count_id(qbit1);
                let new_count_2 = bucket.count_id(qbit2);
                let score_1 = new_count_1 as i32 - init_id_count_1 as i32;
                let score_2 = new_count_2 as i32 - init_id_count_2 as i32;
                assert!(score_1 == 0 || score_2 == 0);
                let score = if score_1 == 0 { score_2 } else { score_1 };
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

fn build_graph(bucket: &mut PauliSet) -> (UnGraph<(), i32>, HashMap<(usize, usize), Chunk>) {
    let mut graph: UnGraph<(), i32> = UnGraph::new_undirected();
    let mut best_chunks: HashMap<(usize, usize), Chunk> = HashMap::new();
    for _ in 0..bucket.n {
        graph.add_node(());
    }
    for qbit1 in 0..bucket.n {
        for qbit2 in (qbit1 + 1)..bucket.n {
            // computing the initial identity count

            let init_id_count = bucket.count_id(qbit1) + bucket.count_id(qbit2);
            let mut max_score = 0;
            let mut best_chunk: Chunk = [None; 3];
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

fn single_synthesis_step_depth(bucket: &mut PauliSet) -> CliffordCircuit {
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

pub fn single_synthesis_step(bucket: &mut PauliSet, metric: &Metric) -> CliffordCircuit {
    return match metric {
        Metric::COUNT => single_synthesis_step_count(bucket),
        Metric::DEPTH => single_synthesis_step_depth(bucket),
    };
}

pub fn pauli_network_synthesis(
    bucket: &mut PauliSet,
    metric: &Metric,
    skip_sort: bool,
) -> CliffordCircuit {
    if bucket.len() == 0 {
        return CliffordCircuit::new(0);
    }

    let nqbits = bucket.n;
    let mut output = CliffordCircuit::new(nqbits);

    loop {
        if !skip_sort {
            bucket.support_size_sort();
        }
        while bucket.support_size(0) <= 1 && bucket.len() > 0 {
            bucket.pop();
        }
        if bucket.len() == 0 {
            break;
        }
        let circuit_piece = single_synthesis_step(bucket, &metric);
        output.extend_with(&circuit_piece);
    }
    return output;
}

#[cfg(test)]
mod greedy_synthesis_tests {
    use super::*;
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
                CliffordGate::CZ(i, j) => {
                    bucket.cz(*i, *j);
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
        let axes = vec!["XX".to_owned(), "ZZ".to_owned(), "YY".to_owned()];
        let mut input = PauliSet::from_slice(&axes);
        let circuit = pauli_network_synthesis(&mut input, &Metric::COUNT, false);
        println!("{circuit:?}");
        assert!(check_circuit(&axes, &circuit))
    }

    #[test]
    fn depth_synthesis() {
        let axes = vec!["XX".to_owned(), "ZZ".to_owned(), "YY".to_owned()];
        let mut input = PauliSet::from_slice(&axes);
        let circuit = pauli_network_synthesis(&mut input, &Metric::DEPTH, false);
        println!("{circuit:?}");
    }

    #[test]
    fn count_synthesis_ss() {
        let axes = vec!["XX".to_owned(), "ZZ".to_owned(), "YY".to_owned()];
        let mut input = PauliSet::from_slice(&axes);
        let circuit = pauli_network_synthesis(&mut input, &Metric::COUNT, true);
        println!("{circuit:?}");
        assert!(check_circuit(&axes, &circuit))
    }

    #[test]
    fn depth_synthesis_ss() {
        let axes = vec!["XX".to_owned(), "ZZ".to_owned(), "YY".to_owned()];
        let mut input = PauliSet::from_slice(&axes);
        let circuit = pauli_network_synthesis(&mut input, &Metric::DEPTH, true);
        println!("{circuit:?}");
    }
}
