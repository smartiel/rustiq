use super::structures::clifford_circuit::{CliffordCircuit, CliffordGate};
use super::structures::metric::Metric;
use super::structures::pauli_set::PauliSet;
use super::synthesis::pauli_network::{
    pauli_network_synthesis, pauli_network_synthesis_no_permutation,
};
use crate::structures::pauli_like::PauliLike;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use std::collections::HashSet;
fn check_circuit(input: &[String], circuit: &CliffordCircuit) {
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
    // println!("Synthesized {} operators", hit_map.len());
    // println!("{:?}", bucket);
    assert!(
        hit_map.len() == input.len(),
        "Synthesized {} operators, expected {}",
        hit_map.len(),
        input.len()
    );
}

#[pyfunction]
/// Single interface function for the 4 greedy algorithms
/// The return value can be easily used on the python side
pub fn greedy_pauli_network(
    operator_sequence: Vec<String>,
    metric: Metric,
    preserve_order: bool,
    check: bool,
) -> Vec<(String, Vec<usize>)> {
    let circuit = if preserve_order {
        pauli_network_synthesis_no_permutation(operator_sequence.clone(), metric)
    } else {
        pauli_network_synthesis(operator_sequence.clone(), metric)
    };
    if check {
        check_circuit(&operator_sequence, &circuit);
    }
    return circuit.gates.iter().map(|gate| gate.to_vec()).collect();
}

#[pymodule]
fn rustiq(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(greedy_pauli_network))?;
    m.add_class::<Metric>()?;
    Ok(())
}
