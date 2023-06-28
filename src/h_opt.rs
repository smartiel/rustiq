use super::circuit::{Circuit, Gate};
use super::greedy_order_preserving::PauliDag;
use super::pauli_set::PauliSet;
use pyo3::prelude::*;

#[pyfunction]
pub fn diagonalization_network(input: Vec<String>) -> Vec<(Vec<(String, Vec<usize>)>, String)> {
    let mut pauli_set = PauliSet::from_slice(&input);
    let mut output = Vec::new();
    let mut rotations = Vec::new();
    if pauli_set.len() == 0 {
        return Vec::new();
    }
    for i in 0..pauli_set.len() {
        let (_, vec) = pauli_set.get_as_vec_bool(i);
        let mut x_support: Vec<usize> = (0..pauli_set.n).filter(|qbit| vec[*qbit]).collect();
        let mut piece = Circuit::new(pauli_set.n);
        if let Some(control) = x_support.pop() {
            for target in x_support {
                pauli_set.cnot(control, target);
            }
            let (_, vec) = pauli_set.get_as_vec_bool(i);
            if vec[pauli_set.n + control] {
                piece.gates.push(Gate::S(control));
                pauli_set.s(control);
            }
            piece.gates.push(Gate::H(control));
            pauli_set.h(control);
        }
        output.push(piece);
        let (_, axis) = pauli_set.get(i);
        rotations.push(axis);
    }
    return output.iter().map(|c| c.to_vec()).zip(rotations).collect();
}

pub fn h_opt(axes: Vec<String>) {
    let mut dag = PauliDag::from_slice(&axes);
}
