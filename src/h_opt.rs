use super::circuit::{Circuit, Gate};
use super::pauli_set::PauliSet;
use pyo3::prelude::*;

pub fn diagonalization_network_pset(pauli_set: &mut PauliSet) -> Vec<(Circuit, String)> {
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
    return output.into_iter().zip(rotations).collect();
}

pub fn h_opt(
    axes: Vec<String>,
) -> (
    Vec<(String, Vec<usize>)>,
    Vec<(Vec<(String, Vec<usize>)>, String)>,
) {
    let reversed_axes: Vec<_> = axes.clone().into_iter().rev().collect();
    let mut reversed_pset = PauliSet::from_slice(&reversed_axes);
    let diag_net = diagonalization_network_pset(&mut reversed_pset);
    let mut diag_net_generators = PauliSet::new(reversed_pset.n);

    for i in 0..reversed_pset.n {
        let z_pstring = (0..reversed_pset.n)
            .map(|j| if j == i { "Z" } else { "I" })
            .collect::<Vec<_>>()
            .join("");
        diag_net_generators.insert(&z_pstring, false);
    }

    for (circuit, _) in diag_net.iter().rev() {
        diag_net_generators.conjugate_with_circuit(&circuit);
    }
    for pauli in axes {
        diag_net_generators.insert(&pauli, false);
    }
    let mut final_diag = diagonalization_network_pset(&mut diag_net_generators);

    let mut first_circuit = Circuit::new(reversed_pset.n);
    for (circuit, _) in final_diag.drain(..reversed_pset.n) {
        first_circuit.gates.extend(circuit.gates.into_iter());
    }
    if let Some((circuit, _)) = final_diag.get_mut(0) {
        first_circuit.gates.extend(circuit.gates.drain(..));
    }

    return (
        first_circuit.to_vec(),
        final_diag
            .into_iter()
            .map(|(c, s)| (c.to_vec(), s))
            .collect(),
    );
}

#[pyfunction]
pub fn diagonalization_network(
    input: Vec<String>,
    optimal: bool,
) -> (
    Vec<(String, Vec<usize>)>,
    Vec<(Vec<(String, Vec<usize>)>, String)>,
) {
    if optimal {
        return h_opt(input);
    }
    let mut pauli_set = PauliSet::from_slice(&input);
    let output = diagonalization_network_pset(&mut pauli_set);
    return (
        Vec::new(),
        output.into_iter().map(|(c, s)| (c.to_vec(), s)).collect(),
    );
}
