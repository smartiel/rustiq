use crate::routines::f2_linalg::*;
use crate::structures::{CliffordCircuit, CliffordGate, GraphState, PauliSet};

fn extract_zx_parts(pset: &PauliSet) -> (Matrix, Matrix) {
    let mut z_part = vec![vec![false; pset.n]; pset.n];
    let mut x_part = vec![vec![false; pset.n]; pset.n];
    for i in 0..pset.n {
        let (_, vec) = pset.get_as_vec_bool(i);
        for j in 0..pset.n {
            z_part[j][i] = vec[j + pset.n];
            x_part[j][i] = vec[j];
        }
    }
    (z_part, x_part)
}

/// Finds the H to inject in order to make the X part of the stabilizers full rank
fn make_x_full_rank(z_part: &mut Matrix, x_part: &mut Matrix) -> CliffordCircuit {
    let mut circuit = CliffordCircuit::new(x_part.len());
    let mut wit = Vec::new();
    for i in 0..x_part.len() {
        wit.push(x_part[i].clone());
        if f2_rank(&wit) < i + 1 {
            wit[i] = z_part[i].clone();
            std::mem::swap(&mut z_part[i], &mut x_part[i]);
            circuit.gates.push(CliffordGate::H(i));
        }
    }
    assert_eq!(f2_rank(&x_part), x_part.len());
    circuit
}

pub fn extract_graph_state(pset: &PauliSet) -> (GraphState, CliffordCircuit) {
    let (mut z_part, mut x_part) = extract_zx_parts(pset);
    let circuit = make_x_full_rank(&mut z_part, &mut x_part);
    let n = z_part.len();
    diagonalize(&mut x_part, &mut z_part, n);
    (GraphState::from_adj(z_part), circuit)
}
