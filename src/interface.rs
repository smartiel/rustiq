use super::structures::{CliffordGate, GraphState, Metric, PauliSet};
use super::synthesis::clifford::codiagonalization::{
    codiagonalize as codiag, codiagonalize_subsetwise,
};
use super::synthesis::clifford::graph_state::{
    synthesize_graph_state, synthesize_stabilizer_state,
};
use super::synthesis::clifford::isometry::isometry_synthesis as iso_synth;
use super::synthesis::pauli_network::{check_circuit, greedy_pauli_network};
use crate::structures::{IsometryTableau, PauliLike};
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pyfunction]
/// Single interface function for the 4 greedy algorithms for Pauli network synthesis
/// The return value can be easily used on the python side
pub fn pauli_network_synthesis(
    operator_sequence: Vec<String>,
    metric: Metric,
    preserve_order: bool,
    check: bool,
    nshuffles: usize,
    skip_sort: bool,
    fix_clifford: bool,
) -> Vec<(String, Vec<usize>)> {
    let mut bucket = PauliSet::from_slice(&operator_sequence);
    let circuit = greedy_pauli_network(
        &mut bucket,
        &metric,
        preserve_order,
        nshuffles,
        skip_sort,
        fix_clifford,
    );
    if check {
        check_circuit(&operator_sequence, &circuit);
    }
    return circuit.gates.iter().map(|gate| gate.to_vec()).collect();
}

#[pyfunction]
/// Single interface function for the 2 algorithms for stabilizer state synthesis
/// The return value can be easily used on the python side
pub fn stabilizer_state_synthesis(
    stabilizers: Vec<String>,
    metric: Metric,
    niter: usize,
) -> Vec<(String, Vec<usize>)> {
    let ps = PauliSet::from_slice(&stabilizers);
    let circuit = synthesize_stabilizer_state(&ps, &metric, niter);
    circuit.gates.iter().map(|gate| gate.to_vec()).collect()
}

#[pyfunction]
/// Single interface function for the 2 algorithms for graph state synthesis
/// The return value can be easily used on the python side
pub fn graph_state_synthesis(
    graph_adj: Vec<Vec<bool>>,
    metric: Metric,
    niter: usize,
) -> Vec<(String, Vec<usize>)> {
    let gs = GraphState::from_adj(graph_adj);
    let mut circuit = synthesize_graph_state(&gs, &metric, niter);
    for i in 0..gs.n {
        circuit.gates.push(CliffordGate::H(i));
    }
    circuit.gates.iter().map(|gate| gate.to_vec()).collect()
}

#[pyfunction]
/// Single interface function for the 2 algorithms for codiagonalization
/// The return value can be easily used on the python side
pub fn codiagonalization(
    paulis: Vec<String>,
    metric: Metric,
    niter: usize,
) -> Vec<(String, Vec<usize>)> {
    let mut pset = PauliSet::from_slice(&paulis);
    let circuit = codiag(&mut pset, &metric, niter);
    return circuit.gates.iter().map(|gate| gate.to_vec()).collect();
}

#[pyfunction]
/// Single interface function for the subset-wise codiagonalization algorithm
/// The return value can be easily used on the python side
pub fn codiagonalization_sswise(paulis: Vec<String>, k: usize) -> Vec<(String, Vec<usize>)> {
    let pset = PauliSet::from_slice(&paulis);
    let circuit = codiagonalize_subsetwise(&pset, k);
    return circuit.gates.iter().map(|gate| gate.to_vec()).collect();
}

#[pyfunction]
/// Single interface function for the isometry synthesis algorithm
/// The return value can be easily used on the python side
pub fn isometry_synthesis(
    logicals: Vec<String>,
    stabilizers: Vec<String>,
    metric: Metric,
    niter: usize,
) -> Vec<(String, Vec<usize>)> {
    let logicals = PauliSet::from_slice(&logicals);
    let stabilizers = PauliSet::from_slice(&stabilizers);
    let isometry = IsometryTableau {
        n: logicals.len() / 2,
        k: stabilizers.len(),
        logicals,
        stabilizers,
    };
    let circuit = iso_synth(&isometry, &metric, niter);
    return circuit.gates.iter().map(|gate| gate.to_vec()).collect();
}

#[pymodule]
fn rustiq(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(pauli_network_synthesis))?;
    m.add_wrapped(wrap_pyfunction!(graph_state_synthesis))?;
    m.add_wrapped(wrap_pyfunction!(stabilizer_state_synthesis))?;
    m.add_wrapped(wrap_pyfunction!(codiagonalization))?;
    m.add_wrapped(wrap_pyfunction!(codiagonalization_sswise))?;
    m.add_wrapped(wrap_pyfunction!(isometry_synthesis))?;
    m.add_class::<Metric>()?;
    Ok(())
}
