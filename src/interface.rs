use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use rustiq_core::routines::rotation_extraction::extract_rotations as extract_rot;
use rustiq_core::routines::rotation_optimization::full_initial_state_propagation;
use rustiq_core::routines::rotation_optimization::zhang_rotation_optimization as zhang_opt;
use rustiq_core::structures::{CliffordGate, GraphState, Metric as RsMetric, PauliSet};
use rustiq_core::structures::{IsometryTableau, Parameter, PauliLike, Tableau};
use rustiq_core::synthesis::clifford::codiagonalization::{
    codiagonalize as codiag, codiagonalize_subsetwise,
};
use rustiq_core::synthesis::clifford::graph_state::{
    synthesize_graph_state, synthesize_stabilizer_state,
};
use rustiq_core::synthesis::clifford::isometry::isometry_synthesis as iso_synth;
use rustiq_core::synthesis::pauli_network::{check_circuit, greedy_pauli_network};

#[pyclass]
#[derive(Clone)]
pub struct Metric {
    data: RsMetric,
}

impl Metric {
    pub fn new(metric_type: &str) -> Self {
        Self {
            data: match metric_type {
                "count" => RsMetric::COUNT,
                "depth" => RsMetric::DEPTH,
                _ => panic!(
                    "Unknown metric type `{}` (must be 'count' or 'depth')",
                    metric_type
                ),
            },
        }
    }

    pub fn unwrap(self) -> RsMetric {
        return self.data;
    }
}

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
        &metric.unwrap(),
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
    let circuit = synthesize_stabilizer_state(&ps, &metric.unwrap(), niter);
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
    let mut circuit = synthesize_graph_state(&gs, &metric.unwrap(), niter);
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
    let circuit = codiag(&mut pset, &metric.unwrap(), niter);
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
    let circuit = iso_synth(&isometry, &metric.unwrap(), niter);
    return circuit.gates.iter().map(|gate| gate.to_vec()).collect();
}

#[pyfunction]
pub fn extract_rotations(
    circuit: Vec<(String, Vec<usize>)>,
    angles: Vec<String>,
    nqubits: usize,
    optimize: bool,
) -> (Vec<(String, String)>, Vec<(bool, String)>) {
    let (rotations, mut clifford) = extract_rot(&circuit, nqubits);
    let angles = angles.into_iter().map(|s| Parameter::from_string(s));
    let mut rotations: Vec<_> = angles
        .zip(rotations)
        .map(|(mut angle, (phase, axis))| {
            if phase {
                angle.flip_sign();
            }
            (axis, angle)
        })
        .collect();
    if optimize {
        let (new_rotations, inverse_final_clifford) = zhang_opt(rotations, nqubits);
        rotations = new_rotations;
        clifford = clifford * inverse_final_clifford.adjoint();
    }
    let rotations = rotations
        .into_iter()
        .map(|(x, y)| (x, y.to_string()))
        .collect();
    return (
        rotations,
        (0..2 * nqubits).map(|i| clifford.logicals.get(i)).collect(),
    );
}

#[pyfunction]
pub fn zhang_rotation_optimization(
    rotations: Vec<(String, String)>,
    nqubits: usize,
) -> (Vec<(String, String)>, Vec<(bool, String)>) {
    let rotations = rotations
        .into_iter()
        .map(|(a, b)| (a, Parameter::from_string(b)))
        .collect();
    let (rotations, inverse_final_clifford) = zhang_opt(rotations, nqubits);
    let rotations = rotations
        .into_iter()
        .map(|(x, y)| (x, y.to_string()))
        .collect();
    let clifford = inverse_final_clifford.adjoint();
    return (
        rotations,
        (0..2 * nqubits).map(|i| clifford.logicals.get(i)).collect(),
    );
}

#[pyfunction]
pub fn initial_state_propagation(
    rotations: Vec<(String, String)>,
) -> (Vec<(String, String)>, Vec<(bool, String)>) {
    let rotations: Vec<_> = rotations
        .into_iter()
        .map(|(a, b)| (a, Parameter::from_string(b)))
        .collect();
    let (rotations, final_clifford) = full_initial_state_propagation(&rotations);
    let rotations = rotations
        .into_iter()
        .map(|(x, y)| (x, y.to_string()))
        .collect();
    return (
        rotations,
        (0..2 * final_clifford.logicals.n)
            .map(|i| final_clifford.logicals.get(i))
            .collect(),
    );
}

#[pyfunction]
pub fn tableau_mul(t1: Vec<(bool, String)>, t2: Vec<(bool, String)>) -> Vec<(bool, String)> {
    let t1 = Tableau::from_operators(&t1);
    let t2 = Tableau::from_operators(&t2);
    let t3 = t1 * t2;
    (0..2 * t3.logicals.n).map(|i| t3.logicals.get(i)).collect()
}

#[pymodule]
fn rustiq(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pauli_network_synthesis, m)?)?;
    m.add_function(wrap_pyfunction!(graph_state_synthesis, m)?)?;
    m.add_function(wrap_pyfunction!(stabilizer_state_synthesis, m)?)?;
    m.add_function(wrap_pyfunction!(codiagonalization, m)?)?;
    m.add_function(wrap_pyfunction!(codiagonalization_sswise, m)?)?;
    m.add_function(wrap_pyfunction!(isometry_synthesis, m)?)?;
    m.add_function(wrap_pyfunction!(extract_rotations, m)?)?;
    m.add_function(wrap_pyfunction!(zhang_rotation_optimization, m)?)?;
    m.add_function(wrap_pyfunction!(initial_state_propagation, m)?)?;
    m.add_function(wrap_pyfunction!(tableau_mul, m)?)?;
    m.add_class::<Metric>()?;
    Ok(())
}
