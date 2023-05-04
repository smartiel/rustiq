use super::commutation_dag::pauli_network_synthesis_no_permutation;
use super::greedy_pauli_network::{pauli_network_synthesis, Metric};
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
#[pyfunction]
/// Single interface function for the 4 greedy algorithms
/// The return value can be easily used on the python side
pub fn greedy_pauli_network(
    operator_sequence: Vec<String>,
    metric: Metric,
    preserve_order: bool,
) -> Vec<(String, Vec<usize>)> {
    let circuit = if preserve_order {
        pauli_network_synthesis_no_permutation(operator_sequence, metric)
    } else {
        pauli_network_synthesis(operator_sequence, metric)
    };
    return circuit.gates.iter().map(|gate| gate.to_vec()).collect();
}

#[pymodule]
fn rustiq(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(greedy_pauli_network))?;
    m.add_class::<Metric>()?;
    Ok(())
}
