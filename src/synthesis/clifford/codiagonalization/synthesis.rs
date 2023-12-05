use crate::structures::{CliffordCircuit, Metric, PauliSet};

use super::count::codiagonalize_count;
use super::depth::codiagonalize_depth;

pub fn codiagonalize(pauli_set: &mut PauliSet, metric: &Metric, niter: usize) -> CliffordCircuit {
    match metric {
        Metric::COUNT => codiagonalize_count(pauli_set, niter),
        Metric::DEPTH => codiagonalize_depth(pauli_set),
    }
}
