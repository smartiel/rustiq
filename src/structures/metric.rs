use super::CliffordCircuit;
use pyo3::prelude::*;
#[pyclass]
#[derive(Debug, Clone)]
pub enum Metric {
    COUNT,
    DEPTH,
}

impl Metric {
    /// Attempts to build a Metric from a string
    pub fn from_string(name: &str) -> Result<Self, String> {
        match name {
            "depth" => Result::Ok(Self::DEPTH),
            "count" => Result::Ok(Self::COUNT),
            &_ => Result::Err(format!("Unknown metric name `{}`", name)),
        }
    }

    pub fn on_circuit(&self, circuit: &CliffordCircuit) -> usize {
        match self {
            Self::DEPTH => circuit.entangling_depth(),
            Self::COUNT => circuit.entangling_count(),
        }
    }
}
