use pyo3::prelude::*;
use rand::Rng;
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CliffordGate {
    CNOT(usize, usize),
    CZ(usize, usize),
    H(usize),
    S(usize),
    Sd(usize),
    SqrtX(usize),
    SqrtXd(usize),
}
impl CliffordGate {
    pub fn dagger(&self) -> Self {
        match self {
            Self::S(i) => Self::Sd(*i),
            Self::SqrtX(i) => Self::SqrtXd(*i),
            Self::Sd(i) => Self::S(*i),
            Self::SqrtXd(i) => Self::SqrtX(*i),
            _ => *self,
        }
    }
    pub fn to_vec(&self) -> (String, Vec<usize>) {
        match self {
            CliffordGate::CNOT(i, j) => ("CNOT".to_owned(), vec![*i, *j]),
            CliffordGate::CZ(i, j) => ("CZ".to_owned(), vec![*i, *j]),
            CliffordGate::H(i) => ("H".to_owned(), vec![*i]),
            CliffordGate::S(i) => ("S".to_owned(), vec![*i]),
            CliffordGate::Sd(i) => ("Sd".to_owned(), vec![*i]),
            CliffordGate::SqrtX(i) => ("SqrtX".to_owned(), vec![*i]),
            CliffordGate::SqrtXd(i) => ("SqrtXd".to_owned(), vec![*i]),
        }
    }
}
#[pyclass]
#[derive(Debug)]
pub struct CliffordCircuit {
    pub nqbits: usize,
    pub gates: Vec<CliffordGate>,
}

impl CliffordCircuit {
    pub fn new(nqbits: usize) -> Self {
        Self {
            nqbits,
            gates: Vec::new(),
        }
    }

    pub fn random(nqubits: usize, ngates: usize) -> Self {
        let mut rng = rand::thread_rng();
        let mut circuit = Self::new(nqubits);
        for _ in 0..ngates {
            if rng.gen_bool(0.5) {
                // CNOT
                let i = rng.gen_range(0..nqubits);
                let mut j = rng.gen_range(0..nqubits);
                while j == i {
                    j = rng.gen_range(0..nqubits);
                }
                circuit.gates.push(CliffordGate::CNOT(i, j));
                continue;
            }
            if rng.gen_bool(0.5) {
                // H
                let i = rng.gen_range(0..nqubits);
                circuit.gates.push(CliffordGate::H(i));
                continue;
            }
            let i = rng.gen_range(0..nqubits);
            circuit.gates.push(CliffordGate::S(i));
        }
        circuit
    }

    pub fn extend_with(&mut self, other: &CliffordCircuit) {
        self.gates.extend_from_slice(&other.gates);
    }
    /// Counts the number of CNOT gates
    pub fn cnot_count(&self) -> usize {
        return self
            .gates
            .iter()
            .filter(|gate| {
                if let CliffordGate::CNOT(_, _) = gate {
                    true
                } else {
                    false
                }
            })
            .count();
    }
    /// Counts the number of CNOT gates
    pub fn entangling_count(&self) -> usize {
        return self
            .gates
            .iter()
            .filter(|gate| match gate {
                CliffordGate::CNOT(_, _) => true,
                CliffordGate::CZ(_, _) => true,
                _ => false,
            })
            .count();
    }
    /// Computes the CNOT depth of the circuit
    pub fn cnot_depth(&self) -> usize {
        let mut depths: Vec<usize> = vec![0; self.nqbits];
        for gate in self.gates.iter() {
            match gate {
                CliffordGate::CNOT(i, j) => {
                    let gate_depth = std::cmp::max(depths[*i], depths[*j]) + 1;
                    depths[*i] = gate_depth;
                    depths[*j] = gate_depth;
                }
                _ => {}
            }
        }
        return *depths.iter().max().unwrap();
    }
    /// Computes the CNOT depth of the circuit
    pub fn entangling_depth(&self) -> usize {
        let mut depths: Vec<usize> = vec![0; self.nqbits];
        for gate in self.gates.iter() {
            match gate {
                CliffordGate::CNOT(i, j) => {
                    let gate_depth = std::cmp::max(depths[*i], depths[*j]) + 1;
                    depths[*i] = gate_depth;
                    depths[*j] = gate_depth;
                }
                CliffordGate::CZ(i, j) => {
                    let gate_depth = std::cmp::max(depths[*i], depths[*j]) + 1;
                    depths[*i] = gate_depth;
                    depths[*j] = gate_depth;
                }
                _ => {}
            }
        }
        return *depths.iter().max().unwrap();
    }
    /// Returns the inverse of the circuit
    pub fn dagger(&self) -> Self {
        let new_gates = self.gates.iter().rev().map(|gate| gate.dagger()).collect();
        return Self {
            nqbits: self.nqbits,
            gates: new_gates,
        };
    }
}
