use pyo3::prelude::*;
#[derive(Debug, Clone, Copy)]
pub enum Gate {
    CNOT(usize, usize),
    H(usize),
    S(usize),
    Sd(usize),
    SqrtX(usize),
    SqrtXd(usize),
}
impl Gate {
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
            Gate::CNOT(i, j) => ("CNOT".to_owned(), vec![*i, *j]),
            Gate::H(i) => ("H".to_owned(), vec![*i]),
            Gate::S(i) => ("S".to_owned(), vec![*i]),
            Gate::Sd(i) => ("Sd".to_owned(), vec![*i]),
            Gate::SqrtX(i) => ("SqrtX".to_owned(), vec![*i]),
            Gate::SqrtXd(i) => ("SqrtXd".to_owned(), vec![*i]),
        }
    }
}
#[pyclass]
#[derive(Debug)]
pub struct Circuit {
    pub nqbits: usize,
    pub gates: Vec<Gate>,
}

impl Circuit {
    pub fn new(nqbits: usize) -> Self {
        Self {
            nqbits,
            gates: Vec::new(),
        }
    }
    pub fn extend_with(&mut self, other: &Circuit) {
        self.gates.extend_from_slice(&other.gates);
    }
    /// Counts the number of CNOT gates
    pub fn cnot_count(&self) -> usize {
        return self
            .gates
            .iter()
            .filter(|gate| {
                if let Gate::CNOT(_, _) = gate {
                    true
                } else {
                    false
                }
            })
            .count();
    }
    /// Computes the CNOT depth of the circuit
    pub fn cnot_depth(&self) -> usize {
        let mut depths: Vec<usize> = vec![0; self.nqbits];
        for gate in self.gates.iter() {
            match gate {
                Gate::CNOT(i, j) => {
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
    pub fn to_vec(&self) -> Vec<(String, Vec<usize>)> {
        self.gates.iter().map(|gate| gate.to_vec()).collect()
    }
}
