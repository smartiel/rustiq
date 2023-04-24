#[derive(Debug, Clone, Copy)]
pub enum Gate {
    CNOT(usize, usize),
    H(usize),
    S(usize),
    SqrtX(usize),
}

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
}
