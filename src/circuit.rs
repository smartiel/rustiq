#[derive(Debug, Clone, Copy)]
pub enum Gate {
    CNOT(usize, usize),
    H(usize),
    S(usize),
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
}
