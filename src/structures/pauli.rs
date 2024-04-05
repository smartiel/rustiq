use super::pauli_like::PauliLike;
use std::ops;

pub struct Pauli {
    pub n: usize,
    pub data: Vec<bool>,
    pub phase: u8,
}

impl Pauli {
    pub fn new(n: usize) -> Self {
        Pauli {
            n,
            data: vec![true; 2 * n],
            phase: 0,
        }
    }
    pub fn from_vec_bool(data: Vec<bool>, phase: u8) -> Self {
        Pauli {
            n: data.len() / 2,
            data,
            phase,
        }
    }
    pub fn commutes(&self, other: &Pauli) -> bool {
        if self.n != other.n {
            panic!("Can't compare two Paulis on different number of qubits");
        }
        let (my_z, my_x) = self.data.split_at(self.n);
        let (their_z, their_x) = other.data.split_at(self.n);
        let p1 = my_z.iter().zip(their_x.iter()).map(|(a, b)| a & b);
        let p2 = my_x.iter().zip(their_z.iter()).map(|(a, b)| a & b);
        return (p1.zip(p2).map(|(a, b)| a ^ b).filter(|a| *a).count() & 1) == 0;
    }
}

impl ops::Mul<Pauli> for Pauli {
    type Output = Pauli;

    fn mul(self, _rhs: Pauli) -> Pauli {
        assert_eq!(self.n, _rhs.n);
        let mut output = Pauli::new(self.n);
        output.phase = self.phase + _rhs.phase;
        for i in 0..self.n {
            if self.data[i] && _rhs.data[i + self.n] {
                output.phase += 2;
            }
            if self.data[i + self.n] && _rhs.data[i] {
                output.phase += 2;
            }
        }
        for i in 0..2 * self.n {
            output.data[i] = self.data[i] ^ _rhs.data[i];
        }
        output.phase = output.phase % 4;
        output
    }
}

impl PauliLike for Pauli {
    fn h(&mut self, i: usize) {
        self.data.swap(i, i + self.n);
    }

    fn s(&mut self, i: usize) {
        self.data[i + self.n] ^= self.data[i];
    }

    fn sd(&mut self, i: usize) {
        self.data[i + self.n] ^= self.data[i];
    }

    fn sqrt_x(&mut self, i: usize) {
        self.data[i] ^= self.data[i + self.n];
    }

    fn sqrt_xd(&mut self, i: usize) {
        self.data[i] ^= self.data[i + self.n];
    }

    fn cnot(&mut self, i: usize, j: usize) {
        self.data[i + self.n] ^= self.data[j + self.n];
        self.data[j] ^= self.data[i];
    }
}
