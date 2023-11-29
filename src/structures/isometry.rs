use super::pauli_like::PauliLike;
use super::pauli_set::PauliSet;

pub struct IsometryTableau {
    n: usize,
    k: usize,
    logicals: PauliSet,
    stabilizers: PauliSet,
}

impl IsometryTableau {
    pub fn new(n: usize, k: usize) -> Self {
        let mut logicals = PauliSet::new(n + k);
        for i in 0..n {
            let mut vecbool = vec![false; 2 * n + 2 * k];
            vecbool[i] = true;
            logicals.insert_vec_bool(&vecbool, false);
        }
        for i in 0..n {
            let mut vecbool = vec![false; 2 * n + 2 * k];
            vecbool[i + n + k] = true;
            logicals.insert_vec_bool(&vecbool, false);
        }
        let mut stabilizers = PauliSet::new(n + k);
        for i in 0..k {
            let mut vecbool = vec![false; 2 * n + 2 * k];
            vecbool[n + i] = true;
            stabilizers.insert_vec_bool(&vecbool, false);
        }
        IsometryTableau {
            n,
            k,
            logicals,
            stabilizers,
        }
    }
}

impl PauliLike for IsometryTableau {
    fn h(&mut self, i: usize) {
        self.logicals.h(i);
        self.stabilizers.h(i);
    }

    fn s(&mut self, i: usize) {
        self.logicals.s(i);
        self.stabilizers.s(i);
    }

    fn sd(&mut self, i: usize) {
        self.logicals.sd(i);
        self.stabilizers.sd(i);
    }

    fn sqrt_x(&mut self, i: usize) {
        self.logicals.sqrt_x(i);
        self.stabilizers.sqrt_x(i);
    }

    fn sqrt_xd(&mut self, i: usize) {
        self.logicals.sqrt_xd(i);
        self.stabilizers.sqrt_xd(i);
    }

    fn cnot(&mut self, i: usize, j: usize) {
        self.logicals.cnot(i, j);
        self.stabilizers.cnot(i, j);
    }
}
