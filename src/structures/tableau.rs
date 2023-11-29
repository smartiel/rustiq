use super::pauli_like::PauliLike;
use super::pauli_set::PauliSet;

pub struct Tableau {
    logicals: PauliSet,
}

impl Tableau {
    pub fn new(n: usize) -> Self {
        let mut logicals = PauliSet::new(n);
        for i in 0..2 * n {
            // fn insert_vec_bool(&mut self, axis: &Vec<bool>, phase: bool) -> usize {
            let mut vecbool = vec![false; 2 * n];
            vecbool[i] = true;
            logicals.insert_vec_bool(&vecbool, false);
        }
        Tableau { logicals }
    }
}

impl PauliLike for Tableau {
    fn h(&mut self, i: usize) {
        self.logicals.h(i);
    }

    fn s(&mut self, i: usize) {
        self.logicals.s(i);
    }

    fn sd(&mut self, i: usize) {
        self.logicals.sd(i);
    }

    fn sqrt_x(&mut self, i: usize) {
        self.logicals.sqrt_x(i);
    }

    fn sqrt_xd(&mut self, i: usize) {
        self.logicals.sqrt_xd(i);
    }

    fn cnot(&mut self, i: usize, j: usize) {
        self.logicals.cnot(i, j);
    }
}
