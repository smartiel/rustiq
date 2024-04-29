use super::pauli_like::PauliLike;
use super::pauli_set::PauliSet;
use super::{CliffordCircuit, IsometryTableau};
use rand::Rng;

fn compute_phase_product_pauli(pset0: &PauliSet, vec: &Vec<bool>) -> bool {
    let mut phase = false;
    for j in 0..2 * pset0.n {
        phase ^= pset0.get_phase(j) & vec[j];
    }
    let mut ifact: u8 = 0;
    for i in 0..pset0.n {
        if vec[i] & vec[i + pset0.n] {
            ifact += 1;
        }
    }
    ifact = ifact % 4;
    for j in 0..pset0.n {
        let mut x: bool = false;
        let mut z: bool = false;
        for i in 0..2 * pset0.n {
            if vec[i] {
                let x1: bool = pset0.get_entry(j, i);
                let z1: bool = pset0.get_entry(j + pset0.n, i);
                let entry = (x1, z1, x, z);
                if LOOKUP_0.contains(&entry) {
                    ifact += 1;
                }
                if LOOKUP_1.contains(&entry) {
                    ifact += 3;
                }
                x ^= x1;
                z ^= z1;
                ifact = ifact % 4;
            }
        }
    }
    return (((ifact % 4) >> 1) != 0) ^ phase;
}

#[derive(Debug, Clone, PartialEq)]
pub struct Tableau {
    pub logicals: PauliSet,
}

impl Tableau {
    /// Allocates a new Tableau representing the identity operator over `n` qubits
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
    /// Build the Tableau corresponding to a Clifford circuit
    pub fn from_circuit(circuit: &CliffordCircuit) -> Self {
        let mut tab = Self::new(circuit.nqbits);
        tab.conjugate_with_circuit(circuit);
        return tab;
    }
    /// Generates a random Tableau (no garantuees, just here for testing)
    pub fn random(n: usize) -> Self {
        let mut rng = rand::thread_rng();
        let mut iso = Self::new(n);
        for _ in 0..(n) * (n) {
            let i = rng.gen::<usize>() % (n);
            loop {
                let j = rng.gen::<usize>() % (n);
                if i == j {
                    continue;
                }
                iso.cnot(i, j);
                break;
            }
            for _ in 0..(n) {
                let gate_i = rng.gen::<u8>() % 3;
                if gate_i == 1 {
                    let q = rng.gen::<usize>() % (n);
                    iso.h(q);
                }
                if gate_i == 2 {
                    let q = rng.gen::<usize>() % (n);
                    iso.s(q);
                }
            }
        }
        return iso;
    }
    /// Build a Tableau from a PauliSet
    pub fn from_operators(logicals: &Vec<(bool, String)>) -> Self {
        if logicals.len() == 0 {
            return Self::new(0);
        }
        let nqbits = logicals[0].1.len();
        let mut pset = PauliSet::new(nqbits);
        for (phase, string) in logicals {
            pset.insert(string, *phase);
        }
        Self { logicals: pset }
    }
    /// Returns the inverse Tableau
    pub fn adjoint(&self) -> Self {
        let mut new_logicals = PauliSet::new(self.logicals.n);
        for i in 0..self.logicals.n {
            let (_, string) = self.logicals.get_inverse_x(i);
            new_logicals.insert(&string, self.logicals.get_phase(i));
        }
        for i in 0..self.logicals.n {
            let (_, string) = self.logicals.get_inverse_z(i);
            new_logicals.insert(&string, self.logicals.get_phase(i + self.logicals.n));
        }
        let prod = self.clone()
            * Tableau {
                logicals: new_logicals.clone(),
            };
        for i in 0..2 * self.logicals.n {
            new_logicals.set_phase(i, new_logicals.get_phase(i) ^ prod.logicals.get_phase(i));
        }
        return Self {
            logicals: new_logicals,
        };
    }

    pub fn get_inverse_z(&self, qbit: usize) -> (bool, String) {
        let (_, string) = self.logicals.get_inverse_z(qbit);
        let mut as_vec_bool = vec![false; 2 * self.logicals.n];
        for qbit in 0..self.logicals.n {
            match string.chars().nth(qbit).unwrap() {
                'X' => {
                    as_vec_bool[qbit] = true;
                }
                'Y' => {
                    as_vec_bool[qbit] = true;
                    as_vec_bool[qbit + self.logicals.n] = true;
                }
                'Z' => {
                    as_vec_bool[qbit + self.logicals.n] = true;
                }
                _ => {}
            }
        }
        let phase = compute_phase_product_pauli(&self.logicals, &as_vec_bool);
        return (phase, string);
    }
    pub fn get_inverse_x(&self, qbit: usize) -> (bool, String) {
        let (_, string) = self.logicals.get_inverse_x(qbit);
        let mut as_vec_bool = vec![false; 2 * self.logicals.n];
        for qbit in 0..self.logicals.n {
            match string.chars().nth(qbit).unwrap() {
                'X' => {
                    as_vec_bool[qbit] = true;
                }
                'Y' => {
                    as_vec_bool[qbit] = true;
                    as_vec_bool[qbit + self.logicals.n] = true;
                }
                'Z' => {
                    as_vec_bool[qbit + self.logicals.n] = true;
                }
                _ => {}
            }
        }
        let phase = compute_phase_product_pauli(&self.logicals, &as_vec_bool);
        return (phase, string);
    }
    /// Lifts the Taleau into an IsometryTableau (k = 0)
    pub fn to_isometry(self) -> IsometryTableau {
        IsometryTableau {
            n: self.logicals.n,
            k: 0,
            stabilizers: PauliSet::new(self.logicals.n),
            logicals: self.logicals,
        }
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

const LOOKUP_0: [(bool, bool, bool, bool); 3] = [
    (false, true, true, true),
    (true, false, false, true),
    (true, true, true, false),
];

const LOOKUP_1: [(bool, bool, bool, bool); 3] = [
    (false, true, true, false),
    (true, false, true, true),
    (true, true, false, true),
];

fn compute_phase_product_single_col(pset0: &PauliSet, pset1: &PauliSet, col: usize) -> u8 {
    let mut ifact = pset1.get_i_factors_single_col(col);
    for j in 0..pset0.n {
        let mut x: bool = false;
        let mut z: bool = false;
        for i in 0..2 * pset0.n {
            if pset1.get_entry(i, col) {
                let x1: bool = pset0.get_entry(j, i);
                let z1: bool = pset0.get_entry(j + pset0.n, i);
                let entry = (x1, z1, x, z);
                if LOOKUP_0.contains(&entry) {
                    ifact += 1;
                }
                if LOOKUP_1.contains(&entry) {
                    ifact += 3;
                }
                x ^= x1;
                z ^= z1;
                ifact = ifact % 4;
            }
        }
    }
    return ifact % 4;
}
fn compute_phase_product(pset0: &PauliSet, pset1: &PauliSet) {
    let mut ifacts = pset1.get_i_factors();
    for k in 0..2 * pset0.n {
        for j in 0..pset0.n {
            let mut x: bool = false;
            let mut z: bool = false;
            for i in 0..2 * pset0.n {
                if pset1.get_entry(i, k) {
                    let x1: bool = pset0.get_entry(j, i);
                    let z1: bool = pset0.get_entry(j + pset0.n, i);
                    let entry = (x1, z1, x, z);
                    if LOOKUP_0.contains(&entry) {
                        ifacts[k] += 1;
                    }
                    if LOOKUP_1.contains(&entry) {
                        ifacts[k] += 3;
                    }
                    x ^= x1;
                    z ^= z1;
                    ifacts[k] = ifacts[k] % 4;
                }
            }
        }
        ifacts[k] = ifacts[k] % 4;
    }
}

impl std::ops::Mul<Tableau> for Tableau {
    type Output = Tableau;
    fn mul(self, rhs: Tableau) -> Self::Output {
        assert_eq!(self.logicals.n, rhs.logicals.n);
        let mut new_tableau = Tableau::new(self.logicals.n);
        for i in 0..2 * self.logicals.n {
            let (mut phase, col) = rhs.logicals.get_as_vec_bool(i);
            for j in 0..2 * self.logicals.n {
                phase ^= self.logicals.get_phase(j) & col[j];
            }
            new_tableau.logicals.set_phase(i, phase);
        }
        let mut ifacts = rhs.logicals.get_i_factors();
        for k in 0..2 * self.logicals.n {
            for j in 0..self.logicals.n {
                let mut x: bool = false;
                let mut z: bool = false;
                for i in 0..2 * self.logicals.n {
                    if rhs.logicals.get_entry(i, k) {
                        let x1: bool = self.logicals.get_entry(j, i);
                        let z1: bool = self.logicals.get_entry(j + self.logicals.n, i);
                        let entry = (x1, z1, x, z);
                        if LOOKUP_0.contains(&entry) {
                            ifacts[k] += 1;
                        }
                        if LOOKUP_1.contains(&entry) {
                            ifacts[k] += 3;
                        }
                        x ^= x1;
                        z ^= z1;
                        ifacts[k] = ifacts[k] % 4;
                    }
                }
            }
            ifacts[k] = ifacts[k] % 4;
        }
        let p: Vec<bool> = ifacts.into_iter().map(|v| 0 != ((v % 4) >> 1)).collect();
        for (i, ph) in p.iter().enumerate() {
            new_tableau
                .logicals
                .set_phase(i, new_tableau.logicals.get_phase(i) ^ ph);
        }

        for i in 0..2 * self.logicals.n {
            for j in 0..2 * self.logicals.n {
                let (_, col) = rhs.logicals.get_as_vec_bool(j);
                new_tableau
                    .logicals
                    .set_raw_entry(i, j, self.logicals.and_row_acc(i, &col));
            }
        }
        new_tableau
    }
}

#[cfg(test)]
mod tests {
    use super::Tableau;

    #[test]
    fn test_mul_adjoint() {
        let t1 = Tableau::random(5);
        let t2 = t1.adjoint();
        let t3 = t1 * t2;
        let t4 = Tableau::new(5);
        assert_eq!(t3, t4);
    }
}
