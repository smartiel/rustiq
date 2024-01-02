use crate::routines::f2_linalg::*;
use crate::structures::{CliffordCircuit, CliffordGate, IsometryTableau};

pub fn extract_abcd(isometry: &IsometryTableau) -> (Matrix, Matrix, Matrix, Matrix) {
    let mut a = vec![vec![false; isometry.n + isometry.k]; isometry.n + isometry.k];
    let mut b = vec![vec![false; isometry.n + isometry.k]; isometry.n + isometry.k];
    for i in 0..isometry.n {
        let (_, vec) = isometry.logicals.get_as_vec_bool(i + isometry.n);
        for j in 0..(isometry.n + isometry.k) {
            a[j][i] = vec[j + isometry.n + isometry.k];
            b[j][i] = vec[j];
        }
    }
    for i in 0..isometry.k {
        let (_, vec) = isometry.stabilizers.get_as_vec_bool(i);
        for j in 0..(isometry.n + isometry.k) {
            a[j][i + isometry.n] = vec[j + isometry.n + isometry.k];
            b[j][i + isometry.n] = vec[j];
        }
    }
    let mut c = vec![vec![false; isometry.n]; isometry.n + isometry.k];
    let mut d = vec![vec![false; isometry.n]; isometry.n + isometry.k];
    for i in 0..isometry.n {
        let (_, vec) = isometry.logicals.get_as_vec_bool(i);
        for j in 0..(isometry.n + isometry.k) {
            c[j][i] = vec[j + isometry.n + isometry.k];
            d[j][i] = vec[j];
        }
    }
    return (a, b, c, d);
}

fn make_b_full_rank(
    a: &mut Matrix,
    b: &mut Matrix,
    c: &mut Matrix,
    d: &mut Matrix,
) -> CliffordCircuit {
    let mut circuit = CliffordCircuit::new(a.len());
    let mut wit = Vec::new();
    for i in 0..b.len() {
        wit.push(b[i].clone());
        if f2_rank(&wit) < i + 1 {
            wit[i] = a[i].clone();
            std::mem::swap(&mut a[i], &mut b[i]);
            std::mem::swap(&mut c[i], &mut d[i]);
            circuit.gates.push(CliffordGate::H(i));
        }
    }
    assert_eq!(f2_rank(&b), b.len());
    circuit
}

pub fn decompose(isometry: &IsometryTableau) -> (Matrix, Matrix, Matrix, CliffordCircuit) {
    let (mut a, mut b, mut c, mut d) = extract_abcd(isometry);
    let piece = make_b_full_rank(&mut a, &mut b, &mut c, &mut d);
    let inv_b = inverse_f2(&b);
    let b_k: Matrix = inv_b.clone().drain(..isometry.n).collect();
    let gn = mult_f2(&a, &inv_b);
    let gk = mult_f2(&inv_b, &d).drain(..isometry.n).collect();
    return (gk, gn, b_k, piece);
}
// Write a test module for this function
#[cfg(test)]
mod tests {
    use super::*;
    use crate::structures::IsometryTableau;

    fn is_symmetric(table: &Matrix) -> bool {
        for i in 0..table.len() {
            for j in 0..table.len() {
                if table[j][i] != table[i][j] {
                    return false;
                }
            }
        }
        return true;
    }
    #[test]
    fn test_decompose() {
        let n = 2; // Number of qubits
        let k = 1; // Size of the extension
        let isometry = IsometryTableau::random(n, k);

        let (gk, gn, b_k, _) = decompose(&isometry);

        assert_eq!(gk.len(), n);
        assert!(gk.iter().all(|row| row.len() == n));
        assert_eq!(gn.len(), k + n);
        assert!(gn.iter().all(|row| row.len() == k + n));

        assert!(is_symmetric(&gk));
        assert!(is_symmetric(&gn));
        assert_eq!(f2_rank(&b_k), b_k.len());
    }
}
