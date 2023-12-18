/// This module implements codiagonalization algorithms based on https://arxiv.org/abs/2007.10515
use super::common::build_table;
use crate::routines::f2_linalg::{f2_rank, rowop};
use crate::structures::{CliffordCircuit, CliffordGate, PauliSet};
use std::collections::HashSet;

fn vec_xor(v1: &mut Vec<bool>, v2: &Vec<bool>) {
    for (a, b) in v1.iter_mut().zip(v2.iter()) {
        *a ^= *b;
    }
}

fn greedy_reduce(
    z_table: &mut Vec<Vec<bool>>,
    x_table: &mut Vec<Vec<bool>>,
    qubits: &Vec<usize>,
) -> (usize, CliffordCircuit) {
    let mut min_weight = z_table.len() + 1;
    let mut min_index = None;
    for i in 0..x_table.first().unwrap().len() {
        let mut count = 0;
        for j in qubits.iter() {
            if x_table[*j][i] || z_table[*j][i] {
                count += 1;
            }
        }
        if count > 0 && count < min_weight {
            min_weight = count;
            min_index = Some(i);
        }
    }
    if let Some(min_index) = min_index {
        let mut circuit = CliffordCircuit::new(z_table.len());
        let mut support = Vec::new();
        for j in qubits.iter() {
            if x_table[*j][min_index] || z_table[*j][min_index] {
                support.push(*j);
            }
            if x_table[*j][min_index] && z_table[*j][min_index] {
                circuit.gates.push(CliffordGate::SqrtX(*j));
                vec_xor(&mut x_table[*j], &z_table[*j]);
            } else if x_table[*j][min_index] && !z_table[*j][min_index] {
                circuit.gates.push(CliffordGate::H(*j));
                std::mem::swap(&mut x_table[*j], &mut z_table[*j]);
            }
        }
        for control in support.iter().skip(1) {
            circuit.gates.push(CliffordGate::CNOT(*control, support[0]));
            rowop(z_table, support[0], *control);
            rowop(x_table, *control, support[0]);
        }

        return (support[0], circuit);
    }
    panic!("Something went wrong, no qbit was diagonalized");
}

fn find_rank_nm1(
    z_part: &mut Vec<Vec<bool>>,
    x_part: &mut Vec<Vec<bool>>,
) -> Option<CliffordCircuit> {
    // We want to turn x_part into a rank k-1 matrix (or at least attempt to)
    // We are allowed to:
    // * swap row i of Z with row i of X (H gate)
    // * xor row i of Z in row i of X (SQRT_X gate)
    // S gates won't change the rank of X (no point in using them)
    for i in 0..3_usize.pow(z_part.len() as u32) {
        let mut circuit = CliffordCircuit::new(0);
        let mut down: usize = i;
        for k in 0..z_part.len() {
            match down % 3 {
                1 => {
                    std::mem::swap(&mut z_part[k], &mut x_part[k]);
                    circuit.gates.push(CliffordGate::H(k))
                }
                2 => {
                    for l in 0..z_part[k].len() {
                        x_part[k][l] ^= z_part[k][l];
                    }
                    circuit.gates.push(CliffordGate::SqrtX(k))
                }
                _ => (),
            }
            down /= 3;
        }
        let rk = f2_rank(&x_part);
        if rk == x_part.len() - 1 {
            return Some(circuit);
        }
        let mut down: usize = i;
        for k in 0..z_part.len() {
            match down % 3 {
                1 => {
                    std::mem::swap(&mut z_part[k], &mut x_part[k]);
                    circuit.gates.push(CliffordGate::H(k))
                }
                2 => {
                    for l in 0..z_part[k].len() {
                        x_part[k][l] ^= z_part[k][l];
                    }
                    circuit.gates.push(CliffordGate::SqrtX(k))
                }
                _ => (),
            }
            down /= 3;
        }
    }
    None
}

fn update_rank_nm1(
    z_table: &mut Vec<Vec<bool>>,
    x_table: &mut Vec<Vec<bool>>,
    subset: &Vec<usize>,
) -> CliffordCircuit {
    for k in subset.iter() {
        let mut target = x_table[*k].clone();
        for j in subset.iter() {
            if *j != *k {
                vec_xor(&mut target, &x_table[*j]);
            }
        }
        if target.iter().filter(|a| **a).count() == 0 {
            let mut output = CliffordCircuit::new(z_table.len());
            for j in subset.iter() {
                if *j != *k {
                    rowop(x_table, *j, *k);
                    rowop(z_table, *k, *j);
                    output.gates.push(CliffordGate::CNOT(*j, *k));
                }
            }
            return output;
        }
    }
    panic!("This should never happen");
}

fn subset_codiag(
    z_table: &mut Vec<Vec<bool>>,
    x_table: &mut Vec<Vec<bool>>,
    subset: &Vec<usize>,
) -> Option<CliffordCircuit> {
    let mut z_part = Vec::new();
    let mut x_part = Vec::new();
    for q in subset.iter() {
        z_part.push(z_table[*q].clone());
    }
    for q in subset.iter() {
        x_part.push(x_table[*q].clone());
    }
    let circuit = find_rank_nm1(&mut z_part, &mut x_part);
    match circuit {
        None => return None,
        Some(circuit) => {
            let mut output_circuit = CliffordCircuit::new(z_table.len());
            for gate in circuit.gates.iter() {
                match gate {
                    CliffordGate::H(i) => {
                        std::mem::swap(&mut z_table[subset[*i]], &mut x_table[subset[*i]]);
                        output_circuit.gates.push(CliffordGate::H(subset[*i]));
                    }
                    CliffordGate::SqrtX(i) => {
                        for l in 0..z_part[*i].len() {
                            x_table[*i][l] ^= z_table[*i][l];
                        }
                        output_circuit.gates.push(CliffordGate::SqrtX(subset[*i]));
                    }
                    _ => (),
                }
            }
            let cnot_circuit = update_rank_nm1(z_table, x_table, subset);
            output_circuit.extend_with(&cnot_circuit);
            return Some(output_circuit);
        }
    }
}

fn is_qbit_trivial(
    z_table: &mut Vec<Vec<bool>>,
    x_table: &mut Vec<Vec<bool>>,
    qbit: usize,
) -> Option<CliffordCircuit> {
    let observed: HashSet<_> = z_table[qbit]
        .iter()
        .zip(x_table[qbit].iter())
        .map(|(a, b)| (*a, *b))
        .collect();
    if observed.len() < 2 {
        if observed.contains(&(true, true)) {
            vec_xor(&mut x_table[qbit], &z_table[qbit]);
            return Some(CliffordCircuit {
                nqbits: z_table.len(),
                gates: vec![CliffordGate::SqrtX(qbit)],
            });
        }
        if observed.contains(&(false, true)) {
            std::mem::swap(&mut x_table[qbit], &mut z_table[qbit]);
            return Some(CliffordCircuit {
                nqbits: z_table.len(),
                gates: vec![CliffordGate::H(qbit)],
            });
        }
        return Some(CliffordCircuit::new(z_table.len()));
    }
    return None;
}

fn enumerate_subsets(qubits: &HashSet<usize>, size: usize) -> Vec<Vec<usize>> {
    if size == 0 {
        return Vec::new();
    }
    let subsets = enumerate_subsets(qubits, size - 1);
    let mut new_subsets = Vec::new();
    for qbit in qubits.iter() {
        for subset in subsets.iter() {
            if subset.iter().filter(|e| **e > *qbit).count() > 0 {
                continue;
            }
            if !subset.contains(qbit) {
                let mut new_subset = subset.clone();
                new_subset.push(*qbit);
                new_subsets.push(new_subset);
            }
        }
    }
    new_subsets
}

fn step(
    z_table: &mut Vec<Vec<bool>>,
    x_table: &mut Vec<Vec<bool>>,
    qubits: &mut HashSet<usize>,
    k: usize,
) -> (usize, CliffordCircuit) {
    // Trivial case
    for qbit in qubits.iter() {
        if let Some(circuit) = is_qbit_trivial(z_table, x_table, *qbit) {
            return (*qbit, circuit);
        }
    }
    for subset_size in 2..k + 1 {
        let subsets = enumerate_subsets(qubits, subset_size);
        for subset in subsets {
            if let Some(circuit) = subset_codiag(z_table, x_table, &subset) {
                for qbit in subset.iter() {
                    if x_table[*qbit].iter().filter(|e| **e).count() == 0 {
                        return (*qbit, circuit);
                    }
                }
            }
        }
    }
    return greedy_reduce(z_table, x_table, &qubits.iter().map(|e| *e).collect());
}

pub fn codiagonalize_subsetwise(pauli_set: &PauliSet, k: usize) -> CliffordCircuit {
    let (mut z_table, mut x_table) = build_table(&pauli_set);
    let mut circuit = CliffordCircuit::new(pauli_set.n);
    let mut qbits: HashSet<usize> = (0..pauli_set.n).collect();
    while qbits.len() > 0 {
        let (qbit, piece) = step(&mut z_table, &mut x_table, &mut qbits, k);
        assert!(qbits.remove(&qbit));
        circuit.extend_with(&piece);
    }
    circuit
}

#[cfg(test)]
mod codiag_subset_tests {

    use super::*;
    use crate::structures::PauliLike;
    #[test]
    fn test_paper_ex() {
        let mut x_table = vec![
            vec![false, false, true, true, true, true, false, false],
            vec![false, false, true, true, true, true, false, false],
        ];
        let mut z_table = vec![
            vec![false, false, false, false, true, true, true, true],
            vec![false, true, false, true, false, true, false, true],
        ];
        let circuit = subset_codiag(&mut z_table, &mut x_table, &vec![0, 1]);
        println!("{:?}", circuit);
    }

    use rand::Rng;
    #[test]
    fn test_paper_ex_random_lc() {
        let mut rng = rand::thread_rng();
        let mut x_table = vec![
            vec![false, false, true, true, true, true, false, false],
            vec![false, false, true, true, true, true, false, false],
        ];
        let mut z_table = vec![
            vec![false, false, false, false, true, true, true, true],
            vec![false, true, false, true, false, true, false, true],
        ];
        for _ in 0..8 {
            for i in 0..2 {
                let pick = rng.gen::<u8>() % 3;
                if pick == 1 {
                    for l in 0..8 {
                        x_table[i][l] ^= z_table[i][l];
                    }
                }
                if pick == 2 {
                    std::mem::swap(&mut x_table[i], &mut z_table[i]);
                }
            }
        }
        let circuit = subset_codiag(&mut z_table, &mut x_table, &vec![0, 1]);
        assert!(matches!(circuit, Some(..)));
    }
    fn random_instance(n: usize, m: usize) -> PauliSet {
        let mut rng = rand::thread_rng();
        let mut pset = PauliSet::new(n);
        for _ in 0..m {
            let mut vec: Vec<bool> = Vec::new();
            for _ in 0..n {
                vec.push(rng.gen::<bool>());
            }
            for _ in 0..n {
                vec.push(false);
            }
            pset.insert_vec_bool(&vec, false);
        }
        for _ in 0..n * n {
            let i = rng.gen::<usize>() % n;
            loop {
                let j = rng.gen::<usize>() % n;
                if j != i {
                    pset.cnot(i, j);
                    let g2 = rng.gen::<bool>();
                    if g2 {
                        pset.h(j);
                    } else {
                        pset.s(j);
                    }
                    break;
                }
            }
            let g1 = rng.gen::<bool>();
            if g1 {
                pset.h(i);
            } else {
                pset.s(i);
            }
        }
        pset
    }

    #[test]
    fn test_thin() {
        for _ in 0..10 {
            let instance = random_instance(50, 20);
            let mut copy_instance = instance.clone();
            let circuit = codiagonalize_subsetwise(&instance, 2);
            copy_instance.conjugate_with_circuit(&circuit);
            for i in 0..instance.len() {
                let (_, vec) = copy_instance.get_as_vec_bool(i);
                for j in 0..instance.n {
                    assert!(!vec[j]);
                }
            }
        }
    }
    #[test]
    fn test_thick() {
        for _ in 0..10 {
            let instance = random_instance(20, 50);
            let mut copy_instance = instance.clone();
            let circuit = codiagonalize_subsetwise(&instance, 2);
            copy_instance.conjugate_with_circuit(&circuit);
            for i in 0..instance.len() {
                let (_, vec) = copy_instance.get_as_vec_bool(i);
                for j in 0..instance.n {
                    assert!(!vec[j]);
                }
            }
        }
    }
}
