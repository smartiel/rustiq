use std::collections::HashSet;

use super::greedy_order_preserving::pauli_network_synthesis_no_permutation;
use super::greedy_pauli_network::pauli_network_synthesis;
use crate::structures::{
    CliffordCircuit, CliffordGate, IsometryTableau, Metric, PauliLike, PauliSet,
};
use crate::synthesis::clifford::isometry::isometry_synthesis;
use rand::thread_rng;
use rand::Rng;

fn permute_input(pset: &mut PauliSet) -> Vec<usize> {
    let mut rng = thread_rng();
    let mut permutation = (0..pset.n).collect::<Vec<usize>>();
    if pset.n <= 1 {
        return permutation;
    }
    for _ in 0..pset.n * pset.n {
        let i = rng.gen_range(0..pset.n - 1);
        let j = rng.gen_range(i + 1..pset.n);
        permutation.swap(i, j);
        pset.swap_qbits(i, j);
    }
    permutation
}

fn permute_circuit(circuit: &CliffordCircuit, permutation: &Vec<usize>) -> CliffordCircuit {
    let mut output = CliffordCircuit::new(circuit.nqbits);
    for gate in circuit.gates.iter() {
        match gate {
            CliffordGate::H(i) => output.gates.push(CliffordGate::H(permutation[*i])),
            CliffordGate::S(i) => output.gates.push(CliffordGate::S(permutation[*i])),
            CliffordGate::Sd(i) => output.gates.push(CliffordGate::Sd(permutation[*i])),
            CliffordGate::SqrtX(i) => output.gates.push(CliffordGate::SqrtX(permutation[*i])),
            CliffordGate::SqrtXd(i) => output.gates.push(CliffordGate::SqrtXd(permutation[*i])),
            CliffordGate::CNOT(a, b) => output
                .gates
                .push(CliffordGate::CNOT(permutation[*a], permutation[*b])),
            CliffordGate::CZ(a, b) => output
                .gates
                .push(CliffordGate::CZ(permutation[*a], permutation[*b])),
        }
    }
    output
}

pub fn check_circuit(input: &[String], circuit: &CliffordCircuit) {
    let mut hit_map: HashSet<usize> = HashSet::new();
    let mut bucket = PauliSet::from_slice(input);
    for i in 0..bucket.len() {
        if bucket.support_size(i) == 1 {
            hit_map.insert(i);
        }
    }
    for gate in circuit.gates.iter() {
        bucket.conjugate_with_gate(&gate);

        for i in 0..bucket.len() {
            if bucket.support_size(i) == 1 {
                hit_map.insert(i);
            }
        }
    }
    assert!(
        hit_map.len() == input.len(),
        "Synthesized {} operators, expected {}",
        hit_map.len(),
        input.len()
    );
}

pub fn greedy_pauli_network(
    operator_sequence: &mut PauliSet,
    metric: &Metric,
    preserve_order: bool,
    nshuffles: usize,
    skip_sort: bool,
    fix_clifford: bool,
) -> CliffordCircuit {
    let synth = if preserve_order {
        pauli_network_synthesis_no_permutation
    } else {
        pauli_network_synthesis
    };
    let mut circuit = synth(&mut operator_sequence.clone(), metric, skip_sort);
    let mut cost = metric.on_circuit(&circuit);
    for _ in 0..nshuffles {
        let mut pset = operator_sequence.clone();
        let permutation = permute_input(&mut pset);
        let new_circuit = synth(&mut pset, metric, skip_sort);
        let new_cost = metric.on_circuit(&new_circuit);
        if new_cost < cost {
            cost = new_cost;
            circuit = permute_circuit(&new_circuit, &permutation);
        }
    }
    if fix_clifford {
        let mut tableau = IsometryTableau::new(circuit.nqbits, 0);
        tableau.conjugate_with_circuit(&circuit.dagger());
        let fix = isometry_synthesis(&mut tableau, &metric, 100);
        circuit.extend_with(&fix);
    }
    return circuit;
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_inifinte_loop_identities() {
        let mut operator_sequence = PauliSet::new(4);
        operator_sequence.insert("XZYX", false);
        operator_sequence.insert("XXIY", false);
        operator_sequence.insert("IIII", false);
        operator_sequence.insert("ZZYI", false);
        operator_sequence.insert("XXZZ", false);
        operator_sequence.insert("ZYZY", false);

        let metric = Metric::DEPTH;
        let preserve_order = true;
        let nshuffles = 0;

        let _result = greedy_pauli_network(
            &mut operator_sequence,
            &metric,
            preserve_order,
            nshuffles,
            false,
            true,
        );
    }
    #[test]
    fn test_fix_clifford() {
        let mut operator_sequence = PauliSet::new(4);
        operator_sequence.insert("XZYX", false);
        operator_sequence.insert("XXIY", false);
        operator_sequence.insert("ZZYI", false);
        operator_sequence.insert("XXZZ", false);
        operator_sequence.insert("ZYZY", false);

        let metric = Metric::COUNT;
        let preserve_order = true;
        let nshuffles = 0;

        let result = greedy_pauli_network(
            &mut operator_sequence,
            &metric,
            preserve_order,
            nshuffles,
            false,
            true,
        );
        check_circuit(
            &[
                "XZYX".to_owned(),
                "XXIY".to_owned(),
                "ZZYI".to_owned(),
                "XXZZ".to_owned(),
                "ZYZY".to_owned(),
            ],
            &result,
        );
        let mut clifford = IsometryTableau::new(4, 0);
        clifford.conjugate_with_circuit(&result);
        for i in 0..clifford.n {
            let mut st = "".to_owned();
            for q in 0..clifford.n {
                if q != i {
                    st += "I";
                } else {
                    st += "X";
                }
            }
            let (_, t) = clifford.logicals.get(i);
            assert_eq!(t, st);
            let mut st = "".to_owned();
            for q in 0..clifford.n {
                if q != i {
                    st += "I";
                } else {
                    st += "Z";
                }
            }
            let (_, t) = clifford.logicals.get(i + clifford.n);
            assert_eq!(t, st);
        }
    }
    #[test]
    fn test_shuffle() {
        let mut operator_sequence = PauliSet::new(4);
        operator_sequence.insert("XZYX", false);
        operator_sequence.insert("XXIY", false);
        operator_sequence.insert("ZZYI", false);
        operator_sequence.insert("XZZZ", false);
        operator_sequence.insert("ZYZY", false);

        let metric = Metric::COUNT;
        let preserve_order = true;
        let nshuffles = 10;

        let result = greedy_pauli_network(
            &mut operator_sequence,
            &metric,
            preserve_order,
            nshuffles,
            false,
            false,
        );
        check_circuit(
            &[
                "XZYX".to_owned(),
                "XXIY".to_owned(),
                "ZZYI".to_owned(),
                "XZZZ".to_owned(),
                "ZYZY".to_owned(),
            ],
            &result,
        );
    }
}
