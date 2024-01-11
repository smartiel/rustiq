use super::greedy_order_preserving::pauli_network_synthesis_no_permutation;
use super::greedy_pauli_network::pauli_network_synthesis;
use crate::structures::PauliLike;
use crate::structures::{CliffordCircuit, CliffordGate, IsometryTableau, Metric, PauliSet};
use crate::synthesis::clifford::isometry::isometry_synthesis;
use rand::thread_rng;
use rand::Rng;

fn permute_input(pset: &mut PauliSet) -> Vec<usize> {
    let mut rng = thread_rng();
    let mut permutation = (0..pset.n).collect::<Vec<usize>>();
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
            CliffordGate::CNOT(i, j) => output
                .gates
                .push(CliffordGate::CNOT(permutation[*i], permutation[*j])),
            CliffordGate::CZ(i, j) => output
                .gates
                .push(CliffordGate::CZ(permutation[*i], permutation[*j])),
        }
    }
    output
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
        let fix = isometry_synthesis(&mut tableau, &Metric::COUNT, 100);
        circuit.extend_with(&fix);
    }
    return circuit;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shuffle() {
        let mut operator_sequence = PauliSet::new(4);
        operator_sequence.insert("XX", false);
        let metric = Metric::COUNT;
        let preserve_order = true;
        let nshuffles = 50;

        let result = greedy_pauli_network(
            &mut operator_sequence,
            &metric,
            preserve_order,
            nshuffles,
            false,
            false,
        );
        assert_eq!(result.gates[0], CliffordGate::CNOT(1, 0));
    }
}
