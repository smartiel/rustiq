use crate::structures::pauli_dag::{build_dag_from_pauli_set, get_front_layer, Dag};
/// Implementation of the Zhang et al algorithm for T-count optimization
/// The algorithm is adapted to work with any angles (not just pi/4 or pi/2)
use crate::structures::{CliffordCircuit, CliffordGate, Parameter, PauliLike, PauliSet, Tableau};
use std::collections::HashSet;

fn update_rot_pi2<T: PauliLike>(axis: &String, k: i32, rest: &mut T, dagger: bool) {
    let mut support: Vec<_> = (0..axis.len())
        .filter(|i| axis.chars().nth(*i) != Some('I'))
        .collect();
    for qbit in support.iter() {
        match axis.chars().nth(*qbit).unwrap() {
            'X' => rest.h(*qbit),
            'Y' => {
                if dagger {
                    rest.sqrt_xd(*qbit)
                } else {
                    rest.sqrt_x(*qbit)
                }
            }
            _ => {}
        }
    }
    let target = support.pop().unwrap();
    for qbit in support.iter() {
        rest.cnot(*qbit, target);
    }
    for _ in 0..k {
        if dagger {
            rest.sd(target)
        } else {
            rest.s(target)
        }
    }
    for qbit in support.iter() {
        rest.cnot(*qbit, target);
    }
    for qbit in support.iter() {
        match axis.chars().nth(*qbit).unwrap() {
            'X' => rest.h(*qbit),
            'Y' => {
                if dagger {
                    rest.sqrt_x(*qbit)
                } else {
                    rest.sqrt_xd(*qbit)
                }
            }
            _ => {}
        }
    }
}

fn zhang_internal(
    rotations: &Vec<(String, Parameter)>,
    nqubits: usize,
    inverse_final_clifford: &mut Tableau,
) -> Vec<(String, Parameter)> {
    let mut bucket = PauliSet::new(nqubits);
    let mut axes = Vec::new();
    let mut future_angles = Vec::new();
    for (axis, angle) in rotations.iter() {
        axes.push(axis.clone());
        future_angles.push(angle.clone());
    }
    let mut rest = PauliSet::from_slice(&axes);
    let mut angles: Vec<Parameter> = Vec::new();
    let mut rot_index = 0;
    while rest.len() > 0 {
        let (phase, axis) = rest.get(0);
        rest.pop();
        let angle = &mut future_angles[rot_index];
        if phase {
            angle.flip_sign();
        }
        let index = bucket.insert(&axis, false);
        let mut merged = false;
        for i in (0..index).rev() {
            if !bucket.commute(index, i) {
                break;
            }
            if bucket.equals(i, index) {
                angles[i] += angle.clone();
                merged = true;
                let (new_angle, mult_pi_2) = angles[i].simplify();
                update_rot_pi2(&axis, mult_pi_2, &mut rest, true);
                update_rot_pi2(&axis, mult_pi_2, inverse_final_clifford, true);
                angles[i] = new_angle;
                if angles[i].is_zero_mod_two_pi() {
                    bucket.set_to_identity(i);
                }
                break;
            }
        }
        if merged {
            bucket.pop_last();
        } else {
            angles.push(angle.clone());
        }
        rot_index += 1;
    }
    let mut output = Vec::new();
    for i in 0..bucket.len() {
        let (phase, pstring) = bucket.get(i);
        assert!(!phase);
        if pstring.chars().any(|c| c != 'I') {
            output.push((pstring, angles[i].clone()));
        }
    }
    output
}

/// Implementation of the Zhang et al algorithm for T-count optimization adapted to work
/// with any angles (even parametrized rotations)
pub fn zhang_rotation_optimization(
    mut rotations: Vec<(String, Parameter)>,
    nqubits: usize,
) -> (Vec<(String, Parameter)>, Tableau) {
    let mut current_size = rotations.len();
    let mut inverse_final_clifford = Tableau::new(nqubits);
    loop {
        let new_rotations = zhang_internal(&rotations, nqubits, &mut inverse_final_clifford);
        if new_rotations.len() == current_size {
            break;
        }
        current_size = new_rotations.len();
        rotations = new_rotations;
    }
    (rotations, inverse_final_clifford)
}

/// Data structure of the initial state propagation
/// It stores a Pauli DAG and a set of marked qubits (i.e. qubits that are not stabilized by |0> anymore)
struct MarkedPauliDag {
    pauli_set: PauliSet,
    final_clifford: Tableau,
    dag: Dag,
    marked: HashSet<usize>,
    output_pauli_set: PauliSet,
    output_rotations: Vec<usize>,
    did_something: bool,
}

impl MarkedPauliDag {
    fn new(pauli_set: PauliSet) -> Self {
        let dag = build_dag_from_pauli_set(&pauli_set);
        Self {
            output_pauli_set: PauliSet::new(pauli_set.n),
            final_clifford: Tableau::new(pauli_set.n),
            pauli_set,
            dag,
            marked: HashSet::new(),
            output_rotations: Vec::new(),
            did_something: false,
        }
    }

    fn get_front_indices(&self) -> Vec<usize> {
        let front_layer = get_front_layer(&self.dag);
        front_layer
            .into_iter()
            .map(|ni| *self.dag.node_weight(ni).unwrap())
            .collect()
    }

    fn get_unmarked_xy(&self, rotation_index: usize) -> Option<usize> {
        for qbit in 0..self.pauli_set.n {
            if self.pauli_set.get_entry(qbit, rotation_index) & !self.marked.contains(&qbit) {
                return Some(qbit);
            }
        }
        return None;
    }
    fn has_unmarked_xy(&self, rotation_index: usize) -> bool {
        for qbit in 0..self.pauli_set.n {
            if self.pauli_set.get_entry(qbit, rotation_index) & !self.marked.contains(&qbit) {
                return true;
            }
        }
        return false;
    }
    fn get_rotation_score(&self, rotation_index: usize) -> usize {
        if self.has_unmarked_xy(rotation_index) {
            return self.pauli_set.support_size(rotation_index) - 1;
        }
        let mut score = 0;
        for qbit in 0..self.pauli_set.n {
            if !self.pauli_set.get_entry(qbit, rotation_index)
                & self
                    .pauli_set
                    .get_entry(qbit + self.pauli_set.n, rotation_index)
                & !self.marked.contains(&qbit)
            {
                score += 1;
            }
        }
        return score;
    }
    /// Optimize a given rotation from the front layer of the PDAG
    fn optimize_rotation(&mut self, rotation_index: usize) {
        // Remove unmarked Z components:
        for qbit in 0..self.pauli_set.n {
            if !self.marked.contains(&qbit) {
                if self
                    .pauli_set
                    .get_entry(qbit + self.pauli_set.n, rotation_index)
                    & !self.pauli_set.get_entry(qbit, rotation_index)
                {
                    self.pauli_set.set_entry(rotation_index, qbit, false, false);
                    self.did_something = true;
                }
            }
        }
        // Fold X components onto a single unmarked qubit (if there is any)
        let mut piece = CliffordCircuit::new(self.pauli_set.n);
        if let Some(control) = self.get_unmarked_xy(rotation_index) {
            for qbit in self.pauli_set.get_support(rotation_index) {
                if qbit == control {
                    continue;
                }
                assert!(
                    self.marked.contains(&qbit) | self.pauli_set.get_entry(qbit, rotation_index)
                );
                if self
                    .pauli_set
                    .get_entry(qbit + self.pauli_set.n, rotation_index)
                {
                    if self.pauli_set.get_entry(qbit, rotation_index) {
                        piece.gates.push(CliffordGate::Sd(qbit));
                    } else {
                        piece.gates.push(CliffordGate::H(qbit));
                    }
                }
                piece.gates.push(CliffordGate::CNOT(control, qbit));
            }
            for qbit in self.pauli_set.get_support(rotation_index) {
                if qbit == control {
                    continue;
                }
                if self
                    .pauli_set
                    .get_entry(qbit + self.pauli_set.n, rotation_index)
                {
                    if self.pauli_set.get_entry(qbit, rotation_index) {
                        piece.gates.push(CliffordGate::S(qbit));
                    } else {
                        piece.gates.push(CliffordGate::H(qbit));
                    }
                }
                self.pauli_set.set_entry(rotation_index, qbit, false, false);
                self.did_something = true;
            }
            self.marked.insert(control);
        }
        if self.pauli_set.support_size(rotation_index) > 0 {
            let (phase, bv) = self.pauli_set.get_as_vec_bool(rotation_index);
            self.output_pauli_set.insert_vec_bool(&bv, phase);
            self.pauli_set.conjugate_with_circuit(&piece.dagger());
            let c = Tableau::from_circuit(&piece);
            self.final_clifford = self.final_clifford.clone() * c;
            self.output_rotations.push(rotation_index);
        }
        self.dag.retain_nodes(|graph, node_index| {
            return *graph.node_weight(node_index).unwrap() != rotation_index;
        });
    }

    /// Attempts to simplify a rotation from the front layer of the PDAG
    /// Returns true if a simplification was made
    fn simplify_once(&mut self) -> bool {
        let front_layer = self.get_front_indices();
        // Looking for the best rotation to remove
        // score(R) = |R| -1 if there is an unmarked X/Y
        // score(R) = # of unmarked Zs otherwise
        let best_candidate = front_layer
            .into_iter()
            .map(|ri| (ri, self.get_rotation_score(ri)))
            .max_by_key(|(_, score)| *score);
        if let Some((rotation_index, score)) = best_candidate {
            if score > 0 {
                self.optimize_rotation(rotation_index);
                return true;
            }
        }
        return false;
    }

    fn pop_rest(&mut self) {
        loop {
            let front_layer = self.get_front_indices();
            if front_layer.len() == 0 {
                break;
            }
            for rotation_index in front_layer.iter() {
                let (phase, bv) = self.pauli_set.get_as_vec_bool(*rotation_index);
                self.output_pauli_set.insert_vec_bool(&bv, phase);
                self.output_rotations.push(*rotation_index);
            }
            self.dag.retain_nodes(|graph, node_index| {
                return !front_layer.contains(graph.node_weight(node_index).unwrap());
            })
        }
    }

    pub fn propagate(mut self) -> (PauliSet, Tableau, Vec<usize>, bool) {
        while (&mut self).simplify_once() {}
        self.pop_rest();
        return (
            self.output_pauli_set,
            self.final_clifford,
            self.output_rotations,
            self.did_something,
        );
    }
}

pub fn full_initial_state_propagation(
    rotations: &Vec<(String, Parameter)>,
) -> (Vec<(String, Parameter)>, Tableau) {
    let axes: Vec<_> = rotations.iter().map(|e| e.0.clone()).collect();
    let mut angles: Vec<_> = rotations.iter().map(|e| e.1.clone()).collect();
    let mut pset = PauliSet::from_slice(&axes);
    let mut final_clifford = Tableau::new(pset.n);
    loop {
        let mpdag = MarkedPauliDag::new(pset);
        let (new_pset, clifford, new_rotations, carry_on) = mpdag.propagate();
        angles = new_rotations
            .into_iter()
            .map(|i| angles[i].clone())
            .collect();
        pset = new_pset;
        final_clifford = final_clifford * clifford;
        if !carry_on {
            break;
        }
    }
    let mut new_rotations = Vec::new();
    for i in 0..pset.len() {
        let (phase, string) = pset.get(i);
        if phase {
            angles[i].flip_sign();
        }
        new_rotations.push((string, angles[i].clone()));
    }
    return (new_rotations, final_clifford);
}
