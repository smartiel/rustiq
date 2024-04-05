/// Implementation of the Zhang et al algorithm for T-count optimization
/// The algorithm is adapted to work with any angles (not just pi/4 or pi/2)
use crate::structures::{Parameter, PauliLike, PauliSet, Tableau};

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

fn zhang_internal_no_simpl(
    rotations: &Vec<(String, Parameter)>,
    nqubits: usize,
) -> Vec<(String, Parameter)> {
    let mut bucket = PauliSet::new(nqubits);
    let mut angles: Vec<Parameter> = Vec::new();
    let mut merge_count = 0;
    for (axis, angle) in rotations.iter() {
        let index = bucket.insert(axis, false);
        let mut merged = false;
        for i in (0..index).rev() {
            if !bucket.commute(index, i) {
                break;
            }
            if bucket.equals(i, index) {
                angles[i] += angle.clone();
                merged = true;
                if angles[i].is_zero_mod_two_pi() {
                    bucket.set_to_identity(i);
                }
                break;
            }
        }
        if merged {
            bucket.pop_last();
            merge_count += 1;
        } else {
            angles.push(angle.clone());
        }
    }
    let mut output = Vec::new();
    for i in 0..bucket.len() {
        let (_, pstring) = bucket.get(i);
        if pstring.chars().any(|c| c != 'I') {
            output.push((pstring, angles[i].clone()));
        }
    }
    println!("# of merges: {merge_count}");
    output
}

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

pub fn full_initial_state_propagation(
    rotations: &Vec<(String, Parameter)>,
) -> Vec<(String, Parameter)> {
    let mut output = Vec::new();
    output
}
