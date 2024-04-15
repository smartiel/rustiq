/// A simple function that expresses a given circuit as a sequence of Pauli rotations
/// followed by a final Clifford operator
use crate::structures::{PauliLike, Tableau};

pub fn extract_rotations(
    circuit: &Vec<(String, Vec<usize>)>,
    nqubits: usize,
) -> (Vec<(bool, String)>, Tableau) {
    let mut clifford = Tableau::new(nqubits);
    let mut rotations = Vec::new();
    for (gate_name, qbits) in circuit.iter() {
        match gate_name.as_str() {
            "CX" => clifford.cnot(qbits[0], qbits[1]),
            "CZ" => clifford.cz(qbits[0], qbits[1]),
            "H" => clifford.h(qbits[0]),
            "S" => clifford.s(qbits[0]),
            "Sd" => clifford.sd(qbits[0]),
            "SqrtX" => clifford.sqrt_x(qbits[0]),
            "SqrtXd" => clifford.sqrt_xd(qbits[0]),
            "X" => {
                clifford.sqrt_x(qbits[0]);
                clifford.sqrt_x(qbits[0])
            }
            "Z" => {
                clifford.s(qbits[0]);
                clifford.s(qbits[0])
            }
            "RZ" => {
                rotations.push(clifford.get_inverse_z(qbits[0]));
            }
            _ => panic!("Unsupported gate {}", gate_name),
        }
    }
    (rotations, clifford)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rotations_to_circuit(rotations: &[(&str, f64)], n: usize) -> Vec<(String, Vec<usize>)> {
        let mut circuit = Vec::new();
        for (axis, _) in rotations {
            for q in 0..n {
                match axis.chars().nth(q).unwrap() {
                    'X' => circuit.push(("H".to_string(), vec![q])),
                    'Y' => circuit.push(("SqrtX".to_string(), vec![q])),
                    _ => {}
                }
            }
            let mut support: Vec<_> = (0..n)
                .filter(|q| axis.chars().nth(*q).unwrap() != 'I')
                .collect();
            let control = support.pop().unwrap();
            for q in support.iter() {
                circuit.push(("CX".to_string(), vec![*q, control]));
            }
            circuit.push(("RZ".to_string(), vec![control]));
            for q in support.iter() {
                circuit.push(("CX".to_string(), vec![*q, control]));
            }
            for q in 0..n {
                match axis.chars().nth(q).unwrap() {
                    'X' => circuit.push(("H".to_string(), vec![q])),
                    'Y' => circuit.push(("SqrtXd".to_string(), vec![q])),
                    _ => {}
                }
            }
        }
        circuit
    }

    #[test]
    fn simple_rotation_extraction() {
        let rotations = [
            ("ZX", 4.512201785772802),
            ("XX", 2.851130732235927),
            ("YZ", 6.202871194609474),
            ("ZY", 5.709017144348731),
            ("XI", 5.20162260375333),
            ("ZZ", 4.0981647318566905),
            ("IX", 1.20103093057112),
            ("XX", 3.748132932778756),
            ("IX", 2.7013063221306455),
            ("YI", 2.085429040429552),
        ];
        let circuit = rotations_to_circuit(&rotations, 2);
        let (new_rotations, clifford) = extract_rotations(&circuit, 2);
        println!("{:?}", new_rotations);
        assert_eq!(clifford, Tableau::new(2));
        for (r1, r2) in rotations.iter().zip(new_rotations.iter()) {
            assert_eq!(r1.0, r2.1);
            assert!(!r2.0);
        }
        println!("{}", clifford.logicals);
        for i in 0..4 {
            assert!(!clifford.logicals.get_phase(i));
        }
    }
}
