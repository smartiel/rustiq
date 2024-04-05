/// A simple function that expresses a given circuit as a sequence of Pauli rotations
/// followed by a final Clifford operator
use crate::structures::{PauliLike, Tableau};

pub fn extract_rotations(
    circuit: &Vec<(String, Vec<usize>)>,
    nqubits: usize,
) -> (Vec<(bool, String)>, Tableau) {
    let mut isometry = Tableau::new(nqubits);
    let mut rotations = Vec::new();
    for (gate_name, qbits) in circuit.iter() {
        match gate_name.as_str() {
            "CX" => isometry.cnot(qbits[0], qbits[1]),
            "H" => isometry.h(qbits[0]),
            "S" => isometry.s(qbits[0]),
            "RZ" => {
                rotations.push(isometry.logicals.get_inverse_z(qbits[0]));
            }
            _ => panic!("Unsupported gate {}", gate_name),
        }
    }
    (rotations, isometry)
}
