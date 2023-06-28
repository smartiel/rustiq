import numpy as np
from .rustiq import diagonalization_network


class Pauli:
    def __init__(self, pstring):
        self.nqbits = len(pstring)
        self.x_vec = np.zeros(self.nqbits, dtype=np.byte)
        self.z_vec = np.zeros(self.nqbits, dtype=np.byte)
        for i, pauli in enumerate(pstring):
            if pauli in "ZY":
                self.z_vec[i] = 1
            if pauli in "YX":
                self.x_vec[i] = 1

    def commutes(self, other_pauli):
        return (
            sum((self.x_vec & other_pauli.z_vec) ^ (self.z_vec & other_pauli.x_vec)) % 2
        ) == 0

    def __repr__(self):
        pstring_out = ""
        for i in range(self.nqbits):
            if self.z_vec[i] & self.x_vec[i]:
                pstring_out += "Y"
            elif self.z_vec[i]:
                pstring_out += "Z"
            elif self.x_vec[i]:
                pstring_out += "X"
            else:
                pstring_out += "I"
        return pstring_out

    def conjugate_with_gate(self, gate):
        if gate[0] == "H":
            self.x_vec[gate[1][0]] ^= self.z_vec[gate[1][0]]
            self.z_vec[gate[1][0]] ^= self.x_vec[gate[1][0]]
            self.x_vec[gate[1][0]] ^= self.z_vec[gate[1][0]]
        elif gate[0] == "S":
            self.z_vec[gate[1][0]] ^= self.x_vec[gate[1][0]]
        elif gate[0] == "CNOT":
            self.z_vec[gate[1][0]] ^= self.z_vec[gate[1][1]]
            self.x_vec[gate[1][1]] ^= self.x_vec[gate[1][0]]
        else:
            raise ValueError(f"Unknown gate {gate[0]}")

    def conjugate_with(self, circuit):
        for gate in circuit:
            self.conjugate_with_gate(gate)


def propagate_x_correction(network, start_index, target_qbit, current_mapping):
    correction = Pauli(
        "".join(
            "I" if qbit != target_qbit else "X"
            for qbit in range(max(current_mapping.values()) + 1)
        )
    )
    # We always anti-commute with the next rotation
    # (this is why we had a H gate in the first place)
    to_flip = [start_index]
    for index, (piece, rotation) in enumerate(network[start_index + 1 :]):
        correction.conjugate_with(piece)
        axis = ["I"] * (max(current_mapping.values()) + 1)
        for qbit, pauli in enumerate(rotation):
            axis[current_mapping[qbit]] = pauli
        rotation = Pauli(axis)
        if not rotation.commutes(correction):
            to_flip.append(index + start_index + 1)
    return to_flip, correction


def simplify_phase_polynomial(phase_polynomial_part, nqbits):
    phase_poly = []
    clifford_correction = []
    for gate in reversed(phase_polynomial_part):
        # Gate can be:
        # - a (diagonal) rotation
        # - a CNOT gate ====> here we need to do a bit of work
        # - a S gate (commutes trivially with any rotation)
        # - a CZ gate (commutes trivially with any rotation)
        if gate[0] == "ROTATION":
            z_vec = np.zeros((nqbits,), dtype=np.byte)
            for i, p in enumerate(gate[1]):
                if p == "Z":
                    z_vec[i] = 1
            phase_poly.append((z_vec, gate[2]))
            continue
        elif gate[0] == "CNOT":
            for rotation in phase_poly:
                rotation[gate[1][0]] ^= rotation[gate[1][1]]
        clifford_correction.append(gate)

    return phase_poly, list(reversed(clifford_correction))


def gadgetize(network):
    nqbits = len(network[0][1])
    mapping = dict(zip(range(nqbits), range(nqbits)))
    corrections = {}
    fresh = nqbits
    measurement_index = 0
    beginning = []
    phase_polynomial_part = []
    end = []
    predicted_nqbits = nqbits
    for piece, _ in network:
        for gate in piece:
            if gate[0] == "H":
                predicted_nqbits += 1

    for index, (circuit_piece, rotation) in enumerate(network):
        for gate, qbits in circuit_piece:
            if gate == "H":
                beginning.append(("H", [fresh]))
                phase_polynomial_part.append(("CZ", [mapping[qbits[0]], fresh]))
                end.append(("H", [mapping[qbits[0]]]))
                end.append(("MEASUREMENT", [mapping[qbits[0]]], measurement_index))
                mapping[qbits[0]] = fresh
                to_flip, _ = propagate_x_correction(network, index, fresh, mapping)
                for rot_index in to_flip:
                    corr = corrections.get(rot_index, [])
                    corr.append(measurement_index)
                    corrections[rot_index] = corr
                measurement_index += 1
                fresh += 1
            else:
                phase_polynomial_part.append((gate, [mapping[qbit] for qbit in qbits]))
        new_rotation = ["I"] * predicted_nqbits
        for qbit, pauli in enumerate(rotation):
            new_rotation[mapping[qbit]] = pauli
        phase_polynomial_part.append(
            ("ROTATION", "".join(new_rotation), corrections.get(index, []))
        )
    return (
        beginning,
        simplify_phase_polynomial(phase_polynomial_part, predicted_nqbits),
        end,
    )


def predict_nbqbits(pauli_sequence, optimal=False, display=False):
    """
    Predict the total number of qubits required in order
    to implement this sequence while gadgetizing every H gate.
    """
    init_circuit, network = diagonalization_network(pauli_sequence, optimal)
    if display:
        print(init_circuit)
        print(network)
    nqbits = len(network[0][1])
    predicted_nqbits = nqbits
    for piece, _ in network:
        for gate in piece:
            if gate[0] == "H":
                predicted_nqbits += 1
    return predicted_nqbits