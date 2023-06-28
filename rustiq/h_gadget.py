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


def gadgetize(network):
    nqbits = len(network[0][1])
    mapping = dict(zip(range(nqbits), range(nqbits)))
    corrections = {}
    fresh = nqbits
    measurement_index = 0
    beginning = []
    phase_polynomial_part = []
    end = []
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
        new_rotation = [
            (pauli, mapping[qbit])
            for qbit, pauli in enumerate(rotation)
            if pauli != "I"
        ]
        phase_polynomial_part.append(
            ("ROTATION", new_rotation, corrections.get(index, []))
        )
    return beginning, phase_polynomial_part, end
