"""
Utility functions for rustiq's circuits
"""

import typing

import numpy as np

_Circuit: typing.TypeAlias = list[tuple[str, list[int]]]


def entangling_depth(circ: _Circuit) -> int:
    """
    Computes the entangling depth of a circuit.
    Only checks the number of qubits of the gate to decide if it is entangling or not.

    Arguments:
        circ: A circuit.

    Returns:
        int: the entangling depth
    """
    depths = {}
    for _, qbits in circ:
        if len(qbits) > 1:
            gate_depth = max(depths.get(qbits[0], 0), depths.get(qbits[1], 0)) + 1
            for qbit in qbits:
                depths[qbit] = gate_depth
    return max(depths.values()) if depths else 0


def entangling_count(circ: _Circuit) -> int:
    """
    Count the number of entangling gates in a circuit.
    Only checks the number of qubits of the gate to decide if it is entangling or not.

    Arguments:
        circ: A circuit.

    Returns:
        int: the entangling count
    """
    return sum(len(qbits) == 2 for _, qbits in circ)


def convert_circuit(circuit: list) -> (list, list):
    """
    Converts a circuit into a CNOT + H + S + RZ(theta) circuit with theta non Clifford (i.e. != kpi/2)
    """
    new_circuit = []
    angles = []
    for gate in circuit:
        if gate[0] in ["CX", "CZ", "H", "S", "Sd", "X", "Z", "SqrtX", "SqrtXd"]:
            new_circuit.append(gate)
            continue
        assert gate[0] == "RZ"
        if isinstance(gate[2], str):
            new_circuit.append(gate[:2])
            angles.append(gate[2])
            continue
        if np.isclose(0.0, gate[2] % (np.pi / 2)) or np.isclose(
            np.pi / 2, gate[2] % (np.pi / 2)
        ):
            mult, offset = np.divmod(gate[2], np.pi / 2)
            if np.isclose(offset, np.pi / 2):
                mult += 1
            for _ in mult:
                new_circuit.append(("S", gate[1]))
            continue
        new_circuit.append(gate[:2])
        angles.append(gate[2])
    return new_circuit, angles
