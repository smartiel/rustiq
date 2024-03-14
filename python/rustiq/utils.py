"""
Utility functions for rustiq's circuits
"""

import typing

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
