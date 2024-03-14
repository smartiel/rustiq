"""
Utility functions for rustiq's circuit format
"""


def entangling_depth(circ):
    """
    Computes the entangling depth of a circuit.
    """
    depths = {}
    for _, qbits in circ:
        if len(qbits) > 1:
            gate_depth = max(depths.get(qbits[0], 0), depths.get(qbits[1], 0)) + 1
            for qbit in qbits:
                depths[qbit] = gate_depth
    return max(depths.values()) if depths else 0


def entangling_count(circ):
    """
    Count the number of entangling gates in a circuit.
    """
    return sum(len(qbits) == 2 for _, qbits in circ)
