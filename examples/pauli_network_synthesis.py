import numpy as np
from rustiq import greedy_pauli_network, Metric

np.random.seed(1)

# Inputs to the synthesis algorithm are list of strings representing Pauli operators
n = 10
m = 100
pauli_sequence = list(
    set("".join(np.random.choice(["I", "X", "Y", "Z"], size=n)) for _ in range(m))
)


def cnot_depth(circ):
    depths = {}
    for gate, qbits in circ:
        if gate == "CNOT":
            gate_depth = max(depths.get(qbits[0], 0), depths.get(qbits[1], 0)) + 1
            for qbit in qbits:
                depths[qbit] = gate_depth
    return max(depths.values())


def cnot_count(circ):
    return sum(gate == "CNOT" for gate, _ in circ)


# Optimizing COUNT, not preserving order
circuit = greedy_pauli_network(pauli_sequence, Metric.COUNT, False, check=True)
print(cnot_count(circuit), cnot_depth(circuit))

# Optimizing DEPTH, not preserving order
circuit = greedy_pauli_network(pauli_sequence, Metric.DEPTH, False, check=True)
print(cnot_count(circuit), cnot_depth(circuit))

# Optimizing COUNT, preserving order
circuit = greedy_pauli_network(pauli_sequence, Metric.COUNT, True, check=True)
print(cnot_count(circuit), cnot_depth(circuit))

# Optimizing DEPTH, preserving order
circuit = greedy_pauli_network(pauli_sequence, Metric.DEPTH, True, check=True)
print(cnot_count(circuit), cnot_depth(circuit))
