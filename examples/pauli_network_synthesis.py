import numpy as np
from rustiq import pauli_network_synthesis, Metric

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
circuit = pauli_network_synthesis(pauli_sequence, Metric.COUNT, False)
print("Count no order:", cnot_count(circuit), cnot_depth(circuit))

# Optimizing DEPTH, not preserving order
circuit = pauli_network_synthesis(pauli_sequence, Metric.DEPTH, False)
print("Depth no order:", cnot_count(circuit), cnot_depth(circuit))

# Optimizing COUNT, preserving order
circuit = pauli_network_synthesis(pauli_sequence, Metric.COUNT, True)
print("Count order   :", cnot_count(circuit), cnot_depth(circuit))

# Optimizing DEPTH, preserving order
circuit = pauli_network_synthesis(pauli_sequence, Metric.DEPTH, True)
print("Depth order   :", cnot_count(circuit), cnot_depth(circuit))


# Pushing the optimization further by randomizing the qbit ordering:
circuit = pauli_network_synthesis(pauli_sequence, Metric.DEPTH, False, nshuffles=40)
print("# With 40 qubit shuffling (40 times slower):")
print("Depth no order:", cnot_count(circuit), cnot_depth(circuit))
