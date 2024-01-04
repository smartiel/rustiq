import numpy as np
from rustiq import pauli_network_synthesis, Metric
from rustiq.utils import entangling_count, entangling_depth

np.random.seed(1)

# Inputs to the synthesis algorithm are list of strings representing Pauli operators
n = 10
m = 100
pauli_sequence = list(
    set("".join(np.random.choice(["I", "X", "Y", "Z"], size=n)) for _ in range(m))
)


# Optimizing COUNT, not preserving order
circuit = pauli_network_synthesis(pauli_sequence, Metric.COUNT, False)
print("Count no order:", entangling_count(circuit), entangling_depth(circuit))

# Optimizing DEPTH, not preserving order
circuit = pauli_network_synthesis(pauli_sequence, Metric.DEPTH, False)
print("Depth no order:", entangling_count(circuit), entangling_depth(circuit))

# Optimizing COUNT, preserving order
circuit = pauli_network_synthesis(pauli_sequence, Metric.COUNT, True)
print("Count order   :", entangling_count(circuit), entangling_depth(circuit))

# Optimizing DEPTH, preserving order
circuit = pauli_network_synthesis(pauli_sequence, Metric.DEPTH, True)
print("Depth order   :", entangling_count(circuit), entangling_depth(circuit))


# Optimizing COUNT, not preserving order, fixing the final clifford operator
circuit = pauli_network_synthesis(
    pauli_sequence, Metric.COUNT, False, fix_clifford=True
)
print("Count no order (fixed):", entangling_count(circuit), entangling_depth(circuit))

# Optimizing DEPTH, not preserving order, fixing the final clifford operator
circuit = pauli_network_synthesis(
    pauli_sequence, Metric.DEPTH, False, fix_clifford=True
)
print("Depth no order (fixed):", entangling_count(circuit), entangling_depth(circuit))

# Optimizing COUNT, preserving order, fixing the final clifford operator
circuit = pauli_network_synthesis(pauli_sequence, Metric.COUNT, True, fix_clifford=True)
print("Count order    (fixed):", entangling_count(circuit), entangling_depth(circuit))

# Optimizing DEPTH, preserving order, fixing the final clifford operator
circuit = pauli_network_synthesis(pauli_sequence, Metric.DEPTH, True, fix_clifford=True)
print("Depth order    (fixed):", entangling_count(circuit), entangling_depth(circuit))

# Pushing the optimization further by randomizing the qbit ordering:
print("# With 40 qubit shuffling (40 times slower):")
circuit = pauli_network_synthesis(pauli_sequence, Metric.DEPTH, False, nshuffles=40)
print("Depth no order:", entangling_count(circuit), entangling_depth(circuit))
