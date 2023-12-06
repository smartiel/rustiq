from rustiq import codiagonalization, Metric
import numpy as np


def ent_depth(circ):
    depths = {}
    for gate, qbits in circ:
        if gate in ["CNOT", "CZ"]:
            gate_depth = max(depths.get(qbits[0], 0), depths.get(qbits[1], 0)) + 1
            for qbit in qbits:
                depths[qbit] = gate_depth
    return max(depths.values() or [0])


def ent_count(circ):
    return sum(gate in ["CNOT", "CZ"] for gate, _ in circ)


PAULI = {
    (False, False): "I",
    (True, False): "Z",
    (False, True): "X",
    (True, True): "Y",
}


def generate_random_commuting(n, m):
    z_part = np.random.random((n, m)) < 0.5
    x_part = np.zeros((n, m), dtype=bool)
    for _ in range(n * n):
        for _ in range(n):
            if np.random.random() < 0.5:
                i = np.random.choice(range(n))
                x_part[i] ^= z_part[i]
            else:
                i = np.random.choice(range(n))
                x_part[i] ^= z_part[i]
                z_part[i] ^= x_part[i]
                x_part[i] ^= z_part[i]
        i, j = np.random.choice(range(n), size=2, replace=False)
        z_part[i] ^= z_part[j]
        x_part[j] ^= x_part[i]
    return [
        "".join(PAULI[(z_part[j][i], x_part[j][i])] for j in range(n)) for i in range(m)
    ]


instance = generate_random_commuting(60, 30)
circuit = codiagonalization(instance, Metric.COUNT)
print(ent_count(circuit), ent_depth(circuit))

circuit = codiagonalization(instance, Metric.DEPTH)
print(ent_count(circuit), ent_depth(circuit))


# Increasing the number of iterations in the syndrome decoding (only for count):
for niter in [1, 10, 20, 50, 100]:
    circuit = codiagonalization(instance, Metric.COUNT, niter)
    print(niter, ent_count(circuit))
