from rustiq import codiagonalization, codiagonalization_sswise, Metric
from rustiq.utils import entangling_count, entangling_depth
import numpy as np


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


instance = generate_random_commuting(10, 6)
circuit = codiagonalization(instance, Metric.COUNT)
print("Count minimization: ", entangling_count(circuit), entangling_depth(circuit))

circuit = codiagonalization(instance, Metric.DEPTH)
print("Depth minimization: ", entangling_count(circuit), entangling_depth(circuit))

circuit = codiagonalization_sswise(instance)
print("Cowtan et al method:", entangling_count(circuit), entangling_depth(circuit))


print("Pushing the syndrome decoding:")
for niter in [1, 10, 20, 50, 100]:
    circuit = codiagonalization(instance, Metric.COUNT, niter)
    print(niter, entangling_count(circuit))
