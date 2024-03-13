# Example: build an encoding circuit for the Steane code
from rustiq import clifford_isometry_synthesis, Metric
from rustiq.utils import entangling_count, entangling_depth

stabilizers = ["IIIXXXX", "IXXIIXX", "XIXIXIX", "IIIZZZZ", "IZZIIZZ", "ZIZIZIZ"]
logical_operators = ["ZZZZZZZ", "XXXXXXX"]

circuit_count = clifford_isometry_synthesis(
    logical_operators, stabilizers, Metric.COUNT, syndrome_iter=1000
)
print(entangling_count(circuit_count), entangling_depth(circuit_count))

circuit_depth = clifford_isometry_synthesis(
    logical_operators, stabilizers, Metric.DEPTH
)
print(entangling_count(circuit_depth), entangling_depth(circuit_depth))
