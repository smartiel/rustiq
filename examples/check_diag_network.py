import numpy as np
from qiskit.circuit import QuantumCircuit

with open("uccsd_0.txt", "r") as fin:
    paulis = [line.strip("\n") for line in fin.readlines()]

# paulis = ["XX", "ZZ"]


def trivial_impl(paulis, angles):
    nqbits = len(paulis[0])
    circuit = QuantumCircuit(nqbits)
    for rotation, angle in zip(paulis, angles):
        for qbit, pauli in enumerate(rotation):
            if pauli == "X":
                circuit.h(qbit)
            elif pauli == "Y":
                circuit.rx(np.pi / 2, qbit)
        support = [qbit for qbit, pauli in enumerate(rotation) if pauli != "I"]
        for control in support[1:]:
            circuit.cx(control, support[0])
        circuit.rz(angle, support[0])
        for control in support[1:]:
            circuit.cx(control, support[0])
        for qbit, pauli in enumerate(rotation):
            if pauli == "X":
                circuit.h(qbit)
            elif pauli == "Y":
                circuit.rx(-np.pi / 2, qbit)
    return circuit


angles = np.random.random((len(paulis),)) * 2 * np.pi

from rustiq.h_gadget import network_to_qiskit, diagonalization_network

init, network = diagonalization_network(paulis, False)
circuit_diag_1 = network_to_qiskit(init, network, angles)

init, network = diagonalization_network(paulis, True)
circuit_diag_2 = network_to_qiskit(init, network, angles)
circuit_std = trivial_impl(paulis, angles)

from qiskit import Aer

backend = Aer.get_backend("statevector_simulator")
statevectors = [
    (backend.run(circuit).result().get_statevector(circuit, decimals=10).data)
    for circuit in (circuit_diag_1, circuit_diag_2, circuit_std)
]

ref = statevectors[2]

for sv in statevectors[:-1]:
    print(np.abs(sv.dot(ref.conj())))
