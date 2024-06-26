"""
A Rust powered Quantum Circuit synthesis library.

Algorithms are based on papers:

- "A graph-state based synthesis framework for Clifford isometries" 
  by T. Goubault de Brugière, S. Martiel, and C. Vuillot
- "Faster and shorter synthesis of Hamiltonian simulation circuits" 
  by T. Goubault de Brugère, and S. Martiel
"""

import typing

import networkx as nx

from .rustiq import (
    Metric,
    graph_state_synthesis as rust_gs,
    stabilizer_state_synthesis as rust_ss,
    pauli_network_synthesis as rust_pauli_network,
    codiagonalization as rust_codiag,
    codiagonalization_sswise as rust_sswise,
    isometry_synthesis as rust_isometry,
    extract_rotations as rust_extract_rotations,
    zhang_rotation_optimization,
)
from .utils import convert_circuit

__all__ = [
    "Metric",
    "pauli_network_synthesis",
    "graph_state_synthesis",
    "stabilizer_state_synthesis",
    "codiagonalization",
    "clifford_isometry_synthesis",
    "clifford_synthesis",
    "extract_rotations",
]

_Circuit: typing.TypeAlias = list[tuple[str, list[int]]]


def pauli_network_synthesis(
    paulis: list[str],
    metric: Metric,
    preserve_order: bool = True,
    nshuffles: int = 0,
    skip_sort: bool = False,
    fix_clifford: bool = False,
    check: bool = False,
) -> _Circuit:
    """
    Synthesize a circuit implementing a Pauli network for a given set of Pauli operators.

    Args:
        paulis (list): List of Pauli operators (as strings).
        metric (:class:`.Metric`): The metric to minimize.
        preserve_order (optional, bool): Whether to preserve the order of the Pauli operators.
          Default is True.
        nshuffles (optional, int): Number of qbit ordering shuffles to perform. Default is 0.
        skip_sort (optional, bool): Whether to skip sorting Pauli operators by hamming weight.
          If set to True, the algorithm might fail to converge (use with caution). Default is False.
        fix_clifford (optional, bool): Whether to fix the final Clifford operator. Default is False.
        check (optional, bool): Whether to check the validity of the circuit. Default is False.

    Returns:
        list: The synthesized circuit.
    """
    if paulis:
        return rust_pauli_network(
            paulis, metric, preserve_order, check, nshuffles, skip_sort, fix_clifford
        )
    return []


def graph_state_synthesis(
    graph: nx.Graph, metric: Metric, syndrome_iter: int = 1
) -> _Circuit:
    """
    Synthesize a circuit preparing a given graph state specified by a networkx graph.

    Args:
        graph (:class:`network.Graph`): A graph (networkx object).
        metric (:class:`.Metric`): The metric to minimize.
        syndrome_iter (optional, int): The number of syndrome decoding iteration used
          in the ISD solver (for `Metric.COUNT` only). Default is 1.
    Returns:
        list: The synthesized circuit.
    """
    syndrome_iter = syndrome_iter or 1
    adj = nx.to_numpy_array(graph) > 0.5
    return rust_gs(adj.tolist(), metric, syndrome_iter)


def stabilizer_state_synthesis(
    stabilizers: list[str], metric: Metric, syndrome_iter: int = 1
) -> _Circuit:
    """
    Synthesize a circuit preparing a given graph state specified by a networkx graph.

    Args:
        stabilizers (list): A list of strings representing the stabilizers.
        metric (:class:`.Metric`): The metric to minimize.
        syndrome_iter (optional, int): The number of syndrome decoding iteration used
          in the ISD solver (for `Metric.COUNT` only). Default is 1.
    Returns:
        list: The synthesized circuit.
    """
    return rust_ss(stabilizers, metric, syndrome_iter)


def codiagonalization(
    paulis: list[str], metric: Metric, syndrome_iter: int = 1
) -> _Circuit:
    """
    Synthesize a circuit preparing a codiagonalization Clifford circuit for a given set of
    pairwise commuting Pauli operators.

    Args:
        paulis (list): List of the Pauli operators to co-diagonalize (as strings).
        metric (:class:`.Metric`): The metric to minimize.
        syndrome_iter (optional, int): The number of syndrome decoding iteration used
          in the ISD solver (for `Metric.COUNT` only). Default is 1.

    Returns:
        list: The synthesized circuit.
    """
    syndrome_iter = syndrome_iter or 1
    return rust_codiag(paulis, metric, syndrome_iter)


def codiagonalization_sswise(paulis: list[str]) -> _Circuit:
    """
    Synthesize a circuit preparing a codiagonalization Clifford circuit for a given set of
    pairwise commuting Pauli operators.

    This is the heuristic of the paper "A Generic Compilation Strategy for the Unitary
    Coupled Cluster Ansatz".
    It is here for comparison purposes with the `codiagonalization` method .

    Args:
        paulis (list): List of the Pauli operators to co-diagonalize (as strings).

    Returns:
        list: The synthesized circuit.
    """
    return rust_sswise(paulis, 2)


def clifford_isometry_synthesis(
    logicals: list[str], stabilisers: list[str], metric: Metric, syndrome_iter: int = 1
) -> _Circuit:
    """
    Synthesize a Clifford circuit implementing a target Clifford isometry specified by
    its logical operators and its stabilizers.
    The tableau is specified by a list of operators (as strings). The first (resp. second) half
    of the list corresponding to the images of `Zi` (resp. `Xi`) operators.

    Args:
        logicals (list): The logical operators.
        stabilisers (list): The stabilizers.
        metric (:class:`.Metric`): The metric to minimize.
        syndrome_iter (optional, int): The number of syndrome decoding iteration used in the
          ISD solver in the metric=`Metric.COUNT` case only. Default is 1.

    Returns:
        list: The synthesized circuit.
    """
    syndrome_iter = syndrome_iter or 1
    return rust_isometry(logicals, stabilisers, metric, syndrome_iter)


def clifford_synthesis(
    logicals: list[str], metric: Metric, syndrome_iter: int = 1
) -> _Circuit:
    """
    Synthesize a Clifford circuit implementing a target Clifford operator specified by its
    tableau.
    The tableau is specified by a list of operators (as strings). The first (resp. second) half
    of the list corresponding to the images of `Zi` (resp. `Xi`) operators.

    Args:
        logicals (list): The tableau.
        metric (:class:`.Metric`): The metric to minimize.
        syndrome_iter (optional, int): The number of syndrome decoding iteration used in
          the ISD solver in the metric=`Metric.COUNT` case only. Default is 1.

    Returns:
        list: The synthesized circuit.
    """
    syndrome_iter = syndrome_iter or 1
    return rust_isometry(logicals, [], metric, syndrome_iter)


def extract_rotations(
    circuit: list[
        typing.Union[tuple[str, list[int]], tuple[str, list[int], typing.Any]]
    ],
    optimize: bool = True,
):
    """
    Converts a circuit into a sequence of Pauli rotations and a final Clifford operator.
    The circuit is specified as a sequence of tuples:
    - ("CX", [q0, q1]) for a CNOT gate
    - ("H", [q0]) for a Hadamard gate
    - ("S", [q0]) for a S gate
    - ("RZ", [q0], angle) for a RZ gate with angle in radian

    This gate set is exactly universal for QC.
    """
    circuit_no_angles, angles = convert_circuit(circuit)
    nqubits = max(max(g[1]) for g in circuit) + 1
    rotations, final_clifford_tableau = rust_extract_rotations(
        circuit_no_angles, [str(a) for a in angles], nqubits, optimize
    )
    return rotations, final_clifford_tableau
