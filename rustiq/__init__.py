import networkx as nx
from .rustiq import (
    Metric,
    graph_state_synthesis as rust_gs,
    pauli_network_synthesis as rust_pauli_network,
    codiagonalization as rust_codiag,
    codiagonalization_sswise as rust_sswise,
    isometry_synthesis as rust_isometry,
)

__all__ = [
    "Metric",
    "graph_state_synthesis",
    "pauli_network_synthesis",
    "codiagonalization",
    "clifford_synthesis",
]


def pauli_network_synthesis(
    paulis,
    metric,
    preserve_order=True,
    nshuffles=0,
    skip_sort=False,
    fix_clifford=False,
    check=False,
):
    """
    Synthesize a circuit implementing a Pauli network for a given set of Pauli operators.

    Args:
        paulis (list): List of Pauli operators (as strings).
        metric (Metric): The metric to minimize.
        preserve_order (optional, bool): Whether to preserve the order of the Pauli operators. Default is True.
        nshuffles (optional, int): Number of qbit ordering shuffles to perform. Default is 0.
        skip_sort (optional, bool): Whether to skip the sorting of the Pauli operators (by hamming weight).
          If set to True, the algorithm might fail to converge (use with caution). Default is False.
        fix_clifford (optional, bool): Whether to fix the final Clifford operator. Default is False.
        check (optional, bool): Whether to check the validity of the circuit. Default is False.

    Returns:
        list: The synthesized circuit.
    """
    return rust_pauli_network(
        paulis, metric, preserve_order, check, nshuffles, skip_sort, fix_clifford
    )


def graph_state_synthesis(graph, metric, syndrome_iter=1):
    """
    Synthesize a circuit preparing a given graph state specified by a networkx graph.

    Args:
        graph (nx.Graph): A graph (networkx object).
        metric (Metric): The metric to minimize.
        syndrome_iter (optional, int): The number of syndrome decoding iteration used in the ISD solver
          (for `Metric.COUNT` only). Default is 1.
    Returns:
        list: The synthesized circuit.
    """
    syndrome_iter = syndrome_iter or 1
    adj = nx.to_numpy_array(graph) < 0.5
    return rust_gs(adj.tolist(), metric, syndrome_iter)


def codiagonalization(paulis, metric, syndrome_iter=1):
    """
    Synthesize a circuit preparing a codiagonalization Clifford circuit for a given set of
    pairwise commuting Pauli operators.

    Args:
        paulis (list): List of the Pauli operators to co-diagonalize (as strings).
        metric (Metric): The metric to minimize.
        syndrome_iter (optional, int): The number of syndrome decoding iteration used in the ISD solver
          (for `Metric.COUNT` only). Default is 1.

    Returns:
        list: The synthesized circuit.
    """
    syndrome_iter = syndrome_iter or 1
    return rust_codiag(paulis, metric, syndrome_iter)


def codiagonalization_sswise(paulis):
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


def clifford_isometry_synthesis(logicals, stabilisers, metric, syndrome_iter=1):
    """
    Synthesize a Clifford circuit implementing a target Clifford isometry specified by its logical operators and its
    stabilizers.
    The tableau is specified by a list of operators (as strings). The first (resp. second) half of the list corresponding
    to the images of `Zi` (resp. `Xi`) operators.

    Args:
        logicals (list): The logical operators.
        stabilisers (list): The stabilizers.
        metric (Metric): The metric to minimize.
        syndrome_iter (optional, int): The number of syndrome decoding iteration used in the ISD solver in the metric=`Metric.COUNT` case only. Default is 1.

    Returns:
        list: The synthesized circuit.
    """
    syndrome_iter = syndrome_iter or 1
    return rust_isometry(logicals, stabilisers, metric, syndrome_iter)


def clifford_synthesis(logicals, metric, syndrome_iter=1):
    """
    Synthesize a Clifford circuit implementing a target Clifford operator specified by its
    tableau.
    The tableau is specified by a list of operators (as strings). The first (resp. second) half of the list corresponding
    to the images of `Zi` (resp. `Xi`) operators.

    Args:
        logicals (list): The tableau.
        metric (Metric): The metric to minimize.
        syndrome_iter (optional, int): The number of syndrome decoding iteration used in the ISD solver in the metric=`Metric.COUNT` case only. Default is 1.

    Returns:
        list: The synthesized circuit.
    """
    syndrome_iter = syndrome_iter or 1
    return rust_isometry(logicals, [], metric, syndrome_iter)
