import networkx as nx
from .rustiq import (
    Metric,
    graph_state_synthesis as rust_gs,
    greedy_pauli_network,
    codiagonalization as rust_codiag,
    codiagonalization_sswise as rust_sswise,
    isometry_synthesis as rust_isometry,
)

__all__ = [
    "Metric",
    "graph_state_synthesis",
    "greedy_pauli_network",
    "codiagonalization",
    "clifford_synthesis",
]


def graph_state_synthesis(graph, metric, syndrome_iter=None):
    """
    Synthesize a circuit preparing a given graph state specified by a networkx graph.

    """
    syndrome_iter = syndrome_iter or 1
    adj = nx.to_numpy_array(graph) < 0.5
    return rust_gs(adj.tolist(), metric, syndrome_iter)


def codiagonalization(paulis, metric, syndrome_iter=None):
    """
    Synthesize a circuit preparing a codiagonalizing a given set of pairwise commuting Pauli operators

    """
    syndrome_iter = syndrome_iter or 1
    return rust_codiag(paulis, metric, syndrome_iter)


def codiagonalization_sswise(paulis, k=2):
    """
    Synthesize a circuit preparing a codiagonalizing a given set of pairwise commuting Pauli operators

    """
    return rust_sswise(paulis, k)


def clifford_synthesis(logicals, metric=Metric.COUNT, syndrome_iter=None):
    """
    Synthesize a circuit preparing a codiagonalizing a given set of pairwise commuting Pauli operators

    """
    syndrome_iter = syndrome_iter or 1
    return rust_isometry(logicals, [], metric, syndrome_iter)
