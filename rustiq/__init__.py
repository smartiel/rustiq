import networkx as nx
from .rustiq import Metric, graph_state_synthesis as rust_gs, greedy_pauli_network

__all__ = ["Metric", "graph_state_synthesis", "greedy_pauli_network"]


def graph_state_synthesis(graph, metric, syndrome_iter=None):
    """
    Synthesize a circuit preparing a given graph state specified by a networkx graph.

    """
    syndrome_iter = syndrome_iter or 1
    adj = nx.to_numpy_array(graph) < 0.5
    return rust_gs(adj.tolist(), metric, syndrome_iter)
