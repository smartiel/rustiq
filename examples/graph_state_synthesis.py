from rustiq import graph_state_synthesis, Metric
from rustiq.utils import entangling_count, entangling_depth
import networkx as nx


graph = nx.generators.erdos_renyi_graph(40, 0.5)

circuit = graph_state_synthesis(graph, Metric.COUNT)
print(entangling_count(circuit), entangling_depth(circuit))

circuit = graph_state_synthesis(graph, Metric.DEPTH)
print(entangling_count(circuit), entangling_depth(circuit))


# Increasing the number of iterations in the syndrome decoding (only for count):
for niter in [1, 10, 20, 50, 100]:
    circuit = graph_state_synthesis(graph, Metric.COUNT, niter)
    print(niter, entangling_count(circuit))
