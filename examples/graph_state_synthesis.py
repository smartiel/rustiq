from rustiq import graph_state_synthesis, Metric
import networkx as nx


def ent_depth(circ):
    depths = {}
    for gate, qbits in circ:
        if gate in ["CNOT", "CZ"]:
            gate_depth = max(depths.get(qbits[0], 0), depths.get(qbits[1], 0)) + 1
            for qbit in qbits:
                depths[qbit] = gate_depth
    return max(depths.values())


def ent_count(circ):
    return sum(gate in ["CNOT", "CZ"] for gate, _ in circ)


graph = nx.generators.erdos_renyi_graph(40, 0.5)

circuit = graph_state_synthesis(graph, Metric.COUNT)
print(ent_count(circuit), ent_depth(circuit))

circuit = graph_state_synthesis(graph, Metric.DEPTH)
print(ent_count(circuit), ent_depth(circuit))


# Increasing the number of iterations in the syndrome decoding (only for count):
for niter in [1, 10, 20, 50, 100]:
    circuit = graph_state_synthesis(graph, Metric.COUNT, niter)
    print(niter, ent_count(circuit))
