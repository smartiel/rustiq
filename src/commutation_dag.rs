/// This module implements a data structure encoding the commutation DAG of a sequence of Pauli operators.
use petgraph::graph::DiGraph;

pub type PauliDag = DiGraph<, ()>