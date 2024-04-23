use petgraph::algo::{bellman_ford, min_spanning_tree};
use petgraph::data::FromElements;
use petgraph::prelude::*;
use petgraph::Graph;
use std::collections::HashMap;

fn convert_shortest_paths(
    all_paths: &Vec<bellman_ford::Paths<NodeIndex, f64>>,
    nodes: &Vec<NodeIndex>,
) -> Vec<Vec<(usize, Vec<usize>)>> {
    let mut my_paths: Vec<Vec<(usize, Vec<usize>)>> = Vec::new();
    for n1 in nodes.iter() {
        my_paths.push(Vec::new());
        for n2 in nodes.iter() {
            let path = get_shortest_path(all_paths, n1, n2);
            my_paths.last_mut().unwrap().push((path.len() - 1, path));
        }
    }
    my_paths
}

fn get_shortest_path(
    all_paths: &Vec<bellman_ford::Paths<NodeIndex, f64>>,
    n1: &NodeIndex,
    n2: &NodeIndex,
) -> Vec<usize> {
    let path = &all_paths[n1.index()];
    let mut spath: Vec<usize> = Vec::new();
    spath.push(n2.index());
    let mut v = path.predecessors[n2.index()];
    while let Some(ni) = v {
        spath.push(ni.index());
        v = path.predecessors[ni.index()];
    }
    return spath;
}

fn build_closure_mst(
    all_paths: &Vec<Vec<(usize, Vec<usize>)>>,
    terminals: &Vec<usize>,
) -> Graph<(), f64, Undirected> {
    let mut g = Graph::new_undirected();
    let nodes: Vec<_> = (0..terminals.len()).map(|_| g.add_node(())).collect();

    for (i1, n1) in nodes.iter().enumerate() {
        for (i2, n2) in nodes.iter().enumerate() {
            if i1 != i2 {
                g.add_edge(*n1, *n2, all_paths[terminals[i1]][terminals[i2]].0 as f64);
            }
        }
    }
    let mst: Graph<(), f64, Undirected> = Graph::from_elements(min_spanning_tree(&g));
    mst
}

#[derive(Debug, Clone)]
pub struct SteinerTree {
    pub graph: UnGraph<(), i32, u32>,
    pub mapping: HashMap<usize, NodeIndex>,
    pub inverse_mapping: HashMap<NodeIndex, usize>,
    pub terminals: Vec<usize>,
}

impl SteinerTree {
    pub fn len(&self) -> usize {
        self.graph.node_count()
    }

    pub fn cnot_cost(&self) -> usize {
        return 2 * self.len() - self.terminals.len() - 1;
    }

    pub fn contains(&self, node: usize) -> bool {
        self.mapping.get(&node).is_some()
    }
    pub fn is_leaf(&self, node: usize) -> bool {
        let node_index = self.mapping.get(&node).unwrap();
        self.graph.neighbors_undirected(*node_index).count() == 1
    }
    pub fn neighbors(&self, node: usize) -> Vec<usize> {
        let mut neighs = Vec::new();
        for nei in self.graph.neighbors_undirected(self.mapping[&node]) {
            neighs.push(self.inverse_mapping[&nei]);
        }
        neighs
    }
    pub fn nodes(&self) -> Vec<usize> {
        self.mapping.keys().map(|k| *k).collect()
    }
    pub fn add_node(&mut self, node: usize) {
        let node_index = self.graph.add_node(());
        self.mapping.insert(node, node_index);
        self.inverse_mapping.insert(node_index, node);
    }

    pub fn add_edge(&mut self, n1: usize, n2: usize) {
        self.graph.add_edge(self.mapping[&n1], self.mapping[&n2], 1);
    }

    pub fn remove_node(&mut self, node: usize) {
        if !self.contains(node) {
            panic!("Node {} is not in the tree", node);
        }
        let last_index: NodeIndex<u32> = NodeIndex::new(self.graph.node_count() - 1);
        let node_index = self.mapping.get(&node).unwrap().clone();
        let i = self.inverse_mapping.get(&last_index).unwrap().clone();
        if i == node {
            self.graph.remove_node(node_index);
            self.mapping.remove(&i);
            self.inverse_mapping.remove(&node_index);
            return;
        }
        self.inverse_mapping.remove(&last_index);
        self.mapping.remove(&node);
        self.mapping.insert(i, node_index);
        self.inverse_mapping.insert(node_index, i);
        self.graph.remove_node(node_index);
    }

    pub fn update_tree(&mut self, new_support: &Vec<usize>, qbits: &[usize; 2]) {
        self.terminals = new_support.clone();
        if self.contains(qbits[0]) && self.contains(qbits[1]) {
            if new_support.contains(&qbits[0]) && new_support.contains(&qbits[1]) {
                return;
            }
            if !new_support.contains(&qbits[0]) && self.is_leaf(qbits[0]) {
                self.remove_node(qbits[0]);
            }
            if !new_support.contains(&qbits[1]) && self.is_leaf(qbits[1]) {
                self.remove_node(qbits[1]);
            }
            return;
        }

        if !self.contains(qbits[0]) && new_support.contains(&qbits[0]) {
            self.add_node(qbits[0]);
            self.add_edge(qbits[0], qbits[1]);
            return;
        }

        if !self.contains(qbits[1]) && new_support.contains(&qbits[1]) {
            self.add_node(qbits[1]);
            self.add_edge(qbits[0], qbits[1]);
            return;
        }
        return;
    }
    pub fn prune_non_terminal_leaves(&mut self) {
        loop {
            let mut popped = false;
            for node in self.graph.node_indices() {
                if !self.terminals.contains(&self.inverse_mapping[&node])
                    && self.graph.neighbors(node).count() == 1
                {
                    self.remove_node(self.inverse_mapping[&node]);
                    popped = true;
                    break;
                }
            }
            if !popped {
                break;
            }
        }
    }
}

fn build_steiner_from_mst(
    closure_mst: &Graph<(), f64, Undirected>,
    terminals: &Vec<usize>,
    all_paths: &Vec<Vec<(usize, Vec<usize>)>>,
) -> SteinerTree {
    // Building the extended mst
    let mut g = Graph::new_undirected();
    let mut nodes = HashMap::new();
    let mut inverse_map = HashMap::new();
    for i in terminals.iter() {
        let nodei = g.add_node(());
        nodes.insert(*i, nodei);
        inverse_map.insert(nodei, *i);
    }
    for edge in closure_mst.raw_edges() {
        let t1 = terminals[edge.source().index()];
        let t2 = terminals[edge.target().index()];
        let (_, path) = &all_paths[t1][t2];
        for i in path.iter() {
            if !nodes.contains_key(i) {
                let nodei = g.add_node(());
                nodes.insert(*i, nodei);
                inverse_map.insert(nodei, *i);
            }
        }
        let mut left = path[0];
        for next_v in path.iter().skip(1) {
            g.add_edge(nodes[&left], nodes[next_v], 1);
            left = *next_v;
        }
    }
    let mst: Graph<(), i32, Undirected> = Graph::from_elements(min_spanning_tree(&g));
    let mut as_st = SteinerTree {
        graph: mst,
        mapping: nodes,
        inverse_mapping: inverse_map,
        terminals: terminals.clone(),
    };
    as_st.prune_non_terminal_leaves();
    return as_st;
}

fn get_steiner_tree(
    conv_paths: &Vec<Vec<(usize, Vec<usize>)>>,
    terminals: &Vec<usize>,
) -> SteinerTree {
    let closure = build_closure_mst(&conv_paths, &terminals);
    return build_steiner_from_mst(&closure, terminals, &conv_paths);
}
#[derive(Debug, Clone)]
pub struct HardwareGraph {
    pub graph: UnGraph<(), f64, u32>,
    shortest_paths: Vec<Vec<(usize, Vec<usize>)>>,
}

impl HardwareGraph {
    pub fn from_couplings(couplings: &Vec<(usize, usize)>) -> Self {
        let mut graph = UnGraph::new_undirected();
        let n = couplings.iter().map(|c| c.0.max(c.1)).max().unwrap() + 1;
        let nodes: Vec<_> = (0..n).map(|_| graph.add_node(())).collect();
        graph.extend_with_edges(couplings.iter().map(|c| (nodes[c.0], nodes[c.1], 1.)));
        let all_paths: Vec<bellman_ford::Paths<NodeIndex, f64>> = nodes
            .iter()
            .map(|node_index| bellman_ford(&graph, *node_index).unwrap())
            .collect();
        let shortest_paths = convert_shortest_paths(&all_paths, &nodes);
        Self {
            graph,
            shortest_paths,
        }
    }
    pub fn len(&self) -> usize {
        return self.graph.node_count();
    }
    pub fn get_steiner_tree(&self, terminals: &Vec<usize>) -> SteinerTree {
        get_steiner_tree(&self.shortest_paths, terminals)
    }
}
