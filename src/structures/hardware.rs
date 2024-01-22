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
            g.add_edge(*n1, *n2, all_paths[terminals[i1]][terminals[i2]].0 as f64);
        }
    }
    let mst: Graph<(), f64, Undirected> = Graph::from_elements(min_spanning_tree(&g));
    mst
}

fn prune_non_terminal_leaves(
    closure_mst: &mut Graph<(), i32, Undirected>,
    terminals: &Vec<usize>,
    mapping: &HashMap<NodeIndex, usize>,
) {
    println!("{:?}", mapping);
    loop {
        let mut popped = false;
        for node in closure_mst.node_indices() {
            if !terminals.contains(&mapping[&node]) && closure_mst.neighbors(node).count() == 1 {
                closure_mst.remove_node(node);
                popped = true;
                break;
            }
        }
        if !popped {
            break;
        }
    }
}

fn build_steiner_from_mst(
    closure_mst: &Graph<(), f64, Undirected>,
    terminals: &Vec<usize>,
    all_paths: &Vec<Vec<(usize, Vec<usize>)>>,
) -> (
    Graph<(), i32, Undirected>,
    HashMap<usize, NodeIndex>,
    HashMap<NodeIndex, usize>,
) {
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
                println!("New node: {} -> {:?}", i, nodes[i]);
            }
        }
        let mut left = path[0];
        for next_v in path.iter().skip(1) {
            println!("Adding edge between {} and {}", left, next_v);
            g.add_edge(nodes[&left], nodes[next_v], 1);
            left = *next_v;
        }
    }
    let mut mst: Graph<(), i32, Undirected> = Graph::from_elements(min_spanning_tree(&g));
    prune_non_terminal_leaves(&mut mst, terminals, &inverse_map);
    return (mst, nodes, inverse_map);
}

fn get_steiner_tree(
    conv_paths: &Vec<Vec<(usize, Vec<usize>)>>,
    terminals: &Vec<usize>,
) -> (
    Graph<(), i32, Undirected>,
    HashMap<usize, NodeIndex>,
    HashMap<NodeIndex, usize>,
) {
    let closure = build_closure_mst(&conv_paths, &vec![1, 2, 4]);
    return build_steiner_from_mst(&closure, terminals, &conv_paths);
}

pub struct SteinerTree {
    pub graph: UnGraph<(), i32, u32>,
    mapping: HashMap<usize, NodeIndex>,
    inverse_mapping: HashMap<NodeIndex, usize>,
    terminals: Vec<usize>,
}

impl SteinerTree {
    pub fn len(&self) -> usize {
        self.graph.node_count()
    }
    pub fn cnot_cost(&self) -> usize {
        return 2 * self.len() - self.terminals.len() - 1;
    }
    pub fn update_tree(&mut self, new_support: &Vec<usize>, qbits: &[usize; 2]) {
        if self.graph.node_weight(self.mapping[&qbits[0]]).is_some()
            && self.graph.node_weight(self.mapping[&qbits[1]]).is_some()
        {
            if new_support.contains(&qbits[0]) && new_support.contains(&qbits[1]) {
                return;
            }
            if !new_support.contains(&qbits[0]) {
                if self.graph.neighbors(self.mapping[&qbits[0]]).count() == 1 {
                    self.graph.remove_node(self.mapping[&qbits[0]]);
                }
            }
            if !new_support.contains(&qbits[1]) {
                if self.graph.neighbors(self.mapping[&qbits[1]]).count() == 1 {
                    self.graph.remove_node(self.mapping[&qbits[1]]);
                }
            }
            return;
        }
        if self.graph.node_weight(self.mapping[&qbits[0]]).is_none()
            && new_support.contains(&qbits[0])
        {
            self.mapping.insert(qbits[0], self.graph.add_node(()));
            self.graph
                .add_edge(self.mapping[&qbits[0]], self.mapping[&qbits[1]], 1);
            return;
        }
        if self.graph.node_weight(self.mapping[&qbits[1]]).is_none()
            && new_support.contains(&qbits[1])
        {
            self.mapping.insert(qbits[1], self.graph.add_node(()));
            self.graph
                .add_edge(self.mapping[&qbits[1]], self.mapping[&qbits[0]], 1);
            return;
        }
        return;
    }
}

pub struct HardwareGraph {
    pub graph: UnGraph<(), f64, u32>,
    shortest_paths: Vec<Vec<(usize, Vec<usize>)>>,
}

impl HardwareGraph {
    pub fn from_couplings(couplings: &Vec<(usize, usize)>) -> Self {
        let mut graph = UnGraph::new_undirected();
        let n = couplings.iter().map(|c| c.0).max().unwrap() + 1;
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
        let (graph, mapping, inverse_mapping) = get_steiner_tree(&self.shortest_paths, terminals);
        SteinerTree {
            graph,
            mapping,
            inverse_mapping,
            terminals: terminals.clone(),
        }
    }
}
