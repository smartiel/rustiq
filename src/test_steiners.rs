use petgraph::algo::{bellman_ford, min_spanning_tree};
use petgraph::data::Build;
use petgraph::data::FromElements;
use petgraph::graph::Node;
use petgraph::prelude::*;
use petgraph::visit::IntoNeighbors;
use petgraph::Graph;
use std::collections::HashMap;
use std::hash::Hash;

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
    graph: &Graph<(), f64, Undirected>,
    conv_paths: &Vec<Vec<(usize, Vec<usize>)>>,
    terminals: &Vec<usize>,
) -> (
    Graph<(), i32, Undirected>,
    HashMap<usize, NodeIndex>,
    HashMap<NodeIndex, usize>,
) {
    let closure = build_closure_mst(&conv_paths, &vec![1, 2, 4]);
    return build_steiner_from_mst(&closure, &vec![1, 2, 4], &conv_paths);
}

fn main() {
    let mut graph = UnGraph::new_undirected();
    let nodes: Vec<_> = (0..6).map(|_| graph.add_node(())).collect();
    graph.add_edge(NodeIndex::new(0), NodeIndex::new(1), 1.);
    graph.add_edge(NodeIndex::new(1), NodeIndex::new(2), 1.);
    graph.add_edge(NodeIndex::new(2), NodeIndex::new(3), 1.);
    graph.add_edge(NodeIndex::new(0), NodeIndex::new(3), 1.);
    graph.add_edge(NodeIndex::new(0), NodeIndex::new(4), 1.);
    graph.add_edge(NodeIndex::new(5), NodeIndex::new(4), 1.);
    graph.remove_node(NodeIndex::new(2));
    println!("{:?}", graph);
}
