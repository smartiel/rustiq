// The example graph:
//
//    +-- b ---- d ---- f
//   /    |      |
//  a     |      |
//   \    |      |
//    +-- c ---- e
//
// Maximum matching: { (a, b), (c, e), (d, f) }

#[cfg(test)]
mod graph_test_test {
    use petgraph::algo::maximum_matching;
    use petgraph::prelude::*;
    #[test]
    fn test_matching() {
        let mut graph: UnGraph<(), ()> = UnGraph::new_undirected();
        let a = graph.add_node(());
        let b = graph.add_node(());
        let c = graph.add_node(());
        let d = graph.add_node(());
        let e = graph.add_node(());
        let f = graph.add_node(());
        graph.extend_with_edges(&[(a, b), (a, c), (b, c), (b, d), (c, e), (d, e), (d, f)]);

        let matching = maximum_matching(&graph);
        assert!(matching.contains_edge(a, b));
        assert!(matching.contains_edge(c, e));
        assert_eq!(matching.mate(d), Some(f));
        assert_eq!(matching.mate(f), Some(d));
    }
}
