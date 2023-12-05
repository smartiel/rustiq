use rustiq::structures::graph_state::GraphState;
use rustiq::structures::metric::Metric;
use rustiq::synthesis::clifford::graph_state::synthesize_graph_state;

fn main() {
    // for n in &[10, 30, 50] {
    //     let graph = GraphState::random(*n);
    //     for niter in &[1, 10, 50, 100] {
    //         let circuit_count = synthesize_graph_state(&graph, &Metric::COUNT, *niter);
    //         println!("{}\t{}\t{}", n, niter, circuit_count.entangling_count());
    //     }
    // }
    for n in &[10, 30, 50, 100] {
        let graph = GraphState::random(*n);
        let circiut_depth = synthesize_graph_state(&graph, &Metric::DEPTH, 0);
        println!("{}\t{}", n, circiut_depth.entangling_depth());
    }
}
