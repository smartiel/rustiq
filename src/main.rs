extern crate clap;
use clap::{command, Arg};
use rustiq::greedy_order_preserving::pauli_network_synthesis_no_permutation;
use rustiq::greedy_pauli_network::{pauli_network_synthesis, Metric};
use rustiq::pauli_set::PauliSet;
use std::env;
use std::fs;
use std::time::Instant;

fn load_pauli_set_from_file(file_path: &str) -> Vec<String> {
    return fs::read_to_string(file_path)
        .expect("File does not exist")
        .split('\n')
        .filter(|s| s.len() > 0)
        .map(|s| s.to_owned())
        .collect::<Vec<String>>();
}

fn main() {
    let mut pset = PauliSet::new(2);
    for _ in 0..(10 * 1024) {
        pset.insert_vec_bool(&vec![false; 4], false);
    }
    let start = Instant::now();
    for _ in 0..1000 {
        pset.cnot(0, 1)
    }
    let duration = start.elapsed();
    println!("{}", duration.as_secs_f64(),);

    // let args = command!()
    //     .version("0.0.1")
    //     .about("Command line interface for rustiq")
    //     .arg(
    //         Arg::new("FILE")
    //             .help("A file containing the set of Pauli operators")
    //             .required(true),
    //     )
    //     .arg(
    //         Arg::new("metric")
    //             .long("metric")
    //             .help("The quantity to minimize [depth/count]")
    //             .default_value("count"),
    //     )
    //     .arg(
    //         Arg::new("onlyinfo")
    //             .long("onlyinfo")
    //             .action(clap::ArgAction::SetTrue)
    //             .help("If set, only display a tuple (synthesis time in seconds, CNOT count, CNOT depth)"),
    //     ).arg(Arg::new("keeporder")
    //     .long("keeporder")
    //     .action(clap::ArgAction::SetTrue)
    //     .help("If set, preserves the operators order (up to allowed commutations)"),)
    //     .get_matches();
    // let file_path = args
    //     .get_one::<String>("FILE")
    //     .expect("This should never fail :)");
    // let pauli_operators = load_pauli_set_from_file(file_path);
    // let metric = Metric::from_string(
    //     args.get_one::<String>("metric")
    //         .expect("This should never fail :)"),
    // )
    // .expect("Unknown metric");
    // let start = Instant::now();
    // let result = if args.get_flag("keeporder") {
    //     pauli_network_synthesis_no_permutation(pauli_operators, metric)
    // } else {
    //     pauli_network_synthesis(pauli_operators, metric)
    // };
    // let duration = start.elapsed();
    // if args.get_flag("onlyinfo") {
    //     println!(
    //         "{},{},{}",
    //         duration.as_secs_f64(),
    //         result.cnot_count(),
    //         result.cnot_depth()
    //     );
    // } else {
    //     for gate in result.gates {
    //         println!("{:?}", gate);
    //     }
    // }
}
