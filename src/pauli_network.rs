use rand::Rng;
use rustiq::structures::metric::Metric;
use rustiq::synthesis::pauli_network::greedy_pauli_network;
use std::time::{Duration, SystemTime};
const PAULIS: [char; 4] = ['I', 'X', 'Y', 'Z'];

fn random_instance(n: usize, m: usize) -> Vec<String> {
    let mut rng = rand::thread_rng();
    let mut output = Vec::new();
    for _ in 0..m {
        let mut s = String::new();
        for _ in 0..n {
            s.push(PAULIS[rng.gen::<usize>() % 4]);
        }
        output.push(s);
    }
    return output;
}

fn main() {
    let instance = random_instance(30, 200);
    let now = SystemTime::now();
    let mut pset = rustiq::structures::PauliSet::from_slice(&instance);
    match now.elapsed() {
        Ok(elapsed) => {
            // it prints '2'
            println!("Time for PauliSet construction:{}", elapsed.as_secs_f64());
        }
        Err(e) => {
            // an error occurred!
            println!("Error: {e:?}");
        }
    }
    let now = SystemTime::now();
    let result = greedy_pauli_network(&mut pset, &Metric::DEPTH, false, 0, false, false);
    match now.elapsed() {
        Ok(elapsed) => {
            // it prints '2'
            println!("Time for synthesis:{}", elapsed.as_secs_f64());
        }
        Err(e) => {
            // an error occurred!
            println!("Error: {e:?}");
        }
    }
    println!("{:?}", result.cnot_count());
}
