/// This module contains the implementation of a struct representing graph states
use super::pauli_like::PauliLike;
use rand::Rng;

#[derive(Clone, Debug)]
pub struct GraphState {
    pub adj: Vec<Vec<bool>>,
    pub n: usize,
}

impl GraphState {
    pub fn new(n: usize) -> Self {
        Self {
            adj: vec![vec![false; n]; n],
            n,
        }
    }
    pub fn random(n: usize) -> Self {
        let mut rng = rand::thread_rng();
        let mut gs = Self::new(n);
        for i in 0..n {
            for j in i..n {
                let entry = rng.gen::<bool>();
                gs.adj[i][j] = entry;
                gs.adj[j][i] = entry;
            }
        }
        return gs;
    }
    pub fn from_adj(adj: Vec<Vec<bool>>) -> Self {
        let n = adj.len();
        for row in adj.iter() {
            assert_eq!(row.len(), n, "Matrix is not square");
        }
        for i in 0..n {
            for j in 0..n {
                assert_eq!(adj[i][j], adj[j][i], "Matrix is not symmetric");
            }
        }
        Self { adj, n }
    }
    pub fn count_ones(&self) -> usize {
        self.adj
            .iter()
            .map(|row| row.iter().filter(|a| **a).count())
            .sum::<usize>()
    }
}
impl PauliLike for GraphState {
    fn h(&mut self, _: usize) {
        panic!("You are not supposed to apply H to a graph state!");
    }

    fn s(&mut self, i: usize) {
        self.adj[i][i] ^= true;
    }

    fn sd(&mut self, i: usize) {
        self.adj[i][i] ^= true;
    }

    fn sqrt_x(&mut self, _: usize) {
        panic!("You are not supposed to apply SQRT_X to a graph state!");
    }

    fn sqrt_xd(&mut self, _: usize) {
        panic!("You are not supposed to apply SQRT_XD to a graph state!");
    }

    fn cnot(&mut self, i: usize, j: usize) {
        for k in 0..self.n {
            self.adj[i][k] ^= self.adj[j][k];
        }
        for k in 0..self.n {
            self.adj[k][i] ^= self.adj[k][j];
        }
    }

    fn cz(&mut self, i: usize, j: usize) {
        self.adj[i][j] ^= true;
        self.adj[j][i] ^= true;
    }
}
