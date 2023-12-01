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
}
impl PauliLike for GraphState {
    fn h(&mut self, i: usize) {
        panic!("You are not supposed to apply H to a graph state!");
    }

    fn s(&mut self, i: usize) {
        self.adj[i][i] ^= true;
    }

    fn sd(&mut self, i: usize) {
        self.adj[i][i] ^= true;
    }

    fn sqrt_x(&mut self, i: usize) {
        panic!("You are not supposed to apply SQRT_X to a graph state!");
    }

    fn sqrt_xd(&mut self, i: usize) {
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
