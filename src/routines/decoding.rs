use rand::seq::SliceRandom;
use rand::thread_rng;

pub fn syndrome_decoding(parities: &Vec<Vec<bool>>, input_target: &Vec<bool>) -> Vec<bool> {
    let mut target = input_target.clone();
    let mut solution = vec![false; parities.len()];
    let mut hweight = target.iter().filter(|a| **a).count();
    loop {
        let mut best_reduction = 0;
        let mut best_index: i32 = -1;
        for (index, parity) in parities.iter().enumerate() {
            let new_hweight = parity
                .iter()
                .zip(target.iter())
                .map(|(a, b)| a ^ b)
                .filter(|a| *a)
                .count();
            if (hweight as i32 - new_hweight as i32) > best_reduction {
                best_reduction = hweight as i32 - new_hweight as i32;
                best_index = index as i32;
            }
        }
        if best_index == -1 {
            break;
        }
        for (a, b) in target.iter_mut().zip(parities[best_index as usize].iter()) {
            *a ^= b;
        }
        solution[best_index as usize] ^= true;
        hweight = target.iter().filter(|a| **a).count();
    }
    let mut true_target = vec![false; target.len()];
    for (i, b) in solution.iter().enumerate() {
        if *b {
            for (x, y) in parities[i].iter().zip(true_target.iter_mut()) {
                *y ^= x;
            }
        }
    }
    assert_eq!(true_target, *input_target);
    return solution;
}

fn colop(parities: &mut Vec<Vec<bool>>, i: usize, j: usize) {
    for k in 0..parities.len() {
        parities[k][j] ^= parities[k][i];
    }
}

fn shuffle_parities(parities: &mut Vec<Vec<bool>>, target: &mut Vec<bool>) -> Vec<usize> {
    let n = parities.first().unwrap().len();
    let mut row_permutation: Vec<usize> = (0..parities.len()).collect();
    row_permutation.shuffle(&mut thread_rng());
    let mut new_parities = Vec::new();
    for j in row_permutation.iter() {
        new_parities.push(parities[*j].clone());
    }
    let mut rank = 0;
    for i in 0..parities.len() {
        let mut pivot = None;
        for j in rank..n {
            if new_parities[i][j] {
                pivot = Some(j);
                break;
            }
        }
        if let Some(pivot) = pivot {
            if pivot != rank {
                colop(&mut new_parities, pivot, rank);
                target[rank] ^= target[pivot];
            }
            for j in 0..n {
                if new_parities[i][j] && j != rank {
                    colop(&mut new_parities, rank, j);
                    target[j] ^= target[rank];
                }
            }
            rank += 1;
            if rank == n {
                break;
            }
        }
    }
    *parities = new_parities;

    row_permutation
}

fn fix_permutation(solution: &Vec<bool>, permutation: &Vec<usize>) -> Vec<bool> {
    let mut new_solution = vec![false; solution.len()];
    for (i, j) in permutation.iter().enumerate() {
        new_solution[*j] = solution[i];
    }
    new_solution
}

pub fn information_set_decoding(
    input_parities: &Vec<Vec<bool>>,
    input_target: &Vec<bool>,
    ntries: usize,
) -> Vec<bool> {
    let mut best_solution: Option<Vec<bool>> = None;
    let mut best_cost = None;
    for _ in 0..ntries {
        let mut parities = input_parities.clone();
        let mut target = input_target.clone();
        let permutation = shuffle_parities(&mut parities, &mut target);
        let solution = fix_permutation(&syndrome_decoding(&parities, &mut target), &permutation);
        let cost = solution.iter().filter(|a| **a).count();
        if let Some(best_cost) = best_cost {
            if best_cost < cost {
                continue;
            }
        }
        best_cost = Some(cost);
        best_solution = Some(solution);
    }
    let solution = best_solution.unwrap();
    let mut true_target = vec![false; input_target.len()];
    for (i, b) in solution.iter().enumerate() {
        if *b {
            for (x, y) in input_parities[i].iter().zip(true_target.iter_mut()) {
                *y ^= x;
            }
        }
    }
    assert_eq!(true_target, *input_target);
    solution
}

#[cfg(test)]
mod decoding_tests {
    use super::*;
    #[test]
    fn test_raw_decoding() {
        let mut parities = vec![
            vec![true, false, false],
            vec![false, true, false],
            vec![false, false, true],
            vec![true, false, true],
        ];
        let mut target = vec![true, true, true];
        let solution = syndrome_decoding(&mut parities, &mut target);
        println!("{:?}", solution);
        assert_eq!(solution, vec![false, true, false, true]);
    }
    #[test]
    fn test_isd() {
        let mut parities = vec![
            vec![true, false, false],
            vec![false, true, false],
            vec![false, false, true],
            vec![true, false, true],
        ];
        let mut target = vec![true, true, true];
        let solution = information_set_decoding(&mut parities, &mut target, 100);
        println!("{:?}", solution);
    }
}
