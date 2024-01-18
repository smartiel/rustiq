use crate::structures::{CliffordCircuit, CliffordGate};
pub type Matrix = Vec<Vec<bool>>;

pub fn xor_vec(a: &mut Vec<bool>, b: &Vec<bool>) {
    for i in 0..a.len() {
        a[i] ^= b[i];
    }
}

pub fn rowop(table: &mut Matrix, i: usize, j: usize) {
    if table.len() != 0 {
        for k in 0..table.first().unwrap().len() {
            table[j][k] ^= table[i][k];
        }
    }
}

pub fn colop(table: &mut Matrix, i: usize, j: usize) {
    for k in 0..table.len() {
        table[k][j] ^= table[k][i];
    }
}

pub fn row_echelon(table: &mut Matrix, k: usize) {
    let mut rank = 0;
    for i in 0..table.first().unwrap().len() {
        let mut pivot = None;
        for j in rank..k {
            if table[j][i] {
                pivot = Some(j);
                break;
            }
        }
        if let Some(pivot) = pivot {
            if pivot != rank {
                rowop(table, pivot, rank);
            }
            for j in 0..table.len() {
                if table[j][i] && j != rank {
                    rowop(table, rank, j);
                }
            }
            rank += 1;
        }
    }
}

pub fn diagonalize(table: &mut Matrix, friend: &mut Matrix, rank: usize) {
    let n = table.first().unwrap().len();
    for i in 0..rank {
        let mut pivot = None;
        for j in i..n {
            if table[i][j] {
                pivot = Some(j);
                break;
            }
        }
        if let Some(pivot) = pivot {
            if pivot != i {
                colop(table, pivot, i);
                colop(friend, pivot, i);
            }
            for j in 0..n {
                if table[i][j] && j != i {
                    colop(table, i, j);
                    colop(friend, i, j);
                }
            }
        } else {
            panic!("This is not gooood!");
        }
    }
}

pub fn f2_rank(table: &Matrix) -> usize {
    let mut rank = 0;
    let mut table = table.clone();
    let nrows = table.len();
    let ncols = table.first().unwrap().len();
    for i in 0..ncols {
        let mut pivot = None;
        for j in rank..nrows {
            if table[j][i] {
                pivot = Some(j);
                break;
            }
        }
        if let Some(pivot) = pivot {
            if pivot != rank {
                rowop(&mut table, pivot, rank);
            }
            for j in rank + 1..nrows {
                if table[j][i] {
                    rowop(&mut table, rank, j);
                }
            }
            rank += 1;
        }
    }

    rank
}

pub fn inverse_f2(table: &Matrix) -> Matrix {
    let n = table.len();
    let mut table = table.clone();
    let mut friend = vec![vec![false; n]; n];
    for i in 0..n {
        friend[i][i] = true;
    }
    for i in 0..n {
        let mut pivot = None;
        for j in i..n {
            if table[j][i] {
                pivot = Some(j);
                break;
            }
        }
        if let Some(pivot) = pivot {
            table.swap(i, pivot);
            friend.swap(i, pivot);
            for j in 0..n {
                if j != i && table[j][i] {
                    rowop(&mut table, i, j);
                    rowop(&mut friend, i, j)
                }
            }
        }
    }
    return friend;
}

pub fn mult_f2(a: &Matrix, b: &Matrix) -> Matrix {
    let (k, l, m) = (a.len(), a.first().unwrap().len(), b.first().unwrap().len());
    assert_eq!(l, b.len());
    let mut result = vec![vec![false; m]; k];
    for i in 0..k {
        for j in 0..m {
            for y in 0..l {
                result[i][j] ^= a[i][y] & b[y][j];
            }
        }
    }
    result
}

pub fn transpose(table: &Matrix) -> Matrix {
    let n = table.len();
    let k = table.first().unwrap().len();
    let mut result = vec![vec![false; n]; k];
    for i in 0..k {
        for j in 0..n {
            result[i][j] = table[j][i];
        }
    }
    result
}

pub fn plu_facto(table: &Matrix) -> (Matrix, Matrix, Matrix) {
    let n = table.len();
    let k = table.first().unwrap().len();
    let mut u = table.clone();
    let mut l = vec![vec![false; k]; n];
    let mut p = vec![vec![false; n]; n];
    for i in 0..n {
        p[i][i] = true;
    }
    for i in 0..k {
        let mut pivot = None;
        for j in i..n {
            if u[j][i] {
                pivot = Some(j);
                break;
            }
        }

        if let Some(pivot) = pivot {
            if i != pivot {
                u.swap(i, pivot);
                p.swap(i, pivot);
                l.swap(i, pivot);
            }
            l[i][i] = true;
            for j in i + 1..n {
                if u[j][i] {
                    l[j][i] = true;
                    rowop(&mut u, i, j);
                }
            }
        }
    }

    return (p, l, u.drain(..k).collect());
}

pub fn lu_facto(table: &Matrix) -> (Matrix, Matrix, Matrix, CliffordCircuit) {
    let mut table = table.clone();
    let (c, ops) = non_zero_leading_principal_minors(&table);
    let mut output = CliffordCircuit::new(table.len());
    for (a, b) in ops.iter() {
        rowop(&mut table, *a, *b);
        output.gates.push(CliffordGate::CNOT(*b, *a));
    }
    let (p, l, u) = plu_facto(&table);
    for i in 0..table.len() {
        assert!(p[i][i]);
    }
    return (l, u, c, output);
}

fn f2_rank_square(matrix: &Matrix) -> usize {
    let matrix: Matrix = matrix
        .iter()
        .map(|v| v.clone().into_iter().take(matrix.len()).collect())
        .collect();
    return f2_rank(&matrix);
}
pub fn print_matrix(matrix: &Matrix) {
    for row in matrix.iter() {
        for elem in row.iter() {
            if *elem {
                print!("1");
            } else {
                print!("0");
            }
        }
        println!("");
    }
}
/// Finds a sequence of row operations that makes sure that all the leading principal
/// minors are 1
/// In other words, makes sure that matrix[:i, :i] is invertible
pub fn non_zero_leading_principal_minors(matrix: &Matrix) -> (Matrix, Vec<(usize, usize)>) {
    let mut piece: Matrix = Vec::new();
    let mut moves = Vec::new();
    for i in 0..f2_rank(&matrix) {
        piece.push(matrix[i].clone());
        let mut current_rk = f2_rank_square(&piece);
        while current_rk == i {
            for k in i + 1..matrix.len() {
                xor_vec(&mut piece[i], &matrix[k]);
                let new_rank = f2_rank_square(&piece);
                if new_rank == i + 1 {
                    moves.push((k, i));
                    current_rk = new_rank;
                    break;
                } else {
                    xor_vec(&mut piece[i], &matrix[k]);
                }
            }
            if current_rk == i + 1 {
                break;
            }
        }
    }
    let mut c = vec![vec![false; matrix.len()]; matrix.len()];
    for i in 0..c.len() {
        c[i][i] = true;
    }
    for (a, b) in moves.iter() {
        rowop(&mut c, *a, *b);
    }
    let m = mult_f2(&c, &matrix);
    for i in 0..piece.len() {
        assert_eq!(piece[i], m[i]);
    }
    (c, moves)
}

pub fn count_ones(matrix: &Matrix) -> usize {
    matrix
        .iter()
        .map(|row| row.iter().filter(|x| **x).count())
        .sum()
}

pub fn count_ones_except_diag(matrix: &Matrix) -> usize {
    let mut count = 0;
    for i in 0..matrix.len() {
        for j in 0..matrix[i].len() {
            if i != j && matrix[i][j] {
                count += 1;
            }
        }
    }
    return count;
}
#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    fn random_skinny(n: usize, m: usize) -> Matrix {
        assert!(m <= n);
        let mut rng = rand::thread_rng();
        let mut matrix = vec![vec![false; m]; n];
        for i in 0..m {
            matrix[i][i] = true;
        }
        for _ in 0..n * n {
            let i = rng.gen_range(0..n);
            let j = rng.gen_range(0..n);
            if i != j {
                rowop(&mut matrix, i, j);
            }
        }
        matrix
    }
    #[test]
    fn test_plu_facto() {
        for _ in 0..10 {
            let n = 40;
            let m = 20;
            let matrix = random_skinny(n, m);
            let (p, l, u) = plu_facto(&matrix);

            // Check if P is a permutation matrix
            assert_eq!(
                p.iter()
                    .map(|r| r.iter().filter(|&&v| v).count())
                    .collect::<Vec<_>>(),
                vec![1; matrix.len()]
            );
            let mut perm = p
                .iter()
                .map(|r| r.iter().position(|&v| v).unwrap())
                .collect::<Vec<_>>();
            perm.sort();
            assert_eq!(perm, (0..n).collect::<Vec<_>>());

            for i in 0..l.len() {
                for j in i + 1..l[i].len() {
                    assert!(!l[i][j], "L should be lower triangular");
                }
            }

            for i in 0..u.len() {
                for j in 0..std::cmp::min(u[i].len(), i) {
                    assert!(!u[i][j], "U should be upper triangular");
                }
            }

            // Check if PA = LU
            let pa = mult_f2(&p, &matrix);
            let lu = mult_f2(&l, &u);
            assert_eq!(pa, lu, "PA should be equal to LU");
        }
    }

    #[test]
    fn test_lu_facto() {
        for _ in 0..10 {
            let n = 40;
            let m = 22;
            let matrix = random_skinny(n, m);
            let (l, u, c, _) = lu_facto(&matrix);

            for i in 0..l.len() {
                for j in i + 1..l[i].len() {
                    assert!(!l[i][j], "L should be lower triangular");
                }
            }

            for i in 0..u.len() {
                for j in 0..std::cmp::min(u[i].len(), i) {
                    assert!(!u[i][j], "U should be upper triangular");
                }
            }
            let lu = mult_f2(&l, &u);
            let clu = mult_f2(&inverse_f2(&c), &lu);
            assert_eq!(clu, matrix, "C^-1 LU should be equal to A");
        }
    }

    fn random_invertible(n: usize) -> Matrix {
        let mut rng = rand::thread_rng();
        let mut matrix = vec![vec![false; n]; n];
        for i in 0..n {
            matrix[i][i] = true;
        }
        for _ in 0..n * n {
            let i = rng.gen_range(0..n);
            let j = rng.gen_range(0..n);
            if i != j {
                rowop(&mut matrix, i, j);
            }
        }
        matrix
    }

    #[test]
    fn test_leading_full_rank() {
        for _ in 0..40 {
            // Test data
            let n = 40;
            let mut test_matrix = random_invertible(n);
            let (_, moves) = non_zero_leading_principal_minors(&test_matrix);
            for (a, b) in moves {
                rowop(&mut test_matrix, a, b);
            }
            for i in 1..n {
                let piece: Matrix = test_matrix.clone().into_iter().take(i).collect();
                assert_eq!(f2_rank_square(&piece), i);
            }
        }
    }
}
