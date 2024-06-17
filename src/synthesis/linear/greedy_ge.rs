use crate::routines::f2_linalg::{Matrix, lu_facto, row_op};



fn select_rows(matrix: &Matrix, offset: usize) -> (usize, usize) {
    let mut rows: Vec<usize> = offset..matrix.len().collect();
    for k in offset..matrix.len() {
        let mut set_0 = Vec::new();
        let mut set_1 = Vec::new();
        for row in rows {
            if matrix[row][k]{
                set_1.push(row);
            } else {
                set_0.push(row);
            }
        }
        if set_1.len() > 1 {
            rows = set_1;
        }else{
            rows = set_0;
        }
        if rows.len() == 2 {
            break;
        }
    }
    return (rows[0], rows[1]);
}


fn greedy_ge_l(matrix: &Matrix) -> Vec<(usize, usize)> {
    let mut matrix = matrix.clone();
    let mut moves = Vec::new();
    loop {
        let (a, b) = select_rows(&matrix, 0);
        moves.push((a, b));
        row_op(&mut matrix, a, b);
        let mut is_identity = true;
        for i in 0..matrix.len() {
            for j in 0..i {
                if matrix[i][j] {
                    is_identity = false;
                    break
                }
            }
            if !is_identity {
                break;
            }
        }
        if !is_identity {
            continue;
        }
        break;
    }
    moves
}

pub fn greedy_ge(matrix: &Matrix) -> Vec<(usize, usize)> {
    let (l, u, _, _) = lu_facto(matrix);
    let mut moves = greedy_ge_l(&l);
    moves.extend(greedy_ge_l(&u));
    moves
} 