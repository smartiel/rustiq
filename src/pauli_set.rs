use std::cmp::max;

const WIDTH: usize = 64;
fn get_stride(index: usize) -> usize {
    return index / WIDTH;
}
fn get_offset(index: usize) -> usize {
    return index % WIDTH;
}
/// A set of Pauli operators (module global phase)
/// Conjugation by Clifford gates are vectorized
/// Warning:
///     Phases are not tracked ! (because I am lazy)
#[derive(Clone, Debug)]
pub struct PauliSet {
    pub n: usize,
    nstrides: usize,
    noperators: usize,
    start_offset: usize,
    /// The X and Z parts of the Pauli operators (in row major)
    /// The X part spans the first `n` rows and the Z part spans the last `n` rows
    data_array: Vec<Vec<u64>>,
}

impl PauliSet {
    /// Allocate an empty set of n-qubit Pauli operators
    pub fn new(n: usize) -> Self {
        Self {
            n,
            nstrides: 0,
            noperators: 0,
            start_offset: 0,
            data_array: vec![Vec::new(); 2 * n],
        }
    }
    // Construction from a list of operators
    pub fn from_slice(n: usize, data: &[String]) -> Self {
        let mut pset = Self::new(n);
        for piece in data {
            pset.insert(piece);
        }
        return pset;
    }
    /// Returns the number of operators stored in the set
    pub fn len(&self) -> usize {
        return self.noperators;
    }
    /// Inserts a new Pauli operator in the set and returns its index
    pub fn insert(&mut self, axis: &str) -> usize {
        let stride_index = get_stride(self.noperators + self.start_offset);
        let offset = get_offset(self.noperators + self.start_offset);
        if stride_index == self.nstrides {
            self.nstrides += 1;
            self.data_array.iter_mut().for_each(|row| row.push(0));
        }
        for (index, pauli) in axis.chars().enumerate() {
            match pauli {
                'Z' => {
                    self.data_array[index + self.n][stride_index] =
                        self.data_array[index + self.n][stride_index] | (1 << offset)
                }
                'X' => {
                    self.data_array[index][stride_index] =
                        self.data_array[index][stride_index] | (1 << offset)
                }
                'Y' => {
                    self.data_array[index][stride_index] =
                        self.data_array[index][stride_index] | (1 << offset);
                    self.data_array[index + self.n][stride_index] =
                        self.data_array[index + self.n][stride_index] | (1 << offset)
                }
                _ => {}
            }
        }
        self.noperators += 1;
        return self.noperators - 1;
    }

    /// Inserts a new Pauli operator described as a vector of bool in the set and returns its index
    pub fn insert_vec_bool(&mut self, axis: &Vec<bool>) -> usize {
        let stride = get_stride(self.noperators + self.start_offset);
        let offset = get_offset(self.noperators + self.start_offset);
        if stride == self.nstrides {
            self.nstrides += 1;
            self.data_array.iter_mut().for_each(|row| row.push(0));
        }
        for (index, value) in axis.iter().enumerate() {
            if *value {
                self.data_array[index][stride] |= 1 << offset;
            }
        }
        self.noperators += 1;
        return self.noperators - 1;
    }
    /// Clears the data of the Pauli set
    pub fn clear(&mut self) {
        for i in 0..2 * self.n {
            for j in 0..self.nstrides {
                self.data_array[i][j] = 0;
            }
        }
        self.noperators = 0;
        self.start_offset = 0;
    }
    /// Pops the first rotation in the set
    pub fn pop(&mut self) {
        let stride = get_stride(self.start_offset);
        let offset = get_offset(self.start_offset);
        for i in 0..2 * self.n {
            self.data_array[i][stride] &= !(1 << offset);
        }
        self.start_offset += 1;
        self.noperators -= 1;
    }
    /// Get the operator at index `operator_index`
    pub fn get(&self, operator_index: usize) -> String {
        let operator_index = operator_index + self.start_offset;
        let mut output = String::new();
        let stride = get_stride(operator_index);
        let offset = get_offset(operator_index);
        for i in 0..self.n {
            match (
                (self.data_array[i][stride] >> offset) & 1,
                (self.data_array[i + self.n][stride] >> offset) & 1,
            ) {
                (1, 0) => {
                    output += "X";
                }
                (0, 1) => {
                    output += "Z";
                }
                (1, 1) => {
                    output += "Y";
                }
                _ => {
                    output += "I";
                }
            }
        }
        output
    }

    /// Get the operator at index `operator_index`
    pub fn get_as_vec_bool(&self, operator_index: usize) -> Vec<bool> {
        let operator_index = operator_index + self.start_offset;
        let mut output = Vec::new();
        let stride = get_stride(operator_index);
        let offset = get_offset(operator_index);
        for i in 0..2 * self.n {
            output.push(((self.data_array[i][stride] >> offset) & 1) != 0);
        }
        output
    }

    /// Returns the support size of the operator (i.e. the number of non-I Pauli term in the operator)
    pub fn support_size(&self, index: usize) -> usize {
        let index = index + self.start_offset;
        let mut count = 0;
        let stride = get_stride(index);
        let offset = get_offset(index);
        for i in 0..self.n {
            if (((self.data_array[i][stride] | self.data_array[i + self.n][stride]) >> offset) & 1)
                != 0
            {
                count += 1;
            }
        }
        return count;
    }

    /*
           Internal methods
    */

    /// XORs row `i` into row `j`
    fn row_op(&mut self, i: usize, j: usize) {
        let (left, right) = self.data_array.split_at_mut(max(i, j));
        let (target_row, source_row) = if i < j {
            (right.get_mut(0).unwrap(), left.get(i).unwrap())
        } else {
            (left.get_mut(j).unwrap(), right.get(0).unwrap())
        };

        for (v1, v2) in source_row.iter().zip(target_row.iter_mut()) {
            *v2 ^= *v1;
        }
    }
    /*
       Gate conjugation
    */

    /// Conjugate the set of rotations via a H gate
    pub fn h(&mut self, i: usize) {
        self.data_array.swap(i, i + self.n);
    }
    /// Conjugate the set of rotations via a S gate
    pub fn s(&mut self, i: usize) {
        self.row_op(i, i + self.n);
    }
    /// Conjugate the set of rotations via a SQRT_X gate
    pub fn sqrt_x(&mut self, i: usize) {
        self.row_op(i + self.n, i);
    }
    /// Conjugate the set of rotations via a CNOT gate
    pub fn cnot(&mut self, i: usize, j: usize) {
        self.row_op(j + self.n, i + self.n);
        self.row_op(i, j);
    }

    /*
       Metrics for synthesis algorithms
    */
    pub fn count_id(&self, qbit: usize) -> usize {
        let mut count: usize = 0;
        for (vx, vz) in self.data_array[qbit]
            .iter()
            .zip(self.data_array[qbit + self.n].iter())
        {
            let value = vx | vz;
            if value == 0 {
                count += 64;
            } else {
                return count + value.trailing_zeros() as usize;
            }
        }
        count
    }
    /// Sorts the set by support size
    pub fn support_size_sort(&mut self) {
        // We first build the "transpose" of data_array (cheaper this way)
        let mut transposed: Vec<Vec<bool>> = (0..self.noperators)
            .map(|i| self.get_as_vec_bool(i))
            .collect();
        transposed.sort_by_key(|vec| {
            (0..self.n)
                .map(|i| if vec[i] | vec[i + self.n] { 1 } else { 0 })
                .sum::<i32>()
        });
        self.clear();
        for axis in transposed {
            self.insert_vec_bool(&axis);
        }
    }
}

#[cfg(test)]
mod pauli_set_tests {
    use super::*;

    #[test]
    fn construction() {
        let pset = PauliSet::new(10);
        assert_eq!(pset.data_array.len(), 20);
        assert_eq!(pset.n, 10);
        assert_eq!(pset.nstrides, 0);
        assert_eq!(pset.noperators, 0);
    }

    #[test]
    fn insertion() {
        let mut pset = PauliSet::new(4);
        pset.insert(&"XYZI");
        assert_eq!(pset.data_array.len(), 8);
        assert_eq!(pset.n, 4);
        assert_eq!(pset.nstrides, 1);
        assert_eq!(pset.noperators, 1);
        assert_eq!(pset.data_array[0][0], 1);
        assert_eq!(pset.data_array[1][0], 1);
        assert_eq!(pset.data_array[2][0], 0);
        assert_eq!(pset.data_array[3][0], 0);

        assert_eq!(pset.data_array[4][0], 0);
        assert_eq!(pset.data_array[5][0], 1);
        assert_eq!(pset.data_array[6][0], 1);
        assert_eq!(pset.data_array[7][0], 0);
    }

    #[test]
    fn get() {
        let mut pset = PauliSet::new(4);
        pset.insert(&"XYZI");
        assert_eq!(pset.get(0), "XYZI");
    }

    #[test]
    fn h_test() {
        let mut pset = PauliSet::new(1);
        pset.insert(&"X");
        pset.insert(&"Z");
        pset.insert(&"Y");
        pset.insert(&"I");

        pset.h(0);

        assert_eq!(pset.get(0), "Z");
        assert_eq!(pset.get(1), "X");
        assert_eq!(pset.get(2), "Y");
        assert_eq!(pset.get(3), "I");
    }

    #[test]
    fn s_test() {
        let mut pset = PauliSet::new(1);
        pset.insert(&"X");
        pset.insert(&"Z");
        pset.insert(&"Y");
        pset.insert(&"I");

        pset.s(0);

        assert_eq!(pset.get(0), "Y");
        assert_eq!(pset.get(1), "Z");
        assert_eq!(pset.get(2), "X");
        assert_eq!(pset.get(3), "I");
    }
    #[test]
    fn cnot_test() {
        let mut pset = PauliSet::new(2);
        pset.insert(&"XX");
        pset.insert(&"ZZ");

        pset.cnot(0, 1);

        assert_eq!(pset.get(0), "XI");
        assert_eq!(pset.get(1), "IZ");
    }
    #[test]
    fn support_size_test() {
        let mut pset = PauliSet::new(4);
        pset.insert(&"XYIZ");
        pset.insert(&"XYII");
        pset.insert(&"IYIZ");
        pset.insert(&"IIII");
        assert_eq!(pset.support_size(0), 3);
        assert_eq!(pset.support_size(1), 2);
        assert_eq!(pset.support_size(2), 2);
        assert_eq!(pset.support_size(3), 0);
    }
    #[test]
    fn count_id() {
        let mut pset = PauliSet::new(4);
        pset.insert(&"IIII");
        pset.insert(&"XIII");
        pset.insert(&"XXII");
        pset.insert(&"XXXI");
        pset.insert(&"XXXX");
        assert_eq!(pset.count_id(0), 1);
        assert_eq!(pset.count_id(1), 2);
        assert_eq!(pset.count_id(2), 3);
        assert_eq!(pset.count_id(3), 4);
    }
    #[test]
    fn sort_test() {
        let mut pset = PauliSet::new(4);
        pset.insert(&"IIII");
        pset.insert(&"XXII");
        pset.insert(&"XXXX");
        pset.insert(&"XIII");
        pset.insert(&"XXXI");
        pset.support_size_sort();
        assert_eq!(pset.get(0), "IIII");
        assert_eq!(pset.get(1), "XIII");
        assert_eq!(pset.get(2), "XXII");
        assert_eq!(pset.get(3), "XXXI");
        assert_eq!(pset.get(4), "XXXX");
    }
    #[test]
    fn pop_test() {
        let mut pset = PauliSet::new(1);
        pset.insert(&"I");
        pset.insert(&"X");
        assert_eq!(pset.noperators, 2);
        pset.pop();
        assert_eq!(pset.noperators, 1);
        assert_eq!(pset.start_offset, 1);
        assert_eq!(pset.get(0), "X");
    }
}
