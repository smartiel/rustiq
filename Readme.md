[![Rust](https://github.com/smartiel/rustiq/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/smartiel/rustiq/actions/workflows/rust.yml)

# Rustiq: A quantum circuit compiler in Rust

This library contains implementations of various quantum circuit synthesis algorithms.
It can handle:
 - Stabilizer and Graph state synthesis
 - Synthesis of codiagonalization circuits for sets of commuting Pauli operators
 - Pauli network synthesis (i.e. synthesis of sequences of Pauli rotations)
 - Clifford operator and Clifford isometry synthesis



## Python library

All the synthesis methods are binded in python and the repo can be installed via pip:

```bash
pip install . 
```

You might need to install rust first.

## Building the project (Rust)

The project is also available as a rust crate:

```bash
cargo build --release
```


## Running the synthesis algorithms

### In python

See the `examples` folder for examples of python usage (the file names are pretty self-explainatory)
