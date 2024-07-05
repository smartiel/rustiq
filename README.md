[![Rust](https://github.com/smartiel/rustiq/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/smartiel/rustiq/actions/workflows/rust.yml)

# Rustiq: A python interface for quantum circuit synthesis routines

This library contains implementations of various quantum circuit synthesis algorithms.
It can handle:
 - Pauli network synthesis (i.e. synthesis of sequences of Pauli rotations)
 - Stabilizer and Graph state synthesis
 - Synthesis of codiagonalization circuits for sets of commuting Pauli operators
 - Clifford operators and Clifford isometries synthesis



## Python library

All the synthesis methods are binded in python and the repo can be installed via pip:

```bash
pip install git+https://github.com/smartiel/rustiq.git@main
```
or
```bash
git clone https://github.com/smartiel/rustiq.git
pip install ./rustiq
```

You might need to install rust first.

## Underlying rust crate

This project is essentially a python wrapper around the rustiq-core crate



## Running the synthesis algorithms


See the `examples` folder for examples of python usage (the file names are pretty self-explainatory).

Most of the synthesis algorithms are binded and wrapped in a nice python interface.


## References and how to cite

The algorithms implemented in this package have been developed in the following papers:

For Pauli network synthesis: 
```bibtex
@misc{debrugière2024faster,
      title={Faster and shorter synthesis of Hamiltonian simulation circuits}, 
      author={Timothée Goubault de Brugière and Simon Martiel},
      year={2024},
      eprint={2404.03280},
      archivePrefix={arXiv},
      primaryClass={quant-ph}
}
```

For graph/stabilizer state synthesis, Clifford operators, Clifford isometries, co-diagonalization:
```bibtex
@misc{debrugière2022graphstate,
      title={A graph-state based synthesis framework for Clifford isometries}, 
      author={Timothée Goubault de Brugière and Simon Martiel and Christophe Vuillot},
      year={2022},
      eprint={2212.06928},
      archivePrefix={arXiv},
      primaryClass={quant-ph}
}
```