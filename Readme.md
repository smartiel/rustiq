[![Rust](https://github.com/smartiel/qrust/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/smartiel/qrust/actions/workflows/rust.yml)

# A quantum circuit compiler in Rust

## Building the project

```bash
cargo build --release
```

## Running the synthesis algorithm

The directory `benchmarks/chem/` contains text files describing set of Pauli operators used in UCCSD Ansätze.
The file names are pretty straightforward to parse.

Once the crate is compiled, you can run the synthesizer as follows:

```bash
./target/release/qrust <FILENAME> 
```

By default, the program prints the output circuit. You can simply ask the program for the running time, the CNOT count, and the CNOT depth
by adding the flag `--onlyinfo` (the circuits might be huge).

```bash
./target/release/qrust <FILENAME> --onlyinfo
```

You can switch between the count and depth minimizing algorithms via the `--metric` option:

```bash
./target/release/qrust <FILENAME> --onlyinfo --metric=depth
```

or

```bash
./target/release/qrust <FILENAME> --onlyinfo --metric=count
```

