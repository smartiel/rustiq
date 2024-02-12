use crate::structures::CliffordGate;

pub type Chunk = [Option<CliffordGate>; 3];
pub const ALL_CHUNKS: [Chunk; 18] = [
    [None, None, Some(CliffordGate::CNOT(0, 1))],
    [None, None, Some(CliffordGate::CNOT(1, 0))],
    [
        None,
        Some(CliffordGate::H(1)),
        Some(CliffordGate::CNOT(0, 1)),
    ],
    [
        None,
        Some(CliffordGate::H(0)),
        Some(CliffordGate::CNOT(1, 0)),
    ],
    [
        None,
        Some(CliffordGate::S(1)),
        Some(CliffordGate::CNOT(0, 1)),
    ],
    [
        None,
        Some(CliffordGate::S(0)),
        Some(CliffordGate::CNOT(1, 0)),
    ],
    [
        Some(CliffordGate::H(0)),
        None,
        Some(CliffordGate::CNOT(0, 1)),
    ],
    [
        Some(CliffordGate::H(1)),
        None,
        Some(CliffordGate::CNOT(1, 0)),
    ],
    [
        Some(CliffordGate::H(0)),
        Some(CliffordGate::H(1)),
        Some(CliffordGate::CNOT(0, 1)),
    ],
    [
        Some(CliffordGate::H(1)),
        Some(CliffordGate::H(0)),
        Some(CliffordGate::CNOT(1, 0)),
    ],
    [
        Some(CliffordGate::H(0)),
        Some(CliffordGate::S(1)),
        Some(CliffordGate::CNOT(0, 1)),
    ],
    [
        Some(CliffordGate::H(1)),
        Some(CliffordGate::S(0)),
        Some(CliffordGate::CNOT(1, 0)),
    ],
    [
        Some(CliffordGate::SqrtX(0)),
        None,
        Some(CliffordGate::CNOT(0, 1)),
    ],
    [
        Some(CliffordGate::SqrtX(1)),
        None,
        Some(CliffordGate::CNOT(1, 0)),
    ],
    [
        Some(CliffordGate::SqrtX(0)),
        Some(CliffordGate::H(1)),
        Some(CliffordGate::CNOT(0, 1)),
    ],
    [
        Some(CliffordGate::SqrtX(1)),
        Some(CliffordGate::H(0)),
        Some(CliffordGate::CNOT(1, 0)),
    ],
    [
        Some(CliffordGate::SqrtX(0)),
        Some(CliffordGate::S(1)),
        Some(CliffordGate::CNOT(0, 1)),
    ],
    [
        Some(CliffordGate::SqrtX(1)),
        Some(CliffordGate::S(0)),
        Some(CliffordGate::CNOT(1, 0)),
    ],
];
