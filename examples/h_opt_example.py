from rustiq.h_gadget import diagonalization_network, gadgetize, predict_nbqbits


with open("uccsd.txt", "r") as fin:
    paulis = [line.strip("\n") for line in fin.readlines()]


def chunkize(paulis, n_h_ancillae, h_optimal=False):
    nqbits = len(paulis[0])
    # print("# of logical qubits:", nqbits)
    # print("# of H ancillae:", n_h_ancillae)
    total = nqbits + n_h_ancillae
    chunks = []
    start = 0
    end = 1
    while True:
        while predict_nbqbits(paulis[start:end], h_optimal) <= total and end <= len(
            paulis
        ):
            end += 1
        chunks.append(paulis[start : end - 1])
        if end == len(paulis) + 1:
            break
        start = end - 1
        end = start + 1
        if start == len(paulis):
            break
    # print("# of phase polynomials:", len(chunks))
    return chunks


print(predict_nbqbits(paulis, False))
print(predict_nbqbits(paulis, True))
for n_anc in range(1, 31):
    print("=" * 40)
    chunks = chunkize(paulis, n_anc)
    print(max([predict_nbqbits(chunk, False) for chunk in chunks if chunk]))
    print(max([predict_nbqbits(chunk, True) for chunk in chunks if chunk]))
    chunks_opt = chunkize(paulis, n_anc, True)
    # print([predict_nbqbits(chunk, False) for chunk in chunks_opt if chunk])
    # print([predict_nbqbits(chunk, True) for chunk in chunks_opt if chunk])
    print(f"{n_anc},{len(chunks)},{len(chunks_opt)}")
