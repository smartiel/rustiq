from rustiq.h_gadget import diagonalization_network, gadgetize

example_ali = ["ZZ", "ZX"]

network = diagonalization_network(example_ali)
begin, (phase_pol, clifford_correction), end = gadgetize(network)
print("Preparation:")
print(begin)
print("Phase pol:")
print(phase_pol)
print("Measurements:")
print(end)
