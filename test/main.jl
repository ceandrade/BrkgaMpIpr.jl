push!(LOAD_PATH, ".")

import BrkgaMpIpr
import TestInstance
import Decoder

chr_size = 100

print("\n\n>> Building instance")
instance = TestInstance.Instance(chr_size)

print("\n\n>> Building BRKGA data")
brkga_data = BrkgaMpIpr.init(chr_size, chr_size)

print("\n\n>> Evolving")
@time BrkgaMpIpr.evolve!(brkga_data, instance, Decoder.decode!)

