t = randn(Float32, 5)
open(io->write(io, t), "../build/arrayio.bin", "w")
println.(t)