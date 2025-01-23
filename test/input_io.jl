include("../julia/FMI_IO.jl")

v = Float32.(collect(1:1000))

traces = [FMI_IO.Trace(length(v), 0, 1, 2, 3, 4, v), FMI_IO.Trace(length(v), 1, 2, 3, 4, 5, v)]
mts = [FMI_IO.MomentTensor((1:6)...)]
gfs = [FMI_IO.GFtrace(length(v), 1, 2, v, v, v, v, v, v)]

FMI_IO.write_c_input("../build/input.bin", traces, mts, gfs)
