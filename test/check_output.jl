include("../julia/FMI_IO.jl")

t = FMI_IO.read_c_result("../build/output.bin")
