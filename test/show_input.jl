include("../julia/FMI_IO.jl")

(trs, mts, gfs) = FMI_IO.read_c_input("../build/input.bin")

println("Traces($(length(trs))):")
for i = eachindex(trs)
    tr = trs[i]
    println(join([string.([i, tr.n_samples, tr.ref_i, tr.pwin_i, tr.pwin_l, tr.swin_i, tr.swin_l]); string.(tr.obs[1:5])], ' '))
end

println("\nMomentTensors($(length(mts))):")
# for i = 1:min(5, length(mts))
for i = 1:length(mts)
    mt = mts[i]
    println(join(string.([mt.m11, mt.m22, mt.m33, mt.m12, mt.m13, mt.m23]), " "))
end

println("\nGF($(length(gfs))):")
for i = 1:length(gfs)
    gf = gfs[i]
    println(join([string.([gf.n_samples, gf.p_i, gf.s_i]); string.(gf.g11[1:5])], " "))
end
