using SeisTools, GLMakie, Statistics, LinearAlgebra
include("../julia/FMI_IO.jl")

(trs, mts, gfs) = FMI_IO.read_c_input("input_cuda.bin");
result_omp = FMI_IO.read_c_result("output_omp.bin");
result_cuda = FMI_IO.read_c_result("output_cuda.bin");

for f = [:pl2, :pshift, :sl2, :sshift, :polarity, :psr]
    println(string(f), " ", all( getfield(result_omp, f) .== getfield(result_cuda, f)) ? "True" : "False")
end

dmis = result_cuda.pl2 - result_omp.pl2

a0 = 0.005

fig = Figure();
ax11 = Axis(fig[1,1]);
hm1 = heatmap!(ax11, result_cuda.pl2[2,:,:], colorrange=(0.0, a0));
Colorbar(fig[1,2], hm1);
ax21 = Axis(fig[2,1]);
hm2 = heatmap!(ax21, result_omp.pl2[2,:,:], colorrange=(0.0, a0));
Colorbar(fig[2,2], hm2);
ax31 = Axis(fig[3,1]);
hm3 = heatmap!(ax31, dmis[2,:,:], colorrange=(-0.0001, 0.0001), colormap=:bwr);
Colorbar(fig[3,2], hm3);
ax41 = Axis(fig[4,1]);
hm4 = heatmap!(ax41, abs.(dmis[2,:,:])./result_omp.pl2[2,:,:], colorrange=(0.0, 0.1),colormap=:binary);
Colorbar(fig[4,2], hm4);
fig
