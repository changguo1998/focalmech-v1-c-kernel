using SeisTools, GLMakie, Statistics, LinearAlgebra, Printf
include("../julia/FMI_IO.jl")

(trs, mts, gfs) = FMI_IO.read_c_input("input_cuda.bin");
result_omp = FMI_IO.read_c_result("output_omp.bin");
result_cuda = FMI_IO.read_c_result("output_cuda.bin");

for f = [:pl2, :pshift, :sl2, :sshift, :polarity, :psr]
    println(string(f), " ", all( getfield(result_omp, f) .== getfield(result_cuda, f)) ? "True" : "False")
end

mis_cuda = result_cuda.psr[2, :, :];
mis_omp = result_omp.psr[2, :, :];
dmis = mis_cuda - mis_omp;
rdmis = 2 .* abs.(dmis) ./ (mis_cuda + mis_omp);

fig = Figure(size=(1200, 900));
ax11 = Axis(fig[1,1]; ylabel="CUDA misfit");
hm1 = heatmap!(ax11, mis_cuda, colormap=[:white, :red]);
text!(ax11, 0, size(mis_cuda, 2), text=
    @sprintf("min: %.2f\nmax: %.2f", minimum(mis_cuda), maximum(mis_cuda)),
    align=(:left, :top));
Colorbar(fig[1,2], hm1);
ax21 = Axis(fig[2,1]; ylabel="OpenMP misfit");
hm2 = heatmap!(ax21, mis_omp, colormap=[:white, :red]);
text!(ax21, 0, size(mis_omp, 2), text=
    @sprintf("min: %.2f\nmax: %.2f", minimum(mis_omp), maximum(mis_omp)),
    align=(:left, :top));
Colorbar(fig[2,2], hm2);
ax31 = Axis(fig[3,1]; ylabel="log(abs(CUDA - OpenMP))");
hm3 = heatmap!(ax31, log.(abs.(dmis)), colormap=[:white, :red]);
text!(ax31, 0, size(dmis, 2), text=
    @sprintf("min: %.2e\nmax: %.2e", minimum(dmis), maximum(dmis)),
    align=(:left, :top));
Colorbar(fig[3,2], hm3);
ax41 = Axis(fig[4,1]; ylabel="relative misfit");
hm4 = heatmap!(ax41, rdmis, colormap=:binary);
text!(ax41, 0, size(rdmis, 2), text=
    @sprintf("min: %.2f\nmax: %.2f", minimum(rdmis), maximum(rdmis)),
    align=(:left, :top));
Colorbar(fig[4,2], hm4);
fig
