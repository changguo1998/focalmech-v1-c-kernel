using SeisTools, CairoMakie, LinearAlgebra, Statistics, DelimitedFiles

include("../julia/FMI_IO.jl")

(trs, mts, gfs) = FMI_IO.read_c_input("input_omp.bin");

fig1 = Figure();
ax1 = Axis(fig1[1,1]);
for i = eachindex(trs)
    w = trs[i].obs;
    lines!(ax1, w ./ std(w) ./ 2.0 .+ i);
    lines!(ax1, [1.0, 1.0].* trs[i].pwin_i, i .+ [-0.5, 0.5], color=:blue);
    lines!(ax1, [1.0, 1.0].* trs[i].swin_i, i .+ [-0.5, 0.5], color=:red);
end
save("traces.png", fig1);

t = FMI_IO.read_c_result("output_omp.bin");

misfit = dropdims(sum(t.pl2 + t.sl2, dims=1), dims=1);

(min_misfit, min_index) = findmin(misfit);

# misfit_resample = misfit[round.(Int, (0:2000)./2000.0.*(size(misfit, 1)-1)).+1, :]

fig2 = Figure();
ax2 = Axis(fig2[1,1]);
hm = heatmap!(ax2, misfit);
save("misfit.png", fig2);

recvs = readdlm("receivers.txt")

fig3 = Figure();
ax3 = Axis(fig3[1,1]);
scatter!(recvs[1,:], recvs[2,:]);
scatter!([0.1], [-0.1], marker=:star5, color=:red)
save("receivers.png", fig3);

imech = min_index[1]
isrc = min_index[2]

mt = mts[imech]

sdr = SeisTools.Source.SDR(SeisTools.Source.MomentTensor(mt.m11, mt.m22, mt.m33, mt.m12, mt.m13, mt.m23))

x_search = -3:1.0:3.0
y_search = -3:1.0:3.0
z_search = 1.0:1.0:5

sourclocs = CartesianIndices((length(z_search), length(y_search), length(x_search))) |> vec
isrc_cart = sourclocs[isrc]

x_result = x_search[isrc_cart[3]]
y_result = y_search[isrc_cart[2]]
z_result = z_search[isrc_cart[1]]

[x_result y_result z_result]
