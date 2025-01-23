using SeisTools, LinearAlgebra, DelimitedFiles

include("../julia/FMI_IO.jl")

function gf(x1, y1, z1, x2, y2, z2, vp, vs, M0, npts, hw)
    gaussw = exp.(-range(-3.0, 3.0, hw * 2 + 1) .^ 2)
    rv = [x2 - x1, y2 - y1, z2 - z1]
    r = norm(rv)
    rn = rv ./ r
    tp = r / vp
    ts = r / vs
    ip = round(Int, tp / dt)
    is = round(Int, ts / dt)
    w = zeros(npts, 6)
    Rp = zeros(6)
    Rp[1] = rn[3] * rn[1] * rn[1]
    Rp[2] = rn[3] * rn[2] * rn[2]
    Rp[3] = rn[3] * rn[3] * rn[3]
    Rp[4] = rn[3] * rn[1] * rn[2]
    Rp[5] = rn[3] * rn[1] * rn[3]
    Rp[6] = rn[3] * rn[2] * rn[3]
    Rs = zeros(6)
    Rs[1] = (0.0 - rn[3] * rn[1])* rn[1]
    Rs[2] = (0.0 - rn[3] * rn[2])* rn[2]
    Rs[3] = (1.0 - rn[3] * rn[3])* rn[3]
    Rs[4] = (0.0 - rn[3] * rn[1])* rn[2]
    Rs[5] = (0.0 - rn[3] * rn[1])* rn[3]
    Rs[6] = (0.0 - rn[3] * rn[2])* rn[3]

    for ic = 1:6
        w[ip:ip+2*hw, ic] .+= (Rp[ic] * M0 / r / vp) .* gaussw
        w[is:is+2*hw, ic] .+= (Rs[ic] * M0 / r / vs) .* gaussw
    end
    return (ip, is, w)
end

n_receiver = 6
npts = 400
dt = 0.01
halfwidth = 20
M0 = 30.0
sdr = SeisTools.Source.SDR(10, 20, 30)
mt = SeisTools.Source.MomentTensor(sdr)
mtvec = collect(mt.values)

receivers = randn(2, n_receiver)
writedlm("receivers.txt", receivers)
source = [0.1, -0.1, 3.2]
vp = 5.0
vs = 3.0

x_search = -3:1.0:3.0
y_search = -3:1.0:3.0
z_search = 1.0:1.0:5

traces = map(1:n_receiver) do ir
    (ip, is, w) = gf(source[1], source[2], source[3], receivers[1, ir], receivers[2, ir], 0.0, vp, vs, M0, npts, halfwidth)
    return FMI_IO.Trace(npts, 0, ip, 2 * halfwidth, is, 2 * halfwidth, w * mtvec .+ randn(npts) .* 0.001)
end

# sourclocs = CartesianIndices((length(z_search), length(y_search), length(x_search))) |> vec

mts = map(vec(CartesianIndices((0:10:350, 0:5:90, -90:10:90)))) do ci
    mt0 = SeisTools.Source.MomentTensor(SeisTools.Source.SDR(ci[1], ci[2], ci[3]))
    return FMI_IO.MomentTensor(mt0.values...)
end

gf_db = map(vec(CartesianIndices((length(z_search), length(y_search), length(x_search), n_receiver)))) do ci
    z = z_search[ci[1]]
    y = y_search[ci[2]]
    x = x_search[ci[3]]
    (ip, is, w) = gf(x, y, z, receivers[1, ci[4]], receivers[2, ci[4]], 0.0, vp, vs, M0, npts, halfwidth)
    return FMI_IO.GFtrace(npts, ip, is, w[:, 1], w[:, 2], w[:, 3], w[:, 4], w[:, 5], w[:, 6])
end

open(io -> FMI_IO.write_c_input(io, traces, mts, gf_db), "../build/input.bin", "w")
open(io -> FMI_IO.write_c_input(io, traces, mts, gf_db), "../build/input_omp.bin", "w")
