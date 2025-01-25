module FMI_IO

import Base: read, write

struct Trace
    n_samples::Int32
    ref_i::Int32
    pwin_i::Int32
    pwin_l::Int32
    swin_i::Int32
    swin_l::Int32
    obs::Vector{Float32}
end

function Trace(n_samples::Integer, ref_i::Integer,
    pwin_i::Integer, pwin_l::Integer,
    swin_i::Integer, swin_l::Integer,
    obs::AbstractVector{<:AbstractFloat})
    return Trace(Int32(n_samples), Int32(ref_i), Int32(pwin_i), Int32(pwin_l),
        Int32(swin_i), Int32(swin_l), Float32.(obs))
end

function read(io::IO, ::Type{Trace})
    n_samples = read(io, Int32)
    ref_i = read(io, Int32)
    pwin_i = Base.read(io, Int32)
    pwin_l = read(io, Int32)
    swin_i = read(io, Int32)
    swin_l = read(io, Int32)
    obs = zeros(Float32, n_samples)
    read!(io, obs)
    return Trace(n_samples, ref_i, pwin_i, pwin_l, swin_i, swin_l, obs)
end

function write(io::IO, tr::Trace)
    # @info "n_samples"
    write(io, tr.n_samples)
    write(io, tr.ref_i)
    write(io, tr.pwin_i)
    write(io, tr.pwin_l)
    write(io, tr.swin_i)
    write(io, tr.swin_l)
    write(io, tr.obs)
end

struct MomentTensor
    m11::Float32
    m22::Float32
    m33::Float32
    m12::Float32
    m13::Float32
    m23::Float32
end

function MomentTensor(m11::AbstractFloat, m22::AbstractFloat, m33::AbstractFloat,
    m12::AbstractFloat, m13::AbstractFloat, m23::AbstractFloat)
    return MomentTensor(Float32(m11), Float32(m22), Float32(m33),
        Float32(m12), Float32(m13), Float32(m23))
end

function read(io::IO, ::Type{MomentTensor})
    m11 = read(io, Float32)
    m22 = read(io, Float32)
    m33 = read(io, Float32)
    m12 = read(io, Float32)
    m13 = read(io, Float32)
    m23 = read(io, Float32)
    return MomentTensor(m11, m22, m33, m12, m13, m23)
end

function write(io::IO, s::MomentTensor)
    write(io, s.m11)
    write(io, s.m22)
    write(io, s.m33)
    write(io, s.m12)
    write(io, s.m13)
    write(io, s.m23)
end

struct GFtrace
    n_samples::Int32
    p_i::Int32
    s_i::Int32
    g11::Vector{Float32}
    g22::Vector{Float32}
    g33::Vector{Float32}
    g12::Vector{Float32}
    g13::Vector{Float32}
    g23::Vector{Float32}
end

function GFtrace(n_samples::Integer, p_i::Integer, s_i::Integer,
    g11::AbstractVector{<:AbstractFloat},
    g22::AbstractVector{<:AbstractFloat},
    g33::AbstractVector{<:AbstractFloat},
    g12::AbstractVector{<:AbstractFloat},
    g13::AbstractVector{<:AbstractFloat},
    g23::AbstractVector{<:AbstractFloat}
)
    return GFtrace(Int32(n_samples), Int32(p_i), Int32(s_i),
        Float32.(g11), Float32.(g22), Float32.(g33), Float32.(g12), Float32.(g13), Float32.(g23))
end

function read(io::IO, ::Type{GFtrace})
    n = read(io, Int32)
    p_i = read(io, Int32)
    s_i = read(io, Int32)
    g11 = zeros(Float32, n)
    g22 = zeros(Float32, n)
    g33 = zeros(Float32, n)
    g12 = zeros(Float32, n)
    g13 = zeros(Float32, n)
    g23 = zeros(Float32, n)
    read!(io, g11)
    read!(io, g22)
    read!(io, g33)
    read!(io, g12)
    read!(io, g13)
    read!(io, g23)
    return GFtrace(n, p_i, s_i, g11, g22, g33, g12, g13, g23)
end

function write(io::IO, gf::GFtrace)
    write(io, gf.n_samples)
    write(io, gf.p_i)
    write(io, gf.s_i)
    write(io, gf.g11)
    write(io, gf.g22)
    write(io, gf.g33)
    write(io, gf.g12)
    write(io, gf.g13)
    write(io, gf.g23)
end

function read_c_input(io::IO)
    ntrs = read(io, Int32)
    @info "ntrs: $ntrs"
    nmts = read(io, Int32)
    @info "nmts: $nmts"
    nsrc = read(io, Int32)
    @info "nsrc: $nsrc"
    ngfs = ntrs * nsrc
    traces = Vector{Trace}(undef, ntrs)
    mts = Vector{MomentTensor}(undef, nmts)
    gf = Vector{GFtrace}(undef, ngfs)

    for i = eachindex(traces)
        traces[i] = read(io, Trace)
    end
    for i = eachindex(mts)
        mts[i] = read(io, MomentTensor)
    end
    for i = eachindex(gf)
        gf[i] = read(io, GFtrace)
    end
    return (traces, mts, gf)
end

read_c_input(fname::AbstractString) = open(read_c_input, fname, "r")

function write_c_input(io::IO, traces::AbstractVector{Trace}, mts::AbstractVector{MomentTensor}, gf::AbstractVector{GFtrace})
    write(io, Int32(length(traces)))
    write(io, Int32(length(mts)))
    write(io, Int32(length(gf)/length(traces)))
    for i = eachindex(traces)
        write(io, traces[i])
    end
    for i = eachindex(mts)
        write(io, mts[i])
    end
    for i = eachindex(gf)
        write(io, gf[i])
    end
    return nothing
end

write_c_input(fname::AbstractString, traces::AbstractVector{Trace}, mts::AbstractVector{MomentTensor}, gf::AbstractVector{GFtrace}) =
    open(io -> write_c_input(io, traces, mts, gf), fname, "w")

function read_c_result(io::IO)
    ntrs = read(io, Int32)
    nmts = read(io, Int32)
    nsrc = read(io, Int32)
    misfit_pl2 = zeros(Float32, ntrs, nmts, nsrc)
    misfit_pshift = zeros(Int32, ntrs, nmts, nsrc)
    misfit_sl2 = zeros(Float32, ntrs, nmts, nsrc)
    misfit_sshift = zeros(Int32, ntrs, nmts, nsrc)
    misfit_pol = zeros(Float32, ntrs, nmts, nsrc)
    misfit_psr = zeros(Float32, ntrs, nmts, nsrc)
    read!(io, misfit_pl2)
    read!(io, misfit_pshift)
    read!(io, misfit_sl2)
    read!(io, misfit_sshift)
    read!(io, misfit_pol)
    read!(io, misfit_psr)
    return (pl2=misfit_pl2, pshift=misfit_pshift, sl2=misfit_sl2, sshift=misfit_sshift, polarity=misfit_pol, psr=misfit_psr)
end

read_c_result(fname::AbstractString) = open(read_c_result, fname, "r")

end
