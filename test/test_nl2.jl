using DelimitedFiles, GLMakie

samples = map(readlines("test_nl2.txt")) do l
    t = split(l, ' ')
    return (parse(Float64, t[1]), parse(Float64, t[3]), parse(Float64, t[5]))
end

(wr, ws) = let
    buf = filter(!isempty, readlines("test_nl2_raw.txt"))
    i = findfirst(contains("="), buf)
    (parse.(Float64, buf[1:i-1]), parse.(Float64, buf[i+1:end]))
end

x = getindex.(samples, 1)
y = getindex.(samples, 2)
dxy = getindex.(samples, 3)

fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, wr)
lines!(ax1, ws)
ax2 = Axis(fig[2, 1])
scatter!(ax2, x, color=:black)
lines!(ax2, y, color=:blue)
scatter!(ax2, dxy, color=:red, markersize=2)
ax3 = Axis(fig[3, 1])
scatter!(ax3, dxy, color=:red, markersize=4)
