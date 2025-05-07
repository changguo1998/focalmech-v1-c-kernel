d = abspath(@__DIR__, "../out/build/x64-Debug")

for f in filter(endswith("exe"), readdir(d))
    cp(joinpath(d, f), joinpath(@__DIR__, f); force=true)
end
