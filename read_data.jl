open("input_data") do file
    lines = readlines(file)
    println(parse(Float64, lines))
end

IJulia.load("input_data")

cfl    = 0.2
tfinal = 1.
nstepmax = 500    
idiag  = 10
md = 2
nd = 2
nx = 100
ny = 100
dimx = 1.0
dimy = 1.0


using JLD
jldopen("input_data.jld", "w") do file
    JLD.@write file cfl
    JLD.@write file tfinal
    JLD.@write file nstepmax    
    JLD.@write file idiag
    JLD.@write file md
    JLD.@write file nd
    JLD.@write file nx
    JLD.@write file ny
    JLD.@write file dimx
    JLD.@write file dimy
end


data = load("input_data.jld")

using JSON
open("input_data.json","w") do f
    JSON.print(f, data)
end

!cat input_data.json

open("input_data.json", "r") do f
    data = read(f, String)
    println(JSON.parse(data))
end


