include("data.jl")
include("model.jl")

function run_model(filename::String)
    d = read_data(filename)
    m = hevtsptw(d)
end
