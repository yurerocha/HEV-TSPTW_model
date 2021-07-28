include("data.jl")
include("model.jl")

function run_model(filename::String)
    d = read_data(filename)
    m = hevtsptw(d)
    # print(m)
    optimize!(m)
    if termination_status(m) == MOI.INFEASIBLE
        println("Infeasible")
    else
        solution_summary(m)
        raw_status(m)
        println("Objective: ", objective_value(m) - sum(m[max_wb[i]] for i in 0:d.ncustomers) + sum(m[min_lmax[i]] for i in 0:d.ncustomers))
    end
end
