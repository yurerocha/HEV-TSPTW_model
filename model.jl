using JuMP, CPLEX

include("data.jl")

function hevtsptw(d::Data)
    m = Model(CPLEX.Optimizer)

    M = 1_000_000
    n = d.ncustomers

    @variable(m, x[i in 0:n, j in 1:n + 1; i != j], binary=true)
    @variable(m, x_cc[i in 0:n, j in 1:n + 1; i != j], binary=true)
    @variable(m, x_b[i in 0:n, j in 1:n + 1; i != j], binary=true)
    @variable(m, x_c[i in 0:n, j in 1:n + 1; i != j], binary=true)
    @variable(m, x_e[i in 0:n, j in 1:n + 1; i != j], binary=true)

    @variable(m, w[0:n + 1])
    @variable(m, l[0:n + 1])

    @variable(m, max_wb[0:n]) # max(wi, bi)
    @variable(m, min_lmax[0:n]) # min(li + tcc_ij*rc, lmax)

    ex_cc = @expression(m, sum((i != j ? x_cc[i, j] * d.costs["cc"][i, j] : 0) for i in 0:n, j in 1:n + 1))
    ex_b = @expression(m, sum((i != j ? x_b[i, j] * d.costs["b"][i, j] : 0) for i in 0:n, j in 1:n + 1))
    ex_c = @expression(m, sum((i != j ? x_c[i, j] * d.costs["c"][i, j] : 0) for i in 0:n, j in 1:n + 1))
    ex_e = @expression(m, sum((i != j ? x_e[i, j] * d.costs["e"][i, j] : 0) for i in 0:n, j in 1:n + 1))

    @objective(m, Min, ex_cc + ex_b + ex_c + ex_e + sum(max_wb[i] for i in 0:n) - sum(min_lmax[i] for i in 0:n))

    @constraint(m, con2[i in 0:n, j in 1:n + 1; i != j], x_cc[i, j] + x_b[i, j] + x_c[i, j] + x_e[i, j] == x[i, j])

    @constraint(m, con3[j in 1:n], sum((i != j ? x[i, j] : 0) for i in 0:n) == 1)
    @constraint(m, con4[i in 1:n], sum((i != j ? x[i, j] : 0) for j in 1:n + 1) == 1)

    @constraint(m, con5, sum(x[0, j] for j in 1:n) == 1)
    @constraint(m, con6, sum(x[i, n + 1] for i in 1:n) == 1)

    @constraint(m, con7, w[0] == 0)

    # Time flow constraints.
    @constraint(m, [i in 0:n], max_wb[i] >= w[i]) # max(wi, bi)
    @constraint(m, [i in 0:n], max_wb[i] >= d.tw_b[i]) # max(wi, bi)
    # @expression(m, tf_ex[i in 0:n, j in 1:n + 1; i != j], max_wb[i] + d.service[i] - w[j])

    @constraint(m, con8a[i in 0:n, j in 1:n + 1; i != j], max_wb[i] + d.service[i] + d.times["cc"][i, j] - w[j] <= (1 - x_cc[i, j]) * M)
    @constraint(m, con8b[i in 0:n, j in 1:n + 1; i != j], max_wb[i] + d.service[i] + d.times["cc"][i, j] - w[j] >= -(1 - x_cc[i, j]) * M)

    @constraint(m, con9a[i in 0:n, j in 1:n + 1; i != j], max_wb[i] + d.service[i] + d.times["b"][i, j] - w[j] <= (1 - x_b[i, j]) * M)
    @constraint(m, con9b[i in 0:n, j in 1:n + 1; i != j], max_wb[i] + d.service[i] + d.times["b"][i, j] - w[j] >= -(1 - x_b[i, j]) * M)

    @constraint(m, con10a[i in 0:n, j in 1:n + 1; i != j], max_wb[i] + d.service[i] + d.times["c"][i, j] - w[j] <= (1 - x_c[i, j]) * M)
    @constraint(m, con10b[i in 0:n, j in 1:n + 1; i != j], max_wb[i] + d.service[i] + d.times["c"][i, j] - w[j] >= -(1 - x_c[i, j]) * M)

    @constraint(m, con11a[i in 0:n, j in 1:n + 1; i != j], max_wb[i] + d.service[i] + d.times["e"][i, j] - w[j] <= (1 - x_e[i, j]) * M)
    @constraint(m, con11b[i in 0:n, j in 1:n + 1; i != j], max_wb[i] + d.service[i] + d.times["e"][i, j] - w[j] >= -(1 - x_e[i, j]) * M)
    
    @constraint(m, con12[j in 1:n + 1], w[j] + d.service[j] <= d.tw_e[j])
    @constraint(m, con13, l[0] == 0)
        
    # Battery charging constraints.
    @constraint(m, [i in 0:n, j in 1:n + 1; i != j], min_lmax[i] <= l[i] + d.times["cc"][i, j] * d.rate_c) # min(li + tcc_ij*rc, lmax).
    @constraint(m, [i in 0:n], min_lmax[i] <= d.l_max) # min(li + tcc_ij*rc, lmax)

    @constraint(m, con14a[i in 0:n, j in 1:n + 1; i != j], min_lmax[i] - l[j] <= (1 - x_cc[i, j]) * M)
    @constraint(m, con14b[i in 0:n, j in 1:n + 1; i != j], min_lmax[i] - l[j] >= -(1 - x_cc[i, j]) * M)

    @constraint(m, con15a[i in 0:n, j in 1:n + 1; i != j], l[i] - d.times["b"][i, j] * d.rate_d - l[j] <= (1 - x_b[i, j]) * M)
    @constraint(m, con15b[i in 0:n, j in 1:n + 1; i != j], l[i] - d.times["b"][i, j] * d.rate_d - l[j] >= -(1 - x_b[i, j]) * M)

    @constraint(m, con16a[i in 0:n, j in 1:n + 1; i != j], l[i] - l[j] <= (1 - x_c[i, j]) * M)
    @constraint(m, con16b[i in 0:n, j in 1:n + 1; i != j], l[i] - l[j] >= -(1 - x_c[i, j]) * M)

    @constraint(m, con17a[i in 0:n, j in 1:n + 1; i != j], l[i] - d.times["e"][i, j] * d.rate_d - l[j] <= (1 - x_e[i, j]) * M)
    @constraint(m, con17b[i in 0:n, j in 1:n + 1; i != j], l[i] - d.times["e"][i, j] * d.rate_d - l[j] >= -(1 - x_e[i, j]) * M)

    @constraint(m, con18[i in 0:n + 1], 0 <= l[i] <= d.l_max)

    return m
end
