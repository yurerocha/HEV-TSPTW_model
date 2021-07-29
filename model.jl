using JuMP, CPLEX

include("data.jl")

function hevtsptw(d::Data)
    m = Model(CPLEX.Optimizer)

    M = 1_000_000_000
    n = d.ncustomers

    @variable(m, x[i in 0:n, j in 1:n + 1; i != j], binary=true)
    @variable(m, x_cc[i in 0:n, j in 1:n + 1; i != j], binary=true)
    @variable(m, x_b[i in 0:n, j in 1:n + 1; i != j], binary=true)
    @variable(m, x_c[i in 0:n, j in 1:n + 1; i != j], binary=true)
    @variable(m, x_e[i in 0:n, j in 1:n + 1; i != j], binary=true)

    @variable(m, w[0:n + 1])
    @variable(m, l[0:n + 1])

    @variable(m, max_wb[i in 0:n, j in 1:n + 1; i != j]) # max(wi, bi)
    @variable(m, min_lmax[i in 0:n, j in 1:n + 1; i != j]) # min(li + tcc_ij*rc, lmax)

    ex_cc = @expression(m, sum((i != j ? x_cc[i, j] * d.costs["cc"][i, j] : 0.0) for i in 0:n, j in 1:n + 1))
    ex_b = @expression(m, sum((i != j ? x_b[i, j] * d.costs["b"][i, j] : 0.0) for i in 0:n, j in 1:n + 1))
    ex_c = @expression(m, sum((i != j ? x_c[i, j] * d.costs["c"][i, j] : 0.0) for i in 0:n, j in 1:n + 1))
    ex_e = @expression(m, sum((i != j ? x_e[i, j] * d.costs["e"][i, j] : 0.0) for i in 0:n, j in 1:n + 1))

    @objective(m, Min, ex_cc + ex_b + ex_c + ex_e + sum((i != j ? max_wb[i, j] - min_lmax[i, j] : 0.0) for i in 0:n, j in 1:n + 1))

    @constraint(m, con2[i in 0:n, j in 1:n + 1; i != j], x_cc[i, j] + x_b[i, j] + x_c[i, j] + x_e[i, j] == x[i, j])

    @constraint(m, con3[j in 1:n], sum((i != j ? x[i, j] : 0.0) for i in 0:n) == 1.0)
    @constraint(m, con4[i in 1:n], sum((i != j ? x[i, j] : 0.0) for j in 1:n + 1) == 1.0)

    @constraint(m, con5, sum(x[0, j] for j in 1:n) == 1.0)
    @constraint(m, con6, sum(x[i, n + 1] for i in 1:n) == 1.0)

    @constraint(m, con7, w[0] == 0)

    # Time flow constraints.
    @constraint(m, [i in 0:n, j in 1:n + 1; i != j], max_wb[i, j] >= w[i]) # max(wi, bi)
    @constraint(m, [i in 0:n, j in 1:n + 1; i != j], max_wb[i, j] >= d.tw_b[i]) # max(wi, bi)
    # @expression(m, tf_ex[i in 0:n, j in 1:n + 1; i != j], max_wb[i, j] + d.service[i] - w[j])

    @constraint(m, con8a[i in 0:n, j in 1:n + 1; i != j], max_wb[i, j] + d.service[i] + d.times["cc"][i, j] - w[j] <= (1 - x_cc[i, j]) * M)
    @constraint(m, con8b[i in 0:n, j in 1:n + 1; i != j], max_wb[i, j] + d.service[i] + d.times["cc"][i, j] - w[j] >= -(1 - x_cc[i, j]) * M)

    @constraint(m, con9a[i in 0:n, j in 1:n + 1; i != j], max_wb[i, j] + d.service[i] + d.times["b"][i, j] - w[j] <= (1 - x_b[i, j]) * M)
    @constraint(m, con9b[i in 0:n, j in 1:n + 1; i != j], max_wb[i, j] + d.service[i] + d.times["b"][i, j] - w[j] >= -(1 - x_b[i, j]) * M)

    @constraint(m, con10a[i in 0:n, j in 1:n + 1; i != j], max_wb[i, j] + d.service[i] + d.times["c"][i, j] - w[j] <= (1 - x_c[i, j]) * M)
    @constraint(m, con10b[i in 0:n, j in 1:n + 1; i != j], max_wb[i, j] + d.service[i] + d.times["c"][i, j] - w[j] >= -(1 - x_c[i, j]) * M)

    @constraint(m, con11a[i in 0:n, j in 1:n + 1; i != j], max_wb[i, j] + d.service[i] + d.times["e"][i, j] - w[j] <= (1 - x_e[i, j]) * M)
    @constraint(m, con11b[i in 0:n, j in 1:n + 1; i != j], max_wb[i, j] + d.service[i] + d.times["e"][i, j] - w[j] >= -(1 - x_e[i, j]) * M)
    
    @constraint(m, con12[j in 1:n + 1], w[j] + d.service[j] <= d.tw_e[j])
    @constraint(m, con13, l[0] == 0)
        
    # Battery charging constraints.
    @constraint(m, [i in 0:n, j in 1:n + 1; i != j], min_lmax[i, j] <= l[i] + d.times["cc"][i, j] * d.rate_c) # min(li + tcc_ij*rc, lmax).
    @constraint(m, [i in 0:n, j in 1:n + 1; i != j], min_lmax[i, j] <= d.l_max) # min(li + tcc_ij*rc, lmax)

    @constraint(m, con14a[i in 0:n, j in 1:n + 1; i != j], min_lmax[i, j] - l[j] <= (1 - x_cc[i, j]) * M)
    @constraint(m, con14b[i in 0:n, j in 1:n + 1; i != j], min_lmax[i, j] - l[j] >= -(1 - x_cc[i, j]) * M)

    @constraint(m, con15a[i in 0:n, j in 1:n + 1; i != j], l[i] - d.times["b"][i, j] * d.rate_d - l[j] <= (1 - x_b[i, j]) * M)
    @constraint(m, con15b[i in 0:n, j in 1:n + 1; i != j], l[i] - d.times["b"][i, j] * d.rate_d - l[j] >= -(1 - x_b[i, j]) * M)

    @constraint(m, con16a[i in 0:n, j in 1:n + 1; i != j], l[i] - l[j] <= (1 - x_c[i, j]) * M)
    @constraint(m, con16b[i in 0:n, j in 1:n + 1; i != j], l[i] - l[j] >= -(1 - x_c[i, j]) * M)

    @constraint(m, con17a[i in 0:n, j in 1:n + 1; i != j], l[i] - d.times["e"][i, j] * d.rate_d - l[j] <= (1 - x_e[i, j]) * M)
    @constraint(m, con17b[i in 0:n, j in 1:n + 1; i != j], l[i] - d.times["e"][i, j] * d.rate_d - l[j] >= -(1 - x_e[i, j]) * M)

    @constraint(m, con18[i in 0:n + 1], 0 <= l[i] <= d.l_max)
    
    # print(m)
    optimize!(m)
    println(termination_status(m))
    if termination_status(m) == MOI.INFEASIBLE
        println("Infeasible")
    else
        # solution_summary(m)
	i = 0
	print(0, " ")
	while true
	    for j in 1:n + 1
	        if i != j && diffValues(value(x[i, j]), 0.0)
		    if j != n + 1
		    	print(j, " ")
		    else
			println(0)
		    end
		    i = j
	        end
	        if i == n + 1
		    break
		end
	    end       
	    if i == n + 1
		break
	    end
	end
        println("Objective: ", objective_value(m))
        println("Objective: ", objective_value(m) + sum((i != j ? value(max_wb[i, j]) - value(min_lmax[i, j]) : 0.0) for i in 0:n, j in 1:n + 1))
    end

    return m
end

function diffValues(a, b)
    return a - b > 0.000001
end
