mutable struct Data
    ncustomers::Int
    service::Dict{Int, Float64} # Service time.
    # Time Windows.
    tw_b::Dict{Int, Float64}
    tw_e::Dict{Int, Float64}
    costs::Dict{String, Dict{Tuple{Int, Int}, Float64}} # Costs for each mode of operation.
    times::Dict{String, Dict{Tuple{Int, Int}, Float64}} # Travel times for each mode of operation.
    l_init::Float64 # Battery initial charging.
    l_max::Float64 # Maximum battery capacity.
    rate_c::Float64 # Charging rate.
    rate_d::Float64 # Discharging rate.

    Data() = new(0, Dict{Int, Float64}(), Dict{Int, Float64}(), Dict{Int, Float64}(), 
                 Dict{String, Dict{Tuple{Int, Int}, Float64}}(), 
                 Dict{String, Dict{Tuple{Int, Int}, Float64}}(), 0, 0, 0, 0)
end

function read_data(filename::String)
    file = readlines(filename)
    d = Data()
    d.ncustomers = parse(Int, file[5])
    nvertices = d.ncustomers + 1

    # Service time and time windows.
    s = map(x -> parse(Float64, x), split(file[8], ", "))
    tw = map(x -> parse(Float64, x), split(file[35], ", "))
    tw_b = []
    tw_e = []
    i = 1
    while i < length(tw)
        push!(tw_b, tw[i])
        push!(tw_e, tw[i + 1])
        i += 2
    end
    push!(tw_b, tw[1])
    push!(tw_e, tw[2])

    for i in 0:nvertices
        d.service[i] = s[i + 1]
        d.tw_b[i] = tw_b[i + 1]
        d.tw_e[i] = tw_e[i + 1]
    end

    # Battery data.
    d.l_init = parse(Float64, file[11])
    d.l_max = parse(Float64, file[14])
    d.rate_c = parse(Float64, file[17]) / 60.0
    d.rate_d = parse(Float64, file[20]) / 60.0

    # Costs and travel times.
    modes = ["c", "cc", "e", "b"]
    t = map(x -> parse(Float64, x), split(file[29], ", "))
    c = map(x -> parse(Float64, x), split(file[32], ", "))
    times = Dict{String, Matrix{Float64}}()
    costs = Dict{String, Matrix{Float64}}()
    i = 1
    for (_, md) in enumerate(modes)
        times[md] = Matrix{Float64}(undef, nvertices, nvertices)
        costs[md] = Matrix{Float64}(undef, nvertices, nvertices)
        for j in 1:nvertices, k in 1:nvertices
            # Put the depot as the last column. For instance: 0 1 2 3 4 5 6 7 8 => 1 2 3 4 5 6 7 8 0.
            k_aux = (k == 1 ? nvertices : k - 1)
            times[md][j, k_aux] = t[i]
            costs[md][j, k_aux] = c[i]
            i += 1
        end
    end

    for (_, md) in enumerate(modes)
        d.costs[md] = Dict{Tuple{Int, Int}, Float64}()
        d.times[md] = Dict{Tuple{Int, Int}, Float64}()
        for i in 0:d.ncustomers, j in 1:d.ncustomers + 1
            d.costs[md][i, j] = costs[md][i + 1, j]
            d.times[md][i, j] = times[md][i + 1, j]
        end
    end

    return d
end
