using ArgParse
using Metaheuristics

include("../models/dde_model.jl")

"""
    parse_command_line()

Function returning the values from command lines

# Returns
- `Dict{String, Any}`: Dictionary including the values form command line
"""
function parse_command_line()
    s = ArgParseSettings()
    
    @add_arg_table s begin
        "--mu", "-m"  # name
        help = "Value for administration timing"  # help
        arg_type = Float64  # typeof
        default = 9.0  # default value

        "--circ"
        help = "Enable circadian rhythm"
        arg_type = Bool
        default = false
        action = :store_true
        
        "--no-circ"
        help = "Disable circadian rhythm"
        arg_type = Bool
        default = false
        action = :store_true
    end

    return parse_args(ARGS, s)
end

"""
    select_states(st::State, bound::Float64, num_vis::Int64)

Select states which have a minimal amount (obj2) within the condition of obj1 less than bound.

# Arguments
- `st::State`: State varible during optimization algorithm in Metaheruistics.jl.
- `bound::Float64`: Upper bound for feasible set in terms of obj1
- `num_vis::Int64`: The selection number of smallest ones.

# Returns
- `smallest_st::Vector{Metaheuristics.xFgh_solution{Vector{Float64}}}`: The chosen ones
"""
function select_states(st::State, bound::Float64, num_vis::Int64, is_several_cycle::Bool)
    # Selection criteria for P₄ is less than `bound`
    if is_several_cycle
        filtered_st = filter(x -> x.f[3] < bound, st.population)
    else
        filtered_st = filter(x -> x.f[1] < bound, st.population)
    end
    sorted_st_pop = sort(filtered_st, by = x -> x.f[2])
    
    # If there is no population, select the minimum usage.
    isempty(sorted_st_pop) ? sorted_st_pop = sort(st.population, by = x -> x.f[2]) : nothing

    # pick the num_vis number results
    smallest_st = sorted_st_pop[1:num_vis]

    return smallest_st
end