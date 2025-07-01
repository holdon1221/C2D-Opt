"""
utils.jl

Utility functions for command-line argument parsing and result selection in the optimization workflow.
Provides:
- Command-line interface parsing for experiment configuration
- Selection of feasible and optimal states from optimization results
"""

using ArgParse
using Metaheuristics

include("../models/dde_model.jl")

# =======================
# Command-Line Argument Parsing
# =======================
"""
    parse_command_line()

Parse command-line arguments for the optimization script.

# Returns
- `Dict{String, Any}`: Dictionary containing parsed values from the command line
"""
function parse_command_line()
    s = ArgParseSettings()
    # Define command-line arguments
    @add_arg_table s begin
        "--mu", "-m"  # Administration timing argument
        help = "Value for administration timing"
        arg_type = Float64
        default = 9.0

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
    # Parse and return arguments as a dictionary
    return parse_args(ARGS, s)
end

# =======================
# State Selection Utility
# =======================
"""
    select_states(st::State, bound::Float64, num_vis::Int64, is_several_cycle::Bool)

Select states with minimal objective 2 (e.g., drug usage) among those with objective 1 below a bound.
Handles both single- and multi-cycle selection modes.

# Arguments
- `st::State`: State variable from Metaheuristics.jl containing optimization population
- `bound::Float64`: Upper bound for feasible set in terms of objective 1 (e.g., constraint)
- `num_vis::Int64`: Number of solutions to select
- `is_several_cycle::Bool`: Whether to use multi-cycle (true) or single-cycle (false) selection

# Returns
- `smallest_st::Vector{Metaheuristics.xFgh_solution{Vector{Float64}}}`: The selected solutions
"""
function select_states(st::State, bound::Float64, num_vis::Int64, is_several_cycle::Bool)
    # Filter population by constraint on P₄ (or other objective)
    if is_several_cycle
        filtered_st = filter(x -> x.f[3] < bound, st.population)
    else
        filtered_st = filter(x -> x.f[1] < bound, st.population)
    end
    # Sort by objective 2 (e.g., drug usage)
    sorted_st_pop = sort(filtered_st, by = x -> x.f[2])
    # If no feasible population, select by minimal usage regardless of constraint
    isempty(sorted_st_pop) ? sorted_st_pop = sort(st.population, by = x -> x.f[2]) : nothing
    # Pick the top num_vis results
    smallest_st = sorted_st_pop[1:num_vis]
    return smallest_st
end