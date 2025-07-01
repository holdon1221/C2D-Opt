"""
save_results.jl

This module provides functions to save the results of the optimization process,
including the maximum P₄ value for a given number of cycles, optimization results,
and dynamics and control bars.

# Functions
- `save_maximum_P4_given_cycles!`: Save the maximum P₄ value for a given number of cycles.
- `save_optimization_results`: Save the optimization results.
- `save_dynamics_and_control_bars`: Save dynamics and control bars.
"""

include("utils.jl")
include("../models/dde_model.jl")

using Dates
using CodecZlib
using Plots
using Serialization

# ============================================================================
# Save Maximum P₄ Given Cycles
# ============================================================================
"""
    save_maximum_P4_given_cycles!(res::State, num_cycle::Float64, p::Vector{Any}, is_circadian::Bool)

Save the maximum P₄ value for a given number of cycles.

# Arguments
- `res::State`: Result state
- `num_cycle::Float64`: Number of cycles for finding maximum P₄
- `p::Vector{Any}`: Parameter vector
- `is_circadian::Bool`: Boolean parameter for considering circadian rhythm

# Key Steps
1. Calculate the candidates for maximum P₄.
2. Iterate over the population and calculate the maximum P₄ for each candidate.
3. Update the result state with the maximum P₄ values.
"""
function save_maximum_P4_given_cycles!(res::State, num_cycle::Float64, p::Vector{Any}, is_circadian::Bool)
    # Calculate the candidates for maximum P₄
    println("Calculating the candidates... Num sol: $(length(res.population))")
    for (i, r) in enumerate(res.population)
        # Copy the parameter value to independently manipulate it
        local_p = copy(p)

        # Calculate the taking dates and breaking dates
        taking_dates = Int64(length(r.x) / 2)
        dates_per_cycle = 28
        breaking_dates = dates_per_cycle - taking_dates

        # Create the control vectors
        u1 = vcat(r.x[1:taking_dates], zeros(breaking_dates))
        u2 = vcat(r.x[(taking_dates+1):end], zeros(breaking_dates))

        # Update the local parameter with the control vectors
        local_p[end-2] = u1
        local_p[end-1] = u2

        # Solve the DDE and calculate the maximum P₄
        sol, _, _ = get_solution(r, p, num_cycle, is_circadian)
        if is_circadian
            _, P₄, _ = calculate_auxiliary(sol, local_p, sol.t)
        else
            _, P₄, _ = calculate_auxiliary_no_circ(sol, local_p, sol.t)
        end

        # Update the result state with the maximum P₄ value
        push!(r.f, maximum(P₄))
    end

    return res
end

# ============================================================================
# Save Optimization Results
# ============================================================================
"""
    save_optimization_results(res::State, μ_received::Float64, last_plot, p::Vector{Any}, is_circadian::Bool)

Save the optimization results.

# Arguments
- `res::State`: Result state
- `μ_received::Float64`: The administration time.
- `last_plot::Plots.Plot`: The last plot drawn during optimization process (nothing possible)
- `p::Vector{Any}`: Parameter values vector
- `is_circadian::Bool`: Boolean parameter for considering circadian rhythm

# Key Steps
1. Create the filename and directory for saving the results.
2. Save the last plot if it exists.
3. Save the maximum P₄ value for 10 cycles.
4. Save the result variable.
5. Select the feasible population and save the dynamics and control bars.
"""
function save_optimization_results(res::State, μ_received::Float64, last_plot, p::Vector{Any}, is_circadian::Bool)
    # Create the filename and directory for saving the results
    date = Dates.format(now(), "yyyymmdd")
    generation = res.iteration  # The last generation
    population = length(res.population) # Population size
    time = μ_received # Taking pill

    filename = "$(date)_Clock$(time)_CCMO{NSGA2}_Gen$(generation)_Pop$(population)" * (is_circadian ? "" : "_NoCirc")
    filename = joinpath("..", "..", "res", filename)
    isdir(joinpath("..", "..", "res")) || mkdir(joinpath("..", "..", "res"))
    isdir(filename) || mkdir(filename)

    # Save the last plot if it exists
    if last_plot !== nothing
        savefig(last_plot, joinpath(filename, "pareto_front.png"))
    end

    # Save the maximum P₄ value for 10 cycles
    save_maximum_P4_given_cycles!(res, 10.0, p, is_circadian)

    # Save the result variable
    serialize(joinpath(filename, "result_variable.jls"), res)

    # Select the feasible population and save the dynamics and control bars
    smallest_five = select_states(res, 3.0, 5, true)
    for (i, candidate) in enumerate(smallest_five)
        save_dynamics_and_control_bars(candidate, p, joinpath(filename, "f$(i)_"), is_circadian)
    end
end

# ============================================================================
# Save Dynamics and Control Bars
# ============================================================================
"""
    save_dynamics_and_control_bars(candidate::Metaheuristics.xFgh_solution, p::Vector{Any}, filename::String, is_circadian::Bool)

Save dynamics and control bars.

# Arguments
- `candidate::Metaheuristics.xFgh_solution`: Drawn result candidate
- `p::Vector{Any}`: Parameter values vector
- `filename::String`: File path
- `is_circadian::Bool`: Boolean parameter for considering circadian rhythm

# Key Steps
1. Copy the parameter value to independently manipulate it.
2. Solve the DDE and calculate the P₄ values.
3. Create the dynamics and control bar plots.
4. Save the plots.
"""
function save_dynamics_and_control_bars(candidate::Metaheuristics.xFgh_solution, p::Vector{Any}, filename::String, is_circadian::Bool)
    # Copy the parameter value to independently manipulate it
    local_p = copy(p)

    # Solve the DDE and calculate the P₄ values
    sol, u1, u2 = get_solution(candidate, local_p, 5.0, is_circadian) # drawing 5 cycles
    if is_circadian
        _, P₄, _ = calculate_auxiliary(sol, local_p, sol.t)
    else
        _, P₄, _ = calculate_auxiliary_no_circ(sol, local_p, sol.t)
    end
    tspan = (sol.t[1], sol.t[end])

    # Create the dynamics and control bar plots
    expressed_amount_u = round(candidate.f[2], digits=1)

    dynamics_FSH = plot(sol, vars=(0, 4), xlims=tspan, xlabel="Time (days)", ylabel="FSH [IU/L]", 
        legend=false, linewidth=4, lc = :dodgerblue,
        y_guidefontcolor=:dodgerblue, y_foreground_color_text=:dodgerblue, 
        y_foreground_color_axis=:dodgerblue, y_foreground_color_border=:dodgerblue, 
        guidefont=15, tickfont=15)
    savefig(dynamics_FSH, string(filename, "sum$(expressed_amount_u)_dynamics_FSH.png"))

    dynamics_P4 = plot(sol.t, P₄, xlims=tspan, ylabel="P₄ [ng/ml]", 
        legend=false, linewidth=4, lc=:crimson,  
        y_guidefontcolor=:crimson, y_foreground_color_text=:crimson, 
        y_foreground_color_axis=:crimson, y_foreground_color_border=:crimson, 
        guidefont=15, tickfont=15)
    savefig(dynamics_P4, string(filename, "sum$(expressed_amount_u)_dynamics_P4.png"))

    # Create the control bar plots
    days = 0:27
    b1 = bar(days, u1, xlabel="Time [days]", ylabel="Exo E₂ [μg]", legend=false, 
        bar_width=0.7, xticks=days, color=:dodgerblue)

    # Add control info to the plot
    total_amount_u1 = sum(u1)
    formatted_total_amount_u1 = round(total_amount_u1, digits=2)
    dose_count_u1 = sum(u1 .!= 0)

    annotate!(maximum(days)*1, maximum(u1)*0.95, text("Total Amount: $formatted_total_amount_u1", 10, :right))
    annotate!(maximum(days)*1, maximum(u1)*0.85, text("Dose Count: $dose_count_u1", 10, :right))

    savefig(b1, string(filename, "sum$(expressed_amount_u)_exo_E2.png"))

    b2 = bar(days, u2, xlabel="Time [days]", ylabel="Exo P₄ [μg]", legend=false, 
        bar_width=0.7, xticks=days, color=:crimson)

    # Add control info to the plot
    total_amount_u2 = sum(u2)
    formatted_total_amount_u2 = round(total_amount_u2, digits=2)
    dose_count_u2 = sum(u2 .!= 0)

    annotate!(maximum(days)*1, maximum(u2) *0.95, text("Total Amount: $formatted_total_amount_u2", 10, :right))
    annotate!(maximum(days)*1, maximum(u2) *0.85, text("Dose Count: $dose_count_u2", 10, :right))

    savefig(b2, string(filename, "sum$(expressed_amount_u)_exo_P4.png"))
end