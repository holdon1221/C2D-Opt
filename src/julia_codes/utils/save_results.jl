include("utils.jl")
include("../models/dde_model.jl")

using Dates
using CodecZlib
using Plots
using Serialization

"""
    save_maximum_P4_given_cycles!(res::State, num_cycle::Float64, p::Vector{Any})

Insert the maximum value of P₄ for all given results.

# Arguments
- `res::State`: Result state
- `num_cycle::Float64`: Number of cycles for finding maximum P₄
- `p::Vector{Any}`: Parameter vector
- `is_circadian:Bool`: Boolean parameter for considering circadian rhythm
"""
function save_maximum_P4_given_cycles!(res::State, num_cycle::Float64, p::Vector{Any}, is_circadian::Bool)
    println("Calculating the candidates... Num sol: $(length(res.population))")
    for (i, r) in enumerate(res.population)
        local_p = copy(p)

        taking_dates = Int64(length(r.x) / 2)
        dates_per_cycle = 28
        breaking_dates = dates_per_cycle - taking_dates

        u1 = vcat(r.x[1:taking_dates], zeros(breaking_dates))
        u2 = vcat(r.x[(taking_dates+1):end], zeros(breaking_dates))

        local_p[end-2] = u1
        local_p[end-1] = u2

        sol, _, _ = get_solution(r, p, num_cycle, is_circadian)
        if is_circadian
            _, P₄, _ = calculate_auxiliary(sol, local_p, sol.t)
        else
            _, P₄, _ = calculate_auxiliary_no_circ(sol, local_p, sol.t)
        end

        push!(r.f, maximum(P₄))
    end

    return res
end

"""
    save_optimization_results()

Function to save the restults in file and visual graph

# Arguments
- `res::State`: Result state
- `μ_received::Float64`: The administration time.
- `last_plot::Plots.Plot`: The last plot drawn during optimization process (nothing possible)
- `p::Vector{Any}`: Parameter values vector
- `is_circadian:Bool`: Boolean parameter for considering circadian rhythm
"""
function save_optimization_results(res::State, μ_received::Float64, last_plot, p::Vector{Any}, is_circadian::Bool)
    date = Dates.format(now(), "yyyymmdd")
    generation = res.iteration  # The last generation
    population = length(res.population) # Population size
    time = μ_received # Taking pill

    filename = "$(date)_Clock$(time)_CCMO{NSGA2}_Gen$(generation)_Pop$(population)" * (is_circadian ? "" : "_NoCirc")
    filename = joinpath("..", "..", "res", filename)
    isdir(joinpath("..", "..", "res")) || mkdir(joinpath("..", "..", "res"))
    isdir(filename) || mkdir(filename)

    if last_plot !== nothing
        savefig(last_plot, joinpath(filename, "pareto_front.png"))
    end

    # Save the maximum P₄ value for 10 cycles
    save_maximum_P4_given_cycles!(res, 10.0, p, is_circadian)

    # Save the result variable.
    serialize(joinpath(filename, "result_variable.jls"), res)

    # Select feasible population
    smallest_five = select_states(res, 3.0, 5, true)

    # Save the result dynamics and control variables
    for (i, candidate) in enumerate(smallest_five)
        save_dynamics_and_control_bars(candidate, p, joinpath(filename, "f$(i)_"), is_circadian)
    end
end

"""
    save_dynamics_and_control_bars(candidate::Metaheuristics.xFgh_solution)

Save dynamics and control bar figures

# Arguments
- `candidate::Metaheuristics.xFgh_solution`: Drawn result candidate
- `p::Vector{Any}`: Parameter values vector
- `filename::String`: File path
- `is_circadian:Bool`: Boolean parameter for considering circadian rhythm
"""
function save_dynamics_and_control_bars(candidate::Metaheuristics.xFgh_solution, p::Vector{Any}, filename::String, is_circadian::Bool)
    # Copy the parameter value in order to independently manipulate it
    local_p = copy(p)

    # Draw dynamics and the optimal controls
    days = 0:27

    # Solve and get P₄ values
    sol, u1, u2 = get_solution(candidate, local_p, 5.0, is_circadian) # drawing 5 cycles
    if is_circadian
        _, P₄, _ = calculate_auxiliary(sol, local_p, sol.t)
    else
        _, P₄, _ = calculate_auxiliary_no_circ(sol, local_p, sol.t)
    end
    tspan = (sol.t[1], sol.t[end])
    
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

    # Optimal control bar
    b1 = bar(days, u1, xlabel="Time [days]", ylabel="Exo E₂ [μg]", legend=false, 
        bar_width=0.7, xticks=days, color=:dodgerblue)

    # Optimal control info
    total_amount_u1 = sum(u1)
    formatted_total_amount_u1 = round(total_amount_u1, digits=2)
    dose_count_u1 = sum(u1 .!= 0)

    annotate!(maximum(days)*1, maximum(u1)*0.95, text("Total Amount: $formatted_total_amount_u1", 10, :right))
    annotate!(maximum(days)*1, maximum(u1)*0.85, text("Dose Count: $dose_count_u1", 10, :right))

    savefig(b1, string(filename, "sum$(expressed_amount_u)_exo_E2.png"))

    # Optimal control bar
    b2 = bar(days, u2, xlabel="Time [days]", ylabel="Exo P₄ [μg]", legend=false, 
        bar_width=0.7, xticks=days, color=:crimson)

    # Optimal control info
    total_amount_u2 = sum(u2)
    formatted_total_amount_u2 = round(total_amount_u2, digits=2)
    dose_count_u2 = sum(u2 .!= 0)

    annotate!(maximum(days)*1, maximum(u2) *0.95, text("Total Amount: $formatted_total_amount_u2", 10, :right))
    annotate!(maximum(days)*1, maximum(u2) *0.85, text("Dose Count: $dose_count_u2", 10, :right))

    savefig(b2, string(filename, "sum$(expressed_amount_u)_exo_P4.png"))
end