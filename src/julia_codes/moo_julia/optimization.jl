"""
optimization.jl

Defines the main optimization routine for running multi-objective optimization (MOO)
on the menstrual cycle model. Handles algorithm configuration, logging, and visualization.
"""

include("cost_function.jl")
include("../visualization/plot_results.jl")

using Metaheuristics
using Plots
using Random


# =======================
# Main MOO Optimization Routine
# =======================
"""
    run_moo_optimization(option::Dict, p::Vector{Any}, is_circadian::Bool)

Run the multi-objective optimization using Metaheuristics.jl and the provided cost function.
Handles configuration, logging, and (optionally) parallel evaluation.

# Arguments
- `option::Dict`: Options for the optimization algorithm (e.g., tolerances, population size, bounds)
- `p::Vector{Any}`: Parameter values vector
- `is_circadian::Bool`: Whether to include circadian rhythm in the model

# Returns
- `res`: Optimization result (Metaheuristics.State)
- `last_plot`: Last plot generated during optimization (for visualization)
"""
function run_moo_optimization(option::Dict, p::Vector{Any}, is_circadian::Bool)
    # Extract options with defaults for the optimizer
    x_tol = get(option, :x_tol, 1e-4)
    f_tol = get(option, :f_tol, 1e-4)
    verbose = get(option, :verbose, true)
    iterations = get(option, :iterations, 800)
    population_size = get(option, :population_size, 400)
    seed = get(option, :seed, 1)
    parallel_evaluation = get(option, :parallel_evaluation, false)
    bounds = get(option, :bounds, [zeros(42), ones(42)])
    store_convergence = get(option, :store_convergence, false)

    # Set up optimization options and algorithm
    options = Options(
        x_tol=x_tol,
        f_tol=f_tol,
        verbose=verbose,
        iterations=iterations,
        f_calls_limit=Inf,
        parallel_evaluation=parallel_evaluation,
        store_convergence=store_convergence, 
        seed=seed,
        rng=Random.MersenneTwister(seed)
    )
    # Create the CCMO optimizer with NSGA2 backend
    opt = CCMO(NSGA2(N=population_size), options=options)
    # # Example for hyperparameter tuning (uncomment if needed)
    # opt = CCMO(NSGA2(N=population_size, η_cr=10, η_m=10, p_m=0.1), options=options)

    # Set up algorithm name and last_plot for logging
    algorithm_name = string(typeof(opt.parameters))
    last_plot = nothing

    # Logger function for visualization during optimization
    logger(st::State) = begin
        A = fvals(st)
        # Select the best state (lowest obj2 among those with obj1 < 3.0)
        smallest_states = select_states(st, 3.0, 1, false)[1]
        u1 = vcat(smallest_states.x[1:21], zeros(7))
        u2 = vcat(smallest_states.x[22:end], zeros(7))
        last_plot = logger_drawer(A, u1, u2, algorithm_name, st.iteration)
    end

    # Choose parallel or serial evaluation based on option
    if parallel_evaluation
        res = optimize(X -> f_parallel(X, p, is_circadian), bounds, opt, logger=logger)
    else
        res = optimize(X -> moo_cost_function(X, p, is_circadian), bounds, opt, logger=logger)
    end
    return res, last_plot
end
