include("cost_function.jl")
include("../visualization/plot_results.jl")

using Metaheuristics
using Plots
using Random

"""
    run_moo_optimization(option::Dict, p::Vector{Any})

Calculate the result using Multi-objective optimization function.

# Arguments
- `option::Dict`: Options for the optimization algorithm. The type is from Metaheuristics.jl
- `p::Vector{Any}`: Parameter values vector
- `is_circadian:Bool`: Boolean parameter for considering circadian rhythm
"""
function run_moo_optimization(option::Dict, p::Vector{Any}, is_circadian::Bool)
    # Get the options
    x_tol = get(option, :x_tol, 1e-4)
    f_tol = get(option, :f_tol, 1e-4)
    verbose = get(option, :verbose, true)
    iterations = get(option, :iterations, 800)
    population_size = get(option, :population_size, 400)
    seed = get(option, :seed, 1)
    parallel_evaluation = get(option, :parallel_evaluation, false)
    bounds = get(option, :bounds, [zeros(42), ones(42)])
    store_convergence = get(option, :store_convergence, false)

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
    # Default case
    opt = CCMO(NSGA2(N=population_size), options=options)
    # # hyp case
    # opt = CCMO(NSGA2(N=population_size, η_cr=10, η_m=10, p_m=0.1), options=options)

    # Set values for visualizations
    algorithm_name = string(typeof(opt.parameters))
    last_plot = nothing

    logger(st::State) = begin
        A = fvals(st)

        smallest_states = select_states(st, 3.0, 1, false)[1]

        u1 = vcat(smallest_states.x[1:21], zeros(7))
        u2 = vcat(smallest_states.x[22:end], zeros(7))

        last_plot = logger_drawer(A, u1, u2, algorithm_name, st.iteration)
    end

    if parallel_evaluation
        res = optimize(X -> f_parallel(X, p, is_circadian), bounds, opt, logger=logger)
    else
        res = optimize(X -> moo_cost_function(X, p, is_circadian), bounds, opt, logger=logger)
    end
    return res, last_plot
end
