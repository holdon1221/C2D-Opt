using Pkg
Pkg.activate(".")
Pkg.instantiate()

include("utils/utils.jl")
include("utils/save_results.jl")

include("models/parameter_loader.jl")
include("moo_julia/optimization.jl")

# Command line μ
parsed_args = parse_command_line()

μ_received = parsed_args["mu"]
μ = μ_received / 24 # to day unit

# Command line for circadian rhythm
if parsed_args["circ"] && parsed_args["no-circ"]
    error("Cannot specify both --circ and --no-circ at the same time.")
end

is_circadian = if parsed_args["circ"]
    true
elseif parsed_args["no-circ"]
    false
else
    true 
end

# Load parameter values and set the initial condition, tspan and lag time
path_para = joinpath("..", "..", "res", "parameter_values.csv")
p = load_parameter(path_para, μ)

# Multi-objective optimization
E₂_upper = 30
P₄_upper = 2000
optimization_options = Dict(
    :x_tol => 1e-4,
    :f_tol => 1e-4,
    :verbose => true,
    :iterations => 10000,
    :population_size => 400,
    :seed => 1,
    :parallel_evaluation => true,
    :bounds => [zeros(42) vcat(E₂_upper * ones(21), P₄_upper * ones(21))]
)
res, last_plot = run_moo_optimization(optimization_options, p, is_circadian)

# Visualize and save the result
save_optimization_results(res, μ_received, last_plot, p, is_circadian)