"""
main_moo.jl

Main entry point for running multi-objective optimization (MOO) on the model.
Handles command-line arguments, parameter loading, optimization, and result saving.
"""

# =======================
# Environment Setup
# =======================
using Pkg
Pkg.activate(".")                # Activate the local Julia environment
Pkg.instantiate()                 # Install missing packages if needed

# =======================
# Include Dependencies
# =======================
include("utils/utils.jl")         # Utility functions (argument parsing, etc.)
include("utils/save_results.jl")  # Functions to save optimization results
include("models/parameter_loader.jl")  # Parameter loading logic
include("moo_julia/optimization.jl")   # MOO algorithm implementation

# =======================
# Parse Command-Line Arguments
# =======================
parsed_args = parse_command_line()  # Custom function to parse CLI args

# Extract and convert μ (mu) argument to per-day units
μ_received = parsed_args["mu"]
μ = μ_received / 24  # Convert μ from per-hour to per-day units

# =======================
# Circadian Rhythm Option Parsing
# =======================
# Ensure mutually exclusive circadian flags
if parsed_args["circ"] && parsed_args["no-circ"]
    error("Cannot specify both --circ and --no-circ at the same time.")
end

# Determine if circadian rhythm should be enabled (default: true)
is_circadian = if parsed_args["circ"]
    true
elseif parsed_args["no-circ"]
    false
else
    true  # Default: circadian rhythm is ON
end

# =======================
# Load Model Parameters
# =======================
# Load parameter values from CSV and set initial conditions
path_para = joinpath("..", "..", "res", "parameter_values.csv")
p = load_parameter(path_para, μ)

# =======================
# Multi-Objective Optimization Setup
# =======================
E₂_upper = 30      # Upper bound for E2 variable
P₄_upper = 2000    # Upper bound for P4 variable
optimization_options = Dict(
    :x_tol => 1e-4,                  # Tolerance for solution variables
    :f_tol => 1e-4,                  # Tolerance for objective function
    :verbose => true,                # Print optimization progress
    :iterations => 10000,            # Maximum number of iterations
    :population_size => 400,         # Population size for evolutionary algorithm
    :seed => 1,                      # Random seed for reproducibility
    :parallel_evaluation => true,    # Enable parallel evaluation
    :bounds => [zeros(42) vcat(E₂_upper * ones(21), P₄_upper * ones(21))] # Variable bounds
)

# Run the multi-objective optimization
res, last_plot = run_moo_optimization(optimization_options, p, is_circadian)

# =======================
# Save and Visualize Results
# =======================
# Save optimization results and plot
save_optimization_results(res, μ_received, last_plot, p, is_circadian)