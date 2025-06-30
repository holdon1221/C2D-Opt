"""
cost_function.jl

Defines the cost functions and parallel evaluation utilities for multi-objective optimization
of the menstrual cycle model. Includes both single and parallel cost function evaluation,
and helper functions for mathematical operations.
"""

include("../models/dde_model.jl")

using DifferentialEquations
using Integrals

# =======================
# Main Multi-Objective Cost Function
# =======================
"""
    moo_cost_function(u::Vector{Float64}, local_p::Vector{Any}, is_circadian::Bool)

Compute the objective and constraint values for a candidate control vector in MOO.

# Arguments
- `u::Vector{Float64}`: Control variable (concatenated E2 and P4 dosing)
- `local_p::Vector{Any}`: Parameter values (copied for thread safety)
- `is_circadian::Bool`: Whether to include circadian rhythm in the model

# Returns
- `fu::Vector{Float64}`: Objective function values ([max P4, total dosing])
- `gu::Vector{Float64}`: Inequality constraint values
- `hu::Vector{Float64}`: Equality constraint values
"""
function moo_cost_function(u::Vector{Float64}, local_p::Vector{Any}, is_circadian::Bool)
    # Split control vector into E2 (u1) and P4 (u2) with 28 days (21+7 zeros)
    u1 = vcat(u[1:21], zeros(7))
    u2 = vcat(u[22:end], zeros(7))

    # Insert control vectors into parameter array for this evaluation
    local_p[end-2] = u1
    local_p[end-1] = u2

    # Set up initial conditions and DDE problem
    x0, tspan, lags = setup_initial_conditions_time_info(local_p)
    if is_circadian
        prob = DDEProblem(menstrual_cycle_dynamics!, x0, history_func, tspan, local_p; constant_lags=lags)
    else
        prob = DDEProblem(menstrual_cycle_dynamics_no_circ!, x0, history_func, tspan, local_p; constant_lags=lags)
    end

    # Solve the DDE
    sol = solve(prob, MethodOfSteps(RK4()), dtmax=0.05)
    # Compute P4 time course using auxiliary function
    if is_circadian
        _, P₄, _ = calculate_auxiliary(sol, local_p, sol.t)
    else
        _, P₄, _ = calculate_auxiliary_no_circ(sol, local_p, sol.t)
    end

    # Main objective: maximum P4 concentration over 3 cycles
    obj1 = maximum(P₄)
    # Secondary objective: total amount of dosing (E2 + P4)
    obj2 = sum(u1)+sum(u2)

    fu = [obj1, obj2]
    gu = [0.0]  # No inequality constraints (placeholder)
    hu = [0.0]  # No equality constraints (placeholder)

    return fu, gu, hu
end

# =======================
# Parallel Cost Function Evaluation
# =======================
"""
    f_parallel(X::Matrix{Float64}, p::Vector{Any}, is_circadian::Bool)

Evaluate the cost function in parallel for a batch of candidate solutions.
Each row of X is a candidate control vector.

# Arguments
- `X::Matrix{Float64}`: Input matrix, each row is a candidate
- `p::Vector{Any}`: Parameter values
- `is_circadian::Bool`: Whether to include circadian rhythm

# Returns
- `fx::Matrix{Float64}`: Objective values for each candidate
- `gx::Matrix{Float64}`: Inequality constraint values
- `hx::Matrix{Float64}`: Equality constraint values
"""
function f_parallel(X::Matrix{Float64}, p::Vector{Any}, is_circadian::Bool)
    N = size(X, 1)
    nobjectives = 2
    fx, gx, hx = zeros(N, nobjectives), zeros(N, 1), zeros(N, 1)
    errors = Threads.Atomic{Int}(0)

    # Launch a thread for each candidate
    handles = Vector{Task}(undef, N)
    for i in 1:N
        handles[i] = Threads.@spawn begin
            try
                local_p = copy(p)
                fx[i, :], gx[i, :], hx[i, :] = moo_cost_function(X[i, :], local_p, is_circadian)
            catch e
                println("Error occurred at X[$i,:] = ", X[i, 1:10])
                println("Error: ", e)
                fx[i, :], gx[i, :], hx[i, :] = NaN, NaN, NaN 
                Threads.atomic_add!(errors, 1)
            end
        end
    end

    # Wait for all threads to finish
    for handle in handles
        wait(handle)
    end

    return fx, gx, hx
end

# =======================
# Helper Functions
# =======================
"""
    σ(x::Float64)

Numerically stable sigmoid function (from NNlib.jl)

# Arguments
- `x::Float64`: Input value
# Returns
- `σ(x)`: Sigmoid of x
"""
function σ(x::Float64)
    t = exp(-abs(x))
    ifelse(x ≥ 0, inv(1 + t), t / (1 + t))
end