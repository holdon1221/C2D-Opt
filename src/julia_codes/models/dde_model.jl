"""
dde_model.jl

Defines the delay differential equation (DDE) model for the menstrual cycle,
including initial condition setup, auxiliary calculations for hormone dynamics, and the DDE system definition.
All functions are documented and commented for clarity and maintainability.
"""

using DifferentialEquations

# =======================
# Initial Condition and Time Setup
# =======================
"""
    setup_initial_conditions_time_info(p)

Set up the initial conditions (history), simulation time span, and delay terms for the DDE model.

# Arguments
- `p::Vector{Any}`: Model parameter variables

# Returns
- `x0::Vector{Float64}`: Initial values (constant history)
- `tspan::Tuple{Float64, Float64}`: Simulation time span
- `lags::Vector{Float64}`: Delay terms
"""
function setup_initial_conditions_time_info(p::Vector{Any})
    # Extract initial state values from the parameter vector
    x0 = Float64[p[end-15:end-3]...]

    # Compute initial values for E₂, P₄, Inh using the auxiliary function
    E₂_init, P₄_init, Inh_init = calculate_auxiliary(x0, p, 0.0)
    push!(x0, E₂_init, P₄_init, Inh_init)

    # Simulate three cycles (84 days); adjust as needed
    tspan = (0.0, 3*28.0)
    # Set delay (lag) for DDE in days
    lags = [1.5]  # 1.5-day delay for hormone dynamics

    return x0, tspan, lags
end

# =======================
# Auxiliary Calculations
# =======================
"""
    calculate_auxiliary(x, p, t; n_cycles=nothing)

Calculate E2, P4, Inh values

# Arguments
- `x::Vector{Float64}`: State variables
- `p::Vector{Any}`: Model parameter values
- `t::Float64`: Current time
- `n_cycles::Union{Nothing, Int}`: Number of cycles to simulate; if `nothing`, simulate indefinitely

# Returns
- `(E2, P4, Inh)`: Calculated values
"""
function calculate_auxiliary(x, p, t; n_cycles=nothing)
    # Extract relevant state variables
    GrF, DomF, Lut₂, Lut₃, Lut₄ = map(i -> x[i,:], [6:7;11:13])
    
    # Extract model parameters
    e_0, e_1, e_2, e_3, p_1, p_2, h_0, h_1, h_2, h_3, θ_1, θ_2, θ_3, θ_4, _, _, _, _, η_E₂, F_E₂, k_a_E₂, k_21_E₂, alpha1_E₂, beta1_E₂, η_P₄, F_P₄, k_a_P₄, k_21_P₄, alpha1_P₄, beta1_P₄, V_c_E₂, V_c_P₄, u1, u2, μ = p[[31:62; end-2:end]]
    
    # Control
    today = floor.(Int, t)
    exo_E₂ = zeros(length(t))
    exo_P₄ = zeros(length(t))

    # Code optimization
    E₂_constant = k_a_E₂ * F_E₂ / V_c_E₂
    N_E₂_constant = (k_21_E₂ - k_a_E₂) / ( (alpha1_E₂ - k_a_E₂) * (beta1_E₂ - k_a_E₂) )
    L_E₂_constant = (k_21_E₂ - alpha1_E₂) / ( (k_a_E₂ - alpha1_E₂) * (beta1_E₂ - alpha1_E₂) )
    M_E₂_constant = (k_21_E₂ - beta1_E₂) / ( (k_a_E₂ - beta1_E₂) * (alpha1_E₂ - beta1_E₂) )
    P₄_constant = k_a_P₄ * F_P₄ / V_c_P₄
    N_P₄_constant = (k_21_P₄ - k_a_P₄) / ( (alpha1_P₄ - k_a_P₄) * (beta1_P₄ - k_a_P₄) )
    L_P₄_constant = (k_21_P₄ - alpha1_P₄) / ( (k_a_P₄ - alpha1_P₄)*(beta1_P₄ - alpha1_P₄) )
    M_P₄_constant = (k_21_P₄ - beta1_P₄) / ( (k_a_P₄ - beta1_P₄)*(alpha1_P₄ - beta1_P₄) )
    current_pill_time = today .+ μ

    # Calculate exogenous hormone contributions
    for i in eachindex(t)
        current_cycle = div(today[i], 28) + 1 # Current cycle number

        if today[i] > 0 # From the 1st cycle (simulation starts from day 0)
            past_days = 0:(today[i] - 1) # past days
            mod_past_days = mod.(past_days, 28) # dosing strategy have a cycle of 28 days
            i_pill_taking_time = past_days .+ μ # vector of taking time
            passed_time_from_pills = t[i] .- i_pill_taking_time # vector of passed time from each administration

            # Past pills' contribution
            exo_E₂[i] = sum(u1[mod_past_days .+ 1] .* E₂_constant .* (N_E₂_constant .* exp.(-k_a_E₂ .* passed_time_from_pills) .+ L_E₂_constant .* exp.(-alpha1_E₂ .* passed_time_from_pills) .+ M_E₂_constant .* exp.(-beta1_E₂ .* passed_time_from_pills))) 

            exo_P₄[i] = sum(u2[mod_past_days .+ 1] .* P₄_constant .* (N_P₄_constant .* exp.(-k_a_P₄ .* passed_time_from_pills) .+ L_P₄_constant .* exp.(-alpha1_P₄ .* passed_time_from_pills) .+ M_P₄_constant .* exp.(-beta1_P₄ .* passed_time_from_pills)))
        end

        # Current pill's contribution
        if (isnothing(n_cycles) || current_cycle <= n_cycles) && (t[i] - today[i] >= μ)
            mod_today = mod(today[i], 28)
            passed_time_from_current_pill = t[i] - current_pill_time[i]

            exo_E₂[i] += u1[mod_today + 1] * E₂_constant * (N_E₂_constant * exp(-k_a_E₂ * passed_time_from_current_pill) + L_E₂_constant * exp(-alpha1_E₂ * passed_time_from_current_pill) + M_E₂_constant * exp(-beta1_E₂ * passed_time_from_current_pill)) 
	    
            exo_P₄[i] += u2[mod_today + 1] * P₄_constant * (N_P₄_constant * exp(-k_a_P₄ * passed_time_from_current_pill) + L_P₄_constant * exp(-alpha1_P₄ * passed_time_from_current_pill) + M_P₄_constant * exp(-beta1_P₄ * passed_time_from_current_pill))            
        end
    end
    
    # Calculate endogenous hormone values
    C_E₂ = @. 1 + θ_1 * cos(2π * (t - θ_2))
    C_P₄ = @. 1 + θ_3 * cos(2π * (t - θ_4))
    
    E₂ = @. (e_0 + e_1 * GrF + e_2 * DomF + e_3 * Lut₄) * C_E₂ + η_E₂ * 1e3 * exo_E₂
    P₄ = @. (p_1 * Lut₃ + p_2 * Lut₄) * C_P₄ + η_P₄ * 1e3 * exo_P₄
    Inh = @. h_0 + h_1 * DomF + h_2 * Lut₂ + h_3 * Lut₃
    
    return length(E₂)==1 ? E₂[1] : E₂,
            length(P₄)==1 ? P₄[1] : P₄,
            length(Inh)==1 ? Inh[1] : Inh
end

# =======================
# Menstrual Cycle Dynamics
# =======================
"""
    menstrual_cycle_dynamics!(dx, x, h, p, t; n_cycles=nothing)

The dynamics for menstrual cycle.

# Arguments
- `dx::Vector{Float64}`: Vector variable containing differential values
- `x::Vector{Float64}`: Current state variables
- `h`: Delay function (returns the past values)
- `p::Vector{Any}`: Model parameter values
- `t::Float64`: Current time
- `n_cycles`: Optional parameter for the number of cycles

# Updates
- Modifies `dx` in place with calculated derivatives
"""
function menstrual_cycle_dynamics!(dx::Vector{Float64}, x::Vector{Float64}, h, p::Vector{Any}, t::Float64; n_cycles=nothing)
    # Extract state variables
    RP_LH, LH, RP_FSH, FSH, RcF, GrF, DomF, Sc_1, Sc_2, Lut_1, Lut_2, Lut_3, Lut_4 = x
    
    # Extract model parameters
    V_0_LH, V_1_LH, Km_LH, Ki_LH_P, k_LH, c_LH_P, c_LH_E, v, α_LH, Ki_FSH_Inh, w, k_FSH, c_FSH_P, c_FSH_E, V_FSH, α_FSH, b, c_1, q, c_2, c_3, c_4, d_1, d_2, k_1, k_2, k_3, k_4, α, γ, _, _, _, _, _, _, h_0, h_1, h_2, h_3, _, _, _, _, θ_5, θ_6, θ_7, θ_8 = p[1:48]
    
    # Control
    if isnothing(n_cycles)
        E2, P4, _ = calculate_auxiliary(x, p, t)
    else
        E2, P4, _ = calculate_auxiliary(x, p, t; n_cycles=n_cycles)
    end
    
    # Delayed state variables
    delayed_state = h(p, t - 1.5)
    DomF_delayed = delayed_state[7]
    Lut_2_delayed = delayed_state[11]
    Lut_3_delayed = delayed_state[12]
    Inh_delayed = h_0 + h_1 * DomF_delayed + h_2 * Lut_2_delayed + h_3 * Lut_3_delayed
    
    # Circadian cycle modulation
    C_LH = 1 + θ_5 * cos(2π * (t - θ_6))
    C_FSH = 1 + θ_7 * cos(2π * (t - θ_8))

    # Code optimization: precompute repeated calculations
    E2_squared = E2 * E2
    E2_eighed = E2_squared * E2_squared * E2_squared * E2_squared
    Km_LH_squared = Km_LH * Km_LH
    Km_LH_eighted = Km_LH_squared * Km_LH_squared * Km_LH_squared * Km_LH_squared
    LH_α = LH^α
    LH_γ = LH^γ

    # Precompute divisors for efficiency
    div1 = 1 / (1 + c_LH_E * E2)
    div2 = 1 / (1 + c_FSH_E * E2_squared)
    div3 = 1 / (1 + P4/q)
    div4 = 1 / v
    
    # Differential equations representing the dynamics
    dx[1] = (V_0_LH + V_1_LH * E2_eighed / (Km_LH_eighted + E2_eighed)) / (1 + (P4 / Ki_LH_P)) * C_LH - k_LH * (1 + c_LH_P * P4) * RP_LH * div1
    dx[2] = div4 * k_LH * (1 + c_LH_P * P4) * RP_LH * div1 - α_LH * LH
    dx[3] = V_FSH / (1 + (Inh_delayed / Ki_FSH_Inh) + P4/w) * C_FSH - k_FSH * (1 + c_FSH_P * P4) * RP_FSH * div2
    dx[4] = div4 * k_FSH * (1 + c_FSH_P * P4) * RP_FSH * div2 - α_FSH * FSH
    dx[5] = (b + c_1 * RcF) * FSH * div3 - c_2 * LH_α * RcF
    dx[6] = c_2 * LH_α * RcF - c_3 * GrF * LH
    dx[7] = c_3 * GrF * LH - c_4 * LH_γ * DomF
    dx[8] = c_4 * LH_γ * DomF - d_1 * Sc_1
    dx[9] = d_1 * Sc_1 - d_2 * Sc_2
    dx[10] = d_2 * Sc_2 - k_1 * Lut_1
    dx[11] = k_1 * Lut_1 - k_2 * Lut_2
    dx[12] = k_2 * Lut_2 - k_3 * Lut_3
    dx[13] = k_3 * Lut_3 - k_4 * Lut_4
end

# ===============================================
# Auxiliary Calculations (No Circadian Rhythm)
# ===============================================
"""
    calculate_auxiliary_no_circ(x, p, t; n_cycles=nothing)

Calculate E2, P4, Inh values without considering circadian rhythm.

# Arguments
- `x::Vector{Float64}`: State variables
- `p::Vector{Any}`: Model parameter values
- `t::Float64`: Current time
- `n_cycles`: Number of cycles (optional)

# Returns
- `(E2, P4, Inh)`: Calculated values
"""
function calculate_auxiliary_no_circ(x, p, t; n_cycles=nothing)
    GrF, DomF, Lut₂, Lut₃, Lut₄ = map(i -> x[i,:], [6:7;11:13])
    e_0, e_1, e_2, e_3, p_1, p_2, h_0, h_1, h_2, h_3, θ_1, θ_2, θ_3, θ_4, _, _, _, _, η_E₂, F_E₂, k_a_E₂, k_21_E₂, alpha1_E₂, beta1_E₂, η_P₄, F_P₄, k_a_P₄, k_21_P₄, alpha1_P₄, beta1_P₄, V_c_E₂, V_c_P₄, u1, u2, μ = p[[31:62; end-2:end]]
    
    # Control
    today = floor.(Int, t)
    exo_E₂ = zeros(length(t))
    exo_P₄ = zeros(length(t))

    # Code optimization
    E₂_constant = k_a_E₂ * F_E₂ / V_c_E₂
    N_E₂_constant = (k_21_E₂ - k_a_E₂) / ( (alpha1_E₂ - k_a_E₂) * (beta1_E₂ - k_a_E₂) )
    L_E₂_constant = (k_21_E₂ - alpha1_E₂) / ( (k_a_E₂ - alpha1_E₂) * (beta1_E₂ - alpha1_E₂) )
    M_E₂_constant = (k_21_E₂ - beta1_E₂) / ( (k_a_E₂ - beta1_E₂) * (alpha1_E₂ - beta1_E₂) )
    P₄_constant = k_a_P₄ * F_P₄ / V_c_P₄
    N_P₄_constant = (k_21_P₄ - k_a_P₄) / ( (alpha1_P₄ - k_a_P₄) * (beta1_P₄ - k_a_P₄) )
    L_P₄_constant = (k_21_P₄ - alpha1_P₄) / ( (k_a_P₄ - alpha1_P₄)*(beta1_P₄ - alpha1_P₄) )
    M_P₄_constant = (k_21_P₄ - beta1_P₄) / ( (k_a_P₄ - beta1_P₄)*(alpha1_P₄ - beta1_P₄) )
    current_pill_time = today .+ μ

    for i in eachindex(t)
        current_cycle = div(today[i], 28) + 1 # Current cycle number

        if today[i] > 0 # From the 1st cycle (simulation starts from day 0)
            past_days = 0:(today[i] - 1) # past days
            mod_past_days = mod.(past_days, 28) # dosing strategy have a cycle of 28 days
            i_pill_taking_time = past_days .+ μ # vector of taking time
            passed_time_from_pills = t[i] .- i_pill_taking_time # vector of passed time from each administration

            # Past pills' contribution
            exo_E₂[i] = sum(u1[mod_past_days .+ 1] .* E₂_constant .* (N_E₂_constant .* exp.(-k_a_E₂ .* passed_time_from_pills) .+ L_E₂_constant .* exp.(-alpha1_E₂ .* passed_time_from_pills) .+ M_E₂_constant .* exp.(-beta1_E₂ .* passed_time_from_pills))) 
            exo_P₄[i] = sum(u2[mod_past_days .+ 1] .* P₄_constant .* (N_P₄_constant .* exp.(-k_a_P₄ .* passed_time_from_pills) .+ L_P₄_constant .* exp.(-alpha1_P₄ .* passed_time_from_pills) .+ M_P₄_constant .* exp.(-beta1_P₄ .* passed_time_from_pills)))
        end

        # Current pill's contribution
        if (isnothing(n_cycles) || current_cycle <= n_cycles) && (t[i] - today[i] >= μ)
            mod_today = mod(today[i], 28)
            passed_time_from_current_pill = t[i] - current_pill_time[i]

            exo_E₂[i] += u1[mod_today + 1] * E₂_constant * (N_E₂_constant * exp(-k_a_E₂ * passed_time_from_current_pill) + L_E₂_constant * exp(-alpha1_E₂ * passed_time_from_current_pill) + M_E₂_constant * exp(-beta1_E₂ * passed_time_from_current_pill)) 
            exo_P₄[i] += u2[mod_today + 1] * P₄_constant * (N_P₄_constant * exp(-k_a_P₄ * passed_time_from_current_pill) + L_P₄_constant * exp(-alpha1_P₄ * passed_time_from_current_pill) + M_P₄_constant * exp(-beta1_P₄ * passed_time_from_current_pill)) 
        end
    end
    
    E₂ = @. (e_0 + e_1 * GrF + e_2 * DomF + e_3 * Lut₄) + η_E₂ * 1e3 * exo_E₂
    P₄ = @. (p_1 * Lut₃ + p_2 * Lut₄) + η_P₄ * 1e3 * exo_P₄
    Inh = @. h_0 + h_1 * DomF + h_2 * Lut₂ + h_3 * Lut₃
    
    return length(E₂)==1 ? E₂[1] : E₂,
            length(P₄)==1 ? P₄[1] : P₄,
            length(Inh)==1 ? Inh[1] : Inh
end

# ===============================================
# Dynamics without Circadian Rhythm
# ===============================================
"""
    menstrual_cycle_dynamics_no_circ!(dx, x, h, p, t; n_cycles=nothing)

The dynamics for menstrual cycle without considering circadian rhythm.

# Arguments
- `dx::Vector{Float64}`: Vector variable containing differential values
- `x::Vector{Float64}`: Current states
- `h`: Delay function (returns the past values)
- `p::Vector{Any}`: Model parameter values
- `t::Float64`: Current time
- `n_cycles`: Optional parameter for the number of cycles

# Updates
- Modifies `dx` in place with calculated derivatives
"""
function menstrual_cycle_dynamics_no_circ!(dx::Vector{Float64}, x::Vector{Float64}, h, p::Vector{Any}, t::Float64; n_cycles=nothing)
    # Extract state variables from the current state vector
    RP_LH, LH, RP_FSH, FSH, RcF, GrF, DomF, Sc_1, Sc_2, Lut_1, Lut_2, Lut_3, Lut_4 = x
    
    # Extract model parameters
    V_0_LH, V_1_LH, Km_LH, Ki_LH_P, k_LH, c_LH_P, c_LH_E, v, α_LH, Ki_FSH_Inh, w, k_FSH, c_FSH_P, c_FSH_E, V_FSH, α_FSH, b, c_1, q, c_2, c_3, c_4, d_1, d_2, k_1, k_2, k_3, k_4, α, γ, _, _, _, _, _, _, h_0, h_1, h_2, h_3, _, _, _, _, θ_5, θ_6, θ_7, θ_8 = p[1:48]

    # Control
    if isnothing(n_cycles)
        E2, P4, _ = calculate_auxiliary_no_circ(x, p, t)
    else
        E2, P4, _ = calculate_auxiliary_no_circ(x, p, t; n_cycles=n_cycles)
    end
    
    # Retrieve delayed state values for feedback regulation
    delayed_state = h(p, t - 1.5)
    DomF_delayed = delayed_state[7]
    Lut_2_delayed = delayed_state[11]
    Lut_3_delayed = delayed_state[12]
    Inh_delayed = h_0 + h_1 * DomF_delayed + h_2 * Lut_2_delayed + h_3 * Lut_3_delayed

    # Precompute constant values for efficiency
    E2_squared = E2 * E2
    E2_eighed = E2_squared * E2_squared * E2_squared * E2_squared
    Km_LH_squared = Km_LH * Km_LH
    Km_LH_eighted = Km_LH_squared * Km_LH_squared * Km_LH_squared * Km_LH_squared
    LH_α = LH^α
    LH_γ = LH^γ

    # Precompute divisors to avoid recalculation
    div1 = 1 / (1 + c_LH_E * E2)
    div2 = 1 / (1 + c_FSH_E * E2_squared)
    div3 = 1 / (1 + P4/q)
    div4 = 1 / v
    
    # Differential equations representing the dynamics of the menstrual cycle
    dx[1] = (V_0_LH + V_1_LH * E2_eighed / (Km_LH_eighted + E2_eighed)) / (1 + (P4 / Ki_LH_P)) - k_LH * (1 + c_LH_P * P4) * RP_LH * div1
    dx[2] = div4 * k_LH * (1 + c_LH_P * P4) * RP_LH * div1 - α_LH * LH
    dx[3] = V_FSH / (1 + (Inh_delayed / Ki_FSH_Inh) + P4/w) - k_FSH * (1 + c_FSH_P * P4) * RP_FSH * div2
    dx[4] = div4 * k_FSH * (1 + c_FSH_P * P4) * RP_FSH * div2 - α_FSH * FSH
    dx[5] = (b + c_1 * RcF) * FSH * div3 - c_2 * LH_α * RcF
    dx[6] = c_2 * LH_α * RcF - c_3 * GrF * LH
    dx[7] = c_3 * GrF * LH - c_4 * LH_γ * DomF
    dx[8] = c_4 * LH_γ * DomF - d_1 * Sc_1
    dx[9] = d_1 * Sc_1 - d_2 * Sc_2
    dx[10] = d_2 * Sc_2 - k_1 * Lut_1
    dx[11] = k_1 * Lut_1 - k_2 * Lut_2
    dx[12] = k_2 * Lut_2 - k_3 * Lut_3
    dx[13] = k_3 * Lut_3 - k_4 * Lut_4
end

# =======================
# History Function
# =======================
"""
    history_func(p, t)

Return the initial history vector for the DDE problem, used by the solver.

# Arguments
- `p::Vector{Any}`: Model parameter values
- `t::Float64`: Current time (ignored; history is constant)

# Returns
- `x0::Vector{Float64}`: Initial state vector for all variables
"""
function history_func(p, t)
    # Retrieve initial conditions for all state variables
    x0, _, _ = setup_initial_conditions_time_info(p)
    return x0
end

# =======================
# Solution Retrieval
# =======================
"""
    get_solution(candidate, p::Vector{Any}, num_cycle::Float64, is_circadian::Bool)

The method getting solution from the result population

# Arguments
- `candidate`: Metaheuristics solution value (vector of f, g, h, x)
- `p::Vector{Any}`: Parameter vector
- `num_cycle::Float64`: Number of cycle for simulation
- `is_circadian:Bool`: Boolean parameter for considering circadian rhythm

# Returns
- `sol`: The DDE solution object
- `u1`: E2 control vector (with zeros for non-dosing days)
- `u2`: P4 control vector (with zeros for non-dosing days)
"""
function get_solution(candidate, p::Vector{Any}, num_cycle::Float64, is_circadian::Bool)
    local_p = copy(p)

    # Split candidate.x into E2 and P4 controls, pad with zeros for non-dosing days
    taking_dates = Int64(length(candidate.x) / 2)
    dates_per_cycle = 28
    breaking_dates = dates_per_cycle - taking_dates
    u1 = vcat(candidate.x[1:taking_dates], zeros(breaking_dates))
    u2 = vcat(candidate.x[(taking_dates+1):end], zeros(breaking_dates))

    local_p[end-2] = u1
    local_p[end-1] = u2

    # Set up initial conditions and DDE problem
    x0, _, lags = setup_initial_conditions_time_info(local_p)
    tspan = (0.0, num_cycle*28.0)
    if is_circadian
        # Use DDE model with circadian modulation
        prob = DDEProblem(menstrual_cycle_dynamics!, x0, history_func, tspan, local_p; constant_lags=lags)
    else
        # Use DDE model without circadian modulation
        prob = DDEProblem(menstrual_cycle_dynamics_no_circ!, x0, history_func, tspan, local_p; constant_lags=lags)
    end
    # Solve the DDE system using RK4 method
    sol = solve(prob, MethodOfSteps(RK4()), dtmax=0.05)

    return sol, u1, u2
end