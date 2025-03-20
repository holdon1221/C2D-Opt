using Plots
using Metaheuristics

"""
    logger_drawer(A::Matrix{Float64}, alg_name::String, num_iter::Int64)

Logger drawer for optimization process

# Arguments
- `A::Matrix{Float64}`: Matrix has the objective function values for all population
- `u1::Vector{Float64}`: The current best strategy for the exogenous E₂ control
- `u2::Vector{Float64}`: The current best strategy for the exogenous P₄ control
- `alg_name::String`: The name of optimization algorithm
- `num_iter::Int64`: Iteration number
"""
function logger_drawer(A::Matrix{Float64}, u1::Vector{Float64}, u2::Vector{Float64}, alg_name::String, num_iter::Int64)
    scatter_graph = scatter(A[:,1], A[:,2], label=alg_name, title="Gen: $(num_iter)", xlabel="J₁", ylabel="J₂")

    days = 0:27
    # Update the bar graph for results
    bar_u1 = bar(days, u1, xlabel="Time [days]", ylabel="Exo E₂ [μg]", legend=false, 
                bar_width=0.7, xticks=days, color=:dodgerblue)
    bar_u2 = bar(days, u2, xlabel="Time [days]", ylabel="Exo P₄ [μg]", legend=false, 
                bar_width=0.7, xticks=days, color=:crimson)
    
    # Create a custom layout
    layout = @layout [a{0.7w} [b; c]]
    combined_plot = plot(scatter_graph, bar_u1, bar_u2, layout=layout, size=(1000, 600))

    display(combined_plot)

    return scatter_graph
end