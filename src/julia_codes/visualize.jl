using Pkg
Pkg.activate(".")

using Plots
using Serialization
using Metaheuristics
using FileIO
using ArgParse
using HDF5

include("models/parameter_loader.jl")
include("utils/save_results.jl")

# Function to select the directory with the smallest dose for a specific clock
function select_smallest_dose_directory(clock_name, dirs)
    smallest_dose = Inf
    best_file_path = nothing

    for dir in dirs # Pick one directory for various population strategies
        dir_path = filter(d -> occursin(clock_name, d), dir)
        
        if length(dir_path) == 1
            file_path = joinpath(dir_path[1], "result_variable.jls")
            res = deserialize(file_path)
            smallest = select_states(res, 3.0, 1, true)[1]
            if smallest.f[2] < smallest_dose
                smallest_dose = smallest.f[2]
                best_file_path = file_path
            end
        end
    end

    if best_file_path === nothing
        println("Warning: No valid directory found for ", clock_name)
    end
    
    return best_file_path
end

# Function to handle the input directories and generate plots/heatmaps
function visualize(base_dir::String, cases::Union{String, Vector{String}}, p::Vector{Any}; start_time=7, end_time=22, plot_type="both")
    # Convert cases to a vector if a single case is provided as a string
    if typeof(cases) == String
        cases = [cases]
    end
    
    # Generate full paths for each case
    directories = [joinpath(base_dir, case) for case in cases]
    
    visualize(directories, p; start_time=start_time, end_time=end_time, plot_type=plot_type)
end

# Overloaded function to handle a single directory input
function visualize(directory::String, p::Vector{Any}; start_time=7, end_time=22, plot_type="both")
    visualize([directory], p; start_time=start_time, end_time=end_time, plot_type=plot_type)
end

# Core function to process directories (modified from previous implementation)
function visualize(directories::Vector{String}, p::Vector{Any}; start_time=7, end_time=22, plot_type="both")
    time_slot = start_time:end_time
    clocks = ["Clock$(i).0" for i in time_slot]
    days = 1:28
    
    e2_doses, p4_doses, sum_e2_doses, sum_p4_doses = Vector{Float64}[], Vector{Float64}[], Float64[], Float64[]
    
    # Load directories
    all_dirs = [readdir(dir, join=true) for dir in directories]
    
    # Determine the result directory
    if length(directories) == 1
        result_dir_path = directories[1]
    else
        folder_suffixes = [splitpath(dir)[end] for dir in directories]
        sort!(folder_suffixes) # The same directory with a different order
        result_dir_name = join(folder_suffixes, "_")
        common_parent = joinpath(splitpath(directories[1])[1:end-1]...)
        result_dir_path = joinpath(common_parent, result_dir_name)
        
        # Create the result directory only if it doesn't exist
        if !isdir(result_dir_path)
            mkdir(result_dir_path)
        else
            println("Directory already exists: ", result_dir_path, " - Using existing directory.")
        end
    end
    
    for i in eachindex(time_slot)
        clock_name = clocks[i]
        file_path = select_smallest_dose_directory(clock_name, all_dirs)
        if file_path !== nothing
            res = deserialize(file_path)
            smallest = select_states(res, 3.0, 1, true)[1]
            
            push!(e2_doses, vcat(smallest.x[1:21], zeros(7)))
            push!(sum_e2_doses, sum(smallest.x[1:21]))
            push!(p4_doses, vcat(smallest.x[22:end], zeros(7)))
            push!(sum_p4_doses, sum(smallest.x[22:end]))
        end
    end
    
    # Generate heatmap data by concatenating doses across clocks
    heatmap_data_e2 = hcat(e2_doses...)
    heatmap_data_p4 = hcat(p4_doses...)
    heatmap_data_all = heatmap_data_p4 .+ heatmap_data_e2
    
    if plot_type == "heatmap" || plot_type == "both"
        # Create heatmap for E2
        h_E₂ = heatmap(
            days, clocks, heatmap_data_e2', 
            xlabel="Time [days]", ylabel="Administration timing [o'clock]",
            color=cgrad(:RdPu), size=(600, 600), 
            xticks=1:28, yticks=((time_slot.-time_slot[1]) .+0.5, time_slot)
        )
        h5write(joinpath(result_dir_path, "Optimal_Strategy.h5"), "/E2_$(start_time)-$(end_time)", heatmap_data_e2)
        savefig(h_E₂, joinpath(result_dir_path, "Optimal_Strategy_E2_$(start_time)-$(end_time).png"))
        
        # Create heatmap for P4
        h_P₄ = heatmap(
            days, clocks, heatmap_data_p4', 
            xlabel="Time [days]", ylabel="Administration timing [o'clock]", 
            color=cgrad(:RdPu), size=(600, 600), 
            xticks=1:28, yticks=((time_slot.-time_slot[1]) .+0.5, time_slot)
        )
        h5write(joinpath(result_dir_path, "Optimal_Strategy.h5"), "/P4_$(start_time)-$(end_time)", heatmap_data_p4)
        savefig(h_P₄, joinpath(result_dir_path, "Optimal_Strategy_P4_$(start_time)-$(end_time).png"))
        
        # Create heatmap for total doses
        h_all = heatmap(
            days, clocks, heatmap_data_all', 
            xlabel="Time [days]", ylabel="Administration timing [o'clock]",
            color=cgrad(:RdPu), size=(600, 600), 
            xticks=1:28, yticks=((time_slot.-time_slot[1]) .+0.5, time_slot)
        )
        h5write(joinpath(result_dir_path, "Optimal_Strategy.h5"), "/E2_P4_$(start_time)-$(end_time)", heatmap_data_all)
        savefig(h_all, joinpath(result_dir_path, "Optimal_Strategy_E2_P4_$(start_time)-$(end_time).png"))
    end

    if plot_type == "plot" || plot_type == "both"
        # Create plot for E2
        p_E₂ = plot(
            time_slot, sum_e2_doses, xlabel="Dosing time [o'clock]", ylabel="Optimal dose of E₂ [μg]", 
            lw=4, xticks=time_slot, color=:dodgerblue, label=false, grid=false
        )
        h5write(joinpath(result_dir_path, "Optimal_Dose.h5"), "/E2_$(start_time)-$(end_time)", sum_e2_doses)
        savefig(p_E₂, joinpath(result_dir_path, "Optimal_Dose_E2_$(start_time)-$(end_time).png"))
        
        # Create plot for P4
        p_P₄ = plot(
            time_slot, sum_p4_doses, xlabel="Dosing time [o'clock]", ylabel="Optimal dose of P₄ [μg]", 
            lw=4, xticks=time_slot, color=:crimson, label=false, grid=false
        )
        h5write(joinpath(result_dir_path, "Optimal_Dose.h5"), "/P4_$(start_time)-$(end_time)", sum_p4_doses)
        savefig(p_P₄, joinpath(result_dir_path, "Optimal_Dose_P4_$(start_time)-$(end_time).png"))
        
        # Create plot for total doses
        p_all = plot(
            time_slot, sum_e2_doses .+ sum_p4_doses, xlabel="Dosing time [o'clock]", ylabel="Optimal total dose [μg]", 
            lw=4, xticks=time_slot, color=:darkolivegreen4, label=false, grid=false
        )
        h5write(joinpath(result_dir_path, "Optimal_Dose.h5"), "/E2_P4_$(start_time)-$(end_time)", sum_e2_doses .+ sum_p4_doses)
        savefig(p_all, joinpath(result_dir_path, "Optimal_Dose_E2_P4_$(start_time)-$(end_time).png"))
    end
    
    println("Results saved in: ", result_dir_path)
end

# Example usage:
# visualize("../../res/obj_sigmoidP4FSH/no_thresholding/tol4/", ["N700", "N800"])
# visualize("../../res/obj_sigmoidP4FSH/no_thresholding/tol4/", "N700")
# visualize("../../res/obj_sigmoidP4FSH/no_thresholding/tol4/N700")

# Main entry point for command line execution
function main()
    s = ArgParseSettings()
    
    @add_arg_table! s begin
        "--base_dir", "-b"
        help = "Base directory containing the cases"
        arg_type = String
        required = true

        "--cases", "-c"
        help = "Comma-separated list of cases (e.g., N700,N800)"
        arg_type = String
        required = true

        "--start_time", "-s"
        help = "Start time for the analysis (default: 7)"
        arg_type = Int
        default = 7

        "--end_time", "-e"
        help = "End time for the analysis (default: 22)"
        arg_type = Int
        default = 22

        "--plot_type", "-p"
        help = "Type of plot to generate (options: heatmap, plot, both) (default: both)"
        arg_type = String
        default = "both"
    end
    
    parsed_args = parse_args(s)

    base_dir = parsed_args["base_dir"]
    cases = split(parsed_args["cases"], ",")
    cases = String.(cases) # Convert SubString to String
    start_time = parsed_args["start_time"]
    end_time = parsed_args["end_time"]
    plot_type = parsed_args["plot_type"]

    # Load parameter values and set the initial condition, tspan and lag time
    path_para = joinpath("..", "..", "res", "parameter_values.csv")
    μ = 0.0
    p = load_parameter(path_para, μ)

    visualize(base_dir, cases, p; start_time=start_time, end_time=end_time, plot_type=plot_type)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end