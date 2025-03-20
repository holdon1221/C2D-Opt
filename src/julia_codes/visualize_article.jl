# This file is for visualizations in the article.
# We have to show the followings:
#   1. Near-optimal dosage in the morning time requires way less amount than the one in the other timings.
#   2. Near-optimal dosing strategy shows the less administration would be enough for birth controls.
#   3. Dynamcis of them
using Pkg
Pkg.activate(".")

using HDF5
using CSV
using DataFrames
using ArgParse
using PyPlot
using PyCall

include("models/parameter_loader.jl")
include("models/dde_model.jl")

function generate_summary_figure(optimal_dose_total::Vector{Float64}, start_time::Int64, end_time::Int64)
    administration_timings = start_time:end_time

    try
        fig, ax = plt.subplots(figsize=(16,6))

        # Highlight the morning time (6 AM - 12 PM)
        ax.axvspan(6, 12, color="#fff68f", alpha=0.6, label="Morning (6 AM - 12 PM)")


        # Highlight the night time (6 PM - 12 AM)
        ax.axvspan(18, 24, color="#D3D3D3", alpha=0.6, label="Night (6 PM - 12 AM)")

        # Draw the total near-optimal dosage for each administration timings
        ax.plot(administration_timings, optimal_dose_total,
                linestyle=:solid, color="#696969", linewidth=5, 
                marker="o", markersize=13, 
                label="Optimal Dose")

        # Set axis labels and ticks
        ax.set_xlabel("Administration timings [o'clock]", fontsize=12)
        ax.set_ylabel("Near-optimal total dosage [μg]", fontsize=12)
        ax.set_xticks([6, 12, 18, 24])

        # Morning text
        ylim = ax.get_ylim()
        ax.text(9, (ylim[1] + ylim[2]) / 2, "Morning", horizontalalignment="center", verticalalignment="center", fontsize=18, color=:black)

        # Night text
        ax.text(21, (ylim[1] + ylim[2]) / 2, "Night", horizontalalignment="center", verticalalignment="center", fontsize=18, color=:black)

        plt.tight_layout()

        return fig

    catch e
        println("Error occurred: ", e)
        return nothing
    end
end

function generate_summary_figure_bar(optimal_strategy_E2::Matrix{Float64}, optimal_strategy_P4::Matrix{Float64})
    try
        sum_optimal_strategy_E2 = sum(optimal_strategy_E2, dims=1)[1,:]
        sum_optimal_strategy_P4 = sum(optimal_strategy_P4, dims=1)[1,:]
        
        fig, ax1 = plt.subplots(figsize=(16,6))

        # Draw the total near-optimal dosage for each administration timings
        bar_width = 0.4
        
        # Visualize all EE strategies
        ax1.bar(collect(1:24) .- bar_width/2, sum_optimal_strategy_E2, color="#87CEFA", label="EE", width=bar_width)
        ax1.set_xlabel("Administration timings [o'clock]")
        ax1.set_ylabel("EE dosage [μg]", color="black")
        ax1.set_xticks(collect(1:24))

        # Visualize all DNG strategies
        ax2 = ax1.twinx()
        ax2.bar(collect(1:24) .+ bar_width/2, sum_optimal_strategy_P4, color="#FFB6C1", label="DNG", width=bar_width)
        # ax2.set_ylim([0, maximum([1, 1.05*maximum(sum_optimal_strategy_P4)])])
        ax2.set_ylabel("DNG dosage [μg]", color="black")

        # Integrate legends in the one legend
        handles1, labels1 = ax1.get_legend_handles_labels()  # EE
        handles2, labels2 = ax2.get_legend_handles_labels()  # DNG

        ax1.legend(vcat(handles1, handles2), vcat(labels1, labels2), loc="upper right")

        plt.tight_layout()

        return fig

    catch e
        println("Error occurred: ", e)
        return nothing
    end
end

function generate_bar_graph_with_input(strategy_E2::Vector{Float64}, strategy_P4::Vector{Float64})
    fig, ax1 = plt.subplots(figsize=(8,3))

    bar_width = 0.4
    
    # Visualize all EE strategies
    ax1.bar(collect(0:27) .- bar_width/2, strategy_E2, color="#87CEFA", label="EE", width=bar_width)
    ax1.set_xlabel("Time [day]")
    ax1.set_ylabel("EE dosage [μg]", color="black")
    ax1.set_xticks(collect(0:27))
    ax1.set_ylim([0,30])

    # Add green shading from day 21 to 27
    ax1.axvspan(21, 27, color="lightgreen", alpha=0.4)

    # Add text "Break" in the center of the green-shaded region
    ylim_EE = ax1.get_ylim()
    plt.text(24, (ylim_EE[1] + ylim_EE[2]) / 2, "Break", horizontalalignment="center", verticalalignment="center", fontsize=12, color="black")

    # Visualize all DNG strategies
    ax2 = ax1.twinx()
    ax2.bar(collect(0:27) .+ bar_width/2, strategy_P4, color="#FFB6C1", label="DNG", width=bar_width)
    ax2.set_ylim([0, maximum([1, 1.05*maximum(strategy_P4)])])
    ax2.set_ylabel("DNG dosage [μg]", color="black")

    # Integrate legends in the one legend
    handles1, labels1 = ax1.get_legend_handles_labels()  # EE
    handles2, labels2 = ax2.get_legend_handles_labels()  # DNG

    ax1.legend(vcat(handles1, handles2), vcat(labels1, labels2), loc="upper right")

    plt.tight_layout()

    return fig
end

function generate_figure_dynamics(optimal_strategy_E2_1::Vector{Float64}, optimal_strategy_P4_1::Vector{Float64}, optimal_strategy_E2_2::Vector{Float64}, optimal_strategy_P4_2::Vector{Float64}, const_strategy_E2_1::Vector{Float64}, const_strategy_P4_1::Vector{Float64}, const_strategy_E2_2::Vector{Float64}, const_strategy_P4_2::Vector{Float64}, clock::Vector{Float64}, num_cycle::Int64, n_cycles=nothing)
    try
        # Upload parameter values
        μ_1 = clock[1]/24.0
        μ_2 = clock[2]/24.0

        path_para = joinpath("..", "..", "res", "parameter_values.csv")
        p_1 = load_parameter(path_para, μ_1)
        p_2 = load_parameter(path_para, μ_2)

        # Set the parameter vector without control
        p_wo_control = copy(p_1)

        # Set the parameter vector with a uniform dose strategy
        p_const_1 = copy(p_1)
        p_const_1[end-2] = const_strategy_E2_1
        p_const_1[end-1] = const_strategy_P4_1

        # Set the parameter vector with a uniform dose strategy
        p_const_2 = copy(p_2)
        p_const_2[end-2] = const_strategy_E2_2
        p_const_2[end-1] = const_strategy_P4_2

        # Set the control variables
        p_optimal_1 = copy(p_1)
        p_optimal_1[end-2] = optimal_strategy_E2_1
        p_optimal_1[end-1] = optimal_strategy_P4_1

        # Set the control variables
        p_optimal_2 = copy(p_2)
        p_optimal_2[end-2] = optimal_strategy_E2_2
        p_optimal_2[end-1] = optimal_strategy_P4_2

        # Set the DDE problems
        tspan = (0, num_cycle*28)
        x0_wo_control, _, lags = setup_initial_conditions_time_info(p_wo_control)

        # if n_cycles is nothing, then there is no stopping cycles
        if isnothing(n_cycles)
            prob_wo_control = DDEProblem(
                menstrual_cycle_dynamics!,
                x0_wo_control, history_func, tspan, p_wo_control; 
                constant_lags=lags
            )
        else
            prob_wo_control = DDEProblem(
                (dx, x, h, p, t) -> menstrual_cycle_dynamics!(dx, x, h, p, t; n_cycles=n_cycles),
                x0_wo_control, history_func, tspan, p_wo_control; 
                constant_lags=lags
            )
        end

        x0_const, _, lags = setup_initial_conditions_time_info(p_const_1)

        if isnothing(n_cycles)
            prob_const_1 = DDEProblem(
                menstrual_cycle_dynamics!,
                x0_const, history_func, tspan, p_const_1; 
                constant_lags=lags
            )
        else
            prob_const_1 = DDEProblem(
                (dx, x, h, p, t) -> menstrual_cycle_dynamics!(dx, x, h, p, t; n_cycles=n_cycles),
                x0_const, history_func, tspan, p_const_1; 
                constant_lags=lags
            )
        end

        x0_const, _, lags = setup_initial_conditions_time_info(p_const_2)

        if isnothing(n_cycles)
            prob_const_2 = DDEProblem(
                menstrual_cycle_dynamics!,
                x0_const, history_func, tspan, p_const_2; 
                constant_lags=lags
            )
        else
            prob_const_2 = DDEProblem(
                (dx, x, h, p, t) -> menstrual_cycle_dynamics!(dx, x, h, p, t; n_cycles=n_cycles),
                x0_const, history_func, tspan, p_const_2; 
                constant_lags=lags
            )
        end

        x0_optimal, _, lags = setup_initial_conditions_time_info(p_optimal_1)

        if isnothing(n_cycles)
            prob_optimal_1 = DDEProblem(
                menstrual_cycle_dynamics!,
                x0_optimal, history_func, tspan, p_optimal_1; 
                constant_lags=lags
            )
        else
            prob_optimal_1 = DDEProblem(
                (dx, x, h, p, t) -> menstrual_cycle_dynamics!(dx, x, h, p, t; n_cycles=n_cycles),
                x0_optimal, history_func, tspan, p_optimal_1; 
                constant_lags=lags
            )
        end

        x0_optimal, _, lags = setup_initial_conditions_time_info(p_optimal_2)

        if isnothing(n_cycles)
            prob_optimal_2 = DDEProblem(
                menstrual_cycle_dynamics!,
                x0_optimal, history_func, tspan, p_optimal_2; 
                constant_lags=lags
            )
        else
            prob_optimal_2 = DDEProblem(
                (dx, x, h, p, t) -> menstrual_cycle_dynamics!(dx, x, h, p, t; n_cycles=n_cycles),
                x0_optimal, history_func, tspan, p_optimal_2; 
                constant_lags=lags
            )
        end

        # Solve the mathematical model
        sol_wo_control = solve(prob_wo_control, MethodOfSteps(RK4()), dtmax=0.05)
        LH_wo_control = sol_wo_control[2,:]
        _, P₄_wo_control, _ = calculate_auxiliary(sol_wo_control, p_wo_control, sol_wo_control.t; n_cycles=7)

        sol_const_1 = solve(prob_const_1, MethodOfSteps(RK4()), dtmax=0.05)
        LH_const_1 = sol_const_1[2,:]
        _, P₄_const_1, _ = calculate_auxiliary(sol_const_1, p_const_1, sol_const_1.t; n_cycles=7)

        sol_const_2 = solve(prob_const_2, MethodOfSteps(RK4()), dtmax=0.05)
        LH_const_2 = sol_const_2[2,:]
        _, P₄_const_2, _ = calculate_auxiliary(sol_const_2, p_const_2, sol_const_2.t; n_cycles=7)

        sol_optimal_1 = solve(prob_optimal_1, MethodOfSteps(RK4()), dtmax=0.05)
        LH_optimal_1 = sol_optimal_1[2,:]
        _, P₄_optimal_1, _ = calculate_auxiliary(sol_optimal_1, p_optimal_1, sol_optimal_1.t; n_cycles=7)

        sol_optimal_2 = solve(prob_optimal_2, MethodOfSteps(RK4()), dtmax=0.05)
        LH_optimal_2 = sol_optimal_2[2,:]
        _, P₄_optimal_2, _ = calculate_auxiliary(sol_optimal_2, p_optimal_2, sol_optimal_2.t; n_cycles=7)

        # Draw them
        fig, axs = plt.subplots(1, 2, figsize=(16, 3))#, sharex="col")

        # Morning P₄ dynamics
        ax1 = axs[1]

        if !isnothing(n_cycles)
            ax1.axvspan(28*n_cycles, 28*num_cycle, color="lightgreen", alpha=0.4)
        end
        for i in 28:28:28*num_cycle
            ax1.axvline(i, linestyle=:dashed, color=:black)
        end

        ax1.plot(sol_wo_control.t, P₄_wo_control, color=:black, linewidth=1, label="No administration")
        ax1.plot(sol_const_1.t, P₄_const_1, color="#1F77B4", linestyle=:dotted, linewidth=1.5, label="Constant")
        ax1.plot(sol_optimal_1.t, P₄_optimal_1, color="#E74C3c", linewidth=3, label="Optimized")
        ax1.axhline(y=3, color=:red, linestyle=:dashed)
        ax1.text(sol_optimal_1.t[1] - 1, 3.1, "Contraceptive threshold", fontsize=10, verticalalignment="bottom", horizontalalignment="left")
        ax1.set_ylabel("P₄ [ng/ml]")

        ax1.set_xticks(14:28:(28*num_cycle-14))
        ax1.set_xticklabels([i for i in 1:num_cycle])
        ax1.set_xlabel("Number of cycles")

        # Night P₄ dynamics
        ax2 = axs[2]

        if !isnothing(n_cycles)
            ax2.axvspan(28*n_cycles, 28*num_cycle, color="lightgreen", alpha=0.4)
        end
        for i in 28:28:28*num_cycle
            ax2.axvline(i, linestyle=:dashed, color=:black)
        end

        ax2.plot(sol_wo_control.t, P₄_wo_control, color=:black, linewidth=1, label="No administration")
        ax2.plot(sol_const_2.t, P₄_const_2, color="#1F77B4", linestyle=:dotted, linewidth=1.5, label="Constant")
        ax2.plot(sol_optimal_2.t, P₄_optimal_2, color="#E74C3c", linewidth=3, label="Optimized")
        ax2.axhline(y=3, color=:red, linestyle=:dashed)
        ax2.text(sol_optimal_2.t[1] - 1, 3.1, "Contraceptive threshold", fontsize=10, verticalalignment="bottom", horizontalalignment="left")
        ax2.set_ylabel("P₄ [ng/ml]")

        ax2.set_xticks(14:28:(28*num_cycle-14))
        ax2.set_xticklabels([i for i in 1:num_cycle])
        ax2.set_xlabel("Number of cycles")

        # add a common legend
        h, l = ax1.get_legend_handles_labels()
        fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 1.00), ncol=3)
        # fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 1.00), ncol=3)

        plt.tight_layout()
        # fig.subplots_adjust(top=0.92)
        fig.subplots_adjust(top=0.85)

        return fig, ((sol_wo_control.t, P₄_wo_control, LH_wo_control), (sol_const_1.t, P₄_const_1, LH_const_1), (sol_const_2.t, P₄_const_2, LH_const_2), (sol_optimal_1.t, P₄_optimal_1, LH_optimal_1), (sol_optimal_2.t, P₄_optimal_2, LH_optimal_2))

    catch e
        println("Error occurred: ", e)
        return nothing
    end
end


function main()
    s = ArgParseSettings()
    
    @add_arg_table! s begin
        "--base_dir", "-b"
        help = "Base directory containing the cases"
        arg_type = String
        required = true

        "--timing", "-t"
        help = "Drawing timings (provide two values)"
        arg_type = Float64
        nargs = 2
        required = true

        "--dpi", "-d"
        help = "The resolution of the outputs"
        arg_type = Int
        default = 300

        "--mode", "-m"
        help = "Set the mode only for saving"
        arg_type = String
        default = "display"
    end

    parsed_args = parse_args(s)

    res_directory = parsed_args["base_dir"]
    raw_data_path = joinpath(res_directory, "raw_data")
    mkpath(raw_data_path)

    dpi = parsed_args["dpi"]
    timing = parsed_args["timing"]
    mode = parsed_args["mode"]
    isdisplay = (mode=="display")
    int_timing = [Int(t) for t in timing]
    
    start_time = 1
    end_time = 24

    optimal_dose_total = h5read(joinpath(res_directory, "Optimal_Dose.h5"), "E2_P4_$(start_time)-$(end_time)")
    optimal_strategy_E2 = h5read(joinpath(res_directory, "Optimal_Strategy.h5"), "E2_$(start_time)-$(end_time)")
    selected_optimal_strategy_E2_1 = optimal_strategy_E2[:,int_timing[1]]
    # selected_optimal_strategy_E2_1[selected_optimal_strategy_E2_1 .< 1] .= 0
    selected_optimal_strategy_E2_2 = optimal_strategy_E2[:,int_timing[2]]
    # selected_optimal_strategy_E2_2[selected_optimal_strategy_E2_2 .< 1] .= 0
    optimal_strategy_P4 = h5read(joinpath(res_directory, "Optimal_Strategy.h5"), "P4_$(start_time)-$(end_time)")
    selected_optimal_strategy_P4_1 = optimal_strategy_P4[:,int_timing[1]]
    # selected_optimal_strategy_P4_1[selected_optimal_strategy_P4_1 .< 1] .= 0
    selected_optimal_strategy_P4_2 = optimal_strategy_P4[:,int_timing[2]]
    # selected_optimal_strategy_P4_2[selected_optimal_strategy_P4_2 .< 1] .= 0

    const_strategy = CSV.read("../../res/Constant_dosing.csv", DataFrame; datarow=1, limit=24)
    const_strategy_E2 = const_strategy[:,2] ./ 1e3 # (ng to μg) 
    selected_const_strategy_E2_1 = vcat(const_strategy_E2[int_timing[1]]*ones(21), zeros(7))
    selected_const_strategy_E2_2 = vcat(const_strategy_E2[int_timing[2]]*ones(21), zeros(7))
    const_strategy_P4 = const_strategy[:,3] ./ 1e3 # (ng to μg)
    selected_const_strategy_P4_1 = vcat(const_strategy_P4[int_timing[1]]*ones(21), zeros(7))
    selected_const_strategy_P4_2 = vcat(const_strategy_P4[int_timing[2]]*ones(21), zeros(7))

    # #
    const_str_bar_fig = generate_bar_graph_with_input(selected_const_strategy_E2_1, selected_const_strategy_P4_1)
    if isnothing(const_str_bar_fig)
        println("Warning: Failed to generate const startegy bars for the morning.")
    elseif isdisplay
        const_str_bar_fig_path = joinpath(res_directory, "fig_const_stg_bar_clock$(int_timing[1])")
        const_str_bar_fig.savefig(const_str_bar_fig_path * ".png", dpi=dpi)
        const_str_bar_fig.savefig(const_str_bar_fig_path * ".eps")
        const_str_bar_fig.savefig(const_str_bar_fig_path * ".tiff", dpi=dpi)
    end

    const_str_bar_fig = generate_bar_graph_with_input(selected_const_strategy_E2_2, selected_const_strategy_P4_2)
    if isnothing(const_str_bar_fig)
        println("Warning: Failed to generate const startegy bars for the night.")
    elseif isdisplay
        const_str_bar_fig_path = joinpath(res_directory, "fig_const_stg_bar_clock$(int_timing[2])")
        const_str_bar_fig.savefig(const_str_bar_fig_path * ".png", dpi=dpi)
        const_str_bar_fig.savefig(const_str_bar_fig_path * ".eps")
        const_str_bar_fig.savefig(const_str_bar_fig_path * ".tiff", dpi=dpi)
    end

    CSV.write(joinpath(raw_data_path, "fig_const_stg_bar_clock$(int_timing[1]).csv"), DataFrame(day=collect(0:27),EE=selected_const_strategy_E2_1, DNG=selected_const_strategy_P4_1))
    CSV.write(joinpath(raw_data_path, "fig_const_stg_bar_clock$(int_timing[2]).csv"), DataFrame(day=collect(0:27),EE=selected_const_strategy_E2_2, DNG=selected_const_strategy_P4_2))

    # #
    str_bar_fig = generate_bar_graph_with_input(selected_optimal_strategy_E2_1, selected_optimal_strategy_P4_1)
    if isnothing(str_bar_fig)
        println("Warning: Failed to generate optimal startegy bars for the morning.")
    elseif isdisplay
        str_bar_fig_path = joinpath(res_directory, "fig_stg_bar_clock$(int_timing[1])")
        str_bar_fig.savefig(str_bar_fig_path * ".png", dpi=dpi)
        str_bar_fig.savefig(str_bar_fig_path * ".eps")
        str_bar_fig.savefig(str_bar_fig_path * ".tiff", dpi=dpi)
    end

    str_bar_fig = generate_bar_graph_with_input(selected_optimal_strategy_E2_2, selected_optimal_strategy_P4_2)
    if isnothing(str_bar_fig)
        println("Warning: Failed to generate optimal startegy bars for the night.")
    elseif isdisplay
        str_bar_fig_path = joinpath(res_directory, "fig_stg_bar_clock$(int_timing[2])")
        str_bar_fig.savefig(str_bar_fig_path * ".png", dpi=dpi)
        str_bar_fig.savefig(str_bar_fig_path * ".eps")
        str_bar_fig.savefig(str_bar_fig_path * ".tiff", dpi=dpi)
    end

    CSV.write(joinpath(raw_data_path, "fig_stg_bar_clock$(int_timing[1]).csv"), DataFrame(day=collect(0:27),EE=selected_optimal_strategy_E2_1, DNG=selected_optimal_strategy_P4_1))
    CSV.write(joinpath(raw_data_path, "fig_stg_bar_clock$(int_timing[2]).csv"), DataFrame(day=collect(0:27),EE=selected_optimal_strategy_E2_2, DNG=selected_optimal_strategy_P4_2))

    # #
    poc_fig, (sol_wo_control, sol_const_1, sol_const_2, sol_optimal_1, sol_optimal_2) = generate_figure_dynamics(selected_optimal_strategy_E2_1, selected_optimal_strategy_P4_1, selected_optimal_strategy_E2_2, selected_optimal_strategy_P4_2, selected_const_strategy_E2_1, selected_const_strategy_P4_1, selected_const_strategy_E2_2, selected_const_strategy_P4_2, timing, 3)
    if isnothing(poc_fig)
        println("Warning: Failed to generate POC figure.")
    elseif isdisplay
        poc_fig_path = joinpath(res_directory, "fig_poc")
        poc_fig.savefig(poc_fig_path * ".png", dpi=dpi)
        poc_fig.savefig(poc_fig_path * ".eps")
        poc_fig.savefig(poc_fig_path * ".tiff", dpi=dpi)
    end

    CSV.write(joinpath(raw_data_path, "fig_poc_wo_control.csv"), DataFrame(time=sol_wo_control[1], P4=sol_wo_control[2], LH=sol_wo_control[3]))
    CSV.write(joinpath(raw_data_path, "fig_poc_const_clock$(int_timing[1]).csv"), DataFrame(time=sol_const_1[1], P4=sol_const_1[2], LH=sol_const_1[3]))
    CSV.write(joinpath(raw_data_path, "fig_poc_const_clock$(int_timing[2]).csv"), DataFrame(time=sol_const_2[1], P4=sol_const_2[2], LH=sol_const_2[3]))
    CSV.write(joinpath(raw_data_path, "fig_poc_optimal_clock$(int_timing[1]).csv"), DataFrame(time=sol_optimal_1[1], P4=sol_optimal_1[2], LH=sol_optimal_1[3]))
    CSV.write(joinpath(raw_data_path, "fig_poc_optimal_clock$(int_timing[2]).csv"), DataFrame(time=sol_optimal_2[1], P4=sol_optimal_2[2], LH=sol_optimal_2[3]))

    sum_fig = generate_summary_figure(optimal_dose_total, start_time, end_time)
    if isnothing(sum_fig)
        println("Warning: Failed to generate summary figure.")
    elseif isdisplay
        sum_fig_path = joinpath(res_directory, "fig_sum")
        sum_fig.savefig(sum_fig_path * ".png", dpi=dpi)
        sum_fig.savefig(sum_fig_path * ".eps")
        sum_fig.savefig(sum_fig_path * ".tiff", dpi=dpi)
    end

    sum_fig_bar = generate_summary_figure_bar(optimal_strategy_E2, optimal_strategy_P4)
    if isnothing(sum_fig_bar)
        println("Warning: Failed to generate summary bar figure.")
    elseif isdisplay
        sum_fig_bar_path = joinpath(res_directory, "fig_sum_bar")
        sum_fig_bar.savefig(sum_fig_bar_path * ".png", dpi=dpi)
        sum_fig_bar.savefig(sum_fig_bar_path * ".eps")
        sum_fig_bar.savefig(sum_fig_bar_path * ".tiff", dpi=dpi)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end