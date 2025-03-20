using Pkg
Pkg.activate(".")

using Metaheuristics
using Serialization
using PyPlot
using HDF5

include("models/parameter_loader.jl")
include("models/dde_model.jl")

include("models/parameter_loader.jl")
include("visualize_article.jl")

function main()
    path_para = joinpath("..", "..", "res", "parameter_values.csv")
    p = load_parameter(path_para, 10.0/24)
    res_directory = joinpath("..", "..", "res", "RBA1.7_0.01", "N400_N400_hyp_N500_N500_hyp")

    initial_res = deserialize(joinpath(res_directory, "initial_example.jls"))

    maxP4_vals = [child.f[1] for child in initial_res.population]
    total_amount = [child.f[2] for child in initial_res.population]

    index_rep = 5
    EE_rep = vcat(initial_res.population[index_rep].x[1:21], zeros(7))
    DNG_rep = vcat(initial_res.population[index_rep].x[22:end], zeros(7))

    fig, ax = plt.subplots(figsize=(8,6))

    # Scatter of initial population
    ax.scatter(total_amount, maxP4_vals)
    ax.scatter(total_amount[index_rep], maxP4_vals[index_rep], color=:red)
    ax.axvline(x=total_amount[index_rep], color=:red, linestyle=:dashed, linewidth=2)
    ax.axhline(y=maxP4_vals[index_rep], color=:red, linestyle=:dashed, linewidth=2)
    ax.set_xlabel("Total amount", fontsize=12)
    ax.set_ylabel("Maximum P₄", fontsize=12)

    plt.tight_layout()

    fig.savefig("tmp.png")

    fig_rep_strategy = generate_bar_graph_with_input(EE_rep, DNG_rep)
    
    fig_rep_strategy.savefig("tmp2.png")

    p_wo_control = copy(p)

    tspan = (0, 3.0*28)
    x0_wo_control, _, lags = setup_initial_conditions_time_info(p_wo_control)
    prob_wo_control = DDEProblem(
                menstrual_cycle_dynamics!,
                x0_wo_control, history_func, tspan, p_wo_control; 
                constant_lags=lags
            )
    sol_wo_control = solve(prob_wo_control, MethodOfSteps(RK4()), dtmax=0.05)
    _, P₄_wo_control, _ = calculate_auxiliary(sol_wo_control, p_wo_control, sol_wo_control.t)

    p[end-2] = EE_rep
    p[end-1] = DNG_rep

    # Set the DDE problems
    x0, _, lags = setup_initial_conditions_time_info(p)

    prob = DDEProblem(menstrual_cycle_dynamics!, x0, history_func, tspan, p; constant_lags=lags)

    sol = solve(prob, MethodOfSteps(RK4()), dtmax=0.05)
    _, P₄, _ = calculate_auxiliary(sol, p, sol.t)

    fig_rep_dynamics, ax_rep_dynamics = plt.subplots(figsize=(8,3))

    for i in 28:28:28*3
        ax_rep_dynamics.axvline(i, linestyle=:dashed, color=:black)
    end

    ax_rep_dynamics.plot(sol_wo_control.t, P₄_wo_control, color=:black, linewidth=1)
    ax_rep_dynamics.plot(sol.t, P₄, color="#E74C3c", linewidth=3)

    ax_rep_dynamics.axhline(y=3, color=:red, linestyle=:dashed)
    ax_rep_dynamics.text(sol.t[1] - 1, 3.1, "Contraceptive threshold", fontsize=10, verticalalignment="bottom", horizontalalignment="left")

    ax_rep_dynamics.set_xticks(14:28:(28*3-14))
    ax_rep_dynamics.set_xticklabels([i for i in 1:3])

    ax_rep_dynamics.set_xlabel("Number of cycles")
    ax_rep_dynamics.set_ylabel("P₄ [ng/ml]")

    # Continue with the rest of the plot as before
    ax_rep_dynamics.plot(sol_wo_control.t, P₄_wo_control, color=:black, linewidth=1)
    ax_rep_dynamics.plot(sol.t, P₄, color="#E74C3c", linewidth=3)

    plt.tight_layout()

    fig_rep_dynamics.savefig("tmp3.png")

    res_400 = deserialize(joinpath("..", "..", "res", "RBA1.7_0.01", "N400", "20240927_Clock10.0_CCMO{NSGA2}_Gen940_Pop400", "result_variable.jls"))

    res_400_maxP4_vals = [child.f[1] for child in res_400.population]
    res_400_total_amount = [child.f[2] for child in res_400.population]
    
    fig_N400, ax_N400 = plt.subplots(figsize=(2,1.5))

    # Scatter of initial population
    ax_N400.scatter(res_400_total_amount, res_400_maxP4_vals, s=2)
    ax_N400.set_xticks([])
    ax_N400.set_yticks([])

    plt.tight_layout()

    fig_N400.savefig("tmp4.png")

    res_400_hyp = deserialize(joinpath("..", "..", "res", "RBA1.7_0.01", "N400_hyp", "20241006_Clock10.0_CCMO{NSGA2}_Gen850_Pop400", "result_variable.jls"))

    res_400_hyp_maxP4_vals = [child.f[1] for child in res_400_hyp.population]
    res_400_hyp_total_amount = [child.f[2] for child in res_400_hyp.population]
    
    fig_N400_hyp, ax_N400_hyp = plt.subplots(figsize=(2,1.5))

    # Scatter of initial population
    ax_N400_hyp.scatter(res_400_hyp_total_amount, res_400_hyp_maxP4_vals, s=2)
    ax_N400_hyp.set_xticks([])
    ax_N400_hyp.set_yticks([])

    plt.tight_layout()

    fig_N400_hyp.savefig("tmp5.png")

    res_500 = deserialize(joinpath("..", "..", "res", "RBA1.7_0.01", "N500", "20241002_Clock10.0_CCMO{NSGA2}_Gen1000_Pop500", "result_variable.jls"))

    res_500_maxP4_vals = [child.f[1] for child in res_500.population]
    res_500_total_amount = [child.f[2] for child in res_500.population]
    
    fig_N500, ax_N500 = plt.subplots(figsize=(2,1.5))

    # Scatter of initial population
    ax_N500.scatter(res_500_total_amount, res_500_maxP4_vals, s=2)
    ax_N500.set_xticks([])
    ax_N500.set_yticks([])

    plt.tight_layout()

    fig_N500.savefig("tmp6.png")

    res_500_hyp = deserialize(joinpath("..", "..", "res", "RBA1.7_0.01", "N500_hyp", "20241008_Clock10.0_CCMO{NSGA2}_Gen800_Pop500", "result_variable.jls"))

    res_500_hyp_maxP4_vals = [child.f[1] for child in res_500_hyp.population]
    res_500_hyp_total_amount = [child.f[2] for child in res_500_hyp.population]
    
    fig_N500_hyp, ax_N500_hyp = plt.subplots(figsize=(2,1.5))

    # Scatter of initial population
    ax_N500_hyp.scatter(res_500_hyp_total_amount, res_500_hyp_maxP4_vals, s=2)
    ax_N500_hyp.set_xticks([])
    ax_N500_hyp.set_yticks([])

    plt.tight_layout()

    fig_N500_hyp.savefig("tmp7.png")

    integrated_population = vcat(res_400.population, res_400_hyp.population, res_500.population, res_500_hyp.population)
    integrated_maxP4_vals = [child.f[1] for child in integrated_population]
    integrated_total_amount = [child.f[2] for child in integrated_population]
    
    integrated_population_morethan3 = [points for points in integrated_population if points.f[1]<3 && points.f[3]>3]
    integrated_morethan3_maxP4_vals = [child.f[1] for child in integrated_population_morethan3]
    integrated_morethan3_total_amount = [child.f[2] for child in integrated_population_morethan3]

    integrated_population_lessthan3 = [points for points in integrated_population if points.f[1]<3 && points.f[3]<3]
    integrated_lessthan3_maxP4_vals = [child.f[1] for child in integrated_population_lessthan3]
    integrated_lessthan3_total_amount = [child.f[2] for child in integrated_population_lessthan3]

    index_min = argmin(integrated_lessthan3_total_amount)

    fig_integrated, ax_integrated = plt.subplots(figsize=(8,6))

    # Scatter of initial population
    ax_integrated.scatter(integrated_total_amount, integrated_maxP4_vals, s=5)
    ax_integrated.scatter(integrated_morethan3_total_amount, integrated_morethan3_maxP4_vals, s=5, c="#c04f15")
    ax_integrated.scatter(integrated_lessthan3_total_amount, integrated_lessthan3_maxP4_vals, s=5, c="#8ed973")
    
    ax_integrated.scatter(integrated_lessthan3_total_amount[index_min], integrated_lessthan3_maxP4_vals[index_min], s=100, c=:magenta, marker="*")

    origin_ylim = ax_integrated.get_ylim()
    ax_integrated.vlines(x=integrated_lessthan3_total_amount[index_min], ymin=0, ymax=integrated_lessthan3_maxP4_vals[index_min], color=:magenta, linestyle=:dashed, linewidth=1.5)
    ax_integrated.set_ylim(origin_ylim)
    
    ax_integrated.axhline(y=3, color=:red, linestyle=:dashed)
    ax_integrated.text(0, 3.1, "Contraceptive threshold", fontsize=10, verticalalignment="bottom", horizontalalignment="left")
    ax_integrated.set_xlabel("Total amount", fontsize=12)
    ax_integrated.set_ylabel("Maximum P₄", fontsize=12)

    plt.tight_layout()

    fig_integrated.savefig("tmp8.png")

    fig_sum, ax_sum = plt.subplots(figsize=(16,6))

    administration_timings = 1:24
    optimal_dose_total = h5read(joinpath(res_directory, "Optimal_Dose.h5"), "E2_P4_1-24")

    # Highlight the morning time (6 AM - 12 PM)
    ax_sum.axvspan(6, 12, color="#fff68f", alpha=0.6, label="Morning (6 AM - 12 PM)")

    # Highlight the night time (6 PM - 12 AM)
    ax_sum.axvspan(18, 24, color="#D3D3D3", alpha=0.6, label="Night (6 PM - 12 AM)")

    # Draw the total near-optimal dosage for each administration timings
    ax_sum.plot(administration_timings, optimal_dose_total,
            linestyle=:solid, color="#696969", linewidth=5, 
            marker="o", markersize=13, 
            label="Optimal Dose")

    # Set axis labels and ticks
    ax_sum.set_xlabel("Administration timings [o'clock]", fontsize=12)
    ax_sum.set_ylabel("Near-optimal total dosage [μg]", fontsize=12)
    ax_sum.set_xticks([6, 12, 18, 24])

    # Morning text
    ylim = ax_sum.get_ylim()
    ax_sum.text(9, (ylim[1] + ylim[2]) / 2, "Morning", horizontalalignment="center", verticalalignment="center", fontsize=18, color=:black)

    # Night text
    ax_sum.text(21, (ylim[1] + ylim[2]) / 2, "Night", horizontalalignment="center", verticalalignment="center", fontsize=18, color=:black)

    ax_sum.scatter(administration_timings[10], optimal_dose_total[10], s=500, c=:magenta, marker="*", zorder=10)

    plt.tight_layout()

    fig_sum.savefig("tmp9.png")

    contdict_p = copy(p)
    contdict_p[end-2] = vcat(integrated_population_morethan3[1].x[1:21], zeros(7))
    contdict_p[end-1] = vcat(integrated_population_morethan3[1].x[22:end], zeros(7))

    tspan = (0, 10.0*28)
    x0_contdict, _, lags = setup_initial_conditions_time_info(contdict_p)
    prob_contdict = DDEProblem(
                menstrual_cycle_dynamics!,
                x0_contdict, history_func, tspan, contdict_p; 
                constant_lags=lags
            )
    sol_contdict = solve(prob_contdict, MethodOfSteps(RK4()), dtmax=0.05)
    _, P₄_contdict, _ = calculate_auxiliary(sol_contdict, contdict_p, sol_contdict.t)

    fig_contdict, ax_contdict = plt.subplots(figsize=(6,1.8))

    ax_contdict.plot(sol_contdict.t, P₄_contdict, color=:black, linewidth=1)

    ax_contdict.axhline(y=3, color=:red, linestyle=:dashed)
    ax_contdict.text(sol_contdict.t[1] - 1, 3.1, "Contraceptive threshold", fontsize=10, verticalalignment="bottom", horizontalalignment="left")

    ax_contdict.set_xticks(14:28:(28*10.0-14))
    ax_contdict.set_xticklabels([i for i in 1:10])

    ax_contdict.set_xlabel("Number of cycles")
    ax_contdict.set_ylabel("P₄ [ng/ml]")

    plt.tight_layout()

    fig_contdict.savefig("tmp10.png")

    min_p = copy(p)
    min_p[end-2] = vcat(integrated_population_lessthan3[index_min].x[1:21], zeros(7))
    min_p[end-1] = vcat(integrated_population_lessthan3[index_min].x[22:end], zeros(7))

    tspan = (0, 10.0*28)
    x0_min, _, lags = setup_initial_conditions_time_info(min_p)
    prob_min = DDEProblem(
                menstrual_cycle_dynamics!,
                x0_min, history_func, tspan, min_p; 
                constant_lags=lags
            )
    sol_min = solve(prob_min, MethodOfSteps(RK4()), dtmax=0.05)
    _, P₄_min, _ = calculate_auxiliary(sol_min, min_p, sol_min.t)

    fig_min, ax_min = plt.subplots(figsize=(6,1.8))

    ax_min.plot(sol_min.t, P₄_min, color=:black, linewidth=1)

    ax_min.axhline(y=3, color=:red, linestyle=:dashed)
    ax_min.text(sol_min.t[1] - 1, 2.5, "Contraceptive threshold", fontsize=10, verticalalignment="bottom", horizontalalignment="left")

    ax_min.set_xticks(14:28:(28*10.0-14))
    ax_min.set_xticklabels([i for i in 1:10])

    ax_min.set_xlabel("Number of cycles")
    ax_min.set_ylabel("P₄ [ng/ml]")

    plt.tight_layout()

    fig_min.savefig("tmp11.png")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end