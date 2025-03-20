Seed: 1 (default), 2024 several times
tol: setting 1e-4 (both x_tol and f_tol) as default

****** From RBA results
Among the cases where P4 stays below 3.0 over three cycles (days 0–84), select the one that uses the smallest amount. Then, for each result, check the maximum P4 value over 10 cycles and record the five results that used the least amount. 

****** hyp
x_tol = 1e-4, f_tol =1e-4, maximum iteration = 10000,  f_call_limits = Inf
N=400 or 500,
seed = 1, rng = Random.MersenneTwister(1),
Algorithm = CCMO(NSGA2())

non-hyp result: Default setting(
    η_cr = 20,
    p_cr = 0.9,
    η_m = 20,
    p_m = 1.0 / D,
)
Because our problem has 42 varibles, p_m = 1.0/42 ~= 0.0238

hyp result: Setting with more variation (
    η_cr = 10,
    p_cr = 0.9,
    η_m = 10,
    p_m = 0.1,
)
