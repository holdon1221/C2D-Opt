## General Settings
- **Seed**: 1 (default)

- **Tolerance settings**:
   Both `x_tol` and `f_tol` were set to `1e-4` to ensure convergence with moderate precision.

---

## Results
- From the results across **three full menstrual cycles (days 0–84)**, we filtered the solutions where the maximum P4 concentration remained **below 3.0 ng/mL** throughout the entire period.
- Among those valid solutions, we identified the **single solution with the minimum total hormone usage.**
- For this subset of minimal-use solutions, we then **evaluated the long-term stability** by checking the **maximum P4 value over 10 cycles (i.e., days 0–280).**
- We selected the **top five solutions** from this extended set based on **lowest total hormone usage.**

---

## Optimization Settings

### Base configuration
- `x_tol=1e-4, f_tol=1e-4`
- **Maximum iterations:** 10_000
- **Function call limit:** Unbounded
- **Population size (N):** 400 or 500
- **Random number generator:** `Random.MersenneTwister(1), seed=1`
- **Algorithm:** `CCMO(NSGA2())`

### Genetic Operator Settings (default)
- Crossover distribution index (`η_cr`) = 20
- Crossover probability (`p_cr`) = 0.9
- Mutation distribution index (`η_m`) = 20
- Mutation probability (`p_m`) = 1/D = 1/42 ≈ 0.0238

### Hyperparameter-Tuned Configuration (more exploration)
- Crossover distribution index (`η_cr`) = 10
- Crossover probability (`p_cr`) = 0.9
- Mutation distribution index (`η_m`) = 10
- Mutation probability (`p_m`) = 0.1

This setting increases search diversity during evolution by setting lower distribution indices (`η_cr, η_m`), which allows offspring to deviate more significantly from parents, and by using a higher mutation probability (`p_m`), encouraging broader exploration and helping to avoid premature convergence.