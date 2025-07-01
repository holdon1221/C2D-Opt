# C2D-Opt: Circadian & Cycle-based Dosing Optimization

This repository contains all the codes used in the manuscript:

> **"Optimizing oral contraceptive timing: Morning intake reduces doses and enhances efficacy"**

We build upon the model proposed by Gavina et al. to simulate and optimize oral contraceptive regimens, integrating circadian rhythm and pharmacokinetics/pharmacodynamics (PK/PD).

---

## 📁 Directory Structure

- `constant/`: MATLAB code for parameter estimation
  - `Cosine fitting/`: Circadian hormone data fitting
  - `Parameter estimation circadian/`: Circadian parameter estimation
  - `Parameter estimation drug PK/`: PK parameter fitting
- `src/julia_codes/`: Julia code for multi-objective optimization
- `res/`: Precomputed results and configuration notes

---

## 🚀 Usage Guide

### A. Parameter Estimation (MATLAB)

1. **Cosine fitting**  
  
  ```matlab
  cd constant/Cosine fitting
  run('cosine_fit.m')
  ``` 

2. **Estimate circadian parameters**

  ```matlab
  cd ../Parameter estimation circadian
  run('circ_par_est.m')
  ```

3. **Estimate PK parameters**

  ```matlab
  cd ../Parameter estimation drug PK
  run('PK_est.m')
  ```

### B. Model Simulation (MATLAB)

  ```matlab
  cd constant
  run('main_model.m')
  ```

### C. Optimization (Julia)

1. **Navigate to:**

  ```cmd
  cd src/julia_codes
  ```

2. **Run:**

  ```cmd
  bash run_optimization_clock.sh [1|2|3|4] [recipient_email]
  ```
This script performs time-of-day optimizations (from 1:00 to 24:00) using a multi-objective NSGA-II algorithm implemented in Metaheuristics.jl.
- `MODE=1`: Run all 24 optimizations sequentially
- `MODE=2`: Parallelize into odd vs even groups
- `MODE=3`: Parallelize into 3 groups based on `mu % 3`
- `MODE=4`: Parallelize into 4 groups based on `mu % 4`

*Note: Email notifications are supported on macOS only.*

### ⚙️ Optimization Configuration
See res/ReadMe.txt for detailed hyperparameters.
Final results reflect a combination of optimization runs with the full range of hyperparameter settings shown below. 

*Note: While results are stable across platforms, strict numerical reproducibility (identical outputs across runs) is currently guaranteed only on Windows due to OS-level floating-point behavior.*

- **Seed**: 1
- **Population**: 400 or 500
- **Algorithm**: CCMO(NSGA2)
- **Stopping criteria**: `x_tol = 1e-4, f_tol = 1e-4, max_iter = 10000, f_call_limits = Inf`
- **Mutation rates**:
  - Default: η_cr = 20, p_cr = 0.9, η_m=20, p_m = 1.0/42
  - Hyperparameter variation: η_cr = 10, p_cr = 0.9, η_m=10, p_m = 0.1

---

## 📦 Dependencies

### MATLAB
- Standard MATLAB (tested on R2023+)

### Julia
- [Julia v1.10.5](https://julialang.org/downloads/oldreleases/#:~:text=v1.10.5%2C%20on%202024%2D08%2D28T15%3A43%3A13Z)
- **Package management:** via `Project.toml` and `Manifest.toml`

To set up the environment:

```bash
cd src/julia_codes
julia --project -e 'using Pkg; Pkg.instantiate()'
```

### 📌 Main packages include:
- [`Metaheuristics.jl v3.3.5`](https://jmejia8.github.io/Metaheuristics.jl/v3.3.5/)
- `CSV.jl v0.10.15`
- `DataFrames.jl v1.7.0`
- `Plots.jl v1.40.9`

And more (see Project.toml / Manifest.toml for full list)

📖 For Julia’s official reproducibility mechanism using Manifest.toml, see [this](https://pkgdocs.julialang.org/v1/environments/#manifests).

---

## 📝 Citation (filled later)
If you use this code, please cite the manuscript:

Optimizing oral contraceptive timing: Morning intake reduces doses and enhances efficacy (preprint or journal info here)

## 📧 Contact (filled later)
For questions, please contact the corresponding author.
