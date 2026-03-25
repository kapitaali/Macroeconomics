# Appendix H: Numerical Algorithms Pseudocode

All pseudocode uses language-agnostic notation. APL equivalents are noted where distinctive.

---

## H.1 Newton–Raphson (Scalar)

```
Input: f, df, x0, tol, max_iter
x ← x0
for n = 1, 2, ..., max_iter:
    x_new ← x - f(x) / df(x)
    if |x_new - x| < tol: return x_new
    x ← x_new
return x  # (or raise: not converged)

APL: nr_step ← {⍵ - (f ⍵) ÷ df ⍵}
     x_star  ← nr_step ⍣ (tol∘>|⊢-nr_step) ⊢ x0
```

## H.2 Newton–Raphson (System)

```
Input: F: Rⁿ→Rⁿ, J_F: Rⁿ→Rⁿˣⁿ (Jacobian), x0, tol
x ← x0
for n = 1, 2, ..., max_iter:
    Δx ← solve(J_F(x), -F(x))   # linear system J_F * Δx = -F
    x_new ← x + Δx
    if ||x_new - x||∞ < tol: return x_new
    x ← x_new

APL: step ← {x - ((-F x) ⌹ J x)}
     x_star ← step ⍣ (tol∘>⌈/|⊢-step) ⊢ x0
```

## H.3 Value Function Iteration (VFI)

```
Input: u(c,s): utility, beta, a_grid (N), y_grid (S), Pi_y (S×S), tol
V ← zeros(N, S)   # initial guess

repeat:
    V_old ← copy(V)
    for each (i, j) in grid:
        c_choices ← feasible consumption at (a_grid[i], y_grid[j])
        a_next   ← (1+r)*(a_grid[i] + w*y_grid[j] - c_choices)
        EV       ← Pi_y[j,:] @ V_interp(a_next)    # expected next-period value
        payoffs  ← u(c_choices) + beta * EV
        V[i,j]   ← max(payoffs)
        policy[i,j] ← argmax(payoffs)
until ||V - V_old||∞ < tol
return V, policy

APL (inner loop): ⌈⌿ u_mat + beta × V_next_mat
     Iteration:  bellman ⍣ (tol∘>⌈/⌈/|⊢-bellman) ⊢ V0
```

## H.4 Howard Policy Improvement (PFI)

```
Input: u, beta, P_transition (N×N), a_grid, tol
policy ← initial_policy   # e.g., consume 90% of wealth

repeat:
    # Step 1: Policy evaluation
    u_policy ← u(a_grid - policy(a_grid))   # N-vector of utilities
    P_c      ← transition_matrix(policy)    # N×N Markov matrix
    V        ← solve((I - beta*P_c), u_policy)   # linear system

    # Step 2: Policy improvement
    policy_new ← argmax over c of {u(c) + beta * V_interp(a')}
until ||policy_new - policy||∞ < tol
return V, policy_new
```

## H.5 Gaussian Quadrature (Gauss–Hermite)

```
# For E[f(Y)] where Y ~ N(mu, sigma²)
nodes, weights ← hermgauss(n)    # n-point GH nodes and weights
y_nodes ← mu + sqrt(2) * sigma * nodes
E_f ← sum(weights * f(y_nodes)) / sqrt(pi)
return E_f

APL: E_f ← (gh_weights +.× f¨ mu+sig×(2*0.5)×gh_nodes) ÷ (○1)*0.5
```

## H.6 Runge–Kutta 4th Order (RK4)

```
# Solve ẋ = F(x, t) with step size h
rk4_step(F, x, t, h):
    k1 ← F(x,           t)
    k2 ← F(x + h/2*k1,  t + h/2)
    k3 ← F(x + h/2*k2,  t + h/2)
    k4 ← F(x + h*k3,    t + h)
    return x + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

APL: rk4 ← {F h x ← ⍵
     k1 ← F x ⋄ k2 ← F x+(h÷2)×k1 ⋄ k3 ← F x+(h÷2)×k2 ⋄ k4 ← F x+h×k3
     x + (h÷6)×k1+2×k2+2×k3+k4}
```

## H.7 Kalman Filter (Prediction + Update)

```
Input: F, H, Q, R (system matrices), a0, P0, y[1..T] (observations)

a ← a0;  P ← P0;  log_lik ← 0

for t = 1, ..., T:
    # Prediction step
    a_pred ← F @ a
    P_pred ← F @ P @ F.T + Q

    # Innovation
    v  ← y[t] - H @ a_pred
    Fv ← H @ P_pred @ H.T + R

    # Log-likelihood contribution
    log_lik -= 0.5 * (log_det(Fv) + v.T @ inv(Fv) @ v + p*log(2π))

    # Kalman gain
    K ← P_pred @ H.T @ inv(Fv)

    # Update step
    a ← a_pred + K @ v
    P ← (I - K @ H) @ P_pred

return log_lik, a, P   # final filtered state

APL: predict ← {a P F Q c ← ⍵ ⋄ (F+.×a)+c,  (F+.×P+.×⍉F)+Q}
     update  ← {a P H R d y ← ⍵
       v←y-(H+.×a)+d ⋄ Fv←(H+.×P+.×⍉H)+R
       K←P+.×(⍉H)+.×⌹Fv ⋄ a+(K+.×v), ((=⍨⍳≢a)-(K+.×H))+.×P}
```

## H.8 Metropolis–Hastings

```
Input: log_posterior, theta0, sigma_proposal, N_draws

theta ← theta0;  log_p ← log_posterior(theta)
chain ← zeros(N_draws, k);  n_accept ← 0

for s = 1, ..., N_draws:
    proposal ← theta + sigma_proposal * randn(k)
    log_p_prop ← log_posterior(proposal)
    log_alpha ← log_p_prop - log_p

    if log(rand()) < log_alpha:
        theta ← proposal;  log_p ← log_p_prop
        n_accept += 1

    chain[s, :] ← theta

return chain, n_accept / N_draws   # chain and acceptance rate

APL: mh_step ← {theta ← ⍵
     prop ← theta + scale × normal_draw k
     (logpost prop) > logpost theta + ⍟?0: prop ⋄ theta}
     chain ← mh_step \ N ⍴ theta0
```

## H.9 Krusell–Smith Outer Loop

```
# Approximate aggregation for HA models with aggregate shocks
Initialize: a0, a1, a2 in log K_{t+1} = a0 + a1*log K_t + a2*log A_t

repeat:
    # 1. Solve household problem given law of motion
    policy(a, y, K, A) ← EGM or VFI on (a,y,K,A) grid

    # 2. Simulate
    for t = 1, ..., T_sim:
        Draw A_t from Tauchen Markov chain
        Update all households: a_{i,t+1} = g*(a_{i,t}, y_{i,t}, K_t, A_t)
        K_{t+1} = mean(a_{i,t+1})   # aggregate capital

    # 3. Update law of motion via OLS regression
    (a0, a1, a2) ← OLS of log K_{t+1} on (1, log K_t, log A_t)
    R² ← compute from regression

until R² > 0.9999 and coefficients converged
return policy, (a0, a1, a2)
```

---

# Appendix I: Software Guides

---

## I.1 Dynare

**Installation.** Download from dynare.org; requires MATLAB (R2019b+) or Octave (6.0+). Add Dynare to path: `addpath('/path/to/dynare/matlab')`.

**Key commands:**

| Command | Purpose |
|---|---|
| `stoch_simul(order=1)` | Solve and simulate the model (1st order) |
| `stoch_simul(order=2)` | Second-order perturbation |
| `estimation(datafile=...)` | Bayesian estimation via MH |
| `identification` | Iskrev identification check |
| `shock_decomposition` | Structural shock decomposition |
| `forecast` | $h$-step-ahead forecasts |

**Standard .mod file blocks:**

```dynare
var y pi r;           // endogenous variables
varexo eps_u eps_r;   // exogenous shocks
parameters beta kappa sigma phi_pi phi_y rho_u rho_r;

// Assign parameter values
beta = 0.99; kappa = 0.15; sigma = 1.0;
phi_pi = 1.5; phi_y = 0.5; rho_u = 0.5; rho_r = 0.8;

model(linear);
  pi = beta*pi(+1) + kappa*y + eps_u;           // NKPC
  y  = y(+1) - sigma*(r - pi(+1) - rho_r*r(-1)); // DIS
  r  = phi_pi*pi + phi_y*y;                      // Taylor rule
end;

initval; pi=0; y=0; r=0; end;
shocks; var eps_u = 0.01^2; var eps_r = 0.01^2; end;

stoch_simul(order=1, irf=20);
```

**Accessing output from Dynare (MATLAB/APL bridge):**
- Decision rule: `oo_.dr.ghx` (state response), `oo_.dr.ghu` (shock response).
- IRFs: `oo_.irfs.y_eps_u` (IRF of $y$ to shock `eps_u`).
- Posterior: `oo_.posterior_mean.parameters`, `oo_.posterior_mode.parameters`.

**Common errors:**
- `NKPC is not a forward-looking equation`: check the `(+1)` notation.
- `Blanchard-Kahn conditions not satisfied`: verify Taylor principle $\phi_\pi > 1$.
- `Singular matrix in QZ`: check for redundant equations or unidentified parameters.

## I.2 Dyalog APL

**Installation.** Download Dyalog APL 18.2+ Unicode edition from dyalog.com. The free Community Edition is suitable for all exercises in this book.

**Book-wide session defaults** (set at the start of every workspace):
```apl
⎕IO ← 0    ⍝ Index origin 0: arrays are zero-indexed
⎕ML ← 1    ⍝ Migration level 1: standard dyadic ⊃, ⊢, ⊣ behavior
```

**APL Primitive Reference for Macroeconomics:**

| Primitive | Syntax | Purpose | Example |
|---|---|---|---|
| `+.×` | `A +.× B` | Matrix multiply | `G0 ⌹ G1` for $\Gamma_0^{-1}\Gamma_1$ |
| `⌹` | `b ⌹ A` | Solve $Ax=b$ or pseudoinverse | OLS: `Y ⌹ X` |
| `⍉` | `⍉A` | Transpose | `⍉Pi` |
| `∘.×` | `u ∘.× v` | Outer product | Covariance matrix |
| `⍤k` | `f⍤k` | Rank-$k$ application | Apply function to each row |
| `⍣` | `f⍣n` or `f⍣≡` | Power: apply $f$ $n$ times or to fixed point | VFI: `bellman⍣≡V0` |
| `/` | `+/v` | Reduce (sum) | `+/v` = sum; `×/v` = product |
| `\` | `+\v` | Scan (cumulative) | `+\v` = running total |
| `⍳n` | `⍳5` → `0 1 2 3 4` | Index generator | Grid: `a_grid ← a_min + h×⍳N` |
| `⍴` | `⍴A` or `n⍴v` | Shape / Reshape | `2 2⍴1 0 0 1` = identity |
| `⊂` | `⊂v` | Enclose (scalar of array) | Nested array element |
| `⊃` | `⊃L` | Disclose (first element) | Extract first of nested |
| `@` | `val @ idx ⊢ arr` | Selective update | Update element at index |
| `⌈⌿` | `⌈⌿M` | Column-wise maximum | Bellman: `⌈⌿payoff_matrix` |
| `⌊⌿` | `⌊⌿M` | Column-wise minimum | |
| `¨` | `f¨v` | Each (map function) | `u¨ c_grid` = utility at each point |

**APL coding patterns for macroeconomics:**

```apl
⍝ OLS / Linear system solve
B_hat ← (⌹ X) +.× Y

⍝ Lyapunov equation: vec(Σ) = (I - Φ⊗A)⁻¹ vec(Q)
kron  ← {⍺ ∘.× ⍵}
I4    ← =⍨ ⍳ 4
vec_Sigma ← (,Q) ⌹ I4 - (⍉Phi) kron A

⍝ Newton–Raphson fixed-point iteration
nr    ← {⍵ - (f ⍵) ÷ df ⍵}
xstar ← nr ⍣ (1e¯10∘>|⊢-nr) ⊢ x0

⍝ VFI Bellman operator
bellman ← {V_in←⍵ ⋄ ⌈⌿ u_mat + beta × V_next_mat V_in}
V_star  ← bellman ⍣ (1e¯6∘>⌈/⌈/|⊢-bellman) ⊢ V0

⍝ Tauchen outer product (core step)
diff  ← (⍉z_grid) ∘.- rho_A × z_grid    ⍝ N_A × N_A matrix
Pi    ← Phi_cdf_upper - Phi_cdf_lower    ⍝ normal CDF evaluated at diff±dz

⍝ IRF computation
irf_h ← {C⍣⍵ +.× D +.× shock} ¨ ⍳ H

⍝ HP filter
hp_filter ← {lam y ← ⍵ ⋄ T←≢y ⋄ I←=⍨⍳T
    D ← ... ⍝ second-difference matrix
    y - ⌹(I + lam × (⍉D)+.×D) +.× y}

⍝ MH chain
mh_step ← {theta scale ← ⍵
    prop ← theta + scale × rnorm ≢theta
    (logpost prop) > (logpost theta)+⍟?0: prop ⋄ theta}
chain ← mh_step ⍣ N ⊢ theta0   ⍝ or {mh_step ⍵ scale}\ N⍴theta0
```

**Connecting APL to Python** via `⎕PY`:
```apl
⎕PY.Import 'numpy as np'
result ← 'np.linalg.eigvals' ⎕PY.Call A    ⍝ eigenvalues of APL matrix A
normal_draws ← 'np.random.randn' ⎕PY.Call N ⍝ N standard normal draws
```

**Connecting APL to FRED** via `HttpCommand`:
```apl
⎕SE.SALT.Load 'HttpCommand'
url ← 'https://api.stlouisfed.org/fred/series/observations?series_id=GDPC1&api_key=YOUR_KEY&file_type=json'
data ← (HttpCommand.Get url).Data
```

**Workspace organization:**
```
#.Models     ← namespace for DSGE/RBC models
#.Data       ← namespace for loaded datasets
#.Util       ← namespace for numerical utilities (RK4, bisection, etc.)
```

Save workspace: `⎕SAVE 'my_dsge'`. Save source to file: `]SAVE Models/nk_model.aplf`.

## I.3 Python

**Core stack:** NumPy (arrays), SciPy (optimization, integration, sparse), matplotlib (plotting), statsmodels (time series), sklearn (ML).

**Key packages for macroeconomics:**

| Package | Purpose | Install |
|---|---|---|
| `quantecon` | Markov chains, VFI, Tauchen | `pip install quantecon` |
| `pydsge` | DSGE model solution | `pip install pydsge` |
| `statsmodels` | VAR, ARIMA, Kalman filter | `pip install statsmodels` |
| `pykalman` | Kalman filter and smoother | `pip install pykalman` |
| `fredapi` | FRED data access | `pip install fredapi` |
| `SALib` | Sobol sensitivity analysis | `pip install SALib` |
| `dynare` | Python-Dynare bridge | See dynare.org |

**Solving linear systems:** `scipy.linalg.solve(A, b)` (LU); `scipy.linalg.lstsq(X, y)` (OLS/QR).

**VAR estimation (statsmodels):**
```python
from statsmodels.tsa.api import VAR
model = VAR(data)
results = model.fit(maxlags=4, ic='bic')
irfs = results.irf(20)
irfs.plot(orth=True)
```

**Publication-quality figures (matplotlib):**
```python
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':11, 'figure.dpi':150,
                     'axes.spines.top':False, 'axes.spines.right':False})
```

## I.4 Julia

**Installation.** Download from julialang.org. Current LTS: Julia 1.10.

**Key packages:**

| Package | Purpose | Add command |
|---|---|---|
| `DifferentialEquations.jl` | ODE/SDE solvers | `]add DifferentialEquations` |
| `QuantEcon.jl` | VFI, Markov chains, Tauchen | `]add QuantEcon` |
| `Optim.jl` | Unconstrained/constrained optimization | `]add Optim` |
| `StateSpaceModels.jl` | Kalman filter, structural models | `]add StateSpaceModels` |
| `GLMNet.jl` | LASSO/Ridge | `]add GLMNet` |
| `ForwardDiff.jl` | Automatic differentiation | `]add ForwardDiff` |
| `LinearAlgebra` | (stdlib) LU, QR, eigen, SVD | (built-in) |
| `GlobalSensitivity.jl` | Sobol indices | `]add GlobalSensitivity` |
| `NLsolve.jl` | Nonlinear system solving | `]add NLsolve` |

**Performance tips:**
- First run is slow (JIT compilation); subsequent runs are fast.
- `@code_warntype f(x)` to check for type instabilities.
- Use `BenchmarkTools.@btime` to measure performance.
- Pre-allocate arrays with `zeros(n)` rather than building with `push!`.

## I.5 R

**Key packages:**

| Package | Purpose |
|---|---|
| `vars` | VAR estimation, IRFs, Granger causality |
| `KFAS` | Kalman filter and smoother (state-space) |
| `gmm` | GMM estimation |
| `nloptr` | Nonlinear optimization |
| `deSolve` | ODE solvers |
| `bvars` | BVAR with Minnesota prior |
| `glmnet` | LASSO/Ridge/elastic net |
| `sensitivity` | Sobol sensitivity analysis |
| `nleqslv` | Nonlinear system solving |
| `fredr` | FRED API access |

**BVAR with Minnesota prior (bvars):**
```r
library(bvars)
data <- ts(cbind(gdp_growth, inflation, rate), frequency=4)
fit  <- bvar(data, lags=4, prior=set_prior(type="minnesota", lambda=0.2))
irfs <- irf(fit, n.ahead=20, identification="chol")
plot(irfs)
```

**Kalman filter (KFAS):**
```r
library(KFAS)
model <- SSModel(y ~ SSMtrend(2, Q=list(NA, NA)), H=NA)
fit   <- fitSSM(model, inits=log(c(0.01, 0.001, 0.05)))
out   <- KFS(fit$model)
```
