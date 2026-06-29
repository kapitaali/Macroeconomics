# Chapter 17: Recursive Methods for Real Business Cycle Models

*Value Function Iteration over the (k, A) State Space*

> *"The RBC model is a quantitative theory machine: feed in calibrated parameters, compute the equilibrium policy functions, simulate the model, and compare the resulting second moments to the data."*
> — Finn Kydland and Edward Prescott

**Cross-reference:** *Principles* Ch. 27 (RBC model: specification, calibration, stylized facts, criticisms) **[P:Ch.27]**

---

## 17.1 The RBC Model as a Recursive Problem

The **Real Business Cycle (RBC)** model of Kydland and Prescott (1982) and Long and Plosser (1983) is the discrete-time stochastic extension of the neoclassical growth model. It adds:

1. **Stochastic TFP:** Technology $A_t$ follows an AR(1) process in logs.
2. **Labor-leisure choice:** Households choose hours worked as well as consumption and saving.
3. **Competitive equilibrium:** Prices are determined in competitive markets; the equilibrium can be decentralized.

The key insight of Kydland and Prescott: the competitive equilibrium of the RBC model is equivalent to the solution of a social planner's problem (by the First Welfare Theorem, since there are no externalities or market failures). This means we can solve the planner's problem using the dynamic programming tools of Chapter 15 and interpret the solution as the competitive equilibrium.

**The state variables** of the RBC model are $(k_t, A_t)$: the current capital stock and the current TFP level. Both are known at the beginning of period $t$; both are inherited from the past. The control variables are $(c_t, n_t)$: consumption and hours worked chosen at the beginning of period $t$.

---

## 17.2 Model Specification

### 17.2.1 Preferences

The representative household maximizes:

$$\mathbb{E}_0\sum_{t=0}^\infty\beta^t u(c_t, 1-n_t), \quad u(c, 1-n) = \frac{[c^\mu(1-n)^{1-\mu}]^{1-\sigma}}{1-\sigma},$$

where $n_t \in [0,1]$ is hours worked, $1-n_t$ is leisure, $\mu \in (0,1)$ is the consumption share in the composite good, and $\sigma > 0$ is the curvature parameter (the coefficient of relative risk aversion over the composite good).

This King-Plosser-Rebelo (1988) utility specification is consistent with balanced growth: with log-utility in the composite good ($\sigma = 1$), the hours-worked ratio $n^*$ is constant along the balanced growth path — consistent with the Kaldor fact of stable labor supply [P:Ch.5.1].

### 17.2.2 Production and Capital Accumulation

$$Y_t = A_t K_t^\alpha n_t^{1-\alpha}, \quad K_{t+1} = (1-\delta)K_t + I_t, \quad Y_t = C_t + I_t.$$

### 17.2.3 TFP Process

$$\ln A_t = \rho_A\ln A_{t-1} + \varepsilon_t^A, \quad \varepsilon_t^A \sim \mathcal{N}(0, \sigma_A^2),$$

with $\rho_A \in (0,1)$ and $\sigma_A > 0$. The TFP shock $\varepsilon_t^A$ is the fundamental driver of business cycles in the RBC model.

### 17.2.4 The Planner's Problem

The social planner maximizes household utility subject to feasibility and the TFP process. The state at the beginning of period $t$ is $\mathbf{s}_t = (k_t, A_t)$. The planner's Bellman equation:

$$V(k, A) = \max_{c, n, k'}\left\{u(c, 1-n) + \beta\mathbb{E}_A[V(k', A')]\right\}$$

subject to:
$$c + k' = AK^\alpha n^{1-\alpha} + (1-\delta)k, \quad c \geq 0, \quad n \in [0,1], \quad k' \geq 0.$$

The expectation $\mathbb{E}_A[V(k', A')]$ is over the next period's TFP realization, conditional on the current $A$.

---

## 17.3 Value Function Iteration over the $(k, A)$ Grid

### 17.3.1 State Space Discretization

**Capital grid:** $\mathcal{K} = \{k_1, \ldots, k_{N_k}\}$ with $k_1 = k_{min} \approx (1-\varepsilon)k^*$ and $k_{N_k} = k_{max} \approx (1+\varepsilon)k^*$ for some $\varepsilon > 0$ (e.g., $\varepsilon = 0.5$ means ±50% around the steady state).

**TFP grid:** $\mathcal{A} = \{A_1, \ldots, A_{N_A}\}$ discretized from the AR(1) using Tauchen (1986) with $N_A$ points and transition matrix $\Pi \in \mathbb{R}^{N_A\times N_A}$.

The value function is represented as an $N_k \times N_A$ matrix $\mathbf{V}$, with $V_{ij} \approx V(k_i, A_j)$.

### 17.3.2 The Bellman Operator in Matrix Form

For each $(k_i, A_j)$:

1. Compute output: $Y_{ij}^n = A_j k_i^\alpha n^{1-\alpha}$ for each feasible $n$.
2. Compute consumption: $c = Y_{ij}^n + (1-\delta)k_i - k'$ for each feasible $(n, k')$.
3. Compute utility: $u(c, 1-n)$ if $c > 0$, else $-\infty$.
4. Compute expected next-period value: $\mathbb{E}_j[V(k')] = \sum_{j'}\Pi_{jj'}V_{k',j'}$.
5. Maximize over $(n, k')$ to get $V^{new}_{ij}$.

**Theorem 17.1 (VFI Convergence for the RBC Model).** Under standard regularity conditions on $u$ (continuous, strictly concave, satisfying Inada conditions), the Bellman operator for the RBC model is a contraction with modulus $\beta$. VFI converges to the unique fixed point.

*Proof.* Same structure as Theorem 15.2, extended to the two-dimensional state space. The key step: Blackwell's monotonicity and discounting conditions are satisfied because $u$ is increasing in $c$, $V$ is evaluated at the same $\beta$ in each iteration, and the constraint set is compact. $\square$

### 17.3.3 Simplification: Analytical Labor Optimality

For the KPR utility specification with $\sigma = 1$ (log in composite):

$$u(c, 1-n) = \mu\ln c + (1-\mu)\ln(1-n),$$

the labor-leisure optimality condition is:

$$\frac{(1-\mu)}{1-n_t} = \frac{\mu}{c_t}(1-\alpha)\frac{Y_t}{n_t} = \mu\frac{(1-\alpha)A_t k_t^\alpha n_t^{-\alpha}}{c_t}.$$

This pins down $n_t$ as a function of $(k_t, A_t, c_t)$. Substituting into the resource constraint, we can eliminate $n$ and solve the remaining $(k, c)$ problem — reducing the dimension of the optimization by one. This computational trick speeds up VFI substantially.

For general $\sigma$, we solve for $n$ numerically at each grid point using the first-order condition above.

---

## 17.4 Calibration: The Kydland–Prescott Method

**Definition 17.1 (Calibration).** **Calibration** is the process of assigning numerical values to model parameters by: (i) matching parameters that have direct empirical counterparts (e.g., $\alpha$ = capital share in national income, $\delta$ = depreciation rate); (ii) setting preference parameters to match long-run ratios (e.g., $\mu$ set so that $n^* = 0.33$, one-third of time at work).

**Standard RBC calibration (quarterly data):**

| Parameter | Value | Source |
|---|---|---|
| $\alpha$ (capital share) | 0.36 | National income accounts |
| $\delta$ (depreciation) | 0.025 per quarter | BEA fixed asset accounts |
| $\beta$ (discount factor) | 0.99 | Matches $r^* \approx 4\%$ annual |
| $\sigma$ (CRRA) | 1.0 | Log utility |
| $\mu$ (consumption share) | 0.357 | Matches $n^* = 1/3$ |
| $\rho_A$ (TFP persistence) | 0.95 | Solow residual AR(1) |
| $\sigma_A$ (TFP volatility) | 0.0072 | Solow residual std dev |

**Steady-state relations used in calibration:**

From the Euler equation at the steady state: $1 = \beta(1 + r^* - \delta) \Rightarrow r^* = 1/\beta - 1 + \delta$.

From the capital-output ratio: $K^*/Y^* = \alpha/(r^*+\delta)$.

From labor market clearing: $(1-\mu)/\mu = (1-\alpha)Y^*/c^* \cdot (1-n^*)/n^*$, pins down $\mu$ given $n^* = 1/3$.

---

## 17.5 The HP Filter and Business Cycle Moments

Macroeconomic data contain trends (long-run growth) and cycles (short-run fluctuations). The **Hodrick–Prescott filter** [P:Ch.6.3] separates them by solving:

$$\min_{\{\tau_t\}} \sum_{t=1}^T(y_t - \tau_t)^2 + \lambda\sum_{t=2}^{T-1}[(\tau_{t+1}-\tau_t) - (\tau_t-\tau_{t-1})]^2,$$

where $y_t$ is log output, $\tau_t$ is the trend component, and $\lambda = 1600$ (for quarterly data). The cycle is $c_t^{HP} = y_t - \tau_t$.

In matrix form, the solution is $\bm\tau = (I + \lambda D'D)^{-1}\mathbf{y}$, where $D$ is the second-difference matrix. In APL, this is:

```apl
⍝ APL — HP filter via matrix inversion
⎕IO←0 ⋄ ⎕ML←1

hp_filter ← {lambda y ← ⍵
    T ← ≢y
    ⍝ Second-difference matrix D (T-2 × T)
    D ← {(T-2) T ⍴ 0}   ⍝ build D properly
    ⍝ For each row i: D[i,i]=1, D[i,i+1]=-2, D[i,i+2]=1
    ⍝ Use outer product and diagonal positioning
    rows ← ⍳T-2
    D ← (T-2) T ⍴ 0
    D[rows; rows]   +← 1
    D[rows; rows+1] +← -2
    D[rows; rows+2] +← 1
    ⍝ Filter: τ = (I + λ D'D)⁻¹ y
    I  ← =⍨⍳T
    tau ← ⌹(I + lambda × (⍉D) +.× D) +.× y
    y - tau}   ⍝ return cycle component
```

**Target second moments (U.S. quarterly data, HP-filtered, 1954–2019):**

| Variable | Std dev (% of Y std) | Correlation with Y | Autocorrelation |
|---|---|---|---|
| Output $Y$ | 1.72% | 1.00 | 0.87 |
| Consumption $C$ | 0.86% (0.50) | 0.88 | 0.88 |
| Investment $I$ | 5.30% (3.08) | 0.93 | 0.87 |
| Hours $n$ | 1.59% (0.93) | 0.88 | 0.88 |
| Real wage $w$ | 0.68% (0.40) | 0.12 | 0.66 |

The RBC model must match these moments — particularly the relative volatility of consumption (less volatile than output) and investment (more volatile), and the high contemporaneous correlation of all variables with output.

---

## 17.6 Full VFI Implementation

```apl
⍝ APL — Full RBC model VFI over (k, A) grid
⎕IO←0 ⋄ ⎕ML←1

⍝ Parameters
alpha←0.36  ⋄  delta←0.025  ⋄  beta←0.99  ⋄  sigma←1  ⋄  mu←0.357
rho_A←0.95  ⋄  sig_A←0.0072

⍝ Steady state
r_star  ← (÷beta) - 1 + delta
k_y_rat ← alpha ÷ r_star + delta
n_star  ← 1÷3               ⍝ calibration target
⍝ With σ=1 (log), y* given by other conditions; normalize A*=1
A_star  ← 1
kstar   ← (alpha×A_star÷r_star+delta)*÷1-alpha    ⍝ from MPK=r*+δ

⍝ Grids
N_k ← 50  ⋄  N_A ← 7
k_grid ← kstar × 0.5 + (⍳N_k) × ÷N_k   ⍝ [0.5k*, 1.5k*]

⍝ Tauchen discretization for log(A)
⍝ (Use ⎕PY to call Python for normal CDF in production)
sigma_z ← sig_A ÷ (1-rho_A*2)*0.5
z_min   ← ¯3×sigma_z  ⋄  z_max ← 3×sigma_z
z_grid  ← z_min + (z_max-z_min) × (⍳N_A) ÷ N_A-1
A_grid  ← *z_grid    ⍝ A = exp(log_A)

⍝ Transition matrix Pi (simplified: diagonal dominant)
⍝ In production: use full Tauchen formula with normal CDF via ⎕PY
Pi ← N_A N_A ⍴ 0.07   ⍝ placeholder; replace with proper Tauchen
Pi[⍳N_A; ⍳N_A] ← 0.60  ⍝ diagonal

⍝ Production: Y(k,n,A) = A*k^alpha*n^(1-alpha)
prod ← {A k n ← ⍵ ⋄ A×(k*alpha)×n*1-alpha}

⍝ Utility: log composite for σ=1
util ← {c n ← ⍵ ⋄ (mu×⍟c) + (1-mu)×⍟1-n}

⍝ Optimal labor: from FOC (1-mu)/c = mu*(1-alpha)*Y/(n*c)
⍝ => n = (1-alpha)*(1-mu)/((1-alpha)*(1-mu) + mu*(1-n)... implicit
⍝ For log utility, n* = 1-mu in frictionless model (simplified)
n_opt ← 1-mu    ⍝ constant labor for log-log (special case)

⍝ Initialize value function
V ← N_k N_A ⍴ 0    ⍝ N_k × N_A matrix

⍝ Bellman operator: one application
bellman_rbc ← {V_in ← ⍵
    V_new ← N_k N_A ⍴ 0
    ⍝ For each (k_i, A_j): find optimal k' ∈ k_grid
    ⍝ Expected next-period value: EV[k'] = V[k',*] +.× Pi[j,*]
    EV ← V_in +.× ⍉Pi    ⍝ N_k×N_A: EV[i,j] = sum_j' Pi[j,j'] V[i,j']
    ⍝ Resource constraint: c = A*k^α*n^(1-α) + (1-δ)*k - k'
    ⍝ For each k_i and A_j and candidate k':
    ⍝ Build payoff tensor using outer operations
    ⍝ (Simplified: evaluate at n=n_opt)
    :For j :In ⍳N_A
        Aj ← A_grid[j]
        :For i :In ⍳N_k
            ki   ← k_grid[i]
            Y_ij ← prod Aj ki n_opt
            x_ij ← Y_ij + (1-delta)×ki    ⍝ cash on hand
            c_kp ← x_ij - k_grid          ⍝ consumption for each k'
            feas ← c_kp > 0               ⍝ feasibility
            u_kp ← feas × util¨ (0⌈c_kp) n_opt
            u_kp +← (1-feas)×¯1e20       ⍝ penalize infeasible
            payoff ← u_kp + beta × EV[;j] ⍝ + discounted EV at each k'
            V_new[i;j] ← ⌈/ payoff
        :EndFor
    :EndFor
    V_new}

⍝ Iterate to convergence
V_star ← bellman_rbc ⍣ (1e¯5∘>⌈/⌈/|⊢-bellman_rbc) ⊢ V
```

```python
import numpy as np
from scipy.interpolate import RegularGridInterpolator

# RBC model VFI
alpha, delta, beta, mu = 0.36, 0.025, 0.99, 0.357
rho_A, sig_A = 0.95, 0.0072

# Steady state
r_star = 1/beta - 1 + delta
kstar  = (alpha / (r_star + delta))**(1/(1-alpha))
n_star = 1/3

# Grids
Nk, NA = 100, 7
k_grid = np.linspace(0.5*kstar, 1.5*kstar, Nk)

# Tauchen discretization for log(A)
from scipy.stats import norm
sig_z = sig_A / np.sqrt(1-rho_A**2)
z_grid = np.linspace(-3*sig_z, 3*sig_z, NA)
A_grid = np.exp(z_grid)
dz = z_grid[1]-z_grid[0]
Pi = np.zeros((NA, NA))
for i in range(NA):
    mu_z = rho_A * z_grid[i]
    Pi[i,0]    = norm.cdf((z_grid[0]+dz/2-mu_z)/sig_A)
    Pi[i,NA-1] = 1 - norm.cdf((z_grid[NA-1]-dz/2-mu_z)/sig_A)
    for j in range(1, NA-1):
        Pi[i,j] = norm.cdf((z_grid[j]+dz/2-mu_z)/sig_A) - norm.cdf((z_grid[j]-dz/2-mu_z)/sig_A)

n_opt = n_star  # constant labor (log-log utility, simplified)

def bellman_rbc(V):
    V_new = np.full((Nk, NA), -1e20)
    c_opt = np.zeros((Nk, NA))
    k_opt = np.zeros((Nk, NA))
    
    # Expected value: EV[i,j] = sum_j' Pi[j,j'] V[i,j']
    EV = V @ Pi.T  # shape: Nk × NA
    V_interp = RegularGridInterpolator((k_grid, A_grid), EV, method='linear', bounds_error=False, fill_value=None)
    
    for j, Aj in enumerate(A_grid):
        for i, ki in enumerate(k_grid):
            Y_ij = Aj * ki**alpha * n_opt**(1-alpha)
            x_ij = Y_ij + (1-delta)*ki
            
            # Vectorize over k' choices
            k_choices = k_grid[k_grid < x_ij]
            if len(k_choices) == 0: continue
            c_choices = x_ij - k_choices
            
            # Utility
            util_c = mu*np.log(c_choices) + (1-mu)*np.log(1-n_opt)
            
            # Expected continuation
            points = np.column_stack([k_choices, np.full(len(k_choices), Aj)])
            EV_kp = V_interp(points)
            
            payoffs = util_c + beta * EV_kp
            idx = np.argmax(payoffs)
            V_new[i, j] = payoffs[idx]
            c_opt[i, j] = c_choices[idx]
            k_opt[i, j] = k_choices[idx]
    
    return V_new, c_opt, k_opt

# Initialize and run VFI
V = np.zeros((Nk, NA))
for iteration in range(300):
    V_new, c_opt, k_opt = bellman_rbc(V)
    diff = np.max(np.abs(V_new - V))
    V = V_new
    if iteration % 20 == 0:
        print(f"Iter {iteration:3d}: max diff = {diff:.2e}")
    if diff < 1e-5:
        print(f"Converged at iteration {iteration}")
        break

# Simulate the model
T_sim = 10000
np.random.seed(42)
eps = np.random.normal(0, sig_A, T_sim)
log_A = np.zeros(T_sim); log_A[0] = 0
for t in range(1, T_sim):
    log_A[t] = rho_A * log_A[t-1] + eps[t]
A_sim = np.exp(log_A)

k_sim = np.zeros(T_sim); k_sim[0] = kstar
Y_sim = np.zeros(T_sim); C_sim = np.zeros(T_sim)

c_interp = RegularGridInterpolator((k_grid, A_grid), c_opt, bounds_error=False, fill_value=None)
for t in range(T_sim-1):
    pt = np.array([[k_sim[t], np.clip(A_sim[t], A_grid[0], A_grid[-1])]])
    C_sim[t] = float(c_interp(pt))
    Y_sim[t] = A_sim[t] * k_sim[t]**alpha * n_opt**(1-alpha)
    I_t = Y_sim[t] - C_sim[t]
    k_sim[t+1] = (1-delta)*k_sim[t] + I_t

# HP filter and compute second moments
from scipy.signal import lfilter

def hp_filter(y, lam=1600):
    T = len(y)
    # Build second difference matrix via solving the linear system
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve
    D = diags([1, -2, 1], [0, 1, 2], shape=(T-2, T)).toarray()
    I_mat = np.eye(T)
    tau = np.linalg.solve(I_mat + lam * D.T @ D, y)
    return y - tau

log_Y = np.log(Y_sim[100:]); log_C = np.log(C_sim[100:])
cyc_Y = hp_filter(log_Y); cyc_C = hp_filter(log_C)

print(f"\nModel second moments (HP-filtered):")
print(f"  Std(Y): {100*np.std(cyc_Y):.3f}%")
print(f"  Std(C)/Std(Y): {np.std(cyc_C)/np.std(cyc_Y):.3f}  (data: 0.50)")
print(f"  Corr(C,Y): {np.corrcoef(cyc_C, cyc_Y)[0,1]:.3f}  (data: 0.88)")
```

```julia
# Julia — RBC model: Tauchen + VFI essentials
using Statistics, Distributions, Interpolations

alpha, delta, beta, mu = 0.36, 0.025, 0.99, 0.357
rho_A, sig_A = 0.95, 0.0072
r_star = 1/beta - 1 + delta
kstar  = (alpha/(r_star+delta))^(1/(1-alpha))
n_opt  = 1/3  # simplified constant labor

Nk, NA = 80, 7
k_grid = range(0.5*kstar, 1.5*kstar, length=Nk) |> collect
sig_z  = sig_A/sqrt(1-rho_A^2)
z_grid = range(-3sig_z, 3sig_z, length=NA) |> collect
A_grid = exp.(z_grid)

# Tauchen transition matrix
d = Normal(0, sig_A); dz = step(range(-3sig_z,3sig_z,length=NA))
Pi = zeros(NA, NA)
for i in 1:NA
    mu_z = rho_A*z_grid[i]
    Pi[i,1]    = cdf(d, z_grid[1]+dz/2-mu_z)
    Pi[i,NA]   = 1 - cdf(d, z_grid[NA]-dz/2-mu_z)
    for j in 2:NA-1
        Pi[i,j] = cdf(d, z_grid[j]+dz/2-mu_z) - cdf(d, z_grid[j]-dz/2-mu_z)
    end
end

util(c, n) = mu*log(max(c,1e-10)) + (1-mu)*log(max(1-n,1e-10))

V = zeros(Nk, NA)
println("Running RBC VFI...")
for iter in 1:200
    EV   = V * Pi'   # Nk × NA: expected continuation
    V_new = fill(-Inf, Nk, NA)
    for j in 1:NA, i in 1:Nk
        Aj, ki = A_grid[j], k_grid[i]
        x = Aj*ki^alpha*n_opt^(1-alpha) + (1-delta)*ki
        best = -Inf
        for ip in 1:Nk
            kp = k_grid[ip]; c = x - kp
            c > 0 || continue
            pay = util(c, n_opt) + beta*EV[ip,j]
            pay > best && (best = pay)
        end
        V_new[i,j] = best
    end
    diff = maximum(abs.(V_new - V)); V .= V_new
    iter % 25 == 0 && println("  iter=$iter, diff=$(round(diff,digits=6))")
    diff < 1e-5 && (println("Converged at iter $iter"); break)
end
println("VFI complete. V range: [$(round(minimum(V),digits=2)), $(round(maximum(V),digits=2))]")
```

```r
# R — RBC calibration and moment matching
alpha<-0.36; delta<-0.025; beta<-0.99; rho_A<-0.95; sig_A<-0.0072
r_star <- 1/beta - 1 + delta
kstar  <- (alpha/(r_star+delta))^(1/(1-alpha))
n_opt  <- 1/3

cat(sprintf("Steady state: k*=%.3f, r*=%.3f (annual %.1f%%)\n",
            kstar, r_star, r_star*4*100))

# Simulate AR(1) TFP and compute moments analytically
T <- 5000; set.seed(1)
log_A <- arima.sim(list(ar=rho_A), T, sd=sig_A)
A_sim <- exp(log_A)

# Simple linearized simulation: y_hat ≈ alpha*k_hat + A_hat
# (first-order approximation, full VFI above for exact)
# IRF to TFP shock
H <- 20
irf_A <- rho_A^(0:(H-1))   # TFP IRF (AR1)
# Output IRF ≈ (1/(1-alpha*beta))*TFP_IRF (simplified)
irf_Y <- irf_A / (1-alpha*beta*rho_A)   # rough approximation

cat("\nTheoretical RBC IRF (output response to 1% TFP shock):\n")
cat(sprintf("  Horizon 0: %.3f, 4: %.3f, 8: %.3f, 20: %.3f\n",
            irf_Y[1], irf_Y[5], irf_Y[9], irf_Y[21]))
cat("\nKey calibration targets vs model:\n")
cat(sprintf("  Capital share α = %.2f  (target: 0.36)\n", alpha))
cat(sprintf("  Annual depreciation = %.1f%%  (target ≈ 10%%)\n", delta*4*100))
cat(sprintf("  Annual interest rate = %.1f%%  (target ≈ 4%%)\n", r_star*4*100))
```

---

## 17.7 Programming Exercises

### Exercise 17.1 (APL — HP Filter)

Implement the HP filter in APL using `⌹` (matrix divide) to solve the linear system $(I + \lambda D'D)\bm\tau = \mathbf{y}$. (a) Construct the second-difference matrix $D$ using outer products and diagonal indexing. (b) Apply to simulated AR(1) GDP data with $\lambda = 1600$. (c) Verify that the cycle component has the correct spectral properties (band-pass behavior around business cycle frequencies 6–32 quarters).

### Exercise 17.2 (Python — Second Moments)

Using the VFI solution from Section 17.6: (a) simulate 10,000 periods of the RBC model; (b) compute HP-filtered second moments for $(Y, C, I, n)$; (c) compare to the data targets in Section 17.5. What are the model's main failures? (The standard RBC model generates too-correlated labor with output and too-high investment volatility relative to data.)

### Exercise 17.3 (Julia — Calibration Sensitivity)

For the RBC calibration, vary $\rho_A \in \{0.7, 0.8, 0.9, 0.95, 0.99\}$ and $\sigma_A \in \{0.005, 0.007, 0.010, 0.015\}$. For each combination: (a) simulate the model for 1,000 periods; (b) compute the output standard deviation and the consumption-to-output relative volatility; (c) identify which combination best matches U.S. data moments. Plot the results as a heat map.

### Exercise 17.4 — Log-Linearization Preview ($\star$)

Log-linearize the RBC equilibrium conditions around the non-stochastic steady state (following Chapter 27's approach): (a) the Euler equation; (b) the labor optimality condition; (c) the resource constraint; (d) the capital accumulation equation. Write the log-linearized system in the state-space form $\hat{\mathbf{y}}_{t+1} = A\hat{\mathbf{y}}_t + B\hat\varepsilon_{t+1}$. (e) Compute the eigenvalues of $A$ and verify the Blanchard–Kahn condition (one stable eigenvalue for one predetermined variable $\hat{k}_t$).

### Exercise 17.5 — Endogenous Labor ($\star\star$)

Extend the VFI implementation to include endogenous labor supply: for each $(k_i, A_j, k'_m)$ triple, solve the household's labor optimality condition $[(1-\mu)/mu] \cdot c / (1-n) = (1-\alpha)Y/n$ jointly with the resource constraint $c + k' = Y + (1-\delta)k$. (a) Implement as an inner optimization loop using scipy's `brentq`. (b) Compute the policy functions $c^*(k,A)$ and $n^*(k,A)$. (c) Show that endogenous labor amplifies the output response to TFP shocks (the labor supply channel).

---

## 17.8 Chapter Summary

**Key results:**

- The **RBC model** is the stochastic extension of the Ramsey growth model with endogenous labor, solved as a social planner's problem using VFI over the $(k, A)$ state space.
- The **state space** is $(k_t, A_t)$; the TFP process is discretized to a Markov chain via Tauchen's method, yielding the $N_A\times N_A$ transition matrix $\Pi$.
- The **Bellman equation** $V(k,A) = \max_{c,n,k'}\{u(c,1-n)+\beta\mathbb{E}[V(k',A')]\}$ is solved by VFI; convergence is guaranteed by Theorem 17.1 (same as Theorem 15.2, extended to 2D).
- **Calibration** uses national accounts (capital share $\alpha$, depreciation $\delta$), asset-return data (discount factor $\beta$), and the Solow residual (TFP parameters $\rho_A$, $\sigma_A$). The preference parameter $\mu$ is set to match a labor supply target.
- The **HP filter** $(I+\lambda D'D)^{-1}$ separates trend from cycle; in APL it is implemented as `⌹(I+lambda×(⍉D)+.×D)+.×y`.
- **Second moments** (standard deviations and correlations of HP-filtered variables) are the primary targets for model evaluation. The standard RBC model matches output variance and consumption relative volatility reasonably well, but overpredicts labor volatility.

**Connections forward:** Chapter 27 (log-linearization) provides an analytic approximation to the VFI policy functions, dramatically reducing computation time. Chapter 28 (Blanchard–Kahn) applies the determinacy framework to verify the RBC model has a unique stable solution. Chapter 37 (Great Recession replication) adds financial frictions to the RBC structure and estimates the augmented model on data.

---

*Next: Part V — Stochastic Methods for Macroeconomic Modeling*
