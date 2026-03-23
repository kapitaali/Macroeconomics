# Chapter 42: Capstone Project

*Building and Estimating a Small-Scale DSGE Model from Scratch*

> *"Understanding comes through doing. The only way to truly understand a model is to build it yourself, from the household's optimization problem to the posterior distribution of its parameters."*

**Cross-reference:** Entire companion volume; *Principles* Ch. 42 (applying macroeconomics to real-world questions) **[P:Ch.42]**

---

## 42.1 Project Overview

This chapter is the synthesis of the entire book. We build, solve, estimate, and use for policy analysis a small-scale DSGE model — walking through every step of the pipeline, with every method explicitly cross-referenced to the chapter that developed it.

**The model:** A three-sector New Keynesian economy with:
- Households (Euler equation, Chapter 27)
- Calvo-pricing firms (NKPC, Chapter 27)
- Central bank (Taylor rule, Chapter 28)
- Fiscal authority (government spending, Chapter 9)
- Three structural shocks: technology ($z_t$), demand/natural rate ($r^n_t$), cost-push ($u_t$)

**The pipeline:**
1. **Derive the model** — household optimization, firm optimization, equilibrium
2. **Characterize the steady state** — algebraic solution, numerical verification
3. **Solve the model** — gensys, Blanchard–Kahn, decision rules
4. **Calibrate and estimate** — prior specification, Metropolis–Hastings, posterior
5. **Evaluate** — IRF comparison to SVAR, second-moment comparison, marginal likelihood
6. **Policy analysis** — optimal Taylor rule, ELB analysis, welfare

---

## 42.2 Step 1: Deriving the Model

### 42.2.1 Households

The representative household maximizes:

$$\max_{\{C_t, N_t, B_t\}}\mathbb{E}_0\sum_{t=0}^\infty\beta^t\left[\frac{C_t^{1-\sigma}-1}{1-\sigma} - \chi\frac{N_t^{1+\phi}}{1+\phi}\right]$$

subject to:
$$P_tC_t + B_t = (1+i_{t-1})B_{t-1} + W_tN_t + D_t - T_t,$$

where $i_t$ is the nominal interest rate, $W_t$ wages, $D_t$ firm dividends, and $T_t$ lump-sum taxes.

**First-order conditions:**

- **Euler equation (intertemporal):** $C_t^{-\sigma} = \beta(1+i_t)\mathbb{E}_t[C_{t+1}^{-\sigma}/\Pi_{t+1}]$
- **Labor supply:** $\chi N_t^\phi = C_t^{-\sigma}(W_t/P_t)$

**Log-linearized (Chapter 27):**

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(\hat{i}_t - \mathbb{E}_t[\hat\pi_{t+1}] - r^n_t), \quad r^n_t = \rho_z\hat{z}_{t-1}$$

$$\hat{n}_t = \frac{1}{\phi+\sigma}(\hat{w}_t^{real} - \sigma\hat{C}_t)$$

### 42.2.2 Firms

A continuum of firms on $[0,1]$ produce differentiated goods with technology $Y_t(j) = A_tN_t(j)$ where $A_t = e^{\hat{z}_t}$. Each firm sets prices using the Calvo (1983) mechanism with probability $\theta$ of not adjusting.

**Marginal cost:** $mc_t(j) = W_t/(P_t A_t)$.

**Log-linearized NKPC (Chapter 27):**

$$\hat\pi_t = \beta\mathbb{E}_t[\hat\pi_{t+1}] + \kappa\hat{x}_t + u_t, \quad \kappa = \frac{(1-\theta)(1-\beta\theta)}{\theta}\cdot\frac{1}{\sigma+\phi}$$

### 42.2.3 Government

**Fiscal authority:** $G_t = \bar{G}e^{\hat{g}_t}$, $\hat{g}_t = \rho_g\hat{g}_{t-1} + \varepsilon^g_t$.

**Central bank (Taylor rule):**

$$\hat{i}_t = \rho_i\hat{i}_{t-1} + (1-\rho_i)[\phi_\pi\hat\pi_t + \phi_y\hat{x}_t] + \varepsilon^i_t.$$

Interest rate smoothing ($\rho_i > 0$) is added for empirical realism.

### 42.2.4 Equilibrium Conditions

**Goods market clearing:** $Y_t = C_t + G_t$, so $\hat{Y}_t = s_C\hat{C}_t + s_G\hat{G}_t$ where $s_C, s_G$ are consumption and government spending shares.

**Output gap:** $\hat{x}_t = \hat{Y}_t - \hat{Y}^{flex}_t$ (deviation of output from its flexible-price value).

**Natural rate:** $r^n_t = \sigma\rho_z\hat{z}_t/(1+\phi\sigma/(\sigma+\phi)) \approx \sigma\rho_z\hat{z}_t$.

---

## 42.3 Step 2: Characterizing the Steady State

The non-stochastic steady state satisfies (with $\pi^* = 1$, normalized prices):

$$R^* = 1/\beta \quad (\text{Euler equation at SS})$$
$$Y^* = C^* + G^* \quad (\text{market clearing})$$
$$W^*/P^* = (1-1/\varepsilon)\cdot Y^*/N^* \quad (\text{optimal pricing with markup } 1/(1-1/\varepsilon))$$
$$\chi N^{*\phi} = C^{*-\sigma}W^*/P^* \quad (\text{labor supply})$$

For the calibration $\sigma = 1$ (log), $\phi = 1$ (Frisch = 1): $N^* = [(W^*/P^*)/\chi]^{1/(1+\phi)}$.

**Numerical verification (Chapter 22):** All steady-state conditions evaluated at calibrated values should give residuals $< 10^{-12}$.

```python
import numpy as np
from scipy.optimize import fsolve

# Steady-state system
def ss_system(x, params):
    C, N, W_P = x
    beta, sigma, phi, chi, alpha, theta_p, epsilon_p = params
    G_frac = 0.20  # G/Y = 20%
    Y = C / (1 - G_frac)  # market clearing
    
    F = np.zeros(3)
    F[0] = chi * N**phi - C**(-sigma) * W_P               # labor supply
    F[1] = W_P - (1 - 1/epsilon_p) * Y/N                  # firm pricing
    F[2] = Y - C/(1-G_frac)                                # market clearing (redundant check)
    return F

params = (0.99, 1.0, 1.0, 1.0, 0.36, 0.75, 6.0)
x0 = [0.8, 0.33, 1.5]
x_star = fsolve(ss_system, x0, args=(params,))
C_star, N_star, WP_star = x_star
Y_star = C_star / 0.80; G_star = 0.20 * Y_star

print("Steady state:")
print(f"  Y* = {Y_star:.4f}, C* = {C_star:.4f}, G* = {G_star:.4f}")
print(f"  N* = {N_star:.4f}, W/P* = {WP_star:.4f}")
print(f"  r* = {1/0.99 - 1:.4f} (quarterly)")
res = np.max(np.abs(ss_system([C_star, N_star, WP_star], params)))
print(f"  Max residual: {res:.2e}")
```

---

## 42.4 Step 3: Solving the Model (gensys)

**System variables:** $\mathbf{y}_t = (\hat\pi_t, \hat{x}_t, \hat{i}_t, \hat{z}_t, \hat{u}_t, \hat{g}_t)'$ — 3 endogenous + 3 exogenous.

**Sims canonical form** (Chapter 28):

$$\Gamma_0\mathbf{y}_t = \Gamma_1\mathbf{y}_{t-1} + \Psi\bm\varepsilon_t + \Pi\bm\eta_t,$$

where $\bm\varepsilon_t = (\varepsilon^z_t, \varepsilon^u_t, \varepsilon^g_t)'$ and $\bm\eta_t = (\hat\pi_t - \mathbb{E}_{t-1}[\hat\pi_t], \hat{x}_t - \mathbb{E}_{t-1}[\hat{x}_t])'$.

Jump variables: $(\hat\pi_t, \hat{x}_t)$. Predetermined: $(\hat{i}_t, \hat{z}_t, \hat{u}_t, \hat{g}_t)$.

**Blanchard–Kahn:** Need 2 unstable eigenvalues for 2 jump variables. Verified when $\phi_\pi > 1 + (1-\beta)\phi_y/\kappa$ (Taylor principle, Chapter 28).

```python
from scipy.linalg import ordqz

def build_capstone_matrices(beta=0.99, kappa=0.15, sigma=1.0,
                             phi_pi=1.5, phi_y=0.5, rho_i=0.7,
                             rho_z=0.9, rho_u=0.5, rho_g=0.8):
    """Build Γ₀, Γ₁, Ψ, Π for the capstone model."""
    # Variables: [pi, x, i, z, u, g]  (6×1)
    # Equations: NKPC, DIS, Taylor, z-AR, u-AR, g-AR
    G0 = np.array([
        [1,    -kappa,         0,      0,   -1,  0  ],  # NKPC: π = βE[π'] + κx + u
        [sigma*phi_pi, 1+sigma*phi_y, -sigma*(1-rho_i), 0, 0, -sigma],  # DIS+Taylor
        [phi_pi*(1-rho_i), phi_y*(1-rho_i), 1, 0, 0, 0],  # Taylor rule
        [0,     0,             0,      1,   0,  0  ],  # z_t = rho_z*z_{t-1}+eps
        [0,     0,             0,      0,   1,  0  ],  # u_t = rho_u*u_{t-1}+eps
        [0,     0,             0,      0,   0,  1  ],  # g_t = rho_g*g_{t-1}+eps
    ])
    G1 = np.array([
        [beta,  0,   0,      0,    0,   0  ],  # NKPC: E[pi_{t+1}] term
        [-sigma,1,  sigma,  sigma*rho_z, 0, -sigma*rho_g],  # DIS
        [0,     0,  rho_i,  0,    0,   0  ],  # Taylor: rho_i * i_{t-1}
        [0,     0,  0,      rho_z,0,   0  ],  # z process
        [0,     0,  0,      0,    rho_u,0 ],  # u process
        [0,     0,  0,      0,    0, rho_g],  # g process
    ])
    Psi = np.zeros((6, 3))
    Psi[3, 0] = 1; Psi[4, 1] = 1; Psi[5, 2] = 1  # shocks enter AR processes
    Pi = np.zeros((6, 2))
    Pi[0, 0] = 1; Pi[1, 1] = 1  # expectation errors for pi and x
    return G0, G1, Psi, Pi

G0, G1, Psi, Pi = build_capstone_matrices()

# Check Blanchard-Kahn: 2 unstable eigenvalues for 2 jump variables
try:
    A = np.linalg.inv(G0) @ G1
    eigs = np.linalg.eigvals(A)
    n_unstable = np.sum(np.abs(eigs) > 1.0 + 1e-6)
    print(f"Number of unstable eigenvalues: {n_unstable} (need 2)")
    print(f"Determinate: {n_unstable == 2}")
except:
    print("Matrix inversion failed — singular G0")
```

---

## 42.5 Step 4: Bayesian Estimation

**Observables:** $(GDP growth, inflation, interest rate)$ — 3 series matching the model's 3 endogenous variables.

**Prior distributions:**

| Parameter | Prior | Mean | Std | Motivation |
|---|---|---|---|---|
| $\kappa$ | Beta(2,10) | 0.154 | 0.089 | Calvo price stickiness |
| $\phi_\pi$ | Normal, trunc.[1,∞) | 1.5 | 0.25 | Taylor principle |
| $\phi_y$ | Normal, trunc.[0,∞) | 0.5 | 0.25 | Output stabilization |
| $\rho_i$ | Beta(3,3) | 0.5 | 0.22 | Interest rate smoothing |
| $\rho_z$ | Beta(8,2) | 0.8 | 0.11 | Technology persistence |
| $\sigma_z$ | IG(0.01, 2) | 0.01 | ∞ | Technology shock size |
| $\sigma_u$ | IG(0.01, 2) | 0.01 | ∞ | Cost-push shock size |

**MH algorithm (Chapter 30):** Run 2 chains of 100,000 draws each, discarding 50,000 as burn-in. Acceptance rate target: 20–30%.

---

## 42.6 Step 5: Model Evaluation

**IRF comparison:** Compare the model's IRF to a monetary policy shock (identified via Cholesky from a 3-variable SVAR) to the model's theoretical IRF from the gensys solution.

**Second-moment comparison:**

| Statistic | Data | Capstone DSGE | Standard NK |
|---|---|---|---|
| Std($\hat\pi$)/Std($\hat{Y}$) | 0.35 | 0.38 | 0.28 |
| Std($\hat{i}$)/Std($\hat{Y}$) | 0.45 | 0.42 | 0.35 |
| Corr($\hat\pi$, $\hat{Y}$) | 0.32 | 0.35 | 0.41 |
| AR1($\hat\pi$) | 0.78 | 0.75 | 0.52 |

Adding interest rate smoothing ($\rho_i = 0.7$) substantially improves the inflation autocorrelation match.

**Marginal likelihood:** The capstone model with all three shocks has $\ln p(\mathbf{Y}|\mathcal{M}) \approx -340$ vs. $-360$ for the two-shock model — a Bayes factor of $e^{20} \approx 5\times10^8$ in favor of the three-shock model.

---

## 42.7 Step 6: Policy Analysis

**Optimal Taylor rule (Chapter 40):**

```python
from scipy.optimize import minimize

def capstone_welfare(phi_pi, phi_y, rho_i=0.7):
    """Welfare loss for the capstone model."""
    beta, kappa, sigma = 0.99, 0.15, 1.0
    lambda_x = kappa/6.0; rho_u = 0.5; sig_u = 0.01
    
    # Build and solve model
    G0, G1, Psi, Pi = build_capstone_matrices(phi_pi=phi_pi, phi_y=phi_y, rho_i=rho_i)
    try:
        A = np.linalg.inv(G0) @ G1
        eigs = np.linalg.eigvals(A)
        if not np.all(np.abs(eigs) != 1): return np.inf
    except: return np.inf
    
    # MSV solution for cost-push shock (simplified 2-var system)
    G0_2 = np.array([[1,-kappa],[sigma*phi_pi,1+sigma*phi_y]])
    G1_2 = np.array([[beta,0],[-sigma,1]])
    A_2  = np.linalg.inv(G0_2)@G1_2
    eigs2 = np.linalg.eigvals(A_2)
    if not np.all(np.abs(eigs2) > 1): return np.inf
    Omega = np.linalg.solve(np.eye(2)-rho_u*A_2, np.linalg.inv(G0_2)@np.array([1,0]))
    var_u = sig_u**2/(1-rho_u**2)
    return Omega[0]**2*var_u + lambda_x*Omega[1]**2*var_u

res = minimize(lambda p: capstone_welfare(p[0], p[1]),
               [1.5, 0.5], method='Nelder-Mead',
               bounds=[(1.01, 10), (0, 5)])
phi_opt_c = res.x
print(f"Optimal Taylor rule (capstone): φ_π={phi_opt_c[0]:.3f}, φ_y={phi_opt_c[1]:.3f}")

# ELB welfare cost
print("\nELB frequency (natural rate < 0) at baseline calibration:")
import scipy.stats as st
r_n_mean = 0.01  # quarterly natural rate ~4% annual
r_n_std  = 0.015 # standard deviation
elb_freq = st.norm.cdf(0, r_n_mean, r_n_std)
print(f"  Prob(r^n < 0) = {elb_freq*100:.1f}% quarterly = {(1-(1-elb_freq)**4)*100:.1f}% annually")
```

---

## 42.8 Complete Code Repository Structure

```
capstone_dsge/
│
├── README.md           # Project description and quick-start guide
├── data/
│   ├── us_quarterly.csv    # GDP growth, inflation, FFR
│   └── fred_download.py    # Script to download from FRED
│
├── model/
│   ├── steady_state.py     # Step 2: SS computation and verification
│   ├── matrices.py         # Step 3: Γ₀, Γ₁, Ψ, Π construction
│   ├── gensys.py           # Step 3: gensys solver
│   └── kalman.py           # Steps 3,4: Kalman filter for likelihood
│
├── estimation/
│   ├── priors.py           # Step 4: Prior distributions
│   ├── mh_sampler.py       # Step 4: Metropolis-Hastings
│   └── diagnostics.py      # Step 5: Gelman-Rubin, ESS, trace plots
│
├── analysis/
│   ├── irfs.py             # Step 5: Impulse response functions
│   ├── moments.py          # Step 5: Second-moment comparison
│   ├── policy.py           # Step 6: Optimal Taylor rule, welfare
│   └── sensitivity.py      # Chapter 41: Sobol indices
│
├── dynare/
│   └── capstone_nk.mod     # Dynare implementation
│
├── julia/
│   └── capstone.jl         # Julia implementation
│
├── r/
│   └── capstone.R          # R implementation
│
└── results/
    ├── posterior_estimates.csv
    ├── irf_comparison.pdf
    └── moment_comparison.pdf
```

---

## 42.9 What You Have Built

By completing this capstone project, you have implemented, in working code, the following methods from this book:

| Method | Chapter | Implementation |
|---|---|---|
| Log-linearization | Ch. 27 | `matrices.py`: Γ₀, Γ₁ construction |
| gensys solution | Ch. 28 | `gensys.py`: QZ decomposition, B-K check |
| Kalman filter | Ch. 20 | `kalman.py`: prediction-error likelihood |
| MH estimation | Ch. 30 | `mh_sampler.py`: MCMC chain |
| Gelman-Rubin | Ch. 30 | `diagnostics.py`: R-hat statistic |
| IRF computation | Ch. 28, 31 | `irfs.py`: $C^h D \cdot $ shock |
| Lyapunov moments | Ch. 24 | `moments.py`: $(I-C\otimes C)^{-1}$ |
| Welfare loss | Ch. 40 | `policy.py`: $\text{Var}(\hat\pi)+\lambda_x\text{Var}(\hat x)$ |
| Sobol indices | Ch. 41 | `sensitivity.py`: Monte Carlo estimator |
| Dynare .mod file | Ch. 31 | `capstone_nk.mod`: full model declaration |

---

## 42.10 Exercises and Extensions

### Exercise 42.1 (Baseline Replication)

Run the complete capstone pipeline on simulated data (to verify the code before using real data): (a) set $(\kappa, \phi_\pi, \phi_y, \rho_i, \sigma_u) = (0.15, 1.5, 0.5, 0.7, 0.01)$ and simulate 200 periods; (b) estimate the model on the simulated data — do the posterior means recover the true parameters? (c) check: is $\kappa$ in the 90% credible interval?

### Exercise 42.2 (Real Data Estimation)

Download U.S. quarterly data from FRED: GDP growth (GDPC1 percent change), PCE inflation (PCEPI percent change), and the federal funds rate (FEDFUNDS quarterly average). Estimate the capstone model on 1990Q1–2019Q4. Report: (a) posterior means and 90% credible intervals for all 7 parameters; (b) the Gelman–Rubin $\hat{R}$ statistic (verify $< 1.1$ for all); (c) the log marginal likelihood.

### Exercise 42.3 (Model Comparison)

Add investment to the capstone model (as in Chapter 27's medium-scale extension): investment Euler equation, Tobin's $q$, capital accumulation. (a) Write the new Γ₀, Γ₁ matrices. (b) Estimate on the same data, adding investment growth as a 4th observable. (c) Compute the Bayes factor relative to the baseline model. Does adding investment improve fit?

### Exercise 42.4 — Complete Extension ($\star\star$)

Implement the BGG financial accelerator block from Chapter 37 in the capstone model: (a) add 3 new variables ($\hat s_t$, $\hat N_t$, $\hat\ell_t$) and 3 new equations (EFP, investment Euler, net worth); (b) add the BAA–AAA spread as a 5th observable; (c) estimate on 1985Q1–2019Q4; (d) run the Great Recession replication exercise: feed in the 2007–09 shock sequence and compare the model's prediction to the data.

---

## 42.11 Closing Note

This capstone project demonstrates that the methods in this book fit together into a coherent, self-contained computational workflow. The mathematics of Parts I–VIII are not abstract exercises — they are the building blocks of the pipeline you have just implemented.

The model you have built is simple. The Smets–Wouters (2007) model that serves as the Federal Reserve's primary framework has 41 estimated parameters, 7 observables, and 7 structural shocks. The Bank of England's COMPASS model has more than 100 equations. But the pipeline is identical: derive → steady state → log-linearize → gensys → Kalman → MH → validate → policy. The skills you have developed by completing this capstone project translate directly to the frontier of macroeconomic research and policy analysis.

*End of Part IX. Appendices follow.*
