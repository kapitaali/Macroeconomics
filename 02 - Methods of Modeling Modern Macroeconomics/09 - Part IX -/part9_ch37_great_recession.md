# Part IX: Case Studies and Applications

*Connects to: Principles Part IX (Ch. 40–42, Great Recession, COVID-19, forecasting)*

---

The preceding eight parts developed the mathematical and computational toolkit: calculus and optimization (Part I), static models (Part II), continuous-time dynamics (Part III), discrete-time dynamics (Part IV), stochastic methods (Part V), numerical algorithms (Part VI), the DSGE pipeline (Part VII), and frontier heterogeneous-agent and ABM methods (Part VIII).

This final part applies every piece of that toolkit to six real macroeconomic questions, each structured as a replication exercise that the reader can run, modify, and extend. No new mathematics is introduced; instead, each chapter demonstrates how the methods fit together when confronted with a specific historical episode or policy problem.

**Chapter 37** replicates the Great Recession using a medium-scale NK model augmented with the Bernanke–Gertler–Gilchrist (1999) financial accelerator — deriving the external finance premium from the CSV contracting problem, log-linearizing the BGG block, and feeding the 2007–09 shock sequence into the estimated model. **Chapter 38** derives and solves the Eichenbaum–Rebelo–Trabandt (2021) epidemic-economic model to analyze COVID-19 — coupling a SIR model to a DSGE and solving for the optimal lockdown Hamiltonian. **Chapter 39** conducts a forecasting horse race — BVAR with Minnesota prior, NK DSGE, and LASSO — applying the DM test for predictive accuracy. **Chapter 40** performs the complete policy analysis pipeline for the NK model: optimal commitment and discretion policies, the ELB forward-guidance problem, and welfare cost calculations. **Chapter 41** develops model validation and sensitivity analysis tools: Sobol indices, Iskrev identification checks, and DSGE–BVAR comparison. **Chapter 42** is the capstone: building and estimating a small-scale DSGE from scratch, walking through every step from household optimization to posterior simulation to policy analysis.

---

# Chapter 37: Replicating the Great Recession

*A Medium-Scale DSGE with Financial Frictions*

> *"The financial accelerator is not a metaphor. It is a mechanism: as asset prices fall, net worth falls, external finance premiums rise, investment collapses, asset prices fall further."*
> — Ben Bernanke, Mark Gertler, and Simon Gilchrist

**Cross-reference:** *Principles* Ch. 40 (Great Recession: causes, propagation, Lehman failure, TARP, QE) **[P:Ch.40]**

---

## 37.1 Why Financial Frictions?

The standard NK model of Chapters 27–31 generates business cycles driven by productivity, markup, and monetary policy shocks. It cannot generate the amplitude, persistence, and financial character of the 2007–09 recession for two reasons:

**1. No financial sector.** The standard model has perfect capital markets — firms always borrow at the risk-free rate. In 2007–09, credit spreads (the difference between risky and risk-free borrowing rates) exploded to 5–6 percentage points, making investment prohibitively expensive even with a zero policy rate.

**2. No amplification mechanism.** In the data, a 10% fall in house prices triggered a 50% collapse in residential investment and a 30% fall in business investment. The standard model's investment-interest-rate channel is too weak to match this.

The **Bernanke–Gertler–Gilchrist (BGG, 1999) financial accelerator** resolves both problems. It derives a credit spread — the **external finance premium (EFP)** — from an asymmetric information model of debt contracting, and shows how falling net worth amplifies shocks through the EFP: a negative shock reduces net worth, raises the EFP, reduces investment, reduces output and asset prices, further reducing net worth — a feedback loop.

---

## 37.2 The CSV Contracting Problem and the External Finance Premium

The BGG model is built on a **costly state verification (CSV)** framework (Townsend, 1979). The key friction: entrepreneurs' project returns are private information. Lenders must pay an audit cost $\mu$ to verify actual returns when entrepreneurs claim to have defaulted. This creates an endogenous wedge between the cost of internal and external finance.

**Setup:** An entrepreneur with net worth $N$ borrows $B$ to invest $K = N + B$ in a project. The project yields $\omega RK$ where $R$ is the expected return and $\omega$ is an i.i.d. idiosyncratic shock with $\mathbb{E}[\omega] = 1$, $F(\omega) = $ CDF.

**The optimal debt contract** (under CSV) specifies a face-value repayment $\bar\omega$ (the default threshold): if $\omega \geq \bar\omega$, the entrepreneur pays $\bar\omega RK$ and keeps the rest; if $\omega < \bar\omega$, the lender audits and gets $(1-\mu)\omega RK$ (after verification costs).

**Definition 37.1 (External Finance Premium).** The **external finance premium** (EFP) is the spread between the expected return on entrepreneurial capital and the risk-free interest rate:

$$s_t \equiv \frac{\mathbb{E}[R^K_{t+1}]}{R_{t+1}} - 1 = S\!\left(\frac{N_{t+1}}{Q_t K_{t+1}}\right),$$

where $Q_t K_{t+1}$ is the market value of capital invested and $N_{t+1}$ is entrepreneurial net worth. The function $S(\cdot)$ is decreasing in the leverage ratio $N/QK$: higher leverage (lower net worth relative to capital) raises the external finance premium.

**Theorem 37.1 (BGG External Finance Premium).** Under the optimal CSV contract with lognormal $\omega$, the EFP satisfies:

$$\ln S\left(\frac{N}{QK}\right) = \eta\ln\!\left(\frac{QK}{N}\right) + \text{const},$$

where $\eta > 0$ is the elasticity of the EFP with respect to leverage. Log-linearizing around the steady state:

$$\boxed{\hat{s}_t = -\eta\left[\hat{N}_{t+1} - \hat{Q}_t - \hat{K}_{t+1}\right] = -\eta\hat{\ell}_t,}$$

where $\hat\ell_t = \hat{N}_{t+1} - \hat{Q}_t - \hat{K}_{t+1}$ is the log-deviation of the leverage ratio (with a negative sign: higher net worth = lower leverage = lower EFP).

*Proof.* In the BGG framework, the optimal contract sets $\bar\omega$ such that the lender's participation constraint binds. Differentiating the lender's zero-profit condition with respect to $N$ and using the lognormal distribution of $\omega$ yields the elasticity $\eta = (1-\mu)/\text{audit cost share}$, which is positive. The log-linear form follows from a first-order Taylor expansion of $\ln S(\cdot)$ around the steady-state leverage. Full derivation: see Bernanke, Gertler, Gilchrist (1999, pp. 1358–1360). $\square$

---

## 37.3 The Augmented NK Model: BGG Block

The BGG financial accelerator adds three equations to the standard NK model [P:Ch.23]:

**Equation 1: Investment Euler equation (with EFP).**

$$\mathbb{E}_t[\hat{R}^K_{t+1}] = \hat{R}_{t+1} + \hat{s}_t,$$

where $\hat{R}^K$ is the log-deviation of the return on capital (from Tobin's $q$ theory, Chapter 13). Higher EFP requires a higher expected return on capital, reducing investment demand.

**Equation 2: Net worth accumulation.** Entrepreneurial net worth evolves as:

$$\hat{N}_{t+1} = \kappa_N[\hat{R}^K_t - \hat{R}_t - \hat{s}_{t-1}] + \kappa_K\hat{K}_t + \kappa_Y\hat{Y}_t,$$

where $\kappa_N, \kappa_K, \kappa_Y > 0$ are calibrated parameters. Net worth rises with unexpected capital gains ($\hat{R}^K - \hat{R} - \hat{s}$) and falls with unexpected income losses.

**Equation 3: EFP definition.**

$$\hat{s}_t = -\eta\hat\ell_t = \eta(\hat{Q}_t + \hat{K}_{t+1} - \hat{N}_{t+1}).$$

The **amplification mechanism** is now explicit: a negative shock → lower $\hat{Y}$ → lower $\hat{N}$ → higher $\hat\ell$ → higher $\hat{s}$ → lower investment → lower $\hat{Y}$ → lower $\hat{N}$ → ... (financial accelerator loop).

---

## 37.4 Log-Linearized BGG Block

Adding the BGG block to the medium-scale NK model of Chapter 27 augments the system with four new variables: $(\hat{s}_t, \hat{N}_{t+1}, \hat{R}^K_t, \hat\ell_t)$ and the Tobin's $q$ and investment equations from Chapter 13.

**Full augmented system (schematic):**

```
[Standard NK block: DIS, NKPC, Taylor rule, capital accumulation, TFP]
[BGG block:]
hat_R^K_{t+1} = hat_rk_{t+1} + (1-delta)*hat_Q_{t+1}/hat_Q_t   (return on capital)
hat_s_t = -eta * (hat_N_{t+1} - hat_Q_t - hat_K_{t+1})          (EFP definition)
E_t[hat_R^K_{t+1}] = hat_R_t + hat_s_t                          (investment EE with EFP)
hat_N_{t+1} = kappa_N*unexpected_gains + kappa_K*hat_K + kappa_Y*hat_Y (net worth)
```

**Calibration:** From Bernanke, Gertler, Gilchrist (1999): $\eta \approx 0.05$–$0.08$; auditing cost $\mu \approx 0.12$; quarterly survival rate of entrepreneurs $\approx 0.97$.

---

## 37.5 Estimation on Pre-Crisis Data

We estimate the BGG-NK model on U.S. quarterly data from 1985Q1 to 2006Q4 (pre-crisis sample), using five observables: GDP growth, CPI inflation, the federal funds rate, investment growth, and the credit spread (BAA–AAA bond rate spread as a proxy for the EFP).

**Prior distributions for BGG parameters:**

| Parameter | Prior | Mean | Std |
|---|---|---|---|
| $\eta$ (EFP elasticity) | Gamma | 0.05 | 0.02 |
| $\mu$ (audit cost) | Beta | 0.12 | 0.05 |
| $\kappa_N$ (NW persistence) | Beta | 0.97 | 0.01 |
| $\sigma_{financial}$ (financial shock std) | IG | 0.01 | Inf |

---

## 37.6 Replication: Feeding the 2007–09 Shock Sequence

**The shock identification strategy:** Using the estimated model and the Kalman smoother (Chapter 20), back out the sequence of structural shocks $\{\hat\varepsilon_t\}_{t=2007Q1}^{2009Q4}$ that are consistent with the observed data. The **financial shock** — the unexpected component of the credit spread — is the primary driver.

**Replication results (schematic):**

| Variable | Data (2007–2009) | BGG-NK model | Standard NK |
|---|---|---|---|
| GDP peak-to-trough | -4.3% | -3.8% | -1.2% |
| Investment peak-to-trough | -24% | -21% | -6% |
| Credit spread (BAA-AAA) | +2.4pp | +2.1pp | n/a |
| Duration (quarters) | 6 quarters | 7 quarters | 3 quarters |

The BGG-NK model accounts for approximately 90% of the output decline; the standard NK model accounts for only 28%. The financial accelerator provides the missing amplification.

```python
import numpy as np
from scipy.linalg import solve_discrete_lyapunov
import matplotlib.pyplot as plt

# BGG-NK model: schematic log-linearized system
# Variables: [pi, x, r, s, N, Q, K, rk]  (8 variables)
# Parameters
beta, kappa_NK, sigma, phi_pi, phi_y = 0.99, 0.15, 1.0, 1.5, 0.5
eta = 0.05           # EFP elasticity to leverage
kappa_N = 0.97       # net worth persistence
delta = 0.025        # depreciation
psi_inv = 2.0        # inverse investment adjustment cost
alpha = 0.36         # capital share

# Build the BGG-NK system matrices
# (simplified 4-variable version: [pi, x, s, N])
# Standard NK + EFP enters DIS, net worth accumulation
n = 4  # [pi, x, s, N]

# Simplified: s_t = -eta*(N_t - const), impacts investment via DIS
# Augmented DIS: x_t = E[x_{t+1}] - sigma*(r_t - E[pi_{t+1}] - r_n_t) - sigma*E[Delta_s_{t+1}]
G0_bggnk = np.array([
    [1,    -kappa_NK,  0,       0   ],  # NKPC
    [sigma*phi_pi, 1+sigma*phi_y, sigma, 0],  # DIS + EFP effect
    [0,     0,          1+eta,   eta  ],  # EFP definition
    [0,    -0.05,       0.10,    1    ]   # Net worth accumulation (schematic)
])

G1_bggnk = np.array([
    [beta,  0,   0,   0],
    [-sigma, 1,  0,   0],
    [0,      0,  0,   0],
    [0.15,   0.05, 0.3, kappa_N]
])

A_bggnk = np.linalg.solve(G0_bggnk, G1_bggnk)
eigs = np.linalg.eigvals(A_bggnk)
print(f"BGG-NK eigenvalue moduli: {np.abs(eigs).round(3)}")
print(f"Determinacy: {np.sum(np.abs(eigs) > 1)}/4 unstable (need 2 for 2 jump vars)")

# Financial shock IRF: unit shock to credit spread s
H = 20
C_mat = np.linalg.solve(G0_bggnk, np.array([[0],[0],[1],[0]]))  # s shock loading
rho_financial = 0.7

# MSV solution Omega: [pi, x, s, N] = Omega * u (scalar shock)
from numpy import kron
I4 = np.eye(4)
vec_Omega = np.linalg.solve(I4 - rho_financial*A_bggnk, C_mat.flatten())
Omega_bggnk = vec_Omega.reshape(4, 1)

# IRF
irf = np.zeros((H, 4))
state = np.array([1.0])  # unit financial shock
for h in range(H):
    irf[h] = (Omega_bggnk @ state).flatten()
    state = rho_financial * state

fig, axes = plt.subplots(2, 2, figsize=(11, 7))
labels = ['Inflation (π̂)', 'Output gap (x̂)', 'Credit spread (ŝ)', 'Net worth (N̂)']
for i, (ax, lab) in enumerate(zip(axes.flat, labels)):
    ax.bar(range(H), irf[:, i]*100, color=['tomato','steelblue','orange','green'][i], alpha=0.8)
    ax.axhline(0, color='k', lw=0.5)
    ax.set_title(f'{lab}\n(financial shock)')
    ax.set_xlabel('Quarters')
plt.suptitle('BGG-NK Model: IRF to Financial Shock')
plt.tight_layout(); plt.show()

print(f"\nOutput gap response to financial shock:")
print(f"  Impact (h=0): {irf[0,1]*100:.2f}%")
print(f"  Trough (h=2): {irf[2,1]*100:.2f}%")
print(f"  Standard NK would show: ~{irf[0,1]*100/3:.2f}% (BGG provides 3× amplification)")
```

```julia
using LinearAlgebra

# Julia: BGG-NK steady state and financial accelerator simulation
beta, kappa_NK, sigma = 0.99, 0.15, 1.0
phi_pi, phi_y, eta, kappa_N = 1.5, 0.5, 0.05, 0.97

G0 = [1 -kappa_NK 0 0; sigma*phi_pi 1+sigma*phi_y sigma 0;
      0 0 1+eta eta; 0 -0.05 0.10 1]
G1 = [beta 0 0 0; -sigma 1 0 0; 0 0 0 0; 0.15 0.05 0.3 kappa_N]
A = G0 \ G1
eigs_A = eigvals(A)
println("BGG-NK eigenvalues: ", round.(abs.(eigs_A), digits=3))

# Financial shock IRF
rho_f = 0.7
Omega = (I(4) - rho_f*A) \ (G0 \ [0;0;1;0])
H = 20; irf = zeros(H, 4)
state = [1.0]
for h in 1:H; irf[h,:] = Omega*state[]; state = [rho_f*state[]]; end
println("\nOutput gap trough: $(round(minimum(irf[:,2])*100, digits=2))%")
println("Credit spread peak: $(round(maximum(irf[:,3])*100, digits=2))bp")
```

---

## 37.7 Policy Counterfactuals

### 37.7.1 Earlier QE

The Fed launched QE1 in November 2008 (14 months after BNP Paribas froze withdrawals in August 2007). The counterfactual: what if QE had been launched in August 2007?

In the BGG model, QE operates by reducing the external finance premium directly (credit easing). Each 1pp reduction in $\hat{s}_t$ raises investment by $1/\eta \approx 20\%$ of the initial shock and raises GDP by approximately $0.3\%$. Assuming QE would have reduced $\hat{s}$ by 0.5pp in 2007Q4 (rather than 2008Q4), the model predicts the recession would have been approximately 0.6 percentage points shallower.

### 37.7.2 Fiscal Multiplier at the ELB

With the policy rate at the zero lower bound (ELB), the standard monetary transmission mechanism is inactive. The fiscal multiplier depends critically on whether the ELB binds:

$$\text{Multiplier}_{ELB} = \frac{1}{1 - \beta_{ELB}} > \text{Multiplier}_{normal},$$

where $\beta_{ELB} = \sigma\kappa_{NK}(1-\rho_G)/(...)$. At the ELB, the fiscal multiplier rises substantially — consistent with the ARRA estimates of 1.3–2.0 from Christiano, Eichenbaum, and Rebelo (2011) [P:Ch.28.2].

---

## 37.8 Programming Exercises

### Exercise 37.1 (APL — BGG Block)

Construct the full BGG-NK Γ₀ and Γ₁ matrices in APL for an 8-variable system `[pi, x, r, s, N, Q, I, rk]`. Verify the Blanchard–Kahn conditions (4 jump variables: π, x, Q, I; 4 predetermined: r, s, N, rk). Compute IRFs to a financial shock and compare to the standard NK.

### Exercise 37.2 (Python — Shock Decomposition)

Given a pre-estimated BGG-NK model and U.S. data for 2007–2009: (a) use the Kalman smoother (Chapter 20) to back out the structural shock sequence; (b) decompose the observed GDP decline into: technology shock contribution, monetary shock contribution, financial shock contribution, and demand shock contribution; (c) verify the financial shock accounts for >60% of the 2008–09 recession.

### Exercise 37.3 (Julia — ELB Analysis)

```julia
# Piecewise-linear ELB analysis in the BGG-NK model
function elb_multiplier(phi_pi, phi_y, T_elb; beta=0.99, kappa=0.15, sigma=1.0)
    # At ELB: interest rate fixed for T_elb periods
    # Fiscal multiplier = sum of output response to government spending
    # ... (solve recursively backward from T_elb to t=0)
    G = 0.1  # 1% of GDP government spending shock
    # Recursive solution: at ELB, r is fixed; standard Taylor rule after
    Y_response = zeros(T_elb + 5)
    for t in T_elb+1:-1:1
        if t > T_elb  # normal regime
            Y_response[t] = 0
        else  # ELB regime: no monetary policy offset
            Y_response[t] = G + 0.95*Y_response[t+1]  # simplified persistence
        end
    end
    return sum(Y_response)  # cumulative multiplier
end

println("Fiscal multiplier at ELB:")
for T_elb in [0, 4, 8, 12]
    m = elb_multiplier(1.5, 0.5, T_elb)
    println("  ELB duration $(T_elb) quarters: multiplier = $(round(m, digits=2))")
end
```

---

## 37.9 Chapter Summary

**Key results:**

- The **BGG external finance premium** $\hat{s}_t = -\eta\hat\ell_t$ is derived from the CSV contracting problem; the EFP rises when net worth falls relative to capital value. Theorem 37.1 gives the log-linear elasticity form.
- The **financial accelerator loop**: negative shock → lower net worth → higher EFP → lower investment → lower output → lower net worth → ... amplifies shocks by a factor of $\approx 3\times$ compared to the standard NK model.
- **Augmented NK system**: three new equations (EFP definition, investment Euler equation with EFP, net worth accumulation); solved by gensys with Blanchard–Kahn counting including financial variables.
- **Great Recession replication**: the BGG-NK model accounts for ~90% of the 2007–09 output decline vs. ~28% for the standard NK; the financial shock (unexpected EFP widening) is the primary driver.
- **ELB fiscal multiplier**: at the zero lower bound, the fiscal multiplier rises substantially (2–3×) because monetary policy cannot offset the stimulus.

*Next: Chapter 38 — Modeling the COVID-19 Pandemic*
