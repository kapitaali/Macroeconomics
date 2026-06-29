# Chapter 29: Perturbation Methods

*Higher-Order Approximations for Nonlinear Dynamics*

> *"The first-order approximation tells you the average. The second-order approximation tells you about the risk. You need both."*

**Cross-reference:** *Principles* Ch. 39 (future methods: non-Gaussian risks, rare disasters); Ch. 12 (investment under uncertainty, real options) **[P:Ch.39, P:Ch.12]**

---

## 29.1 Why Go Beyond First Order?

Chapter 27's log-linearization is accurate to $O(\varepsilon^2)$ and is the workhorse for most policy questions. But three economically important phenomena are invisible at first order:

**1. Precautionary saving.** The Euler equation under uncertainty gives $\mathbb{E}_t[u'(c_{t+1})] = u'(c_t)/\beta(1+r)$. By Jensen's inequality, $\mathbb{E}[u'(c)] > u'(\mathbb{E}[c])$ when $u'' < 0$ (concavity) and $u''' > 0$ (prudence, positive third derivative). Households save more than certainty-equivalent income — the precautionary saving motive. This is a second-order effect ($\propto \sigma_c^2$) invisible in the first-order approximation.

**2. Risk premia.** The equity premium is $\mathbb{E}[R_{t+1}^e - R^f] = -\text{Cov}(M_{t+1}, R_{t+1}^e) / \mathbb{E}[M_{t+1}]$. The covariance between the SDF $M_{t+1}$ and equity returns is a second-order quantity — the product of two first-order deviations. Asset pricing requires at least second-order approximations.

**3. ELB dynamics and welfare comparisons.** The welfare loss from inflation volatility is $\mathcal{L} \approx \frac{\kappa}{2}\text{Var}(\hat\pi) + \frac{\lambda_x}{2}\text{Var}(\hat{x})$. Variances are second-order terms. Comparing welfare across policy regimes requires second-order accuracy.

---

## 29.2 The Perturbation Approach

**Definition 29.1 (Perturbation Parameter).** Perturbation methods introduce a scalar $\varepsilon \in [0,1]$ that scales all shocks: the model's stochastic equations replace $\varepsilon_t \to \varepsilon\,\varepsilon_t$. When $\varepsilon = 0$: the deterministic steady state. When $\varepsilon = 1$: the full stochastic model. The solution is expanded as a power series in $\varepsilon$ and the state deviations.

The $k$-th order perturbation approximates the policy functions:

$$c(k_t, A_t; \varepsilon) = c^{SS} + c_k(k_t - k^*) + c_A(A_t - A^*) + \frac{1}{2}c_{kk}(k_t-k^*)^2 + c_{kA}(k_t-k^*)(A_t-A^*) + \frac{1}{2}c_{AA}\sigma^2 + \ldots$$

The $\frac{1}{2}c_{AA}\sigma^2$ term is the **precautionary term** — nonzero only at second order and proportional to the variance of shocks.

---

## 29.3 Setting Up the Second-Order System

For the scalar RBC model with state $(k_t, A_t)$ and control $c_t$, the equilibrium conditions are:

$$F(k_{t+1}, c_{t+1}, k_t, c_t, A_{t+1}, A_t; \varepsilon) = 0,$$

where $F$ collects the Euler equation, resource constraint, and TFP process. The solution at order $\varepsilon^0$: the steady state $(k^*, c^*, A^*)$. The solution at order $\varepsilon^1$: the log-linearized policy functions from Chapter 27. The solution at order $\varepsilon^2$: the second-order correction.

**The second-order expansion.** Writing the policy function as:

$$c_t = g(k_t, A_t; \varepsilon) = g_0 + g_x\mathbf{x}_t + \frac{1}{2}\mathbf{x}_t'G_{xx}\mathbf{x}_t + \frac{1}{2}g_{\varepsilon\varepsilon}\varepsilon^2 + \varepsilon g_{x\varepsilon}\mathbf{x}_t + O(\varepsilon^3),$$

where $\mathbf{x}_t = (k_t - k^*, A_t - A^*)'$ is the state deviation vector. The new terms at second order:

- $\mathbf{x}_t'G_{xx}\mathbf{x}_t/2$: interaction and squared terms (state-dependent risk adjustment).
- $g_{\varepsilon\varepsilon}\varepsilon^2/2$: **stochastic steady state correction** — the constant shift in the level of the policy function due to the presence of risk.

**Theorem 29.1 (Breakdown of Certainty Equivalence at Second Order).** At first order, the policy function takes the form $c_t^{(1)} = c^{SS} + g_x^{(1)}\mathbf{x}_t$ — independent of the variance $\varepsilon^2\sigma^2$. At second order, the stochastic correction $g_{\varepsilon\varepsilon}\varepsilon^2/2$ appears: the policy function depends on the variance of shocks. This is the **breakdown of certainty equivalence**: agents behave differently under uncertainty than under certainty, even near the steady state.

*Proof.* At first order: $\mathbb{E}[F_c] = F_c(\mathbb{E}[c])$ (by linearity). The solution depends only on $\mathbb{E}[\mathbf{x}_{t+1}|\mathbf{x}_t] = A\mathbf{x}_t$ and does not involve $\text{Var}[\mathbf{x}_{t+1}|\mathbf{x}_t] = \varepsilon^2Q$. At second order: the Euler equation involves $\mathbb{E}_t[F_{cc}(c_{t+1}-c^*)^2]$ — the variance of future consumption appears. The coefficient $g_{\varepsilon\varepsilon}$ is determined by solving the second-order system, which depends on $\text{Var}[\varepsilon_{t+1}]$. $\square$

---

## 29.4 Computing the Second-Order Solution

The second-order coefficient matrices $G_{xx}$ and $g_{\varepsilon\varepsilon}$ are found by differentiating the equilibrium conditions twice with respect to states and the perturbation parameter.

### 29.4.1 The System to Be Solved

Taking the second derivative of $\mathbb{E}_t[F(\cdot;\varepsilon)] = 0$ with respect to $(\mathbf{x}_t, \varepsilon)$ at $(\mathbf{x}^*, 0)$:

$$\underbrace{A_{y'y'}}_{n\times n^2}\text{vec}(G_{xx}) + \underbrace{B_{y'y}}_{n\times n^2}(G_x\otimes G_x) + \underbrace{C_{yy}}_{n\times n^2}(I\otimes I) = 0,$$

where $A_{y'y'}$, $B_{y'y'}$, $C_{yy}$ are collections of second partial derivatives of $F$ evaluated at the steady state — computable analytically or by automatic differentiation (ForwardDiff.jl, JAX).

This is a **linear system** in the unknown $\text{vec}(G_{xx})$, solvable by standard methods. The stochastic correction:

$$g_{\varepsilon\varepsilon} = -(A_{y'y'} + B_{y'}\Phi)\backslash(A_{y'\varepsilon\varepsilon} + B_{y\varepsilon\varepsilon} + \text{tr}[G_{xx}\Sigma]),$$

where $\Sigma = \text{Cov}(\varepsilon_t)$.

---

## 29.5 Pruning: Avoiding Explosive Paths

A practical problem with higher-order perturbation: the second-order terms $(k_t - k^*)^2$ can grow without bound even when the linear part is stable, because squared deviations can be explosive even if deviations themselves are stable. **Pruning** (Kim, Kim, Schaumburg, and Sims, 2008) addresses this by separating the first- and second-order components.

**Definition 29.2 (Pruning).** The pruned second-order approximation maintains two separate state vectors:

$$\mathbf{x}_t^{(1)} = A\mathbf{x}_{t-1}^{(1)} + B\varepsilon_t \quad \text{(first-order component)},$$
$$\mathbf{x}_t^{(2)} = A\mathbf{x}_{t-1}^{(2)} + \frac{1}{2}\mathbf{x}_{t-1}^{(1)'} G_{xx}\mathbf{x}_{t-1}^{(1)} + \frac{1}{2}g_{\varepsilon\varepsilon}\sigma^2 \quad \text{(second-order correction)}.$$

The pruned approximation: $c_t \approx c^{SS} + g_x\mathbf{x}_t^{(1)} + \mathbf{x}_t^{(1)'}G_{cx}\mathbf{x}_t^{(1)}/2 + g_{\varepsilon\varepsilon}\sigma^2/2$.

**Why pruning works:** By separately tracking the first-order state, we prevent the second-order correction from feeding back into the first-order dynamics. The first-order component is stable (by the Blanchard–Kahn conditions); the second-order correction inherits this stability because it is driven by the square of the stable first-order component, which remains bounded.

---

## 29.6 Worked Example: Precautionary Saving in the RBC Model

*Cross-reference: Principles Ch. 11.5 (buffer-stock saving, precautionary motive)* **[P:Ch.11.5]**

**Setup:** Standard RBC calibration ($\alpha = 0.36$, $\delta = 0.025$, $\beta = 0.99$, $\sigma = 2$, $\rho_A = 0.95$, $\sigma_A = 0.0072$).

**First-order solution:** At first order, the consumption policy function is $\hat{c}_t = g_k^{(1)}\hat{k}_t + g_A^{(1)}\hat{A}_t$. The stochastic steady state coincides with the deterministic steady state: $\mathbb{E}[\hat{c}] = 0$.

**Second-order correction:** The precautionary saving term shifts the stochastic steady state below the deterministic one: $\hat{c}^{SS,2} = g_{\varepsilon\varepsilon}\sigma_A^2/2 < 0$. Households facing income uncertainty save more than under certainty — they hold a buffer stock.

**Quantitative magnitude:** For the standard calibration, $g_{\varepsilon\varepsilon}\sigma_A^2/2 \approx -0.5\sigma_A^2\cdot\sigma = -0.5\times(0.0072)^2\times 2 \approx -5\times10^{-5}$, i.e., approximately 0.005% of steady-state consumption. Small but nonzero — and critical for welfare analysis.

```python
import numpy as np
from scipy.linalg import solve_discrete_lyapunov

# RBC model: first-order solution and precautionary saving computation
alpha, delta, beta, sigma_crra = 0.36, 0.025, 0.99, 2.0
rho_A, sig_A = 0.95, 0.0072

# Steady state
r_ss = 1/beta - 1 + delta
k_ss = (alpha/r_ss)**(1/(1-alpha))
c_ss = k_ss**alpha - delta*k_ss
y_ss = k_ss**alpha

print(f"Steady state: k*={k_ss:.4f}, c*={c_ss:.4f}, y*={y_ss:.4f}")

# First-order policy functions (from log-linearization and gensys):
# From solving the linearized RBC, coefficients approx:
g_k = (1 - delta*(1 - alpha)) / (1 - beta*(1-delta))  # rough estimate
g_A = 1.0 / (1 - alpha*beta*rho_A)  # output response to TFP

print(f"\nFirst-order policy: ĉ_t ≈ {g_k:.3f}*k̂_t + {g_A:.3f}*Â_t")

# Second-order stochastic correction:
# From Euler eq: E[u'(c_{t+1})] = u'(c^SS) / [beta*(1+r)]
# Under uncertainty: E[(c^SS + hat_c)^{-sigma}] ≈ (c^SS)^{-sigma}(1 + sigma*(sigma+1)/2 * Var(hat_c)/c^SS^2)
# Stochastic steady state: hat_c^SS_2 = -sigma/2 * Var(hat_c) / c^SS^2 * (something)

# Variance of consumption: from Lyapunov equation on linearized system
# State: hat_k
# hat_k_{t+1} = (1-delta)*hat_k_t + delta/y_ss * g_k * hat_k_t + delta*g_A*hat_A_t + ... 
# Approximate: Var(hat_c) = g_A^2 * Var(hat_A)
var_A = sig_A**2 / (1 - rho_A**2)
var_c1 = g_A**2 * var_A  # rough first-order variance

# Precautionary saving term (second-order stochastic SS shift)
# From Kimball prudence: savings increase proportional to var(c)*prudence
prudence = sigma_crra * (sigma_crra + 1)  # CRRA: u'''(c) > 0
precautionary_correction = -0.5 * prudence * var_c1  # approximate
print(f"\nPrecautionary saving term: {precautionary_correction:.6f} (≈ {100*precautionary_correction:.4f}% of SS consumption)")
print(f"Var(ĉ) ≈ {var_c1:.6f}, std(ĉ) ≈ {np.sqrt(var_c1)*100:.2f}%")

# Compare 1st vs 2nd order moments
print(f"\nFirst order: E[ĉ] = 0 (no precautionary correction)")
print(f"Second order: E[ĉ] ≈ {precautionary_correction:.6f} (precautionary motive)")
```

```julia
using LinearAlgebra, ForwardDiff

println("Second-order perturbation: Dynare approach")
println("In Dynare: add 'order=2;' to stoch_simul command")
println("Comparison of 1st vs 2nd order impulse responses:")

# For illustration: compare IRF shapes qualitatively
alpha, delta, beta, rho_A, sig_A = 0.36, 0.025, 0.99, 0.95, 0.0072
sigma_crra = 2.0

H = 20
irf_1st = zeros(H); irf_2nd = zeros(H)

# First-order IRF: exponential decay
for h in 1:H
    irf_1st[h] = rho_A^(h-1) * sig_A
end

# Second-order correction: includes (shock)^2 term
# Rough approximation: adds concavity correction
for h in 1:H
    base = rho_A^(h-1)
    corr = -0.5 * sigma_crra * base^2 * sig_A^2  # second-order concavity term
    irf_2nd[h] = base * sig_A + corr
end

println("TFP IRF comparison (1st vs 2nd order, first 5 periods):")
for h in 1:5
    println("  h=$h: 1st-order=$(round(irf_1st[h],digits=6)), 2nd-order=$(round(irf_2nd[h],digits=6)), diff=$(round(irf_2nd[h]-irf_1st[h],digits=8))")
end
println("  (2nd-order correction is O(σ²) ≈ $(round(sig_A^2,digits=8)) per period)")
```

---

## 29.7 Programming Exercises

### Exercise 29.1 (APL — Variance of HP-Filtered Output)

Using the first-order solution from Chapter 28, compute the theoretical variance of HP-filtered output. (a) From the decision rule $\hat{Y}_t = C_Y\hat{k}_t + D_Y\hat{A}_t$, compute $\text{Var}(\hat{Y})$ using the Lyapunov equation: `Sigma_y = (⌹ I - C kron C) +.× vec_D_Sigma_D`. (b) Apply the HP filter transfer function (in frequency domain: multiply spectrum by HP weight function) to convert GDP variance to HP-filtered variance. (c) Verify this matches the simulated standard deviation from Chapter 26.

### Exercise 29.2 (Python — Dynare Order=2 Replication)

Run the RBC model in Dynare with `order=1` and `order=2`. (a) Compare the simulated moments (HP-filtered standard deviations) from both orders. (b) Extract the `oo_.dr.ghs2` vector (the stochastic steady state correction) and verify it is negative for consumption (precautionary saving). (c) Plot the impulse responses from both orders and identify at which horizons the second-order correction is largest.

### Exercise 29.3 (Julia — Risk Premium via 2nd Order)

The equity premium in the RBC model requires second-order accuracy. (a) Compute $\text{Cov}_t(\hat{M}_{t+1}, \hat{R}_{t+1}^K)$ using the second-order solution, where $\hat{M}_{t+1} = -\sigma\hat{c}_{t+1}$ and $\hat{R}_{t+1}^K = \alpha\hat{Y}_{t+1} - \hat{K}_t + (1-\delta)$. (b) The equity premium is $EP = -\text{Cov}(\hat{M}_{t+1},\hat{R}_{t+1}^K)\times(1+R^*)$. (c) For $\sigma = 2$ and the standard calibration, what is the annual equity premium? How does it compare to the data ($\approx 6\%$)?

### Exercise 29.4 — Pruning Stability ($\star$)

Show analytically that the pruned second-order state vector remains bounded: (a) the first-order state $\mathbf{x}_t^{(1)}$ converges in mean square because $\rho(C) < 1$ (from Blanchard–Kahn); (b) the second-order correction $\mathbf{x}_t^{(2)}$ is driven by $(\mathbf{x}_{t-1}^{(1)})^2$, which has bounded variance $\text{Var}[(\mathbf{x}^{(1)})^2] = 2(\text{Var}[\mathbf{x}^{(1)}])^2$; (c) the second-order state therefore has bounded variance. Compare to the unpruned case where $\mathbf{x}_t^{(2)}$ can be driven by its own lagged square — potentially explosive.

---

## 29.8 Chapter Summary

**Key results:**

- First-order perturbation satisfies **certainty equivalence** — the stochastic solution equals the deterministic one. Second-order perturbation breaks this: the stochastic steady state is shifted by a term proportional to $\sigma^2$ (shock variance).
- The **second-order policy function** has the form $c_t = c^{SS} + g_x\mathbf{x}_t + \frac{1}{2}\mathbf{x}_t'G_{xx}\mathbf{x}_t + \frac{1}{2}g_{\varepsilon\varepsilon}\sigma^2$, with a **precautionary saving term** $g_{\varepsilon\varepsilon}\sigma^2/2 < 0$ for CRRA households.
- The **second-order coefficient matrices** $G_{xx}$ are computed by solving a linear system involving second derivatives of equilibrium conditions (via ForwardDiff or analytical expressions).
- **Pruning** (Kim et al., 2008) separates first- and second-order state components to prevent explosive paths while retaining second-order accuracy.
- **Risk premia** and **welfare comparisons** across policy regimes require second-order accuracy; first-order models give equity premia of exactly zero and cannot rank policies by welfare.
- In Dynare: `order=2` and `order=3` activate higher-order perturbation; `pruning` option implements the Kim et al. scheme.

*Next: Chapter 30 — Bayesian Estimation of DSGE Models*
