# Chapter 41: Model Validation and Sensitivity Analysis

*Sobol Indices, Identification, and DSGE–BVAR Comparison*

> *"A model that fits every observed fact is likely fitting noise. A model that fits only a few key facts may be capturing the true structure."*

**Cross-reference:** *Principles* Appendix B (model evaluation: Bayes factors, DSGE vs. BVAR); Ch. 39 (epistemic humility in macroeconomics) **[P:AppB, P:Ch.39]**

---

## 41.1 What Does It Mean to Validate a Macroeconomic Model?

The Lucas critique (Chapter 18) tells us that reduced-form parameters change with policy — so reduced-form fit is not evidence of structural validity. But structural identification requires assumptions that are themselves difficult to test. This creates the fundamental tension in macroeconomic model evaluation: the models that are most useful for policy are least testable.

Despite this, model validation is both possible and necessary. We validate along three dimensions:

1. **Internal consistency:** Does the model satisfy its own theoretical restrictions (Blanchard–Kahn conditions, positive wealth effects, etc.)?
2. **External fit:** Do the model's predictions for observable moments match the data? Do the model's IRFs align with SVAR estimates?
3. **Robustness:** Are the model's policy conclusions stable to reasonable changes in parameterization?

This chapter develops formal tools for each dimension.

---

## 41.2 Moment Matching and DSGE vs. VAR Comparison

The standard approach to DSGE validation is **moment matching**: comparing model-implied second moments (standard deviations, autocorrelations, cross-correlations) to their data counterparts, all HP-filtered.

**Definition 41.1 (Spectral Score).** The **spectral score** of a DSGE model is:

$$S(\mathcal{M}) = \sum_{j}\left[\frac{\hat{m}_j^{model} - m_j^{data}}{s_j^{data}}\right]^2,$$

where $\hat{m}_j^{model}$ are model-implied moments, $m_j^{data}$ are data moments, and $s_j^{data}$ are their standard errors. Lower score = better fit.

**DSGE–BVAR distance (Sims, 2002):** The Sims test projects the DSGE model into the BVAR space and measures the KL divergence between the DSGE-implied VAR and the unrestricted BVAR. A significant KL divergence indicates the DSGE model is missing important features of the data.

---

## 41.3 Bayesian Model Comparison: Marginal Likelihood

The **marginal likelihood** $p(\mathbf{Y}|\mathcal{M})$ is the gold standard for Bayesian model comparison (Chapter 30). For two models:

$$\ln BF_{12} = \ln p(\mathbf{Y}|\mathcal{M}_1) - \ln p(\mathbf{Y}|\mathcal{M}_2).$$

From the Jeffreys (1961) scale: $BF_{12} > 10$ is "strong" evidence; $BF_{12} > 100$ is "decisive."

**The Laplace approximation** (Chapter 30): $\ln p(\mathbf{Y}|\mathcal{M}) \approx \ell(\hat\theta^{MAP}) + \frac{k}{2}\ln(2\pi) - \frac{1}{2}\ln|\hat{H}|$, where $k$ is the number of parameters. This penalizes model complexity (larger $k$ → lower marginal likelihood, other things equal).

**Practical implementation in Dynare:** The `identification` command computes the Fisher information matrix and checks for rank deficiency. The `estimation` command reports `oo_.MarginalDensity.LaplaceApproximation` — the log marginal likelihood for Bayes factor computation.

---

## 41.4 Global Sensitivity Analysis: Sobol Indices

**Local sensitivity analysis** (partial derivatives of outputs w.r.t. parameters) is efficient but only measures sensitivity in a small neighborhood of the calibrated parameter vector. **Global sensitivity analysis** (GSA) quantifies how much of the model's output variance can be attributed to each parameter, over the entire feasible parameter space.

**Definition 41.2 (Sobol Variance Decomposition).** For a model output $Y = f(\bm\theta)$ with independent parameters $\theta_1, \ldots, \theta_k$:

$$\text{Var}(Y) = \sum_i V_i + \sum_{i<j}V_{ij} + \cdots + V_{12\cdots k},$$

where $V_i = \text{Var}_{\theta_i}(\mathbb{E}_{\bm\theta_{-i}}[Y|\theta_i])$ is the first-order Sobol index variance (main effect of $\theta_i$).

**Definition 41.3 (First-Order Sobol Index).** The **first-order Sobol index** for parameter $\theta_i$:

$$S_i = \frac{V_i}{\text{Var}(Y)} = \frac{\text{Var}_{\theta_i}(\mathbb{E}_{\bm\theta_{-i}}[Y|\theta_i])}{\text{Var}(Y)}.$$

$S_i$ measures the fraction of total output variance attributable to $\theta_i$ alone (marginalizing over all other parameters). $\sum_i S_i \leq 1$, with equality when parameters are non-interacting.

**Theorem 41.1 (Monte Carlo Estimator for Sobol Indices).** The Saltelli (2002) estimator for $S_i$:

$$\hat{S}_i = \frac{1}{N}\sum_{n=1}^N f(\mathbf{B}_n)\left[f(\mathbf{A}^{(i)}_n) - f(\mathbf{A}_n)\right] \bigg/ \hat{V},$$

where $\mathbf{A}, \mathbf{B}$ are two independent $N\times k$ sample matrices, $\mathbf{A}^{(i)}_n$ is matrix $\mathbf{A}$ with column $i$ replaced by column $i$ of $\mathbf{B}$, and $\hat{V} = \text{Var}_n[f(\mathbf{A}_n)]$.

This requires $N(k+2)$ model evaluations — feasible for NK models (each evaluation ≈ 1 ms) but expensive for large DSGE models (each evaluation ≈ 1 s).

In APL, the Sobol estimate uses outer products:

```apl
⍝ APL — Sobol first-order index via Monte Carlo
⎕IO←0 ⋄ ⎕ML←1

sobol_S1 ← {F A B i ← ⍵
    N ← ≢A
    ⍝ A^{(i)}: replace column i of A with column i of B
    A_i ← A
    A_i[;i] ← B[;i]
    ⍝ Evaluate F at A and A^(i)
    fA   ← F¨ ⊂¨↓A      ⍝ F evaluated at each row of A
    fA_i ← F¨ ⊂¨↓A_i    ⍝ F evaluated at each row of A^(i)
    fB   ← F¨ ⊂¨↓B
    ⍝ Saltelli estimator: (1/N) * sum(fB * (fA_i - fA)) / Var(fA)
    V_hat ← (÷N) × +/ (fA - (+/fA)÷N)*2    ⍝ sample variance of fA
    S_i   ← (÷N×V_hat) × fB +.× fA_i - fA
    S_i}
```

---

## 41.5 Identification Analysis: Iskrev (2010)

**Definition 41.4 (Local Identification).** A DSGE model with parameters $\bm\theta$ is **locally identified** at $\bm\theta_0$ if there is a neighborhood $U(\bm\theta_0)$ such that no other $\bm\theta \in U(\bm\theta_0) \setminus \{\bm\theta_0\}$ generates the same population moments (or likelihood).

**Theorem 41.2 (Iskrev Identification Condition).** The DSGE model is locally identified at $\bm\theta_0$ if and only if the **Jacobian matrix** $J(\bm\theta_0) = \partial m(\bm\theta_0)/\partial\bm\theta'$ has full column rank, where $m(\bm\theta)$ is the vector of model-implied moments (spectral density or autocovariances).

*Proof.* If $J$ has full column rank $k$ (the number of parameters), then for small $\Delta\bm\theta \neq 0$: $m(\bm\theta_0 + \Delta\bm\theta) \approx m(\bm\theta_0) + J\Delta\bm\theta \neq m(\bm\theta_0)$ (since $J\Delta\bm\theta \neq 0$ by full rank). Hence nearby parameters generate different moments — local identification. Conversely, if $J$ is rank-deficient, there exists $\Delta\bm\theta \neq 0$ with $J\Delta\bm\theta = 0$ — a direction of parameter variation that leaves moments unchanged (lack of identification). $\square$

**Practical test:** Compute $J$ numerically (column-by-column finite differences of $m(\bm\theta)$ w.r.t. each $\theta_i$) and check `rank(J) == k`. The numerical rank is computed via the SVD: rank equals the number of singular values above a threshold $\tau = \sigma_1\cdot\varepsilon_M\cdot\max(p,k)$ where $\sigma_1$ is the largest singular value and $p$ is the number of moments.

In Dynare, the `identification` command implements the Iskrev (2010) test and reports which parameters are locally unidentified and why (showing which parameter combinations are degenerate).

---

## 41.6 Worked Example: Smets–Wouters Sensitivity Analysis

*Cross-reference: Principles Ch. 27 (Smets–Wouters calibration and estimation)* **[P:Ch.27]**

We conduct a global sensitivity analysis of the NK welfare loss $\mathcal{L} = \text{Var}(\hat\pi) + \lambda_x\text{Var}(\hat{x})$ to the five key NK parameters $(\beta, \kappa, \sigma, \phi_\pi, \phi_y)$.

```python
import numpy as np
from scipy.stats import uniform

def nk_welfare(theta):
    """NK welfare loss as a function of parameters."""
    beta, kappa, sigma, phi_pi, phi_y = theta
    lambda_x = kappa / 6.0
    rho_u, sigma_u = 0.5, 0.01
    
    G0 = np.array([[1, -kappa], [sigma*phi_pi, 1+sigma*phi_y]])
    G1 = np.array([[beta, 0], [-sigma, 1]])
    try:
        A = np.linalg.inv(G0) @ G1
    except: return np.nan
    eigs = np.linalg.eigvals(A)
    if not np.all(np.abs(eigs) > 1): return np.nan
    Omega = np.linalg.solve(np.eye(2) - rho_u*A, np.linalg.inv(G0) @ np.array([1,0]))
    var_u = sigma_u**2/(1-rho_u**2)
    return Omega[0]**2*var_u + lambda_x*Omega[1]**2*var_u

# Saltelli Sobol estimator
np.random.seed(42); N = 2000; k = 5
# Parameter ranges: beta in [0.97,0.999], kappa in [0.05,0.3], sigma in [0.5,3],
#                   phi_pi in [1.1,4], phi_y in [0,1.5]
ranges = [(0.97, 0.999), (0.05, 0.30), (0.50, 3.00), (1.10, 4.00), (0.00, 1.50)]

def sample_params(N):
    A = np.column_stack([np.random.uniform(lo, hi, N) for lo, hi in ranges])
    return A

A_mat = sample_params(N); B_mat = sample_params(N)
fA = np.array([nk_welfare(A_mat[n]) for n in range(N)])
fB = np.array([nk_welfare(B_mat[n]) for n in range(N)])
# Remove NaN entries
valid = ~(np.isnan(fA) | np.isnan(fB))
fA, fB = fA[valid], fB[valid]; A_mat, B_mat = A_mat[valid], B_mat[valid]
N_valid = valid.sum()

V_hat = np.var(fA)
param_names = ['β', 'κ', 'σ', 'φ_π', 'φ_y']
S1 = np.zeros(k); ST = np.zeros(k)

for i in range(k):
    A_i = A_mat.copy(); A_i[:, i] = B_mat[:, i]
    fA_i = np.array([nk_welfare(A_i[n]) for n in range(N_valid)])
    invalid_i = np.isnan(fA_i)
    fA_i[invalid_i] = fA[invalid_i]  # replace NaN with fA value
    S1[i] = np.mean(fB * (fA_i - fA)) / V_hat
    ST[i] = np.mean(fA * (fA - fA_i)) / V_hat

print("Sobol sensitivity indices for NK welfare loss:")
print(f"{'Parameter':<8} {'S1 (main)':<12} {'ST (total)':<12} {'Interpretation'}")
for i, name in enumerate(param_names):
    print(f"{name:<8} {S1[i]:<12.4f} {ST[i]:<12.4f}  {'important' if ST[i]>0.1 else 'minor'}")
print(f"\nSum of S1: {sum(S1):.4f} (≤ 1; gap = interaction effects)")
```

```julia
using GlobalSensitivity, Statistics

function nk_welfare_jl(θ)
    beta, kappa, sigma, phi_pi, phi_y = θ
    lx = kappa/6.0; rho_u = 0.5; sig_u = 0.01
    G0 = [1 -kappa; sigma*phi_pi 1+sigma*phi_y]; G1=[beta 0;-sigma 1]
    det(G0) ≈ 0 && return NaN
    A = inv(G0)*G1
    all(abs.(eigvals(A)).>1) || return NaN
    Omega = (I(2)-rho_u*A)\(inv(G0)*[1;0])
    vu = sig_u^2/(1-rho_u^2)
    Omega[1]^2*vu + lx*Omega[2]^2*vu
end

# Parameter bounds
lb = [0.97, 0.05, 0.5, 1.1, 0.0]
ub = [0.999, 0.30, 3.0, 4.0, 1.5]

# Sobol analysis using GlobalSensitivity.jl
res = gsa(nk_welfare_jl, Sobol(), [[lb[i],ub[i]] for i in 1:5]; N=2000)
println("Sobol S1 (main effects): ", round.(res.S1, digits=4))
println("Sobol ST (total effects): ", round.(res.ST, digits=4))
```

---

## 41.7 Programming Exercises

### Exercise 41.1 (APL — Local Sensitivity)

For the NK model with baseline $\bm\theta_0 = (\beta_0, \kappa_0, \sigma_0, \phi_{\pi,0}, \phi_{y,0}) = (0.99, 0.15, 1.0, 1.5, 0.5)$: (a) compute the Jacobian $J_{ij} = \partial m_i(\bm\theta_0)/\partial\theta_j$ using central finite differences with $h = 10^{-5}$; (b) compute the condition number $\kappa(J)$; (c) identify which parameters are poorly identified (singular values near zero); (d) verify that $\phi_\pi$ and $\phi_y$ are separately identified but that $\kappa$ and $\lambda_x = \kappa/\varepsilon$ are collinear (not separately identified without observing price dispersion).

### Exercise 41.2 (Python — Morris Screening)

Before running the full Sobol analysis (which requires $N(k+2)$ evaluations), use the **Morris method** (elementary effects) to screen out unimportant parameters with only $r(k+1)$ evaluations (typically $r = 10$): (a) generate $r$ Morris trajectories in the $k$-dimensional parameter space; (b) compute the elementary effect $EE_{ij} = [f(\bm\theta + \Delta\mathbf{e}_j) - f(\bm\theta)]/\Delta$ for parameter $j$ at design point $i$; (c) rank parameters by $\bar{\mu}^* = r^{-1}\sum_i |EE_{ij}|$ and flag those with $\bar\mu^* < 0.01 \times \bar\mu^*_{max}$ as unimportant; (d) confirm these match the low-$S_1$ parameters from the Sobol analysis.

### Exercise 41.3 (Julia — Iskrev Identification Check)

Implement the Iskrev (2010) identification test for the 3-parameter NK model ($\kappa, \phi_\pi, \rho_u$): (a) compute the model-implied autocovariance function $m(\bm\theta) = (\gamma_\pi(0), \gamma_\pi(1), \gamma_x(0), \gamma_x(1), \gamma_{\pi x}(0))$; (b) compute $J = \partial m/\partial(\kappa,\phi_\pi,\rho_u)'$ by finite differences; (c) compute the SVD of $J$ and check rank; (d) if rank-deficient, identify the null space vector (which parameter combination is unidentified?).

### Exercise 41.4 — DSGE vs. BVAR Comparison ($\star$)

The **DSGE-BVAR distance** (Sims, 2002) measures how much the DSGE model is forced away from the data by its structural restrictions. (a) Estimate a 3-variable BVAR(4) on $(GDP, \pi, i)$ data. (b) Compute the log marginal likelihood $\ln p(\mathbf{Y}|BVAR)$. (c) Compute $\ln p(\mathbf{Y}|NK-DSGE)$ using the Laplace approximation. (d) The Bayes factor $BF = p(\mathbf{Y}|NK)/p(\mathbf{Y}|BVAR)$ — if $\ln BF < -2$, the BVAR is significantly better; the DSGE is missing important features of the data.

---

## 41.8 Chapter Summary

**Key results:**

- **Moment matching**: spectral score $S(\mathcal{M}) = \sum_j[(\hat{m}_j^{model}-m_j^{data})/s_j^{data}]^2$ — lower is better; compare to VAR-implied moments.
- **Marginal likelihood** for Bayes factors: $\ln BF_{12} = \ln p(\mathbf{Y}|\mathcal{M}_1) - \ln p(\mathbf{Y}|\mathcal{M}_2)$; Laplace approximation penalizes model complexity.
- **Sobol first-order indices** $S_i = V_i/\text{Var}(Y)$ decompose output variance into parameter contributions; Saltelli estimator (Theorem 41.1) requires $N(k+2)$ evaluations.
- **Iskrev identification** (Theorem 41.2): model is locally identified iff the Jacobian $\partial m/\partial\bm\theta'$ has full column rank; checked via SVD.
- For the NK model: $\phi_\pi$ and $\kappa$ drive most of the welfare loss variance ($S_1 > 0.3$); $\beta$ is nearly invariant ($S_1 \approx 0$) over reasonable ranges.
- In APL: Sobol estimator is `(÷N×V_hat) × fB +.× fA_i - fA`; Jacobian is `{(welfare theta+h×⍵) - welfare theta-h×⍵) ÷ 2×h}¨ identity_columns`.

*Next: Chapter 42 — Capstone Project*
