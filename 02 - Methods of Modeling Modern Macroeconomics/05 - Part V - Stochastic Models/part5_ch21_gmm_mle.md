# Chapter 21: GMM and Maximum Likelihood

*Estimating Structural Macroeconomic Parameters*

> *"GMM is honest about what it assumes. MLE is efficient about what it gains. The choice between them is a trade-off between robustness and precision."*

**Cross-reference:** *Principles* Appendix B (calibration vs. estimation, structural DSGE estimation); Ch. 27 (calibration of the RBC model) **[P:AppB, P:Ch.27]**

---

## 21.1 The Estimation Problem in Structural Macroeconomics

Macroeconomic models contain parameters — the CRRA coefficient $\sigma$, the Calvo price-setting probability $\theta$, the Taylor rule coefficients $(\phi_\pi, \phi_y)$, the AR(1) persistence $\rho_A$ — that cannot be read directly from data. They must be **estimated** from observable moments of the data.

Three estimation strategies are used in modern macroeconomics:

1. **Calibration** (Chapter 17): set parameters equal to values from microeconomic evidence, national accounts, or prior literature. Simple and transparent, but ignores the information in the aggregate data.

2. **GMM** (this chapter): choose parameters to make the model's implied moments match the data moments. Robust to distributional misspecification; requires only moment conditions.

3. **Maximum likelihood / Bayesian** (this chapter and Chapter 30): evaluate the model's full probability density for the observed data and maximize it. Efficient (uses all information in the data) but requires correct distributional specification.

This chapter develops GMM and MLE in full generality, with applications to the consumption Euler equation (GMM) and state-space DSGE models (MLE via the Kalman filter from Chapter 20).

---

## 21.2 The GMM Framework

### 21.2.1 Moment Conditions

**Definition 21.1 (Moment Conditions).** A model with parameter vector $\bm\theta \in \mathbb{R}^k$ implies **moment conditions**:

$$\mathbb{E}[f(\mathbf{x}_t, \bm\theta)] = \mathbf{0},$$

where $f: \mathbb{R}^d\times\mathbb{R}^k \to \mathbb{R}^q$ is a vector-valued function of the data $\mathbf{x}_t$ and parameters. There are $q$ moment conditions and $k$ parameters; the model is:
- **Just-identified** if $q = k$.
- **Over-identified** if $q > k$ (more moment conditions than parameters — testable).
- **Under-identified** if $q < k$ (not estimable from moment conditions alone).

**Example (Consumption Euler equation):** The Euler equation $u'(c_t) = \beta(1+r_{t+1})\mathbb{E}_t[u'(c_{t+1})]$ implies, with CRRA utility:

$$\mathbb{E}_t\left[\beta(1+r_{t+1})\left(\frac{c_{t+1}}{c_t}\right)^{-\sigma} - 1\right] = 0.$$

This is an orthogonality condition: the forecast error $e_{t+1} \equiv \beta(1+r_{t+1})(c_{t+1}/c_t)^{-\sigma} - 1$ should be uncorrelated with any variable in the time-$t$ information set $\mathcal{F}_t$. Using instruments $z_t \in \mathcal{F}_t$ (e.g., $z_t = (1, c_{t-1}/c_{t-2}, r_t)'$):

$$\mathbb{E}[f_t(\bm\theta)] = \mathbb{E}[e_{t+1}(\bm\theta)\cdot z_t] = \mathbf{0}.$$

With $z_t$ having $q$ elements, we have $q$ moment conditions for $k = 2$ parameters $(\beta, \sigma)$.

### 21.2.2 The GMM Estimator

**Definition 21.2 (GMM Estimator).** Let $\bar{f}_T(\bm\theta) = T^{-1}\sum_{t=1}^T f(\mathbf{x}_t, \bm\theta)$ be the sample analogue of the moment conditions. The **GMM estimator** minimizes the weighted quadratic form:

$$\hat{\bm\theta}_{GMM} = \arg\min_{\bm\theta}\;\bar{f}_T(\bm\theta)'\,W\,\bar{f}_T(\bm\theta),$$

where $W$ is a $q\times q$ positive definite weight matrix.

For just-identified models ($q = k$), GMM sets $\bar{f}_T(\hat{\bm\theta}) = \mathbf{0}$ and the weight matrix is irrelevant.

**Theorem 21.1 (Consistency of GMM).** Under standard regularity conditions (identification, stationarity, uniform convergence):

$$\hat{\bm\theta}_{GMM} \xrightarrow{p} \bm\theta_0 \quad \text{as } T\to\infty.$$

*Proof sketch.* The population objective $\bar{f}(\bm\theta)'W\bar{f}(\bm\theta)$ has a unique minimum at $\bm\theta_0$ (identification). Uniform convergence of $\bar{f}_T(\bm\theta) \to \mathbb{E}[f(\mathbf{x}_t,\bm\theta)]$ and continuity of the objective give consistency by the argmax continuous mapping theorem. $\square$

**Theorem 21.2 (Asymptotic Normality of GMM).** Under additional smoothness conditions:

$$\sqrt{T}(\hat{\bm\theta}_{GMM} - \bm\theta_0) \xrightarrow{d} \mathcal{N}(\mathbf{0}, V_{GMM}),$$

where the asymptotic variance is:

$$V_{GMM} = (D'WD)^{-1}(D'W\,S\,WD)(D'WD)^{-1},$$

with $D = \partial\mathbb{E}[f_t(\bm\theta_0)]/\partial\bm\theta'$ (the derivative of the moment conditions) and $S = \sum_{j=-\infty}^\infty\text{Cov}(f_t, f_{t-j})$ (the long-run covariance of the moment conditions, estimated by Newey–West).

### 21.2.3 Optimal GMM and the Efficient Weight Matrix

**Theorem 21.3 (Efficiency of Optimal GMM).** The asymptotic variance $V_{GMM}$ is minimized (in the matrix sense) when $W = S^{-1}$, the inverse of the long-run covariance matrix. The **efficient GMM** estimator with $W^* = S^{-1}$ achieves:

$$V_{GMM}^* = (D'S^{-1}D)^{-1},$$

which is the semiparametric efficiency bound — no consistent and asymptotically normal estimator based only on the moment conditions $\mathbb{E}[f_t(\bm\theta)] = \mathbf{0}$ can achieve lower asymptotic variance.

*Proof.* For any weight matrix $W$: $V_{GMM} - V_{GMM}^* = (D'WD)^{-1}D'W(S - D(D'S^{-1}D)^{-1}D')WD(D'WD)^{-1} \geq 0$ by the positive semi-definiteness of $S - D(D'S^{-1}D)^{-1}D'$ (a Schur complement). $\square$

**Practical two-step GMM:**

1. **Step 1:** Estimate $\hat{\bm\theta}^{(1)}$ using $W = I$ (identity weight matrix).
2. Compute the Newey–West estimator of $S$: $\hat{S} = \hat\Gamma_0 + \sum_{j=1}^{B}(1-j/(B+1))(\hat\Gamma_j + \hat\Gamma_j')$, where $\hat\Gamma_j = T^{-1}\sum_{t=j+1}^T\hat f_t\hat f_{t-j}'$.
3. **Step 2:** Re-estimate $\hat{\bm\theta}^{(2)}$ using $W = \hat{S}^{-1}$.

---

## 21.3 Application: Estimating the Euler Equation

**Moment conditions for the consumption Euler equation:**

With $f_{t+1}(\beta,\sigma) = \beta(1+r_{t+1})(c_{t+1}/c_t)^{-\sigma} - 1$ and instruments $z_t = (1, \Delta\ln c_{t-1}, r_t)'$:

$$\mathbb{E}[f_{t+1}(\beta,\sigma)\cdot z_t] = \mathbf{0}.$$

This gives $q = 3$ moment conditions for $k = 2$ parameters — an over-identified system. The $J$-test for over-identifying restrictions:

$$J = T\cdot\bar{f}_T(\hat{\bm\theta})'\hat{S}^{-1}\bar{f}_T(\hat{\bm\theta}) \xrightarrow{d} \chi^2(q-k) = \chi^2(1).$$

The $J$-statistic tests whether the over-identifying restrictions are satisfied. Rejection suggests model misspecification (the Euler equation moment conditions are inconsistent with the data given the instruments).

In APL, the GMM objective `g'Wg` is:

```apl
⍝ APL — GMM objective and gradient
⎕IO←0 ⋄ ⎕ML←1

⍝ Given: moment conditions g_bar (q-vector), weight matrix W (q×q)
gmm_obj ← {g W ← ⍵ ⋄ g +.× W +.× g}    ⍝ g'Wg

⍝ Moment conditions: f_t(β, σ) = β(1+r_{t+1})(c_{t+1}/c_t)^{-σ} - 1
⍝ Instruments Z: T×q matrix
moment_cond ← {beta sigma c_growth r_next Z ← ⍵
    e ← beta × (1+r_next) × c_growth*(-sigma)) - 1   ⍝ Euler error (T-vector)
    g_bar ← (÷≢e) × +⌿ Z × e ∘.× ⍬    ⍝ gbar = mean(Z × e)
    ⍝ = (1/T) Z'e  [q-vector]
    (÷T) × (⍉Z) +.× e}

⍝ Two-step GMM
⍝ Step 1: W = I
⍝ Step 2: W = S^{-1} (Newey-West)
newey_west ← {e Z B ← ⍵
    T ← ≢e
    e_z ← Z × e ∘.× ⍬    ⍝ T×q: Z_t × e_t
    ⍝ Gamma_0
    G0 ← (T⍴1) {(⍉⍺)+.×⍺} e_z  ⍝ E[e_z e_z']... need matrix version
    ⍝ APL: G0 = (1/T) × (⍉e_z) +.× e_z
    G0 ← (÷T) × (⍉e_z) +.× e_z
    ⍝ Add Bartlett-kernel lags
    S ← G0
    :For j :In 1+⍳B
        Gj ← (÷T) × (⍉j↓e_z) +.× (j↓⍵e_z)    ⍝ Gamma_j
        S +← (1-j÷B+1)×Gj + ⍉Gj
    :EndFor
    S}
```

```python
import numpy as np
from scipy.optimize import minimize

# GMM estimation of consumption Euler equation parameters (β, σ)
np.random.seed(42)
T = 300
true_beta, true_sigma = 0.97, 1.5

# Simulate data from the model
r = 0.015 + 0.005*np.random.randn(T)  # quarterly interest rate
c_growth = (true_beta*(1+r))**(1/true_sigma) * np.exp(0.02*np.random.randn(T))

# Instruments: constant, lagged c_growth, lagged r
Z = np.column_stack([np.ones(T-1), c_growth[:-1], r[:-1]])
c_g = c_growth[1:]  # consumption growth in t+1
r_t = r[1:]         # interest rate in t+1

def moment_conditions(params, c_g, r_t, Z):
    beta, sigma = params
    e = beta * (1+r_t) * c_g**(-sigma) - 1
    g_bar = Z.T @ e / len(e)   # q-vector of sample moments
    return g_bar, e

def gmm_objective(params, W, c_g, r_t, Z):
    g, _ = moment_conditions(params, c_g, r_t, Z)
    return g @ W @ g

def newey_west(e, Z, B=4):
    T = len(e); q = Z.shape[1]
    ez = Z * e[:,None]  # T×q
    S = ez.T @ ez / T
    for j in range(1, B+1):
        Gj = ez[j:].T @ ez[:-j] / T
        S += (1 - j/(B+1)) * (Gj + Gj.T)
    return S

# Step 1: W = I
W1 = np.eye(Z.shape[1])
res1 = minimize(gmm_objective, [0.96, 1.0], args=(W1, c_g, r_t, Z),
                method='Nelder-Mead', options={'xatol':1e-8})
beta1, sigma1 = res1.x

# Step 2: W = S^{-1}
_, e1 = moment_conditions(res1.x, c_g, r_t, Z)
S_hat = newey_west(e1, Z)
W2 = np.linalg.inv(S_hat)
res2 = minimize(gmm_objective, res1.x, args=(W2, c_g, r_t, Z),
                method='Nelder-Mead', options={'xatol':1e-10})
beta2, sigma2 = res2.x

print(f"True:  β={true_beta:.3f}, σ={true_sigma:.3f}")
print(f"Step1: β={beta1:.3f}, σ={sigma1:.3f}")
print(f"Step2: β={beta2:.3f}, σ={sigma2:.3f}")

# J-test for over-identification
g_hat, _ = moment_conditions(res2.x, c_g, r_t, Z)
J_stat = len(c_g) * g_hat @ W2 @ g_hat
q, k = Z.shape[1], 2
print(f"\nJ-statistic = {J_stat:.3f}, df = {q-k}, p-value ≈ {1-0.95:.3f} CV at 5%")
```

```julia
using LinearAlgebra, Optim

# GMM for Euler equation
function gmm_euler(params, c_g, r_next, Z, W)
    beta, sigma = params
    e = @. beta * (1+r_next) * c_g^(-sigma) - 1
    g = Z' * e / length(e)
    return g' * W * g
end

# Simulate data
T=200; Random.seed!(42)
beta_t, sigma_t = 0.97, 1.5
r = 0.015 .+ 0.005.*randn(T); cg = (beta_t.*(1.+r)).^(1/sigma_t).*exp.(0.02.*randn(T))
Z = hcat(ones(T-1), cg[1:end-1], r[1:end-1]); cg_use=cg[2:end]; r_use=r[2:end]

W1 = I(size(Z,2))
res1 = optimize(p->gmm_euler(p,cg_use,r_use,Z,W1), [0.95,1.0], NelderMead())
println("Step 1: β=$(round(Optim.minimizer(res1)[1],digits=3)), σ=$(round(Optim.minimizer(res1)[2],digits=3))")
```

```r
# R GMM using gmm package
library(gmm)

set.seed(42); T<-200; beta_t<-0.97; sigma_t<-1.5
r<-0.015+0.005*rnorm(T); cg<-(beta_t*(1+r))^(1/sigma_t)*exp(0.02*rnorm(T))

# Moment conditions: g(θ,x) = Z_t × [β(1+r_{t+1})(c_{t+1}/c_t)^{-σ} - 1]
dat <- data.frame(cg=cg[-1], r=r[-1], cg_lag=cg[-T], r_lag=r[-T])
g_euler <- function(theta, x) {
  beta<-theta[1]; sigma<-theta[2]
  e <- beta*(1+x$r)*(x$cg^(-sigma)) - 1
  cbind(e, e*x$cg_lag, e*x$r_lag)
}
fit <- gmm(g_euler, dat, c(beta=0.96, sigma=1.0), method="BFGS")
print(summary(fit))
```

---

## 21.4 Maximum Likelihood Estimation

**Definition 21.3 (Likelihood Function).** For independent and identically distributed observations $\{x_t\}_{t=1}^T$ with density $f(x; \bm\theta)$, the **likelihood function** is:

$$L(\bm\theta) = \prod_{t=1}^T f(x_t; \bm\theta), \quad \ell(\bm\theta) = \sum_{t=1}^T\ln f(x_t; \bm\theta).$$

The **MLE** is $\hat{\bm\theta}_{MLE} = \arg\max_{\bm\theta}\ell(\bm\theta)$.

**Theorem 21.4 (Asymptotic Theory of MLE).** Under regularity conditions (identifiability, compactness, continuity, dominated convergence):

$$\sqrt{T}(\hat{\bm\theta}_{MLE} - \bm\theta_0) \xrightarrow{d} \mathcal{N}(\mathbf{0}, I(\bm\theta_0)^{-1}),$$

where $I(\bm\theta) = -\mathbb{E}[\partial^2\ln f(x;\bm\theta)/\partial\bm\theta\partial\bm\theta'] = \mathbb{E}[(\partial\ln f/\partial\bm\theta)(\partial\ln f/\partial\bm\theta)']$ is the **Fisher information matrix**.

MLE is **asymptotically efficient**: no consistent and asymptotically normal estimator achieves a smaller asymptotic variance than $I(\bm\theta_0)^{-1}$ (Cramér–Rao lower bound).

### 21.4.1 The Score Function and Information Matrix

**Definition 21.4 (Score).** The **score function** is the gradient of the log-likelihood:

$$s(\bm\theta; x) = \frac{\partial\ln f(x;\bm\theta)}{\partial\bm\theta}.$$

At the true parameter $\bm\theta_0$: $\mathbb{E}[s(\bm\theta_0; x)] = \mathbf{0}$ (score has zero mean) and $\text{Var}[s(\bm\theta_0; x)] = I(\bm\theta_0)$ (variance equals the Fisher information).

The MLE FOC: $\partial\ell(\hat{\bm\theta})/\partial\bm\theta = \mathbf{0}$ sets the sample score to zero.

### 21.4.2 The Delta Method

The **delta method** computes standard errors for nonlinear transformations of estimators.

**Theorem 21.5 (Delta Method).** If $\sqrt{T}(\hat{\bm\theta} - \bm\theta_0) \xrightarrow{d} \mathcal{N}(\mathbf{0}, V)$ and $g: \mathbb{R}^k\to\mathbb{R}^m$ is differentiable at $\bm\theta_0$, then:

$$\sqrt{T}(g(\hat{\bm\theta}) - g(\bm\theta_0)) \xrightarrow{d} \mathcal{N}\left(\mathbf{0},\; \nabla g(\bm\theta_0)\, V\, \nabla g(\bm\theta_0)'\right).$$

*Application:* If $\hat\sigma = 1/\hat\rho$ (EIS estimated from the inverse of the CRRA), the delta method gives $\text{SE}(\hat\sigma) \approx |{-1/\hat\rho^2}|\cdot\text{SE}(\hat\rho)$.

---

## 21.5 Structural vs. Reduced-Form Estimation: The Lucas Critique

**Definition 21.5 (Lucas Critique).** The **Lucas critique** [P:Ch.16.1] states that the parameters of reduced-form econometric models — consumption functions, investment equations, Phillips curves — are not structural: they depend on the policy rule in place and will change when policy changes. Only structural parameters (preferences, technology) are invariant to policy changes.

**Implication for GMM:** Estimating the reduced-form Euler equation via GMM yields consistent estimates of $(\beta, \sigma)$ as long as the regression structure is correctly specified. But reduced-form VAR coefficients (slopes of IS curves, Phillips curves estimated in the data) are contaminated by the policy rule and will change when the policy rule changes — exactly the Lucas critique.

**Solution:** Estimate structural parameters (the deep parameters of preferences and technology) either by:
1. **GMM using structural moment conditions** — orthogonality conditions derived directly from the structural model (as in Section 21.3).
2. **Full-information MLE** — evaluate the likelihood of the entire DSGE model (using the Kalman filter from Chapter 20) rather than single equations.

---

## 21.6 Programming Exercises

### Exercise 21.1 (APL — GMM Objective)

Implement the two-step GMM estimator in APL. (a) The moment conditions function `g_bar ← {(÷T)×(⍉Z)+.×e}` where `e` is the vector of Euler errors. (b) The GMM objective `Q ← {g W ← ⍵ ⋄ g+.×W+.×g}`. (c) The Newey–West estimator of $S$ using a loop with `+←` accumulation. (d) The two-step estimator via `⍣≡` Newton iterations on the gradient of $Q$.

### Exercise 21.2 (Python — MLE vs. GMM Comparison)

For the CRRA Euler equation, compare MLE (assuming log-normal consumption growth: $\ln(c_{t+1}/c_t) \sim \mathcal{N}(\mu_c, \sigma_c^2)$) to two-step GMM. (a) Derive the MLE estimating equations. (b) Simulate 100 datasets of length $T = 200$ from the model. (c) For each dataset, compute the MLE, Step-1 GMM, and Step-2 GMM estimates. (d) Compare bias, variance, and coverage of 95% confidence intervals across the 100 replications.

### Exercise 21.3 (Julia — J-Test Power)

```julia
# Power of the J-test: does it detect model misspecification?
using Distributions, Statistics, LinearAlgebra, Optim

function simulate_and_test(T, beta_t, sigma_t, spec_error=0.0)
    # Generate data, possibly with misspecification
    r = 0.015 .+ 0.005.*randn(T)
    c_growth = (beta_t.*(1.+r)).^(1/sigma_t) .* exp.(0.02.*randn(T))
    
    # Add specification error: habit persistence
    if spec_error > 0
        c_growth .*= (1 .+ spec_error.*[0; c_growth[1:end-1]])
    end
    
    # GMM with 3 instruments
    Z = hcat(ones(T-1), c_growth[1:end-1], r[1:end-1])
    cg = c_growth[2:end]; rv = r[2:end]
    
    # Step 1
    obj1(p) = (Z'*(p[1]*(1.+rv).*cg.^(-p[2]).-1)/length(cg))' * (Z'*(p[1]*(1.+rv).*cg.^(-p[2]).-1)/length(cg))
    res1 = optimize(obj1, [0.96, 1.0], NelderMead())
    p1 = Optim.minimizer(res1)
    
    # Newey-West S
    e1 = p1[1].*(1.+rv).*cg.^(-p1[2]).-1
    ez = Z .* e1; T_=length(e1); q=size(Z,2)
    S = ez'*ez/T_
    for j in 1:4; S += (1-j/5)*(ez[j+1:end,:]'*ez[1:end-j,:]/T_ + ez[1:end-j,:]'*ez[j+1:end,:]/T_); end
    W2 = inv(S)
    
    # Step 2
    obj2(p) = (Z'*(p[1]*(1.+rv).*cg.^(-p[2]).-1)/length(cg))' * W2 * (Z'*(p[1]*(1.+rv).*cg.^(-p[2]).-1)/length(cg))
    res2 = optimize(obj2, p1, NelderMead())
    p2 = Optim.minimizer(res2)
    
    # J-test
    g = Z'*(p2[1].*(1.+rv).*cg.^(-p2[2]).-1)/length(cg)
    J = length(cg) * g'*W2*g
    return J > quantile(Chisq(q-2), 0.95)  # reject at 5%?
end

# Size (no misspecification): rejection rate should be ~5%
rejections_size = mean([simulate_and_test(200, 0.97, 1.5, 0.0) for _ in 1:500])
# Power (with misspecification)
rejections_power = mean([simulate_and_test(200, 0.97, 1.5, 0.3) for _ in 1:500])
println("J-test size:  $(round(rejections_size*100,digits=1))% (target: 5%)")
println("J-test power: $(round(rejections_power*100,digits=1))% (with habit misspec)")
```

### Exercise 21.4 — Cramér–Rao Bound ($\star$)

For the AR(1) model $y_t = \rho y_{t-1} + \varepsilon_t$ with $\varepsilon_t \sim \mathcal{N}(0, \sigma^2)$: (a) compute the log-likelihood $\ell(\rho, \sigma^2; \{y_t\})$; (b) compute the score functions $\partial\ell/\partial\rho$ and $\partial\ell/\partial\sigma^2$; (c) compute the Fisher information matrix $I(\rho, \sigma^2)$; (d) show the Cramér–Rao lower bound for $\text{Var}(\hat\rho)$ is $(1-\rho^2)/T$; (e) verify by simulation that the OLS estimate $\hat\rho = \sum y_t y_{t-1}/\sum y_{t-1}^2$ achieves this bound.

---

## 21.7 Chapter Summary

**Key results:**

- The **GMM estimator** $\hat{\bm\theta}_{GMM}$ minimizes $\bar{f}_T(\bm\theta)'W\bar{f}_T(\bm\theta)$; consistent (Theorem 21.1) and asymptotically normal (Theorem 21.2) under weak regularity conditions.
- **Optimal GMM** uses $W^* = S^{-1}$ (inverse long-run covariance), achieving the semiparametric efficiency bound $V^* = (D'S^{-1}D)^{-1}$ (Theorem 21.3).
- The **two-step GMM** procedure: estimate with $W = I$, compute Newey–West $\hat{S}$, re-estimate with $W = \hat{S}^{-1}$.
- The **$J$-test** $J = T\bar{f}'S^{-1}\bar{f} \to \chi^2(q-k)$ tests over-identifying restrictions — model misspecification.
- The **MLE** achieves the Cramér–Rao bound: asymptotic variance $= I(\bm\theta_0)^{-1}$ (Theorem 21.4). More efficient than GMM but requires correct distributional specification.
- The **delta method** (Theorem 21.5): SE of $g(\hat\bm\theta) \approx |\nabla g(\hat\bm\theta)|\cdot\text{SE}(\hat\bm\theta)$.
- The **Lucas critique**: only structural parameters (from Euler equations, not reduced-form regressions) are policy-invariant; GMM on structural moment conditions sidesteps the critique.
- In APL: GMM objective is `g +.× W +.× g`; Newey–West is accumulated via `S +← (1-j÷B+1)×Gj + ⍉Gj`; MLE maximization uses `⍣≡` Newton steps on the score.

---

*Next: Part VI — Numerical Methods for Solving Macroeconomic Models*
