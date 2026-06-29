# Chapter 30: Bayesian Estimation of DSGE Models

*MCMC Methods and the Metropolis–Hastings Algorithm*

> *"Bayesian estimation takes the model seriously as a data-generating process, while acknowledging that we are uncertain about its parameters."*
> — Frank Smets and Raf Wouters

**Cross-reference:** *Principles* Appendix B (Bayesian estimation, Metropolis–Hastings); Ch. 27 (Smets–Wouters calibration and estimation) **[P:AppB, P:Ch.27]**

---

## 30.1 Why Bayesian Estimation?

Chapter 21 developed GMM and MLE for structural parameter estimation. For DSGE models, **Bayesian estimation** — combining the likelihood with prior distributions over parameters — has become the dominant approach for three reasons:

**1. Identification.** DSGE models typically have many parameters but limited observables. Many parameter combinations are nearly observationally equivalent — different parameter values produce similar time series. The prior provides regularization, guiding the posterior toward economically plausible regions of the parameter space.

**2. Incorporating prior information.** Microeconomic studies provide independent evidence about preference and technology parameters (e.g., labor supply elasticities from household data, price rigidity from firm surveys). The Bayesian prior formalizes this information.

**3. Model comparison.** The Bayesian framework provides a principled method for comparing non-nested models via the **Bayes factor** — the ratio of marginal likelihoods. This is particularly useful when comparing different DSGE specifications.

The price: Bayesian inference requires specifying priors and running MCMC chains that can take hours to converge. For large models (Smets–Wouters with 41 estimated parameters), each likelihood evaluation takes seconds, and the chain needs $10^5$–$10^6$ draws. In practice, this is done using Dynare's optimized implementation.

---

## 30.2 Bayesian Inference: The Posterior Distribution

**Bayes' theorem for parameters:**

$$p(\bm\theta | \mathbf{Y}_{1:T}) = \frac{p(\mathbf{Y}_{1:T} | \bm\theta)\cdot p(\bm\theta)}{p(\mathbf{Y}_{1:T})},$$

where:
- $p(\bm\theta|\mathbf{Y})$ — **posterior**: our beliefs about $\bm\theta$ after seeing the data.
- $p(\mathbf{Y}|\bm\theta) = \mathcal{L}(\bm\theta)$ — **likelihood**: probability of the data given parameters (from Kalman filter, Chapter 20).
- $p(\bm\theta)$ — **prior**: beliefs before seeing the data.
- $p(\mathbf{Y}) = \int p(\mathbf{Y}|\bm\theta)p(\bm\theta)d\bm\theta$ — **marginal likelihood**: normalizing constant.

The posterior combines the prior with the likelihood. When the likelihood is very informative, the posterior is concentrated near the MLE regardless of the prior. When the data are sparse, the posterior is closer to the prior.

**Definition 30.1 (Posterior Mode and Posterior Mean).** The **posterior mode** $\hat{\bm\theta}^{MAP}$ (maximum a posteriori) maximizes $\ln p(\bm\theta|\mathbf{Y}) = \ell(\bm\theta) + \ln p(\bm\theta) + \text{const}$. The **posterior mean** $\bar{\bm\theta} = \mathbb{E}[\bm\theta|\mathbf{Y}] = \int\bm\theta\, p(\bm\theta|\mathbf{Y})d\bm\theta$ is the Bayes estimator under squared error loss.

Finding the posterior mode is an optimization problem (Chapter 24); computing the posterior mean requires integrating over the posterior distribution — which is analytically intractable for DSGE models. **MCMC methods** approximate the posterior by sampling.

---

## 30.3 Prior Distributions for DSGE Parameters

The choice of prior reflects both economic beliefs and practical considerations.

**Standard choices:**

| Parameter | Typical Prior | Support | Motivation |
|---|---|---|---|
| Calvo probability $\theta$ | Beta($\alpha,\beta$) | $[0,1]$ | Fraction, must be in $[0,1]$ |
| CRRA $\sigma$ | Gamma($\alpha,\beta$) | $(0,\infty)$ | Must be positive |
| AR persistence $\rho$ | Beta($\alpha,\beta$) | $[0,1]$ | Must be in $[0,1]$ for stationarity |
| Taylor rule $\phi_\pi$ | Normal (truncated) | $[1,\infty)$ | Determinacy requires $\phi_\pi > 1$ |
| Shock std dev $\sigma_i$ | Inverse-Gamma | $(0,\infty)$ | Scale parameter, must be positive |
| Investment adj. cost $\psi$ | Normal (positive) | $(0,\infty)$ | Must be positive |

**Definition 30.2 (Conjugate Prior).** A prior $p(\bm\theta)$ is **conjugate** to the likelihood $p(\mathbf{Y}|\bm\theta)$ if the posterior $p(\bm\theta|\mathbf{Y})$ is in the same distributional family as the prior. Conjugate priors allow analytical posterior updates.

For DSGE models, the full likelihood is nonlinear in parameters (through the Kalman filter), so no conjugate priors exist. MCMC is necessary.

---

## 30.4 The Metropolis–Hastings Algorithm

The **Metropolis–Hastings (MH) algorithm** generates samples from the posterior distribution $p(\bm\theta|\mathbf{Y})$ without computing the normalizing constant $p(\mathbf{Y})$.

**Algorithm 30.1 (Metropolis–Hastings).**

Initialize: $\bm\theta^{(0)}$ (e.g., the posterior mode from optimization). Scale matrix $\Sigma_q$ (proposal covariance, typically tuned to give 20–30% acceptance rate).

For $s = 1, 2, \ldots, N$:

1. **Propose:** Draw $\bm\theta^* \sim q(\bm\theta | \bm\theta^{(s-1)}) = \mathcal{N}(\bm\theta^{(s-1)}, \Sigma_q)$ (symmetric random walk proposal).

2. **Compute acceptance ratio:**
$$\alpha(\bm\theta^{(s-1)}, \bm\theta^*) = \min\!\left\{1,\; \frac{p(\bm\theta^*|\mathbf{Y})}{p(\bm\theta^{(s-1)}|\mathbf{Y})}\right\} = \min\!\left\{1,\; \frac{\mathcal{L}(\bm\theta^*)p(\bm\theta^*)}{\mathcal{L}(\bm\theta^{(s-1)})p(\bm\theta^{(s-1)})}\right\}.$$

3. **Accept/Reject:** Draw $u \sim \text{Uniform}[0,1]$.
   - If $u \leq \alpha$: $\bm\theta^{(s)} = \bm\theta^*$ (accept).
   - If $u > \alpha$: $\bm\theta^{(s)} = \bm\theta^{(s-1)}$ (reject, stay).

**Theorem 30.1 (MH Invariance).** The MH chain with target $\pi(\bm\theta) \propto \mathcal{L}(\bm\theta)p(\bm\theta)$ has $\pi$ as its stationary distribution. Starting from any initial point $\bm\theta^{(0)}$, the chain converges in distribution to $\pi$ under mild regularity conditions (irreducibility, aperiodicity).

*Proof sketch.* The **detailed balance condition** ensures $\pi$ is stationary. For any two states $\bm\theta, \bm\theta^*$:

$$\pi(\bm\theta)q(\bm\theta^*|\bm\theta)\min\!\left\{1,\frac{\pi(\bm\theta^*)}{\pi(\bm\theta)}\right\} = \pi(\bm\theta^*)q(\bm\theta|\bm\theta^*)\min\!\left\{1,\frac{\pi(\bm\theta)}{\pi(\bm\theta^*)}\right\}.$$

This is verified by considering the two cases $\pi(\bm\theta^*) \geq \pi(\bm\theta)$ and $\pi(\bm\theta^*) < \pi(\bm\theta)$. Detailed balance $\Rightarrow$ stationarity. Irreducibility and aperiodicity $\Rightarrow$ convergence. $\square$

In APL, the MH chain is a scan over draws:

```apl
⍝ APL — Metropolis-Hastings chain
⎕IO←0 ⋄ ⎕ML←1

⍝ Scalar log-posterior function (for illustration)
log_post ← {theta ← ⍵
    ⍝ log-likelihood (Kalman filter result — evaluated via ⎕PY in production)
    ⍝ + log-prior
    ll ← -0.5×(theta-1.5)*2    ⍝ placeholder: Gaussian log-likelihood
    lp ← -0.5×(theta-1.0)*2    ⍝ Normal(1, 1) prior
    ll + lp}

⍝ One MH step: theta is current draw, scale is proposal std dev
mh_step ← {theta scale ← ⍵
    proposal ← theta + scale × (2○(?0)×(-2×⍟?0)*0.5)  ⍝ normal draw (Box-Muller)
    log_ratio ← (log_post proposal) - log_post theta
    accept ← log_ratio > ⍟?0    ⍝ log(U) < log_ratio ⟺ U < ratio
    accept ⊃ proposal theta}    ⍝ return proposal if accepted, else theta

⍝ Run chain: N_draws iterations
N_draws ← 10000
scale   ← 0.5     ⍝ proposal standard deviation
theta0  ← 0       ⍝ starting value

⍝ Collect chain via scan (each step takes [theta, scale] and returns new theta)
chain ← {mh_step ⍵ scale}\ N_draws ⍴ theta0

⍝ Compute posterior statistics
burnin  ← 2000
samples ← burnin↓chain
post_mean ← (+/samples)÷≢samples
post_std  ← ((+/(samples-post_mean)*2)÷≢samples)*0.5
accept_rate ← (+/chain[1:]≠chain[:-1])÷N_draws-1

post_mean post_std accept_rate    ⍝ display results
```

---

## 30.5 MCMC Diagnostics

Running an MCMC chain is only half the task — the other half is verifying that the chain has **converged** to the target distribution and that the draws provide reliable posterior estimates.

### 30.5.1 Trace Plots and Mixing

The trace plot shows the parameter value at each iteration. A well-mixing chain moves freely around the parameter space, not getting stuck in one region. Poor mixing (high autocorrelation in the chain) suggests the proposal variance $\Sigma_q$ is too small.

**Definition 30.3 (Effective Sample Size).** The **effective sample size** accounts for autocorrelation in the chain:

$$N_{eff} = \frac{N}{1 + 2\sum_{j=1}^\infty\hat\rho_j},$$

where $\hat\rho_j$ is the sample autocorrelation at lag $j$. If the chain has high autocorrelation, $N_{eff} \ll N$ — 10,000 correlated draws may be worth only 100 independent draws.

### 30.5.2 The Gelman–Rubin Statistic

**Definition 30.4 (Gelman–Rubin $\hat{R}$ Statistic).** Run $m \geq 2$ independent chains from different starting points, each of length $n$. The $\hat{R}$ statistic compares within-chain variance to between-chain variance:

$$\hat{R} = \sqrt{\frac{\hat{V}}{W}},$$

where $W = \frac{1}{m}\sum_j S_j^2$ (average within-chain variance), $B = \frac{n}{m-1}\sum_j(\bar\theta_j - \bar{\bar\theta})^2$ (between-chain variance), and $\hat{V} = (1-1/n)W + B/n$ (pooled variance estimate).

$\hat{R} \approx 1.0$ indicates convergence (within- and between-chain variances agree). The standard threshold: convergence when $\hat{R} < 1.1$ for all parameters.

---

## 30.6 The DSGE Likelihood via the Kalman Filter

The log-likelihood for the DSGE model with parameters $\bm\theta$ is evaluated using the Kalman filter (Chapter 20):

$$\ell(\bm\theta) = -\frac{Tp}{2}\ln(2\pi) - \frac{1}{2}\sum_{t=1}^T\left[\ln|\mathbf{F}_t(\bm\theta)| + \mathbf{v}_t(\bm\theta)'\mathbf{F}_t(\bm\theta)^{-1}\mathbf{v}_t(\bm\theta)\right],$$

where:
1. From $\bm\theta$: construct the DSGE model matrices $\Gamma_0(\bm\theta)$, $\Gamma_1(\bm\theta)$, $\Psi(\bm\theta)$.
2. Solve using gensys (Chapter 28): get decision rules $C(\bm\theta)$, $D(\bm\theta)$.
3. Map to state-space form: $F(\bm\theta)$, $H(\bm\theta)$, $Q(\bm\theta)$, $R(\bm\theta)$.
4. Run Kalman filter: get innovations $\mathbf{v}_t(\bm\theta)$ and innovation variances $\mathbf{F}_t(\bm\theta)$.
5. Sum the log-likelihood contributions.

Each MH iteration requires one full Kalman filter evaluation — the bottleneck of the estimation.

---

## 30.7 Model Comparison: Marginal Likelihood and Bayes Factors

**Definition 30.5 (Marginal Likelihood).** The **marginal likelihood** of model $\mathcal{M}$ is:

$$p(\mathbf{Y}|\mathcal{M}) = \int p(\mathbf{Y}|\bm\theta, \mathcal{M})p(\bm\theta|\mathcal{M})d\bm\theta.$$

The **Bayes factor** for model $\mathcal{M}_1$ vs. $\mathcal{M}_2$:

$$BF_{12} = \frac{p(\mathbf{Y}|\mathcal{M}_1)}{p(\mathbf{Y}|\mathcal{M}_2)}.$$

If $BF_{12} > 1$: $\mathcal{M}_1$ is supported by the data; $BF_{12} > 10$ is "strong" evidence.

**The Laplace approximation** for the marginal likelihood:

$$\ln p(\mathbf{Y}|\mathcal{M}) \approx \ell(\hat{\bm\theta}^{MAP}) + \ln p(\hat{\bm\theta}^{MAP}|\mathcal{M}) + \frac{k}{2}\ln(2\pi) - \frac{1}{2}\ln|\hat{H}|,$$

where $\hat{H}$ is the negative Hessian of the log-posterior at the mode and $k$ is the number of parameters. The Laplace approximation is used in Dynare as the default marginal likelihood estimate.

---

## 30.8 Worked Example: MH for a 3-Parameter NK Model

```python
import numpy as np
from scipy.stats import norm, gamma, beta as beta_dist
from scipy.linalg import solve_discrete_lyapunov

# Simplified 3-parameter NK model
# Parameters: theta = [kappa, phi_pi, sigma_u]
# Prior: kappa ~ Beta(2,10), phi_pi ~ N(1.5, 0.25) truncated [1,∞), sigma_u ~ IG(0.1, 2)

def log_prior(theta):
    kappa, phi_pi, sigma_u = theta
    if kappa <= 0 or kappa >= 1 or phi_pi <= 1 or sigma_u <= 0:
        return -np.inf
    lp = beta_dist.logpdf(kappa, 2, 10)
    lp += norm.logpdf(phi_pi, 1.5, 0.5)
    lp += -2*np.log(sigma_u) - 0.1/sigma_u**2  # IG(0.1, 2) log-density
    return lp

def dsge_kalman_loglik(theta, Y_obs):
    """Evaluate NK model log-likelihood via Kalman filter."""
    kappa, phi_pi, sigma_u = theta
    beta_NK, sigma_IS = 0.99, 1.0
    phi_y = 0.5
    
    # Build model matrices
    G0 = np.array([[1.0, -kappa], [sigma_IS*phi_pi, 1+sigma_IS*phi_y]])
    G1 = np.array([[beta_NK, 0.0], [-sigma_IS, 1.0]])
    
    A = np.linalg.solve(G0, G1)
    eigs = np.linalg.eigvals(A)
    if not np.all(np.abs(eigs) > 1.0):
        return -np.inf  # indeterminate
    
    # MSV solution: [pi; x] = Omega * u_t, u_t = rho_u * u_{t-1} + eps_t
    rho_u = 0.7; Psi = np.linalg.solve(G0, np.array([1.0, 0.0]))
    Omega = np.linalg.solve(np.eye(2) - rho_u*A, Psi)
    
    # State-space: alpha_t = [pi_t, x_t, u_t]
    F = np.zeros((3,3)); F[:2,:2] = A; F[2,2] = rho_u
    F[:2, 2] = Psi
    H_obs = np.array([[1,0,0], [0,0,0]])  # only pi observed, x latent
    H_obs = np.eye(2)[:, :2] @ np.eye(2)  # simplify: observe pi only
    H_obs = np.array([[1,0,0]])  # observe pi only
    Q = np.zeros((3,3)); Q[2,2] = sigma_u**2
    R = np.array([[0.01]])  # small measurement error
    
    # Kalman filter
    T = len(Y_obs); m = F.shape[0]; p = H_obs.shape[0]
    a = np.zeros(m); P = np.eye(m) * 10.0  # diffuse initialization
    log_lik = 0.0
    
    for t in range(T):
        a_pred = F @ a
        P_pred = F @ P @ F.T + Q
        v = Y_obs[t:t+1] - H_obs @ a_pred
        Fv = H_obs @ P_pred @ H_obs.T + R
        K = P_pred @ H_obs.T @ np.linalg.inv(Fv)
        
        sign, ld = np.linalg.slogdet(Fv)
        log_lik -= 0.5*(ld + v @ np.linalg.solve(Fv, v) + p*np.log(2*np.pi))
        
        a = a_pred + K @ v
        P = (np.eye(m) - K @ H_obs) @ P_pred
    
    return float(log_lik)

def log_posterior(theta, Y_obs):
    lp = log_prior(theta)
    if not np.isfinite(lp): return -np.inf
    ll = dsge_kalman_loglik(theta, Y_obs)
    return lp + ll if np.isfinite(ll) else -np.inf

# Generate synthetic data
np.random.seed(42)
T_data = 200; kappa_t, phi_pi_t, sigma_u_t = 0.15, 1.5, 0.01
u = np.zeros(T_data)
for t in range(1,T_data): u[t] = 0.7*u[t-1] + sigma_u_t*np.random.randn()

# Simulate pi: approximate first-order solution
beta_NK, sigma_IS, phi_y = 0.99, 1.0, 0.5
G0_t = np.array([[1,-kappa_t],[sigma_IS*phi_pi_t,1+sigma_IS*phi_y]])
G1_t = np.array([[beta_NK,0],[-sigma_IS,1]]); A_t = np.linalg.inv(G0_t)@G1_t
Omega_t = np.linalg.solve(np.eye(2)-0.7*A_t, np.linalg.solve(G0_t,np.array([1,0])))
pi_data = Omega_t[0]*u + 0.003*np.random.randn(T_data)

# MH sampling
N_draws = 5000; burnin = 1000
theta_chain = np.zeros((N_draws, 3))
theta_chain[0] = [0.12, 1.6, 0.012]  # starting point
Sigma_prop = np.diag([0.01**2, 0.05**2, 0.002**2])  # proposal covariance
L_prop = np.linalg.cholesky(Sigma_prop)

lp_curr = log_posterior(theta_chain[0], pi_data)
n_accept = 0

for s in range(1, N_draws):
    proposal = theta_chain[s-1] + L_prop @ np.random.randn(3)
    lp_prop = log_posterior(proposal, pi_data)
    
    log_alpha = lp_prop - lp_curr
    if np.log(np.random.rand()) < log_alpha:
        theta_chain[s] = proposal
        lp_curr = lp_prop
        n_accept += 1
    else:
        theta_chain[s] = theta_chain[s-1]

accept_rate = n_accept / N_draws
posterior = theta_chain[burnin:]
print(f"Acceptance rate: {accept_rate*100:.1f}% (target: 20-30%)")
print(f"\nPosterior estimates (true values in brackets):")
params = ['kappa', 'phi_pi', 'sigma_u']
true_vals = [kappa_t, phi_pi_t, sigma_u_t]
for i, (p, tv) in enumerate(zip(params, true_vals)):
    pm = np.mean(posterior[:,i]); ps = np.std(posterior[:,i])
    ci = np.percentile(posterior[:,i], [5, 95])
    print(f"  {p}: mean={pm:.4f}, std={ps:.4f}, 90%CI=[{ci[0]:.4f},{ci[1]:.4f}]  (true={tv})")
```

---

## 30.9 Programming Exercises

### Exercise 30.1 (APL — Gelman–Rubin Statistic)

Run two independent MH chains from different starting points ($\theta_0^{(1)} = 0.5$ and $\theta_0^{(2)} = 3.0$) for the simple Gaussian posterior. Compute the Gelman–Rubin $\hat{R}$ statistic: (a) within-chain variance $W = (1/m)\sum_j S_j^2$ using `(+/S_j_sq)÷m`; (b) between-chain variance $B = n/(m-1)\sum_j(\bar\theta_j - \bar\theta)^2$; (c) $\hat{R} = \sqrt{\hat{V}/W}$. Verify $\hat{R} \to 1$ as the chain length increases and the two chains converge.

### Exercise 30.2 (Python — Adaptive MH)

Implement **adaptive MH** (Haario, Saksman, Tamminen, 2001): after a warmup of $s_0$ draws, update the proposal covariance as $\Sigma_q^{(s)} = (2.38)^2/k \cdot \text{Cov}(\theta^{(1)},\ldots,\theta^{(s)}) + \varepsilon I$ where $k$ is the dimension. (a) Show that this targets the optimal proposal covariance for Gaussian posteriors. (b) Compare acceptance rates between fixed and adaptive proposals for the 3-parameter NK model. (c) Verify the adaptive chain converges to the same posterior as the fixed chain.

### Exercise 30.3 (Julia — Marginal Likelihood via Laplace)

```julia
using LinearAlgebra, ForwardDiff, Optim

function laplace_marginal_lik(log_posterior, theta_mode)
    # Compute negative Hessian at mode
    H = ForwardDiff.hessian(theta -> -log_posterior(theta), theta_mode)
    k = length(theta_mode)
    
    # Laplace approximation
    lml = log_posterior(theta_mode) + (k/2)*log(2π) - 0.5*log(det(H))
    return lml
end

# Example: simple 2-parameter model
log_post(θ) = -0.5*((θ[1]-1.5)^2/0.25 + (θ[2]-0.15)^2/0.01)  # Gaussian posterior
θ_mode = [1.5, 0.15]
lml = laplace_marginal_lik(log_post, θ_mode)
println("Log marginal likelihood (Laplace): ", round(lml, digits=4))
println("Exact (Gaussian): ", -log(2π) - 0.5*log(0.25*0.01))
```

### Exercise 30.4 — Prior Sensitivity ($\star$)

For the 3-parameter NK model, investigate prior sensitivity. (a) Run the MH estimator with three different priors for $\kappa$: informative Beta(5,20) (mean 0.20), weakly informative Beta(2,10) (mean 0.17), and very flat Beta(1,1). (b) Report the posterior means and standard deviations for all three priors. (c) Compute the Bayes factor between the informative and flat prior specifications. (d) Interpret: is the data informative enough to overcome the prior, or does the posterior depend strongly on the prior specification?

---

## 30.10 Chapter Summary

**Key results:**

- The **Bayesian posterior** $p(\bm\theta|\mathbf{Y}) \propto \mathcal{L}(\bm\theta)p(\bm\theta)$ combines the Kalman filter likelihood with prior distributions over model parameters.
- **Standard priors** for DSGE: Beta for fractions ($\theta$, $\rho$); Gamma/Inverse-Gamma for positive parameters ($\sigma$, $\psi$, shock std devs); Normal (truncated) for Taylor rule coefficients.
- The **Metropolis–Hastings algorithm** generates posterior draws: propose $\bm\theta^* \sim q(\bm\theta|\bm\theta^{(s-1)})$; accept with probability $\min\{1, \pi(\bm\theta^*)/\pi(\bm\theta^{(s-1)})\}$. Invariance to the target distribution proved via detailed balance (Theorem 30.1).
- **Convergence diagnostics**: trace plots, effective sample size $N_{eff} = N/(1+2\sum\hat\rho_j)$, and Gelman–Rubin $\hat{R} < 1.1$.
- The **Laplace approximation** to the marginal likelihood: $\ln p(\mathbf{Y}|\mathcal{M}) \approx \ell(\hat{\bm\theta}^{MAP}) + \ln p(\hat{\bm\theta}^{MAP}) + \frac{k}{2}\ln(2\pi) - \frac{1}{2}\ln|\hat{H}|$.
- In APL: MH chain as a scan `{mh_step ⍵ scale}\ N_draws ⍴ theta0`; acceptance check via `log_ratio > ⍟?0` (accept when $\log U < \log\alpha$).

*Next: Chapter 31 — Implementing DSGE Models: Workflows in Dynare, APL, Python, Julia, and R*
