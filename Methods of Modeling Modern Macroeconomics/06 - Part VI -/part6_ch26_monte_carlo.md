# Chapter 26: Monte Carlo Methods

*Simulating Macroeconomic Models Under Uncertainty*

> *"Monte Carlo is the great equalizer: it makes hard problems easy, at the cost of making easy problems slow."*

**Cross-reference:** *Principles* Ch. 27 (RBC model simulation, second moments, HP filter); Ch. 39 (future methods: particle filters, ML in macro) **[P:Ch.27, P:Ch.39]**

---

## 26.1 From Theory to Simulation

A macroeconomic model is a data-generating process: given parameters and shock sequences, it produces sequences of output, inflation, employment, and other observables. **Simulation** is the process of running the model forward in time, drawing shocks from their specified distributions, and computing the resulting model paths.

Simulation serves four purposes in quantitative macroeconomics:

1. **Second-moment evaluation:** Compute model-implied standard deviations, autocorrelations, and cross-correlations to compare against data (Chapter 17).
2. **Impulse response generation:** Average the response to a single shock across many simulations to estimate the IRF with uncertainty bands.
3. **Policy evaluation:** Simulate the model under different policy rules to compare welfare.
4. **Likelihood approximation:** For nonlinear models where the Kalman filter is unavailable, the **particle filter** approximates the likelihood using simulated particles.

This chapter develops the simulation pipeline from the ground up: random number generation, the Tauchen discretization (formalizing the derivation sketched in Chapter 5), the full RBC simulation workflow, bootstrap confidence bands, and the particle filter.

---

## 26.2 Pseudorandom Number Generation

**Definition 26.1 (Pseudorandom Number Generator).** A **pseudorandom number generator (PRNG)** is a deterministic algorithm that produces a sequence $\{u_t\}$ that is statistically indistinguishable from i.i.d. Uniform$[0,1]$ to standard tests.

The standard PRNG in modern scientific computing is the **Mersenne Twister** (period $2^{19937}-1$, passes all standard statistical tests). All languages use it by default: `numpy.random`, Julia's `MersenneTwister()`, R's `set.seed()`.

**Inverse CDF method:** Any distribution with invertible CDF $F$ can be sampled via $X = F^{-1}(U)$ where $U \sim \text{Uniform}[0,1]$:

- Exponential($\lambda$): $X = -\ln(1-U)/\lambda$
- Geometric($p$): $X = \lceil\ln(U)/\ln(1-p)\rceil$

**Box–Muller transform:** Two independent $U_1, U_2 \sim \text{Uniform}[0,1]$ generate two standard normals:

$$Z_1 = \sqrt{-2\ln U_1}\cos(2\pi U_2), \qquad Z_2 = \sqrt{-2\ln U_1}\sin(2\pi U_2).$$

*Proof.* The joint density of $(Z_1, Z_2)$ factors as $\phi(z_1)\phi(z_2)$, confirming independence and standard normality. The transformation uses the fact that $-2\ln U_1 \sim \chi^2(2)$ and $2\pi U_2 \sim \text{Uniform}[0, 2\pi]$ independently. $\square$

In practice, `np.random.randn()`, `randn(n)` in Julia, and `rnorm(n)` in R implement efficient normal sampling algorithms (typically Ziggurat method, faster than Box–Muller).

---

## 26.3 The Tauchen (1986) Discretization: Full Derivation

The Tauchen method approximates a continuous AR(1) with a discrete Markov chain. We derive it formally here, connecting to the sketch in Chapter 5.

**Setup:** $z_{t+1} = \rho z_t + \varepsilon_{t+1}$, $\varepsilon_t \sim \mathcal{N}(0, \sigma_\varepsilon^2)$. We discretize to an $N$-point grid $\{z_1 < z_2 < \cdots < z_N\}$.

**Definition 26.2 (Tauchen Grid).** Set:
- $z_N = m\sigma_z$ and $z_1 = -z_N$ where $\sigma_z = \sigma_\varepsilon/\sqrt{1-\rho^2}$ (unconditional std dev) and $m = 3$ (covers $\pm3$ standard deviations).
- Equal spacing: $\Delta z = (z_N - z_1)/(N-1)$, so $z_j = z_1 + (j-1)\Delta z$.

**Definition 26.3 (Tauchen Transition Probabilities).** The transition probability from state $i$ to state $j$:

$$\Pi_{ij} = P(z_{t+1} = z_j | z_t = z_i) \approx P\!\left(\frac{z_j - \Delta z/2 - \rho z_i}{\sigma_\varepsilon} < Z < \frac{z_j + \Delta z/2 - \rho z_i}{\sigma_\varepsilon}\right),$$

where $Z \sim \mathcal{N}(0,1)$. With $\Phi$ the standard normal CDF:

$$\Pi_{ij} = \begin{cases} \Phi\!\left(\frac{z_1 + \Delta z/2 - \rho z_i}{\sigma_\varepsilon}\right) & j = 1 \text{ (lower boundary)} \\ \Phi\!\left(\frac{z_N - \Delta z/2 - \rho z_i}{\sigma_\varepsilon}\right) & j = N \text{ (upper boundary)} \\ \Phi\!\left(\frac{z_j + \Delta z/2 - \rho z_i}{\sigma_\varepsilon}\right) - \Phi\!\left(\frac{z_j - \Delta z/2 - \rho z_i}{\sigma_\varepsilon}\right) & \text{otherwise}\end{cases}$$

**Theorem 26.1 (Tauchen Approximation Quality).** As $N \to \infty$ with fixed $m$, the discrete Markov chain converges in distribution to the continuous AR(1). The mean and variance of the discrete chain converge to $0$ and $\sigma_z^2$ respectively. The autocorrelation converges to $\rho$.

*Proof sketch.* The $N$-point chain has mean $\bar{z} = \sum_j z_j\pi_j$ where $\bm\pi$ is the stationary distribution. As $N\to\infty$, the grid covers the support of $\mathcal{N}(0,\sigma_z^2)$ and the discrete probabilities converge to the normal probabilities, so $\bar{z}\to0$ and $\text{Var}(z)\to\sigma_z^2$. Autocorrelation $\to\rho$ follows from the transition probability construction. $\square$

In APL, the Tauchen matrix is built via the `∘.-` outer product — the computational core:

```apl
⍝ APL — Tauchen (1986) discretization
⎕IO←0 ⋄ ⎕ML←1

tauchen ← {rho sig_eps N m ← ⍵
    sig_z ← sig_eps ÷ (1-rho*2)*0.5         ⍝ unconditional std dev
    z_max ← m × sig_z  ⋄  z_min ← -z_max
    z     ← z_min + (z_max-z_min) × (⍳N) ÷ N-1   ⍝ uniform grid
    dz    ← (z_max-z_min) ÷ N-1             ⍝ grid spacing

    ⍝ Conditional means: μ_{ij} = (z_j ± dz/2 - ρ*z_i) / σ_ε
    ⍝ Outer product: diff[i,j] = z_j - ρ*z_i  (N×N)
    diff  ← (⍉z) ∘.- rho × z               ⍝ N×N difference matrix

    ⍝ Normal CDF (approximation; use ⎕PY for scipy.stats.norm.cdf in production)
    Phi ← {0.5 × 1 + 2○⍵ ÷ 2*0.5}         ⍝ erf-based approximation

    ⍝ Transition probabilities: Φ((z_j+dz/2-ρz_i)/σ) - Φ((z_j-dz/2-ρz_i)/σ)
    P ← (Phi (diff + dz÷2) ÷ sig_eps) - (Phi (diff - dz÷2) ÷ sig_eps)

    ⍝ Boundary corrections (absorb tail probability)
    P[;0]   ← Phi (z[0] + dz÷2 - rho × z) ÷ sig_eps
    P[;N-1] ← 1 - +⌿ P[;⍳N-1]

    ⍝ Normalize rows (should already sum to 1 up to floating point)
    P ÷ +/P    ⍝ row-normalize
    z P}       ⍝ return grid and transition matrix

⍝ Test: TFP process for RBC model
rho_A←0.95  ⋄  sig_A←0.0072  ⋄  N←7
z_grid Pi_A ← tauchen rho_A sig_A N 3

z_grid    ⍝ 7-point grid
+/Pi_A    ⍝ row sums: all should be 1
```

---

## 26.4 Simulating Model Economies

**Algorithm 26.1 (Model Simulation Pipeline).**

Input: Policy functions $c^*(k, A)$, $n^*(k, A)$ (from VFI or log-linearization), TFP transition matrix $\Pi$, initial condition $(k_0, A_0)$, simulation length $T$.

1. **Draw shocks:** Simulate the Markov chain $\{A_t\}_{t=0}^T$ using $\Pi$.
2. **Compute decisions:** At each $t$, evaluate $c_t = c^*(k_t, A_t)$.
3. **Update state:** $k_{t+1} = (1-\delta)k_t + I_t$ where $I_t = Y_t - C_t$.
4. **Aggregate:** Collect $(Y_t, C_t, I_t, n_t, k_t)$ for $t = 0, \ldots, T$.
5. **Compute moments:** Apply HP filter; compute standard deviations, autocorrelations, cross-correlations.

**Markov chain simulation:** Given the transition matrix $\Pi$ and initial state $j_0$, simulate a length-$T$ chain:

```apl
⍝ APL — Markov chain simulation
⎕IO←0 ⋄ ⎕ML←1

⍝ Cumulative distribution rows of Pi (for inverse CDF sampling)
Pi_cdf ← +\⍤1 ⊢ Pi_A           ⍝ row-wise cumulative sum

⍝ One step: given current state idx, draw next state
markov_step ← {idx ← ⍵
    u ← ?0                      ⍝ uniform draw in [0,1)
    +/ u > Pi_cdf[idx;]}        ⍝ count how many CDF values u exceeds

⍝ Simulate T steps from state 3 (middle of grid)
T ← 1000  ⋄  j0 ← 3
chain ← markov_step \ T ⍴ j0   ⍝ scan: apply markov_step T times
z_sim ← z_grid[chain]           ⍝ realized z values
A_sim ← *z_sim                  ⍝ TFP levels (exp of log-TFP)
```

---

## 26.5 Bootstrap Methods

The **bootstrap** generates confidence intervals and standard errors for statistics that have no analytical distribution theory.

**Algorithm 26.2 (Residual Bootstrap for VAR IRFs).**

1. Estimate the VAR to obtain $\hat{A}$, residuals $\{\hat{\mathbf{e}}_t\}$, and IRF estimate $\widehat{\text{IRF}}(h)$.
2. For $b = 1, \ldots, B$:
   a. Draw $T$ residuals with replacement: $\mathbf{e}_t^{(b)} \sim \{\hat{\mathbf{e}}_1, \ldots, \hat{\mathbf{e}}_T\}$.
   b. Reconstruct bootstrap data: $\mathbf{y}_t^{(b)} = \hat{A}\mathbf{y}_{t-1}^{(b)} + \mathbf{e}_t^{(b)}$.
   c. Re-estimate the VAR and compute $\widehat{\text{IRF}}^{(b)}(h)$.
3. Report the 16th and 84th percentiles of $\{\widehat{\text{IRF}}^{(b)}(h)\}$ as the 68% confidence band.

**Theorem 26.2 (Bootstrap Consistency).** Under weak regularity conditions (stationarity, ergodicity, finite moments), the bootstrap distribution of $\sqrt{T}(\widehat{\text{IRF}}^{(b)} - \widehat{\text{IRF}})$ consistently estimates the distribution of $\sqrt{T}(\widehat{\text{IRF}} - \text{IRF})$. That is, bootstrap confidence intervals have asymptotically correct coverage.

*Proof sketch.* The bootstrap resamples residuals that converge (by ergodicity) to the true i.i.d. residual distribution. By the functional CLT, the bootstrap sample path distribution converges to the same Gaussian limit as the true distribution. See Gonçalves and Kilian (2004) for the full proof in the VAR context. $\square$

---

## 26.6 The Particle Filter for Nonlinear Models

When the state-space model is **nonlinear** — as in the DSGE model with occasionally-binding constraints, or the Markov-switching model — the Kalman filter (Chapter 20) is no longer optimal. The **particle filter** (Sequential Monte Carlo) approximates the distribution $p(\bm\alpha_t | \mathbf{y}_1, \ldots, \mathbf{y}_t)$ using a cloud of weighted particles.

**Definition 26.4 (Particle Filter).** A **particle filter** represents the distribution $p_t(\bm\alpha_t | \mathbf{Y}_{1:t})$ by $N$ weighted particles $\{(\bm\alpha_t^{(i)}, w_t^{(i)})\}_{i=1}^N$ where $\sum_i w_t^{(i)} = 1$.

**Algorithm 26.3 (Bootstrap Particle Filter / Sequential Importance Resampling).**

Initialize: Draw $\bm\alpha_0^{(i)} \sim p(\bm\alpha_0)$, set $w_0^{(i)} = 1/N$ for all $i$.

For $t = 1, \ldots, T$:

**Step 1 — Propagate:** Draw $\bm\alpha_t^{(i)} \sim p(\bm\alpha_t | \bm\alpha_{t-1}^{(i)})$ (transition density).

**Step 2 — Weight:** Compute importance weights:
$$\tilde{w}_t^{(i)} = w_{t-1}^{(i)} \cdot p(\mathbf{y}_t | \bm\alpha_t^{(i)}) \quad (\text{likelihood of observation given particle}),$$
then normalize: $w_t^{(i)} = \tilde{w}_t^{(i)} / \sum_j\tilde{w}_t^{(j)}$.

**Step 3 — Resample:** When the effective sample size $N_\text{eff} = 1/\sum_i(w_t^{(i)})^2 < N/2$, resample $N$ particles from $\{(\bm\alpha_t^{(i)}, w_t^{(i)})\}$ with replacement (systematic resampling); reset all weights to $1/N$.

**Log-likelihood contribution at $t$:**
$$\ell_t = \ln\left(\frac{1}{N}\sum_i\tilde{w}_t^{(i)}\right) \approx \ln p(\mathbf{y}_t | \mathbf{Y}_{1:t-1}).$$

The total log-likelihood $\mathcal{L} = \sum_t\ell_t$ approximates the Kalman filter log-likelihood for nonlinear models.

**Connecting to the Kalman filter:** For linear Gaussian models, the particle filter converges to the Kalman filter as $N\to\infty$: the particle distribution converges to the Gaussian filtering distribution. The Kalman filter is exact because the Gaussian is fully characterized by mean and variance; the particle filter is exact in the limit for any distribution.

```python
import numpy as np

def particle_filter(y, f_transition, f_likelihood, a0_sampler, N=1000):
    """Bootstrap particle filter for nonlinear state-space model.
    
    f_transition(alpha): sample next state given current state
    f_likelihood(y_obs, alpha): log p(y | alpha)
    a0_sampler(): sample initial state
    """
    T = len(y)
    # Initialize particles
    particles = np.array([a0_sampler() for _ in range(N)])
    weights = np.ones(N) / N
    log_lik = 0.0
    
    filtered_means = np.zeros(T)
    
    for t in range(T):
        # Propagate
        particles = np.array([f_transition(a) for a in particles])
        
        # Weight
        log_w = np.array([f_likelihood(y[t], a) for a in particles])
        log_w -= log_w.max()  # numerical stability
        w_unnorm = np.exp(log_w) * weights
        
        # Log-likelihood contribution
        log_lik += np.log(w_unnorm.mean())
        
        # Normalize
        weights = w_unnorm / w_unnorm.sum()
        
        # Filtered mean
        filtered_means[t] = np.dot(weights, particles)
        
        # Resample if needed
        N_eff = 1.0 / np.sum(weights**2)
        if N_eff < N/2:
            indices = np.random.choice(N, N, p=weights)
            particles = particles[indices]
            weights = np.ones(N) / N
    
    return log_lik, filtered_means

# Example: nonlinear state-space model
# State: alpha_t = 0.8*alpha_{t-1} + eta_t, eta ~ N(0, 1)
# Obs:   y_t = alpha_t^2/20 + eps_t, eps ~ N(0, 1)
np.random.seed(42); T = 100
alpha_true = np.zeros(T)
for t in range(1,T): alpha_true[t] = 0.8*alpha_true[t-1] + np.random.randn()
y_obs = alpha_true**2/20 + np.random.randn(T)

from scipy.stats import norm
f_trans = lambda a: 0.8*a + np.random.randn()
f_lik   = lambda y, a: norm.logpdf(y, a**2/20, 1.0)
a0_sampler = lambda: np.random.randn()

log_lik, filtered = particle_filter(y_obs, f_trans, f_lik, a0_sampler, N=2000)
print(f"Particle filter log-likelihood: {log_lik:.2f}")

import matplotlib.pyplot as plt
plt.figure(figsize=(10,4))
plt.plot(alpha_true, 'b-', label='True state', alpha=0.7)
plt.plot(filtered, 'r-', label='PF filtered mean')
plt.legend(); plt.title('Particle Filter: Nonlinear State Estimation'); plt.show()
```

---

## 26.7 Worked Example: RBC Second Moments via Simulation

*Cross-reference: Principles Ch. 27.3 (calibration and second moments)* **[P:Ch.27.3]**

Using the log-linearized RBC model solution from Chapter 28 (previewed here), simulate 10,000 periods and compute HP-filtered second moments.

```python
import numpy as np
from scipy.stats import norm

# RBC model parameters and Tauchen discretization (from Ch.17/Ch.5)
alpha, delta, beta, sigma_crra = 0.36, 0.025, 0.99, 1.0
rho_A, sig_A = 0.95, 0.0072

# Tauchen discretization: 7-point grid for log(A)
N_A = 7; m = 3
sig_z = sig_A / np.sqrt(1-rho_A**2)
z_grid = np.linspace(-m*sig_z, m*sig_z, N_A)
A_grid = np.exp(z_grid)
dz = z_grid[1]-z_grid[0]

Pi = np.zeros((N_A, N_A))
for i in range(N_A):
    mu_z = rho_A*z_grid[i]
    Pi[i,0] = norm.cdf((z_grid[0]+dz/2-mu_z)/sig_A)
    Pi[i,-1] = 1-norm.cdf((z_grid[-1]-dz/2-mu_z)/sig_A)
    for j in range(1, N_A-1):
        Pi[i,j] = norm.cdf((z_grid[j]+dz/2-mu_z)/sig_A) - norm.cdf((z_grid[j]-dz/2-mu_z)/sig_A)
Pi = Pi / Pi.sum(axis=1, keepdims=True)  # normalize

# Simulate Markov chain
np.random.seed(42); T_sim = 12000
Pi_cdf = np.cumsum(Pi, axis=1)
A_idx = [N_A//2]  # start at middle state
for t in range(T_sim-1):
    u = np.random.rand()
    A_idx.append(int(np.searchsorted(Pi_cdf[A_idx[-1]], u)))
A_sim = A_grid[A_idx]

# Log-linearized RBC decisions (approximate; see Ch.28 for exact computation)
# Approximation: output ~ A * k_ss^alpha * n_ss^(1-alpha) * (1 + hat_a)
# where hat_a = log(A/A_ss) and hat_y ~ (1/(1-alpha*beta*rho_A)) * hat_a
kstar = (alpha/(1/beta-1+delta))**(1/(1-alpha))
Yss = kstar**alpha; Css = Yss - delta*kstar; Iss = delta*kstar

psi_y = 1/(1-alpha*beta*rho_A)  # output multiplier w.r.t. TFP (approximate)
psi_c = psi_y*(1-delta*alpha/(1/beta-1+delta))
psi_i = (psi_y - psi_c*Css/Yss)/Iss*Yss  # rough investment multiplier

hat_a = np.log(A_sim)  # log-deviation of TFP
Y_sim = Yss * np.exp(psi_y * hat_a)
C_sim = Css * np.exp(psi_c * hat_a)
I_sim = Iss * np.exp(psi_i * hat_a)

# HP filter
def hp_filter(x, lam=1600):
    T = len(x)
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve
    e = np.ones(T)
    D2 = diags([e[:-2], -2*e[:-1], e], [0,1,2], shape=(T-2,T))
    I = diags(e)
    tau = spsolve((I + lam*D2.T@D2).toarray(), x)
    return x - tau

burnin = 2000
log_Y = np.log(Y_sim[burnin:]); log_C = np.log(C_sim[burnin:]); log_I = np.log(I_sim[burnin:])
cyc_Y = hp_filter(log_Y); cyc_C = hp_filter(log_C); cyc_I = hp_filter(log_I)

def stats(x, y=None):
    std_x = np.std(x)*100
    if y is None: return std_x, np.corrcoef(x[:-1],x[1:])[0,1]
    return std_x, np.corrcoef(x,y)[0,1]

std_Y, ac_Y = stats(cyc_Y)
std_C, cor_CY = stats(cyc_C, cyc_Y)
std_I, cor_IY = stats(cyc_I, cyc_Y)

print(f"\nRBC Model Second Moments (HP-filtered, λ=1600):")
print(f"{'Variable':<12} {'Std(%)':>8} {'Rel.Std':>8} {'Corr(Y)':>8}")
print(f"{'Output Y':<12} {std_Y:>8.3f} {1.000:>8.3f} {1.000:>8.3f}")
print(f"{'Consump. C':<12} {std_C:>8.3f} {std_C/std_Y:>8.3f} {cor_CY:>8.3f}")
print(f"{'Invest. I':<12} {std_I:>8.3f} {std_I/std_Y:>8.3f} {cor_IY:>8.3f}")
print(f"\nData targets:  Std(C)/Std(Y)≈0.50, Std(I)/Std(Y)≈3.08")
print(f"               Corr(C,Y)≈0.88, Corr(I,Y)≈0.93")
```

---

## 26.8 Programming Exercises

### Exercise 26.1 (APL — Tauchen Discretization)

Implement the complete Tauchen algorithm in APL as described in Section 26.3. (a) The outer product `diff ← (⍉z) ∘.- rho × z` generates the $N\times N$ matrix of conditional means in one expression. (b) Using `⎕PY` to call `scipy.stats.norm.cdf` for the normal CDF, compute the full transition matrix. (c) Verify: row sums equal 1; the stationary distribution (eigenvector of $\Pi'$) matches $\mathcal{N}(0, \sigma_z^2)$ approximately.

### Exercise 26.2 (Python — Rouwenhorst Method)

The **Rouwenhorst (1995) method** is more accurate than Tauchen for highly persistent processes ($\rho > 0.9$). For $N=2$: $\Pi = \begin{pmatrix}p & 1-p \\ 1-q & q\end{pmatrix}$, $p = q = (1+\rho)/2$. For $N > 2$: recurse using $(N-1)$-point matrix. (a) Implement Rouwenhorst for $N \in \{3, 5, 7, 9\}$ and $\rho = 0.95$. (b) Compare the implied autocorrelation $\rho_1^{model}$ to the target $\rho = 0.95$ for both Tauchen and Rouwenhorst. (c) Show Rouwenhorst matches $\rho$ exactly by construction.

### Exercise 26.3 (Julia — Bootstrap IRF Bands)

```julia
using LinearAlgebra, Statistics

function bootstrap_irf(Y, p, H, B=500, identification=:cholesky)
    T, n = size(Y)
    # OLS VAR(p)
    X = hcat([Y[p-j+1:T-j,:] for j in 1:p]...)
    Ytgt = Y[p+1:end,:]
    B_hat = (X'X) \ (X'Ytgt)  # np×n
    resid = Ytgt - X*B_hat
    Sigma = resid'*resid/(T-p-n*p)
    
    A1 = B_hat[1:n,:]'  # companion (first lag only for simplicity)
    P = cholesky(Sigma).L  # Cholesky identification
    
    # Point IRF
    irf_point = zeros(H, n, n)
    for h in 1:H, j in 1:n
        shock = P[:,j]
        irf_point[h,:,j] = (A1^(h-1))*shock
    end
    
    # Bootstrap
    irf_boot = zeros(B, H, n, n)
    for b in 1:B
        # Resample residuals
        idx = rand(1:size(resid,1), size(resid,1))
        resid_b = resid[idx,:]
        # Reconstruct data
        Y_b = zeros(T, n); Y_b[1:p,:] = Y[1:p,:]
        for t in p+1:T
            Xb = vcat([Y_b[t-j,:] for j in 1:p]...)
            Y_b[t,:] = B_hat'*Xb + resid_b[t-p,:]
        end
        # Re-estimate
        Xb = hcat([Y_b[p-j+1:T-j,:] for j in 1:p]...)
        Ytb = Y_b[p+1:end,:]
        Bb = (Xb'Xb)\(Xb'Ytb)
        rb = Ytb - Xb*Bb
        Sb = rb'*rb/(T-p-n*p)
        Pb = cholesky(Sb).L
        Ab = Bb[1:n,:]'
        for h in 1:H, j in 1:n
            irf_boot[b,h,:,j] = (Ab^(h-1))*Pb[:,j]
        end
    end
    
    return irf_point, irf_boot
end

# Simulate a 2-variable VAR and compute bootstrap bands
Random.seed!(42); T=200; n=2
A_true = [0.8 0.1; 0.0 0.7]; L_true = [0.01 0; 0.005 0.008]
Y = zeros(T,n)
for t in 2:T; Y[t,:] = A_true*Y[t-1,:] + L_true*randn(2); end

irf, irf_b = bootstrap_irf(Y, 1, 20, 200)
# 68% confidence bands
lo = dropdims(mapslices(x->quantile(x,0.16), irf_b[:,1,:,1], dims=1), dims=1)
hi = dropdims(mapslices(x->quantile(x,0.84), irf_b[:,1,:,1], dims=1), dims=1)
println("Bootstrap 68% band width at h=1: $(round(hi[1]-lo[1], digits=5))")
```

### Exercise 26.4 — Particle Filter vs. Kalman ($\star$)

For the linear Gaussian state-space model: $\alpha_{t+1} = 0.9\alpha_t + \eta_t$, $\eta_t \sim \mathcal{N}(0,1)$; $y_t = \alpha_t + \varepsilon_t$, $\varepsilon_t \sim \mathcal{N}(0,0.5)$. (a) Run the Kalman filter (exactly) and the particle filter ($N = 100, 500, 1000, 5000$) on the same 200-period dataset. (b) Compare the filtered means and log-likelihoods. (c) Show that the particle filter converges to the Kalman filter as $N\to\infty$. (d) Compute the RMSE of the filtered state relative to the true state for each method and $N$.

---

## 26.9 Chapter Summary

**Key results:**

- **Pseudorandom generation:** Box–Muller transform generates normals from uniforms (proved via change-of-variables); Mersenne Twister is the standard PRNG.
- **Tauchen (1986) discretization** (Theorem 26.1): the $N$-point Markov chain with transition matrix $\Pi_{ij} = \Phi((z_j+\Delta z/2-\rho z_i)/\sigma_\varepsilon) - \Phi((z_j-\Delta z/2-\rho z_i)/\sigma_\varepsilon)$ converges to the continuous AR(1) as $N\to\infty$. In APL: the core is `diff ← (⍉z) ∘.- rho × z` (one outer product).
- **Simulation pipeline:** Draw Markov chain → compute decisions via policy functions → compute aggregate time series → HP-filter → compute second moments.
- **Bootstrap consistency** (Theorem 26.2): residual bootstrap confidence bands for VAR IRFs have asymptotically correct coverage under stationarity and ergodicity.
- **Particle filter:** Represents $p(\alpha_t|Y_{1:t})$ by $N$ weighted particles $\{(\alpha_t^{(i)}, w_t^{(i)})\}$; propagate → weight → resample. Log-likelihood $\sum_t\ln(\frac{1}{N}\sum_i\tilde{w}_t^{(i)})$ approximates the Kalman likelihood for nonlinear models.
- In APL: Markov chain simulation uses cumulative row sums and `⍸` (interval index); the HP filter uses `⌹` on the banded penalty matrix; all second moments are computed via `+.×` and standard deviation operators.

---

*Next: Part VII — Computational Methods for DSGE Models*
