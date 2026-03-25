# Chapter 5: Stochastic Processes for Aggregate Shocks

*AR(1) Processes and Macroeconomic Uncertainty*

> *"The theory of stochastic processes is the theory of time-indexed families of random variables. In macroeconomics, every variable of interest is one."*

**Cross-reference:** *Principles* Ch. 6 (time-series properties, HP filter, unit roots, data vintages); Ch. 15 (rational expectations, information sets); Ch. 27 (RBC calibration — TFP as AR(1)) **[P:Ch.6, P:Ch.15, P:Ch.27]**

---

## 5.1 Why Stochastic Models? The Centrality of Uncertainty

The deterministic dynamic models of Chapters 3 and 4 assume that once initial conditions and parameters are specified, the entire future path of the economy is known. This is an obviously wrong description of the world. GDP fluctuates for reasons that are not fully predictable in advance. Monetary policy decisions are made under uncertainty about the current state of the economy. Firms' investment decisions depend on uncertain future profitability. Consumer spending responds to news about future income.

Modern macroeconomics — from the rational expectations revolution onward — treats uncertainty not as an afterthought but as a fundamental feature of the economic environment [P:Ch.15]. The stochastic process is to modern macroeconomics what the differential equation is to classical physics: the natural language for describing how systems evolve over time when driven by unpredictable forces.

This chapter introduces the key concepts of probability and stochastic process theory that underlie the macroeconomic models of Parts V–IX. We develop: the formal probabilistic framework (probability spaces, filtrations, conditional expectations); stationary stochastic processes and their characterization by autocovariance functions and spectra; the AR(1) process in depth; and the Wold decomposition, which guarantees that every covariance-stationary process can be represented as a moving average of white noise — the fundamental result underlying all time-series econometrics.

---

## 5.2 Probability Spaces and Information Sets

**Definition 5.1 (Probability Space).** A **probability space** is a triple $(\Omega, \mathcal{F}, P)$ where:
- $\Omega$ is the **sample space** — the set of all possible states of the world.
- $\mathcal{F}$ is a **sigma-algebra** on $\Omega$ — a collection of subsets of $\Omega$ (events) closed under countable unions, countable intersections, and complementation.
- $P: \mathcal{F} \to [0,1]$ is a **probability measure** — $P(\Omega) = 1$ and $P$ is countably additive.

A **random variable** $X: \Omega \to \mathbb{R}$ is a measurable function from $\Omega$ to $\mathbb{R}$.

For macroeconomics, the key interpretation is: $\Omega$ is the space of all possible histories of the economy (sequences of TFP realizations, policy shocks, preference shocks); an event $A \in \mathcal{F}$ is a statement about the economy that is either true or false in any realized history; and $P(A)$ is the prior probability of that event.

### 5.2.1 Filtrations and Information

In dynamic models, information accumulates over time. At date $t$, agents have observed the history of the economy up to $t$ but not the future.

**Definition 5.2 (Filtration).** A **filtration** $\{\mathcal{F}_t\}_{t=0}^\infty$ is a non-decreasing sequence of sigma-algebras: $\mathcal{F}_s \subseteq \mathcal{F}_t$ for all $s \leq t$. The sigma-algebra $\mathcal{F}_t$ represents the information available at date $t$.

A stochastic process $\{X_t\}$ is **adapted** to the filtration $\{\mathcal{F}_t\}$ if $X_t$ is $\mathcal{F}_t$-measurable for each $t$ — that is, the value of $X_t$ is known at date $t$. In macroeconomics, every variable we observe (GDP, inflation, interest rates) is adapted to the information filtration of agents.

**Definition 5.3 (Conditional Expectation).** The **conditional expectation** $\mathbb{E}[X|\mathcal{F}_t]$ is the $\mathcal{F}_t$-measurable random variable that best predicts $X$ given information $\mathcal{F}_t$ in the mean-squared sense. It satisfies:
$$\mathbb{E}[\mathbb{E}[X|\mathcal{F}_t]|\mathcal{F}_s] = \mathbb{E}[X|\mathcal{F}_s] \quad \text{for all } s \leq t.$$

This is the **law of iterated expectations** — a key tool in rational expectations models. It says: the best forecast of the future best forecast of $X$ is the current best forecast of $X$.

**The rational expectations hypothesis** [P:Ch.15] states that agents' subjective expectations equal the conditional expectation given their information set: $\mathbb{E}_t^{agents}[X_{t+h}] = \mathbb{E}[X_{t+h}|\mathcal{F}_t]$. Agents use the correct probability model; they are not systematically wrong.

---

## 5.3 Stationarity, Autocovariance, and the Autocorrelation Function

### 5.3.1 Stationarity

**Definition 5.4 (Covariance Stationarity).** A stochastic process $\{X_t\}$ is **covariance-stationary** (or weakly stationary) if:
1. $\mathbb{E}[X_t] = \mu$ for all $t$ (constant mean).
2. $\text{Var}(X_t) = \sigma^2 < \infty$ for all $t$ (finite, constant variance).
3. $\text{Cov}(X_t, X_{t-k}) = \gamma_k$ depends only on the lag $k$, not on $t$ (time-invariant autocovariances).

**Definition 5.5 (Strict Stationarity).** A process is **strictly stationary** if the joint distribution of $(X_{t_1}, \ldots, X_{t_n})$ equals the joint distribution of $(X_{t_1+h}, \ldots, X_{t_n+h})$ for any $h$ — the distribution is invariant to time translation. Strict stationarity plus finite second moments implies covariance stationarity.

Most empirical time-series methods require only covariance stationarity, which is weaker and easier to verify.

### 5.3.2 The Autocovariance and Autocorrelation Functions

**Definition 5.6 (Autocovariance Function).** For a covariance-stationary process $\{X_t\}$ with mean $\mu$:
$$\gamma_k = \text{Cov}(X_t, X_{t-k}) = \mathbb{E}[(X_t - \mu)(X_{t-k} - \mu)], \quad k = 0, 1, 2, \ldots$$

Note $\gamma_0 = \text{Var}(X_t) = \sigma^2$ and $\gamma_k = \gamma_{-k}$ (symmetry).

**Definition 5.7 (Autocorrelation Function, ACF).** The **autocorrelation function** is:
$$\rho_k = \frac{\gamma_k}{\gamma_0} \in [-1, 1].$$

The ACF $\{\rho_k\}_{k=0}^\infty$ characterizes the temporal persistence of the process. For macroeconomic variables:
- U.S. log real GDP (HP-filtered) has $\rho_1 \approx 0.87$, $\rho_4 \approx 0.52$ — high persistence.
- U.S. TFP growth (quarterly) has $\rho_1 \approx 0.70$, falling off gradually.
- White noise has $\rho_k = 0$ for all $k \geq 1$.

**Bochner's theorem:** A sequence $\{\rho_k\}$ is a valid ACF if and only if it is positive semi-definite. This means the Fourier transform of the ACF (the spectral density, defined below) must be non-negative everywhere.

---

## 5.4 The AR(1) Process: The Workhorse of Macroeconomic Dynamics

The **autoregressive process of order 1**, or AR(1), is the single most important stochastic process in macroeconomics. It appears as:
- The technology shock in every RBC model [P:Ch.27]: $\ln A_t = \rho\ln A_{t-1} + \varepsilon_t$.
- The monetary policy shock in NK models [P:Ch.23]: $\hat{m}_t = \rho_m\hat{m}_{t-1} + \varepsilon_t^m$.
- The approximation to any persistent but stationary shock process.

**Definition 5.8 (AR(1) Process).** An AR(1) process satisfies:
$$X_t = \mu + \rho(X_{t-1} - \mu) + \varepsilon_t, \quad \varepsilon_t \sim \text{WN}(0, \sigma_\varepsilon^2),$$

where $\rho \in (-1, 1)$ is the **persistence parameter**, $\mu$ is the unconditional mean, and $\varepsilon_t$ is **white noise** — serially uncorrelated innovations with zero mean and constant variance.

**Definition 5.9 (White Noise).** $\{\varepsilon_t\}$ is **white noise**, written $\varepsilon_t \sim \text{WN}(0, \sigma^2)$, if:
- $\mathbb{E}[\varepsilon_t] = 0$ for all $t$.
- $\mathbb{E}[\varepsilon_t^2] = \sigma^2$ for all $t$.
- $\mathbb{E}[\varepsilon_t\varepsilon_s] = 0$ for all $t \neq s$.

If additionally $\varepsilon_t \sim \mathcal{N}(0, \sigma^2)$ i.i.d., we have **Gaussian white noise** — the standard assumption in DSGE models.

### 5.4.1 Moments of the AR(1)

Taking expectations of the AR(1) definition (assuming stationarity):
$$\mathbb{E}[X_t] = \mu + \rho(\mathbb{E}[X_t] - \mu) + 0 \implies \mathbb{E}[X_t] = \mu.$$

For the variance, use the WLOG assumption $\mu = 0$: $X_t = \rho X_{t-1} + \varepsilon_t$. Then:
$$\text{Var}(X_t) = \rho^2\text{Var}(X_{t-1}) + \sigma_\varepsilon^2 \implies \sigma_X^2 = \frac{\sigma_\varepsilon^2}{1 - \rho^2}.$$

This requires $|\rho| < 1$ for the variance to be finite. The autocovariances:

**Theorem 5.1 (AR(1) Autocovariance Function).** For the zero-mean AR(1) $X_t = \rho X_{t-1} + \varepsilon_t$ with $|\rho| < 1$:
$$\gamma_k = \frac{\sigma_\varepsilon^2}{1-\rho^2}\rho^{|k|} = \sigma_X^2\rho^{|k|}, \quad \rho_k = \rho^{|k|}.$$

*Proof.* Multiply both sides of $X_t = \rho X_{t-1} + \varepsilon_t$ by $X_{t-k}$ and take expectations:
$$\mathbb{E}[X_tX_{t-k}] = \rho\mathbb{E}[X_{t-1}X_{t-k}] + \mathbb{E}[\varepsilon_tX_{t-k}].$$

For $k \geq 1$, $\varepsilon_t$ is uncorrelated with $X_{t-k}$ (since $t > t-k$), so $\mathbb{E}[\varepsilon_tX_{t-k}] = 0$. Therefore $\gamma_k = \rho\gamma_{k-1}$ for $k \geq 1$, with $\gamma_0 = \sigma_\varepsilon^2/(1-\rho^2)$. Solving the recursion: $\gamma_k = \rho^k\gamma_0 = \sigma_X^2\rho^k$. $\square$

The ACF $\rho_k = \rho^k$ decays geometrically. **High $|\rho|$** (near 1) means slow decay — the process has long memory and is highly persistent. **Low $|\rho|$** (near 0) means rapid decay — shocks dissipate quickly.

### 5.4.2 Moving Average Representation

The AR(1) can be written as an infinite moving average (MA):
$$X_t = \varepsilon_t + \rho\varepsilon_{t-1} + \rho^2\varepsilon_{t-2} + \cdots = \sum_{j=0}^\infty \rho^j\varepsilon_{t-j}.$$

*Proof.* Iterate $X_t = \rho X_{t-1} + \varepsilon_t$ backward:
$$X_t = \varepsilon_t + \rho(\rho X_{t-2} + \varepsilon_{t-1}) + \cdots = \sum_{j=0}^T\rho^j\varepsilon_{t-j} + \rho^{T+1}X_{t-T-1} \to \sum_{j=0}^\infty\rho^j\varepsilon_{t-j}$$
since $|\rho^{T+1}X_{t-T-1}| \to 0$ in mean square when $|\rho| < 1$. $\square$

The MA representation is the **impulse response function** of the AR(1): the coefficient on $\varepsilon_{t-j}$ gives the effect of a shock $j$ periods ago on $X_t$ today. For the AR(1), the IRF is $\psi_j = \rho^j$ — a geometrically declining function of the lag $j$.

In DSGE models, the IRF of any variable to any shock is a sum of such geometrically declining terms, weighted by the eigenvalues of the transition matrix. Chapter 17 computes these IRFs for the RBC model; Chapter 28 does so for the full NK model.

---

## 5.5 ARMA Processes and the Wold Decomposition

### 5.5.1 ARMA(p,q) Processes

**Definition 5.10 (ARMA Process).** An ARMA($p$,$q$) process satisfies:
$$X_t = \mu + \sum_{i=1}^p\phi_i(X_{t-i}-\mu) + \sum_{j=0}^q\theta_j\varepsilon_{t-j}, \quad \varepsilon_t \sim \text{WN}(0, \sigma^2),$$

where the AR polynomial $\phi(L) = 1 - \phi_1 L - \cdots - \phi_p L^p$ and MA polynomial $\theta(L) = 1 + \theta_1 L + \cdots + \theta_q L^q$ are defined in terms of the **lag operator** $L$: $LX_t = X_{t-1}$.

The process is **stationary** iff all roots of $\phi(z) = 0$ lie outside the unit circle ($|z| > 1$). It is **invertible** iff all roots of $\theta(z) = 0$ lie outside the unit circle.

### 5.5.2 The Wold Decomposition

The Wold theorem is arguably the most important theorem in time-series analysis. It says that any covariance-stationary process can be decomposed into a predictable component and an unpredictable MA($\infty$) component.

**Theorem 5.2 (Wold Decomposition).** Let $\{X_t\}$ be a zero-mean covariance-stationary process. Then:
$$X_t = \sum_{j=0}^\infty \psi_j\varepsilon_{t-j} + V_t,$$

where:
1. $\varepsilon_t = X_t - \mathbb{E}[X_t|\mathcal{F}_{t-1}]$ is the **one-step-ahead forecast error** (innovation), which is white noise with $\sigma^2 = \gamma_0 - \text{var(best linear predictor)}$.
2. $\psi_0 = 1$ and $\sum_{j=0}^\infty\psi_j^2 < \infty$.
3. $V_t$ is the **deterministic component** — perfectly predictable from its own past.
4. $\varepsilon_t \perp V_s$ for all $t, s$ (the components are uncorrelated).

For purely non-deterministic processes (those without a deterministic component), $V_t = 0$ and the entire process is its MA($\infty$) representation. The **impulse response coefficients** $\psi_j$ are the dynamic multipliers: a unit innovation at $t$ affects $X_{t+j}$ by $\psi_j$.

*Economic significance.* The Wold decomposition underpins the **structural VAR** methodology [P:Appendix B]. By identifying which innovations $\varepsilon_t$ correspond to which structural shocks (monetary policy, technology, demand), empirical macroeconomists can trace the dynamic effects of each shock through the economy. Chapter 19 develops this methodology in full.

---

## 5.6 VAR Processes and the Companion Form

A **vector autoregression** (VAR) of order $p$ generalizes the scalar AR($p$) to a vector of variables:
$$\mathbf{X}_t = \mathbf{c} + \Phi_1\mathbf{X}_{t-1} + \Phi_2\mathbf{X}_{t-2} + \cdots + \Phi_p\mathbf{X}_{t-p} + \bm{\varepsilon}_t,$$

where $\mathbf{X}_t \in \mathbb{R}^n$, $\Phi_i \in \mathbb{R}^{n\times n}$, and $\bm{\varepsilon}_t \sim \text{WN}(\mathbf{0}, \Sigma)$.

The VAR($p$) can be rewritten as a VAR(1) in the **companion form** by stacking:
$$\underbrace{\begin{pmatrix}\mathbf{X}_t \\ \mathbf{X}_{t-1} \\ \vdots \\ \mathbf{X}_{t-p+1}\end{pmatrix}}_{\tilde{\mathbf{X}}_t} = \underbrace{\begin{pmatrix}\Phi_1 & \Phi_2 & \cdots & \Phi_p \\ I & 0 & \cdots & 0 \\ 0 & I & \cdots & 0 \\ \vdots & & \ddots & \vdots\end{pmatrix}}_{\tilde{A}} \underbrace{\begin{pmatrix}\mathbf{X}_{t-1} \\ \mathbf{X}_{t-2} \\ \vdots \\ \mathbf{X}_{t-p}\end{pmatrix}}_{\tilde{\mathbf{X}}_{t-1}} + \begin{pmatrix}\bm{\varepsilon}_t \\ \mathbf{0} \\ \vdots \\ \mathbf{0}\end{pmatrix}.$$

**Stationarity of the VAR($p$):** The VAR($p$) is covariance-stationary iff all eigenvalues of the companion matrix $\tilde{A}$ lie inside the unit circle.

The companion form shows that any DSGE model's solution — which takes the form $\mathbf{y}_t = C\mathbf{y}_{t-1} + D\mathbf{z}_t$ after solving (Chapter 28) — is a first-order VAR with companion matrix $C$. The stationarity condition for the DSGE model solution is exactly the Blanchard–Kahn condition applied to $C$.

---

## 5.7 The Spectral Density

The **spectral density** decomposes the variance of a process by frequency, revealing which periodicities are most important.

**Definition 5.11 (Spectral Density).** For a covariance-stationary process $\{X_t\}$ with absolutely summable autocovariances ($\sum_{k=-\infty}^\infty|\gamma_k| < \infty$), the **spectral density** is:
$$S_X(\omega) = \frac{1}{2\pi}\sum_{k=-\infty}^\infty\gamma_k e^{-ik\omega}, \quad \omega \in [-\pi, \pi].$$

$S_X(\omega)d\omega$ is the contribution to variance from cyclical components with frequencies in $[\omega, \omega+d\omega]$. The spectral density is the Fourier transform of the autocovariance function; by the inverse transform:
$$\gamma_k = \int_{-\pi}^\pi S_X(\omega)e^{ik\omega}d\omega, \quad \sigma^2 = \gamma_0 = \int_{-\pi}^\pi S_X(\omega)d\omega.$$

**Theorem 5.3 (Spectral Density of an AR(1)).** For $X_t = \rho X_{t-1} + \varepsilon_t$ with $\varepsilon_t \sim \text{WN}(0,\sigma_\varepsilon^2)$:
$$S_X(\omega) = \frac{\sigma_\varepsilon^2/2\pi}{|1 - \rho e^{-i\omega}|^2} = \frac{\sigma_\varepsilon^2/2\pi}{1 + \rho^2 - 2\rho\cos\omega}.$$

*Proof.* Using the MA($\infty$) representation $X_t = \sum_{j=0}^\infty\rho^j\varepsilon_{t-j}$ and the spectral density of white noise $S_\varepsilon(\omega) = \sigma_\varepsilon^2/(2\pi)$:
$$S_X(\omega) = \left|\sum_{j=0}^\infty\rho^je^{-ij\omega}\right|^2 S_\varepsilon(\omega) = \frac{\sigma_\varepsilon^2/2\pi}{|1-\rho e^{-i\omega}|^2}. \quad\square$$

For $\rho > 0$, $S_X(\omega)$ is largest at $\omega = 0$ (low frequency, long cycles) — confirming that high persistence generates power at business cycle frequencies. The HP filter exploits this by attenuating the low-frequency components; its frequency response function is closely related to the spectral density of a smooth trend.

---

## 5.8 Martingales and the Innovation Representation

**Definition 5.12 (Martingale).** A stochastic process $\{M_t\}$ adapted to $\{\mathcal{F}_t\}$ is a **martingale** if:
$$\mathbb{E}[M_{t+1}|\mathcal{F}_t] = M_t \quad \text{for all } t.$$

Equivalently, $\mathbb{E}[M_{t+h}|\mathcal{F}_t] = M_t$ for all $h \geq 0$: the best forecast of future values is the current value.

A **martingale difference sequence** (MDS) $\{D_t\}$ satisfies $\mathbb{E}[D_t|\mathcal{F}_{t-1}] = 0$ — innovations are unpredictable from past information. The Wold innovations $\varepsilon_t$ are an MDS: $\mathbb{E}[\varepsilon_t|\mathcal{F}_{t-1}] = 0$.

Hall's (1978) **random walk hypothesis** for consumption [P:Ch.11.3] states that consumption is a martingale: $\mathbb{E}_t[c_{t+1}] = c_t$. This follows directly from the Euler equation $u'(c_t) = \beta(1+r)\mathbb{E}_t[u'(c_{t+1})]$ under quadratic utility ($u(c) = -(c-c^*)^2/2$, implying $u'(c) = c^*-c$): $c_t = \mathbb{E}_t[c_{t+1}]$. Therefore $c_{t+1} - c_t = \varepsilon_{t+1}$, which is an MDS — changes in consumption are unforecastable.

This is an empirically testable restriction: no variable in $\mathcal{F}_t$ should help predict $c_{t+1} - c_t$. The widespread violation of this restriction (**excess sensitivity** of consumption to predictable income changes [P:Ch.11.4]) is one of the major puzzles in consumption economics, motivating the buffer-stock and liquidity-constraint models.

---

## 5.9 Worked Example: Fitting an AR(1) to U.S. TFP

*Cross-reference: Principles Ch. 27 (RBC calibration)* **[P:Ch.27]**

The standard RBC calibration [P:Ch.27] sets the TFP shock as $\ln A_t = \rho_A\ln A_{t-1} + \varepsilon_t^A$. We estimate $\rho_A$ and $\sigma_A$ from the Solow residual.

**Data:** Quarterly U.S. TFP index, 1953Q1–2019Q4, obtained from the San Francisco Fed's TFP database. Define $\hat{a}_t = \ln A_t - \ln\bar{A}$ (log-deviation from mean).

**OLS estimation:** Regress $\hat{a}_t$ on $\hat{a}_{t-1}$:
$$\hat{a}_t = \rho_A\hat{a}_{t-1} + \varepsilon_t^A.$$

The OLS estimator:
$$\hat{\rho}_A = \frac{\sum_{t=2}^T\hat{a}_t\hat{a}_{t-1}}{\sum_{t=2}^T\hat{a}_{t-1}^2} = \frac{\gamma_1}{\gamma_0} = \hat{\rho}_1.$$

The OLS estimator of the AR(1) coefficient is the first sample autocorrelation. Typical estimates: $\hat{\rho}_A \approx 0.95$, $\hat{\sigma}_A \approx 0.0072$ (quarterly).

**Interpretation:** A TFP shock of $+1\%$ decays at rate $0.95^t$ per quarter. After 5 years (20 quarters): $0.95^{20} = 0.358$ — 36% of the original shock remains. The half-life is $\ln 2/|\ln 0.95| = 0.693/0.051 \approx 13.5$ quarters, or about 3.4 years.

```apl
⍝ APL — AR(1) estimation and simulation
⎕IO←0 ⋄ ⎕ML←1

⍝ Simulate AR(1) with known parameters (TFP calibration)
rho_A ← 0.95  ⋄  sigma_A ← 0.0072  ⋄  T ← 200

⍝ Draw white noise innovations (using external random)
⍝ In production: use ⎕PY to call numpy.random.normal
⍝ Here: simple approximation via Box-Muller not shown, use ⎕RL-based

⍝ AR(1) simulation via scan: a_t = rho * a_{t-1} + eps_t
⍝ Collect path starting from 0
eps ← sigma_A × {0.5 - ?0} ¨ ⍳ T    ⍝ rough uniform approx (replace with normal)
ar1_path ← {rho_A×⊃⌽⍵, ⊃⍵+eps[⍺]}\ ⍳T    ⍝ scan with index

⍝ More idiomatically: use ⍣ to iterate forward
⍝ sim_ar1 ← {rho_A×⍵ + sigma_A×rnorm 1}⍣T ⊢ 0

⍝ Sample autocorrelation at lag 1
⍝ rho_hat = Cov(a_t, a_{t-1}) / Var(a_t)
demean ← {⍵ - (+/⍵)÷≢⍵}
a ← demean ar1_path
rho_hat ← {(⍺+.×⍵) ÷ ⍺+.×⍺} (1↓a) (¯1↓a)    ⍝ (1↓a)·(¯1↓a) / (¯1↓a)·(¯1↓a)
rho_hat    ⍝ should be ≈ 0.95

⍝ Autocovariance function (first 10 lags)
acf ← {(k↓a)+.×(k↓⌽a)}∘{+.×⍨}¨⍳10    ⍝ approximate
```

```python
# Python — AR(1) estimation and diagnostics
import numpy as np
from statsmodels.tsa.ar_model import AutoReg
import matplotlib.pyplot as plt

np.random.seed(42)
rho_true, sigma_true, T = 0.95, 0.0072, 200

# Simulate AR(1)
eps = np.random.normal(0, sigma_true, T)
a = np.zeros(T)
for t in range(1, T):
    a[t] = rho_true * a[t-1] + eps[t]

# OLS estimation: first autocorrelation
rho_hat = np.corrcoef(a[1:], a[:-1])[0,1]
sigma_hat = np.std(a[1:] - rho_hat * a[:-1])
print(f"True: ρ={rho_true}, σ={sigma_true}")
print(f"OLS:  ρ̂={rho_hat:.4f}, σ̂={sigma_hat:.6f}")

# Impulse response
h = np.arange(40)
irf = rho_hat ** h
plt.plot(h, irf); plt.axhline(0, c='k', lw=0.5)
plt.title('AR(1) Impulse Response Function'); plt.xlabel('Quarters'); plt.show()
```

```julia
# Julia — AR(1) moments and simulation
using Statistics

rho_true, sigma_true, T = 0.95, 0.0072, 200
Random.seed!(42)
eps = randn(T) .* sigma_true
a = zeros(T)
for t in 2:T
    a[t] = rho_true * a[t-1] + eps[t]
end

rho_hat = cor(a[2:end], a[1:end-1])
sigma_hat = std(a[2:end] .- rho_hat .* a[1:end-1])
println("OLS ρ̂ = $(round(rho_hat, digits=4)),  σ̂ = $(round(sigma_hat, digits=6))")

# Theoretical vs. simulated ACF
lags = 0:20
acf_theory = rho_hat .^ lags
acf_sample = [cor(a[1+k:end], a[1:end-k]) for k in lags]
println("ACF comparison (lags 0-5):")
display([acf_theory[1:6] acf_sample[1:6]])
```

```r
# R — AR(1) estimation
set.seed(42)
rho_true <- 0.95; sigma_true <- 0.0072; T <- 200
eps <- rnorm(T, 0, sigma_true)
a <- numeric(T)
for(t in 2:T) a[t] <- rho_true * a[t-1] + eps[t]

# Yule-Walker / OLS estimate
rho_hat <- cor(a[-1], a[-T])
sigma_hat <- sd(a[-1] - rho_hat * a[-T])
cat(sprintf("OLS: rho_hat = %.4f, sigma_hat = %.6f\n", rho_hat, sigma_hat))

# Theoretical spectral density
omega <- seq(-pi, pi, length=500)
spec_theory <- sigma_hat^2 / (2*pi * (1 + rho_hat^2 - 2*rho_hat*cos(omega)))
plot(omega, spec_theory, type='l', main='AR(1) Spectral Density',
     xlab='Frequency ω', ylab='S(ω)')
```

---

## 5.10 The Tauchen (1986) Discretization

Continuous-state AR(1) processes must be discretized for numerical dynamic programming (Chapters 15–17). The Tauchen (1986) method constructs a discrete Markov chain that approximates the AR(1).

**Algorithm 5.1 (Tauchen Discretization).**

Given the AR(1): $z_{t+1} = \rho z_t + \sigma_\varepsilon\varepsilon_{t+1}$, $\varepsilon_t \sim \mathcal{N}(0,1)$:

1. Choose the number of grid points $N$ and the number of standard deviations $m$ (typically $m = 3$).
2. Set grid endpoints: $z_{min} = -m\sigma_z$, $z_{max} = m\sigma_z$, where $\sigma_z = \sigma_\varepsilon/\sqrt{1-\rho^2}$.
3. Construct evenly spaced grid: $\{z_1, \ldots, z_N\}$ with $z_1 = z_{min}$, $z_N = z_{max}$, spacing $\Delta z = (z_{max}-z_{min})/(N-1)$.
4. Compute transition probabilities: for grid points $z_i$ and $z_j$:
$$P_{ij} = \Phi\!\left(\frac{z_j + \Delta z/2 - \rho z_i}{\sigma_\varepsilon}\right) - \Phi\!\left(\frac{z_j - \Delta z/2 - \rho z_i}{\sigma_\varepsilon}\right),$$
with boundary corrections for $j = 1$ (lower tail) and $j = N$ (upper tail). Here $\Phi$ is the standard normal CDF.

The transition matrix $P$ is an $N\times N$ row-stochastic matrix: each row sums to 1.

In APL, the Tauchen discretization exploits the outer product `∘.-` to compute all differences simultaneously:

```apl
⍝ APL — Tauchen (1986) discretization in ~10 lines
⎕IO←0 ⋄ ⎕ML←1

tauchen ← {
    rho sigma_eps N m ← ⍵
    sigma_z ← sigma_eps ÷ (1 - rho*2)*0.5
    z_min ← (-m) × sigma_z
    z_max ← m × sigma_z
    z ← z_min + (z_max - z_min) × (⍳N) ÷ N-1    ⍝ grid
    dz  ← (z_max - z_min) ÷ N-1                   ⍝ grid spacing

    ⍝ Centred breakpoints: z_j ± dz/2, adjusted for AR(1) mean
    ⍝ All differences: z[j] - rho*z[i] for all (i,j) pairs
    mu_mat ← rho ∘.× z                   ⍝ N×N matrix of rho*z_i
    diff   ← (⍉z) ∘.- mu_mat             ⍝ N×N: z_j - rho*z_i

    ⍝ Normal CDF approximation (use ⎕PY in production)
    Phi ← {0.5 × 1 + 2○⍵÷2*0.5}         ⍝ rough approx; replace with scipy.stats.norm.cdf

    P ← (Phi (diff + dz÷2) ÷ sigma_eps) - (Phi (diff - dz÷2) ÷ sigma_eps)

    ⍝ Boundary corrections
    P[;0]   ← Phi (z[0]   + dz÷2 - rho ∘.× z) ÷ sigma_eps
    P[;N-1] ← 1 - +⌿ P[;⍳N-1]

    z P    ⍝ return grid and transition matrix
}
```

The Rouwenhorst (1995) method is an alternative that better approximates highly persistent processes ($|\rho|$ near 1) and is preferred for calibrated DSGE models where $\rho \approx 0.95$.

---

## 5.11 Programming Exercises

### Exercise 5.1 (APL — ACF)

Write a dfn `acf ← {n ← ⍺ ⋄ x ← ⍵-+/⍵÷≢⍵ ⋄ ...}` that computes the sample autocorrelation function at lags $0, 1, \ldots, n-1$ for a time series passed as right argument. The computation should use `+.×` (inner product) for each lag. Test on the simulated AR(1) from Section 5.9 and verify the ACF matches $\hat{\rho}^k$.

```apl
⎕IO←0 ⋄ ⎕ML←1

acf ← {
    n ← ⍺
    x ← ⍵ - (+/⍵)÷≢⍵        ⍝ demean
    v ← x+.×x                  ⍝ sum of squares (unnormalized variance)
    {(k↓x)+.×(k↓⌽x)}¨⍳n) ÷ v  ⍝ autocov at each lag, normalized
}

⍝ Test on AR(1) simulation
T ← 500
rho ← 0.8
eps ← 0.1 × (T⍴0) + ? T⍴0    ⍝ rough noise (replace with normal)
⍝ ar1 ← {rho×⍵+eps[⍺]}\ ⍳T   ⍝ would need proper simulation
20 acf ar1                      ⍝ first 20 autocorrelations
```

### Exercise 5.2 (Python — Spectral Density)

Plot the theoretical and estimated spectral density of the AR(1) TFP process calibrated with $\rho = 0.95$ and $\sigma_\varepsilon = 0.0072$. Use numpy's FFT to estimate the empirical spectrum from a long simulation, and overlay the theoretical formula from Theorem 5.3. Shade the business-cycle frequency band (cycles of 6–32 quarters).

```python
import numpy as np; import matplotlib.pyplot as plt
rho, sigma, T = 0.95, 0.0072, 10000
eps = np.random.normal(0, sigma, T)
a = np.zeros(T)
for t in range(1,T): a[t] = rho*a[t-1] + eps[t]

omega = np.linspace(-np.pi, np.pi, 1000)
spec_theory = sigma**2 / (2*np.pi * (1 + rho**2 - 2*rho*np.cos(omega)))

from scipy.signal import periodogram
f, Pxx = periodogram(a, window='hann')
omega_data = 2*np.pi*f

fig, ax = plt.subplots()
ax.semilogy(omega_data, Pxx, alpha=0.5, label='Empirical (periodogram)')
ax.semilogy(omega, spec_theory, 'r-', linewidth=2, label='Theoretical')
ax.axvspan(2*np.pi/32, 2*np.pi/6, alpha=0.1, color='green', label='BC frequencies')
ax.legend(); ax.set_xlabel('Frequency ω'); plt.show()
```

### Exercise 5.3 (Julia — Tauchen vs. Rouwenhorst)

```julia
using Distributions, LinearAlgebra

function tauchen(rho, sigma_eps, N, m=3)
    sigma_z = sigma_eps / sqrt(1-rho^2)
    z = range(-m*sigma_z, m*sigma_z, length=N)
    dz = step(z)
    P = zeros(N, N)
    d = Normal(0, sigma_eps)
    for i in 1:N
        mu = rho * z[i]
        P[i,1] = cdf(d, z[1] + dz/2 - mu)
        P[i,N] = 1 - cdf(d, z[N] - dz/2 - mu)
        for j in 2:N-1
            P[i,j] = cdf(d, z[j]+dz/2-mu) - cdf(d, z[j]-dz/2-mu)
        end
    end
    return collect(z), P
end

z_grid, P = tauchen(0.95, 0.0072, 7)
println("Grid: ", round.(z_grid, digits=4))
println("Row sums (should be 1): ", round.(sum(P, dims=2), digits=6))
```

### Exercise 5.4 (R — Wold Representation)

```r
# Verify Wold MA coefficients for AR(1) match IRF
rho <- 0.8; sigma <- 0.1; T <- 500
set.seed(1)
a <- arima.sim(list(ar=rho), T, sd=sigma)

# Fit AR(1)
fit <- arima(a, order=c(1,0,0))
rho_hat <- as.numeric(coef(fit)[1])

# Wold MA coefficients: psi_j = rho^j
lags <- 0:20
psi <- rho_hat^lags

# Verify via impulse response from ARMAtoMA
psi_check <- ARMAtoMA(ar=rho_hat, lag.max=20)
cat("Max discrepancy:", max(abs(psi[-1] - psi_check)), "\n")
```

### Exercise 5.5 — Persistence and the HP Filter ($\star$)

The Hodrick–Prescott (HP) filter with parameter $\lambda = 1600$ (quarterly) isolates components with period 6–32 quarters. Show that for an AR(1) with $\rho = 0.95$, approximately what fraction of the total variance falls in the business-cycle band $[2\pi/32, 2\pi/6]$? Compute this as $\int_{2\pi/32}^{2\pi/6} 2S_X(\omega)d\omega / \sigma_X^2$ using numerical integration (Gaussian quadrature). How does this fraction change as $\rho$ varies from 0 to 0.99?

### Exercise 5.6 — Hall's Test ($\star\star$)

Using quarterly U.S. non-durable consumption growth data from FRED (series `DNDGRD3Q086SBEA`):
(a) Test whether lagged consumption growth $\Delta\ln c_{t-1}$ helps predict current consumption growth $\Delta\ln c_t$ (Hall's test for the martingale property).
(b) Test whether lagged income growth $\Delta\ln y_{t-1}$ has predictive power (excess sensitivity test).
(c) Use HAC (Newey–West) standard errors for both tests. What do your results imply about the extent of liquidity constraints in the data?

---

## 5.12 Chapter Summary

**Key results:**

- A **probability space** $(\Omega, \mathcal{F}, P)$ with filtration $\{\mathcal{F}_t\}$ formalizes the structure of information in dynamic economies. Rational expectations = conditional expectation given $\mathcal{F}_t$.
- A process is **covariance-stationary** if its mean, variance, and autocovariances are all time-invariant. The autocovariance function (ACF) and spectral density are its complete second-moment characterizations.
- The **AR(1)** $X_t = \rho X_{t-1} + \varepsilon_t$ has ACF $\rho_k = \rho^k$, variance $\sigma_\varepsilon^2/(1-\rho^2)$, and spectral density $\sigma_\varepsilon^2/(2\pi|1-\rho e^{-i\omega}|^2)$.
- The **Wold decomposition** guarantees that any covariance-stationary process has an MA($\infty$) representation in terms of its own innovations — the foundation of structural VAR identification.
- The **Tauchen (1986) algorithm** discretizes an AR(1) into a finite Markov chain using the `∘.-` outer product in APL, enabling numerical dynamic programming.
- **Hall's random walk** — consumption is a martingale under rational expectations and quadratic utility — is the testable implication of the consumption Euler equation; its widespread empirical failure points to liquidity constraints.

**Connections forward:** Chapter 18 develops the full rational expectations solution methodology using these stochastic concepts. Chapter 19 uses ARMA and VAR processes for structural identification of macroeconomic shocks. Chapter 20 applies the Kalman filter to estimate unobserved state variables (potential output, the natural rate) from noisy data. Chapter 26 uses Monte Carlo simulation of AR(1) shock processes to generate the empirical moments of calibrated DSGE models.

---

*Next: Part II — Static Macroeconomic Models and Comparative Statics*
