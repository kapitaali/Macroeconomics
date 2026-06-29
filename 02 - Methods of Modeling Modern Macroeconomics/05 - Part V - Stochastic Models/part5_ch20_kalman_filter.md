# Chapter 20: The Kalman Filter

*Estimating Unobserved Components (Potential GDP, NAIRU)*

> *"The Kalman filter is the Swiss Army knife of state estimation: it handles missing data, real-time updating, and likelihood evaluation — all in a unified recursive algorithm."*

**Cross-reference:** *Principles* Ch. 6 (potential GDP and the output gap as unobserved); Appendix B (DSGE estimation via MLE with Kalman filter) **[P:Ch.6, P:AppB]**

---

## 20.1 The State-Space Form

Many macroeconomic concepts are not directly observed: potential output, the natural rate of unemployment (NAIRU), the natural rate of interest, and long-run inflation expectations are all inherently unobservable quantities inferred from observed data. The **state-space model** provides the unified statistical framework for this inference; the **Kalman filter** is the optimal algorithm for doing so recursively in real time.

**Definition 20.1 (Linear Gaussian State-Space Model).** The **state-space model** consists of two equations:

$$\mathbf{y}_t = H\bm\alpha_t + \mathbf{d} + \bm\varepsilon_t, \quad \bm\varepsilon_t \sim \mathcal{N}(\mathbf{0}, R) \quad \text{(measurement equation)}$$
$$\bm\alpha_{t+1} = F\bm\alpha_t + \mathbf{c} + \bm\eta_t, \quad \bm\eta_t \sim \mathcal{N}(\mathbf{0}, Q) \quad \text{(transition equation)}$$

where:
- $\mathbf{y}_t \in \mathbb{R}^p$ is the **observation vector** (observed variables: GDP, unemployment, interest rates).
- $\bm\alpha_t \in \mathbb{R}^m$ is the **state vector** (unobserved variables: potential output, NAIRU, natural rate).
- $H \in \mathbb{R}^{p\times m}$ is the **observation matrix** (links states to observables).
- $F \in \mathbb{R}^{m\times m}$ is the **transition matrix** (law of motion of the state).
- $R \in \mathbb{R}^{p\times p}$ is the **measurement noise covariance**.
- $Q \in \mathbb{R}^{m\times m}$ is the **state noise covariance** (process noise).
- $\mathbf{d}$ and $\mathbf{c}$ are constant terms (intercepts).

**Assumption:** $\bm\varepsilon_t$ and $\bm\eta_t$ are mutually uncorrelated at all leads and lags; the initial state $\bm\alpha_1 \sim \mathcal{N}(\mathbf{a}_1, P_1)$.

The state-space form is extremely general: every ARMA model, every DSGE model (after log-linearization), and every unobserved components (UC) model can be written in this form. Chapter 28 will show that the log-linearized DSGE model reduces to exactly this state-space form, with the Kalman filter providing the likelihood.

---

## 20.2 The Kalman Filter Algorithm

The Kalman filter solves the **linear filtering problem**: given observations $\mathbf{y}_1, \ldots, \mathbf{y}_t$, what is the best (minimum mean squared error) estimate of the unobserved state $\bm\alpha_t$?

**Definition 20.2 (Filtered Estimate).** The **filtered mean** $\mathbf{a}_{t|t} \equiv \mathbb{E}[\bm\alpha_t | \mathbf{y}_1, \ldots, \mathbf{y}_t]$ and **filtered covariance** $P_{t|t} \equiv \text{Var}[\bm\alpha_t | \mathbf{y}_1, \ldots, \mathbf{y}_t]$ are the posterior mean and variance of the state conditional on all information up to time $t$.

The **predicted mean** $\mathbf{a}_{t|t-1} \equiv \mathbb{E}[\bm\alpha_t|\mathbf{y}_1,\ldots,\mathbf{y}_{t-1}]$ and **predicted covariance** $P_{t|t-1} \equiv \text{Var}[\bm\alpha_t|\mathbf{y}_1,\ldots,\mathbf{y}_{t-1}]$ are the prior to the $t$-th observation.

**Algorithm 20.1 (Kalman Filter).**

**Initialization:** $\mathbf{a}_{1|0} = \mathbf{a}_1$, $P_{1|0} = P_1$.

For $t = 1, 2, \ldots, T$:

**Prediction step** (propagate the prior through the transition equation):
$$\mathbf{a}_{t|t-1} = F\mathbf{a}_{t-1|t-1} + \mathbf{c}$$
$$P_{t|t-1} = FP_{t-1|t-1}F' + Q$$

**Prediction error and variance:**
$$\mathbf{v}_t = \mathbf{y}_t - H\mathbf{a}_{t|t-1} - \mathbf{d} \quad \text{(innovation)}$$
$$F_t = HP_{t|t-1}H' + R \quad \text{(innovation variance)}$$

**Kalman gain:**
$$K_t = P_{t|t-1}H'F_t^{-1}$$

**Update step** (incorporate the new observation):
$$\mathbf{a}_{t|t} = \mathbf{a}_{t|t-1} + K_t\mathbf{v}_t$$
$$P_{t|t} = (I - K_tH)P_{t|t-1}$$

In APL, the prediction and update steps are single matrix expressions:

```apl
⍝ APL — Kalman filter: prediction and update steps
⎕IO←0 ⋄ ⎕ML←1

⍝ Prediction step: propagate state and covariance
⍝ a_pred = F × a + c
⍝ P_pred = F × P × F' + Q
predict ← {F Q c a P ← ⍵
    a_new ← (F +.× a) + c
    P_new ← (F +.× P +.× ⍉F) + Q
    a_new P_new}

⍝ Update step: incorporate new observation
⍝ v = y - H×a - d  (innovation)
⍝ Fvar = H × P × H' + R  (innovation variance)
⍝ K = P × H' × Fvar^{-1}  (Kalman gain)
⍝ a_upd = a + K×v
⍝ P_upd = (I - K×H) × P
update ← {H R d y a P ← ⍵
    v    ← y - (H +.× a) + d        ⍝ innovation
    Fvar ← (H +.× P +.× ⍉H) + R    ⍝ innovation variance
    K    ← P +.× (⍉H) +.× ⌹Fvar    ⍝ Kalman gain: P H' F^{-1}
    n    ← ≢a
    I_n  ← =⍨⍳n
    a_upd ← a + K +.× v
    P_upd ← (I_n - K +.× H) +.× P
    a_upd P_upd v Fvar}              ⍝ return updated state, cov, innovation, innov var
```

The full Kalman filter iterates these two steps over $T$ observations, storing the filtered estimates $\{\mathbf{a}_{t|t}, P_{t|t}\}$ and the innovations $\{v_t, F_t\}$.

---

## 20.3 Derivation of the Kalman Gain as MMSE Estimator

**Theorem 20.1 (Optimality of the Kalman Gain).** The Kalman gain $K_t = P_{t|t-1}H'F_t^{-1}$ minimizes the mean squared error of the updated state estimate among all linear functions of $(\mathbf{a}_{t|t-1}, \mathbf{y}_t)$.

*Proof.* The general linear update: $\hat{\bm\alpha}_{t|t} = \mathbf{a}_{t|t-1} + G\mathbf{v}_t$ for some $m\times p$ matrix $G$. The update error:

$$\bm\alpha_t - \hat{\bm\alpha}_{t|t} = (\bm\alpha_t - \mathbf{a}_{t|t-1}) - G\mathbf{v}_t.$$

Since $\mathbf{v}_t = H(\bm\alpha_t - \mathbf{a}_{t|t-1}) + \bm\varepsilon_t$:

$$\bm\alpha_t - \hat{\bm\alpha}_{t|t} = (I - GH)(\bm\alpha_t - \mathbf{a}_{t|t-1}) - G\bm\varepsilon_t.$$

The MSE matrix:

$$\text{MSE}(G) = (I-GH)P_{t|t-1}(I-GH)' + GRG'.$$

Taking the derivative with respect to $G$ and setting to zero:

$$\frac{\partial\text{tr}[\text{MSE}(G)]}{\partial G} = -2P_{t|t-1}H'(I-GH)' + 2GR = 0.$$

Solving: $P_{t|t-1}H' = G(HP_{t|t-1}H' + R) = GF_t$, so $G^* = P_{t|t-1}H'F_t^{-1} = K_t$. $\square$

**The Joseph form for numerical stability:** The update equation $P_{t|t} = (I-K_tH)P_{t|t-1}$ can become non-positive-definite due to numerical errors. The numerically stable **Joseph form** is:

$$P_{t|t} = (I-K_tH)P_{t|t-1}(I-K_tH)' + K_tRK_t'.$$

This is algebraically equivalent but guaranteed to remain symmetric positive semi-definite.

---

## 20.4 The Kalman Smoother

The Kalman filter produces **filtered** estimates $\mathbf{a}_{t|t}$ — based on information up to time $t$. The **Kalman smoother** produces **smoothed** estimates $\mathbf{a}_{t|T}$ — based on all $T$ observations. Smoothed estimates are generally more accurate and are required for computing the log-likelihood efficiently.

**Algorithm 20.2 (Rauch–Tung–Striebel Smoother).**

Run the Kalman filter forward to obtain $\{\mathbf{a}_{t|t}, P_{t|t}\}$ for $t = 1, \ldots, T$.

Then run backward for $t = T-1, T-2, \ldots, 1$:

$$L_t = P_{t|t}F'P_{t+1|t}^{-1} \quad \text{(smoother gain)}$$
$$\mathbf{a}_{t|T} = \mathbf{a}_{t|t} + L_t(\mathbf{a}_{t+1|T} - \mathbf{a}_{t+1|t})$$
$$P_{t|T} = P_{t|t} + L_t(P_{t+1|T} - P_{t+1|t})L_t'$$

with terminal condition $\mathbf{a}_{T|T}$ and $P_{T|T}$ from the filter.

**Interpretation:** The smoother revises the filtered estimate at each date using information from all subsequent observations. For the output gap, the smoother tells us what we would have estimated for 2006Q2 if we had known the entire 2007–2012 sequence of GDP data — a retrospective revision.

---

## 20.5 Maximum Likelihood Estimation via the Kalman Filter

The Kalman filter generates the **prediction error decomposition** of the log-likelihood — the key formula connecting the filter to parameter estimation.

**Theorem 20.2 (Prediction Error Decomposition of Log-Likelihood).** Under the Gaussian assumption, the log-likelihood of the state-space model parameters $\bm\theta = (F, H, Q, R, \mathbf{a}_1, P_1)$ given observations $\{\mathbf{y}_t\}_{t=1}^T$ is:

$$\mathcal{L}(\bm\theta) = -\frac{Tp}{2}\ln(2\pi) - \frac{1}{2}\sum_{t=1}^T\left[\ln|\mathbf{F}_t| + \mathbf{v}_t'\mathbf{F}_t^{-1}\mathbf{v}_t\right],$$

where $\mathbf{v}_t = \mathbf{y}_t - H\mathbf{a}_{t|t-1} - \mathbf{d}$ and $\mathbf{F}_t = HP_{t|t-1}H' + R$ are the prediction errors and their covariances generated by the Kalman filter.

*Proof sketch.* By the chain rule of conditional probability:
$$p(\mathbf{y}_1, \ldots, \mathbf{y}_T) = \prod_{t=1}^T p(\mathbf{y}_t | \mathbf{y}_{t-1}, \ldots, \mathbf{y}_1).$$
Each conditional $p(\mathbf{y}_t | \cdot)$ is Gaussian with mean $H\mathbf{a}_{t|t-1} + \mathbf{d}$ (the predicted observation) and variance $\mathbf{F}_t$. Taking logs gives the stated formula. $\square$

**Practical estimation:** MLE maximizes $\mathcal{L}(\bm\theta)$ over model parameters $\bm\theta$. The gradient $\partial\mathcal{L}/\partial\bm\theta$ can be computed analytically (or by automatic differentiation) and used with a quasi-Newton optimizer (BFGS, Chapter 24). The Kalman filter evaluates $\mathcal{L}$ at each candidate $\bm\theta$ in $O(Tm^3)$ time — fast enough for models with state dimension $m \lesssim 50$.

---

## 20.6 Macroeconomic Applications

### 20.6.1 Estimating the Output Gap

The simplest unobserved components (UC) model for the output gap:

$$y_t = \tau_t + c_t \quad \text{(output = trend + cycle)}$$
$$\tau_{t+1} = \tau_t + g_t + \eta_t^\tau, \quad g_{t+1} = g_t + \eta_t^g \quad \text{(random walk trend + drift)}$$
$$c_{t+1} = \phi_1 c_t + \phi_2 c_{t-1} + \eta_t^c \quad \text{(AR(2) cycle)}$$

State vector: $\bm\alpha_t = (\tau_t, g_t, c_t, c_{t-1})'$. Observation: $y_t = \tau_t + c_t$.

The four unknowns are the signal-to-noise ratios $\{q_\tau = \sigma_\tau^2/\sigma_c^2, q_g = \sigma_g^2/\sigma_c^2\}$ and AR coefficients $(\phi_1, \phi_2)$. MLE via the Kalman filter estimates these from GDP data.

### 20.6.2 The Laubach–Williams Natural Rate Model

Laubach and Williams (2003) estimate the natural rate of interest $r^n_t$ as an unobserved state in a system that links $r^n_t$ to the output gap and inflation:

$$x_t = \mathbb{E}[x_{t+1}|t] - \sigma(r_t - r^n_t) + \varepsilon_t^x \quad \text{(IS curve)}$$
$$\pi_t = \pi_{t-1} + \alpha x_{t-1} + \varepsilon_t^\pi \quad \text{(backward-looking PC)}$$
$$r^n_t = c + g_t + z_t \quad \text{(unobserved natural rate)}$$

The system is cast in state-space form with unobserved $r^n_t$ and estimated by MLE. The LW estimates show $r^n_t$ declining from approximately 3.5% in the early 1980s to near 0% by 2015 — consistent with the secular stagnation narrative of *Principles* Ch. 39 [P:Ch.39.4].

---

## 20.7 Worked Example: Output Gap via UC Model

```python
import numpy as np
from scipy.optimize import minimize

# Simulate U.S.-like GDP data for demonstration
np.random.seed(42)
T = 200
tau = np.zeros(T); g = np.zeros(T); c = np.zeros(T)
g[0] = 0.005  # quarterly growth rate
sig_tau, sig_g, sig_c = 0.002, 0.001, 0.010
phi1, phi2 = 1.5, -0.6

for t in range(1, T):
    g[t]   = g[t-1] + sig_g * np.random.randn()
    tau[t] = tau[t-1] + g[t-1] + sig_tau * np.random.randn()
    c[t]   = (phi1*c[t-1] + phi2*(c[t-2] if t>1 else 0)
              + sig_c * np.random.randn())

y = tau + c + 0.003*np.random.randn(T)  # observed GDP (with tiny meas. noise)

def kalman_filter(y, F, H, Q, R, a0, P0):
    T, p = y.shape if y.ndim > 1 else (len(y), 1)
    m = len(a0)
    y = y.reshape(T, -1)
    
    a = a0.copy(); P = P0.copy()
    log_lik = 0.0
    a_filter = np.zeros((T, m)); P_filter = np.zeros((T, m, m))
    
    for t in range(T):
        # Predict
        a_pred = F @ a
        P_pred = F @ P @ F.T + Q
        # Innovation
        v = y[t] - H @ a_pred
        Fv = H @ P_pred @ H.T + R
        # Log-likelihood contribution
        sign, logdet = np.linalg.slogdet(Fv)
        log_lik -= 0.5 * (logdet + v.T @ np.linalg.solve(Fv, v) + p*np.log(2*np.pi))
        # Update
        K = P_pred @ H.T @ np.linalg.inv(Fv)
        a = a_pred + K @ v
        P = (np.eye(m) - K @ H) @ P_pred
        a_filter[t] = a; P_filter[t] = P
    
    return log_lik, a_filter, P_filter

# Set up UC state-space model
# State: [tau, g, c, c_lag]'
m = 4
F = np.array([[1, 1, 0, 0],
              [0, 1, 0, 0],
              [0, 0, phi1, phi2],
              [0, 0, 1, 0]])
H = np.array([[1, 0, 1, 0]])  # y = tau + c
R = np.array([[1e-6]])  # nearly zero measurement error

def neg_log_lik(params):
    sq, sg, sc = np.exp(params[:3])  # log-parameterize for positivity
    Q = np.diag([sq**2, sg**2, sc**2, 0])
    a0 = np.array([y[0], 0.005, 0, 0])
    P0 = np.diag([1.0, 0.01, 0.1, 0.1])
    ll, _, _ = kalman_filter(y, F, H, Q, R, a0, P0)
    return -ll

# Estimate parameters
x0 = np.log([0.002, 0.001, 0.010])
res = minimize(neg_log_lik, x0, method='Nelder-Mead',
               options={'xatol':1e-6, 'fatol':1e-6, 'maxiter':2000})
sig_hat = np.exp(res.x)
print(f"Estimated σ_τ={sig_hat[0]:.4f}, σ_g={sig_hat[1]:.4f}, σ_c={sig_hat[2]:.4f}")
print(f"True:     σ_τ={sig_tau:.4f}, σ_g={sig_g:.4f}, σ_c={sig_c:.4f}")

# Extract filtered output gap
Q_hat = np.diag([sig_hat[0]**2, sig_hat[1]**2, sig_hat[2]**2, 0])
a0 = np.array([y[0], 0.005, 0, 0]); P0 = np.diag([1.0,0.01,0.1,0.1])
_, a_f, _ = kalman_filter(y, F, H, Q_hat, R, a0, P0)
output_gap = a_f[:,2]  # cycle component

import matplotlib.pyplot as plt
fig,(ax1,ax2) = plt.subplots(2,1,figsize=(12,7),sharex=True)
ax1.plot(y,'b-',alpha=0.6,label='Observed GDP'); ax1.plot(a_f[:,0],'r-',lw=2,label='Filtered trend')
ax1.set_title('UC Model: Trend and Cycle Decomposition'); ax1.legend()
ax2.plot(output_gap,'g-',label='Output gap (filtered)'); ax2.plot(c,'k--',alpha=0.5,label='True cycle')
ax2.axhline(0,color='k',lw=0.5); ax2.set_title('Output Gap'); ax2.legend()
plt.tight_layout(); plt.show()
```

```julia
using LinearAlgebra, Optim

function kalman_filter(y, F, H, Q, R, a0, P0)
    T = length(y); m = length(a0); p = size(H, 1)
    a = copy(a0); P = copy(P0); log_lik = 0.0
    a_store = zeros(T, m)
    for t in 1:T
        a_pred = F * a; P_pred = F * P * F' + Q
        v = [y[t]] - H * a_pred
        Fv = H * P_pred * H' + R
        log_lik -= 0.5*(log(det(Fv)) + dot(v, Fv\v) + p*log(2π))
        K = P_pred * H' * inv(Fv)
        a = a_pred + K * v; P = (I(m) - K*H) * P_pred
        a_store[t,:] = a
    end
    log_lik, a_store
end

println("Kalman filter implemented. Use Optim.jl to maximize log-likelihood over parameters.")
```

```r
# R — Kalman filter via KFAS
library(KFAS)

# Simple local level model (trend extraction)
T <- 200; set.seed(42)
y <- cumsum(rnorm(T, 0.005, 0.002)) + rnorm(T, 0, 0.01)  # trend + noise

model <- SSModel(y ~ SSMtrend(2, Q=list(matrix(NA),matrix(NA))),
                 H=matrix(NA))
fit <- fitSSM(model, inits=c(log(0.001), log(0.001), log(0.01)))
out <- KFS(fit$model)

plot(out$alphahat[,1], type='l', col='red', main='Filtered trend (KFAS)')
lines(y, col='blue', alpha=0.5)
legend('topleft', c('Observed','Filtered trend'), col=c('blue','red'), lty=1)
```

---

## 20.8 Programming Exercises

### Exercise 20.1 (APL — Full Kalman Filter Loop)

Implement the complete Kalman filter in APL as a dfn that takes `(y F H Q R a0 P0)` and returns the filtered states and log-likelihood. Use `{predict_step ⍵}` and `{update_step ⍵}` as helper dfns called within the main loop. The log-likelihood accumulation should use `ln_det F_t + v'F^{-1}v` at each step. Verify on simulated data from a known state-space model.

### Exercise 20.2 (Python — Laubach–Williams Replication)

Implement a simplified Laubach–Williams (2003) model:

$$x_t = \alpha_1 x_{t-1} + \alpha_2(r_{t-1} - r^n_{t-1}) + \varepsilon_t^x$$
$$\pi_t = \beta\pi_{t-1} + \gamma x_{t-1} + \varepsilon_t^\pi$$  
$$r^n_t = c \cdot g_t, \quad g_{t+1} = g_t + \varepsilon_t^g$$

Cast in state-space form with state $(\hat x_t, \pi_t, r^n_t, g_t)'$ and observations $(x_t, \pi_t)'$. Estimate parameters by MLE via the Kalman filter. Plot the estimated natural rate path over 1960–2019 and compare to the LW published estimates.

### Exercise 20.3 (Julia — Kalman Smoother)

Extend the Kalman filter implementation to include the RTS smoother (backward pass). Compare: (a) filtered output gap $c_{t|t}$ (using only information up to $t$); (b) smoothed output gap $c_{t|T}$ (using all information). Show that the smoothed estimates have lower uncertainty (smaller $P_{t|T}$ diagonal elements) than the filtered estimates.

### Exercise 20.4 — DSGE Likelihood ($\star$)

The log-linearized NK model from Chapter 18 can be written in state-space form with state $\bm\alpha_t = (\hat x_t, \hat\pi_t, r^n_t)'$, observation $\mathbf{y}_t = (\hat x_t + e_t^x, \hat\pi_t + e_t^\pi)'$ (output gap and inflation observed with measurement error). (a) Write the $F$, $H$, $Q$, $R$ matrices explicitly in terms of the NK model parameters. (b) Simulate 100 periods of data from the model with $\phi_\pi = 1.5$. (c) Use the Kalman filter to compute $\mathcal{L}(\phi_\pi)$ for $\phi_\pi \in [0.5, 3.0]$ and plot the likelihood surface. Verify the MLE is near the true value.

---

## 20.9 Chapter Summary

**Key results:**

- The **state-space model** has a measurement equation $\mathbf{y}_t = H\bm\alpha_t + \bm\varepsilon_t$ and transition equation $\bm\alpha_{t+1} = F\bm\alpha_t + \bm\eta_t$; virtually all dynamic models reduce to this form.
- The **Kalman filter** is a two-step recursion: predict step ($\mathbf{a}_{t|t-1} = F\mathbf{a}_{t-1|t-1}$, $P_{t|t-1} = FPF'+Q$) and update step ($\mathbf{a}_{t|t} = \mathbf{a}_{t|t-1}+K_t\mathbf{v}_t$, $P_{t|t} = (I-K_tH)P_{t|t-1}$).
- The Kalman gain $K_t = P_{t|t-1}H'F_t^{-1}$ is the **MMSE linear estimator** (proved as Theorem 20.1).
- The **RTS smoother** revises filtered estimates backward in time, producing more accurate estimates $\mathbf{a}_{t|T}$ using all $T$ observations.
- The **prediction error log-likelihood** $\mathcal{L} = -\frac{1}{2}\sum_t[\ln|F_t| + v_t'F_t^{-1}v_t]$ enables MLE estimation of state-space model parameters via the Kalman filter.
- In APL: `P_pred ← (F+.×P+.×⍉F)+Q` and `K ← P+.×(⍉H)+.×⌹Fv` — the prediction and gain steps as single expressions; the full filter as `{update predict ⍵}⍣T⊢init_state`.

*Next: Chapter 21 — GMM and Maximum Likelihood*
