# Appendix E: Differential Equations Cheat Sheet

---

## E.1 First-Order ODEs

**General form:** $\dot{x} = f(x, t)$.

| Type | Form | Solution method |
|---|---|---|
| Separable | $\dot{x} = g(x)h(t)$ | Separate: $\int dx/g(x) = \int h(t)dt$ |
| Linear | $\dot{x} + p(t)x = q(t)$ | Integrating factor: $\mu = e^{\int p\,dt}$ |
| Bernoulli | $\dot{x} + p(t)x = q(t)x^n$ | Substitute $v=x^{1-n}$ → linear |
| Autonomous | $\dot{x} = f(x)$ | Phase line analysis; $x^*$ where $f(x^*)=0$ |

**Linear first-order (constant coefficients):** $\dot{x} = ax + b$.
- Steady state: $x^* = -b/a$ (if $a\neq0$).
- General solution: $x(t) = (x_0 - x^*)e^{at} + x^*$.
- Stable iff $a < 0$.

**Stability:** Linearize $\dot{x} = f(x)$ at $x^*$: $\dot{u} \approx f'(x^*)u$. Stable iff $f'(x^*) < 0$.

## E.2 Second-Order ODEs

**Linear, constant coefficients:** $\ddot{x} + p\dot{x} + qx = r$.

**Characteristic equation:** $\lambda^2 + p\lambda + q = 0$, roots $\lambda_{1,2} = (-p \pm \sqrt{p^2-4q})/2$.

| Case | Roots | General solution |
|---|---|---|
| Distinct real ($p^2>4q$) | $\lambda_1\neq\lambda_2\in\mathbb{R}$ | $C_1e^{\lambda_1 t}+C_2e^{\lambda_2 t}+x^*$ |
| Repeated ($p^2=4q$) | $\lambda = -p/2$ | $(C_1+C_2 t)e^{\lambda t}+x^*$ |
| Complex ($p^2<4q$) | $\lambda = \alpha\pm i\beta$ | $e^{\alpha t}(C_1\cos\beta t + C_2\sin\beta t)+x^*$ |

where $\alpha = -p/2$, $\beta = \sqrt{4q-p^2}/2$, $x^* = r/q$ (particular solution).

**Stability:** All trajectories → $x^*$ iff $\text{Re}(\lambda_i) < 0$ for both roots.

## E.3 Systems of ODEs

**Linear system:** $\dot{\mathbf{x}} = A\mathbf{x} + \mathbf{b}$.
- Steady state: $\mathbf{x}^* = -A^{-1}\mathbf{b}$ (if $A$ invertible).
- General solution: $\mathbf{x}(t) = e^{At}(\mathbf{x}_0 - \mathbf{x}^*) + \mathbf{x}^*$.
- Stable iff all eigenvalues of $A$ have negative real parts.
- Saddle point iff $\det(A) < 0$ (eigenvalues of mixed sign).

**Phase portrait classification** (2D system, eigenvalues $\lambda_1, \lambda_2$):

| Eigenvalue type | $\lambda_1, \lambda_2 < 0$ | $\lambda_1 < 0 < \lambda_2$ | $\lambda_1, \lambda_2 > 0$ | Complex $\alpha\pm i\beta$, $\alpha < 0$ |
|---|---|---|---|---|
| Type | Stable node | Saddle | Unstable node | Stable spiral |

## E.4 Pontryagin Conditions Summary

For $\max_{u(t)}\int_0^\infty e^{-\rho t}F(x,u)dt$ s.t. $\dot{x} = f(x,u)$:

| Condition | Formula | Interpretation |
|---|---|---|
| Current-value Hamiltonian | $\mathcal{H} = F(x,u) + \mu f(x,u)$ | utility + shadow value × rate of change |
| Optimality (FOC) | $\partial\mathcal{H}/\partial u = 0$ | Maximize $\mathcal{H}$ over $u$ |
| Costate equation | $\dot{\mu} = \rho\mu - \partial\mathcal{H}/\partial x$ | Dynamics of shadow price |
| State equation | $\dot{x} = \partial\mathcal{H}/\partial\mu = f(x,u)$ | Law of motion |
| Transversality | $\lim_{t\to\infty}e^{-\rho t}\mu(t)x(t) = 0$ | No asymptotic rents |

**Key costate interpretation:** $\mu(t) = \partial V/\partial x$ (shadow price of state = marginal value of relaxing the constraint).

---

# Appendix F: Difference Equations Cheat Sheet

---

## F.1 First-Order Difference Equations

**Linear:** $x_{t+1} = ax_t + b$.
- Steady state: $x^* = b/(1-a)$ ($a\neq1$).
- General solution: $x_t = a^t(x_0-x^*) + x^*$.
- Stable iff $|a| < 1$.

**Stability and behavior:**

| Condition | Behavior |
|---|---|
| $0 < a < 1$ | Monotone convergence to $x^*$ |
| $-1 < a < 0$ | Oscillatory (damped) convergence |
| $a > 1$ | Monotone divergence |
| $a < -1$ | Oscillatory divergence |

**Forward-looking equation:** $x_t = a\mathbb{E}_t[x_{t+1}] + b z_t$.
- If $|a|<1$: unique bounded solution $x_t = b\sum_{j=0}^\infty a^j\mathbb{E}_t[z_{t+j}]$.
- If $|a|>1$: MSV solution $x_t = b/(1-a\rho_z)z_t$ plus possible sunspots.

## F.2 Second-Order Difference Equations

**Linear:** $x_{t+2} + px_{t+1} + qx_t = c$.

**Characteristic equation:** $\lambda^2 + p\lambda + q = 0$, roots $\lambda_{1,2} = (-p\pm\sqrt{p^2-4q})/2$.

**Stability conditions:** $|a| < 1+q$ AND $q < 1$ (Schur–Cohn conditions — necessary and sufficient).

| Root type | Discriminant | Solution form |
|---|---|---|
| Real, distinct | $p^2 > 4q$ | $C_1\lambda_1^t + C_2\lambda_2^t + x^*$ |
| Real, equal | $p^2 = 4q$ | $(C_1+C_2t)\lambda^t + x^*$ |
| Complex | $p^2 < 4q$, $r=\sqrt{q}$, $\theta=\arccos(-p/(2r))$ | $r^t(C_1\cos\theta t+C_2\sin\theta t)+x^*$ |

**Oscillation period** (complex case): $T = 2\pi/\theta$ periods.

## F.3 Systems of Difference Equations

**Linear system:** $\mathbf{x}_{t+1} = A\mathbf{x}_t + \mathbf{b}$.
- Steady state: $\mathbf{x}^* = (I-A)^{-1}\mathbf{b}$.
- General solution: $\mathbf{x}_t = A^t(\mathbf{x}_0-\mathbf{x}^*)+\mathbf{x}^*$.
- Stable iff all eigenvalues of $A$ satisfy $|\lambda_i| < 1$.
- IRF at horizon $h$: $A^{h-1}\mathbf{b}$ (response to unit impulse).

**Companion form.** Convert $p$-th order scalar system to first-order system with companion matrix:

$$A = \begin{pmatrix}a_1 & a_2 & \cdots & a_p \\ 1 & 0 & \cdots & 0 \\ 0 & 1 & \cdots & 0 \\ \vdots & & \ddots & \vdots \end{pmatrix}.$$

Eigenvalues of $A$ are the roots of the characteristic polynomial.

## F.4 Blanchard–Kahn Counting Rule

For the linearized DSGE $\Gamma_0\mathbf{y}_t = \Gamma_1\mathbf{y}_{t-1}+\Psi\mathbf{z}_t+\Pi\bm\eta_t$:

| Unstable eigenvalues | Jump variables ($n_f$) | Outcome |
|---|---|---|
| $= n_f$ | Match | **Unique bounded solution** ✓ |
| $< n_f$ | Too few | **Indeterminate** (sunspots possible) |
| $> n_f$ | Too many | **No bounded solution** |

## F.5 MSV Solution Procedure Summary

For $\mathbf{y}_t = A\mathbb{E}_t[\mathbf{y}_{t+1}] + C\mathbf{z}_t$ with $\mathbf{z}_t = \Phi\mathbf{z}_{t-1} + \bm\varepsilon_t$:

1. **Conjecture:** $\mathbf{y}_t = \Omega\mathbf{z}_t$.
2. **Substitute:** $\Omega\mathbf{z}_t = A\Omega\Phi\mathbf{z}_t + C\mathbf{z}_t$.
3. **Sylvester equation:** $\Omega - A\Omega\Phi = C$.
4. **Vectorize:** $\text{vec}(\Omega) = (I-\Phi'\otimes A)^{-1}\text{vec}(C)$.
5. **Verify:** All eigenvalues of $A$ outside unit circle iff unique bounded solution.

---

# Appendix G: Time Series Concepts Summary

---

## G.1 Stationarity

**Strict stationarity.** $\{y_t\}$ is strictly stationary if the joint distribution of $(y_{t_1},\ldots,y_{t_k})$ equals that of $(y_{t_1+h},\ldots,y_{t_k+h})$ for all $k$, $h$, and time indices.

**Covariance (weak) stationarity.** $\mathbb{E}[y_t] = \mu$ (constant), $\text{Cov}(y_t, y_{t-j}) = \gamma_j$ (depends only on lag $j$).

**Key distinction:** Strict $\Rightarrow$ weak stationarity (for finite variance processes). Weak $\not\Rightarrow$ strict.

**Integrated processes.** $y_t \sim I(d)$: $d$-th difference is stationary. Most macro levels: $I(1)$; growth rates: $I(0)$.

## G.2 ARMA Representations

**AR($p$):** $y_t = \sum_{j=1}^p\phi_jy_{t-j} + \varepsilon_t$. Stationary iff all roots of $1-\phi_1z-\cdots-\phi_pz^p$ lie outside the unit circle.

**MA($q$):** $y_t = \sum_{j=0}^q\theta_j\varepsilon_{t-j}$. Always stationary.

**ARMA($p,q$):** Combines both. Autocovariance: $\gamma_j = \text{Cov}(y_t, y_{t-j})$ depends on $j$ only.

**Wold decomposition.** Any zero-mean $I(0)$ process: $y_t = \sum_{j=0}^\infty\psi_j\varepsilon_{t-j}$ (MA($\infty$)) with $\sum\psi_j^2<\infty$. Long-run multiplier: $\psi(1) = \sum\psi_j$.

## G.3 VAR Models

**VAR($p$):** $\mathbf{y}_t = \mathbf{c}+\sum_{j=1}^pA_j\mathbf{y}_{t-j}+\mathbf{e}_t$.

**OLS estimation:** $\hat{B} = (X'X)^{-1}X'Y$ — APL: `B ← (⌹ X) +.× Y`.

**Lag selection criteria:**

| Criterion | Formula |
|---|---|
| AIC | $\ln|\hat\Sigma_p| + 2n^2p/T$ |
| BIC | $\ln|\hat\Sigma_p| + n^2p\ln T/T$ |
| HQ | $\ln|\hat\Sigma_p| + 2n^2p\ln\ln T/T$ |

BIC consistent; AIC better for forecasting.

**Granger causality.** $x_t$ Granger-causes $y_t$ iff $A_{yx,j}\neq0$ for some lag $j$ — test via F-test on the restricted vs. unrestricted VAR equations.

## G.4 Structural VAR Identification

**Identification schemes** (additional restrictions beyond $\Sigma = B_0^{-1}(B_0^{-1})'$):

| Method | Restrictions | Advantage |
|---|---|---|
| Cholesky | Lower triangular $B_0$ | Simple; just-identified |
| Long-run | $A(1)$ matrix restrictions | Structural economic restrictions |
| Sign restrictions | IRF signs only | Robust; set-identified |
| External instrument | Proxy $z_t$ correlated with one shock | No ordering assumption |

**IRF:** Response of variable $i$ to shock $j$ at horizon $h$: $(A^h_{comp}\cdot B_0^{-1})_{ij}$.

**FEVD:** Fraction of $h$-step MSE of variable $i$ due to shock $j$: $\sum_{k=0}^{h-1}(\Psi_k B_0^{-1})_{ij}^2/\text{MSE}_i(h)$.

## G.5 Cointegration

**Engle–Granger two-step:**
1. OLS: $\hat\beta = (X'X)^{-1}X'Y$ (super-consistent, converges at rate $T$).
2. ADF test on residuals $\hat{u}_t = y_t - \hat\beta x_t$.

**Johansen trace statistic:** $\Lambda_{trace}(r) = -T\sum_{i=r+1}^n\ln(1-\hat\lambda_i)$.

**Error correction (VECM):** $\Delta\mathbf{y}_t = \alpha\beta'\mathbf{y}_{t-1}+\sum_{j=1}^{p-1}\Gamma_j\Delta\mathbf{y}_{t-j}+\bm\varepsilon_t$, where $\beta$ = cointegrating vectors, $\alpha$ = adjustment speeds.

## G.6 Key Test Statistics

| Test | Statistic | Distribution under $H_0$ | Purpose |
|---|---|---|---|
| ADF | $t_{\hat\gamma}$ in augmented regression | DF distribution (non-standard) | Unit root |
| Johansen trace | $-T\sum\ln(1-\hat\lambda_i)$ | Non-standard | Cointegration rank |
| Diebold–Mariano | $\bar{d}/\sqrt{\hat V/N}$ | $\mathcal{N}(0,1)$ | Forecast accuracy |
| Ljung–Box | $T(T+2)\sum_{j=1}^m\hat\rho_j^2/(T-j)$ | $\chi^2(m)$ | Serial correlation |
| Granger causality | $F$-statistic on excluded lags | $F(p,T-2p-1)$ | Granger causality |
