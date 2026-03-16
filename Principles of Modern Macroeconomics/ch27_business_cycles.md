# Chapter 27 — Business Cycles: Booms, Recessions, and Stabilization

---

## 27.1 Empirical Regularities

Kaldor–Prescott stylized facts in HP-filtered U.S. quarterly data (1960–2019): $\sigma_y \approx 1.5\%$, $\sigma_c \approx 0.9\%$, $\sigma_i \approx 5.5\%$. Correlations: $\mathrm{corr}(y,c) \approx 0.84$, $\mathrm{corr}(y,i) \approx 0.92$. Investment is two to three times more volatile than output; consumption is less volatile; unemployment is countercyclical and lags output.

---

## 27.2 The Real Business Cycle Model

Kydland and Prescott (1982): a calibrated neoclassical model with AR(1) technology shocks can match these second moments. Stochastic technology:

$$Y_t = A_t F(K_t, N_t), \quad \ln A_t = \rho_A \ln A_{t-1} + \epsilon_t^A,\; \epsilon_t^A \sim \mathcal{N}(0,\sigma_A^2).$$

Household utility:

$$\mathbb{E}_0\sum_{t=0}^\infty \beta^t\!\left[\frac{c_t^{1-\sigma}-1}{1-\sigma} + \chi\frac{(1-n_t)^{1-\eta}-1}{1-\eta}\right].$$

With $\rho_A = 0.95$, $\sigma_A = 0.7\%$ and standard parameter values, the model generates $\sigma_y \approx 1.4\%$, $\mathrm{corr}(y,n) \approx 0.95$ — close to the data.

---

## 27.3 Criticisms and the New Keynesian Response

RBC criticisms: (i) the required TFP shock variance is large; (ii) the Frisch elasticity must exceed 2 to match employment fluctuations — far above microeconomic estimates; (iii) the model implies wages are procyclical and flexible. The New Keynesian model adds monopolistic competition and Calvo price stickiness to match business cycle facts with smaller shocks and a role for monetary policy.

---

## 27.4 Structural VAR Identification

A $k$-variable SVAR:

$$\mathbf{y}_t = \mathbf{c} + \sum_{j=1}^p A_j \mathbf{y}_{t-j} + \mathbf{u}_t, \quad \mathbf{u}_t = B\bm{\varepsilon}_t,\; \bm{\varepsilon}_t \sim \mathcal{N}(\mathbf{0}, I).$$

The $k^2$ elements of $B$ are not identified from $\Sigma = BB'$ ($k(k+1)/2$ equations). Identification strategies: Cholesky ordering, sign restrictions, Romer–Romer narrative dates, and external instruments (Stock and Watson, 2012).

---

*Next: Chapter 28 — Fiscal Policy in Practice*
