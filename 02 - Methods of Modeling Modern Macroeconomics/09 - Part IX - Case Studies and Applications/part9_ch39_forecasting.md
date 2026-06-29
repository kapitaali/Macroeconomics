# Chapter 39: Forecasting Inflation and GDP

*Time Series Models vs. DSGE: A Forecasting Horse Race*

> *"All models are wrong, but some are useful. The right question is not which model is true, but which model forecasts better — and under what conditions."*
> — George Box (paraphrased)

**Cross-reference:** *Principles* Ch. 39 (future of macroeconomics: forecasting, ML, model uncertainty); Appendix B (DSGE forecasting, real-time data) **[P:Ch.39, P:AppB]**

---

## 39.1 The Forecasting Problem

Macroeconomic forecasting — predicting GDP growth, inflation, and unemployment — is one of the most practically important applications of the methods in this book. Central banks use forecasts to set policy rates; firms use them to make investment decisions; governments use them to project revenues and plan spending.

**The formal forecasting problem:** Given observations $\mathbf{Y}_{1:T} = \{\mathbf{y}_1, \ldots, \mathbf{y}_T\}$ up to time $T$, produce forecasts $\hat{\mathbf{y}}_{T+h|T}$ for horizons $h = 1, 2, \ldots, H$.

**Competing approaches:**
1. **Reduced-form time series:** BVAR with Minnesota prior, time-varying parameter VAR, factor models.
2. **DSGE models:** The Kalman filter prediction step from Chapter 20 generates $h$-step-ahead forecasts.
3. **Machine learning:** LASSO, Ridge, random forests, neural networks applied to large datasets (FRED-MD: 130 monthly series).

This chapter conducts a horse race: all three approaches on the same data, evaluated by the same metrics.

---

## 39.2 Forecast Evaluation Metrics

**Definition 39.1 (Root Mean Squared Error).** For a sequence of $h$-step-ahead forecasts $\hat{y}_{T+h|T}$ evaluated over $T = T_0, \ldots, T_1$:

$$RMSE_h = \sqrt{\frac{1}{T_1-T_0+1}\sum_{T=T_0}^{T_1}(y_{T+h} - \hat{y}_{T+h|T})^2}.$$

**Definition 39.2 (Continuous Ranked Probability Score).** For density forecasts $\hat{F}_{T+h|T}(\cdot)$:

$$CRPS = \int_{-\infty}^\infty[\hat{F}_{T+h|T}(z) - \mathbf{1}(y_{T+h} \leq z)]^2\,dz.$$

The CRPS reduces to RMSE for point forecasts. It penalizes both poor location (bias) and poor spread (overconfidence or underconfidence). Lower is better.

---

## 39.3 BVAR with Minnesota Prior

**Definition 39.3 (Minnesota Prior).** The **Minnesota prior** (Doan, Litterman, and Sims, 1984) for a VAR($p$) with $n$ variables imposes:

- Own lags: coefficient on lag $j$ of variable $i$ in equation $i$ has prior $\mathcal{N}(\delta_i, \sigma^2_{ii}/j^2)$, where $\delta_i = 1$ for I(1) variables (random walk prior) and $\delta_i = 0$ for I(0).
- Cross-variable lags: coefficient on lag $j$ of variable $l \neq i$ in equation $i$ has prior $\mathcal{N}(0, \lambda^2\sigma^2_{ii}/(j^2\sigma^2_{ll}))$, where $\lambda$ is the overall tightness.
- Intercepts: loose prior centered at zero.

The Minnesota prior shrinks coefficients toward the random walk ($\delta_i = 1$) or zero ($\delta_i = 0$), with shrinkage increasing at longer lags. This regularizes the VAR for forecasting.

**Theorem 39.1 (Minnesota Prior as Augmented Regression).** The Minnesota prior posterior for the VAR coefficients $B$ is equivalent to OLS on augmented data $(\mathbf{Y}_{aug}, \mathbf{X}_{aug})$ where dummy observations encode the prior beliefs:

$$\hat{B}_{Minnesota} = (X_{aug}'X_{aug})^{-1}X_{aug}'Y_{aug}.$$

In APL: `B_MN ← (⌹ X_aug) +.× Y_aug` — one expression for the BVAR estimator.

*Proof.* The Minnesota prior is conjugate (Normal–Wishart) for the VAR with fixed covariance. The posterior mode equals the OLS estimator on the augmented system where the dummy observations implement the prior beliefs as data. Kadiyala and Karlsson (1997) show this equivalence formally. $\square$

**Constructing the dummy observations:** For the standard Minnesota prior with hyperparameters $(\lambda, \delta, \sigma^2_i)$:

```python
import numpy as np

def minnesota_dummies(Y, p, lambda_tight=0.2, delta=1.0):
    """
    Generate Minnesota prior dummy observations.
    Y: T×n data matrix, p: lag order, lambda_tight: shrinkage, delta: prior mean on own lags.
    Returns Y_dum, X_dum to be prepended to actual data for BVAR.
    """
    T, n = Y.shape
    sigma = np.std(Y, axis=0)  # marginal standard deviations
    
    dummies = []
    # Own-lag dummies: one dummy per variable per lag
    for j in range(1, p+1):
        for i in range(n):
            y_d = np.zeros(n); x_d = np.zeros(n*p + 1)
            y_d[i] = sigma[i] * j / lambda_tight
            x_d[i + (j-1)*n] = sigma[i] * j / lambda_tight
            dummies.append((y_d, x_d))
    
    # Sum-of-coefficients dummies
    y_d = delta * np.mean(Y[:5], axis=0)
    x_d = np.zeros(n*p+1)
    for j in range(p): x_d[j*n:(j+1)*n] = delta * np.mean(Y[:5], axis=0)
    dummies.append((y_d, x_d))
    
    Y_dum = np.array([d[0] for d in dummies])
    X_dum = np.array([d[1] for d in dummies])
    return Y_dum, X_dum
```

---

## 39.4 DSGE Forecasting via Kalman Filter

From the gensys solution (Chapter 28), the NK-DSGE has the state-space form:

$$\bm\alpha_{t+1} = F\bm\alpha_t + G\bm\varepsilon_{t+1}, \qquad \mathbf{y}_t = H\bm\alpha_t + \mathbf{e}_t,$$

with $\bm\varepsilon_t \sim \mathcal{N}(0, Q)$ and $\mathbf{e}_t \sim \mathcal{N}(0, R)$.

The **$h$-step-ahead forecast** from the Kalman filter prediction step:

$$\hat{\mathbf{y}}_{T+h|T} = H F^h \hat{\bm\alpha}_{T|T},$$

$$\text{MSE}(\hat{\mathbf{y}}_{T+h|T}) = H F^h P_{T|T} (F^h)' H' + \sum_{j=1}^h HF^{j-1}GQG'(F^{j-1})'H' + R.$$

The DSGE forecasts are computed in real time: at each evaluation date $T$, the state $\hat{\bm\alpha}_{T|T}$ is updated by the Kalman filter using all available data, then projected forward $h$ steps.

---

## 39.5 Machine Learning: LASSO for Macro Forecasting

With $K = 130$ predictor variables (FRED-MD), standard OLS has an $n/K$ degrees-of-freedom problem. **LASSO** (Tibshirani, 1996) adds an $\ell_1$ penalty to shrink many coefficients to exactly zero:

$$\hat{\bm\beta}^{LASSO} = \arg\min_{\bm\beta}\;\frac{1}{T}\|\mathbf{y} - \mathbf{X}\bm\beta\|_2^2 + \lambda_{LASSO}\|\bm\beta\|_1.$$

**Theorem 39.2 (LASSO Solution via Coordinate Descent).** The LASSO objective has a component-wise closed-form update. For each coordinate $j$:

$$\hat\beta_j \leftarrow \mathcal{S}\!\left(\frac{1}{T}\mathbf{x}_j'(\mathbf{y} - \mathbf{X}_{-j}\hat{\bm\beta}_{-j}),\; \lambda_{LASSO}\right),$$

where $\mathcal{S}(z, \lambda) = \text{sign}(z)\max(|z|-\lambda, 0)$ is the **soft-thresholding operator**. Cycling through all $j$ until convergence gives the LASSO solution.

*Proof.* The $j$-th component subproblem (with all other coefficients fixed) is $\min_{\beta_j}\frac{1}{T}(r_j - x_j\beta_j)^2 + \lambda|\beta_j|$, where $r_j = y - X_{-j}\hat\beta_{-j}$ is the partial residual. This is the Lasso regression of $r_j$ on $x_j$, with solution $\hat\beta_j = \mathcal{S}(x_j'r_j/T, \lambda)$ by the KKT conditions for the $|\beta_j|$ term. $\square$

In APL, the soft-threshold operator and coordinate descent step:

```apl
⍝ APL — LASSO coordinate descent
⎕IO←0 ⋄ ⎕ML←1

soft_threshold ← {z lam ← ⍵
    (×z) × 0⌈(|z)-lam}    ⍝ sign(z)*max(|z|-λ, 0)

⍝ One full coordinate descent pass over all K predictors
cd_pass ← {X y lam beta ← ⍵
    :For j :In ⍳ ≢beta
        r_j  ← y - (X[;~j]) +.× beta[~j]  ⍝ partial residual
        z_j  ← (X[;j] +.× r_j) ÷ ≢y       ⍝ univariate OLS coefficient
        beta[j] ← soft_threshold z_j lam    ⍝ soft-threshold
    :EndFor
    beta}

⍝ LASSO: iterate until convergence
lasso ← {X y lam ← ⍵
    K ← (⍴X)[1]
    beta0 ← K ⍴ 0
    cd_pass∘(X y lam) ⍣ (1e¯6∘>⌈/|⊢-cd_pass∘(X y lam)) ⊢ beta0}
```

---

## 39.6 Forecast Combination

No single model dominates across all variables, horizons, and sample periods. **Forecast combination** pools multiple models to reduce forecast risk.

**Simple average:** $\hat{y}^{pool}_{T+h|T} = \frac{1}{M}\sum_{m=1}^M\hat{y}^m_{T+h|T}$.

**Optimal linear pool:** Minimize $\sum_T(y_{T+h}-\sum_m w_m\hat{y}^m)^2$ s.t. $w_m \geq 0$, $\sum w_m = 1$ — a constrained OLS problem.

**Bayesian Model Averaging (BMA):** Weight each model by its marginal likelihood (posterior probability). For $M$ models: $w_m \propto p(\mathbf{Y}|\mathcal{M}_m)\cdot p(\mathcal{M}_m)$.

---

## 39.7 The Diebold–Mariano Test

**Definition 39.4 (Diebold–Mariano Test).** The **DM test** compares the predictive accuracy of two models, $m = 1$ and $m = 2$. Define the differential loss $d_T = L(e^1_T) - L(e^2_T)$, where $e^m_T = y_{T+h} - \hat{y}^m_{T+h|T}$ and $L(\cdot)$ is the loss function (e.g., $L(e) = e^2$ for MSE). The test statistic:

$$DM = \frac{\bar{d}}{\sqrt{\hat{V}(\bar{d})/N}} \xrightarrow{d} \mathcal{N}(0,1) \quad \text{under } H_0: \mathbb{E}[d_T] = 0.$$

where $\bar{d} = N^{-1}\sum_T d_T$ and $\hat{V}(\bar{d})$ is the Newey–West HAC variance estimator.

$H_0$: equal predictive accuracy. $H_1$: one model is significantly better.

---

## 39.8 Worked Example: U.S. Inflation and GDP Forecasting

**Data:** U.S. quarterly GDP growth and CPI inflation, 1960Q1–2019Q4. Pseudo-real-time evaluation: estimate on rolling 20-year windows, forecast 1–8 quarters ahead.

```python
import numpy as np
from sklearn.linear_model import LassoCV
from sklearn.preprocessing import StandardScaler

np.random.seed(42)
T_total = 200   # Simulated quarterly data
T_train = 80    # Initial training window

# Simulate macro data (AR1 for illustration)
y_gdp   = np.zeros(T_total); y_gdp[0] = 2.0
y_infl  = np.zeros(T_total); y_infl[0] = 2.0
for t in range(1, T_total):
    y_gdp[t]  = 0.5*y_gdp[t-1]  + np.random.normal(0, 1.5)
    y_infl[t] = 0.6*y_infl[t-1] + 0.2*y_gdp[t-1] + np.random.normal(0, 0.5)

# Forecasting horse race: BVAR vs. LASSO vs. naive AR
H_horizon = 4  # 4-step ahead (1-year)
n_eval = T_total - T_train - H_horizon

rmse = {'AR1': [], 'BVAR': [], 'LASSO': []}

for T in range(T_train, T_total - H_horizon):
    y_train = y_gdp[:T]; y_target = y_gdp[T + H_horizon]
    
    # AR(1) forecast
    ar_coef = np.polyfit(y_train[:-1], y_train[1:], 1)
    forecast_ar = y_train[-1]
    for _ in range(H_horizon): forecast_ar = ar_coef[0]*forecast_ar + ar_coef[1]
    rmse['AR1'].append((y_target - forecast_ar)**2)
    
    # BVAR (simplified: VAR(2) with Minnesota prior via augmented OLS)
    p_bvar = 2; n_vars = 2
    Y = np.column_stack([y_gdp[:T], y_infl[:T]])
    X_bvar = np.column_stack([Y[:-p_bvar], Y[1:T-p_bvar+1], np.ones(T-p_bvar)])
    Y_dep  = Y[p_bvar:]
    
    # Add Minnesota prior dummies (simplified)
    sigma_vec = np.std(Y, axis=0)
    lam = 0.2
    Y_dum_rows = []; X_dum_rows = []
    for j in range(1, p_bvar+1):
        for i in range(n_vars):
            y_d = np.zeros(n_vars); y_d[i] = sigma_vec[i]*j/lam
            x_d = np.zeros(n_vars*p_bvar+1); x_d[i+(j-1)*n_vars] = sigma_vec[i]*j/lam
            Y_dum_rows.append(y_d); X_dum_rows.append(x_d)
    
    Y_aug = np.vstack([np.array(Y_dum_rows), Y_dep])
    X_aug = np.vstack([np.array(X_dum_rows), X_bvar])
    
    B_bvar = np.linalg.lstsq(X_aug, Y_aug, rcond=None)[0]
    state = np.concatenate([Y[-1], Y[-2], [1]])
    for _ in range(H_horizon): state = np.concatenate([B_bvar.T @ state, state[:n_vars], [1]])
    forecast_bvar = state[0]  # GDP forecast
    rmse['BVAR'].append((y_target - forecast_bvar)**2)
    
    # LASSO (use lagged variables as predictors)
    X_lasso = np.column_stack([y_train[j:T-H_horizon+j] for j in range(8)])
    y_lasso = y_train[H_horizon:]
    if len(y_lasso) > 20:
        scaler = StandardScaler()
        X_sc = scaler.fit_transform(X_lasso)
        lasso = LassoCV(cv=5, max_iter=1000).fit(X_sc, y_lasso)
        x_new = scaler.transform(y_train[-8:].reshape(1,-1))
        forecast_lasso = float(lasso.predict(x_new))
    else:
        forecast_lasso = y_train[-1]
    rmse['LASSO'].append((y_target - forecast_lasso)**2)

# Compute RMSE
for model, losses in rmse.items():
    print(f"{model}: RMSE = {np.sqrt(np.mean(losses)):.4f}")

# DM test: BVAR vs. AR1
d = np.array(rmse['AR1']) - np.array(rmse['BVAR'])
d_bar = np.mean(d); n = len(d)
# Newey-West variance with 4 lags
V_nw = np.var(d)/n
for lag in range(1, 5):
    V_nw += 2*(1-lag/5)*np.cov(d[lag:], d[:-lag])[0,1]/n
DM_stat = d_bar / np.sqrt(V_nw)
from scipy.stats import norm
p_value = 2*(1-norm.cdf(abs(DM_stat)))
print(f"\nDM test (BVAR vs. AR1): DM = {DM_stat:.3f}, p-value = {p_value:.4f}")
print(f"{'BVAR significantly better' if p_value < 0.05 else 'No significant difference'} at 5% level")
```

```julia
using Statistics, GLMNet

# Julia: LASSO for macro forecasting
function lasso_forecast(y, T_train, H)
    X_lags = hcat([y[j:T_train-H+j-1] for j in 1:8]...)
    y_target = y[H+1:T_train]
    path = glmnet(X_lags, y_target)
    # Cross-validate to select lambda
    cv = glmnetcv(X_lags, y_target)
    lambda_opt = cv.lambda[argmin(cv.meanloss)]
    # Predict
    x_new = reshape(y[T_train-7:T_train], 1, 8)
    pred = GLMNet.predict(path, x_new, outtype=:response)
    idx = argmin(abs.(path.lambda .- lambda_opt))
    return pred[1, idx]
end

println("Macro forecasting in Julia: use GLMNet.jl for LASSO")
println("Key finding: LASSO often beats AR(1) at 4-8 quarter horizon due to variable selection")
println("BVAR beats LASSO at short horizons due to Minnesota prior regularization")
println("DSGE beats both when model is correctly specified but loses during crises")
```

---

## 39.9 Programming Exercises

### Exercise 39.1 (APL — LASSO Coordinate Descent)

Implement the LASSO coordinate descent from Section 39.5 in APL: (a) `soft_threshold ← {(×⍺)×0⌈(|⍺)-⍵}` as a dfn taking $(z, \lambda)$; (b) one full coordinate descent pass over $K$ predictors; (c) iterate until convergence using `⍣≡`; (d) test on simulated data with $T = 100$, $K = 50$ predictors, half of which are relevant. Verify LASSO selects approximately the right 25 predictors.

### Exercise 39.2 (Python — Forecasting Horse Race)

Run the full horse race on FRED-MD data (or simulated data with the same structure): (a) three models: AR(4), BVAR(2) with Minnesota prior, LASSO with 8 lags of 20 macro variables; (b) evaluation period: 1990Q1–2019Q4, rolling 40-quarter estimation windows; (c) target: 4-quarter GDP growth and 4-quarter inflation; (d) report RMSE ratios relative to AR(4) and DM test p-values.

### Exercise 39.3 (Julia — DSGE Real-Time Forecasting)

```julia
# Real-time DSGE forecasting via Kalman filter
using LinearAlgebra

function dsge_forecast(A, D, F, H_obs, Q, R, y_obs, h)
    """h-step-ahead forecast from DSGE Kalman filter."""
    T, p = size(y_obs,1), size(H_obs,1)
    m = size(A, 1)
    
    # Filter (Chapter 20)
    alpha = zeros(m); P = 10*I(m)
    for t in 1:T
        # Predict
        alpha_pred = A * alpha
        P_pred = A * P * A' + D * Q * D'
        # Update
        v = y_obs[t,:] - H_obs * alpha_pred
        Fv = H_obs * P_pred * H_obs' + R
        K = P_pred * H_obs' * inv(Fv)
        alpha = alpha_pred + K * v
        P = (I(m) - K*H_obs) * P_pred
    end
    
    # h-step forecast: alpha_{T+h|T} = A^h * alpha_{T|T}
    alpha_fcast = A^h * alpha
    return H_obs * alpha_fcast
end

println("DSGE h-step forecast: H_obs * A^h * alpha_filtered")
println("Error grows with h due to uncertainty accumulation:")
for h in [1, 4, 8]
    uncertainty_factor = h^0.5  # approx: error ~ sqrt(h) for stable systems
    println("  h=$h: relative uncertainty ≈ $(round(uncertainty_factor, digits=2))x baseline")
end
```

### Exercise 39.4 — Model Uncertainty ($\star$)

Implement Bayesian Model Averaging (BMA) for the three-model forecast combination: (a) compute the log marginal likelihood $\ln p(\mathbf{Y}|\mathcal{M}_m)$ for each model using the Laplace approximation; (b) compute BMA weights $w_m \propto p(\mathbf{Y}|\mathcal{M}_m)p(\mathcal{M}_m)$ with equal model priors; (c) compare BMA weights to the ex-post optimal combination weights (OLS on the forecasts). Do the BMA weights select the model that actually forecast best?

---

## 39.10 Chapter Summary

**Key results:**

- **BVAR Minnesota prior**: prior $\mathcal{N}(\delta_i,\sigma^2_{ii}/j^2)$ on own lags, $\mathcal{N}(0,\lambda^2\sigma^2_{ii}/(j^2\sigma^2_{ll}))$ on cross-lags; equivalent to OLS on augmented data (Theorem 39.1), computed as `B_MN ← (⌹ X_aug) +.× Y_aug`.
- **DSGE forecasting**: $h$-step forecast $\hat{\mathbf{y}}_{T+h|T} = HF^h\hat{\bm\alpha}_{T|T}$ from the Kalman filter; valid iff the model is correctly specified (fails during structural breaks).
- **LASSO coordinate descent** (Theorem 39.2): soft-threshold update $\hat\beta_j = \mathcal{S}(x_j'r_j/T, \lambda)$; in APL: `{soft_threshold (X[;j]+.×partial_residual÷T) lam}⍣≡`.
- **DM test** (Theorem implicit): $DM = \bar{d}/\sqrt{\hat{V}(\bar{d})/N} \to \mathcal{N}(0,1)$; tests equal predictive accuracy; uses Newey–West variance.
- **Findings from the literature**: BVAR dominates AR at short horizons (1–2 quarters); LASSO dominates BVAR when many predictors are relevant (large-data settings); DSGE dominates during normal times but fails during crises; combination always helps.

*Next: Chapter 40 — Policy Analysis with a New Keynesian Model*
