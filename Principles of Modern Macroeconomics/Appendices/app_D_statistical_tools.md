# Appendix D — Statistical Tools and Data Analysis

This appendix introduces the statistical and econometric tools most frequently encountered in macroeconomic research. It presupposes familiarity with basic probability theory (random variables, expectations, variances, conditional distributions) and elementary linear regression (OLS, $t$-statistics, $R^2$). The goal is to develop the tools needed to read empirical macroeconomics papers and to interpret quantitative claims critically.

---

## D.1 Time-Series Fundamentals

### Stationarity and Autocovariance

A time series $\{x_t\}_{t=-\infty}^\infty$ is **covariance-stationary** if its mean, variance, and autocovariances are all finite and time-invariant:

$$\mathbb{E}[x_t] = \mu \quad \text{for all } t$$
$$\mathrm{Var}(x_t) = \sigma^2 < \infty \quad \text{for all } t$$
$$\mathrm{Cov}(x_t, x_{t-k}) = \gamma_k \quad \text{depends only on } k, \text{ not } t.$$

The **autocorrelation function (ACF)** is $\rho_k = \gamma_k / \gamma_0$. The **partial autocorrelation function (PACF)** is the correlation between $x_t$ and $x_{t-k}$ after removing the linear dependence explained by $x_{t-1}, \ldots, x_{t-k+1}$.

An AR(1) process $x_t = \phi x_{t-1} + \epsilon_t$, $|\phi| < 1$, $\epsilon_t \sim \text{WN}(0,\sigma^2)$ is stationary with:
$$\mu = 0, \quad \sigma^2_x = \sigma^2/(1-\phi^2), \quad \gamma_k = \phi^k\sigma^2/(1-\phi^2), \quad \rho_k = \phi^k.$$

### Unit Root Processes

A process has a **unit root** if $x_t = x_{t-1} + \epsilon_t$ (random walk) or more generally $\phi(L)x_t = \epsilon_t$ where $\phi(z)$ has a root on the unit circle. Unit-root processes are non-stationary: their variance grows without bound and shocks have permanent effects. The **Augmented Dickey–Fuller (ADF)** test for $H_0: \phi = 1$ against $H_1: |\phi| < 1$ estimates:

$$\Delta x_t = \mu + \gamma t + (\rho - 1)x_{t-1} + \sum_{j=1}^p c_j \Delta x_{t-j} + u_t$$

and tests $H_0: \rho = 1$ using a non-standard distribution (critical values from Dickey and Fuller, 1979). The lags $p$ control for serial correlation in $u_t$ and are chosen by AIC or BIC.

The **Phillips–Perron (PP)** test is a non-parametric alternative that corrects for serial correlation in $u_t$ without adding lags. Both tests have low power in small samples — they often fail to reject the unit root null even when the true process is stationary with a root close to one.

### Cointegration

Two or more unit-root series are **cointegrated** if a linear combination of them is stationary. Cointegration implies a long-run stable relationship: the series share a common stochastic trend and deviations from the long-run relationship are transitory.

If $y_t \sim I(1)$ (integrated of order 1) and $x_t \sim I(1)$ and $y_t - \beta x_t \sim I(0)$ (stationary), then $y_t$ and $x_t$ are cointegrated with cointegrating vector $(1, -\beta)$. Economically, the money demand relation $\ln(M/P) = a + b\ln Y - ci$ may be a cointegrating relationship: each component is non-stationary, but their specific linear combination is stationary, representing a stable long-run equilibrium.

The **Engle–Granger two-step** procedure: (1) regress $y_t$ on $x_t$ by OLS to obtain $\hat{\beta}$; (2) test whether the residuals $\hat{u}_t = y_t - \hat{\beta}x_t$ are stationary using ADF. The **Johansen procedure** tests for cointegration in a multivariate system using maximum likelihood and identifies all cointegrating vectors.

---

## D.2 Regression with Time-Series Data

### OLS with Serially Correlated Errors

When the error term $u_t$ in $y_t = \mathbf{x}_t'\boldsymbol{\beta} + u_t$ is serially correlated — which is the rule rather than the exception in time-series data — OLS remains unbiased and consistent (under stationarity and ergodicity) but its standard errors are invalid. The **Newey–West** estimator provides heteroskedasticity and autocorrelation consistent (HAC) standard errors:

$$\hat{V}_{NW} = (X'X)^{-1}\hat{S}(X'X)^{-1},$$

where $\hat{S} = \hat{\Gamma}_0 + \sum_{j=1}^q w_j(\hat{\Gamma}_j + \hat{\Gamma}_j')$ and $w_j = 1 - j/(q+1)$ (Bartlett weights). The bandwidth $q$ is typically set to $\lfloor 4(T/100)^{2/9} \rfloor$ or similar data-dependent rule.

### Instrumental Variables in Time Series

The IV estimator $\hat{\boldsymbol{\beta}}_{IV} = (Z'X)^{-1}Z'y$ requires: (i) the instrument $Z$ is correlated with the endogenous regressor $X$ (relevance: first stage $F > 10$); (ii) $Z$ is uncorrelated with the error $u$ (exclusion). In time series, finding valid instruments is challenging because all lagged variables are potentially endogenous. Narrative instruments (Romer–Romer monetary policy shocks, Ramey defense news) are constructed to be exogenous from first principles.

### Regression Discontinuity in Macroeconomics

The sharp RD estimator identifies the local average treatment effect at the threshold $c$:

$$\tau_{RD} = \lim_{x\downarrow c}\mathbb{E}[y | x] - \lim_{x\uparrow c}\mathbb{E}[y | x],$$

estimated by local polynomial regression in a bandwidth $h$ around $c$. Bandwidth selection balances bias (choosing too wide a window includes observations far from the threshold where linearity may fail) against variance (too narrow a window has few observations). The optimal MSE bandwidth (Imbens and Kalyanaraman, 2012) is widely used.

---

## D.3 Spectral Analysis

Spectral analysis decomposes the variance of a time series across frequencies, identifying which periodicities are most important. The **spectral density** of a stationary process is:

$$S_x(\omega) = \frac{1}{2\pi}\sum_{k=-\infty}^\infty \gamma_k e^{-ik\omega}, \quad \omega \in [-\pi, \pi].$$

The spectral density $S_x(\omega)$ gives the contribution to the total variance of components with frequency $\omega$ (period $2\pi/\omega$). For an AR(1) with $\phi > 0$: the spectrum is highest at low frequencies (long cycles), reflecting the positive autocorrelation. For white noise: the spectrum is flat at $\sigma^2/(2\pi)$ — all frequencies contribute equally.

In macroeconomics, spectral analysis identifies the "business cycle frequencies": Hodrick and Prescott (1997) define the business cycle as fluctuations with periods of 6–32 quarters (frequencies corresponding to cycles of 1.5–8 years), motivating the HP filter which removes components outside this range.

The **band-pass filter** (Baxter and King, 1999) directly extracts components within a specified frequency band, using a moving average filter designed in the frequency domain. It is more flexible than the HP filter and avoids the endpoint bias of the two-sided HP filter.

---

## D.4 Decomposing Trend and Cycle

Beyond the HP filter (developed in Chapter 6), several approaches to decomposing macroeconomic series into trend and cyclical components have been proposed.

### Beveridge–Nelson Decomposition

The **Beveridge–Nelson (BN) decomposition** (Beveridge and Nelson, 1981) decomposes an integrated series into a random walk trend and a stationary cyclical component:

$$x_t = \tau_t + c_t,$$

where $\tau_t$ is the "permanent component" — the value the series would take in the long run if the disturbance process ended immediately:

$$\tau_t = \lim_{h\to\infty}(\mathbb{E}_t[x_{t+h}] - h\cdot\hat{\mu}),$$

and $c_t = x_t - \tau_t$ is the transitory "cyclical" component. Under the BN decomposition, all cycles are contemporaneously correlated with trend innovations (because they share the same shock structure), which gives very smooth cycles and jagged trends — the opposite of what one might naively expect.

### Unobserved Components Models

An **unobserved components (UC) model** treats the trend and cycle as separate stochastic processes and estimates both simultaneously by the Kalman filter:

$$x_t = \tau_t + c_t$$
$$\tau_t = \mu + \tau_{t-1} + \eta_t, \quad \eta_t \sim \mathcal{N}(0, \sigma_\eta^2)$$
$$c_t = \phi_1 c_{t-1} + \phi_2 c_{t-2} + \epsilon_t, \quad \epsilon_t \sim \mathcal{N}(0, \sigma_\epsilon^2).$$

The signal-to-noise ratio $q = \sigma_\eta^2/\sigma_\epsilon^2$ governs the smoothness of the trend: small $q$ implies a smooth trend (the HP filter with large $\lambda$ is a special case), while large $q$ implies a volatile trend. Watson (1986) and Clark (1987) estimated UC models for U.S. real GDP, finding that the business cycle component explains a substantial share of output variance.

---

## D.5 Panel Data Econometrics

A **balanced panel** consists of observations on $N$ cross-sectional units (countries, firms, households) over $T$ time periods. The general linear panel model:

$$y_{it} = \alpha_i + \lambda_t + \mathbf{x}_{it}'\boldsymbol{\beta} + \epsilon_{it}.$$

**Fixed effects (FE) estimation** treats $\alpha_i$ as free parameters to be estimated (by including $N-1$ dummy variables or by within-transformation). FE controls for all time-invariant heterogeneity across units, so identification comes from within-unit variation over time. The FE estimator is consistent as $T\to\infty$ for fixed $N$ but suffers from **Nickell bias** when $T$ is small: including lags of $y_{it}$ as regressors induces downward bias of order $1/T$.

**Random effects (RE) estimation** treats $\alpha_i$ as drawn from a distribution and estimates $\boldsymbol{\beta}$ by generalized least squares. RE is more efficient than FE if $\alpha_i$ is uncorrelated with $\mathbf{x}_{it}$, but inconsistent if not. The **Hausman test** distinguishes FE and RE: under $H_0$ (no correlation between $\alpha_i$ and $\mathbf{x}_{it}$), both are consistent but RE is efficient; under $H_1$, FE is consistent but RE is not. The Hausman statistic:

$$H = (\hat{\boldsymbol{\beta}}_{FE} - \hat{\boldsymbol{\beta}}_{RE})'[\hat{V}_{FE} - \hat{V}_{RE}]^{-1}(\hat{\boldsymbol{\beta}}_{FE} - \hat{\boldsymbol{\beta}}_{RE}) \sim \chi^2_k.$$

**Two-way fixed effects (TWFE)** includes both unit fixed effects $\alpha_i$ and time fixed effects $\lambda_t$. With recent treatments (staggered treatment timing), the TWFE estimator can produce sign-reversed estimates due to negative weighting of late-treated units (Callaway and Sant'Anna, 2021; Goodman-Bacon, 2021). Modern DiD methods address this by estimating group-time average treatment effects and aggregating appropriately.

---

## D.6 Measuring Uncertainty and Risk

Macroeconomic models increasingly incorporate uncertainty as a state variable, requiring empirical measures. Three approaches are standard.

**Cross-sectional dispersion** measures disagreement among forecasters or across firms. The standard deviation of forecasts in the Survey of Professional Forecasters (SPF) or the Federal Reserve's Greenbook is used as a proxy for aggregate uncertainty. Higher dispersion suggests that agents disagree about the state of the economy, consistent with higher effective uncertainty.

**Realized volatility** measures the ex-post variance of a variable computed from high-frequency data. The VIX (CBOE Volatility Index) is the market-implied standard deviation of S&P 500 returns over the next 30 days, derived from options prices. Bloom (2009) uses the VIX as a proxy for macroeconomic uncertainty in his study of investment responses to uncertainty shocks.

**Text-based measures** use natural language processing applied to newspapers, central bank communications, or firm earnings calls to construct indices of economic policy uncertainty. Baker, Bloom, and Davis (2016) construct the **Economic Policy Uncertainty (EPU)** index from: the frequency of newspaper articles mentioning economic uncertainty and policy; the number of expiring tax provisions; and the dispersion of economic forecasts. The EPU index has become a widely-used empirical proxy for policy uncertainty.

---

## D.7 Bayesian Statistics for Macroeconomics

Bayesian inference updates a prior belief $p(\theta)$ about parameters using observed data $y$ to obtain a posterior belief $p(\theta|y) \propto p(y|\theta) \cdot p(\theta)$.

**Prior distributions** for macroeconomic parameters are typically chosen to: (i) restrict parameters to theoretically plausible ranges (e.g., $\beta \in (0,1)$: Beta distribution; $\sigma > 0$: Inverse-Gamma distribution); (ii) center near values from microeconomic studies; (iii) be relatively diffuse to let the data speak.

The **Metropolis–Hastings (MH) algorithm** is the workhorse for sampling from posterior distributions that cannot be evaluated analytically. At iteration $s$:

1. Propose $\theta^* \sim q(\theta^*|\theta^{(s-1)})$ from a proposal distribution.
2. Compute acceptance probability $\alpha = \min\!\left(1, \frac{p(\theta^*|y)\cdot q(\theta^{(s-1)}|\theta^*)}{p(\theta^{(s-1)}|y)\cdot q(\theta^*|\theta^{(s-1)})}\right)$.
3. Set $\theta^{(s)} = \theta^*$ with probability $\alpha$, else $\theta^{(s)} = \theta^{(s-1)}$.

After a burn-in period, the chain converges to samples from the posterior $p(\theta|y)$.

**Model comparison** in a Bayesian framework uses the **marginal likelihood** (or evidence):

$$p(y|M) = \int p(y|\theta, M)\,p(\theta|M)\,\mathrm{d}\theta.$$

The **Bayes factor** for model $M_1$ against $M_2$ is $BF_{12} = p(y|M_1)/p(y|M_2)$. A Bayes factor above 10 is considered strong evidence for $M_1$. The **Deviance Information Criterion (DIC)** and **Watanabe-Akaike Information Criterion (WAIC)** are alternatives to the Bayes factor for model comparison that are less sensitive to prior specification.

---

*See also Appendix B (Research Methods) for applications of these tools in macroeconomic research.*
