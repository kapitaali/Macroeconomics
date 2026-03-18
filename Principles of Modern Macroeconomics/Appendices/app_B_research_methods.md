# Appendix B — Research Methods in Macroeconomics

This appendix provides a concise introduction to the empirical and theoretical methods used in macroeconomic research. It is intended as a bridge between the substantive material of the main text and the primary literature, equipping the reader with enough methodological vocabulary to read working papers and understand the claims being made.

---

## B.1 The Structure of Macroeconomic Research

Macroeconomic research proceeds along two tracks that are constantly in dialogue but use different tools.

**Theoretical research** constructs and analyzes formal models — systems of equations describing the behavior of optimizing agents in equilibrium. The goal is to identify mechanisms: to show that a particular set of assumptions is sufficient to generate a particular pattern of outcomes, or that a given policy instrument has certain effects. Theoretical papers establish existence and characterization results — does an equilibrium exist? is it unique? is it stable? — and derive comparative statics — how do equilibrium outcomes change when parameters shift?

**Empirical research** confronts theoretical predictions with data. It attempts to measure the parameters that theoretical models require, test whether models' qualitative and quantitative predictions are borne out in the data, and identify the causal effects of specific shocks or policy changes. Empirical macroeconomics is harder than empirical microeconomics for a fundamental reason: the "data-generating process" for macroeconomic aggregates is a single, unrepeatable history of the world. We cannot randomly assign recessions to countries or run controlled experiments on monetary policy.

Modern macroeconomic research typically combines both tracks: a theoretical model motivates an empirical exercise that estimates its key parameters; the empirical results are then used to calibrate or reject the model; and discrepancies motivate refinements of the theory. This iterative process — model, estimate, revise — is the engine of progress in the field.

---

## B.2 Calibration

**Definition (Calibration).** Calibration is the practice of assigning values to model parameters based on empirical observations or micro-level studies, rather than estimating them jointly from aggregate time-series data. A calibrated model is then simulated and its predictions compared to observed data moments (means, variances, correlations) to evaluate its adequacy.

Calibration was introduced to macroeconomics by Kydland and Prescott (1982) in the context of the real business cycle model. Their argument: some parameters (the capital share $\alpha$, the depreciation rate $\delta$) are well-identified from national accounts and can be measured directly from data; estimating them jointly with the remaining parameters from a short aggregate time series introduces noise and potential bias. Calibrating from independent data sources uses information more efficiently.

A typical calibration proceeds as follows:

1. Identify the parameters of the model: $\{\beta, \sigma, \alpha, \delta, n, g, \rho_A, \sigma_A, \ldots\}$.
2. Assign values to well-identified parameters from microeconomic or accounting data: $\alpha \approx 1/3$ from labor's share, $\delta \approx 0.025$ from investment data, $n \approx 0.01$ from population statistics, etc.
3. Set $\beta$ to match the average real interest rate in the data: $r = (1-\beta)/\beta \approx 0.04$ implies $\beta \approx 0.99$ at quarterly frequency.
4. Set $\sigma$ from microeconomic estimates of the EIS or from data on the response of consumption to interest rate changes.
5. Set shock parameters $(\rho_A, \sigma_A)$ to match properties of the Solow residual.
6. Simulate the model and compute second moments (standard deviations, correlations, autocorrelations) of simulated series.
7. Compare simulated moments to empirical moments; evaluate the model's success in matching the data.

The calibration approach has been criticized for its informality: different researchers calibrating the same model can produce different results by making different choices in steps 2–4. The Bayesian estimation approach (Section B.5) addresses this by estimating all parameters simultaneously from the data.

---

## B.3 Vector Autoregressions (VARs)

The **vector autoregression (VAR)** is the principal reduced-form tool for characterizing the joint time-series behavior of macroeconomic variables and for identifying the dynamic effects of structural shocks.

A $k$-variable VAR of order $p$ is:

$$\mathbf{y}_t = \mathbf{c} + \sum_{j=1}^p A_j\,\mathbf{y}_{t-j} + \mathbf{u}_t, \quad \mathbf{u}_t \sim \mathcal{N}(\mathbf{0}, \Sigma),$$

where $\mathbf{y}_t$ is a $k\times 1$ vector of endogenous variables (e.g., output, inflation, the interest rate), $\mathbf{c}$ is a vector of constants, $A_j$ are $k\times k$ coefficient matrices, and $\mathbf{u}_t$ is a vector of reduced-form residuals with contemporaneous covariance matrix $\Sigma$. Each equation of the VAR is an OLS regression of one variable on lags of all variables; OLS is efficient under the null of no serial correlation in $\mathbf{u}_t$.

The reduced-form residuals $\mathbf{u}_t$ are generally correlated across equations (because structural shocks affect multiple variables simultaneously). To recover the underlying structural shocks, we write $\mathbf{u}_t = B\bm{\varepsilon}_t$, where $\bm{\varepsilon}_t \sim \mathcal{N}(\mathbf{0}, I_k)$ are orthogonal structural shocks and $B$ is a $k\times k$ matrix capturing their contemporaneous effects. The data pin down $\Sigma = BB'$, providing $k(k+1)/2$ equations in the $k^2$ elements of $B$. To identify $B$ uniquely, $k(k-1)/2$ additional restrictions must be imposed.

### Identification Strategies

**Cholesky (recursive) identification.** Impose a lower triangular structure on $B$: the ordering of variables in $\mathbf{y}_t$ determines which variables respond to each shock contemporaneously. Variable 1 is affected only by shock 1; variable 2 by shocks 1 and 2; and so on. This is appropriate when there is a clear theoretical justification for the recursive structure (e.g., prices do not respond to output within a quarter, but output responds to price surprises).

**Sign restrictions.** Impose restrictions on the signs of impulse response functions rather than on the elements of $B$ directly. For example, identify a monetary policy shock as one that raises the interest rate, reduces output, and reduces inflation on impact, without specifying the magnitude of these responses. Sign restrictions typically identify a set of admissible decompositions rather than a unique one, and the results are summarized by the set of impulse responses consistent with the restrictions.

**External instruments (proxy SVARs).** Use a variable that is correlated with the structural shock of interest but uncorrelated with all other structural shocks as an instrument to identify the shock (Stock and Watson, 2012; Mertens and Ravn, 2013). For monetary policy shocks, Romer and Romer (2004) use policy intentions extracted from Federal Reserve Greenbook forecasts; for fiscal shocks, Ramey (2011) uses defense spending news.

**Long-run restrictions.** Impose the theoretical restriction that certain shocks have no long-run effect on certain variables (Blanchard and Quah, 1989). For example, identify aggregate demand shocks as shocks with no long-run effect on output (only supply shocks affect long-run output), and aggregate supply shocks as shocks with a permanent output effect.

### Impulse Response Functions

The **impulse response function (IRF)** traces the dynamic response of each variable in $\mathbf{y}_t$ to a unit innovation in one structural shock, holding all other shocks at zero:

$$\mathrm{IRF}_{ij}(h) = \frac{\partial y_{i,t+h}}{\partial \varepsilon_{j,t}}, \quad h = 0, 1, 2, \ldots$$

IRFs are the primary tool for evaluating model predictions against VAR evidence: a well-specified DSGE model should produce theoretical IRFs that fall within the confidence bands of the empirical IRFs from the VAR.

---

## B.4 Identification of Causal Effects

The fundamental challenge in empirical macroeconomics is distinguishing correlation from causation. Does fiscal expansion cause output to rise, or do governments spend more when output rises? Does monetary tightening cause inflation to fall, or does the central bank raise rates when it observes rising inflation? Observational data alone cannot answer these questions without additional identifying assumptions.

### Instrumental Variables

An **instrumental variable (IV)** $Z_t$ identifies the causal effect of a treatment variable $X_t$ on outcome $Y_t$ by exploiting variation in $X_t$ that is exogenous — uncorrelated with the error term in the structural equation. The IV estimator consistently estimates the structural parameter if: (i) **relevance**: $Z_t$ is correlated with $X_t$; and (ii) **exclusion**: $Z_t$ affects $Y_t$ only through $X_t$, not through any other channel.

The IV estimate of the fiscal multiplier in Blanchard and Perotti (2002) exploits the fact that fiscal spending within a quarter is largely predetermined by the budget process, providing quasi-exogenous variation. The IV estimate of the effects of institutions on growth in Acemoglu, Johnson, and Robinson (2001) exploits settler mortality as an instrument for institutional quality.

### Regression Discontinuity

The **regression discontinuity (RD)** design identifies causal effects by exploiting discontinuous jumps in the probability of treatment at a known threshold in a running variable. Close to the threshold, treatment assignment is approximately as-good-as-random, providing local identification of the causal effect. RD designs have been used in macroeconomics to study, for instance, the effects of crossing EU fiscal rules (the 3% deficit threshold) on bond spreads, or the effects of crossing income thresholds on unemployment insurance eligibility.

### Difference-in-Differences

The **difference-in-differences (DiD)** estimator compares the change in outcomes over time between a treatment group (exposed to a policy change) and a control group (not exposed), using the control group's time series to construct the counterfactual. Validity requires the **parallel trends assumption**: in the absence of treatment, the treatment and control groups would have evolved identically. DiD is widely used in subnational fiscal multiplier studies (comparing counties that receive different amounts of federal stimulus spending) and in studying the effects of state-level minimum wage changes.

---

## B.5 Structural Estimation of DSGE Models

The calibration approach of Section B.2 sets parameters from independent information but does not formally test whether the model is consistent with the aggregate time-series data. **Structural estimation** estimates model parameters jointly from aggregate data, allowing formal model comparison and hypothesis testing.

### Maximum Likelihood Estimation

The **likelihood function** of a linearized DSGE model can be evaluated using the **Kalman filter**, which recursively computes the probability of the observed data given the model's state-space representation. Let $\mathbf{y}_t$ be the vector of observed variables and $\mathbf{s}_t$ the vector of state variables (capital, productivity, etc.). The linearized DSGE implies:

$$\mathbf{s}_{t+1} = A(\theta)\,\mathbf{s}_t + B(\theta)\,\bm{\varepsilon}_{t+1}$$
$$\mathbf{y}_t = C(\theta)\,\mathbf{s}_t + D(\theta)\,\bm{\eta}_t,$$

where $\theta$ is the vector of structural parameters and $\bm{\eta}_t$ is measurement error. The Kalman filter evaluates $\mathcal{L}(\theta) = p(\mathbf{y}_1, \ldots, \mathbf{y}_T | \theta)$ recursively, and the MLE is $\hat{\theta} = \arg\max_\theta \mathcal{L}(\theta)$.

### Bayesian Estimation

Most modern DSGE model estimation uses a **Bayesian** approach, combining a likelihood function $\mathcal{L}(\theta | \mathbf{y})$ with a prior distribution $p(\theta)$ over parameters to obtain the posterior distribution:

$$p(\theta | \mathbf{y}) \propto \mathcal{L}(\theta | \mathbf{y}) \cdot p(\theta).$$

The posterior summarizes all available information about parameters — what the data say, weighted by prior beliefs. Priors are typically calibrated from microeconomic studies (a prior on the capital share centered at 1/3) or from the theoretical restrictions of the model (parameters must lie in ranges consistent with existence of a bounded equilibrium). The posterior is sampled using Markov Chain Monte Carlo (MCMC) methods, typically the Metropolis–Hastings algorithm.

Bayesian estimation of DSGE models has several advantages: it incorporates prior information in a disciplined way; it naturally handles identification challenges by placing priors over regions where the likelihood is flat; and it produces full posterior distributions over parameters and model predictions rather than point estimates. The primary reference is Smets and Wouters (2007), whose estimated New Keynesian DSGE model for the euro area and United States has become a benchmark.

---

## B.6 Panel Data Methods in Growth and Development

Cross-country growth regressions typically use panel data with country fixed effects and time fixed effects. The canonical specification:

$$g_{it} = \alpha_i + \lambda_t + \beta \ln y_{i,t-1} + \mathbf{x}_{it}'\boldsymbol{\gamma} + \epsilon_{it},$$

where $g_{it}$ is the growth rate of per-capita income in country $i$ at time $t$, $\alpha_i$ is a country fixed effect, $\lambda_t$ is a time fixed effect (controlling for global growth fluctuations), $\ln y_{i,t-1}$ is initial log income (testing conditional convergence), and $\mathbf{x}_{it}$ is a vector of growth determinants.

Country fixed effects absorb all time-invariant country characteristics — geography, legal tradition, colonial history. The identifying variation is therefore purely within-country over time, which is often too limited to precisely identify the effects of slow-moving institutional variables. This motivated the use of long-difference specifications (comparing initial and final observations only) and instrumental variable approaches to identify the effects of institutions and policies on growth.

The **System GMM** estimator (Blundell and Bond, 1998) instruments the differenced equation with lagged levels and instruments the levels equation with lagged differences, exploiting the full dynamic structure of the panel. It is designed to handle the Nickell bias in fixed effects estimators when the time dimension is short, and is widely used in dynamic panel growth regressions.

---

## B.7 Local Projection Impulse Responses

An alternative to VAR impulse responses that has gained considerable popularity is the **local projection** method of Jordà (2005). Instead of estimating a single VAR and computing multi-step impulse responses by iteration, local projections estimate a separate regression for each horizon $h$:

$$y_{i,t+h} = \alpha_h + \beta_h\,\text{shock}_t + \mathbf{z}_t'\boldsymbol{\gamma}_h + u_{i,t+h},$$

where the coefficient $\beta_h$ directly estimates the impulse response at horizon $h$. Local projections are robust to misspecification: if the true VAR lag order is misspecified, the iterated VAR impulse responses will be biased, but the local projection impulse responses remain consistent (Montiel Olea and Plagborg-Møller, 2021). They also more naturally accommodate nonlinearities and state-dependence by interacting the shock with indicators of the economic state (recession vs. expansion; ELB vs. normal times).

The trade-off: local projections are less efficient than correctly specified VARs (they use less of the dynamic structure of the data) and can be noisy at long horizons. For many applications in macroeconomics where the correct lag order is uncertain, local projections are the preferred approach.

---

*See also Appendix D (Statistical Tools and Data Analysis) for the statistical foundations of the methods described here.*
