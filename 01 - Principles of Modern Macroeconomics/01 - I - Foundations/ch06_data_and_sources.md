# Chapter 6 — Macroeconomic Data and Sources: How We Know What We Know

> *"The first step of wisdom is to know the facts."*
> — Carl von Clausewitz

---

The theories and models of the previous chapters have empirical content only to the extent that their predictions can be compared with data. But macroeconomic data do not fall from the sky fully formed — they are constructed by statistical agencies from surveys, administrative records, and financial reports, according to methodological conventions that change over time, embed numerous approximations, and are never perfectly transparent. A macroeconomist who does not understand how the data are constructed cannot properly evaluate the evidence for or against any empirical claim. A policy analyst who ignores data revisions may draw incorrect conclusions about the state of the economy. And a researcher who applies time-series methods to non-stationary data without testing for unit roots risks estimating spurious relationships that will not replicate. This chapter describes the architecture of the major data systems, the statistical properties of macroeconomic time series, and the econometric tools most commonly used to extract structure from them.

---

## 6.1 The Architecture of National Statistical Systems

### International Harmonization

The national accounts described in Chapter 4 are constructed according to internationally harmonized standards, the most important of which is the **System of National Accounts 2008 (SNA 2008)**, a joint publication of the United Nations, the IMF, the World Bank, the OECD, and Eurostat. The SNA 2008 provides the conceptual definitions, accounting rules, and classification systems that underlie the national accounts of virtually every country in the world. Without this harmonization, international comparisons of GDP — which form the empirical basis of growth research, development policy, and international macroeconomics — would be impossible.

The SNA 2008 is supplemented by several specialized manuals. The **Balance of Payments and International Investment Position Manual, Sixth Edition (BPM6)** governs the measurement of cross-border transactions and external asset positions [Ch. 26]. The **Government Finance Statistics Manual 2014 (GFSM 2014)** covers the public sector accounts, providing definitions of fiscal balance concepts consistent with the national accounts. The **Monetary and Financial Statistics Manual (MFSM)** covers financial sector statistics. Together, these manuals define the concepts, boundaries, and classification systems that allow the accounts of different sectors and different countries to be combined consistently.

### Data Revisions and Vintages

A feature of national accounts data that is poorly appreciated outside the specialist community is that all macroeconomic data are subject to revision — sometimes large revision. The preliminary estimate of U.S. GDP — released approximately 30 days after the end of the reference quarter — is based on incomplete source data: the Advance Estimate uses only about 65% of the source data that will eventually be available. It is revised again at 60 days (Second Estimate), 90 days (Third Estimate), and then in a series of annual revisions over the following three years as more complete survey data become available. Benchmark revisions, occurring roughly every five years, can revise GDP estimates back for decades and sometimes alter the measured depth of recessions and recoveries substantially.

**Definition (Data Vintage).** A **data vintage** is the version of a time series that was available at a specific historical date, before subsequent revisions. Real-time data refers to the vintage available to policymakers and agents at the time decisions were made; revised data refers to the latest available estimates incorporating all subsequent information.

The distinction between real-time and revised data matters enormously for evaluating policy decisions and estimating policy rules. Croushore (2011) documents that U.S. preliminary GDP estimates are revised by a mean absolute amount of approximately 1.2 percentage points at annual rates — a substantial fraction of the typical business cycle movement. This means that in real time, policymakers often genuinely do not know whether the economy has entered recession, whether the output gap is positive or negative, or whether inflation is above or below target. Orphanides (2001) demonstrates that the U.S. Federal Reserve's apparently poor performance in the 1970s largely disappears when the Taylor rule is estimated using real-time data rather than revised data: the Fed was not ignoring inflation; it was acting on badly mismeasured data with large real-time output gap errors. The Federal Reserve Bank of Philadelphia's **Real-Time Data Set for Macroeconomists (RTDSM)** archives all U.S. macro data vintages since 1965, enabling researchers to construct the information sets that policymakers actually had available when making decisions.

### Major U.S. Data Sources

The principal U.S. statistical agencies and their primary releases:

**Bureau of Economic Analysis (BEA)**: the National Income and Product Accounts (GDP, consumption, investment, government spending, net exports); the International Transactions Accounts (balance of payments); the Industry Economic Accounts (input-output tables); and the State and Metropolitan Area Accounts. The BEA publishes the quarterly GDP release and the annual revision each July.

**Bureau of Labor Statistics (BLS)**: the Consumer Price Index (CPI); the Producer Price Index (PPI); the Current Employment Statistics (the monthly payroll survey); the Current Population Survey (the household survey yielding the unemployment rate); and the Job Openings and Labor Turnover Survey (JOLTS, which provides data on vacancies essential for the Beveridge curve analysis of Chapter 31).

**Federal Reserve Board and the Federal Reserve Banks**: the Flow of Funds (now Financial Accounts of the United States), which tracks financial stocks and flows across sectors; the H.6 money stock measures; and the Federal Reserve Bank of St. Louis's FRED database, which aggregates over 800,000 time series from 100+ sources into a single freely accessible interface.

---

## 6.2 Index Number Theory

### The Aggregation Problem

Since macroeconomic data aggregate heterogeneous quantities that cannot be added directly — how do you add automobiles to apples or software to steel? — aggregation requires the use of price or quantity index numbers. The theory of index numbers asks: what properties should a good price or quantity index satisfy, and which practical index formulas best satisfy them?

A **price index** $P(\mathbf{p}_t, \mathbf{q}_t;\, \mathbf{p}_0, \mathbf{q}_0)$ summarizes the change in the price level from a reference period 0 to period $t$. Fisher (1922) proposed a set of **axiomatic tests**:

1. **Identity**: $P(\mathbf{p}, \mathbf{q};\, \mathbf{p}, \mathbf{q}) = 1$.
2. **Proportionality**: $P(\lambda\mathbf{p}, \mathbf{q};\, \mathbf{p}_0, \mathbf{q}_0) = \lambda$ — if all prices scale by $\lambda$, the index scales by $\lambda$.
3. **Time reversal**: $P(\mathbf{p}_0;\, \mathbf{p}_t) = 1 / P(\mathbf{p}_t;\, \mathbf{p}_0)$ — reversing the comparison inverts the index.
4. **Transitivity**: $P_{0,2} = P_{0,1} \times P_{1,2}$ — the index for a two-period span equals the product of sub-period indices.

No commonly used formula satisfies all four simultaneously, creating genuine tradeoffs in index design.

### The Main Index Formulas and Their Properties

The **Laspeyres index** uses base-period quantities as weights: $P^L = \sum p_t q_0 / \sum p_0 q_0$. It satisfies (2) but violates (3) and (4), and is known to exhibit **substitution bias** — it overestimates the true cost of living because it does not allow for consumers substituting toward goods whose relative prices have fallen. The **Paasche index** uses current-period quantities: $P^P = \sum p_t q_t / \sum p_0 q_t$. It exhibits the opposite bias (tending to understate true price growth) and also violates (3).

The **Fisher ideal index** — the geometric mean of Laspeyres and Paasche, $P^F = \sqrt{P^L \cdot P^P}$ — satisfies (3) and is the best second-order approximation to the true cost-of-living index, but violates (4) (it is not transitive across more than two periods). The **Törnqvist index** — a log-change weighted by average expenditure shares — satisfies (4) and provides an excellent approximation to the true index under flexible functional form assumptions, but it is not an exact price index.

National statistical agencies make different choices. The U.S. BEA uses **chain-weighted Fisher indices** for real GDP: in each period, the Fisher index compares consecutive periods, and the chain is built by multiplying consecutive period-over-period indices. This eliminates the base-year dependence of the old fixed-weight Laspeyres real GDP but at the cost that chain-weighted components no longer add up to chain-weighted GDP (the "non-additivity" problem). The U.S. BLS uses a **Laspeyres-type index** for the CPI, with periodic rebasing and geometric mean substitution within categories (to partially address substitution bias), but the CPI is still generally understood to overstate true inflation by approximately 0.3–0.5 percentage points annually due to residual substitution bias, outlet bias, and quality change.

### Quality Change and Hedonic Adjustment

A particularly difficult measurement problem for price indices is **quality change**: when a product improves in quality, not all of a price increase represents inflation — part of it represents more value per unit. The BLS uses **hedonic price indices** for goods with rapidly changing quality — computers, mobile phones, software — in which the price of a product is regressed on its measurable characteristics (processor speed, storage, screen resolution) to extract the pure price change holding quality constant. Hedonic adjustment has been credited with explaining much of the measured productivity slowdown in the 1970s and early 1980s: nominal prices rose less in quality-adjusted terms than in raw-price terms, so real output growth was higher than conventional deflators suggested.

---

## 6.3 Time-Series Properties of Macroeconomic Data

### Stationarity and Non-Stationarity

Before estimating any relationship between macroeconomic variables, researchers must understand the statistical properties of those variables. The most important distinction is between **stationary** and **non-stationary** time series.

**Definition (Stationarity).** A time series $\{x_t\}$ is **covariance-stationary** if its mean, variance, and autocovariances are all finite and time-invariant: $\mathbb{E}[x_t] = \mu$, $\mathrm{Var}(x_t) = \sigma^2 < \infty$, and $\mathrm{Cov}(x_t, x_{t-k}) = \gamma_k$ depending only on the lag $k$ and not on $t$.

Most macroeconomic levels — log real GDP, the price level, nominal money — are not stationary; their means and variances grow over time. Two competing descriptions of this non-stationarity have fundamentally different implications. The **trend-stationary model** $y_t = \alpha + \beta t + u_t$ says that output fluctuates around a deterministic trend; shocks have temporary effects and the economy returns to its trend after disturbances. The **unit-root (difference-stationary) model** $\Delta y_t = \mu + \epsilon_t$ says that changes in output are stationary but the level is not; shocks have permanent effects and the economy never returns to its pre-shock trend. Nelson and Plosser (1982) argued that U.S. macroeconomic series are best described as unit-root processes, a finding with profound implications for business cycle theory: if shocks are permanent, business cycles are not deviations from a stable trend but represent the accumulation of permanent shifts.

**Definition (Unit Root).** A time series $\{x_t\}$ has a **unit root** if $x_t = x_{t-1} + \epsilon_t$ (or more generally, if the autoregressive polynomial has a root equal to one). Shocks permanently shift the level of $x_t$ for all future periods — they accumulate without mean-reversion. The standard test is the **Augmented Dickey-Fuller (ADF) test**: regress $\Delta x_t = \alpha + \beta t + \rho x_{t-1} + \sum_j \gamma_j \Delta x_{t-j} + \epsilon_t$ and test $H_0: \rho = 0$ (unit root) against $H_1: \rho < 0$ (stationarity). Critical values are non-standard (due to the non-standard distribution of the OLS estimator when $\rho = 0$) and must be taken from the Dickey-Fuller tables rather than the standard normal.

### The Hodrick-Prescott Filter

The most widely used method to decompose a macroeconomic series into trend and cycle components is the **Hodrick-Prescott (HP) filter**. Given log-output $\{y_t\}_{t=1}^T$, the HP filter solves:

$$\min_{\{\tau_t\}} \sum_{t=1}^T (y_t - \tau_t)^2 + \lambda \sum_{t=2}^{T-1} \bigl[(\tau_{t+1} - \tau_t) - (\tau_t - \tau_{t-1})\bigr]^2.$$

The first term penalizes deviations of the trend from the data; the second penalizes acceleration in the trend's growth rate (its second difference). The smoothing parameter $\lambda$ governs the trade-off: large $\lambda$ forces the trend to be nearly linear; small $\lambda$ allows the trend to track the data closely. The conventional value for quarterly data is $\lambda = 1{,}600$, chosen to yield trend cycles with the periodicity of typical business cycles.

The matrix solution is $\hat{\boldsymbol{\tau}} = (I + \lambda K'K)^{-1}\mathbf{y}$, where $K$ is the second-difference matrix. The cyclical component $\hat{c}_t = y_t - \hat{\tau}_t$ is then used to compute the business cycle statistics in Chapter 27.

Hamilton (2018) argues that the HP filter has serious statistical pathologies: it introduces spurious cyclical patterns by design, particularly near the endpoints of the sample, and the extracted cycles have spectral properties inconsistent with the raw data. He proposes instead regressing $y_{t+h}$ on $\{y_t, y_{t-1}, y_{t-2}, y_{t-3}\}$ and using the residuals as the cycle measure. The debate over detrending remains active and reflects genuine uncertainty about whether business cycles are deviations from a deterministic trend or permanent shifts in an evolving stochastic trend — a question with substantive implications for whether stabilization policy can meaningfully smooth them.

### Cointegration

When two or more non-stationary series share a common stochastic trend, they are said to be **cointegrated**. Formally, $x_t \sim I(1)$ and $z_t \sim I(1)$ are cointegrated if there exists a vector $\beta$ such that $\beta' (x_t, z_t)'$ is $I(0)$ — stationary despite the individual series being non-stationary. Cointegration represents a long-run equilibrium relationship: though the series individually wander as random walks, they wander together and never drift arbitrarily far apart.

The **Engle-Granger theorem** (Engle and Granger, 1987) states that cointegrated variables have an **error correction representation**: $\Delta x_t = \alpha(\beta' x_{t-1}) + \text{lagged differences} + \epsilon_t$, where $\alpha < 0$ is the error correction coefficient governing how fast the system returns to its long-run relationship after a shock. Economically important cointegrating relationships include: money demand (money, prices, and income); purchasing power parity (domestic and foreign price levels and the nominal exchange rate); and the term structure of interest rates (long and short rates).

---

## 6.4 Cross-Country Data and Panel Methods

### International Data Sources

The study of long-run growth and development requires cross-country data spanning many decades. The **Penn World Tables (PWT)**, maintained by Feenstra, Inklaar, and Timmer (2015) and regularly updated, provide GDP, capital stocks, and factor inputs for 180+ countries from 1950, expressed in a common unit using **purchasing power parity (PPP)** exchange rates. PPP conversion removes differences in price levels: a dollar of GDP in India buys more in terms of real goods and services than a dollar of GDP in the United States, and market exchange rates do not correct for this.

**Definition (Purchasing Power Parity Exchange Rate).** A **PPP exchange rate** is the number of units of domestic currency required to buy the same quantity of goods and services that one unit of a reference currency buys in the reference country. PPP-adjusted GDP measures what a country's output could purchase if priced at a common international set of prices, enabling genuine cross-country comparisons of real living standards. For non-traded services (healthcare, haircuts, education), PPP rates and market rates can diverge by a factor of three or more in comparisons between rich and poor countries.

The **Maddison Project Database** (Bolt et al., 2018) extends historical coverage back to 1820 for many countries and further for a smaller set, enabling the study of growth across two centuries. A striking feature: the vast divergence in living standards between rich and poor countries that characterizes the modern world is largely a product of the past two centuries. Before the Industrial Revolution, per-capita income differences across countries were modest by modern standards — the ratio of richest to poorest countries was perhaps 4:1 in 1800 compared to roughly 50:1 today.

### Panel Regression Methods

The standard empirical tool for cross-country growth analysis is the **panel regression**:

$$y_{it} = \alpha_i + \lambda_t + \mathbf{x}_{it}'\boldsymbol{\beta} + \epsilon_{it},$$

where $\alpha_i$ are country fixed effects (absorbing all time-invariant country characteristics — geography, legal tradition, culture, initial conditions), $\lambda_t$ are time fixed effects (absorbing global trends common to all countries in each period), $\mathbf{x}_{it}$ is a vector of time-varying country characteristics, and $\boldsymbol{\beta}$ is the vector of structural parameters.

The fixed effects estimator controls for omitted country-level variables that are constant within countries, substantially reducing omitted variable bias relative to cross-sectional OLS. The **Hausman test** determines whether the random effects or fixed effects estimator is appropriate: under the null of no correlation between $\alpha_i$ and $\mathbf{x}_{it}$, both are consistent but random effects is more efficient; under the alternative (correlation present), only fixed effects is consistent.

### Instrumental Variables in Growth Regressions

A persistent challenge in empirical growth research is endogeneity: institutions, investment, trade, and schooling are all simultaneously determined with income, making OLS estimates of their effects on growth difficult to interpret causally. **Instrumental variables (IV) estimation** addresses this by finding variables $\mathbf{z}_{it}$ (instruments) that affect the outcome only through the endogenous regressor (the **exclusion restriction**) and are sufficiently correlated with the endogenous regressor (the **relevance condition**).

Acemoglu, Johnson, and Robinson (2001) use early settler mortality rates as an instrument for institutional quality: where Europeans died from tropical disease, they set up extractive institutions; where they survived, they set up inclusive settler colonies. Mortality rates are plausibly exogenous to current income (they reflect historical geography, not current policy) and strongly correlated with current institutional quality — a compelling instrument. Their IV estimates of the causal effect of institutions on income per capita are substantially larger than OLS estimates, suggesting that the OLS coefficient on institutions is downward biased because poor institutions tend to accompany other poverty-causing factors.

---

## Chapter Summary

- National accounts are constructed according to internationally harmonized **SNA 2008** standards. Preliminary estimates are revised substantially — U.S. GDP preliminary estimates carry mean absolute revisions of approximately 1.2 percentage points — making real-time versus revised data distinctions critical for policy analysis and historical evaluation.

- **Index number theory** identifies the tradeoffs among Laspeyres (substitution-biased upward), Paasche (substitution-biased downward), Fisher ideal (best cost-of-living approximation, non-transitive), and Törnqvist (transitive, excellent approximation) formulas. The U.S. uses chain-weighted Fisher indices for GDP and a Laspeyres-type index for CPI; residual substitution and quality-change biases mean CPI likely overstates true inflation by approximately 0.3–0.5 percentage points annually.

- Most macroeconomic levels are **non-stationary** (unit root processes), as documented by Nelson and Plosser (1982). The ADF test evaluates $H_0: \rho = 0$; critical values follow a non-standard distribution. The HP filter ($\lambda = 1{,}600$ for quarterly data) extracts the cyclical component but has statistical pathologies near sample endpoints (Hamilton, 2018). **Cointegration** captures long-run equilibrium relationships among I(1) variables; cointegrated systems have an error-correction representation.

- **PPP-adjusted GDP** (Penn World Tables) enables cross-country comparisons by removing price level differences; **Maddison Project** data extend historical coverage to 1820+. Panel regressions with country and time fixed effects control for time-invariant omitted variables; IV estimation (Acemoglu et al.'s settler mortality instrument) addresses endogeneity in growth regressions.

---

*Next: Chapter 7 — The Aggregate Demand–Aggregate Supply Model*
