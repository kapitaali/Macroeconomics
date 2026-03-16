# Chapter 6 — Macroeconomic Data and Sources: How We Know What We Know

> *"The first step of wisdom is to know the facts."*
> — Carl von Clausewitz

---

The theories and models of the previous chapters have empirical content only to the extent that their predictions can be compared with data. But macroeconomic data do not fall from the sky fully formed — they are constructed by statistical agencies from surveys, administrative records, and financial reports, according to methodological conventions that change over time, are never perfectly transparent, and embed numerous approximations and judgments. A macroeconomist who does not understand how the data are constructed cannot properly evaluate the evidence for or against any empirical claim. This chapter describes the architecture of the major data systems, the statistical properties of macroeconomic time series, and the econometric tools most commonly used to extract structure from them.

---

## 6.1 The Architecture of National Statistical Systems

The national accounts described in Chapter 4 are constructed according to internationally harmonized standards, the most important of which is the **System of National Accounts 2008 (SNA 2008)**, a joint publication of the United Nations, the IMF, the World Bank, the OECD, and Eurostat. The SNA 2008 provides the conceptual definitions, accounting rules, and classification systems that underlie the national accounts of virtually every country in the world. Without this harmonization, international comparisons of GDP would be impossible.

The SNA 2008 is supplemented by several specialized manuals. The **Balance of Payments and International Investment Position Manual, Sixth Edition (BPM6)** governs the measurement of cross-border transactions and external asset positions. The **Government Finance Statistics Manual 2014 (GFSM 2014)** covers the public sector accounts. The **Monetary and Financial Statistics Manual (MFSM)** covers the financial sector. Together, these manuals define the concepts, boundaries, and classification systems that allow the accounts of different sectors and different countries to be combined consistently.

### Data Revisions and Vintages

A feature of national accounts data that is poorly appreciated outside of the specialist community is that all macroeconomic data are subject to revision. The preliminary estimate of U.S. GDP — released approximately 30 days after the end of the reference quarter — is based on incomplete source data. It is revised again at 60 days (second estimate), 90 days (third estimate), and then in a series of annual revisions over the following three years as more complete survey data become available. Benchmark revisions, which occur every five years or so, can revise GDP estimates back for decades.

**Definition (Data Vintage).** A **data vintage** is the version of a time series that was available at a specific historical date, before subsequent revisions. Real-time data refers to the vintage available to policymakers and agents at the time decisions were made; revised data refers to the latest available estimates incorporating all subsequent information.

The distinction between real-time and revised data matters enormously for evaluating policy decisions and for estimating policy rules. Croushore (2011) documents that U.S. preliminary GDP estimates are revised by a mean absolute amount of approximately 1.2 percentage points at annual rates. This means that in real time, policymakers often do not know whether the economy is in recession or whether inflation is above or below target. Orphanides (2001) shows that the U.S. Federal Reserve's apparently poor performance in the 1970s largely disappears when the Taylor rule is estimated using real-time data rather than revised data — the Fed was not ignoring inflation; it was acting on badly mismeasured data. The Federal Reserve Bank of Philadelphia's **Real-Time Data Set for Macroeconomists (RTDSM)** archives all U.S. macro data vintages since 1965, enabling this kind of analysis.

---

## 6.2 Index Number Theory

Since macroeconomic data aggregate heterogeneous quantities that cannot be added directly — how do you sum apples and automobiles? — aggregation requires the use of price or quantity index numbers. The theory of index numbers asks: what properties should a good price or quantity index satisfy, and which practical index formulas best satisfy them?

A **price index** $P(\mathbf{p}_t, \mathbf{q}_t;\, \mathbf{p}_0, \mathbf{q}_0)$ summarizes the change in the price level from period 0 to period $t$, holding some reference quantity vector constant. Fisher (1922) proposed a set of **axiomatic tests** that a price index should satisfy:

1. **Identity**: $P(\mathbf{p}, \mathbf{q};\, \mathbf{p}, \mathbf{q}) = 1$ — an index comparing a period to itself equals one.
2. **Proportionality**: $P(\lambda\mathbf{p}, \mathbf{q};\, \mathbf{p}_0, \mathbf{q}_0) = \lambda$ — if all prices scale by $\lambda$, the index scales by $\lambda$.
3. **Time reversal**: $P(\mathbf{p}_0;\, \mathbf{p}_t) = 1 / P(\mathbf{p}_t;\, \mathbf{p}_0)$ — reversing the comparison inverts the index.
4. **Transitivity**: $P_{0,2} = P_{0,1} \times P_{1,2}$ — the index for a two-period span equals the product of the sub-period indices.

No commonly used index formula satisfies all these properties simultaneously. The Laspeyres index (base-period weights) satisfies (2) but violates (3) and (4). The Paasche index (current-period weights) violates (3). The Fisher ideal index — the geometric mean of Laspeyres and Paasche — satisfies (3) and is the best second-order approximation to the true cost-of-living index but violates (4). The Törnqvist index satisfies (4) but is only approximately equal to the true index. National statistical agencies make different choices: the U.S. BEA uses chain-weighted Fisher indices for real GDP; the U.S. BLS uses a Laspeyres-type index for the CPI, with periodic rebasing to reduce but not eliminate substitution bias.

---

## 6.3 Time-Series Properties of Macroeconomic Data

Before estimating any relationship between macroeconomic variables, we need to understand the statistical properties of those variables. The most important distinction is between **stationary** and **non-stationary** time series.

**Definition (Stationarity).** A time series $\{x_t\}$ is **covariance-stationary** (or weakly stationary) if its mean, variance, and autocovariances are all finite and time-invariant: $\mathbb{E}[x_t] = \mu$ for all $t$, $\text{Var}(x_t) = \sigma^2 < \infty$ for all $t$, and $\text{Cov}(x_t, x_{t-k}) = \gamma_k$ depends only on the lag $k$, not on $t$.

Most macroeconomic time series are not stationary. Log real GDP $y_t = \ln Y_t$ drifts upward over time; its mean changes with every observation. Two competing descriptions of this non-stationarity have different implications. The **trend-stationary** model:

$$y_t = \alpha + \beta t + u_t, \quad u_t \text{ stationary},$$

says that $y_t$ fluctuates around a deterministic trend $\alpha + \beta t$; shocks have temporary effects and the economy returns to the trend after disturbances. The **unit-root (difference-stationary)** model:

$$\Delta y_t = \mu + \epsilon_t, \quad \epsilon_t \sim \text{WN}(0, \sigma^2),$$

says that the *changes* in $y_t$ are stationary but the *level* is not; shocks have permanent effects and the economy never returns to its pre-shock trend. Nelson and Plosser (1982) argued that U.S. macroeconomic series are best described as unit-root processes — a finding with important implications for how we interpret business cycles.

**Definition (Unit Root).** A time series $\{x_t\}$ has a **unit root** if it follows $x_t = \rho x_{t-1} + \epsilon_t$ with $\rho = 1$. In this case, a shock $\epsilon_t$ permanently shifts the level of $x_t$ for all future periods — shocks accumulate without dying out. The Dickey–Fuller test is the standard procedure for testing $\rho = 1$ against $\rho < 1$.

### The Hodrick–Prescott Filter

The most widely used method to decompose a macroeconomic series into trend and cycle components is the **Hodrick–Prescott (HP) filter** (Hodrick and Prescott, 1997). Given the log-output series $\{y_t\}_{t=1}^T$, the HP filter extracts a trend component $\{\tau_t\}$ by solving:

$$\min_{\{\tau_t\}} \sum_{t=1}^T (y_t - \tau_t)^2 + \lambda \sum_{t=2}^{T-1} \bigl[(\tau_{t+1} - \tau_t) - (\tau_t - \tau_{t-1})\bigr]^2.$$

The first term penalizes deviations of the trend from the data; the second penalizes changes in the growth rate of the trend (its second difference). The smoothing parameter $\lambda$ governs the trade-off: large $\lambda$ forces the trend to be nearly linear; small $\lambda$ allows the trend to track the data closely. The conventional value for quarterly data is $\lambda = 1{,}600$.

The solution in matrix form is $\hat{\boldsymbol{\tau}} = (I + \lambda K'K)^{-1}\mathbf{y}$, where $K$ is the second-difference matrix. The **cyclical component** is then $\hat{c}_t = y_t - \hat{\tau}_t$. Hamilton (2018) argues that the HP filter has serious statistical problems — it introduces spurious cyclical patterns, particularly at the endpoints of the sample — and proposes instead a simple regression-based filter. The debate over detrending methods remains active and reflects genuine uncertainty about whether business cycles are deviations from a deterministic trend or permanent shifts in an evolving stochastic trend.

---

## 6.4 Cross-Country Data and Panel Methods

The study of long-run growth and development requires cross-country data spanning many decades. The **Penn World Tables (PWT)**, maintained by Feenstra, Inklaar, and Timmer (2015), provide GDP, capital stocks, and factor inputs for 180+ countries from 1950, expressed in a common unit using **purchasing power parity (PPP)** exchange rates. PPP conversion removes differences in price levels across countries: a dollar of GDP in India buys more in terms of real goods and services than a dollar of GDP in the United States, and PPP adjustment corrects for this.

**Definition (Purchasing Power Parity).** **Purchasing power parity (PPP)** is an exchange rate conversion factor that equalizes the price of a reference basket of goods across countries. PPP-adjusted GDP is GDP expressed in terms of what it can buy in a reference country (usually the United States), allowing meaningful cross-country comparisons of real living standards. Market exchange rates reflect traded goods prices and financial flows; PPP rates better reflect the true cost of living, particularly for non-traded services.

The **Maddison Project Database** (Bolt et al., 2018) extends historical coverage back to 1820 for many countries and to 1000 CE for a smaller set, enabling the study of growth over the very long run. Using these data, one can see that the current vast divergence in living standards between rich and poor countries is largely a product of the past two centuries: before the Industrial Revolution, per-capita income differences across countries were modest by modern standards.

The standard empirical tool for cross-country growth analysis is the **panel regression**:

$$y_{it} = \alpha_i + \lambda_t + \mathbf{x}_{it}'\boldsymbol{\beta} + \epsilon_{it},$$

where $\alpha_i$ is a country fixed effect (absorbing all time-invariant country characteristics — geography, legal tradition, culture), $\lambda_t$ is a time fixed effect (absorbing global trends common to all countries), $\mathbf{x}_{it}$ is a vector of time-varying country characteristics (saving rate, human capital, institutions), and $\boldsymbol{\beta}$ is the vector of parameters of interest. The fixed effects control for omitted variables that are constant within countries, substantially reducing omitted variable bias.

---

*Next: Chapter 7 — The Aggregate Demand–Aggregate Supply Model*
