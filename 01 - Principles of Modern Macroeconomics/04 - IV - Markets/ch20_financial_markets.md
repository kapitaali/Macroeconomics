# Chapter 20 — Financial Markets: Stocks, Bonds, and Asset Prices

> *"Markets can remain irrational longer than you can remain solvent."*
> — Attributed to John Maynard Keynes

---

Financial markets mediate the allocation of resources across time and across states of the world. They allow households to smooth consumption over their lifetimes, firms to raise capital for investment, governments to borrow against future tax revenues, and all agents to insure against idiosyncratic risks. The prices set in financial markets — interest rates, equity prices, exchange rates — are the primary conduit through which monetary policy affects real economic activity, and the dynamics of these prices are central to understanding both the transmission of policy and the origins of financial crises.

This chapter develops the theory of asset pricing from first principles, using the **stochastic discount factor (SDF)** as the unifying framework. The SDF — also called the pricing kernel — prices all assets simultaneously in a way that is consistent with the intertemporal optimality conditions of households. Understanding the SDF framework requires integrating the consumption theory of Chapter 11 with the theory of financial markets, and doing so reveals both the power and the limits of the representative-agent approach to financial economics.

---

## 20.1 The No-Arbitrage Principle and the Stochastic Discount Factor

The foundation of asset pricing theory is the principle of no arbitrage: in equilibrium, there are no riskless profit opportunities. If two assets with identical payoffs had different prices, arbitrageurs would buy the cheaper and sell the more expensive, eliminating the price difference. In a frictionless market, no-arbitrage is almost tautological; in the presence of transaction costs and limits to arbitrage, it holds approximately.

No-arbitrage imposes a fundamental restriction on asset prices: there must exist a positive random variable $M_{t+1}$ such that the price of any asset is the expectation of $M_{t+1}$ times the asset's payoff. This is the **stochastic discount factor** (SDF), also called the **pricing kernel** or the **state-price density**.

**Definition (Stochastic Discount Factor).** The **stochastic discount factor** $M_{t+1}$ is a positive random variable such that the time-$t$ price $p_t$ of any asset with time-$(t+1)$ payoff $(p_{t+1} + d_{t+1})$ satisfies:

$$p_t = \mathbb{E}_t\bigl[M_{t+1}(p_{t+1} + d_{t+1})\bigr].$$

In terms of the gross return $R_{t+1} = (p_{t+1}+d_{t+1})/p_t$, this becomes:

$$\mathbb{E}_t[M_{t+1}\,R_{t+1}] = 1 \quad \text{for every asset.}$$

The SDF must price the risk-free asset as well: if $R_{t+1}^f = 1+r_t^f$ is the risk-free gross return (known at date $t$), then $\mathbb{E}_t[M_{t+1}]\cdot R^f = 1$, so $R^f = 1/\mathbb{E}_t[M_{t+1}]$.

The SDF contains all information about risk pricing. Subtracting the risk-free pricing equation from the general pricing equation:

$$\mathbb{E}_t[R_{t+1}^j - R_{t+1}^f] = -\frac{\mathrm{Cov}_t(M_{t+1},\, R_{t+1}^j)}{\mathbb{E}_t[M_{t+1}]}.$$

The **risk premium** on asset $j$ equals minus the covariance of the SDF with the asset's return, divided by the expected SDF. Assets that covary negatively with the SDF — those that pay poorly in states where $M_{t+1}$ is high — must offer a positive risk premium to attract investors. The SDF is high in states where the marginal utility of consumption is high — recessions, financial crises, disasters — which is precisely why assets that crash in recessions (equities, junk bonds) carry positive risk premia.

### The Representative-Agent SDF

In the consumption-CAPM (CCAPM) with a representative agent having CRRA utility $u(c) = c^{1-\sigma}/(1-\sigma)$, the optimality condition for portfolio choice gives the SDF directly:

$$M_{t+1} = \beta\frac{u'(c_{t+1})}{u'(c_t)} = \beta\left(\frac{c_{t+1}}{c_t}\right)^{-\sigma}.$$

The SDF equals the intertemporal marginal rate of substitution (IMRS): the rate at which the household is willing to trade current for future consumption. A high SDF ($M_{t+1}$ large) means consumption grew slowly or fell — a bad state — so the marginal utility of future consumption is high relative to current consumption. Assets that pay well in bad states (when $M_{t+1}$ is high) act as insurance and therefore command a negative risk premium or zero premium; assets that pay poorly in bad states carry a positive premium.

---

## 20.2 The Term Structure of Interest Rates

The **term structure of interest rates** (or yield curve) describes the relationship between the yields to maturity of government bonds and their time to maturity. Understanding the term structure is essential for monetary policy transmission, since most investment and consumption decisions depend on medium- and long-term interest rates rather than the overnight policy rate.

**Definition (Yield to Maturity).** The **yield to maturity** $i_t^n$ of an $n$-period zero-coupon bond is the constant discount rate that equates the bond's current price $p_t^n$ to its face value of \$1 at maturity:

$$p_t^n = \frac{1}{(1+i_t^n)^n}.$$

In continuous time, $p_t^n = e^{-i_t^n\cdot n}$.

### The Pure Expectations Theory

The **pure expectations theory** (also called the expectations hypothesis) holds that the long-term yield is a simple average of expected future short-term yields:

$$1 + i_t^n = \bigl[(1+\mathbb{E}_t[i_t^1])(1+\mathbb{E}_t[i_{t+1}^1])\cdots(1+\mathbb{E}_t[i_{t+n-1}^1])\bigr]^{1/n},$$

or in log approximation:

$$i_t^n = \frac{1}{n}\sum_{k=0}^{n-1}\mathbb{E}_t[i_{t+k}^1].$$

The intuition: investors are indifferent between holding an $n$-period bond or rolling over short-term bonds for $n$ periods, so the long rate must equal the expected average of short rates. Under the pure expectations theory, the yield curve is upward sloping when short rates are expected to rise and downward sloping (inverted) when they are expected to fall.

### Term Premia

The pure expectations theory is empirically rejected: even after controlling for expected future short rates, longer-maturity bonds tend to offer higher yields. The difference is the **term premium** $\phi_t^n$:

$$i_t^n = \frac{1}{n}\sum_{k=0}^{n-1}\mathbb{E}_t[i_{t+k}^1] + \phi_t^n.$$

The term premium is compensation for bearing **duration risk** — the risk that the value of a long-term bond falls when interest rates unexpectedly rise. An investor who holds a 10-year bond but may need to sell it before maturity faces interest rate risk; the term premium compensates for this.

In the SDF framework, the term premium is:

$$\phi_t^n = -\frac{1}{n}\mathrm{Cov}_t\!\left(\sum_{k=0}^{n-1}M_{t,t+k},\; i_{t+n}^1 - \mathbb{E}_t[i_{t+n}^1]\right) \times \text{(adjustment factor)},$$

which is positive when interest rate surprises are negatively correlated with the SDF — when interest rates tend to fall in recessions (when $M$ is high), long-term bonds act as insurance, commanding a lower term premium; when rates rise in recessions, long-term bonds provide poor insurance and the premium is positive.

### The Yield Curve as a Recession Predictor

The yield spread between long-term (10-year) and short-term (3-month) Treasury yields is one of the most reliable recession predictors in U.S. data (Estrella and Hardouvelis, 1991; Rudebusch and Williams, 2009). An inverted yield curve — when short-term rates exceed long-term rates — has preceded every U.S. recession since 1960, with a lead time of 6–18 months. The mechanism: yield curve inversion typically occurs when the central bank has raised short-term rates significantly above the long-term neutral rate, signaling tight monetary policy that is likely to generate a slowdown.

---

## 20.3 Equity Markets and the CAPM

### The Capital Asset Pricing Model

The **Capital Asset Pricing Model (CAPM)** (Sharpe, 1964; Lintner, 1965) is a special case of the SDF framework derived under specific assumptions: one period, mean-variance preferences, and a representative investor. It predicts that the expected excess return on any asset is proportional to its **beta** — the covariance of the asset's return with the market portfolio divided by the variance of the market:

$$\mathbb{E}[R^j] - R^f = \beta_j\bigl[\mathbb{E}[R^m] - R^f\bigr],$$

$$\beta_j = \frac{\mathrm{Cov}(R^j, R^m)}{\mathrm{Var}(R^m)}.$$

Assets with $\beta_j > 1$ are more volatile than the market and offer higher expected returns; assets with $\beta_j < 0$ (such as gold, which tends to rise in downturns) offer expected returns below the risk-free rate because they provide insurance.

The CAPM derives from the SDF framework by noting that when the SDF is linear in the market return, $M = a + bR^m$, the pricing equation yields the CAPM. This linearity holds under mean-variance preferences.

### The Equity Premium Puzzle

The most famous empirical failure of the CCAPM is the **equity premium puzzle** (Mehra and Prescott, 1985). The annual premium of U.S. equity returns over the risk-free rate has averaged approximately 6–7% since 1890. The CCAPM predicts a premium:

$$\mathbb{E}[R^e - R^f] = \sigma\cdot\mathrm{Cov}\!\left(\Delta\ln c,\, R^e\right)\approx \sigma\cdot\sigma_c\cdot\sigma_e\cdot\rho_{c,e},$$

where $\sigma_c \approx 3\%$ is consumption growth volatility, $\sigma_e \approx 16\%$ is equity return volatility, and $\rho_{c,e} \approx 0.2$ is their correlation. With $\sigma = 1$ (log utility), the predicted premium is $1 \times 0.03 \times 0.16 \times 0.2 = 0.001$ — a mere 0.1%, or 60 times smaller than the data. To match the observed premium, the model requires $\sigma \approx 50$ — a coefficient of relative risk aversion so large that it implies households would pay 50% of their wealth to avoid a 1% gamble on consumption, which is wildly at odds with microeconomic evidence.

Proposed resolutions to the equity premium puzzle include: **habit formation** (Campbell and Cochrane, 1999: when past consumption creates a reference point, effective risk aversion is much higher than the coefficient $\sigma$ in the utility function); **rare disasters** (Barro, 2006: infrequent but catastrophic consumption events — wars, depressions — generate a large precautionary premium even with moderate $\sigma$); **heterogeneous agents** (Mankiw, 1986: idiosyncratic income risk that cannot be fully insured raises the effective equity premium); and **limited stock market participation** (only a subset of households hold equities, so the relevant consumption is stockholder consumption, which is more volatile).

### Gordon Growth Model and Stock Market Valuation

For a dividend-paying stock with dividends $D_t$ growing at constant rate $g$, the **Gordon (1962) growth model** gives:

$$P_t = \frac{D_{t+1}}{i^e - g},$$

where $i^e$ is the required equity return (risk-free rate plus equity risk premium). This formula has several important implications for understanding stock market dynamics. A reduction in the risk-free rate (e.g., from QE) raises stock valuations by reducing $i^e - g$. A permanent increase in expected dividend growth $g$ (e.g., from technological innovation) also raises valuations. But a reduction in risk appetite — a rise in the equity risk premium — reduces valuations.

Shiller's **CAPE ratio** (Cyclically Adjusted Price-Earnings ratio, using 10-year average real earnings) provides an empirical implementation: it is a measure of $P/D$ adjusted for the business cycle. High CAPE predicts low subsequent 10-year real returns on equities, consistent with mean-reversion in valuations driven by time-varying risk premia rather than constant discount rates.

---

## 20.4 Efficient Markets and Empirical Anomalies

The **Efficient Market Hypothesis (EMH)** (Fama, 1970) asserts that asset prices fully reflect all available information. It has three forms of increasing strength.

**Definition (Efficient Market Hypothesis, Three Forms).** A market is:
- **Weak-form efficient** if prices reflect all information in historical prices (technical analysis cannot generate excess returns).
- **Semi-strong-form efficient** if prices reflect all publicly available information (fundamental analysis cannot generate persistent excess returns).
- **Strong-form efficient** if prices reflect all information including private information (even insiders cannot generate excess returns).

Strong-form efficiency is clearly violated (insider trading generates abnormal returns until detection). Semi-strong efficiency is inconsistent with several documented **anomalies** — patterns in returns that cannot be explained by risk compensation under standard models. The most robust anomalies include:

**Momentum** (Jegadeesh and Titman, 1993): stocks that have outperformed over the past 3–12 months continue to outperform over the next 3–12 months. This is inconsistent with both weak-form and semi-strong-form efficiency.

**Value premium** (Fama and French, 1992): stocks with high book-to-market ratios (value stocks) earn higher average returns than growth stocks, even after controlling for beta. Fama and French interpret this as a risk premium (value stocks are riskier in bad times); behavioralists attribute it to investor overreaction.

**Post-earnings announcement drift** (Ball and Brown, 1968): stock prices continue to drift in the direction of an earnings surprise for weeks after the announcement, inconsistent with immediate incorporation of information.

**Low-volatility anomaly** (Baker et al., 2011): low-volatility stocks earn higher risk-adjusted returns than high-volatility stocks — the opposite of what the CAPM predicts.

---

## 20.5 Financial Markets and the Macroeconomy

Financial market conditions affect the macroeconomy through several channels beyond the simple interest rate channel of IS–LM.

### The Wealth Effect

Household consumption depends partly on financial wealth (stocks, bonds, housing). The **wealth effect** on consumption is:

$$\frac{\partial C}{\partial W} = \frac{r}{1+r}\times\text{MPC out of wealth}.$$

From the life-cycle model (Chapter 11), the MPC out of wealth is approximately $r/(r+\delta_{age})$, where $\delta_{age}$ is the rate at which households plan to decumulate. Empirically, the marginal propensity to consume out of stock market wealth is approximately 3–5 cents per dollar (Case, Quigley, and Shiller, 2005); out of housing wealth, it is approximately 5–8 cents.

### The Balance Sheet Channel

When asset prices fall, the net worth of borrowers falls, raising the **external finance premium** (the wedge between the cost of external finance and the risk-free rate), which reduces investment — the financial accelerator of Chapter 24. The balance sheet channel works through:

$$\frac{\partial I}{\partial N_W} = \frac{\partial I}{\partial \rho}\cdot\frac{\partial \rho}{\partial N_W},$$

where $N_W$ is borrower net worth and $\rho$ is the external finance premium. A 10% fall in asset prices reduces net worth, raises the premium, and reduces investment. With leverage, the net worth effect is amplified: a 10% fall in asset values can eliminate 100% of equity for a firm with 90% debt financing.

### The Macro-Finance Interface

The integration of financial economics with macroeconomics — pricing assets with SDFs derived from macroeconomic models, and modeling macroeconomic dynamics with financial sector balance sheets as state variables — is one of the most active research frontiers in economics. Models in which the intermediary sector's leverage and equity are the SDF (He and Krishnamurthy, 2013; Brunnermeier and Sannikov, 2014) can simultaneously explain asset pricing patterns (time-varying risk premia, leverage cycles) and macroeconomic dynamics (amplification, slow recoveries). These models are surveyed in Chapter 39.

---

*Next: Chapter 21 — International Trade and Exchange Rates*
