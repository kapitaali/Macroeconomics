# Chapter 20 — Financial Markets: Stocks, Bonds, and Asset Prices

---

## 20.1 The Term Structure

Under the pure expectations theory, the $n$-period yield is the average of expected future short rates:

$$i_t^n = \frac{1}{n}\sum_{k=0}^{n-1}\mathbb{E}_t[i_{t+k}^1] + \phi_t^n,$$

where $\phi_t^n$ is the term premium for bearing duration risk. The yield curve inverts — a historically reliable recession predictor — when markets expect large future rate cuts.

---

## 20.2 The Stochastic Discount Factor

The SDF $M_{t+1}$ prices all assets via:

$$p_t = \mathbb{E}_t[M_{t+1}(p_{t+1} + d_{t+1})].$$

In the CCAPM with CRRA utility:

$$M_{t+1} = \beta\!\left(\frac{c_{t+1}}{c_t}\right)^{-\sigma}.$$

For any asset with return $R_{t+1}$:

$$\mathbb{E}_t[M_{t+1} R_{t+1}] = 1.$$

The risk premium on asset $j$:

$$\mathbb{E}_t[R_{t+1}^j] - R_t^f = -\frac{\mathrm{Cov}_t(M_{t+1},\, R_{t+1}^j)}{\mathbb{E}_t[M_{t+1}]}.$$

Assets that pay poorly in recessions (when $M_{t+1}$ is high) must offer a positive premium.

**The equity premium puzzle** (Mehra and Prescott, 1985): U.S. equity returns have historically exceeded risk-free returns by ~6% p.a., but the CCAPM with $\sigma \approx 1$ implies a premium below 0.5% — requiring $\sigma \approx 50$ to match the data.

---

## 20.3 Efficient Markets and Anomalies

The EMH (Fama, 1970) holds that prices fully reflect all available information. Documented anomalies: momentum, value premium, post-earnings announcement drift, short-run reversals. Behavioral finance attributes these to cognitive biases; rational finance attributes them to risk premia omitted from simple models.

---

## 20.4 Stock Prices and Monetary Policy

The Gordon growth model: $P_t^{eq} = D_{t+1}/(i^e - g)$. Monetary policy affects valuations through three channels: (i) lower $i^e$ raises the discount factor; (ii) fiscal policy may raise expected $g$; (iii) supply shocks affect dividends $D$. Shiller's CAPE ratio (10-year cyclically adjusted PE) forecasts subsequent real returns — consistent with time-varying risk premia rather than the pure expectations theory.

---

*Next: Chapter 21 — International Trade and Exchange Rates*
