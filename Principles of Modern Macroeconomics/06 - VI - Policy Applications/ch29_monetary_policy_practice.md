# Chapter 29 — Monetary Policy in Practice: Interest Rates and Quantitative Easing

---

Chapter 23 developed the theory of optimal monetary policy in the New Keynesian framework — the ideal central bank, solving a welfare optimization problem with perfect information and full commitment. Real-world monetary policy involves none of these luxuries: information is incomplete and revised continuously, models are uncertain and contested, and commitment mechanisms are imperfect. This chapter examines how monetary policy is actually conducted, what we know about its quantitative effects, and how the toolkit expanded after 2008 when conventional interest rate policy reached its limits.

---

## 29.1 The Taylor Rule in Practice and its Variants

The Taylor (1993) rule described in Chapter 23 has become the benchmark against which actual monetary policy is evaluated. Empirically, the fed funds rate tracks the Taylor rule reasonably well in "normal" times but deviates significantly at turning points and at the ELB. Orphanides (2001) showed that the apparent deviations of 1970s Fed policy from the Taylor principle disappear when real-time data (the data available to policymakers at the time) are used instead of revised data — the Fed may have been following the Taylor principle relative to its own incomplete information.

The **neutral rate estimation problem** is a major practical challenge. The Taylor rule requires the neutral real rate $r^n$, but $r^n$ is not directly observable — it must be estimated from the data. Estimates diverge widely. The Laubach–Williams (2003) state-space model estimates $r^n$ declined from approximately 3.5% in the early 1990s to approximately $-0.5\%$ by 2015–16, a staggering 4-percentage-point decline. If $r^n$ is near zero or negative, the appropriate policy rate under the Taylor rule is also near zero even with inflation at target and output at potential — implying that conventional monetary policy has essentially no "room" at all in the new normal.

---

## 29.2 Inflation Targeting: Design and Empirical Performance

Inflation targeting was introduced by New Zealand in 1990. The framework has four standard features: (i) a public commitment to a specific inflation target (typically 2% CPI); (ii) central bank operational (instrument) independence; (iii) transparency through quarterly Inflation Reports, regular press conferences, and published forecasts; and (iv) accountability through parliamentary testimony and open letters when the target is missed.

**Credibility and the sacrifice ratio**: one of the main claimed benefits of inflation targeting is enhanced credibility — markets and wage-setters believe the central bank will maintain low inflation, reducing the sacrifice ratio for any necessary disinflation. The evidence: Ball and Sheridan (2005) found that inflation fell similarly in OECD targeting and non-targeting countries during the 1990s, challenging the hypothesis that adoption per se reduced inflation. Rose (2007) found targeting reduced inflation level and variance, but identifying causation is difficult because countries often adopted targeting precisely when inflation was high and they were committing to credible disinflation.

**Average Inflation Targeting (AIT)**: the Federal Reserve's August 2020 framework revision, discussed in Chapter 23. Key features: (i) outcome-based forward guidance; (ii) a "make-up" strategy allowing temporary above-target inflation following below-target periods; (iii) an asymmetric labor market goal (maximum employment as a broad, inclusive objective with no fixed numerical target). AIT was designed to raise inflation expectations and reduce the probability of the ELB binding persistently, but its practical implementation during the 2021–22 inflation surge — when delayed tightening allowed inflation to reach 8% — raised questions about whether the framework had been interpreted too loosely.

---

## 29.3 Quantitative Easing: Mechanisms and Evidence

When the policy rate is at the ELB, the central bank's standard tool is exhausted. **Quantitative easing (QE)** extends the toolkit by having the central bank purchase large quantities of longer-term assets.

### Preferred Habitat Theory and Portfolio Balance

The standard theoretical framework for QE's interest rate effects is the **preferred habitat model** of Vayanos and Vila (2009). In this model, investors have preferred maturities (pension funds want 30-year bonds; money market funds want overnight); they will deviate from their preferred maturities only if sufficiently compensated by term premia. When the central bank removes long-duration bonds from the market (by buying them), preferred-habitat investors must hold more short-duration bonds, bidding down long-duration yields relative to short-duration yields — reducing the term premium $\phi_t^n$.

The term premium effect can be written as:

$$\Delta\phi_t^n \approx -\frac{\Delta Q_t}{D_t},$$

where $\Delta Q_t$ is the change in central bank holdings of long-duration bonds and $D_t$ is the dollar duration of the market portfolio. Larger purchases relative to outstanding supply generate larger term premium reductions.

### Empirical Evidence

**Federal Reserve QE1** (November 2008–March 2009, $\$600$ billion in MBS and agency debt; $\$300$ billion in Treasuries): Gagnon et al. (2011) estimate the effect of the full program on 10-year Treasury yields at approximately $-90$ to $-100$ basis points using event studies and time-series regressions. The effect on 30-year MBS yields was similar, reducing mortgage rates and supporting the housing market.

**Federal Reserve QE2** (November 2010–June 2011, $\$600$ billion in Treasuries): D'Amico and King (2013) estimate a term premium reduction of approximately $-26$ bps per $\$100$ billion, consistent with the portfolio balance channel. Hamilton and Wu (2012) estimate similar magnitudes using a term structure model.

**ECB PSPP (Quantitative Easing)** (2015–2018, €2.6 trillion): Altavilla, Carboni, and Motto (2015) estimate 10-year sovereign yield reductions of approximately $-0.5$ to $-0.7$ pp across eurozone countries in the first year, with larger effects for peripheral countries (Italy, Spain) where sovereign spreads were also compressed.

The **transmission from long rates to the real economy** is the critical and less certain step. Lower long-term yields reduce borrowing costs for firms and households, stimulate investment and durable goods purchases, depreciate the currency (stimulating net exports), and raise asset prices (through the wealth effect). The evidence suggests meaningful but modest effects on real activity: Weale and Wieladek (2016) estimate that the Bank of England's QE programs raised UK GDP by approximately 1.5% and inflation by approximately 1.25 pp relative to the no-QE counterfactual.

---

## 29.4 Forward Guidance: Credibility and the Paradoxes

Forward guidance is theoretically the most powerful unconventional tool: since the NK IS curve depends on the full expected future path of real rates, a credible commitment to keep rates low for longer can reduce current long rates without any immediate balance sheet expansion. The NK forward solution:

$$\hat{x}_t = -\sigma\sum_{k=0}^\infty\mathbb{E}_t[i_{t+k} - \pi_{t+k+1} - r_{t+k}^n].$$

A 50-basis-point reduction in expected short rates over the next 5 years reduces the sum by $0.005 \times 20 = 0.10$, generating a 10% output gap increase with $\sigma = 1$ — a large effect.

The **forward-guidance puzzle** (Del Negro, Giannoni, and Patterson, 2015): in the standard NK model, forward guidance is implausibly powerful. A credible announcement that rates will be zero for the next 5 years instead of 2 generates enormous effects on current output and inflation — far larger than VAR evidence supports. The behavioral NK model (Chapter 15) with cognitive discounting $M < 1$ attenuates these effects by making agents discount future policy commitments.

**Credibility constraints** also limit forward guidance. Date-based guidance ("rates will remain at zero through mid-2015") is straightforward to communicate but creates a credibility problem: if the economy recovers faster than expected, the central bank faces pressure to raise rates before the stated date, undermining the commitment. Outcome-based guidance ("rates will remain at zero until unemployment falls below 6.5%") avoids this problem by conditioning on economic outcomes, but is harder to communicate and leaves more uncertainty about the rate path.

---

## 29.5 Monetary Policy and Financial Stability

A long-running debate concerns whether central banks should use interest rate policy to "lean against" financial imbalances — to raise rates above the level dictated by the Taylor rule when credit growth and asset prices are elevated, to reduce the probability of a financial crisis. The debate involves the trade-off between higher current output and inflation variance (the cost of leaning) against a lower probability of future financial crises (the benefit).

**The case against leaning** (Svensson, 2017): the costs of using interest rates to deflate credit and asset bubbles are concentrated in the near term (lower output and employment) while the benefits are diffuse and uncertain (reduced but not eliminated crisis probability, with a long and uncertain lag). Macroprudential tools (countercyclical capital buffers, LTV limits, leverage constraints) are better targeted at financial stability risks without the broad macroeconomic side effects.

**The case for leaning** (Borio and White, 2004): financial imbalances build slowly and are not well captured by current inflation and output gaps; a central bank that ignores them may face much larger eventual instability costs. The 2008 crisis demonstrated that the "clean up after" strategy (wait for the bubble to burst, then cut rates) has limits when the ELB binds and the recession is severe.

The emerging consensus (IMF, 2015; Svensson, 2017 rebuttal): macroprudential policy should be the primary tool for financial stability; monetary policy should incorporate financial stability considerations only when macroprudential tools are unavailable or exhausted. This consensus is embodied in the "two-pillar" frameworks adopted by many central banks, which maintain separate inflation and financial stability mandates.

---

*Next: Chapter 30 — Inflation and Deflation*
