# Chapter 30 — Inflation and Deflation: Causes, Consequences, and Cures

> *"Inflation is the one form of taxation that can be imposed without legislation."*
> — Milton Friedman, *Money and Economic Development*, 1973

---

Inflation and deflation are the two failure modes of the price system as a social institution. Persistent inflation erodes the real value of monetary contracts, distorts relative prices, and — when severe enough — destroys the monetary economy as a coordination device. Persistent deflation raises real debt burdens, produces a self-reinforcing contraction through the debt-deflation spiral, and traps monetary policy at the effective lower bound where it cannot provide further stimulus. Both represent departures from the price stability that allows money to function efficiently as a unit of account, medium of exchange, and store of value. This chapter examines the causes, costs, and remedies for each, drawing on both the theoretical frameworks of earlier chapters and the empirical record from the major inflationary and deflationary episodes of the past century.

---

## 30.1 The Costs of Inflation

The costs of inflation differ importantly depending on whether it is anticipated or unanticipated, and on whether it is moderate or extreme. A complete welfare analysis requires distinguishing these cases.

### Anticipated, Steady Inflation: Shoe-Leather and Distortions

Under fully anticipated, constant inflation, rational agents adjust all nominal contracts — wages, rents, bond coupons — to reflect the expected price path, leaving real variables unaffected in the absence of nominal rigidities. The only welfare costs under these ideal conditions are the **shoe-leather cost** and the **tax interaction effect**.

The shoe-leather cost arises because inflation raises the nominal interest rate ($i = r + \pi$ by the Fisher equation), which raises the opportunity cost of holding non-interest-bearing money balances. Households and firms economize on money balances by making more frequent transactions with financial intermediaries — metaphorically, wearing out their shoes on more trips to the bank. The welfare triangle in the money demand diagram gives:

$$W_{SL} \approx \frac{1}{2}|\eta_i|\left(\frac{M}{P}\right)\pi^2,$$

where $|\eta_i|$ is the interest semi-elasticity of money demand. Lucas (2000) calibrates this cost using the Baumol-Tobin model [Ch. 14] and estimates a welfare cost of approximately 1% of GDP for a shift from 0% to 10% annual inflation — significant, but not catastrophic.

The **tax interaction effect** (Fischer, 1981; Feldstein, 1997) is potentially larger. Most tax systems use nominal rather than real quantities as the tax base: capital gains are taxed on nominal appreciation; interest income is taxed at face value. Inflation therefore raises the effective tax rate on capital income even without any change in statutory rates. At 5% inflation and a 7% nominal return, the real return is 2% but taxes are levied on the full 7% nominal return, dramatically raising the effective tax wedge on saving and investment. Feldstein (1997) estimates this inflation-tax interaction raises the cost of a move from 2% to 0% inflation to approximately 1% of GDP annually — making low-inflation policy look substantially more beneficial than the shoe-leather cost alone implies.

### Unanticipated Inflation: Redistribution and Risk Premia

Unanticipated inflation transfers real wealth from creditors to debtors: fixed-rate nominal obligations decline in real value when inflation is higher than expected. Holders of long-term nominal bonds suffer capital losses; borrowers on fixed-rate mortgages gain. Wage earners locked into nominal contracts lose real wages until renegotiation; firms gain at workers' expense.

These redistributive effects are not in themselves efficiency losses (one party's loss is another's gain), but they distort ex ante decisions in two important ways. First, the risk of unexpected inflation forces lenders to demand an **inflation risk premium**, raising long-term nominal rates above what expected inflation and the real rate would imply and shortening the maturity structure of financial markets — reducing the availability of long-term finance for investment. Second, in the presence of debt contracts, unexpected inflation redistributes from households and firms with net nominal asset positions to those with net nominal liabilities, with macroeconomic effects depending on which group has the higher MPC — an asymmetry at the center of the HANK model [Ch. 25, Ch. 39].

### High Inflation: Relative Price Distortions and Coordination Failure

As inflation rises above approximately 10–20% annually, additional costs accumulate. **Relative price variability** increases: at high inflation, firms adjust prices at heterogeneous frequencies, generating large idiosyncratic price movements that obscure genuine changes in relative demand and misallocate resources. Barro (1976) documents the positive correlation between inflation and relative price variability across countries and time periods; this coordination cost is largely absent in low-inflation environments where price changes carry genuine information about relative scarcity.

At very high inflation, **planning horizons collapse**. Contracts shorten or index to the consumer price level or to foreign currencies. Long-term investment projects become unviable because the real cash flows over their lifespan cannot be reliably estimated. Financial intermediation shrinks or dollarizes. The economy retreats toward shorter production cycles and greater self-sufficiency, both of which reduce the gains from specialization.

---

## 30.2 Hyperinflation: When the Monetary Economy Breaks Down

**Definition (Hyperinflation).** Following Cagan (1956), **hyperinflation** is conventionally defined as a monthly inflation rate exceeding 50%, approximately equivalent to 13,000% annually. The most severe episodes far exceeded this threshold: German inflation in October 1923 reached $29,500\%$ per month; Hungarian inflation in July 1946 reached $41.9$ quintillion percent per month — the most extreme hyperinflation in recorded history.

Hyperinflation represents the collapse of the inflation tax revenue function. As inflation accelerates, agents economize aggressively on money balances — shifting to barter, foreign currencies, and commodity substitutes — causing real money demand to collapse. The **Cagan model** captures this dynamic:

$$\ln(M_t/P_t) = -\alpha\pi_t^e + \text{const}, \quad \pi_t^e = \pi_{t-1}^e + \lambda(\pi_t - \pi_{t-1}^e),$$

where $\alpha > 0$ is the interest (inflation expectation) semi-elasticity of money demand and $\lambda \in (0,1)$ governs the speed of expectations adjustment. Cagan estimated $\alpha \approx 0.3$–$5.5$ across seven post-WWI and post-WWII hyperinflations, with faster adjustment characterizing more severe episodes. As $\pi^e$ rises, real money demand $M/P$ falls, which means the seigniorage revenue $\pi \cdot M/P$ collapses despite the accelerating rate of money creation — the government falls off the right-hand side of the seigniorage Laffer curve.

### The Fiscal Root Cause

Hyperinflations universally involve governments that have lost access to conventional tax revenue and bond markets and are financing their entire primary deficit through money creation. The deeper cause is therefore fiscal, not monetary. Sargent (1982) examines the successful stabilizations of four European hyperinflations (Austria, Germany, Hungary, Poland) in the 1920s and documents that in every case, the hyperinflation ended abruptly — not gradually — coinciding with a credible fiscal reform: a large reduction in the primary deficit, often involving constitutional limits on central bank credit to the government, combined with the formation of a new government with a broader revenue base and credibility. Monetary tightening alone, without fiscal consolidation, did not stabilize any of them. The fiscal theory of the price level (Woodford, 2001; Cochrane, 2023) formalizes this: the price level is jointly determined by monetary and fiscal policy, and inflation ultimately reflects the government's ability to back its nominal liabilities with future fiscal surpluses.

---

## 30.3 Deflation and the Debt-Deflation Spiral

Deflation — a sustained decline in the general price level — is in several respects more dangerous than moderate inflation, particularly when the economy carries high nominal debt. Three mechanisms generate severe macroeconomic consequences.

### The Fisher Debt-Deflation Mechanism

Irving Fisher (1933), writing from the vantage point of the Great Depression's worst year, described the self-reinforcing mechanism through which unexpected deflation can cause economic collapse:

1. A negative shock depresses the price level.
2. Falling prices raise the real value of nominal debt: $D_{real} = D_{nominal}/P$ rises as $P$ falls.
3. Distressed debtors cut consumption and investment, and liquidate assets to reduce debt.
4. Asset liquidation depresses asset prices, further eroding collateral values.
5. Reduced collateral raises external finance premia, tightening credit conditions for all borrowers.
6. Tighter credit and lower investment reduce aggregate demand further, pushing prices down still more — completing the spiral.

The mechanism is self-reinforcing in a way that positive supply shocks are not. In a highly leveraged economy — as the U.S. was in 1929, with corporate debt-to-equity ratios near historical highs and household installment credit having tripled during the 1920s — even moderate deflation can trigger the spiral. The U.S. price level fell approximately 25% between 1929 and 1933, raising real debt burdens by approximately one-third for borrowers with fixed-rate obligations and contributing directly to the approximately 40% of U.S. banks that failed between 1930 and 1933.

### The ELB Deflationary Trap

When the effective lower bound binds and price growth turns negative, the real interest rate rises automatically:

$$r_t = i_t^{ELB} - \pi_t > i_t^{ELB} \quad \text{when } \pi_t < 0.$$

As prices fall and deflationary expectations become entrenched, the expected real rate can become severely elevated even with the nominal rate at zero. The NK-IS equation at the ELB:

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(0 - \mathbb{E}_t[\pi_{t+1}] - r_t^n)$$

with $\mathbb{E}_t[\pi_{t+1}] < 0$ implies $\hat{x}_t < 0$ is self-reinforcing: falling output reduces inflation further, which raises the real rate, which reduces output further. This is the **deflationary trap** that Japan experienced after 1998 and that policymakers feared during the Great Recession.

Eggertsson and Woodford (2003) proved that the optimal policy response to the deflationary trap is a credible commitment to keep rates at zero even after the economy recovers and inflation returns to target — accepting above-target inflation in the future to raise inflation expectations today and reduce the current real rate. This is the theoretical basis for AIT and for the "lower for longer" forward guidance strategies discussed in Chapter 29.

### The Paradox of Thrift

A third deflationary mechanism operates through aggregate demand: when households simultaneously increase saving (to rebuild balance sheets after an adverse wealth shock), their individual prudence reduces aggregate income, which reduces aggregate saving — the paradox of thrift identified by Keynes and formalized by Samuelson. At the macroeconomic level, an increase in the desired saving rate is self-defeating when monetary policy cannot lower the real rate to restore investment-saving balance: the equilibrating mechanism fails, and output falls instead of the real rate, reducing income until saving matches investment at a lower level of activity. This mechanism was prominently at work in 2008–10, when household deleveraging in the U.S., UK, and eurozone simultaneously reduced consumption demand across multiple large economies at the same time, preventing the normal trade-channel adjustment.

---

## 30.4 Disinflation: The Volcker Episode and the Sacrifice Ratio

The most dramatic disinflation in modern economic history — and the most important test of the Phillips curve framework — was the **Volcker disinflation** of 1979–1983.

Paul Volcker became Fed chairman in August 1979, inheriting CPI inflation of 14.8% (Q3 1979). He shifted operating procedures from targeting the federal funds rate to targeting non-borrowed reserves, allowing the funds rate to reach 20% in June 1980 and again in mid-1981 — the highest nominal rates in U.S. history. The results: two recessions (a sharp but brief contraction in early 1980 and a deeper recession in 1981–82, during which unemployment reached 10.8%); and a permanent reduction in inflation to approximately 3% by end-1983.

The **sacrifice ratio** for the Volcker disinflation, computing cumulative output gap over 1980–83 divided by the total inflation reduction, is approximately 1.8 — each percentage point of permanent disinflation required approximately 1.8 percentage points of cumulative lost output.

Ball (1994) systematically computed sacrifice ratios across 28 disinflation episodes in OECD countries from 1960–1990, finding values ranging from approximately 0.5 to 4.0. Several patterns emerge:

- **Speed matters**: faster disinflations have lower sacrifice ratios. The NKPC implies that credible, rapid disinflation shifts inflation expectations down quickly, reducing inflation before large output gaps have accumulated. Gradual disinflation allows adaptive expectations to maintain elevated inflation expectations for longer, requiring a sustained output gap.
- **Central bank credibility and independence**: countries with more independent central banks achieved disinflation at lower sacrifice ratios, consistent with the credibility theory's prediction that institutionally credible commitments reduce the inflation expectations' response lag.
- **Labor market flexibility**: countries with more responsive wages (lower nominal wage rigidity, more decentralized bargaining) had lower sacrifice ratios, consistent with the WS-PS model [Ch. 19] — when nominal wages adjust quickly to unemployment, real wages fall in recessions more rapidly, reducing the inflation pressure without requiring as large an output gap.

---

## 30.5 The 2021–22 Inflation Surge

The inflation surge of 2021–22 was the most significant test of the New Keynesian framework since the Volcker era. U.S. CPI inflation reached 9.1% year-on-year in June 2022; euro area HICP inflation reached 10.6% in October 2022. The episode tested central banks, models, and frameworks in ways that will shape the field for years.

### Multiple, Interacting Causes

The surge did not have a single cause but reflected the simultaneous interaction of several shocks.

**Supply chain disruptions**: COVID-19 caused profound disruptions to global supply chains — semiconductor shortages, shipping container imbalances, port congestion, and just-in-time logistics failures — that raised the relative prices of durable goods dramatically. Semiconductor shortages alone contributed an estimated 2–3 percentage points to vehicle price inflation, which was a major driver of core goods CPI in 2021.

**Sectoral demand imbalance**: the fiscal transfers of 2020–21 (Chapter 41) — approximately $\$5$ trillion in the United States alone — created extraordinary demand for goods at precisely the time when goods supply chains were disrupted, while services demand remained depressed by pandemic caution. This sectoral imbalance meant that demand and supply were simultaneously misallocated, preventing normal price equilibration.

**Energy and commodity shocks**: Russia's invasion of Ukraine in February 2022 caused an abrupt supply shock to European natural gas markets and global wheat and fertilizer markets, transmitting through cost-push channels into broad inflation [NKPC cost-push shock $u_t$].

**Delayed monetary response**: the Fed's AIT framework, combined with widespread forecasts that inflation was "transitory" (driven by identifiable supply bottlenecks that would self-resolve), led to delayed rate increases. Liftoff did not occur until March 2022, by which time annual inflation had exceeded 7% for several months and was still rising.

### The Tightening Cycle and Lessons

The subsequent tightening cycle — 525 basis points of rate increases between March 2022 and July 2023, the fastest pace since Volcker — succeeded in reducing inflation substantially without triggering the severe recession that many forecasters predicted. U.S. CPI inflation fell from 9.1% in June 2022 to approximately 3% by mid-2023, while unemployment remained below 4%. This "soft landing" outcome was better than historical precedent suggested was likely for such rapid tightening, and it raised productive debates about whether the Phillips curve had flattened further, whether supply-side resolution of bottlenecks was doing much of the disinflationary work, or whether the 2024 deceleration would prove temporary.

The episode's lessons for monetary policy frameworks are still being absorbed: whether AIT was the right framework in the first place; whether the "transitory" error reflected model failure, communication failure, or forecast failure; and whether the neutral rate $r^n$ had shifted upward from its post-GFC lows as fiscal expansions and investment-led demand for green transition financing raised the natural rate of interest.

---

## Chapter Summary

- **Anticipated inflation** costs arise primarily through shoe-leather costs $\approx \frac{1}{2}|\eta_i|(M/P)\pi^2$ and tax interaction effects on capital income — both significant but below 1% of GDP at moderate inflation rates.

- **Unanticipated inflation** redistributes from creditors to debtors and requires inflation risk premia in long-term markets; high inflation generates relative price distortions and planning horizon collapse.

- **Hyperinflation** (>50%/month) collapses real money demand and seigniorage revenue simultaneously; it is a fiscal phenomenon requiring coordinated monetary and fiscal reform to end, not monetary tightening alone.

- **Deflation** is dangerous through three mechanisms: the Fisher debt-deflation spiral (real debt rises as prices fall, forcing asset liquidation and tighter credit); the ELB deflationary trap ($r = i^{ELB} - \pi$ rises as $\pi$ falls below zero); and the paradox of thrift (simultaneous deleveraging across sectors reduces income rather than raising the saving rate).

- The **Volcker disinflation** achieved an 11-percentage-point inflation reduction at a sacrifice ratio of approximately 1.8, with lower ratios in faster, more credible, and institutionally independent disinflations (Ball, 1994).

- The **2021–22 inflation surge** reflected interacting supply chain, sectoral, energy, and monetary policy shocks; the subsequent rapid tightening achieved a partial soft landing, generating ongoing debates about AIT design, neutral rate shifts, and the role of supply-side resolution.

---

*Next: Chapter 31 — Unemployment: Types, Policies, and Social Impact*
