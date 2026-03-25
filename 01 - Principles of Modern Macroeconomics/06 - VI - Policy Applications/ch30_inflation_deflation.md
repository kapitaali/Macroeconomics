# Chapter 30 — Inflation and Deflation: Causes, Consequences, and Cures

---

Inflation and deflation are the two failure modes of the price system. Persistent inflation erodes the purchasing power of money, distorts economic decisions, and — when severe enough — destroys the monetary economy itself. Persistent deflation raises real debt burdens, depresses demand through the paradox of postponement (consumers delay purchases anticipating lower future prices), and traps monetary policy at the ELB. This chapter examines the causes, consequences, and remedies for both phenomena.

---

## 30.1 The Costs of Inflation

### Steady, Anticipated Inflation

Under fully anticipated, constant inflation, economic agents can adjust all nominal contracts (wages, rents, bond coupons) to reflect the expected price path, leaving real variables unaffected. The only welfare cost is the **shoe-leather cost**: the cost of economizing on money holdings when the nominal interest rate (the cost of holding money) is high. Since $i = r + \pi$, higher inflation raises $i$, inducing households to hold less money and make more trips to the bank — hence "shoe-leather." The welfare triangle in the money demand diagram:

$$W_{SL} \approx \frac{1}{2}|\eta_i|\frac{M}{P}\pi^2,$$

where $|\eta_i|$ is the interest semi-elasticity of money demand. Lucas (2000) estimates this cost at approximately 1.0% of GDP for a move from 0% to 10% annual inflation — significant but not catastrophic.

A second cost of steady inflation is the **tax interaction effect** (Fischer, 1981): because nominal capital gains and interest income are taxed at the same rate as real income in most tax systems, inflation raises the effective tax rate on capital income even without any change in tax law. With inflation at 5% and a nominal return of 7%, the real return is 2% but taxes are levied on the full 7%, generating a distortion that reduces the incentive to save and invest.

### Unanticipated Inflation

Unanticipated inflation — departures of actual from expected inflation — has distributional consequences that are absent under anticipated inflation. Debtors gain (the real value of their fixed-rate nominal debt falls) and creditors lose. Holders of long-term nominal bonds lose; borrowers on fixed-rate mortgages gain. Workers with nominal wage contracts locked in at low rates lose real wages until renegotiation; firms gain at workers' expense.

These redistributive effects are not in themselves welfare losses (they are transfers), but they distort ex ante decisions: the risk of unexpected inflation makes lenders demand a premium to extend nominal credit (the **inflation risk premium**), shortening the maturity structure of debt markets and reducing investment.

### High Inflation and Hyperinflation

As inflation rises above approximately 10–20% annually, additional costs emerge. **Relative price variability** increases: in high inflation environments, firms adjust prices at different frequencies, generating idiosyncratic price movements that obscure relative price signals and misallocate resources (Barro, 1976). **Economic planning horizons shorten**: contracts in high-inflation economies are indexed to short maturities or to foreign currencies, fragmenting financial markets and reducing investment with long gestation periods.

**Hyperinflation** — Cagan's (1956) definition of monthly inflation exceeding 50% — represents the collapse of the monetary economy. Under hyperinflation, agents shift to barter, foreign currencies, or commodity money; the tax base of the inflation tax collapses as money demand falls to near zero; seigniorage revenue falls below the primary deficit even as the printing press runs at full speed.

The Cagan model of hyperinflation:

$$\ln(M_t/P_t) = -\alpha\pi_t^e + \text{const}, \quad \pi_t^e = \pi_{t-1}^e + \lambda(\pi_t - \pi_{t-1}^e),$$

where $\alpha > 0$ is the interest semi-elasticity of money demand and $\lambda \in (0,1)$ is the adaptive expectations adjustment coefficient. Cagan fitted this model to post-WWI German, Hungarian, Polish, and other hyperinflations and found $\alpha \approx 0.3$–$5.5$ across countries, with faster expectations adjustment in more severe episodes.

The fiscal root cause: hyperinflations universally involve a government that has lost access to conventional tax revenue (through war, occupation, or political collapse) and is financing its primary deficit entirely through seigniorage. The fiscal theory of inflation (Woodford, 2001; Cochrane, 2023) formalizes this: the price level is determined jointly by the fiscal and monetary policy regime, not by monetary policy alone. Stabilization of hyperinflations requires both monetary tightening (stopping money creation) and fiscal consolidation (closing the primary deficit) simultaneously — Sargent's (1982) lesson from the successful European stabilizations of the 1920s.

---

## 30.2 The Costs of Deflation

Deflation — a sustained decline in the general price level — is in several respects more dangerous than moderate inflation. Three mechanisms generate severe macroeconomic consequences.

### The Fisher Debt-Deflation Spiral

Irving Fisher (1933) described the mechanism by which unexpected deflation can generate a self-reinforcing economic contraction. The chain of events:

1. A negative demand shock causes the price level to fall.
2. Falling prices raise the real value of nominal debt: $D_{real} = D_{nominal}/P$ rises as $P$ falls.
3. Higher real debt burdens distress borrowers: households reduce consumption to service debt; firms reduce investment and increase asset sales.
4. Fire sales of assets depress asset prices further, reducing collateral values.
5. Credit conditions tighten (higher EFP), further reducing borrowing and investment.
6. Lower investment and consumption reduce aggregate demand, pushing prices down further — completing the spiral.

This mechanism was central to the Great Depression: the U.S. price level fell by approximately 25% between 1929 and 1933, raising the real value of private debt dramatically and contributing to the wave of bank failures and corporate bankruptcies. In a typical firm with debt-to-equity ratio $D/E$, a 10% fall in the price level raises the real debt burden by approximately $10\%\times D$ in absolute terms — wiping out equity entirely for a firm with leverage $D/E = 10$.

### The ELB Deflationary Trap

When the ELB binds, deflation raises the real interest rate automatically: $r_t = i_t^{ELB} - \pi_t > i_t^{ELB}$. As prices fall ($\pi_t < 0$), the real rate rises, further depressing investment and consumption — a self-reinforcing deflationary trap. The NK model at the ELB:

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i^{ELB} - \mathbb{E}_t[\pi_{t+1}] - r^n),$$

with $i^{ELB} = 0$ and $\mathbb{E}_t[\pi_{t+1}] < 0$ (deflationary expectations): $\hat{x}_t < 0$ is amplified by deflationary expectations, which raise the expected real rate further.

Eggertsson and Woodford (2003) showed that the optimal response to this trap is a credible commitment to keep rates at zero even after the economy recovers and inflation returns to target — "lower for longer" — thereby raising inflation expectations and reducing the real rate at the ELB. This is the theoretical foundation of AIT and forward guidance.

---

## 30.3 Disinflation: The Volcker Episode and its Lessons

The most dramatic disinflation in modern economic history was the **Volcker disinflation** of 1979–83. Paul Volcker became Fed chairman in August 1979, inheriting an inflation rate of 14.8% (Q3 1979). He shifted the Fed's operating procedure from targeting the federal funds rate to targeting non-borrowed reserves, allowing the funds rate to spike to 20% in June 1980 and again to 19% in mid-1981. The result: U.S. CPI inflation fell to 3.2% by the end of 1983, at the cost of two recessions — the brief but sharp 1980 recession and the deep 1981–82 recession in which unemployment reached 10.8%.

The sacrifice ratio for the Volcker disinflation: the cumulative output gap over 1980–83 was approximately 20 percentage points of GDP (using the Congressional Budget Office potential output estimates), and inflation fell by approximately 11 percentage points — giving a sacrifice ratio of approximately 1.8.

Ball (1994) systematically documents sacrifice ratios across 28 disinflation episodes in OECD countries during 1960–1990, finding values from approximately 0.5 to 4.0. Several patterns emerge:

- **Speed of disinflation**: faster disinflations have lower sacrifice ratios, consistent with the NKPC's forward-looking component (credible fast disinflation shifts expectations down quickly, reducing inflation before large output gaps are required).
- **Credibility**: countries with stronger central bank credibility and independence had lower sacrifice ratios.
- **Labor market institutions**: countries with more flexible wages and wages more responsive to unemployment had lower sacrifice ratios (consistent with the WS-PS framework of Chapter 19).

---

## 30.4 The 2021–22 Inflation Surge

The inflation surge of 2021–22 was the most significant test of the New Keynesian framework since the Volcker era. U.S. CPI inflation reached 9.1% year-on-year in June 2022, its highest since 1981. Euro area HICP inflation reached 10.6% in October 2022. The causes were multiple and interacting.

**Supply chain disruptions**: the COVID pandemic caused sharp disruptions in global supply chains — semiconductor shortages, container shipping bottlenecks, port congestion — that raised the relative prices of goods (especially durables) significantly. These are captured as cost-push shocks $u_t$ in the NKPC.

**Demand stimulus**: the massive fiscal transfers of 2020–21 (Chapter 41) created a large positive output gap in the goods sector even as the services sector remained depressed. The sectoral imbalance prevented normal price-adjustment mechanisms from equilibrating quickly.

**Energy price shocks**: the Russian invasion of Ukraine in February 2022 drove natural gas and oil prices sharply higher, generating a classic adverse supply shock.

**Delayed monetary response**: the Fed's AIT framework, combined with the view that inflation was "transitory," led to delayed rate increases — liftoff did not occur until March 2022. The subsequent tightening cycle was the fastest since 1980: 525 bps of rate increases over 16 months (March 2022–July 2023). Inflation subsequently fell — to approximately 3% by mid-2023 — but the episode raised fundamental questions about whether AIT was too permissive and whether the natural rate had shifted upward from its post-GFC lows.

---

*Next: Chapter 31 — Unemployment: Types, Policies, and Social Impact*
