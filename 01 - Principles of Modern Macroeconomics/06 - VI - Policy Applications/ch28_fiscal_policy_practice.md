# Chapter 28 — Fiscal Policy in Practice: Government Spending and Taxation

---

Chapter 22 developed the theoretical foundations of fiscal policy — the debt dynamics, Ricardian equivalence, and the fiscal multiplier. This chapter confronts those theories with evidence, examining what we actually know about how fiscal policy works and under what conditions it is effective. The gap between theory and evidence is large and illuminating: theoretical multipliers range from zero (complete Ricardian equivalence) to several times one (liquidity trap), and empirical multipliers exhibit comparably wide variation depending on the economic environment, the composition of the fiscal impulse, and the monetary policy response.

---

## 28.1 Automatic Versus Discretionary Fiscal Policy

Fiscal stabilization operates through two distinct mechanisms that must be carefully distinguished.

**Definition (Automatic Stabilizers).** **Automatic stabilizers** are changes in tax revenues and government transfers that occur automatically as the business cycle evolves, without any legislative action, cushioning income fluctuations. When GDP falls: income tax revenues fall (reducing the tax drain on household income), payroll taxes fall, corporate taxes fall, and unemployment insurance payments rise. This automatic fiscal relaxation partially offsets the demand decline. The reverse occurs in expansions.

**Definition (Discretionary Fiscal Policy).** **Discretionary fiscal policy** consists of deliberate legislative changes to spending programs or tax rates designed to affect aggregate demand. Stimulus packages, emergency relief legislation, infrastructure programs, and targeted tax cuts are all discretionary. Discretionary policy requires legislative action, which introduces lags: the **recognition lag** (time to identify that stabilization is needed), the **decision lag** (time to design and pass legislation), and the **implementation lag** (time for spending to actually flow into the economy).

The **cyclically adjusted primary surplus** $\hat{s}_t$ strips out the automatic stabilizer component from the observed primary surplus:

$$\hat{s}_t = s_t - \varepsilon^s \hat{x}_t,$$

where $\varepsilon^s$ is the **budget semi-elasticity** — the percentage point change in the primary surplus per percentage point of output gap. For the United States, $\varepsilon^s \approx 0.5$: a 2% negative output gap automatically reduces the primary surplus by about 1% of GDP. Changes in $\hat{s}_t$ represent discretionary policy; changes in $s_t - \hat{s}_t$ represent automatic stabilization.

---

## 28.2 Empirical Multiplier Estimates

The fiscal multiplier $\mu_G = \partial Y/\partial G$ cannot be estimated by OLS because fiscal policy is endogenous: policymakers tend to increase spending during recessions, generating a negative bias. The major identification strategies produce the following estimates.

### Narrative Identification

Romer and Romer (2010) construct narrative tax shocks from presidential tax proposals and Congressional Budget Office estimates, classifying changes as exogenous (motivated by ideology or structural concerns) versus endogenous (responding to the business cycle). Their estimated tax multiplier: a 1% of GDP exogenous tax cut raises output by approximately 3% over three years — a surprisingly large multiplier, substantially above what standard models predict.

Ramey (2011) uses military spending news (anticipated defense buildups identified from newspaper archives) to construct a government spending shock series. Her estimated spending multiplier: approximately 0.6–0.8 over a 5-year horizon, less than one. The contrast with the Romer-Romer tax multiplier reflects in part the different composition of fiscal impulses and in part the different identification assumptions.

### Cross-Sectional Evidence

Nakamura and Steinsson (2014) exploit cross-state variation in U.S. defense spending to estimate a **local multiplier** — the effect of a state-level spending increase on state-level output. Because monetary policy cannot respond to state-level shocks (it responds to national conditions), the local multiplier captures the direct effect of fiscal expansion without monetary policy offset. Their estimate: approximately 1.5 for the local multiplier. However, converting local to aggregate multipliers requires accounting for general equilibrium effects (demand spillovers across states reduce the aggregate multiplier relative to the local estimate).

Suarez Serrato and Wingender (2016) use discontinuities in the federal spending formula for states (based on demographic thresholds) to identify the causal effect of federal transfers on state income. Their estimates: approximately 1.7–2.1 — consistent with the Nakamura-Steinsson range.

### State-Dependence

A robust finding across multiple methodologies is that fiscal multipliers are larger during recessions than during expansions.

Auerbach and Gorodnichenko (2012) use a smooth-transition VAR that allows the multiplier to vary with the business cycle state. Their estimates: approximately 2.5 during recessions versus approximately 0.4 during expansions. The mechanism: in recessions, the economy operates below potential, the marginal propensity to consume of the unemployed and underemployed is high, monetary policy may be at or near the ELB (preventing crowding out), and the automatic stabilizer mechanism is in full operation.

---

## 28.3 The Composition of Fiscal Policy

The multiplier depends not only on the overall fiscal stance but on the composition of the fiscal impulse. Different spending and tax instruments have different multipliers for both theoretical and empirical reasons.

**Infrastructure investment** has a higher multiplier than government consumption (salaries of public employees) because: (i) it is more likely to be directed at constrained households and domestic supply chains; (ii) it raises potential output (the supply-side effects are positive), reducing the inflation pressure that would otherwise trigger monetary policy tightening; (iii) it crowds in private investment through network externalities and complementarities.

**Transfer payments to low-income households** have among the highest multipliers in the short run, because low-income households have high MPCs. Johnson, Parker, and Souleles (2006) estimate that 2001 U.S. tax rebates generated MPCs of approximately 0.25 in the first quarter and 0.60 over two quarters — consistent with a large fraction of liquidity-constrained households.

**Tax cuts for high-income households** tend to have lower short-run multipliers, since wealthy households have low MPCs. However, permanent corporate tax cuts that affect investment incentives may have larger medium-run supply-side effects through the user cost of capital channel.

**Government purchases vs. transfers**: the spending multiplier exceeds the transfer multiplier because every dollar of government purchases creates a direct dollar of demand, while a dollar of transfers is only partly consumed (the remainder is saved by non-constrained households).

---

## 28.4 Fiscal Multipliers at the Effective Lower Bound

The most important environment for fiscal policy is the ELB, where monetary policy cannot offset inflationary pressure from fiscal expansion through interest rate increases. The theoretical prediction (Christiano, Eichenbaum, and Rebelo, 2011) is that the ELB fiscal multiplier substantially exceeds one.

The mechanism in the NK model: a fiscal expansion raises the output gap, which raises inflation expectations, which reduces the real interest rate (since the nominal rate is pinned at the ELB), which stimulates private consumption and investment — a self-reinforcing loop that amplifies the initial fiscal impulse. The ELB multiplier:

$$\mu_G^{ELB} = \frac{1 - \frac{\sigma(1-\delta_G)\phi_y}{1-\rho}}{1 - \frac{\sigma\kappa(1-\delta_G)}{(1-\rho)(1-\beta\rho)} - \frac{\sigma(1-\delta_G)\phi_y}{1-\rho}},$$

where $\rho$ is the persistence of the ELB spell and $\delta_G$ the fraction of spending that is non-productive (transfers to non-constrained households). For standard parameters, this multiplier is approximately 1.5–3.0 — substantially above one.

Empirically, the evidence on ELB multipliers is limited by the infrequency of the event. The 2009 U.S. ARRA provides one observation: Romer and Bernstein (2009) projected a spending multiplier of 1.57 and a tax multiplier of 0.99 for the ELB environment. Subsequent structural VAR analyses (Cogan, Cwik, Taylor, and Wieland, 2010; Christiansen and Eichenbaum, 2012) produced a wide range of estimates (0.5 to 2.5), with the disagreement reflecting different identifying assumptions about the persistence of the ELB and the monetary policy reaction.

---

## 28.5 Fiscal Consolidation and Austerity

In the aftermath of the 2008–09 recession, many governments undertook significant fiscal consolidation — reducing spending and/or raising taxes to stabilize debt-to-GDP ratios. The debate about the macroeconomic effects of fiscal consolidation in a weak economy was contentious.

**Expansionary austerity**: Alesina and Ardagna (2010) documented historical episodes of "expansionary fiscal contractions" — cases where fiscal consolidations coincided with economic expansions. They attributed this to confidence effects (consolidation reduces the sovereign risk premium, lowering borrowing costs) and the composition of consolidation (spending cuts outperforming tax increases). Krugman and others argued this was selection bias: consolidations during expansions (when the economy is not constrained) naturally coincide with growth, while consolidations during recessions are contractionary.

**Blanchard and Leigh (2013)** conducted a systematic empirical test using the IMF's fiscal multiplier assumptions to construct a measure of forecast errors. They found: countries where the IMF assumed larger fiscal consolidations had larger output forecast misses — their growth was worse than the IMF projected. The implied multiplier from the forecast error regression: approximately 1.5, far above the pre-crisis IMF assumption of 0.5. This was a significant moment in the empirical fiscal policy literature, suggesting that austerity in post-crisis European countries was more costly than policymakers had anticipated.

---

*Next: Chapter 29 — Monetary Policy in Practice*
