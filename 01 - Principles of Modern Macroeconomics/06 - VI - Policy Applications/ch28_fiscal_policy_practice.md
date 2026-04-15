# Chapter 28 — Fiscal Policy in Practice: Government Spending and Taxation

> *"Practical men, who believe themselves to be quite exempt from any intellectual influence, are usually the slaves of some defunct economist."*
> — John Maynard Keynes, *The General Theory*, 1936

---

Chapter 22 developed the theoretical foundations of fiscal policy: the government's intertemporal budget constraint, the conditions under which Ricardian equivalence holds, and the analytical fiscal multiplier from the IS–LM and New Keynesian frameworks. This chapter confronts those theories with the empirical evidence. The gap between theory and evidence is both large and informative: theoretical multipliers range from zero under full Ricardian equivalence to several times unity at the effective lower bound, and empirical estimates span an almost equally wide range depending on the economic environment, the composition of the fiscal impulse, the monetary policy response, and the identification strategy. Understanding where along this range the truth lies — and why it depends so critically on context — is essential for evaluating any actual fiscal policy decision.

---

## 28.1 Automatic Versus Discretionary Fiscal Policy

Fiscal stabilization operates through two distinct mechanisms that must be carefully separated both conceptually and empirically.

**Definition (Automatic Stabilizers).** **Automatic stabilizers** are the changes in government revenues and expenditures that occur automatically as output and employment fluctuate, without any legislative action, providing counter-cyclical fiscal impulses that cushion income fluctuations. When GDP falls below potential: income tax revenues fall because taxable incomes fall; corporate tax revenues fall because profits contract; payroll tax contributions decline with the wage bill; and transfer payments — particularly unemployment insurance — rise as more workers become eligible. These automatic changes relax the fiscal stance precisely when private demand is deficient, without the legislative lags that afflict discretionary policy.

**Definition (Discretionary Fiscal Policy).** **Discretionary fiscal policy** consists of deliberate legislative changes to tax rates, spending programs, or transfer systems, designed to affect aggregate demand or supply. Fiscal stimulus packages, infrastructure appropriations, and emergency tax cuts are all discretionary. Three lags limit the stabilizing effectiveness of discretionary fiscal policy: the **recognition lag** (identifying that stabilization is needed, which itself requires several months of data), the **decision lag** (the legislative process, which in the U.S. averages 6–18 months for major legislation), and the **implementation lag** (the delay between legislative passage and actual spending flows, particularly for infrastructure investment where design, contracting, and construction phases may extend over several years).

### Measuring the Structural Fiscal Position

The **cyclically adjusted primary surplus** $\hat{s}_t$ removes the automatic stabilizer component from the observed fiscal balance to reveal the discretionary fiscal stance:

$$\hat{s}_t = s_t - \varepsilon^s \hat{x}_t,$$

where $\varepsilon^s$ is the **budget semi-elasticity** — the change in the primary surplus as a share of GDP per unit change in the output gap. For the United States, $\varepsilon^s \approx 0.5$: a 2-percentage-point negative output gap automatically reduces the primary surplus by approximately 1% of GDP through revenue shortfalls and spending increases. For European economies with stronger automatic stabilizers, $\varepsilon^s$ ranges from 0.4 to 0.6. Changes in $\hat{s}_t$ represent discretionary policy; changes in $s_t - \hat{s}_t$ represent the automatic stabilizer contribution — the component that requires no deliberate policy action.

---

## 28.2 Empirical Multiplier Estimates

The fiscal multiplier $\mu_G = dY/dG$ cannot be estimated by OLS because fiscal policy is endogenous: governments tend to increase spending during recessions and reduce it during booms, inducing a negative correlation between spending and output that biases multiplier estimates toward zero or negative values. Several identification strategies have been developed to overcome this endogeneity problem.

### Narrative Identification

Romer and Romer (2010) construct a narrative measure of exogenous tax changes, classifying each postwar U.S. tax change as either endogenous (responding to the state of the economy) or exogenous (motivated by a desire to raise long-run growth, reduce the deficit for structural reasons, or implement an ideological program), then using the exogenous tax changes as instruments. Their finding: a 1-percentage-point-of-GDP exogenous tax cut raises real GDP by approximately 3% over three years. This is a surprisingly large multiplier, substantially above what New Keynesian models calibrated to normal conditions predict, suggesting either significant non-linearities or model misspecification.

Ramey (2011) uses the anticipation of defense spending increases — identified from newspaper archives covering conflicts and congressional appropriations — to construct government spending shocks that are exogenous to current economic conditions. Her estimated spending multiplier: approximately 0.6–0.8 over a five-year horizon, substantially below the Romer-Romer tax multiplier. The contrast between these two estimates reflects both the different identification assumptions and the genuine difference between the multipliers for spending (which directly adds to aggregate demand but may crowd out private spending) and for tax cuts (which work through the consumption Euler equation and may be larger when households are liquidity-constrained).

### Cross-Sectional and Regional Evidence

Nakamura and Steinsson (2014) exploit cross-state variation in U.S. defense spending contracts to estimate a **local fiscal multiplier**. The key identification insight: U.S. monetary policy responds to national, not state-level, conditions. A defense buildup that primarily benefits Texas (because of its concentration of defense contractors) generates a fiscal stimulus for Texas without triggering a monetary policy response specific to Texas — the ELB or the Taylor rule responds to the national cycle, not the state cycle. This setting therefore identifies the fiscal multiplier holding monetary policy constant. Their point estimate: approximately 1.5 for the local spending multiplier.

The local multiplier, however, overstates the aggregate multiplier because local spending spillovers are not captured (spending in Texas generates income that is partly spent on goods from other states, but this general equilibrium effect does not appear in the Texas GDP figures). Accounting for demand leakages across states reduces the implied aggregate multiplier to approximately 1.0–1.2.

Suarez Serrato and Wingender (2016) use a different source of regional variation: discontinuities in federal formula grants to states and local governments, driven by demographic thresholds that determine funding eligibility. Their estimates of the income multiplier for federal transfers: approximately 1.7–2.1 — consistent with the Nakamura-Steinsson range and suggesting significant short-run multipliers for federal spending at the local level.

### State-Dependence: When the Multiplier Is Large

A robust finding across multiple identification strategies is that fiscal multipliers vary substantially with the state of the business cycle, with the monetary policy stance, and with the degree of household financial constraint.

Auerbach and Gorodnichenko (2012) use a smooth-transition VAR that allows the fiscal multiplier to vary continuously with a measure of the business cycle state (the unemployment gap). Their estimates: approximately 2.5 during recessions versus approximately 0.4 during expansions. The mechanisms are theoretically coherent: in recessions, (i) the economy operates below potential, so fiscal expansion meets underutilized resources rather than crowding out private activity; (ii) the fraction of liquidity-constrained households is higher; (iii) monetary policy is more likely to be at or near the ELB, preventing interest rate increases that would otherwise crowd out private investment; and (iv) the automatic stabilizer feedback is in full operation, reducing the net fiscal cost of each dollar of stimulus.

---

## 28.3 The Composition of Fiscal Policy

The fiscal multiplier depends not only on the magnitude of the fiscal impulse but critically on its composition. Different spending and tax instruments work through different channels with different speeds, magnitudes, and distributional effects.

### Spending Multipliers by Category

**Public infrastructure investment** tends to have the highest spending multiplier among government expenditure categories. The reasons are several: infrastructure investment is directed at domestic contractors and construction workers who are more likely to face liquidity constraints than average households; it generates positive supply-side effects by raising potential output and the marginal product of private capital (through network externalities and complementarities between public and private capital); and it does not crowd out private investment through competition for the same capital goods as government consumption.

**Government consumption** (salaries of public-sector workers, purchases of goods and services) has an intermediate multiplier. It generates direct demand but may crowd out private activity in markets for skilled workers and specific goods, and may substitute for private consumption in ways that infrastructure does not.

**Transfer payments** (unemployment insurance, food stamps, cash transfers) have a high short-run multiplier when directed at liquidity-constrained households, because such households have high MPCs and spend the transfers quickly. Johnson, Parker, and Souleles (2006) estimate that households spent approximately 25 cents of each dollar from the 2001 U.S. tax rebate in the first quarter and approximately 60 cents over two quarters — evidence of significant short-run stimulus from targeted transfers.

**Tax cuts for high-income households** have lower short-run multipliers, since wealthy households save more (low MPC) and are less likely to face liquidity constraints. However, permanent corporate tax cuts affecting investment incentives may have larger medium-run supply-side effects through the user cost of capital channel and TFP growth, effects that are not captured in the multiplier's short-run output focus.

### The Transfer versus Spending Distinction

An important but often overlooked distinction: every dollar of government **purchases** creates a direct dollar of aggregate demand (the government buys goods and services, generating income for the sellers). A dollar of **transfers** creates only $MPC \times \$1$ of direct demand, because the recipient may save part of it. This implies:

$$\mu_G > \mu_T > 0,$$

where $\mu_G$ is the government purchases multiplier and $\mu_T$ is the transfer multiplier. The spending multiplier always exceeds the transfer multiplier under standard assumptions. For large transfers to liquidity-constrained households, however, the difference narrows because $MPC \to 1$ for constrained households. This compositional arithmetic is central to designing stimulus packages efficiently.

---

## 28.4 Fiscal Multipliers at the Effective Lower Bound

The most consequential environment for fiscal policy is the ELB, where monetary policy cannot offset the inflationary pressure from fiscal expansion through interest rate increases. The theoretical prediction from the New Keynesian model (Christiano, Eichenbaum, and Rebelo, 2011) is that the ELB multiplier substantially exceeds one.

### The NK Mechanism

Under normal conditions, a fiscal expansion raises the output gap and inflation, prompting the central bank to raise the nominal interest rate under the Taylor rule. Higher rates crowd out private investment and durable goods consumption, partially offsetting the fiscal stimulus. At the ELB, this monetary offset is absent: the nominal rate cannot rise, so the real rate falls as inflation expectations increase, further stimulating private demand. The self-reinforcing loop amplifies the initial fiscal impulse substantially.

The ELB spending multiplier in the calibrated NK model:

$$\mu_G^{ELB} = \frac{1}{1 - \frac{\sigma\kappa(1-\delta_G)}{(1-\rho)(1-\beta\rho)} - \frac{\sigma(1-\delta_G)\phi_y}{1-\rho}},$$

where $\rho$ is the persistence of the ELB spell, $\delta_G$ is the fraction of spending on transfers to non-constrained households, and $\kappa$ is the NKPC slope. For standard parameter values and a two-year ELB spell, this yields multipliers of approximately 1.5–3.0 — substantially exceeding the normal-times estimate.

### The ARRA Natural Experiment

The 2009 American Recovery and Reinvestment Act (ARRA), at approximately $\$787$ billion over two years, provides the most studied ELB fiscal experiment. Wilson (2012) uses cross-state variation in ARRA spending to estimate a local multiplier of approximately 1.5 for the spending components. Chodorow-Reich et al. (2012) find that each dollar of Medicaid grants to states (which had to be matched by state spending) raised employment substantially more than non-Medicaid grants — consistent with the higher multiplier for spending categories with low leakage and high employment intensity.

The disagreement about the ARRA's aggregate effects largely reflects disagreement about how to attribute post-ARRA output growth — whether the U.S. recovery after 2009 was fast or slow relative to the appropriate counterfactual, and whether the counterfactual should be defined by DSGE models or by the historical experience of post-financial-crisis recoveries (Reinhart and Rogoff, 2009, suggest the latter are systematically slower).

---

## 28.5 Fiscal Consolidation and Austerity

Following the 2008–09 recession, many governments undertook significant fiscal consolidation — cutting spending and/or raising taxes to stabilize rising debt-to-GDP ratios. The macroeconomic debate about the effects of this consolidation was one of the most contentious in recent economic history.

### The Expansionary Austerity Hypothesis

Alesina and Ardagna (2010) surveyed historical fiscal consolidations across OECD countries and identified episodes of "expansionary fiscal contractions" — cases where significant consolidation coincided with economic expansion rather than contraction. They attributed this to two channels: confidence effects (consolidation reduces sovereign risk premia, lowering long-term borrowing costs for the private sector) and composition effects (consolidations dominated by spending cuts outperform those dominated by tax increases, because spending cuts are more credible and have smaller demand effects).

### The IMF Forecast Errors Study

Blanchard and Leigh (2013) conducted a systematic empirical test that became one of the most influential applied macro papers of the post-crisis decade. Their method: using the IMF's own fiscal multiplier assumptions at the time (roughly 0.5), they constructed a measure of the unanticipated fiscal consolidation across European countries. They then regressed IMF forecast errors (actual growth minus IMF forecast) against these unanticipated consolidations. The key finding: countries with larger-than-expected consolidations consistently grew less than the IMF projected. The implied fiscal multiplier from the regression: approximately 1.5 — three times the IMF's assumption. The interpretation: European austerity after 2010 was substantially more contractionary than policymakers expected, partly because the consolidations occurred in a low-growth, ELB environment and partly because simultaneous consolidations across multiple countries eliminated the trade-channel offset that would normally allow a consolidating economy to expand net exports.

The Blanchard-Leigh study contributed to a significant shift in IMF thinking about the appropriate pace of fiscal consolidation in weak economic conditions — a shift formalized in the IMF's 2012 World Economic Outlook, which explicitly revised upward its estimate of fiscal multipliers for depressed economies.

---

## Chapter Summary

- **Automatic stabilizers** cushion income fluctuations without legislative lags; the cyclically adjusted primary surplus $\hat{s}_t = s_t - \varepsilon^s\hat{x}_t$ isolates the discretionary fiscal stance from the automatic cycle. For the U.S., $\varepsilon^s \approx 0.5$.

- **Empirical multiplier estimates** range from approximately 0.6 (Ramey's defense spending) to 3 (Romer-Romer tax changes), with the range reflecting genuine heterogeneity across instruments, environments, and identification strategies.

- The **local multiplier** approach (Nakamura-Steinsson) circumvents the monetary policy endogeneity problem; estimated local multipliers of approximately 1.5 translate to aggregate multipliers of approximately 1.0–1.2 after accounting for demand spillovers.

- **Multipliers are state-dependent**: approximately 2.5 during recessions versus 0.4 during expansions (Auerbach-Gorodnichenko), driven by slack, liquidity constraints, and the monetary offset mechanism.

- **Composition matters**: infrastructure investment and targeted transfers to liquidity-constrained households have higher short-run multipliers than broad tax cuts or government consumption; spending always exceeds the transfer multiplier.

- The **ELB multiplier** exceeds one in NK models and the ARRA evidence; the **Blanchard-Leigh study** found European austerity was three times more contractionary than the IMF's pre-crisis multiplier assumptions implied, contributing to an important revision in consensus views.

---

*Next: Chapter 29 — Monetary Policy in Practice: Interest Rates and Quantitative Easing*
