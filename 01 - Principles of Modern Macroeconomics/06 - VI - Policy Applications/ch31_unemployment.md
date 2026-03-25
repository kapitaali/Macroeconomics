# Chapter 31 — Unemployment: Types, Policies, and Social Impact

---

Unemployment is the most visible and personally damaging failure of the macroeconomy. A person who is unemployed loses income, skills, social connections, and psychological well-being — consequences that extend far beyond the immediate economic loss. Yet unemployment is not a monolithic phenomenon: it has multiple distinct causes, responds to different policy instruments, and carries different welfare implications depending on whether it is frictional, structural, or cyclical. This chapter develops a comprehensive taxonomy of unemployment, examines the major policy tools for reducing it, and discusses the social and psychological dimensions that aggregate statistics alone cannot capture.

---

## 31.1 A Taxonomy of Unemployment

**Definition (Frictional Unemployment).** **Frictional unemployment** is the unemployment that arises because matching workers to jobs takes time, even when aggregate demand is adequate and skills are appropriate. Workers who have just entered the labor force, who have voluntarily quit to seek better opportunities, or who have been laid off from a declining firm all contribute to frictional unemployment while they search. Frictional unemployment is partially efficient: some search time is socially valuable because it leads to better worker-job matches, increasing productivity.

**Definition (Structural Unemployment).** **Structural unemployment** arises from persistent mismatches between workers' skills or locations and employers' needs. Technology displacement (automation of routine tasks), geographic concentration of industry, trade-induced deindustrialization, and the mismatch between college and non-college labor markets all generate structural unemployment. Unlike frictional unemployment, structural unemployment does not resolve through faster job search — it requires retraining, relocation, or fundamental shifts in the labor market.

**Definition (Cyclical Unemployment).** **Cyclical unemployment** is the excess of actual unemployment above the natural rate due to a shortfall of aggregate demand — the component of $u_t$ in $u_t - u^*$. It is measured by the gap in the Beveridge curve (the vacancy rate relative to the unemployment rate) and by Okun's Law ($Y_t < \bar{Y}_t \implies u_t > u^*$). Unlike structural unemployment, cyclical unemployment should respond to macroeconomic stabilization policy.

The sum of frictional and structural unemployment constitutes the natural rate $u^*$. In practice, distinguishing the three types in real time is difficult: structural unemployment often appears first as cyclical, and protracted cyclical recessions can generate structural unemployment through hysteresis.

---

## 31.2 Hysteresis in Unemployment

**Definition (Hysteresis).** **Hysteresis** in the labor market is the mechanism by which cyclical unemployment becomes structural: a prolonged period of high unemployment raises the natural rate $u^*$ itself, so that full employment corresponds to a permanently higher unemployment rate than existed before the recession. Formally:

$$\dot{u}^* = \gamma(u_t - u^*_t), \quad \gamma > 0.$$

If $\gamma > 0$, a temporary negative shock that raises $u_t$ above $u^*$ will raise $u^*$ over time, even after the shock dissipates. The natural rate is not a fixed anchor but a slowly-moving variable that tracks the actual unemployment rate.

Hysteresis operates through several mechanisms. **Human capital atrophy**: workers who are unemployed for extended periods lose skills and work habits, reducing their productivity and attractiveness to employers. The relationship between unemployment duration and reemployment wages documents this: wages fall by approximately 1–3% for each month of unemployment, with the effect accelerating after 6 months (Jacobson, LaLonde, and Sullivan, 1993). **Employer stigma**: employers interpret long unemployment spells as negative signals, screening out long-term unemployed applicants even when the original unemployment was not the worker's fault (Kroft, Lange, and Notowidigdo, 2013). **Labor force withdrawal**: prolonged discouragement leads workers to exit the labor force permanently, shrinking the potential labor supply.

Blanchard and Summers (1986) documented strong hysteresis in European unemployment following the 1970s oil shocks: unemployment in the UK, France, Germany, and Spain rose sharply in 1975–76 and in 1981–82, and then remained stubbornly elevated for a decade rather than returning to pre-shock levels. This persistence was far greater than in the United States and was attributed to insider-outsider labor market dynamics: employed workers (insiders) set wages to protect their own employment without regard for the unemployed (outsiders), preventing the wage flexibility that would restore full employment.

The implication for stabilization policy: if hysteresis is significant, the cost-benefit analysis of aggressive counter-cyclical policy shifts dramatically. A temporary recession that raises unemployment by 2 pp for 2 years is more costly than it appears if it permanently raises the natural rate by 0.5 pp — generating a permanent reduction in potential output. The present value of this permanent output loss far exceeds the transitory output loss, making aggressive stabilization policy look much more attractive.

---

## 31.3 Unemployment Insurance: Theory and Optimal Design

Unemployment insurance (UI) serves two macroeconomic functions: it provides consumption insurance to displaced workers (the **insurance motive**) and it provides an automatic stabilizer that sustains aggregate demand during recessions. Against these benefits, UI reduces the incentive for unemployed workers to search intensively (**moral hazard**).

### The Baily–Chetty Formula

The optimal replacement rate $b/w$ balances the insurance benefit against the search-effort cost. Baily (1978) and Chetty (2006) derived a sufficient statistics formula:

$$\frac{b^*}{w} = \frac{\gamma_c}{\gamma_c + \varepsilon_{1-e,\,b}},$$

where:
- $\gamma_c = -(\Delta\ln c_{unemployed} - \Delta\ln c_{employed})/\Delta(b/w)$ is the **consumption drop risk aversion** — how much utility is lost from the consumption fall during unemployment, per unit of reduction in the replacement rate.
- $\varepsilon_{1-e,\,b} = \partial\ln(1-e)/\partial\ln b$ is the **moral hazard elasticity** — the percentage increase in unemployment duration per percentage increase in benefits.

The formula says the optimal replacement rate is higher when: (i) consumption drops severely during unemployment (high $\gamma_c$, high insurance value); (ii) workers have low search elasticity (low $\varepsilon_{1-e,b}$, low moral hazard cost).

Chetty (2008) estimates $\varepsilon_{1-e,b} \approx 0.5$: a 10% increase in UI benefits increases unemployment duration by approximately 5%. The estimated consumption drop during unemployment (approximately 10–20% for households exhausting UI) implies $\gamma_c \approx 1$–$2$. Plugging in: $b^*/w \approx 1/(1+0.5) = 0.67$ with $\gamma_c = 1$ — a replacement rate of approximately 67%, close to the U.S. UI replacement rate for median earners under normal conditions, and substantially below the 100%+ replacement rate provided temporarily during COVID.

### UI as Automatic Stabilizer

Beyond the individual insurance function, UI acts as a powerful automatic stabilizer because UI recipients have very high MPCs (they are liquidity-constrained and face urgent spending needs). Gruber (1997) estimates that UI prevents approximately 30% of the consumption fall that would otherwise occur during unemployment spells — a significant smoothing effect. Aggregate demand simulations (Chodorow-Reich et al., 2012) suggest that the U.S. extended UI benefits during the 2008–09 recession prevented approximately 440,000 additional job losses by sustaining consumer spending in high-unemployment states.

---

## 31.4 Active Labor Market Policies

**Passive** labor market policies (UI, early retirement, disability insurance) provide income replacement but do not actively facilitate reemployment. **Active** labor market policies (ALMPs) aim to reduce the structural natural rate by improving the matching between workers and employers.

The major ALMPs and their evidence base:

**Job search assistance and activation**: requiring UI recipients to actively search and providing placement services. Card, Kluve, and Weber (2018) meta-analyze 200 randomized evaluations and find that job search assistance programs have the highest average impact on employment rates, with effects typically appearing quickly (within 6 months) and persisting 2–3 years.

**Training programs**: retraining unemployed workers for growing sectors. The evidence is mixed: short-term training (a few weeks) has minimal effects; longer-term vocational training shows positive effects on wages and employment 2–5 years out, but the benefits depend critically on program quality and labor market conditions. OECD (2022) estimates a 5–10 percentage point increase in employment probability from high-quality training programs.

**Wage subsidies**: subsidizing employment of difficult-to-place workers. The main risk is **deadweight loss** (the employer would have hired the worker anyway) and **displacement** (the subsidized hire replaces an unsubsidized worker). Evidence from well-designed programs suggests subsidies can have positive net employment effects when targeted at long-term unemployed workers facing stigma.

**Geographic mobility assistance**: helping workers relocate to regions with better employment prospects. This policy is theoretically appealing (it addresses geographic mismatch) but has limited take-up in practice due to housing costs, family ties, and cultural barriers — the same frictions that limit OCA adjustment in the eurozone.

---

## 31.5 The Social and Psychological Costs of Unemployment

The welfare cost of unemployment extends far beyond the income loss that economic models typically capture. A substantial empirical literature — drawing on life satisfaction surveys, health data, and administrative records — documents large and persistent non-income effects of unemployment.

**Subjective well-being**: Clark and Oswald (1994) found that unemployment reduces self-reported life satisfaction by approximately 0.4 standard deviations — a larger effect than the equivalent income loss alone would predict, suggesting psychological costs beyond the wage effect. Unemployment appears to affect identity, social status, daily structure, and social connections in ways that income replacement cannot compensate for.

**Health effects**: Sullivan and von Wachter (2009) find that male workers who are displaced in mass layoffs (to separate the effect of unemployment from factors that may have caused the job loss) experience a 50–100% increase in mortality rates in the first year after job loss, with mortality remaining elevated for 15–20 years. Unemployment increases smoking, alcohol consumption, cardiovascular disease, and depression.

**Family and social effects**: unemployment increases divorce rates, reduces marriage rates, delays childbearing, and reduces children's educational attainment (through parental stress, income, and residential instability). These effects compound across generations, with children of unemployed parents having lower educational attainment and earnings as adults.

**Crime**: Raphael and Winter-Ebmer (2001) find that a 1-pp increase in the unemployment rate is associated with a 2–4% increase in property crime rates, consistent with the theory that crime is partly an economic substitution for legal employment when legal wages fall below criminal returns.

These social costs of unemployment — far exceeding the income losses that appear in standard welfare calculations — provide a strong argument for aggressive macroeconomic stabilization policy and well-designed ALMPs. They also highlight a limitation of models that evaluate welfare solely through GDP or consumption aggregates.

---

*Next: Chapter 32 — Open Economy Macroeconomics: Exchange Rates and Global Imbalances*
