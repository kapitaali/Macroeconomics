# Chapter 31 — Unemployment: Types, Policies, and Social Impact

> *"Of all the aspects of social misery nothing is so heartbreaking as unemployment."*
> — Jane Addams, *Democracy and Social Ethics*, 1902

---

Unemployment is the most visible and personally damaging failure of the macroeconomy. A person who loses their job loses income, skills, social connections, and — as the empirical literature documents with increasing precision — psychological wellbeing and physical health in ways that extend far beyond the immediate financial loss. Yet unemployment is not a single phenomenon. It has multiple distinct causes, responds to fundamentally different policy instruments, and carries different welfare implications depending on whether it is frictional, structural, or cyclical. This chapter develops a systematic taxonomy of unemployment, examines the major policies for reducing each type, and discusses the social and psychological costs that aggregate statistics alone cannot capture.

---

## 31.1 A Taxonomy of Unemployment

**Definition (Frictional Unemployment).** **Frictional unemployment** is the unemployment that exists because matching workers to appropriate jobs takes time, even when aggregate demand is adequate and the aggregate supply of skills broadly matches the aggregate demand for skills. Workers entering the labor market for the first time, voluntarily quitting to find better-matched positions, or recently laid off from firms in declining sectors all contribute to frictional unemployment while they search for new jobs. Frictional unemployment is partially efficient: some search time generates better worker-job matches that raise long-run productivity; a world of zero frictional unemployment would be one in which workers were instantly reallocated regardless of match quality.

**Definition (Structural Unemployment).** **Structural unemployment** arises from persistent mismatches between the characteristics of unemployed workers — their skills, location, or occupational experience — and the characteristics demanded by available vacancies. Technology displacement (the automation of routine cognitive and manual tasks), trade-induced deindustrialization, geographic concentration of growth in skill-intensive metropolitan areas, and the collapse of middle-skill employment have all generated structural unemployment that cannot be resolved by faster job search alone.

**Definition (Cyclical Unemployment).** **Cyclical unemployment** is the component of observed unemployment above the natural rate attributable to insufficient aggregate demand — the quantity $u_t - u^*$ when $Y_t < \bar{Y}_t$. It is directly linked to the output gap through Okun's Law [Ch. 3]: $u_t - u^* \approx -\psi(Y_t/\bar{Y}_t - 1)$ with $\psi \approx 0.5$. Unlike structural unemployment, cyclical unemployment should respond to demand management — monetary accommodation, fiscal stimulus — that restores aggregate demand to potential.

The sum of frictional and structural unemployment constitutes the **natural rate** $u^*$. In practice, distinguishing the three types in real time is both essential for policy design and extremely difficult empirically. Structural unemployment that has persisted long enough can appear cyclical; cyclical recessions that are deep and prolonged can create structural unemployment through hysteresis — the mechanism by which temporary cyclical shocks become permanent features of the labor market.

---

## 31.2 Hysteresis in Unemployment

**Definition (Hysteresis).** **Hysteresis** in the labor market refers to the process by which a period of high cyclical unemployment raises the natural rate $u^*$ itself, so that after the shock has passed, full employment corresponds to a permanently higher unemployment rate than existed before. The formal condition:

$$\dot{u}^* = \gamma(u_t - u^*_t), \quad \gamma > 0.$$

When $\gamma > 0$, a transitory recession that raises $u_t$ above $u^*$ gradually raises $u^*$ as well, even after the original negative shock has dissipated. The natural rate becomes a slowly-moving endogenous variable that tracks the actual unemployment rate with a lag, rather than an exogenous anchor to which the economy returns.

### Mechanisms

Hysteresis operates through three main channels.

**Human capital atrophy**: extended unemployment spells reduce the skills and work habits of displaced workers. Employers observe this deterioration — or rationally anticipate it — and prefer workers with more recent work experience, reducing the reemployment probability of long-term unemployed workers. The empirical relationship between unemployment duration and subsequent wages confirms this: each month of unemployment reduces subsequent wages by approximately 1–3%, with the effect accelerating after six months (Jacobson, LaLonde, and Sullivan, 1993).

**Employer stigma**: even when the original job loss was entirely the employer's fault (a plant closure, a mass layoff), employers use unemployment spells as a signal of worker quality. Kroft, Lange, and Notowidigdo (2013) document this through a large-scale audit study: equivalent resumes with longer unemployment gaps receive significantly fewer callbacks, even when the callback rate is controlled for the business cycle state at the time of the gap. The stigma effect is stronger in labor markets with more applicants per vacancy, suggesting it reflects statistical discrimination rather than pure taste-based discrimination.

**Labor force withdrawal**: prolonged unemployment leads discouraged workers to exit the labor force permanently. Workforce exit is particularly pronounced among older workers and those with weak labor force attachment; once they exit, labor force re-entry rates are very low. This mechanism shrinks the effective labor supply permanently, raising the natural rate and reducing potential output.

### Blanchard and Summers on European Unemployment

Blanchard and Summers (1986) documented severe hysteresis in European labor markets following the 1970s oil shocks. UK, French, German, and Spanish unemployment all rose sharply in 1975–76 (following the first oil shock) and again in 1981–82 (following the second shock) and then remained persistently elevated for a decade or more — far longer than in the U.S., where unemployment showed much less persistence. They attributed European hysteresis to **insider-outsider dynamics**: employed workers (insiders) bargain over wages to protect their own employment without regard for the unemployed (outsiders). When employment falls, insiders bargain for wages that sustain full employment of the remaining insiders — preventing the wage flexibility that would otherwise reduce unemployment toward its pre-shock level.

### Policy Implications of Hysteresis

If hysteresis is significant, the welfare calculus of counter-cyclical policy shifts dramatically. Consider a recession that raises unemployment by 2 percentage points for two years. Without hysteresis, the welfare cost is the integral of output gap losses over those two years. With hysteresis (say, $\gamma = 0.1$), the natural rate rises by approximately 0.4 percentage points permanently — generating a permanent stream of lost potential output whose present value far exceeds the transitory loss. This implies that aggressive stabilization policy — aiming for rapid return to full employment — has much larger expected benefits than conventional estimates suggest. It also implies that premature withdrawal of stimulus (fiscal consolidation before full recovery) is much more costly when hysteresis is present.

---

## 31.3 Unemployment Insurance: Theory and Optimal Design

Unemployment insurance (UI) serves two macroeconomic functions: it provides consumption insurance to displaced workers (smoothing consumption during the unemployment spell) and it sustains aggregate demand during recessions through its automatic stabilizer properties. Against these benefits, UI reduces the incentive for unemployed workers to search intensively — the **moral hazard** that is the central tradeoff in UI design.

### The Baily-Chetty Formula

The optimal replacement rate $b^*/w$ (the ratio of the benefit level to the pre-unemployment wage) balances the marginal insurance benefit against the marginal moral hazard cost. Baily (1978) and Chetty (2006) derived a sufficient statistics formula that expresses the optimum in terms of directly estimable quantities:

$$\frac{b^*}{w} = \frac{\gamma_c}{\gamma_c + \varepsilon_{1-e,\,b}},$$

where:

- $\gamma_c$ is the **consumption smoothing benefit** — a measure of the concavity of utility relative to the consumption drop during unemployment. It is higher when households face larger consumption drops during unemployment and when risk aversion is greater.
- $\varepsilon_{1-e,\,b} = \partial\ln(1-e)/\partial\ln b$ is the **moral hazard elasticity** — the percentage change in unemployment duration per percentage increase in UI benefits.

Chetty (2008) estimates $\varepsilon_{1-e,b} \approx 0.5$ from U.S. variation in UI generosity across states: a 10% increase in UI benefits extends the average unemployment spell by approximately 5%. The estimated consumption drop during unemployment (approximately 10–20% for households exhausting UI) implies $\gamma_c$ in the range 1–2. Substituting: $b^*/w \approx \gamma_c/(\gamma_c + 0.5) \approx 0.67$ for $\gamma_c = 1$ — a 67% replacement rate, close to the U.S. UI replacement rate for median earners under normal circumstances, and consistent with the interpretation that the U.S. UI system is roughly appropriately designed for average workers in average conditions, but may be suboptimal for workers with high risk aversion (young, low-wealth) or during severe recessions when the consumption drop is larger.

### UI as Automatic Stabilizer

Beyond the individual insurance motive, UI is one of the most powerful automatic stabilizers because UI recipients are overwhelmingly liquidity-constrained and have MPCs close to one — they spend their UI checks rapidly out of necessity. Gruber (1997) estimates that UI prevents approximately 30% of the consumption fall that would otherwise occur during unemployment spells, stabilizing aggregate consumption demand.

Chodorow-Reich et al. (2012) exploit cross-state variation in UI generosity during the Great Recession to estimate the employment multiplier of extended UI benefits. Their finding: states that received more UI extension weeks had significantly less employment decline, not more — contradicting the moral hazard view that UI raises unemployment by reducing search intensity. The aggregate demand effect of UI (sustaining consumption of the unemployed, who spend the transfers) dominated the moral hazard effect in the severe recession environment.

---

## 31.4 Active Labor Market Policies

**Passive** labor market policies (UI, disability insurance, early retirement) provide income replacement but do not actively facilitate reemployment or reduce structural unemployment. **Active labor market policies (ALMPs)** aim to improve the matching process, raise the employability of displaced workers, and reduce the structural natural rate directly.

### Job Search Assistance and Activation

Job search assistance — including counseling, interview preparation, referrals to vacancies, and mandatory reporting requirements that activate displaced workers — is the highest bang-for-buck ALMP in the empirical literature. Card, Kluve, and Weber (2018) conduct a comprehensive meta-analysis of approximately 200 randomized evaluations and find that job search assistance has the largest average impact on employment rates of any ALMP category, with effects appearing quickly (within six months) and persisting for two to three years. The mechanism is straightforward: search assistance reduces the mismatch friction in the Mortensen-Pissarides model [Ch. 13] directly, increasing the matching function's efficiency parameter $s$ without requiring extensive resource investment.

### Retraining and Skills Programs

Retraining displaced workers for occupations in growing sectors addresses structural unemployment directly. The evidence is more mixed than for job search assistance but increasingly positive for high-quality, longer-duration programs. Short-term training courses (a few weeks) typically show minimal effects on wages and employment because they cannot substantially alter workers' skill stocks. Longer-term vocational programs — particularly those with employer partnerships and on-the-job training components — show wage effects of 5–15% and employment probability gains of 5–10 percentage points at horizons of two to five years (OECD, 2022).

The key institutional design features for effective training programs: strong linkages to actual employer skill demands (programs designed with employer input are substantially more effective than programs designed without it); a combination of classroom and workplace training; and targeting workers whose previous occupational skills have been rendered obsolete by technology or trade shocks (who have the most to gain from retraining).

### Wage Subsidies

Employment subsidies reduce the cost of hiring targeted workers — typically long-term unemployed, older workers facing age discrimination, or workers with disabilities. The primary risks are **deadweight loss** (the employer would have hired the worker without the subsidy) and **displacement** (the subsidized hire replaces a non-subsidized worker who would otherwise have been hired). Well-evaluated subsidy programs show positive net employment effects when targeted at workers facing specific barriers to employment (stigma from long unemployment spells, employer discrimination) that create a gap between the worker's productivity and the employer's willingness to hire.

### Geographic Mobility and Place-Based Policies

Geographic mismatch — unemployed workers in depressed regions, vacancies in growing regions — is a significant component of structural unemployment. The theoretical solution is geographic mobility: subsidize workers to relocate to higher-opportunity areas. In practice, take-up of mobility assistance programs is low because workers face housing costs, family ties, and community attachments that create strong preference for remaining in place. European evidence on relocation subsidies (Germany, Sweden) suggests modest effects on mobility, substantially smaller than the implied welfare gains from relocation.

The growing alternative is **place-based policy**: rather than bringing workers to jobs, bring investment to workers. Enterprise zones, regional development grants, and targeted infrastructure investment in economically lagging areas all attempt to create employment opportunities in place. Evidence on effectiveness is mixed: some enterprise zone programs show local employment gains with minimal displacement of activity from nearby areas; others show high deadweight loss. The most successful place-based programs combine tax incentives with complementary investments in infrastructure, education, and public services that address multiple barriers to private investment simultaneously.

---

## 31.5 The Social and Psychological Costs of Unemployment

Aggregate unemployment statistics capture the labor market failure in macroeconomic terms — lost output, foregone tax revenue, higher transfer spending. They systematically understate the total welfare cost because they miss the large and persistent non-income effects of unemployment on individuals.

### Subjective Wellbeing

Clark and Oswald (1994), using British Household Panel Survey data, found that unemployment reduces self-reported life satisfaction by approximately 0.4 standard deviations — an effect substantially larger than the equivalent income loss alone would predict. Subsequent research using panel data in multiple countries has confirmed this finding and shown that unemployment has a "pure" psychological cost beyond the income channel: controlling for income, unemployed individuals have significantly lower life satisfaction than employed individuals, and the gap persists even after re-employment for workers who experienced long spells of unemployment (scarring of wellbeing, not just of wages).

The mechanisms: unemployment disrupts daily structure and social connection (Jahoda, 1982), damages occupational identity and self-efficacy, and exposes workers to social stigma and loss of status in ways that income replacement through UI cannot compensate for. The implication for cost-benefit analyses of stabilization policy is significant: the welfare cost of a one-percentage-point increase in unemployment is substantially greater than the income loss captured in output gap measures. Estimates that include wellbeing effects roughly double the implied welfare cost of recessions.

### Health Effects

Sullivan and von Wachter (2009) exploit mass layoffs in Pennsylvania manufacturing — which generate truly exogenous unemployment unrelated to the displaced workers' own characteristics — and track mortality rates for 20 years afterward. Their finding: workers displaced in mass layoffs experience a 50–100% increase in annual mortality probability in the year immediately following job loss, declining gradually but remaining elevated for 15–20 years. The cumulative mortality effect implies that job displacement reduces life expectancy by approximately 1–1.5 years. Mechanisms include stress-related cardiovascular disease, increases in smoking and alcohol consumption, reduced access to employer-provided health insurance, and depression.

### Family and Intergenerational Effects

Unemployment reduces marriage rates and increases divorce rates; it delays childbearing and, among households that remain intact, reduces children's educational attainment through parental stress, income instability, and residential instability. Oreopoulos, Page, and Stevens (2008) find that children of workers displaced in mass layoffs have lower earnings as adults — approximately 9% lower annual earnings — consistent with significant intergenerational transmission of the unemployment shock. This intergenerational channel substantially raises the long-run social cost of unemployment, extending the welfare loss of recessions beyond the generation directly affected.

### Crime and Community Disintegration

Raphael and Winter-Ebmer (2001) find that a one-percentage-point increase in the unemployment rate is associated with approximately a 2–4% increase in property crime rates, consistent with the economic theory that crime is partly an income-substitution response when legal employment income falls. Beyond individual crime, high-unemployment communities experience reduced social trust, higher rates of political disengagement, and deteriorating local institutions — effects that compound over years and are difficult to reverse even when employment recovers.

---

## Chapter Summary

- The taxonomy of unemployment — **frictional** (matching takes time), **structural** (skills/location mismatch), **cyclical** (demand deficiency) — implies fundamentally different policy responses; distinguishing them in real time is one of applied macroeconomics' hardest problems.

- **Hysteresis** $\dot{u}^* = \gamma(u - u^*)$ transmits cyclical unemployment into permanent natural rate increases through human capital atrophy, employer stigma, and labor force withdrawal. Blanchard-Summers (1986) documented severe hysteresis in European labor markets; its presence raises the welfare benefits of aggressive stabilization substantially.

- The **Baily-Chetty formula** $b^*/w = \gamma_c/(\gamma_c + \varepsilon_{1-e,b})$ characterizes optimal UI replacement rates. With $\varepsilon_{1-e,b} \approx 0.5$ and $\gamma_c \approx 1$, the optimum is approximately 67% — close to U.S. practice. In severe recessions, the automatic stabilizer effect of UI dominates the moral hazard cost.

- **Active labor market policies** vary in effectiveness: job search assistance has the highest and most robust returns; longer-duration retraining shows significant positive effects at 2–5 year horizons; wage subsidies succeed when targeted at workers facing specific hiring barriers; geographic mobility programs have low take-up and are increasingly complemented by place-based policy.

- The **social and psychological costs** of unemployment — subjective wellbeing losses (Clark-Oswald), mortality increases (Sullivan-von Wachter), intergenerational earnings losses (Oreopoulos et al.), and crime — substantially exceed the income losses captured in GDP-based welfare measures, providing a strong quantitative case for aggressive stabilization policy.

---

*Next: Chapter 32 — Open Economy Macroeconomics: Exchange Rates and Global Imbalances*
