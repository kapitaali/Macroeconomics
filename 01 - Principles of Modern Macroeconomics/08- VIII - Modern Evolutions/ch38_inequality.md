# Chapter 38 — Inequality and Macroeconomics: Distributional Effects of Policy

> *"The issue is never whether there will be trade-offs. The question is always: who bears them?"*
> — Arthur Okun, *Equality and Efficiency: The Big Tradeoff*, 1975

---

The standard macroeconomic models developed in this book typically feature a representative household. By construction, these models produce no distributional outcomes: there is only one household, so there is nothing to distribute. This abstraction is analytically convenient but empirically and normatively costly. Income and wealth are distributed very unevenly in every country; the macroeconomic consequences of any policy depend on who receives the gains and who bears the costs; and the distributional consequences of macroeconomic events generate the political economy that shapes future policy. This chapter develops the macroeconomics of inequality — measurement tools, the long-run facts, the distributional consequences of monetary and fiscal policy, and the connections between trade, technology, and political backlash.

---

## 38.1 Measuring Inequality: Tools and Their Properties

### The Gini Coefficient

**Definition (Gini Coefficient).** The **Gini coefficient** $G$ measures inequality in a distribution of income or wealth, ranging from 0 (perfect equality) to 1 (perfect inequality). For $n$ individuals with incomes $y_1 \leq y_2 \leq \cdots \leq y_n$ and mean $\mu$:

$$G = \frac{2}{n^2\mu}\sum_{i=1}^n i\,y_i - \frac{n+1}{n}.$$

Geometrically, $G$ is twice the area between the **Lorenz curve** — the share of cumulative income held by the bottom $p$ fraction of the population — and the 45-degree line of perfect equality: $G = 1 - 2\int_0^1 L(p)\,dp$.

The Gini has several useful properties: it satisfies the Pigou-Dalton transfer principle (a transfer from rich to poor reduces $G$), it is scale-invariant, and it is symmetric. Its main limitation: two distributions can have identical Gini coefficients but very different shapes. The Gini is relatively insensitive to changes at the very top of the distribution, which is where inequality has changed most dramatically since 1980. This motivates using **top income shares** (the share of national income accruing to the top 1%, 0.1%, or 0.01%) as complementary measures, more sensitive to the extreme right tail.

### The Theil Index and Decomposability

**Definition (Theil Index).** The **Theil index** is:

$$T = \frac{1}{n}\sum_{i=1}^n \frac{y_i}{\mu}\ln\!\left(\frac{y_i}{\mu}\right).$$

The Theil index decomposes additively: $T = T^{within} + T^{between}$, where $T^{within}$ is inequality within demographic groups (by education, race, region) and $T^{between}$ is inequality between groups. This decomposability makes it valuable for analyzing whether rising overall inequality reflects widening gaps between groups or widening distributions within previously homogeneous groups.

### Stochastic Dominance

When Lorenz curve $A$ lies entirely above Lorenz curve $B$ for all $p \in [0,1]$, distribution $A$ **Lorenz dominates** $B$ — $A$ is unambiguously more equal than $B$ according to every inequality measure consistent with the Pigou-Dalton principle. Lorenz dominance provides a partial ordering of distributions that is more robust than any single summary statistic but is not a complete ordering (many distribution pairs have crossing Lorenz curves, where neither dominates the other).

---

## 38.2 The Long-Run Facts of Inequality

### The Great Compression and the Great Divergence

Income inequality in the United States followed a U-shaped pattern over the twentieth century. The **Great Compression** (Goldin and Margo, 1992): inequality fell sharply in the 1940s as union coverage expanded from approximately 15% to 35% of the workforce, the minimum wage was introduced and then raised substantially, and wartime wage controls compressed the wage structure. The ratio of 90th to 10th percentile wages fell by approximately 25% between 1940 and 1950 — a historically unprecedented compression.

The **Great Divergence** (Krugman, 2007): from approximately 1978, inequality rose dramatically. The top 1% income share in the United States fell from approximately 20% in 1928 to approximately 8% in 1978 and then rose back to approximately 20% by 2012 (Piketty and Saez, 2003). Real wages for workers in the bottom quartile stagnated or declined in real terms over 1979–2019. The college premium (wages of college graduates relative to high school graduates) rose by approximately 30% over the same period.

The major proximate causes of the Great Divergence: skill-biased technological change (computers and automation raised the relative demand for skills); globalization and the China shock (import competition reduced wages for low-skill manufacturing workers — Autor, Dorn, and Hanson, 2013); declining unionization (union wage premium is approximately 15–20%, and union coverage fell from 35% to 10%); and the erosion of the real minimum wage (the federal minimum wage, adjusted for inflation, fell by approximately 25% between 1968 and 2000). Each factor is quantitatively significant; debates about relative magnitudes continue.

### Piketty's $r > g$ Framework

Piketty (2014) proposed a simple framework for understanding long-run wealth inequality. The capital income share:

$$\alpha = r\cdot\beta,$$

where $r$ is the return on capital and $\beta = K/Y$ is the wealth-to-income ratio. From the Solow model, $\beta^* = s/(n+g)$: lower growth rates raise the capital-income ratio. If $r > g$ — the **Piketty condition** — wealth tends to concentrate: wealthy individuals earning return $r$ on capital with some fraction saved will see their wealth grow at rate $r > g$, outpacing the overall economy.

The Piketty condition has been contested on several grounds. First, in the Ramsey model $r = \rho + \sigma g$: impatient agents with $\rho > 0$ still consume enough to prevent wealth from diverging as dramatically as Piketty implies. Second, the relevant comparison is the after-tax return $(1-\tau_K)r$ vs. $g$: progressive capital and inheritance taxes can stabilize wealth inequality even when the pre-tax $r > g$. Third, empirical evidence on whether wealth inequality has actually followed Piketty's predicted trajectory is contested: Saez and Zucman (2016) find dramatic concentration of the top 0.1% wealth share; Kopczuk and Schoar document more modest trends, with significant measurement uncertainty from offshore wealth concealment and trusts.

---

## 38.3 Distributional Consequences of Monetary Policy

Standard analyses of monetary policy focus on aggregate effects. But monetary policy also redistributes income and wealth, with consequences for both equity and the political sustainability of central bank independence.

### The Auclert Redistribution Framework

Auclert (2019) identifies three channels through which monetary policy redistributes:

**The Unhedged Interest Rate Exposure (URE) channel**: households whose financial liabilities exceed their financial assets (net debtors) gain when rates fall, because their debt service costs fall while their interest income falls less. Net creditors (more assets than liabilities) lose. The aggregate consumption effect:

$$\mathrm{d}C = \int MPC_h \cdot URE_h \cdot (-\mathrm{d}r)\,\mathrm{d}h,$$

where $URE_h = \text{financial assets}_h - \text{financial liabilities}_h$. Because net debtors tend to have higher MPCs than net creditors (debtors are typically less wealthy, with fewer precautionary buffers), the URE channel amplifies the aggregate consumption effect of rate cuts.

**The Fisher channel**: unexpected inflation reduces the real value of all nominal contracts. Debtors holding fixed-rate mortgages or nominal bonds gain; creditors lose. With household debt at approximately 80% of GDP in the United States, a 1-percentage-point inflation surprise redistributes approximately 0.8% of GDP from nominal creditors to debtors.

**The income heterogeneity channel**: expansionary monetary policy raises employment and wages. Because labor income constitutes the majority of income for lower-income households and capital income dominates for wealthy households, the employment channel of monetary policy is distributionally progressive — its gains accrue disproportionately to lower-income households. Kaplan, Moll, and Violante (2018) find this is quantitatively the most important channel in HANK models.

### Quantitative Easing and the Wealth Distribution

QE operates primarily through asset price appreciation: lower discount rates raise the present value of future cash flows, boosting equity and real estate prices. Since wealth is highly concentrated (the top 1% own approximately 40% of financial assets in the United States), QE's asset price effects disproportionately benefit wealthy households — a regressive distributional consequence. The McKinsey Global Institute (2013) estimated that U.S. QE increased financial wealth by approximately $4.6 trillion for the top quintile versus modest increases for the lower quintiles.

However, the distributional assessment of QE must weigh asset price effects against employment effects. If QE prevented a depression that would have devastated low-income workers, the employment channel may dominate the asset price channel in welfare terms. This counterfactual is inherently difficult to establish but is central to a complete evaluation. The HANK framework [Ch. 39] suggests that the employment-income channel is large relative to the intertemporal substitution channel that is QE's primary mechanism in RANK models — potentially partially offsetting the asset-price inequality effect.

---

## 38.4 Fiscal Policy and Redistribution

### The Size and Mechanics of Redistribution

Fiscal policy is the primary instrument of redistribution in advanced economies. The combination of progressive income taxes, means-tested transfers, and universal provision of education and healthcare generates systematic redistribution from higher to lower income households. In OECD countries, the market income Gini (pre-tax, pre-transfer) is approximately 0.45–0.50; the disposable income Gini (post-tax, post-transfer) is approximately 0.25–0.35 in Nordic countries and approximately 0.38 in the United States. Taxes and transfers reduce the Gini by approximately 30–40% — a large effect.

### Tax Incidence and Distributional Consequences

The distributional impact of a tax depends on who ultimately bears its burden, not who legally pays it. The incidence depends on the relative elasticities of supply and demand:

**Payroll taxes**: legally split between employers and employees. In the long run, labor supply is relatively inelastic (workers must work) and the labor demand curve is relatively elastic. The incidence falls primarily on workers through lower take-home wages, making payroll taxes more regressive than their legal structure suggests.

**Corporate income taxes**: contested incidence. In a closed economy, the incidence falls primarily on capital (shareholders and capital owners bear the tax through lower after-tax returns). In an open economy with mobile international capital, the burden may shift to labor through reduced investment and wages. Empirical estimates of the capital-labor split range from 20/80 to 70/30, reflecting the sensitivity to assumptions about capital mobility.

**Value-added taxes**: borne by consumers in proportion to their consumption. Since the consumption share of income declines with income (wealthy households save more), VAT is regressive: lower-income households pay a higher fraction of their income in VAT than wealthy households. However, VATs combined with targeted transfers to low-income households can be made distributionally neutral.

### Fiscal Consolidation and Distribution

Austerity programs — cutting spending and raising taxes — have distributional consequences that depend critically on composition. Cuts to social programs (UI, food assistance, housing vouchers, healthcare) affect primarily lower-income households; cuts to public employment affect a broader range. The 2010–15 European austerity programs — concentrated in the periphery (Greece, Spain, Portugal, Ireland) — were associated with sharp increases in measured disposable income inequality even as market income inequality rose only modestly, as the transfer programs that had previously compressed disposable income were cut back.

---

## 38.5 Trade, Technology, and the Political Economy of Inequality

### The China Shock

China's accession to the WTO in 2001 and the subsequent surge in Chinese manufacturing exports created a natural experiment. Autor, Dorn, and Hanson (2013) exploit cross-regional variation in exposure to Chinese import competition — regions with industry mixes more exposed to Chinese imports experienced sharper declines in manufacturing employment — to estimate causal effects: a $1,000 per worker increase in Chinese import exposure reduced wages by approximately $\$312$ and employment by approximately 0.6 percentage points annually. The estimated total: approximately 2–2.4 million U.S. manufacturing job losses attributable to Chinese competition, with the effects concentrated geographically in the Midwest and South.

The key political economy asymmetry: the gains from trade liberalization are diffuse (slightly lower prices for all consumers — perhaps $\$250$/year per household from Chinese imports) while the losses are concentrated (severe wage and employment losses for specific workers in specific places, with effects lasting more than a decade). The politically visible concentrated losses generate more intense opposition than the diffuse gains generate support — a fundamental asymmetry that makes trade liberalization politically fragile even when it generates positive net aggregate welfare.

### Technology and Distributional Consequences

As Chapter 36 develops, automation of routine tasks has simultaneously increased labor demand for high-skill non-routine cognitive workers and decreased demand for middle-skill routine workers — the polarization hypothesis. Acemoglu and Restrepo's task-based framework implies that the distributional consequences of technology depend critically on whether reinstatement (new tasks for displaced workers) keeps pace with displacement. The evidence from the 1980s–2010s suggests that reinstatement has been insufficient to prevent polarization.

The distributional consequences of AI — which threatens non-routine cognitive tasks — are potentially larger: if AI displaces knowledge workers who previously benefited from technological change, the political economy of technology-driven inequality becomes more widespread and more politically volatile.

### Inequality and Political Macroeconomics

The macroeconomic consequences of rising inequality feed back into policy choices through political economy. High inequality reduces social trust, increases political polarization, and generates pressure for protectionist and redistributive policies that may reduce economic efficiency while addressing distributional concerns. Mian, Straub, and Sufi (2021) argue that rising inequality has been a primary driver of the secular decline in interest rates: wealthy households save more (lower $r$) while workers borrow to sustain consumption (higher leverage); the resulting imbalance requires ever-lower interest rates to maintain demand. This perspective integrates inequality analysis directly into the macroeconomic modeling of secular stagnation.

---

## Chapter Summary

- **Gini coefficient** $G = 2\int_0^1(p - L(p))dp$ summarizes inequality; Theil index $T$ decomposes within/between groups. Top income shares are more sensitive to extreme inequality. Lorenz dominance provides a robust partial ordering.

- The **Great Compression** (1940s: union expansion, minimum wage, wartime controls) reduced U.S. inequality dramatically; the **Great Divergence** (1978–present: skill-biased technology, China shock, union decline, minimum wage erosion) reversed it. The top 1% income share returned to its 1928 level by 2012. **Piketty's $r > g$ framework** explains wealth concentration dynamics; contested on Ramsey model grounds and empirically.

- **Monetary policy distributes** through three Auclert channels: URE (net debtors gain from rate cuts), Fisher (debtors gain from inflation surprises), and income heterogeneity (expansions benefit labor income disproportionately). QE's asset price channel benefits the wealthy; the employment channel benefits the poor; the net effect is contested.

- **Fiscal redistribution** reduces the Gini by 30–40% in advanced economies. Tax incidence diverges from legal structure: payroll taxes fall primarily on workers; VAT is regressive; corporate tax incidence depends on capital mobility. Austerity disproportionately cuts benefits to lower-income households.

- The **China shock** (Autor-Dorn-Hanson) caused approximately 2 million U.S. manufacturing job losses with geographically concentrated long-lasting effects. The political economy asymmetry — diffuse gains, concentrated losses — makes trade liberalization politically fragile. Rising inequality may feed back into lower $r^n$ through the Mian-Straub-Sufi saving mechanism.

---

*Next: Chapter 39 — The Future of Macroeconomics*
