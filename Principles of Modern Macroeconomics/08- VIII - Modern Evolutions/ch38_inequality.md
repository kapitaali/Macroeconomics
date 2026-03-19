# Chapter 38 — Inequality and Macroeconomics: Distributional Effects of Policy

---

The standard macroeconomic models developed in this book — from the Solow growth model to the New Keynesian three-equation system — typically feature a representative household. By assumption, these models produce no distributional outcomes: there is only one household, so there is nothing to distribute. This abstraction is analytically convenient but empirically and normatively problematic. In the real world, income and wealth are distributed very unevenly; the macroeconomic consequences of any policy depend on who receives the gains and who bears the costs; and the distributional consequences of macroeconomic events are themselves the subject of intense political conflict.

This chapter develops the macroeconomics of inequality in three steps. First, measurement: the major statistical tools for characterizing inequality and their properties. Second, the long-run facts of inequality and what drives them. Third, the distributional consequences of the major macroeconomic policy instruments — monetary policy, fiscal policy, and trade — showing that even policies evaluated primarily on aggregate efficiency grounds have substantial distributional implications that feed back into political economy and hence into future policy choices.

---

## 38.1 Measuring Inequality: Tools and Their Properties

**Definition (Gini Coefficient).** The **Gini coefficient** $G$ is a scalar summary of the inequality in a distribution of income (or wealth), ranging from 0 (perfect equality: everyone has the same income) to 1 (perfect inequality: all income accrues to one person). For a distribution with $n$ individuals with incomes $y_1 \leq y_2 \leq \cdots \leq y_n$ and mean $\mu$:

$$G = \frac{2}{n^2\mu}\sum_{i=1}^n i\,y_i - \frac{n+1}{n} = 1 - \frac{2}{n}\sum_{i=1}^n\!\left(1 - \frac{\sum_{j=1}^i y_j}{\mu n}\right).$$

Geometrically, $G$ is twice the area between the **Lorenz curve** (the share of cumulative income held by the bottom $x\%$ of the population) and the 45-degree line of perfect equality.

The Gini coefficient has useful properties but also important limitations. It summarizes the entire distribution in one number, which is convenient for cross-country comparison but loses information about the distribution's shape. In particular, the Gini is relatively insensitive to changes at the very top of the distribution — two distributions can have the same Gini but very different top income shares. This motivates the use of **top income shares** (the share of national income accruing to the top 1%, 0.1%, or 0.01%) as complementary measures, which are more sensitive to changes in extreme inequality.

**Definition (Theil Index).** The **Theil index** is a decomposable inequality measure derived from information theory:

$$T = \frac{1}{n}\sum_{i=1}^n \frac{y_i}{\mu}\ln\!\left(\frac{y_i}{\mu}\right).$$

Unlike the Gini, the Theil index can be decomposed into within-group and between-group components: $T = T^{within} + T^{between}$. This decomposability makes it valuable for analyzing the sources of inequality — how much of total inequality is attributable to differences between demographic groups versus differences within groups.

**Definition (Lorenz Curve and Stochastic Dominance).** The **Lorenz curve** $L(p)$ gives the share of total income received by the bottom $p$ fraction of the population. If distribution $A$'s Lorenz curve lies entirely above distribution $B$'s for all $p \in [0,1]$, then $A$ **Lorenz dominates** $B$ — $A$ is unambiguously more equal than $B$ according to all inequality measures that are symmetric, monotone, and satisfy the Pigou–Dalton transfer principle (a transfer from rich to poor reduces inequality).

---

## 38.2 The Long-Run Facts of Inequality

### The Great Compression and the Great Divergence

Income inequality in the United States and other advanced economies followed a striking U-shaped pattern over the twentieth century. The **Great Compression** (Goldin and Margo, 1992) refers to the sharp decline in income inequality during the 1940s: the ratio of the 90th to the 10th percentile of wages fell dramatically as union coverage expanded, minimum wages were introduced, and wartime wage controls compressed the wage structure. Inequality continued to decline modestly through the 1960s.

The **Great Divergence** (Krugman, 2007) refers to the sustained increase in inequality from approximately 1978 onward. The top 1% income share in the United States, which had fallen from approximately 20% in the late 1920s to approximately 8% by 1978, rose back to approximately 20% by 2012 (Piketty and Saez, 2003). The college premium (wages of college-educated workers relative to non-college) rose substantially. Real wages for workers in the bottom quartile of the wage distribution stagnated or declined in real terms over 1979–2019.

### Piketty's $r > g$ Framework

Piketty (2014) proposed a simple theoretical framework for understanding long-run wealth inequality. In the long run, the share of income accruing to capital in national income is:

$$\alpha = r\cdot\beta,$$

where $r$ is the return on capital and $\beta = K/Y$ is the capital-income ratio (wealth-to-income ratio). From the Solow model, $\beta^* = s/(n+g)$: lower growth rates (lower $g$) generate higher capital-income ratios. Piketty's key argument: if the return on capital $r$ exceeds the growth rate $g$ — the **$r > g$ condition** — wealth tends to concentrate over time. Wealthy individuals earn a return $r$ on their capital, and if they save a fraction of this return, their wealth grows at rate $r > g$, faster than the overall economy. Inequality of wealth rises without bound unless offset by taxes on capital or inheritance.

The $r > g$ condition has been the subject of significant criticism. First, in the Ramsey model, $r = \rho + \sigma g$ — when $r > g$, impatient agents with $\rho > 0$ still consume enough to prevent the wealth distribution from diverging. Second, the relevant comparison is not $r$ vs. $g$ but the after-tax return on wealth $r(1-\tau_K)$ vs. $g$: progressive capital and inheritance taxes can stabilize wealth inequality even when $r > g$ pre-tax. Third, the empirical evidence on whether wealth inequality is actually rising toward Piketty's predicted concentration is mixed: while top income shares have risen, the share of the top 1% in *wealth* has risen more modestly (Saez and Zucman, 2016 vs. Kopczuk, 2015 — with large measurement uncertainty due to tax avoidance and offshore wealth).

---

## 38.3 Distributional Consequences of Monetary Policy

Standard analyses of monetary policy focus on aggregate effects: how does a rate change affect output and inflation? But monetary policy also redistributes income and wealth across groups. Understanding these distributional effects matters both intrinsically (for equity) and instrumentally (for understanding the political constraints on central bank independence).

### The Auclert Redistribution Framework

Auclert (2019) develops a systematic framework for analyzing monetary policy redistribution. He identifies three channels:

**The unhedged interest rate exposure (URE) channel**: households with more financial liabilities than assets (**net debtors**) gain when interest rates fall, because their interest payments fall while their asset returns remain unchanged. Net creditors (more assets than liabilities) lose. The distributional impact equals:

$$\mathrm{dC}_h = MPC_h\cdot\mathrm{URE}_h\cdot(-\mathrm{d}r),$$

where $\mathrm{URE}_h = \text{financial assets}_h - \text{financial liabilities}_h$ and $MPC_h$ is household $h$'s marginal propensity to consume. When a rate cut redistributes from creditors (low MPC, wealthy) to debtors (high MPC, poorer), the aggregate consumption effect is amplified by the redistribution.

**The Fisher channel**: unexpected inflation reduces the real value of nominal debt (benefiting debtors) and the real value of nominal assets (hurting creditors). With household debt at approximately 80% of GDP in the United States, a 1% inflation surprise redistributes approximately 0.8% of GDP from creditors to debtors.

**The income heterogeneity channel**: changes in wages and employment caused by monetary policy affect different households differently. Labor income constitutes the majority of income for poorer households; capital income is more important for wealthy households. An expansionary monetary shock that raises employment disproportionately benefits lower-income households through the labor market.

### Quantitative Easing and Wealth Inequality

QE operates primarily through the portfolio balance channel — raising equity and bond prices. Since wealth is highly concentrated (the top 1% hold approximately 40% of financial assets in the United States), QE's asset price effects disproportionately benefit wealthy households. McKinsey Global Institute estimated that U.S. QE increased the financial wealth of the top quintile by approximately $4.6 trillion more than it increased the wealth of the bottom four quintiles combined (2013).

However, the distributional assessment of QE must weigh its asset-price effects against its employment effects. If QE prevented a depression that would have devastated low-income workers, the distributional calculus may favor the poor even with unequal asset price effects. This counterfactual question — what would have happened without QE — is central to a complete distributional assessment and is inherently difficult to answer.

---

## 38.4 Fiscal Policy and Redistribution

Fiscal policy is the primary instrument of redistribution in advanced economies. The combination of progressive income taxes, means-tested transfers, and public provision of education and healthcare generates a systematic redistribution from higher to lower income households.

**The size of redistribution**: in OECD countries, the ratio of the market income Gini to the disposable income Gini (after taxes and transfers) is approximately 0.45 in the United States (relatively less redistribution) and 0.25 in Nordic countries (more redistribution). Taxes and transfers reduce the Gini by approximately 30–40% in most advanced economies.

**The incidence of taxes**: the distributional impact of a tax depends on who ultimately bears its burden — not necessarily who legally pays it. Payroll taxes (legally split between employers and employees) are likely borne primarily by workers in the long run, since labor supply is inelastic and the incidence falls on the inelastic side. Value-added taxes (sales taxes) are borne by consumers, making them regressive if the consumption share of income declines with income. Corporate income taxes are contested: theory suggests they are borne partly by capital (reducing the after-tax return) and partly by labor (through lower investment and wages), with the shares depending on the openness of the economy.

**Fiscal consolidation and distributional consequences**: austerity programs typically combine spending cuts and tax increases. The distributional impact depends critically on which spending is cut and which taxes are raised. Cuts to social programs (healthcare, housing assistance, food stamps) affect primarily lower-income households; cuts to public sector employment affect a broader range but with a concentration in middle-income households. The 2010–15 European austerity programs were associated with sharp increases in income inequality in the periphery countries (Greece, Spain, Portugal, Ireland) even as aggregate inequality in the EU remained roughly stable.

---

## 38.5 Trade, Technology, and the Political Economy of Inequality

The interaction between trade, technology, and inequality has generated intense political controversy in the 2010s–2020s, as the populations of advanced economies that benefited least from globalization turned to populist political movements. Understanding the macroeconomics underlying this political economy is essential.

The Stolper–Samuelson theorem (developed in Chapter 21) predicts that trade liberalization in a capital-abundant country reduces the real wage of the abundant factor (capital) is benefited and the scarce factor (low-skilled labor) is harmed. The China shock — China's WTO accession in 2001 and the resulting surge in Chinese manufacturing exports — provides a natural experiment. Autor, Dorn, and Hanson (2013) use cross-regional variation in exposure to Chinese import competition to estimate causal effects: regions more exposed to Chinese imports experienced persistently higher unemployment, lower wages, and greater social distress, with the negative effects lasting more than a decade. The estimated wage reduction for affected workers is approximately 3–7% of earnings, with employment losses of approximately 2–2.4 million U.S. manufacturing jobs attributable directly to Chinese competition.

The key political economy insight: the gains from trade liberalization are diffuse (slightly lower prices for all consumers) while the losses are concentrated (large wage and employment losses for specific workers in specific places). The politically visible concentrated losses generate more intense opposition than the diffuse gains generate support — a fundamental asymmetry that makes trade liberalization politically fragile even when it is economically beneficial in aggregate.

---

*Next: Chapter 39 — The Future of Macroeconomics*
