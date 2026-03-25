# Chapter 25 — Households and Demographics: Saving, Aging, and Consumption

> *"Demographics are destiny."*
> — Auguste Comte (attributed)

---

The household sector is the ultimate source of both the labor and the saving that drives capital accumulation, and the demographic structure of the household sector — the age distribution of the population — determines many of macroeconomics' most important long-run trends. An aging population saves less in aggregate, demands more from pension and healthcare systems, and exerts pressure on public finances that will shape fiscal policy for decades. The massive demographic transition currently underway in virtually every advanced economy — the retirement of the post-war baby boom generation — is arguably the most predictable and consequential macroeconomic event of the twenty-first century.

This chapter examines the economics of household saving and demographics in three parts: the overlapping-generations (OLG) model, which provides the rigorous framework for analyzing the macroeconomic implications of demographic change; the life-cycle saving profile and its aggregate implications; and the heterogeneous-agent framework (HANK models) that has transformed our understanding of how monetary and fiscal policy operates in a world with large wealth inequality.

---

## 25.1 The Diamond Overlapping-Generations Model

The Ramsey–Cass–Koopmans model of Chapter 5 assumes an infinitely-lived representative household. This is methodologically convenient but economically problematic: it abstracts from the central fact that households have finite lives, that different generations coexist simultaneously, and that intergenerational transfers — social security, bequests, public debt — involve distributional choices among generations with fundamentally different interests.

The **Diamond (1965) overlapping-generations (OLG) model** addresses these limitations by explicitly modeling the coexistence of multiple generations. In the simplest two-period version:

- At each date $t$, a new "young" generation of mass 1 is born.
- Young households work in period $t$, earning wage $w_t$, and choose how much to save.
- Old households (born at $t-1$) retire in period $t$, consuming their savings plus capital income.
- Each generation lives exactly two periods; there is no bequest motive.

**Definition (Overlapping Generations Equilibrium).** An **OLG equilibrium** is a sequence of wages, interest rates, and capital stocks $\{w_t, r_t, k_t\}_{t=0}^\infty$ such that:
1. Young households maximize two-period utility taking prices as given.
2. Firms maximize profits.
3. Capital market clears: $k_{t+1} = s_t$ (the young's saving becomes the next period's capital).

### Household Optimization

Young household utility over consumption when young $c_t^y$ and old $c_{t+1}^o$:

$$\max_{c_t^y, c_{t+1}^o}\; u(c_t^y) + \beta u(c_{t+1}^o) \quad \text{s.t.} \quad c_t^y + s_t = w_t, \quad c_{t+1}^o = (1+r_{t+1})s_t.$$

First-order condition: $u'(c_t^y) = \beta(1+r_{t+1})u'(c_{t+1}^o)$ — the Euler equation. The saving function $s_t = s(w_t, r_{t+1})$ is positively related to $w_t$ (higher wage income → more lifetime wealth → more saving) and ambiguously related to $r_{t+1}$ (substitution effect: higher return makes future consumption cheaper, encouraging saving; income effect: need less saving to achieve given future consumption, discouraging saving). The standard assumption is that the substitution effect dominates, so $\partial s/\partial r > 0$.

### Capital Market Equilibrium and Dynamics

Capital market equilibrium: $k_{t+1} = s(w_t, r_{t+1})$. Using the competitive factor pricing relations $w_t = f(k_t) - k_t f'(k_t)$ and $r_t = f'(k_t) - \delta$ (where $f(k)$ is output per worker), this becomes:

$$k_{t+1} = s\bigl(f(k_t) - k_t f'(k_t),\; f'(k_{t+1}) - \delta\bigr).$$

This is a first-order implicit difference equation in $k_t$. With Cobb–Douglas $f(k) = k^\alpha$: $w_t = (1-\alpha)k_t^\alpha$ and $r_t = \alpha k_t^{\alpha-1} - \delta$. The saving function under log utility is $s_t = \beta w_t/(1+\beta)$, giving the explicit dynamics:

$$k_{t+1} = \frac{\beta}{1+\beta}(1-\alpha)k_t^\alpha.$$

This is a well-behaved map with a unique stable steady state $k^*$ satisfying $k^* = \frac{\beta(1-\alpha)}{1+\beta}(k^*)^\alpha$.

### Dynamic Inefficiency: A Key OLG Implication

The most important macroeconomic implication of the OLG structure — absent from the Ramsey model — is that the competitive equilibrium can be **dynamically inefficient**: the steady-state capital stock may exceed the Golden Rule level, $k^* > k^{GR}$.

**Definition (Dynamic Inefficiency in OLG).** An OLG economy is **dynamically inefficient** if the steady-state real interest rate is below the growth rate: $r^* < n + g$. In this case, the economy has overaccumulated capital: each generation is saving too much for retirement, driving the capital stock above the level that maximizes steady-state consumption. A policy that reduces saving — such as a pay-as-you-go pension system that transfers from young to old — can make every generation simultaneously better off by reducing the capital stock toward the Golden Rule.

Dynamic inefficiency arises because OLG households do not internalize the effect of their saving on future wages and interest rates, and because there is no mechanism (such as the bequest motive) to correct this externality. Whether actual economies are dynamically inefficient is an empirical question: Abel et al. (1989) found that aggregate investment in the United States is below aggregate capital income ($I < rK$) for all years in their sample, consistent with dynamic efficiency. However, with heterogeneous assets and incomplete markets, the test is more nuanced.

### Social Security and Intergenerational Transfers

A **pay-as-you-go (PAYG) social security system** takes contributions from the current young and distributes them immediately to the current old (as opposed to a fully-funded system that accumulates the contributions and returns them to the same generation in retirement). In an OLG model, PAYG social security:

1. Reduces the young generation's saving (by reducing the need for private retirement saving).
2. Transfers resources from young to old (directly, through the contribution-benefit mechanism).
3. Reduces the capital stock (since the young save less).

If the economy is dynamically inefficient ($r < n+g$), PAYG social security can improve welfare: by reducing overaccumulated capital, it moves the economy toward the Golden Rule and all generations can benefit simultaneously. If the economy is dynamically efficient ($r > n+g$), PAYG social security reduces the capital stock from an already-insufficient level, and the welfare gains of current retirees come at the expense of future generations.

---

## 25.2 Population Aging and the Macroeconomy

The world is aging rapidly. The old-age dependency ratio (population over 65 divided by population 15–64) in the OECD is projected to rise from approximately 0.28 (2020) to 0.53 (2060). Japan, South Korea, and Southern European countries face even more dramatic aging. This demographic transition has profound macroeconomic implications.

### Life-Cycle Saving and the Aggregate Saving Rate

Modigliani's life-cycle hypothesis (Chapter 11) implies that the aggregate saving rate depends on the age distribution of the population. Workers in their prime earning years (40–55) save more than the young (who borrow for education and housing) or the old (who dissave to fund retirement):

$$S_{agg} = \int_0^\infty s(age)\cdot f(age)\,\mathrm{d}(age),$$

where $s(age)$ is the age-specific saving rate and $f(age)$ is the population density. As the baby boom generation retires, $f(age)$ shifts toward the high-dissaving segment, reducing $S_{agg}$.

The macroeconomic implications: (i) lower household saving → lower domestic investment → potentially lower long-run growth (Solow model); (ii) lower private saving financed by public pension liabilities → fiscal pressure as dependency ratios rise; (iii) current account deterioration as domestic saving falls relative to investment — or capital outflows to younger economies if investment also falls.

Whether aging reduces interest rates is contested. The "secular stagnation" hypothesis (Summers, 2014; Rachel and Smith, 2017) argues that demographic aging, combined with rising inequality and declining investment demand, has reduced the natural rate of interest $r^n$ persistently below zero, trapping advanced economies at the ELB. The evidence: the real interest rate on government bonds has been declining in advanced economies since the 1980s, from roughly 3–4% to near zero or below by the 2010s.

### The Fiscal Consequences of Aging

Higher old-age dependency ratios directly increase mandatory government spending on pensions and healthcare. Projections by the IMF and OECD suggest that pension and healthcare spending in advanced economies will rise by 4–8 percentage points of GDP over 2020–2060 due to aging alone, absent policy reforms. The resulting fiscal pressures create a dilemma: raise contribution rates (taxing the shrinking working-age population more), cut benefits (imposing losses on current and future retirees), raise the retirement age (extending working lives), or accept higher debt accumulation (shifting the burden to future generations through the debt dynamics equation).

The intergenerational equity dimension of these choices is fundamental and involves normative judgments about the distribution of welfare across generations — one of the "ethical considerations" developed in Appendix I.

---

## 25.3 Heterogeneous Households and the HANK Framework

The representative-agent models of earlier chapters assume all households are identical — same income, wealth, and preferences. This is obviously false. In the United States, the top 1% of households own approximately 30% of total wealth, and the bottom 50% hold near-zero liquid assets. This wealth heterogeneity is not merely an equity concern: it profoundly affects the aggregate consumption response to monetary and fiscal policy, because households at different wealth levels have very different marginal propensities to consume (MPCs).

**Definition (Hand-to-Mouth Households).** A household is **hand-to-mouth** if it consumes essentially all of its current income — it holds near-zero liquid wealth and its MPC is approximately one. Hand-to-mouth households may not literally have zero wealth: they may hold illiquid wealth (housing, pension funds) that is difficult to access on short notice. The **wealthy hand-to-mouth** (Kaplan and Violante, 2014) are households with substantial illiquid wealth but minimal liquid savings, who behave hand-to-mouth with respect to income changes because liquidating illiquid wealth is costly. This category includes many homeowners and pension savers.

In the United States, approximately 30–40% of households are classified as poor or wealthy hand-to-mouth based on balance sheet data (Kaplan, Violante, and Weidner, 2014). Their aggregate MPC is approximately 0.5–0.8 in response to temporary income transfers, compared to near-zero for wealthy households who smooth consumption over the full lifecycle.

### The HANK Model

The **Heterogeneous-Agent New Keynesian (HANK) model** (Kaplan, Moll, and Violante, 2018) combines the Aiyagari–Bewley incomplete-markets household problem with New Keynesian nominal rigidities to generate macroeconomic dynamics that account for wealth heterogeneity. In this framework, each household $h$ faces idiosyncratic labor income risk $y_t^h$ and solves:

$$\max_{\{c_t^h, a_{t+1}^h\}}\; \mathbb{E}_0\sum_{t=0}^\infty\beta^t u(c_t^h) \quad \text{s.t.} \quad c_t^h + a_{t+1}^h = (1+r_t)a_t^h + w_t n_t^h - T_t^h,$$

subject to a borrowing constraint $a_{t+1}^h \geq \underline{a}$. The idiosyncratic income shocks $y_t^h$ are uninsurable (no complete asset markets): households can only self-insure by accumulating liquid assets. The stationary distribution of household wealth $\mathcal{F}(a)$ is the key state variable, tracking the cross-sectional distribution of assets.

The aggregate consumption response to monetary policy in HANK:

$$\hat{C}_t = \int \hat{c}(h,\, \hat{a}_t^h,\, \hat{y}_t^h)\,\mathcal{F}(\mathrm{d}h).$$

The critical insight of Kaplan, Moll, and Violante (2018): in HANK, the **direct effect** of a monetary policy rate cut on consumption (the intertemporal substitution effect, operating through the Euler equation) is quantitatively small, because the intertemporal elasticity of substitution is low and most households are not on the margin of the Euler equation. The **indirect effect** — the income channel, operating through higher employment and wages — is much larger and dominates. This reverses the RANK conclusion: monetary policy works primarily through the labor market, not through interest rate effects on saving.

Fagereng, Holm, Moll, and Natvik (2021) provide quasi-experimental evidence: Norwegian lottery winnings of random size identify the MPC out of unexpected income across the wealth distribution. Their estimated average MPC over four years is approximately 0.5, declining from near 1 for the lowest-wealth households to near 0 for the highest-wealth group — consistent with the HANK model's heterogeneous MPC structure.

---

*Next: Chapter 26 — The Foreign Sector: Balance of Payments and Capital Flows*
