# Chapter 33: Income Inequality and Macroeconomic Dynamics

*Pareto Laws, Gini Coefficients, and the URE Framework*

> *"The Gini coefficient is a number between 0 and 1. Understanding where it comes from, what drives it, and how policy affects it — that requires a model."*

**Cross-reference:** *Principles* Ch. 38 (inequality, Gini coefficient, Piketty's $r > g$, distributional consequences of policy); Ch. 25 (HANK, heterogeneous MPC, redistribution channel) **[P:Ch.38, P:Ch.25]**

---

## 33.1 Why Inequality Matters for Macroeconomics

Chapter 32 built the tools for solving heterogeneous-agent models. This chapter uses them to analyze inequality — not as a welfare concern in isolation, but as a determinant of aggregate dynamics. Three mechanisms connect inequality to macroeconomic fluctuations:

**1. Heterogeneous MPC.** Poor households have high MPCs (near 1); wealthy households have low MPCs (near 0). A redistribution from rich to poor raises aggregate consumption. The size of the fiscal multiplier depends on the wealth distribution.

**2. Distributional channel of monetary policy.** A rate cut benefits net debtors (who spend more) and hurts net savers (who spend less). The aggregate effect depends on the distribution of net financial wealth — the subject of Auclert's (2019) URE framework developed in Section 33.4.

**3. Inequality and savings.** Higher inequality concentrates wealth in households with lower MPCs, reducing aggregate demand — a potential explanation for secular stagnation [P:Ch.39.4].

---

## 33.2 Pareto Distributions and Power Laws

**Definition 33.1 (Pareto Distribution).** A random variable $X$ has a **Pareto distribution** with scale $x_m > 0$ and shape $\zeta > 0$ if:

$$P(X > x) = \left(\frac{x_m}{x}\right)^\zeta \quad \text{for } x \geq x_m.$$

The PDF: $f(x) = \zeta x_m^\zeta / x^{\zeta+1}$. The mean: $\mathbb{E}[X] = \zeta x_m/(\zeta-1)$ (finite iff $\zeta > 1$). The variance: finite iff $\zeta > 2$.

**Power law tails in wealth data:** Empirically, the upper tail of the wealth distribution follows a Pareto law with $\zeta \approx 1.5$ for the U.S. The top 1% share $s_1$ and the Pareto exponent $\zeta$ are related by $s_1 \approx 0.01^{1-1/\zeta}$.

**Theorem 33.1 (Gini Coefficient for a Pareto Distribution).** For $X \sim \text{Pareto}(x_m, \zeta)$ with $\zeta > 1$:

$$\text{Gini} = \frac{1}{2\zeta - 1}.$$

*Proof.* The Lorenz curve $L(p)$ gives the share of total income held by the bottom fraction $p$. For Pareto: $L(p) = 1 - (1-p)^{1-1/\zeta}$. The Gini coefficient:

$$G = 1 - 2\int_0^1 L(p)\,dp = 1 - 2\int_0^1[1-(1-p)^{1-1/\zeta}]\,dp.$$

Setting $u = 1-p$: $\int_0^1(1-p)^{1-1/\zeta}dp = \int_0^1 u^{1-1/\zeta}du = 1/(2-1/\zeta) = \zeta/(2\zeta-1)$.

So $G = 1 - 2[1 - \zeta/(2\zeta-1)] = 1 - 2(\zeta-1)/(2\zeta-1) = 1/(2\zeta-1)$. $\square$

**Applications:** For U.S. wealth ($\zeta \approx 1.5$): $G \approx 1/(3-1) = 0.5$ for the Pareto tail alone. The overall Gini (including the non-Pareto middle of the distribution) is higher, approximately 0.87.

---

## 33.3 The Gini Coefficient: General Formula

For an arbitrary distribution (not necessarily Pareto), the Gini coefficient is computed from the Lorenz curve.

**Definition 33.2 (Lorenz Curve).** The **Lorenz curve** $L(p)$ maps the cumulative population share $p$ to the cumulative wealth (or income) share held by the bottom $p$ fraction:

$$L(p) = \frac{\int_0^p F^{-1}(q)\,dq}{\int_0^1 F^{-1}(q)\,dq} = \frac{\int_0^{F^{-1}(p)}x\,dF(x)}{\mathbb{E}[X]},$$

where $F^{-1}(p)$ is the $p$-th quantile of the wealth distribution.

**Definition 33.3 (Gini Coefficient).** The **Gini coefficient** is twice the area between the 45-degree line and the Lorenz curve:

$$G = 1 - 2\int_0^1 L(p)\,dp = \frac{2\text{Cov}(X, F(X))}{\mathbb{E}[X]}.$$

**Theorem 33.2 (Gini as Covariance).** For a positive random variable $X$ with CDF $F$:

$$G = \frac{2}{\mathbb{E}[X]}\text{Cov}(X, F(X)).$$

*Proof.* $\text{Cov}(X, F(X)) = \mathbb{E}[XF(X)] - \mathbb{E}[X]\mathbb{E}[F(X)]$. Since $\mathbb{E}[F(X)] = 1/2$ (by the probability integral transform), $\text{Cov}(X,F(X)) = \mathbb{E}[XF(X)] - \mathbb{E}[X]/2$. Rewriting the Gini: $G = (1/\mathbb{E}[X])\mathbb{E}[X(2F(X)-1)] = (2/\mathbb{E}[X])(\mathbb{E}[XF(X)] - \mathbb{E}[X]/2) = (2/\mathbb{E}[X])\text{Cov}(X,F(X))$. $\square$

**Sample formula:** For wealth observations $a_1 \leq a_2 \leq \cdots \leq a_n$ (sorted):

$$G = \frac{2}{n\bar{a}}\sum_{i=1}^n i\cdot a_i - \frac{n+1}{n}.$$

In APL: `G ← (2÷n×a_bar) × +/ (⍳n) × a_sorted) - (n+1)÷n`

---

## 33.4 The Auclert (2019) URE Framework

*Cross-reference: Principles Ch. 25.3 (HANK redistribution channel)* **[P:Ch.25.3]**

Auclert (2019) provides a powerful framework for decomposing the aggregate consumption response to a monetary policy shock into components that depend on the wealth distribution. The key concept is **Unhedged Interest Rate Exposure (URE)**.

**Definition 33.4 (Unhedged Interest Rate Exposure).** Each household's **URE** is:

$$\text{URE}_i = \underbrace{c_i - y_i}_{\text{net borrower/saver}} + \underbrace{a_i^{maturing} - l_i^{maturing}}_{\text{maturing assets minus liabilities}},$$

where $a_i^{maturing}$ and $l_i^{maturing}$ are assets and liabilities maturing in the current period (that will be rolled over at the new interest rate).

**Intuition:** A household with positive URE (net saver) benefits from a rate increase; a household with negative URE (net borrower) is hurt. A rate cut redistributes purchasing power from high-URE to low-URE households.

**Theorem 33.3 (Auclert URE Sufficient Statistic).** The aggregate consumption response to a unit change in the real interest rate $dr$ is:

$$d\bar{C} = \underbrace{-\frac{1}{\sigma}\bar{C}dr}_{\text{intertemporal substitution}} + \underbrace{\text{Cov}\!\left(\text{MPC}_i,\;\text{URE}_i\right)\cdot dr}_{\text{redistribution channel}},$$

where $\sigma$ is the EIS, $\bar{C}$ is aggregate consumption, and $\text{MPC}_i$ is household $i$'s marginal propensity to consume.

*Proof sketch.* Each household's consumption responds to a rate change through: (1) the Euler equation effect ($-\dot{c}/c = (r-\rho)/\sigma$, so $dc_i/dr \propto -c_i/\sigma$); and (2) the income effect from interest rate changes on assets and liabilities at maturity ($dc_i/dr \propto \text{MPC}_i \cdot \text{URE}_i$). Aggregating: $d\bar{C} = -\bar{C}/\sigma \cdot dr + \mathbb{E}[\text{MPC}_i\cdot\text{URE}_i] \cdot dr$. Since $\mathbb{E}[\text{MPC}_i]=\bar{\text{MPC}}$ and $\mathbb{E}[\text{URE}_i]=0$ (economy-wide, interest payments are transfers between agents), $\mathbb{E}[\text{MPC}_i\cdot\text{URE}_i] = \text{Cov}(\text{MPC}_i, \text{URE}_i)$. $\square$

**Key insight:** The redistribution channel $\text{Cov}(\text{MPC},\text{URE})$ is **positive** if poor households (high MPC) are net borrowers (negative URE) and wealthy households (low MPC) are net savers (positive URE). In this case, a rate cut redistributes to high-MPC households, amplifying the aggregate consumption response beyond what the representative-agent Euler equation would predict.

**Calibration of the redistribution channel:** For the U.S.:
- $\text{Cov}(\text{MPC}, \text{URE}) \approx -0.07$ (Auclert 2019 estimate).
- Average MPC $\approx 0.25$.
- So the redistribution channel contributes $-0.07 \cdot dr$ to aggregate consumption.
- For a 1pp rate cut ($dr = -0.01$): redistribution channel contributes $+0.07\%$ of aggregate consumption — substantial relative to the total effect.

---

## 33.5 Decomposing the Gini by Income Source

For understanding the contribution of different income sources (labor, capital, transfers) to total inequality, the **Gini decomposition by source** is useful.

**Theorem 33.4 (Gini Decomposition).** If total income $Y_i = Y_i^L + Y_i^K + Y_i^T$ (labor + capital + transfers), the Gini of total income:

$$G_Y = \sum_k \frac{\mathbb{E}[Y_i^k]}{\mathbb{E}[Y_i]} \cdot G_{Y^k} \cdot \rho_k,$$

where $G_{Y^k}$ is the concentration coefficient of source $k$ and $\rho_k$ is the rank correlation between source $k$ and total income.

---

## 33.6 Worked Example: Gini and URE in a Calibrated HANK Model

*Cross-reference: Principles Ch. 25.3, Ch. 38* **[P:Ch.25.3, P:Ch.38]**

```python
import numpy as np
from scipy.interpolate import interp1d

# Using the Aiyagari wealth distribution from Chapter 32
# Here: simplified version with known distribution
np.random.seed(42)
N = 50000

# Generate a wealth distribution with Pareto upper tail
# Mix: 60% near zero (hand-to-mouth), 40% Pareto
frac_htm = 0.45  # hand-to-mouth fraction
n_htm = int(frac_htm * N)
n_pareto = N - n_htm

wealth_htm = np.random.exponential(0.5, n_htm)  # near-zero wealth
zeta_pareto = 1.5  # Pareto exponent
wealth_pareto = np.random.pareto(zeta_pareto, n_pareto) * 10.0 + 5.0

wealth = np.concatenate([wealth_htm, wealth_pareto])
wealth = np.sort(wealth)

# Gini coefficient
n = len(wealth); mean_w = np.mean(wealth)
G_sample = (2/(n*mean_w)) * np.sum((np.arange(1, n+1))*wealth) - (n+1)/n
print(f"Wealth Gini: {G_sample:.3f}  (data: ~0.87)")

# Pareto tail fit
threshold = np.percentile(wealth, 90)  # fit Pareto to top 10%
tail = wealth[wealth > threshold]
# MLE for Pareto: ζ_hat = n / sum(log(x/x_m))
zeta_hat = len(tail) / np.sum(np.log(tail/threshold))
print(f"Pareto exponent (top 10%): ζ = {zeta_hat:.3f}")
print(f"Implied Gini for Pareto tail: {1/(2*zeta_hat-1):.3f}")

# Lorenz curve
cumshare_pop = np.linspace(0, 1, n)
cumshare_wealth = np.cumsum(wealth) / np.sum(wealth)

# Top shares
for pct in [1, 5, 10, 50]:
    share = np.sum(wealth[int((1-pct/100)*n):]) / np.sum(wealth)
    print(f"Top {pct}% wealth share: {share*100:.1f}%")

# URE computation
# Assign MPC inversely proportional to wealth (hand-to-mouth MPC=1, wealthy MPC→0)
mpc = np.maximum(0.05, 1 - 0.9*cumshare_wealth)  # declining in wealth rank

# URE: net borrowers (bottom 60%) negative; net savers (top 40%) positive
ure_median = 0.0
ure = (cumshare_pop - 0.6) * 2  # roughly: bottom 60% negative URE

# Auclert redistribution channel
cov_mpc_ure = np.cov(mpc, ure)[0,1]
sigma_eis = 1.0; C_bar = np.mean(wealth) * 0.1  # consumption ≈ 10% of wealth
dr = -0.01  # 100bp rate cut
dC_euler = -C_bar/sigma_eis * dr
dC_redist = cov_mpc_ure * dr
print(f"\nResponse to 100bp rate cut:")
print(f"  Intertemporal substitution: {dC_euler*100:.3f}% of C")
print(f"  Redistribution (URE) channel: {dC_redist*100:.3f}% of C")
print(f"  Total dC: {(dC_euler+dC_redist)*100:.3f}% of C")
print(f"  Cov(MPC, URE): {cov_mpc_ure:.4f}")
```

```julia
using Statistics, LinearAlgebra

# Gini coefficient
function gini(wealth)
    n = length(wealth); a = sort(wealth); μ = mean(a)
    (2/(n*μ))*sum(i*a[i] for i in 1:n) - (n+1)/n
end

# Pareto fit
function fit_pareto_mle(x, x_min)
    n = count(x .>= x_min); tail = x[x .>= x_min]
    n / sum(log.(tail./x_min))
end

np = 50000
wealth = sort(vcat(rand(Exponential(0.5), Int(0.45*np)),
                   rand(Pareto(1.5), Int(0.55*np)).*10.0.+5.0))
println("Wealth Gini: $(round(gini(wealth), digits=3))")

thresh = quantile(wealth, 0.90)
zeta = fit_pareto_mle(wealth, thresh)
println("Pareto ζ (top 10%): $(round(zeta,digits=3))")
println("Pareto Gini: $(round(1/(2*zeta-1), digits=3))")
```

```r
# R — Gini coefficient and Lorenz curve
wealth <- sort(c(rexp(22500, 2), rpareto(27500, 1.5)*10+5))
n <- length(wealth); mu <- mean(wealth)
G <- (2/(n*mu))*sum((1:n)*wealth) - (n+1)/n
cat(sprintf("Gini: %.3f\n", G))

# Lorenz curve plot
cum_pop <- (1:n)/n; cum_wealth <- cumsum(wealth)/sum(wealth)
plot(cum_pop, cum_wealth, type='l', col='blue',
     xlab='Population share', ylab='Wealth share', main='Lorenz Curve')
abline(0,1, lty=2, col='red'); legend('topleft', c('Lorenz','Equality'), col=c('blue','red'), lty=1:2)
```

---

## 33.7 Programming Exercises

### Exercise 33.1 (APL — Sample Gini)

Write a dfn `gini_coef ← {a ← a[⍋a] ⋄ n←≢a ⋄ mu←(+/a)÷n ⋄ (2÷n×mu)×+/(⍳n)×a) - (n+1)÷n}`. (a) Verify it equals `(2×mu_mean)+/` on the sorted wealth distribution from Chapter 32's simulation. (b) Compute the Gini as a function of the borrowing limit $\underline{a}$: tighter constraints force more households to be hand-to-mouth, raising the Gini.

### Exercise 33.2 (Python — Piketty $r > g$ and Wealth Concentration)

Piketty's inequality dynamic: $\dot{w}_i = (r-g)w_i + y_i - c_i$ implies that when $r > g$, wealth grows faster than the economy, concentrating at the top. (a) Simulate 1000 households for 100 periods with $r = 0.05$ and $g \in \{0.01, 0.02, 0.03, 0.04, 0.05\}$. (b) Compute the Gini and top-1% share for each $g$. (c) Plot the relationship between $r-g$ and the Gini. (d) Calibrate to U.S. data: what $r-g$ matches the observed top-1% share of 38%?

### Exercise 33.3 (Julia — Decomposing the Rate Cut Effect)

```julia
using Distributions, Statistics

# Simulate heterogeneous households with different MPC and URE
N = 10000; np_rng = 42
Random.seed!(np_rng)

# Wealth distribution (lognormal lower tail + Pareto upper)
a_grid = sort(vcat(rand(LogNormal(0, 1), N÷2), rand(Pareto(1.5), N÷2).*5.0.+2.0))
n = length(a_grid)

# MPC: declining in wealth rank (Carroll-type buffer stock)
rank = (1:n)/n
mpc = max.(0.02, 1.0 .- 0.95.*rank.^0.5)

# URE: negative for borrowers (bottom 55%), positive for savers
ure = (rank .- 0.55) .* 2.0 .* std(a_grid)

# Decompose aggregate dC for dr = -1% (rate cut)
dr = -0.01; sigma_eis = 1.0; C_bar = mean(mpc .* a_grid)

# Aggregate effect components
sub_effect = -C_bar/sigma_eis * dr
redist_effect = cov(mpc, ure) * dr
income_effect = mean(mpc .* a_grid .* dr)  # direct income from interest

println("Rate cut decomposition (Auclert 2019):")
println("  Intertemporal substitution: $(round(sub_effect,digits=4))")
println("  Redistribution (URE):       $(round(redist_effect,digits=4))")
println("  Total:                       $(round(sub_effect+redist_effect,digits=4))")
println("  Cov(MPC, URE): $(round(cov(mpc,ure),digits=4))")
println("  Mean MPC: $(round(mean(mpc),digits=3))")
```

### Exercise 33.4 — Decomposing Wealth Inequality ($\star$)

In the Aiyagari model from Chapter 32, decompose the wealth Gini into: (a) the contribution of income heterogeneity (different $y_i$ realizations); (b) the contribution of precautionary saving (households accumulating buffers at different rates); (c) the contribution of the borrowing constraint (households bunched at $a = \underline{a}$). Use the Gini decomposition by source (Theorem 33.4). Compare to the observed U.S. decomposition: income inequality explains about 30% of wealth inequality; the rest comes from savings behavior heterogeneity and the binding borrowing constraint.

---

## 33.8 Chapter Summary

**Key results:**

- **Pareto distribution** $P(X>x) = (x_m/x)^\zeta$ has Gini $= 1/(2\zeta-1)$ (Theorem 33.1); U.S. wealth data has $\zeta \approx 1.5$, giving Pareto Gini $\approx 0.5$ (upper tail only).
- **Gini coefficient** $G = 2\text{Cov}(X, F(X))/\mathbb{E}[X]$ (Theorem 33.2); sample formula: $G = (2/(n\bar{a}))\sum_i i\cdot a_{(i)} - (n+1)/n$.
- The **Auclert URE framework** decomposes the aggregate consumption response to a rate change into: (1) intertemporal substitution $(-C/\sigma)\cdot dr$; and (2) redistribution $\text{Cov}(\text{MPC},\text{URE})\cdot dr$ (Theorem 33.3). The redistribution channel amplifies policy if poor households are net borrowers.
- For U.S. data: $\text{Cov}(\text{MPC},\text{URE}) > 0$, so rate cuts are amplified by redistribution from low-MPC savers to high-MPC borrowers.
- In APL: Gini is `(2÷n×mu)×(+/(⍳n)×a_sorted) - (n+1)÷n`; Lorenz curve as cumulative sum divided by total `+\a_sorted ÷ +/a_sorted`.

*Next: Chapter 34 — Network Models for Financial Contagion and Systemic Risk*
