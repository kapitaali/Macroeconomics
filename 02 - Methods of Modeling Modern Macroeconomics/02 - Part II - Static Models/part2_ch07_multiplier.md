# Chapter 7: The Multiplier Effect in Closed and Open Economies

*Algebraic Derivation and Extensions*

> *"The multiplier is not one number but a family of numbers, indexed by the model's assumptions about money markets, expectations, and openness."*

**Cross-reference:** *Principles* Ch. 8 (Keynesian cross and all multiplier types); Ch. 26 (open economy, import leakages); Ch. 28 (fiscal policy in practice, ARRA) **[P:Ch.8, P:Ch.26, P:Ch.28]**

---

## 7.1 The Keynesian Cross as a Fixed-Point Problem

*Principles* Chapter 8 introduced the Keynesian cross by means of a diagram: draw the 45-degree line and the expenditure function; equilibrium is their intersection. Here we reformulate this as a fixed-point problem in $Y$, which yields both the existence and uniqueness of equilibrium and the precise multiplier formula from a single piece of mathematics.

The goods-market equilibrium condition is:

$$Y = \mathcal{E}(Y; G, T, \bar{I}), \quad \mathcal{E}(Y) = a + b(Y-T) + \bar{I} + G,$$

where $\mathcal{E}(Y)$ is aggregate planned expenditure as a function of income $Y$ (treating $G$, $T$, and $\bar{I}$ as parameters). Define the **excess demand function**:

$$\phi(Y) \equiv \mathcal{E}(Y) - Y = a + b(Y-T) + \bar{I} + G - Y = \bar{A} - (1-b)Y,$$

where $\bar{A} = a - bT + \bar{I} + G$ is autonomous expenditure. Equilibrium requires $\phi(Y^*) = 0$, i.e., $\bar{A} = (1-b)Y^*$, giving $Y^* = \bar{A}/(1-b) = \kappa_G\bar{A}$.

**The contraction mapping argument:** Consider the iteration $Y_{n+1} = \mathcal{E}(Y_n)$. This defines a sequence $(Y_0, Y_1, Y_2, \ldots)$ starting from any $Y_0$.

**Definition 7.1 (Contraction Mapping).** A function $f: \mathbb{R} \to \mathbb{R}$ is a **contraction** on $\mathbb{R}$ if there exists $\lambda \in [0,1)$ such that $|f(x) - f(y)| \leq \lambda|x - y|$ for all $x, y$.

**Theorem 7.1 (Banach Fixed-Point Theorem, Scalar Version).** If $f$ is a contraction on a complete metric space, then $f$ has a unique fixed point $x^*$ and the iteration $x_{n+1} = f(x_n)$ converges to $x^*$ from any starting point $x_0$.

*Application to the Keynesian cross.* The expenditure function $\mathcal{E}(Y) = \bar{A} + bY$ is a contraction with Lipschitz constant $b = \mathcal{E}'(Y) \in (0,1)$ — precisely the MPC. Since $b < 1$, the iteration $Y_{n+1} = \mathcal{E}(Y_n) = \bar{A} + bY_n$ converges to the unique fixed point $Y^* = \bar{A}/(1-b)$ from any starting income $Y_0$.

The **speed of convergence**: $|Y_n - Y^*| = b^n|Y_0 - Y^*|$. Each "round" of the multiplier process — firms produce, workers earn income, households spend fraction $b$, firms produce again — shrinks the gap from equilibrium by factor $b$. After $n$ rounds, the remaining gap is $b^n$ of the original.

**Definition 7.2 (The Keynesian Multiplier).** The **Keynesian (government spending) multiplier** is:

$$\kappa_G = \frac{\partial Y^*}{\partial G} = \frac{1}{1-b}.$$

It equals the sum of the geometric series $\sum_{n=0}^\infty b^n = 1/(1-b)$: the first-round effect (1 dollar of spending), plus the second-round effect ($b$ dollars of induced consumption), plus the third round ($b^2$), and so on.

---

## 7.2 All Six Multipliers: Complete Algebraic Derivation

The Keynesian cross with proportional income tax rate $t$ (so tax revenue is $T = tY + T_0$ with lump-sum component $T_0$) and transfer payments $TR$:

$$Y = a + b(Y - tY - T_0 + TR) + \bar{I} + G.$$

Collecting terms:
$$Y[1 - b(1-t)] = \underbrace{a - bT_0 + bTR + \bar{I} + G}_{\bar{A}^{eff}}.$$

**Definition 7.3 (Effective Multiplier Denominator).** The expression $1 - b(1-t)$ is the **effective leakage rate**: the fraction of each dollar of income that "leaks out" of the spending stream through saving $(1-b)$ and taxes on marginal income $(bt)$.

$$Y^* = \frac{\bar{A}^{eff}}{1 - b(1-t)} \equiv \kappa^{eff}\bar{A}^{eff}.$$

**Multiplier 1 — Government spending multiplier:**
$$\mu_G = \frac{\partial Y^*}{\partial G} = \frac{1}{1-b(1-t)}.$$

**Multiplier 2 — Lump-sum tax multiplier:**
$$\mu_{T_0} = \frac{\partial Y^*}{\partial T_0} = \frac{-b}{1-b(1-t)}.$$

**Multiplier 3 — Proportional tax rate multiplier:**

$$\mu_t = \frac{\partial Y^*}{\partial t} = \frac{-bY^*}{[1-b(1-t)]^2} \cdot [1-b(1-t)] - \frac{\partial Y^*}{\partial t}\cdot0 = \frac{-bY^*}{1-b(1-t)} < 0.$$

A higher tax rate reduces the effective multiplier and, at any given income level, reduces equilibrium income. The effect depends on $Y^*$ itself (because the tax base is income), making this a nonlinear comparative static.

**Multiplier 4 — Transfer payment multiplier:**
$$\mu_{TR} = \frac{\partial Y^*}{\partial TR} = \frac{b}{1-b(1-t)} = b\cdot\mu_G.$$

Transfers are less stimulative than direct spending: a dollar of transfers raises disposable income by one dollar, of which only fraction $b$ is spent. A dollar of government purchases creates a full dollar of first-round demand.

**Multiplier 5 — The Balanced-Budget Multiplier (Haavelmo's theorem):**

Suppose $\Delta G = \Delta T_0 = \Delta B$ (spending and lump-sum taxes increase equally):

$$\Delta Y^* = \mu_G\Delta B + \mu_{T_0}\Delta B = \frac{1}{1-b(1-t)}\Delta B - \frac{b}{1-b(1-t)}\Delta B = \frac{1-b}{1-b(1-t)}\Delta B.$$

**Theorem 7.2 (Haavelmo's Theorem, General Form).** When government spending and lump-sum taxes both increase by $\Delta B$, output rises by:

$$\Delta Y^* = \frac{1-b}{1-b(1-t)}\Delta B.$$

With $t = 0$ (lump-sum taxes only), this reduces to $\Delta Y^* = \frac{1-b}{1-b}\Delta B = \Delta B$: the balanced budget multiplier is exactly 1.

*Proof (general case, $t = 0$):*

$\Delta Y^* = \kappa_G\Delta B + \mu_{T_0}\Delta B = \frac{1}{1-b}\Delta B - \frac{b}{1-b}\Delta B = \frac{1-b}{1-b}\Delta B = \Delta B.$ $\square$

The economic intuition: government spends every dollar of tax revenue, while households would have saved fraction $(1-b)$. The government therefore contributes more to aggregate demand per dollar taxed than households would have. The net gain equals exactly the saved fraction: $\Delta Y^* = 1 - |\mu_{T_0}| = 1 - b/(1-b) \cdot (1-b) = 1$... but wait — with $t = 0$, $\mu_G = 1/(1-b)$ and $|\mu_{T_0}| = b/(1-b)$, so $\Delta Y^* = [1/(1-b) - b/(1-b)]\Delta B = \Delta B$. The balanced budget multiplier is 1, independent of $b$. A remarkable result.

**Multiplier 6 — Investment multiplier:**
$$\mu_{\bar{I}} = \frac{\partial Y^*}{\partial\bar{I}} = \frac{1}{1-b(1-t)} = \mu_G.$$

An autonomous increase in investment — an animal spirits boom [P:Ch.15.4] — has the same multiplied effect on output as an equivalent increase in government spending.

---

## 7.3 The Open-Economy Multiplier and Import Leakages

In an open economy, some of each round of induced spending falls on imports rather than domestically produced goods [P:Ch.26.5]. This **import leakage** reduces the multiplier.

With import function $IM = m_0 + m_Y Y$ (where $m_Y > 0$ is the marginal propensity to import), the open-economy equilibrium:

$$Y = a + b(Y-T) + \bar{I} + G + X - m_0 - m_Y Y,$$

where $X$ is autonomous exports. Collecting:

$$Y[1 - b(1-t) + m_Y] = \bar{A}^{open}.$$

**Multiplier 7 — Open-economy government spending multiplier:**

$$\mu_G^{open} = \frac{1}{1-b(1-t)+m_Y} < \mu_G^{closed}.$$

The import leakage $m_Y$ in the denominator reduces the multiplier: income generated by fiscal expansion partly flows abroad as import demand, weakening the domestic income circuit.

**Definition 7.4 (Propensity to Spend on Domestic Output).** Define $\tilde{b} = b(1-t) - m_Y$ as the **net marginal propensity to spend on domestic output** — the fraction of each additional dollar of income that re-enters the domestic spending stream. The open-economy multiplier is simply:

$$\mu_G^{open} = \frac{1}{1-\tilde{b}}.$$

This unifies all six multipliers: they all take the form $1/(1-\tilde{b})$ where $\tilde{b}$ accounts for whatever leakages are present (saving, taxes, imports, and in the IS–LM context, the interest rate feedback via crowding out).

**Table of leakages and their effects on multipliers:**

| Leakage source | Reduces denominator by | Effect on $\mu_G$ |
|---|---|---|
| Saving (MPS) | $+(1-b)$ | Reduces multiplier |
| Proportional tax | $+bt$ | Reduces multiplier |
| Import propensity | $+m_Y$ | Reduces multiplier |
| Interest crowding | $+b_r k/h$ | Reduces multiplier (IS–LM) |
| Ricardian saving | $+b\cdot 0 = 0$? | If full Ricardian equiv., $\to 0$ |

---

## 7.4 The Propensity to Spend: Sensitivity Analysis

The multiplier $\kappa = 1/(1-b)$ is extremely sensitive to the MPC $b$ near $b = 1$:

$$\frac{d\kappa}{db} = \frac{1}{(1-b)^2} \to \infty \text{ as } b\to 1.$$

At $b = 0.75$: $\kappa = 4$. At $b = 0.80$: $\kappa = 5$. At $b = 0.90$: $\kappa = 10$. At $b = 0.95$: $\kappa = 20$.

This sensitivity has important implications for policy analysis. The MPC is not a structural constant — it varies across:

1. **Households** by income and wealth: liquidity-constrained households have $b \approx 1$ (they spend every dollar of income); wealthy households have $b \approx 0.3$–$0.5$ (Kaplan et al., 2018, HANK model [P:Ch.25.3]).
2. **Type of fiscal transfer**: lump-sum cash transfers have high MPCs for constrained recipients; corporate tax cuts have low MPCs if they flow to wealthy shareholders.
3. **Expectations about permanence**: permanent income changes have higher MPCs than transitory ones (PIH, [P:Ch.11.2]).

**The elasticity of the multiplier with respect to the MPC:**

$$\varepsilon_{\kappa,b} = \frac{d\kappa}{db}\cdot\frac{b}{\kappa} = \frac{b}{(1-b)^2}\cdot\frac{1-b}{1} = \frac{b}{1-b}.$$

At $b = 0.75$: $\varepsilon = 3$. A 1% increase in the MPC raises the multiplier by 3%. This high elasticity means that the policy effectiveness of fiscal stimulus depends critically on who receives it — and whether recipients are constrained.

---

## 7.5 The Geometric Series Derivation: Round-by-Round Accounting

The multiplier formula $\kappa_G = 1/(1-b) = \sum_{n=0}^\infty b^n$ connects to a concrete narrative about how spending circulates through the economy.

**Round 0:** Government increases spending by $\Delta G = 1$. Output rises by 1.

**Round 1:** The dollar of output becomes income for workers and owners. They spend fraction $b$: induced consumption $\Delta C_1 = b$. Output rises by an additional $b$.

**Round 2:** The additional $b$ of output becomes income. Fraction $b$ is spent: $\Delta C_2 = b^2$. Output rises by $b^2$.

**Round $n$:** $\Delta C_n = b^n$.

**Total:** $\sum_{n=0}^\infty b^n = 1/(1-b) = \kappa_G$.

The sum converges because $b < 1$: each round the increment is smaller by factor $b$, and the infinite sum is finite. If $b = 1$ (no saving), the sum diverges — a $1 spending becomes infinite output. If $b = 0$ (all saving), the sum is 1 — the multiplier equals one because no induced spending occurs.

In APL, the round-by-round accumulation is a natural scan operation:

```apl
⎕IO←0 ⋄ ⎕ML←1

b ← 0.75
n_rounds ← 20

⍝ Round-by-round increments: b^0, b^1, ..., b^(n-1)
increments ← b * ⍳ n_rounds      ⍝ geometric sequence

⍝ Cumulative sum = partial sum of multiplier series
cumulative ← +\ increments        ⍝ scan: running total

⍝ Theoretical limit
kappa_G ← ÷ 1 - b                 ⍝ = 4.0

⍝ Convergence check: how many rounds to reach 99% of limit?
pct_of_limit ← cumulative ÷ kappa_G
n_99pct ← +/ pct_of_limit < 0.99  ⍝ count rounds below 99%
n_99pct    ⍝ ≈ 14 rounds for b=0.75
```

---

## 7.6 Worked Example: The 2009 ARRA Multiplier

*Cross-reference: Principles Ch. 28.2 (ARRA empirical evidence)* **[P:Ch.28.2]**

The American Recovery and Reinvestment Act (2009) provided $787 billion in stimulus over 10 years. A naive Keynesian cross calculation with $b = 0.75$ and no leakages gives $\kappa_G = 4$ — the stimulus should have raised GDP by over $3 trillion. The actual impact was far smaller. Why?

**Calibrating with leakages:**

Let $b = 0.75$, $t = 0.28$ (effective marginal federal plus state income tax rate), $m_Y = 0.12$ (import propensity), $b_r = 1.5$ (investment-interest sensitivity), $h = 4$ (interest semi-elasticity), $k = 0.5$ (income elasticity of money demand).

**IS–LM multiplier (closed):**
$$\mu_G^{closed} = \frac{h}{h(1-b)+b_r k} = \frac{4}{4(0.25)+1.5\times0.5} = \frac{4}{1+0.75} = \frac{4}{1.75} = 2.29.$$

**Adding proportional tax:**
$$\mu_G^{tax} = \frac{1}{1-b(1-t)} = \frac{1}{1-0.75\times0.72} = \frac{1}{0.46} = 2.17.$$

**Adding imports (IS–LM with open economy):**
$$\mu_G^{open,\,IS-LM} = \frac{h}{h(1-b(1-t)+m_Y)+b_r k(1-b(1-t)+m_Y)/(1)}.$$

With the combined denominator $(1-b(1-t)+m_Y) = 0.46 + 0.12 = 0.58$ replacing $(1-b)$:

$$\mu_G^{full} = \frac{h}{h\cdot0.58 + b_r\cdot k} = \frac{4}{2.32+0.75} = \frac{4}{3.07} = 1.30.$$

**At the ELB (estimated):** With $h\to\infty$ (zero crowding out), the ELB multiplier:
$$\mu_G^{ELB} = \frac{1}{1-b(1-t)+m_Y} = \frac{1}{0.58} = 1.72.$$

This is consistent with the Romer–Bernstein (2009) projection of 1.57 and the Nakamura–Steinsson (2014) cross-state estimate of approximately 1.5–2.0 [P:Ch.28.2].

The lesson: the naive multiplier of 4 reflects a closed economy at the liquidity trap with no taxes or imports. Each realistic modification reduces it. Fiscal policy remains stimulative, but the magnitude depends critically on the model calibration.

```apl
⍝ APL — ARRA multiplier under successive leakage assumptions
⎕IO←0 ⋄ ⎕ML←1

b←0.75  ⋄  t←0.28  ⋄  m_Y←0.12  ⋄  b_r←1.5  ⋄  h←4  ⋄  k←0.5

⍝ Progressive multiplier calculations
mu_naive   ← ÷ 1-b                          ⍝ no leakages: 4.0
mu_tax     ← ÷ 1-b×1-t                      ⍝ + proportional tax: 2.17
mu_open    ← ÷ 1-b×1-t)+m_Y                 ⍝ + imports: 1.72 (ELB version)
mu_islm    ← h÷(h×1-b×1-t))+b_r×k          ⍝ + crowding out (closed)
mu_full    ← h÷(h×(1-b×1-t)+m_Y))+b_r×k   ⍝ + imports + crowding out

⍝ Display as a table
labels ← 'Naive' 'Tax' 'Open(ELB)' 'ISLM(closed)' 'Full'
values ← mu_naive mu_tax mu_open mu_islm mu_full
↑ labels ,¨ ⍕¨ 2 ⍕¨ values    ⍝ formatted table
```

---

## 7.7 Programming Exercises

### Exercise 7.1 (APL — Full Multiplier Family)

Write a dfn `all_multipliers ← {b t m_Y ← ⍵ ⋄ ...}` that returns a 6-element vector $(\mu_G, \mu_{T_0}, \mu_{TR}, \mu_{balanced}, \mu_{open}, \mu_{\bar{I}})$ for parameters $(b, t, m_Y)$. Verify: (a) $|\mu_{T_0}| < \mu_G$ always; (b) $\mu_{balanced} = (1-b)/(1-b(1-t)) \leq 1$; (c) $\mu_{TR} = b\cdot\mu_G$.

### Exercise 7.2 (Python — MPC Sensitivity)

```python
import numpy as np, matplotlib.pyplot as plt

b_vals = np.linspace(0.01, 0.99, 500)
# Three multipliers as functions of b (t=0.25, m_Y=0.10)
t, m_Y = 0.25, 0.10
mu_closed = 1 / (1-b_vals)
mu_tax    = 1 / (1-b_vals*(1-t))
mu_open   = 1 / (1-b_vals*(1-t)+m_Y)

fig, ax = plt.subplots()
ax.plot(b_vals, mu_closed, label='Closed, no tax')
ax.plot(b_vals, mu_tax,    label='Closed, with tax t=0.25')
ax.plot(b_vals, mu_open,   label='Open, t=0.25, m_Y=0.10')
ax.axvline(0.75, linestyle='--', color='gray', alpha=0.5, label='b=0.75')
ax.set_ylim(0, 15); ax.set_xlabel('MPC (b)'); ax.set_ylabel('Multiplier')
ax.legend(); plt.title('Keynesian Multiplier vs. MPC'); plt.show()
```

### Exercise 7.3 (Julia — Round-by-Round Convergence)

```julia
b = 0.80; n = 30
increments = b .^ (0:n-1)
cumulative = cumsum(increments)
limit = 1/(1-b)

println("Rounds to 90% of limit: ", findfirst(cumulative ./ limit .>= 0.90))
println("Rounds to 99% of limit: ", findfirst(cumulative ./ limit .>= 0.99))
println("Rounds to 99.9%:        ", findfirst(cumulative ./ limit .>= 0.999))
# Verify: at b=0.8, convergence is slower than b=0.75
```

### Exercise 7.4 (R — ARRA Calibration Sweep)

```r
# Sweep over (b, t, m_Y) to find range of plausible ARRA multipliers
b_vals <- seq(0.5, 0.9, 0.05)
t_vals <- seq(0.20, 0.35, 0.05)
m_Y    <- 0.12

# ELB multiplier: no crowding out
mu_grid <- outer(b_vals, t_vals, function(b, t) 1/(1-b*(1-t)+m_Y))
rownames(mu_grid) <- paste0("b=", b_vals)
colnames(mu_grid) <- paste0("t=", t_vals)
round(mu_grid, 2)
```

### Exercise 7.5 — Haavelmo in the IS–LM Model ($\star$)

Prove that the balanced-budget multiplier in the IS–LM model (with $\Delta G = \Delta T_0 = \Delta B$) equals:

$$\mu_{BB}^{IS-LM} = \frac{(1-b)h}{h(1-b)+b_r k}.$$

Show that this is strictly less than 1 (unlike the Keynesian cross balanced-budget multiplier of 1) because the crowding-out effect also applies to the spending expansion. What is the limit as $h\to\infty$ (liquidity trap)? As $h\to 0$ (classical case)?

### Exercise 7.6 — Negative Multiplier? ($\star\star$)

The "expansionary austerity" hypothesis (Alesina and Ardagna, 2010) suggests that fiscal consolidation can raise output through confidence effects [P:Ch.28.5]. Formally, model confidence as $\bar{I}(G) = \bar{I}_0 - \gamma G$ where $\gamma > 0$ (higher government spending reduces private investment via an uncertainty channel). (a) Derive the multiplier $\partial Y^*/\partial G$ when confidence effects are present. (b) For what $\gamma$ does the multiplier turn negative (contractionary fiscal expansion)? (c) Calibrate with $b = 0.75$, $t = 0.25$, $\gamma = 0.5$: does the multiplier change sign? Interpret.

---

## 7.8 Chapter Summary

**Key results:**

- The Keynesian cross equilibrium is a **fixed point** of the expenditure function $\mathcal{E}(Y)$; the contraction mapping theorem guarantees existence and uniqueness when $b < 1$.
- The **six multipliers** for the closed economy with proportional taxes:
  - Spending: $\mu_G = 1/[1-b(1-t)]$
  - Lump-sum tax: $\mu_{T_0} = -b/[1-b(1-t)]$
  - Transfer: $\mu_{TR} = b/[1-b(1-t)]$
  - Balanced budget: $\mu_{BB} = (1-b)/[1-b(1-t)] \leq 1$
  - Investment: $\mu_{\bar{I}} = \mu_G$
  - Open economy: $\mu_G^{open} = 1/[1-b(1-t)+m_Y]$
- **Haavelmo's theorem**: with lump-sum taxes and $t = 0$, the balanced-budget multiplier is exactly 1, independent of the MPC.
- Every multiplier takes the form $1/(1-\tilde{b})$ where $\tilde{b}$ is the net marginal propensity to spend on domestic output, unifying all cases.
- In APL: the round-by-round accumulation is `+\ b * ⍳ n_rounds` (scan of geometric sequence); the full multiplier table is generated via `∘.f` outer product over parameter grids.

**Connections forward:** Chapter 9 derives the ELB multiplier — the analogous formula for the dynamic NK model at the zero lower bound, where crowding out is absent and forward-looking inflation expectations amplify the fiscal effect beyond $1/(1-\tilde{b})$.

---

*Next: Chapter 8 — The AD–AS Model in Equation Form*
