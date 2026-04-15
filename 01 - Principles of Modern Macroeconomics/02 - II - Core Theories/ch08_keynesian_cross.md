# Chapter 8 — The Keynesian Cross and Multiplier Effect: Understanding Fiscal Policy

> *"The fundamental psychological law ... is that men are disposed to increase their consumption as their income increases, but not by as much as the increase in their income."*
> — Keynes, *The General Theory*, 1936

---

The AS–AD model of Chapter 7 tells us that a rightward shift in aggregate demand raises short-run output. But by how much? And through what mechanism? The Keynesian cross is the simplest model that answers these questions with precision. It abstracts from the interest rate, from supply-side constraints, and from forward-looking behavior — limitations we will address in subsequent chapters — but its simplicity makes the logic of the **multiplier** transparent in a way that more complex models can obscure. The multiplier is the proposition that an initial increase in autonomous spending generates a larger total increase in income. Understanding what determines its size, what limits it, and under what conditions it collapses toward zero is essential for evaluating any real-world fiscal policy.

The chapter proceeds from the basic two-sector model to successively richer environments: proportional taxes (introducing automatic stabilizers), the open economy (introducing import leakages), and finally an account of the model's structural limitations and the extensions that address them.

---

## 8.1 The Basic Model

### Setup

The Keynesian cross rests on three assumptions: (i) the price level is fixed, so we are unambiguously in the short run with nominal rigidity; (ii) the economy has spare productive capacity, so that demand-determined output is feasible below the supply constraint; and (iii) there is a linear relationship between income and consumption.

The **consumption function** is:

$$C = a + b(Y - T), \quad a > 0,\; b \in (0,1),$$

where $a$ is **autonomous consumption** (the amount households would consume if disposable income were zero, financed by asset drawdowns or borrowing), $b$ is the **marginal propensity to consume**, and $Y - T$ is disposable income after lump-sum taxes $T$.

**Definition (Marginal Propensity to Consume).** The **MPC** $b = \partial C / \partial(Y - T)$ is the fraction of each additional dollar of disposable income spent on consumption. It satisfies $b \in (0,1)$: households spend more when income rises but save the remainder. The Keynesian assumption $0 < b < 1$ is the foundation of the multiplier mechanism; without it — if $b = 1$ (full spending) or $b = 0$ (full saving) — the multiplier is either infinite or one.

Investment is exogenous at $\bar{I}$, government spending at $\bar{G}$, and taxes at $\bar{T}$. Goods-market equilibrium (actual output equals desired expenditure):

$$Y = C + I + G = a + b(Y - \bar{T}) + \bar{I} + \bar{G}.$$

Solving for equilibrium income:

$$Y^* = \frac{1}{1-b}\bigl(a - b\bar{T} + \bar{I} + \bar{G}\bigr) = \kappa\,\bar{A},$$

where $\kappa \equiv 1/(1-b)$ is the **Keynesian multiplier** and $\bar{A} = a - b\bar{T} + \bar{I} + \bar{G}$ is total autonomous expenditure.

### The Round-by-Round Logic

The geometric series derivation makes the mechanism transparent. An initial increase $\Delta G$ in government spending creates $\Delta G$ of new income for construction workers, government contractors, or civil servants in the first round. They spend fraction $b$ on consumption, creating $b\Delta G$ of income for retailers, restaurants, and other consumer-goods producers in the second round. Those producers spend $b$ of their new income, creating $b^2\Delta G$ in the third round — and so on through infinitely many rounds. Summing:

$$\Delta Y^* = \Delta G(1 + b + b^2 + b^3 + \cdots) = \frac{\Delta G}{1-b} = \kappa\,\Delta G.$$

With $b = 0.75$, $\kappa = 4$: each dollar of government spending raises equilibrium income by four dollars. With $b = 0.5$, $\kappa = 2$. The multiplier exists because the circular flow of income means every dollar of spending becomes someone else's income.

---

## 8.2 Three Multipliers and Their Implications

### The Government Spending Multiplier

$$\kappa_G = \frac{\partial Y^*}{\partial \bar{G}} = \frac{1}{1-b}.$$

A one-dollar increase in government purchases raises equilibrium income by $1/(1-b)$ dollars, exceeding one because each round of induced consumption creates further income.

### The Tax Multiplier

$$\kappa_T = \frac{\partial Y^*}{\partial \bar{T}} = \frac{-b}{1-b}.$$

A one-dollar tax cut raises income by $b/(1-b)$ dollars. This is smaller in absolute value than $\kappa_G$ by the factor $b$: a tax cut delivers income to households who then spend fraction $b$ in the first round, whereas government spending delivers a full dollar of direct demand in the first round. Equivalently, the tax multiplier is the spending multiplier multiplied by the MPC, because the fiscal instrument works through household consumption rather than directly.

This asymmetry has policy implications. If the goal is maximum short-run stimulus per dollar of fiscal cost, government spending dominates tax cuts under the Keynesian cross logic. If the goal is to minimize the government's direct role in the economy while providing stimulus, tax cuts are preferred at the cost of a smaller multiplier.

### The Balanced Budget Multiplier: Haavelmo's Theorem

Suppose government spending and taxes both increase by the same amount: $\Delta\bar{G} = \Delta\bar{T} = \Delta B$. The combined effect:

$$\Delta Y^* = \kappa_G \Delta B + \kappa_T \Delta B = \frac{1}{1-b}\Delta B - \frac{b}{1-b}\Delta B = \frac{1-b}{1-b}\Delta B = \Delta B.$$

**Definition (Haavelmo's Theorem).** The **balanced budget multiplier** is exactly one, regardless of the MPC: a balanced fiscal expansion raises income one-for-one with the spending increase. The intuition: the government spends every dollar it collects in taxes ($\Delta G = \Delta T$), while the private sector would have saved fraction $1-b$ of the tax payment. The net demand stimulus equals the amount the private sector would have saved, which is $(1-b)\Delta T$. This times the spending multiplier $1/(1-b)$ gives $\Delta T = \Delta B$. Even a strictly balanced fiscal expansion raises income — a counterintuitive but logically rigorous result.

---

## 8.3 Proportional Taxes and Automatic Stabilizers

### The Proportional Tax Multiplier

The lump-sum tax of Section 8.1 is a simplification. In practice, taxes are approximately proportional to income: $T = tY$, where $t \in (0,1)$ is the marginal income tax rate. Substituting:

$$Y = a + b(Y - tY) + \bar{I} + \bar{G} = a + b(1-t)Y + \bar{I} + \bar{G}.$$

Solving:

$$Y^* = \frac{a + \bar{I} + \bar{G}}{1 - b(1-t)}, \qquad \kappa_G^{prop} = \frac{1}{1 - b(1-t)}.$$

Since $b(1-t) < b$, the proportional-tax multiplier is smaller than the lump-sum-tax multiplier. Each round of additional income generates tax revenue that leaks out of the consumption stream, dampening the multiplier.

### Automatic Stabilizers

This leakage is exactly what makes proportional taxation an **automatic stabilizer**.

**Definition (Automatic Stabilizer).** An **automatic stabilizer** is any fiscal mechanism that reduces the amplitude of output fluctuations in response to demand shocks without requiring discretionary policy changes. Proportional income taxes automatically reduce tax revenue when income falls (supporting household consumption) and increase revenue when income rises (dampening the boom). No legislative action is required; the stabilization is built into the tax code.

The effective multiplier with automatic stabilizers is $\kappa_G^{prop} = 1/[1-b(1-t)]$. With $b = 0.8$ and $t = 0.3$: $\kappa_G^{prop} = 1/[1-0.56] = 2.3$ compared to $\kappa_G = 5$ without taxes. Automatic stabilizers substantially reduce the multiplier — and hence the volatility of income in response to any given shock to autonomous expenditure.

Other important automatic stabilizers include unemployment insurance (UI payments rise automatically when unemployment increases, sustaining the income and consumption of displaced workers), means-tested welfare benefits (which rise automatically when incomes fall), and corporate tax collections (which fall in recessions when profits are compressed). The strength of automatic stabilizers differs substantially across countries: European economies, with larger government sectors and more progressive tax schedules, have stronger automatic stabilization than the United States. Blanchard and Perotti (2002) estimate that the EU's automatic stabilizers are roughly twice as powerful as those in the U.S. in terms of the fraction of a GDP shock they offset.

---

## 8.4 The Open Economy Multiplier

### Import Leakages

In an open economy, some of each round of induced consumption falls on imported goods rather than domestic output, reducing the fiscal multiplier further. With import function $M = m_0 + m_Y Y$, where $m_Y \in (0,1)$ is the **marginal propensity to import**, goods-market equilibrium becomes:

$$Y = a + b(Y - \bar{T}) + \bar{I} + \bar{G} + \bar{X} - m_0 - m_Y Y.$$

Solving for equilibrium output:

$$Y^* = \frac{1}{1-b+m_Y}\bigl(a - b\bar{T} + \bar{I} + \bar{G} + \bar{X} - m_0\bigr), \qquad \kappa_G^{open} = \frac{1}{1-b+m_Y}.$$

The open-economy multiplier $\kappa_G^{open} < \kappa_G^{closed}$ because each round of income now generates both domestic consumption $b$ and imports $m_Y$; only the domestic consumption component circulates within the economy.

### Quantitative Importance

For a small open economy with a high propensity to import ($m_Y = 0.3$) and a moderate MPC ($b = 0.75$): $\kappa_G^{open} = 1/(0.25 + 0.30) = 1.82$, compared to the closed-economy $\kappa_G = 4.0$. International leakages more than halve the multiplier. For very small, very open economies — such as Luxembourg, Singapore, or Ireland — the import propensity may approach 0.5 or higher, and fiscal multipliers approach unity or below.

### The Foreign Trade Multiplier

The open-economy framework also allows analysis of the **foreign trade multiplier**: the effect on domestic output of an increase in foreign income $Y^*$, which raises exports $\bar{X}$. A 1-dollar increase in $\bar{X}$ raises domestic output by $1/(1-b+m_Y)$ — the same formula as the government spending multiplier. This implies that domestic fiscal expansion in one country partially stimulates the economies of its trading partners through higher imports, and conversely that foreign fiscal expansion spills over to the domestic economy through higher exports. This **fiscal externality** motivates international policy coordination in open-economy environments.

---

## 8.5 The Limitations of the Keynesian Cross

The Keynesian cross is analytically tractable and pedagogically essential, but it rests on four assumptions that are too strong to be the basis for realistic policy evaluation.

### Fixed Prices

The model assumes a fixed price level, treating output as demand-determined with no supply constraints. But a demand stimulus eventually encounters supply-side constraints as the economy approaches potential output, raising prices and inducing the AS-AD dynamics of Chapter 7. Rising prices reduce real money balances, raise interest rates, and reduce private spending — partially crowding out the fiscal stimulus even before the IS-LM mechanism operates.

### Exogenous Investment

Treating investment as exogenous ignores the fundamental feedback from income to financial markets to investment. When income rises, money demand rises, driving up interest rates (in the IS-LM model of Chapter 9), which discourages investment. This **crowding out** reduces the effective multiplier below $\kappa_G$. The degree of crowding out depends on the interest sensitivity of investment and money demand — parameters that vary across monetary policy regimes and business cycle states.

### Static Expectations

The Keynesian cross assumes that households consume a fixed fraction of current income without any consideration of the future. Forward-looking households, however, will recognize that a temporary government spending increase must eventually be reversed — either through higher future taxes or spending cuts — and will save more today to meet those future tax obligations. In the extreme case of **Ricardian equivalence** [Ch. 22], forward-looking households with no liquidity constraints and perfect foresight fully anticipate the future tax increase and reduce consumption dollar-for-dollar today, driving the fiscal multiplier to zero regardless of $b$.

### Absence of Intertemporal Optimization

The behavioral consumption function $C = a + b(Y-T)$ is assumed rather than derived. Chapter 11 shows that when the consumption function is derived from intertemporal optimization, the MPC depends on whether income changes are permanent or transitory (permanent changes have higher MPCs than transitory changes under the permanent income hypothesis), on whether households face liquidity constraints (constrained households have MPCs near one, unconstrained households have MPCs determined by their Euler equations), and on the degree of precautionary saving. The Keynesian MPC $b$ is an amalgam of these effects and will shift when the underlying economic environment changes — a version of the Lucas critique.

Each of these extensions enriches the framework and explains why empirical estimates of fiscal multipliers in Chapter 28 range from near zero to above two: the realized multiplier depends on which forces dominate in the specific economic environment under consideration.

---

## Chapter Summary

- The **Keynesian cross** assumes fixed prices and exogenous investment; equilibrium income is $Y^* = [1/(1-b)](a - b\bar{T} + \bar{I} + \bar{G})$ with multiplier $\kappa = 1/(1-b)$.

- **Three multipliers**: government spending $\kappa_G = 1/(1-b)$; tax $\kappa_T = -b/(1-b)$ (smaller in absolute value by factor $b$); balanced budget = 1 (Haavelmo's theorem, independent of $b$).

- **Proportional taxes** reduce the multiplier to $\kappa_G^{prop} = 1/[1-b(1-t)]$ and create **automatic stabilizers** — counter-cyclical fiscal impulses requiring no legislative action. European automatic stabilizers are roughly twice as powerful as those in the U.S.

- **Open economy multiplier** $\kappa_G^{open} = 1/(1-b+m_Y)$ is smaller than the closed-economy multiplier by the import propensity $m_Y$; for small, open economies, multipliers may approach unity. Fiscal expansions create positive demand spillovers to trading partners through higher imports.

- **Four limitations** motivate subsequent models: fixed prices (→ AS-AD, Ch. 7), exogenous investment (→ IS-LM, Ch. 9), static expectations (→ rational expectations, Ch. 15–16), and the assumed consumption function (→ consumption theory, Ch. 11, Ricardian equivalence, Ch. 22).

---

*Next: Chapter 9 — The IS–LM Model*
