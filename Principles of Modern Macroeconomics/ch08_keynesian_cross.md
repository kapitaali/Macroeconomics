# Chapter 8 — The Keynesian Cross and Multiplier Effect: Understanding Fiscal Policy

> *"The fundamental psychological law ... is that men are disposed to increase their consumption as their income increases, but not by as much as the increase in their income."*
> — Keynes, *The General Theory*, 1936

---

The AS–AD model of Chapter 7 tells us that a rightward shift in aggregate demand raises short-run output. But by how much? And through what mechanism? The Keynesian cross is the simplest model that answers these questions. It is not the most realistic model we will encounter — it abstracts from the interest rate, from supply-side constraints, and from forward-looking behavior — but its simplicity is a virtue for the purpose of understanding the logic of the **multiplier**, one of the most important concepts in macroeconomics.

The multiplier is the proposition that an initial increase in autonomous spending (say, government expenditure) generates a larger total increase in income. The reason is the circular flow: the government's spending becomes someone's income; that person spends part of it, creating more income for others; those others spend part of their new income; and so on. The total increase in income is the sum of this infinite series of rounds. Understanding the multiplier — what determines its size, what limits it, and when it does not operate at all — is essential for evaluating fiscal policy.

---

## 8.1 The Basic Model

The Keynesian cross assumes: (i) the price level is fixed (we are in the short run with nominal rigidity); (ii) the economy has spare productive capacity, so that output is demand-determined; and (iii) there is a simple linear relationship between income and consumption. These assumptions make the model tractable at the cost of ignoring important feedback mechanisms that the IS–LM model restores.

The **consumption function** is:

$$C = a + b(Y - T), \quad a > 0,\; b \in (0,1),$$

where $a$ is **autonomous consumption** (the amount households would consume even if income were zero, financed by wealth drawdown), $b$ is the **marginal propensity to consume (MPC)**, and $Y - T$ is disposable income after lump-sum taxes $T$.

**Definition (Marginal Propensity to Consume).** The **MPC** is the fraction of each additional dollar of disposable income that households choose to spend on consumption, holding all other factors constant: $b = \partial C / \partial(Y-T)$. It satisfies $b \in (0,1)$: households spend more when income rises but do not spend the full increment (the remainder is saved). The Keynesian assumption $0 < b < 1$ is the foundation of the multiplier mechanism.

Investment is exogenous at $\bar{I}$, government spending at $\bar{G}$, and taxes at $\bar{T}$. Goods-market equilibrium requires that actual output equals desired expenditure:

$$Y = C + I + G = a + b(Y - \bar{T}) + \bar{I} + \bar{G}.$$

Solving for $Y$:

$$Y^* = \frac{1}{1-b}\bigl(a - b\bar{T} + \bar{I} + \bar{G}\bigr).$$

The factor $\kappa \equiv 1/(1-b)$ is the **Keynesian multiplier**. It amplifies any change in the autonomous components of expenditure. To understand why, note that the $n$-th round of spending generates additional income $b^n \Delta G$, where $\Delta G$ is the initial increase in government spending. Summing the geometric series:

$$\Delta Y^* = \Delta G \cdot (1 + b + b^2 + \cdots) = \Delta G \cdot \frac{1}{1-b} = \kappa\, \Delta G.$$

With $b = 0.75$, $\kappa = 4$: each dollar of government spending raises equilibrium income by four dollars. With $b = 0.5$, $\kappa = 2$.

---

## 8.2 Three Multipliers and Their Implications

### The Government Spending Multiplier

$$\kappa_G = \frac{\partial Y^*}{\partial \bar{G}} = \frac{1}{1-b}.$$

A one-dollar increase in government spending raises equilibrium income by $1/(1-b)$ dollars. This exceeds one because each round of induced consumption generates further income.

### The Tax Multiplier

$$\kappa_T = \frac{\partial Y^*}{\partial \bar{T}} = \frac{-b}{1-b}.$$

A one-dollar tax cut raises income by $b/(1-b)$ dollars. This is smaller in absolute value than $\kappa_G$ by a factor of $b$, because a tax cut delivers income to households who then spend fraction $b$ in the first round, whereas government spending delivers a full dollar of expenditure in the first round. The tax multiplier is always less effective than the spending multiplier at stimulus, a comparison with important policy implications.

### The Balanced Budget Multiplier (Haavelmo's Theorem)

Suppose government spending and taxes both increase by the same amount: $\Delta\bar{G} = \Delta\bar{T} = \Delta B$. The combined effect:

$$\Delta Y^* = \kappa_G \Delta B + \kappa_T \Delta B = \frac{1}{1-b}\Delta B + \frac{-b}{1-b}\Delta B = \frac{1-b}{1-b}\Delta B = \Delta B.$$

The balanced budget multiplier is exactly one — the income increase equals the spending increase — regardless of the MPC. This is **Haavelmo's theorem** (Haavelmo, 1945). The intuition: the government's spending directly raises income by $\kappa_G \Delta B$, while the tax increase reduces private spending by $|\kappa_T|\Delta B = b\kappa_G\Delta B$. The net effect is $(1-b)\kappa_G\Delta B = \Delta B$. Even a balanced fiscal expansion raises income, because the government spends every dollar it taxes while households would have saved fraction $1-b$.

---

## 8.3 Proportional Taxes and Automatic Stabilizers

The lump-sum tax of Section 8.1 is a simplification. In practice, taxes are proportional to income: $T = tY$, where $t \in (0,1)$ is the marginal tax rate. Substituting:

$$Y^* = \frac{a + \bar{I} + \bar{G}}{1 - b(1-t)}, \qquad \kappa_G^{prop} = \frac{1}{1 - b(1-t)}.$$

Since $b(1-t) < b$, the multiplier with a proportional tax is smaller than with a lump-sum tax. The intuition is that proportional taxation "leaks" some of each round of additional income to taxes rather than back into consumption, dampening the multiplier.

This leakage is exactly what makes proportional taxation a powerful **automatic stabilizer**.

**Definition (Automatic Stabilizer).** An **automatic stabilizer** is a fiscal mechanism that reduces the size of output fluctuations in response to shocks, without requiring any discretionary policy action by the government. Proportional income taxes are an automatic stabilizer: when income falls in a recession, tax revenues fall automatically, partially offsetting the income decline and supporting consumption. Conversely, in a boom, rising tax revenues automatically reduce disposable income, dampening the expansion. Automatic stabilizers reduce the effective multiplier but also reduce the volatility of income around its equilibrium level.

Other important automatic stabilizers include unemployment insurance (when workers lose jobs, UI payments replace part of lost income), and means-tested welfare benefits (which increase automatically when incomes fall). The strength of automatic stabilizers differs substantially across countries: economies with larger government sectors and more progressive tax systems have stronger automatic stabilization.

---

## 8.4 The Open Economy Multiplier

In an open economy, some of each round of induced consumption falls on imported goods rather than domestic output. With import function $M = m_0 + m_Y Y$, where $m_Y \in (0,1)$ is the **marginal propensity to import**, goods-market equilibrium becomes:

$$Y = a + b(Y - \bar{T}) + \bar{I} + \bar{G} + \bar{X} - m_0 - m_Y Y.$$

Solving:

$$\kappa_G^{open} = \frac{1}{1 - b + m_Y}.$$

Since $m_Y > 0$, the open-economy multiplier is smaller than the closed-economy multiplier. Part of the income stimulus "leaks" abroad through higher imports rather than circulating domestically. For a small open economy with a high propensity to import (say $m_Y = 0.3$) and a moderate MPC ($b = 0.75$), $\kappa_G^{open} = 1/(0.25 + 0.3) = 1.8$ compared to the closed-economy $\kappa_G = 4$. International leakages substantially reduce the effectiveness of fiscal stimulus in open economies.

---

## 8.5 The Limitations of the Keynesian Cross and What Comes Next

The Keynesian cross is a useful pedagogical tool, but it makes four assumptions that are too strong to be left unchallenged.

First, it assumes prices are fixed. But a demand stimulus eventually raises prices as the economy approaches potential output, which crowds out some of the initial stimulus. This motivates the AS–AD model (Chapter 7) and eventually the distinction between demand-pull and cost-push inflation (Chapter 30).

Second, it treats investment as exogenous. In reality, investment depends on the real interest rate, which changes when income changes. Adding this feedback is the contribution of the IS–LM model (Chapter 9), which generates the interest-rate crowding out that makes the fiscal multiplier smaller than $\kappa_G$ suggests.

Third, it models expectations as static. If households are forward-looking (Chapter 11) and expect a temporary fiscal stimulus to be reversed, they will smooth their consumption response across time, reducing the contemporaneous MPC and hence the multiplier. Under Ricardian equivalence (Chapter 22), forward-looking households fully anticipate the future tax increase implied by deficit spending and reduce consumption dollar-for-dollar today, driving the fiscal multiplier to zero.

Fourth, the behavioral consumption function is assumed rather than derived. Chapter 11 derives the consumption function from intertemporal optimization, showing that the MPC depends on whether income changes are perceived as permanent or transitory, on whether households face liquidity constraints, and on the degree of precautionary saving.

Each of these extensions enriches the basic framework, and together they explain why empirical estimates of fiscal multipliers vary from near zero (when all four extensions point toward small multipliers) to above two (at the effective lower bound on interest rates, when crowding out is absent and forward guidance prevents households from discounting fiscal stimulus). The range of estimates reflects genuine uncertainty about which model best describes the economy in any given situation.

---

*Next: Chapter 9 — The IS–LM Model*
