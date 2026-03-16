# Chapter 9 — The IS–LM Model: Interest Rates and Output in the Short Run

> *"Mr Keynes' system is a general equilibrium system. Quantities and prices are all determined together."*
> — J.R. Hicks, *Mr. Keynes and the 'Classics'*, 1937

---

The Keynesian cross assumes that investment is exogenous — fixed at $\bar{I}$ regardless of what happens in financial markets. This is obviously wrong: investment by firms depends critically on the real cost of borrowing, and the cost of borrowing depends on what is happening in money markets and the broader financial system. The IS–LM model, developed by Hicks (1937) as a formalization of Keynes's *General Theory*, corrects this by adding the money market to the goods market and solving for both output and the interest rate simultaneously. The result is a model of short-run macroeconomic equilibrium that remains central to policy analysis despite its limitations.

The name IS–LM refers to the two curves whose intersection defines equilibrium. The **IS curve** (Investment = Saving) is the locus of interest rate and output combinations consistent with goods-market equilibrium. The **LM curve** (Liquidity preference = Money supply) is the locus consistent with money-market equilibrium. Their intersection determines the unique combination of output and the interest rate at which both markets clear simultaneously.

---

## 9.1 The IS Curve: Goods-Market Equilibrium

The IS curve summarizes goods-market equilibrium as a relationship between the interest rate and output. To derive it, take the goods-market equilibrium condition:

$$Y = C(Y - T) + I(r) + G,$$

where $C(\cdot)$ is the consumption function with $C' \in (0,1)$ and $I(r)$ is the investment function with $I'(r) < 0$ — investment decreases when the real interest rate $r$ rises, because borrowing is more expensive and the hurdle rate for investment projects rises.

Totally differentiating with respect to $r$ while holding $G$ and $T$ fixed:

$$\mathrm{d}Y = C'\,\mathrm{d}Y + I'(r)\,\mathrm{d}r \implies \mathrm{d}Y(1-C') = I'(r)\,\mathrm{d}r.$$

The slope of the IS curve in $(Y, r)$ space:

$$\left.\frac{\mathrm{d}r}{\mathrm{d}Y}\right|_{IS} = \frac{1-C'}{I'(r)} < 0,$$

since $1-C' > 0$ and $I' < 0$. The IS curve is **downward sloping**: a higher real interest rate reduces investment; lower investment reduces output through the Keynesian multiplier. The steeper the IS curve (in absolute slope terms), the less sensitive is investment to the interest rate — either because firms' investment decisions are insensitive to borrowing costs, or because the multiplier is small.

**Shifts in the IS Curve.** Any change in autonomous expenditure shifts the IS curve horizontally. The horizontal shift from a change in $G$:

$$\left.\frac{\mathrm{d}Y}{\mathrm{d}G}\right|_{IS \text{ shift}} = \frac{1}{1-C'} = \kappa_G.$$

A fiscal expansion shifts IS to the right by $\kappa_G \Delta G$ at every interest rate. Conversely, a tax increase shifts IS to the left, and a decline in consumer confidence (a fall in autonomous consumption $a$) also shifts IS left.

---

## 9.2 The LM Curve: Money-Market Equilibrium

The LM curve describes equilibrium in the market for money. Money demand depends on two things: income (higher income means more transactions, requiring more money) and the nominal interest rate (which measures the opportunity cost of holding money rather than interest-bearing bonds). The real money demand function $L(Y, i)$ satisfies $L_Y > 0$ and $L_i < 0$.

**Definition (Money Demand).** Real money demand $L(Y,i)$ is the quantity of real money balances (nominal money holdings deflated by the price level) that households and firms wish to hold, as a function of real income $Y$ and the nominal interest rate $i$. The demand for money is increasing in income because higher income means more transactions requiring money as a medium of exchange, and decreasing in the nominal interest rate because a higher $i$ raises the cost of holding non-interest-bearing money rather than bonds.

Money-market equilibrium requires real money demand to equal the real money supply:

$$\frac{M}{P} = L(Y, i).$$

For a given real money supply $M/P$, the LM curve is the set of $(Y, i)$ pairs satisfying this condition. Totally differentiating:

$$0 = L_Y\,\mathrm{d}Y + L_i\,\mathrm{d}i \implies \left.\frac{\mathrm{d}i}{\mathrm{d}Y}\right|_{LM} = -\frac{L_Y}{L_i} > 0.$$

The LM curve is **upward sloping**: when income rises, money demand rises; to restore equilibrium at a fixed money supply, the interest rate must rise to reduce money demand back to its original level.

Two special cases of the LM curve deserve definition, because they represent the extremes of monetary policy effectiveness.

**Definition (Liquidity Trap).** The economy is in a **liquidity trap** when the nominal interest rate is at its lower bound (approximately zero) and money demand becomes perfectly elastic with respect to $i$ — that is, $L_i \to -\infty$. In this case the LM curve is horizontal, and any increase in the money supply is absorbed by money demand without affecting the interest rate or output. Monetary policy is completely ineffective in a liquidity trap.

**Definition (Classical Case).** The **classical case** occurs when money demand is completely insensitive to the interest rate: $L_i = 0$. In this case the LM curve is vertical at $Y = PL^{-1}(M/P)$: the level of income is entirely determined by the money supply, independently of fiscal policy. In the classical case, fiscal expansions raise the interest rate one-for-one and crowd out investment completely, so the fiscal multiplier is zero.

---

## 9.3 IS–LM Equilibrium and Policy Multipliers

Using the linear specifications $Y = \bar{A} - b_r\, r$ (IS, where $\bar{A}$ is autonomous expenditure and $b_r > 0$ measures investment interest sensitivity) and $M/P = kY - h\,i$ (LM, where $k > 0$ and $h > 0$), the equilibrium $(Y^*, i^*)$ satisfies both equations simultaneously. Using $i = r + \pi^e$ with fixed expected inflation:

$$Y^* = \frac{h\bar{A} + b_r(M/P) + b_r h\pi^e}{h + b_r k}.$$

### The Fiscal Multiplier and Crowding Out

An increase $\Delta G$ raises $\bar{A}$ by $\Delta G$. The IS–LM fiscal multiplier:

$$\frac{\Delta Y^*}{\Delta G} = \frac{h}{h + b_r k} < \frac{1}{1-C'} = \kappa_G.$$

This is strictly less than the Keynesian cross multiplier because of **crowding out**.

**Definition (Crowding Out).** **Crowding out** occurs when a fiscal expansion raises income, which raises money demand, which raises the interest rate, which reduces private investment. The increase in government spending is partially offset by the decline in private investment. The degree of crowding out depends on:

- **Investment interest sensitivity $b_r$**: if $b_r = 0$, there is no crowding out (IS is vertical; the multiplier equals $\kappa_G$).
- **Money demand interest sensitivity $h$**: if $h \to \infty$ (liquidity trap), a small interest rate rise generates unlimited money demand, so the interest rate effectively does not rise — no crowding out, and the multiplier approaches $\kappa_G$.

### The Monetary Multiplier

An increase in the real money supply $\Delta(M/P)$:

$$\frac{\Delta Y^*}{\Delta(M/P)} = \frac{b_r}{h + b_r k}.$$

Monetary policy works by shifting LM right (reducing the interest rate), which stimulates investment and raises output. The multiplier is zero in two cases: the liquidity trap ($h \to \infty$, the denominator is infinite) and the case of investment interest insensitivity ($b_r = 0$, the numerator is zero). In the liquidity trap, monetary stimulus fills the economy with money that no one lends, and nothing happens to investment or output.

---

## 9.4 The Mundell–Fleming Model: IS–LM in the Open Economy

In an open economy, a third equilibrium condition is required: the balance of payments must be in equilibrium. The **BP curve** is the locus of $(Y, i)$ consistent with balance-of-payments equilibrium:

$$NX(Y, Y^*, \bar{e}) + KA(i - i^*) = 0,$$

where $i^*$ is the world interest rate, $KA$ is the net capital inflow (increasing in $i - i^*$ as domestic assets become more attractive relative to foreign assets), $Y^*$ is foreign income, and $\bar{e}$ is the fixed exchange rate.

**Definition (Capital Mobility).** **Capital mobility** refers to the ease with which financial capital can flow across borders in response to interest rate differentials. Under **perfect capital mobility**, $KA' \to \infty$: any interest rate above $i^*$ immediately attracts unlimited foreign capital; any rate below $i^*$ immediately triggers capital outflows. Under perfect capital mobility, uncovered interest parity holds continuously: $i = i^*$, and the BP curve is horizontal.

The **Mundell–Fleming trilemma** is one of the most important propositions in international macroeconomics:

**Definition (The Trilemma).** A country cannot simultaneously maintain all three of: (i) a fixed exchange rate; (ii) perfect capital mobility; and (iii) independent monetary policy. It can achieve any two, but not all three at once.

The trilemma arises because: if the exchange rate is fixed and capital flows freely, the domestic interest rate must equal the world rate ($i = i^*$) at all times. If the central bank tries to set $i \neq i^*$, capital will flow in or out until reserves are exhausted or until the peg breaks. Hence the central bank cannot independently set both the exchange rate and the interest rate.

Under a fixed exchange rate with perfect capital mobility, fiscal policy is fully effective (the BP curve prevents crowding out by keeping $i = i^*$), but monetary policy is powerless (any attempt to set $i$ below $i^*$ loses reserves immediately). Under a flexible exchange rate with perfect capital mobility, monetary policy works through the exchange rate channel (a rate cut depreciates the currency, improving net exports), while fiscal policy crowds out net exports through exchange rate appreciation.

---

## 9.5 Critical Assessment of the IS–LM Model

The IS–LM model has been the backbone of macroeconomic policy analysis for nearly nine decades, but it has important limitations that motivate the more sophisticated models developed in subsequent chapters.

The model is **static**: it characterizes a contemporaneous equilibrium without specifying how the economy moves over time. The LM curve will shift when the price level changes (via $M/P$), but the model has no theory of price dynamics. The AS–AD framework of Chapter 7 provides partial dynamics; the full New Keynesian model of Chapters 9 and 23 provides a properly dynamic system.

The model uses **reduced-form behavioral equations** — the consumption and investment functions, the money demand function — rather than deriving behavior from intertemporal optimization. This means it is potentially subject to the Lucas critique: the parameters may shift when policy regimes change. Part III addresses this by deriving each component from first principles.

The model's money demand function proved empirically unstable in the 1970s and 1980s as financial innovation altered the velocity of money, rendering the LM curve difficult to identify. This is one reason modern central banks target the interest rate directly (the Taylor rule) rather than the money supply, making the IS curve alone the primary policy transmission mechanism.

Despite these limitations, IS–LM provides an indispensable conceptual framework: the distinction between crowding out and accommodation, the conditions for monetary versus fiscal policy effectiveness, and the logic of the trilemma all emerge naturally from it and survive in more sophisticated models.

---

*Next: Chapter 10 — The Phillips Curve*
