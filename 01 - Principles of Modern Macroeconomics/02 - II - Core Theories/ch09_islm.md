# Chapter 9 — The IS–LM Model: Interest Rates and Output in the Short Run

> *"Mr Keynes' system is a general equilibrium system. Quantities and prices are all determined together."*
> — J.R. Hicks, *Mr. Keynes and the 'Classics'*, 1937

---

The Keynesian cross assumes that investment is exogenous — fixed at $\bar{I}$ regardless of what happens in financial markets. This is obviously unsatisfactory: investment by firms depends critically on the real cost of borrowing, and the cost of borrowing depends on what is happening in money markets. The IS–LM model, developed by Hicks (1937) as a formalization of Keynes's *General Theory*, corrects this by adding the money market to the goods market and solving for both output and the interest rate simultaneously. The model produced the single most important analytical toolkit in macroeconomic policy analysis for nearly five decades, and its central insights — crowding out, the monetary transmission mechanism, the trilemma — survive in the more sophisticated models that followed.

The name IS–LM refers to the two curves whose intersection defines equilibrium. The **IS curve** is the locus of interest rate-output combinations consistent with goods-market equilibrium. The **LM curve** is the locus consistent with money-market equilibrium. Their intersection determines the unique pair $(Y^*, i^*)$ at which both markets clear simultaneously.

---

## 9.1 The IS Curve: Goods-Market Equilibrium

### Derivation

The IS curve summarizes goods-market equilibrium as a relationship between the interest rate and output. Begin from the equilibrium condition:

$$Y = C(Y - T) + I(r) + G,$$

where $C' \in (0,1)$ and $I'(r) < 0$: investment falls when the real cost of borrowing rises. Totally differentiating with respect to $r$ while holding $G$ and $T$ fixed:

$$\mathrm{d}Y = C'\,\mathrm{d}Y + I'(r)\,\mathrm{d}r \implies \mathrm{d}Y(1-C') = I'(r)\,\mathrm{d}r.$$

The slope of the IS curve in $(Y, r)$ space:

$$\left.\frac{\mathrm{d}r}{\mathrm{d}Y}\right|_{IS} = \frac{1-C'}{I'(r)} < 0.$$

The IS curve is **downward sloping**: higher interest rates reduce investment; lower investment reduces income through the Keynesian multiplier; lower income further reduces consumption — driving output down to maintain goods-market balance.

### The Slope and Its Economic Interpretation

The steepness of the IS curve has direct policy implications. The IS curve is steeper (in absolute value) when:
- **Investment is interest-inelastic** ($|I'|$ small): firms' investment plans barely respond to borrowing costs, so a large interest rate change is required to close any given output gap. Steep IS implies strong fiscal multipliers (interest rates do not rise much to crowd out spending) and weak monetary multipliers (rate cuts barely stimulate investment).
- **The MPC is small** ($C'$ small): each round of induced consumption is small, making the income-multiplying mechanism weak. A given interest rate reduction generates less total output increase.

Conversely, a flat IS curve characterizes an economy where investment is highly sensitive to interest rates and the MPC is high — conditions under which monetary policy is powerful and fiscal policy's interest-rate crowding out substantially offsets its demand effect.

### Shifts in the IS Curve

Any change in autonomous expenditure shifts the IS curve horizontally. A change in $G$ shifts IS by the Keynesian multiplier $\kappa_G = 1/(1-C')$ times the change in $G$, at every interest rate. Fiscal expansions shift IS right; fiscal contractions shift IS left. Changes in consumer or investor confidence (shifts in $a$ or $\bar{I}$), in foreign income (affecting exports in an open economy), or in taxes ($\Delta T$ shifts IS by $\kappa_T\Delta T$) all shift IS without changing its slope.

---

## 9.2 The LM Curve: Money-Market Equilibrium

### Derivation

The LM curve describes equilibrium in the market for real money balances. Money demand $L(Y, i)$ satisfies $L_Y > 0$ (more income means more transactions requiring money) and $L_i < 0$ (a higher nominal interest rate raises the opportunity cost of holding non-interest-bearing money, reducing demand). Money-market equilibrium:

$$\frac{M}{P} = L(Y, i).$$

**Definition (Money Demand).** Real money demand $L(Y,i)$ is the quantity of real money balances that households and firms wish to hold as a function of real income $Y$ and the nominal interest rate $i$. The income elasticity is positive (the transactions motive) and the interest elasticity is negative (the portfolio motive and opportunity cost). The Baumol-Tobin model [Ch. 14] provides microfoundations: $L^{BT} = \sqrt{bY/(2i)}$, predicting income elasticity $+1/2$ and interest elasticity $-1/2$.

Totally differentiating:

$$0 = L_Y\,\mathrm{d}Y + L_i\,\mathrm{d}i \implies \left.\frac{\mathrm{d}i}{\mathrm{d}Y}\right|_{LM} = -\frac{L_Y}{L_i} > 0.$$

The LM curve is **upward sloping**: when income rises, money demand rises; to restore equilibrium with a fixed money supply, the interest rate must rise to reduce money demand back to the supply level.

### Special Cases: Liquidity Trap and Classical Case

**Definition (Liquidity Trap).** The economy is in a **liquidity trap** when the nominal interest rate is at or near its lower bound (approximately zero) and money demand becomes perfectly elastic with respect to $i$: $L_i \to -\infty$. In this case the LM curve is horizontal: any increase in the money supply is absorbed by money demand without affecting the interest rate (since even a tiny rate decline generates unlimited demand for money), so monetary policy loses its power to stimulate investment and output. Keynes identified this as a theoretical possibility; Japan after 1998 and the global economy after 2008 demonstrated it as an empirical reality.

**Definition (Classical Case).** The **classical case** occurs when money demand is completely insensitive to the interest rate: $L_i = 0$. The LM curve is vertical at $Y = (M/P)/k$ (with $L = kY$): income is entirely pinned by the money supply, independently of fiscal policy. In this case, any fiscal expansion raises the interest rate one-for-one (crowding out investment completely) and the fiscal multiplier is zero. This is the original quantity-theory result in IS-LM clothing.

### Shifts in the LM Curve

The real money supply $M/P$ shifts the LM curve: an increase in $M$ (monetary expansion) shifts LM right (at any income level, the interest rate falls to re-equate money demand with the higher supply). A rise in $P$ reduces $M/P$ and shifts LM left (contractionary). Any shock to money demand that alters $L_Y$ or $L_i$ also shifts LM — a channel that proved empirically important when financial innovation in the 1980s destabilized money demand and rendered LM targeting unreliable.

---

## 9.3 IS–LM Equilibrium and Policy Multipliers

### Solving the System

Using linear IS and LM: $i = i_0 - (1/b_r)Y$ (IS, rearranged) and $i = (k/h)Y - M/(Ph)$ (LM), the equilibrium interest rate and output satisfy both simultaneously. Setting the two expressions equal and solving for $Y^*$:

$$Y^* = \frac{h\bar{A} + b_r M/P}{h + b_r k},$$

where $\bar{A}$ is autonomous expenditure, $b_r = |I'| > 0$ is investment interest sensitivity, $k = L_Y > 0$, and $h = -L_i > 0$.

### The Fiscal Multiplier and Crowding Out

An increase $\Delta G$ raises $\bar{A}$ by $\Delta G$. The IS–LM fiscal multiplier:

$$\mu_G^{IS-LM} = \frac{\partial Y^*}{\partial G} = \frac{h}{h + b_r k} < \frac{1}{1-C'} = \kappa_G.$$

The IS-LM multiplier is strictly less than the Keynesian cross multiplier because of **crowding out**.

**Definition (Crowding Out).** **Crowding out** occurs when a fiscal expansion raises income, which raises money demand, which raises the interest rate, which reduces private investment. The net increase in output is smaller than the direct fiscal injection because private investment has partially retreated. The degree of crowding out depends on:
- **Investment interest sensitivity $b_r$**: when $b_r \to 0$, there is no crowding out (investment doesn't respond to rates), and $\mu_G^{IS-LM} \to \kappa_G$.
- **Money demand interest sensitivity $h$**: when $h \to \infty$ (liquidity trap), the interest rate is pinned and money demand absorbs any excess, preventing crowding out, so $\mu_G^{IS-LM} \to \kappa_G$.
- **Money demand income sensitivity $k$**: when $k \to 0$, income growth barely raises money demand and hence barely raises the rate, so crowding out is minimal.

The fiscal multiplier is largest (approaching $\kappa_G$) when the LM curve is flat (liquidity trap or nearly so) and smallest when the LM curve is steep (classical case).

### The Monetary Multiplier

An increase in the real money supply $\Delta(M/P)$:

$$\mu_M = \frac{\partial Y^*}{\partial(M/P)} = \frac{b_r}{h + b_r k}.$$

Monetary easing shifts LM right, reducing the interest rate, which stimulates investment and raises output. The monetary multiplier is zero when $b_r = 0$ (investment insensitive to rates — IS is vertical) or when $h \to \infty$ (liquidity trap — LM is horizontal, rates cannot fall further). These two conditions exactly characterize the cases in which monetary policy is powerless: Keynesian "pushing on a string" (liquidity trap) and classical "investment determined by animal spirits" (vertical IS).

---

## 9.4 The Mundell–Fleming Model: IS–LM in the Open Economy

### Adding the Balance of Payments

In an open economy, a third equilibrium condition joins IS and LM: the **balance of payments** must be in equilibrium. The **BP curve** is the locus of $(Y, i)$ combinations consistent with a zero overall BOP balance:

$$NX(Y, Y^*, \bar{e}) + KA(i - i^*) = 0,$$

where $NX$ is the trade balance (decreasing in $Y$ through higher imports, increasing in $Y^*$ through higher exports, and dependent on the exchange rate $\bar{e}$), and $KA$ is the net capital inflow (increasing in $i - i^*$ as domestic bonds become more attractive).

**Definition (Capital Mobility).** **Capital mobility** refers to the ease with which financial capital flows across borders in response to interest differentials. Under **perfect capital mobility**, any interest rate above $i^*$ immediately attracts unlimited foreign capital (preserving uncovered interest parity $i = i^*$), making the BP curve horizontal. Under zero capital mobility, $KA = 0$ and the BP curve is vertical (only the trade balance matters for BOP equilibrium).

### The Trilemma

**Definition (The Impossible Trinity).** A country cannot simultaneously maintain: (i) a fixed exchange rate; (ii) perfect capital mobility; and (iii) independent monetary policy. At most two of the three are achievable simultaneously.

The logic: with a fixed exchange rate and perfect capital mobility, $i = i^*$ must hold at all times (any deviation triggers unlimited capital flows that exhaust reserves). The central bank therefore cannot independently set $i \neq i^*$ — it must accommodate whatever rate is required to maintain the peg. Monetary policy is entirely subordinated to the exchange rate target.

The trilemma implies three regime corners. **Bretton Woods** (1944–1971): fixed exchange rates plus capital controls (no capital mobility) — preserving monetary autonomy. **Gold standard**: fixed exchange rates plus free capital flows — sacrificing monetary autonomy. **Post-1973 floating**: monetary autonomy plus free capital flows — floating exchange rates. In practice, most economies operate at interior combinations: managed floats with partial capital controls, or currency union membership that sacrifices monetary policy but maintains fiscal flexibility.

### Policy Effectiveness Under Different Regimes

Under **fixed exchange rates with perfect capital mobility**: fiscal policy is fully effective ($i$ is pinned at $i^*$, preventing crowding out), but monetary policy is completely powerless (any monetary easing causes reserves outflows that reverse the expansion until $i$ returns to $i^*$). Under **floating exchange rates with perfect capital mobility**: monetary policy is effective through both the interest rate channel and the exchange rate channel (a rate cut depreciates the currency, boosting net exports — adding to the standard investment channel), while fiscal policy crowds out net exports through currency appreciation — leaving output approximately unchanged in the textbook "perfect mobility" case.

---

## 9.5 Empirical Evidence and the IS-LM Model's Performance

### What the Empirics Show

The IS-LM framework generates several testable predictions. Romer and Romer (1989) identify monetary policy shocks using narrative Fed documents; they find that monetary contractions raise the federal funds rate and persistently reduce output, with the peak output effect approximately 18 months after the shock — consistent with the IS-LM mechanism in which higher rates reduce investment, which then reduces income through multiplier dynamics.

Blanchard and Perotti (2002) use timing restrictions to identify fiscal shocks, finding spending multipliers of approximately 0.9 and tax multipliers of approximately 1.3 over several quarters — broadly consistent with IS-LM logic, though the tax multiplier exceeds the spending multiplier contrary to simple IS-LM predictions (possibly reflecting anticipatory effects and supply-side channels).

The empirical stability of money demand — the foundation of the LM curve — deteriorated sharply in the 1980s as financial innovation (money market mutual funds, NOW accounts, sweep accounts) made the velocity of M1 unstable. This is why modern central banks target the interest rate (the Taylor rule) rather than the money stock: the IS curve alone provides the transmission mechanism, and the LM curve is rendered implicit. The "New Keynesian" IS-LM removes the LM curve entirely, replacing it with a central bank reaction function.

---

## 9.6 Critical Assessment of the IS–LM Model

The IS–LM model has been the backbone of macroeconomic policy analysis for nearly nine decades — a remarkable longevity that reflects both its analytical clarity and its pedagogical power. Its limitations are well-known and have motivated successive generations of more sophisticated models.

The model is **static**: it characterizes a contemporaneous equilibrium without specifying the dynamics of adjustment. The LM curve shifts when the price level changes, but the model offers no theory of price dynamics — that requires the AS-AD framework of Chapter 7 and ultimately the NK three-equation system of Chapter 23.

The model uses **reduced-form behavioral equations** — the consumption function, the investment function, and money demand — rather than deriving behavior from intertemporal optimization. This exposes it to the Lucas critique: the estimated parameters may shift when the policy regime changes, making the model unreliable for policy evaluation across regime changes. Part III addresses this by deriving each behavioral relationship from microfoundations.

Money demand proved **empirically unstable** in the 1980s and 1990s, rendering LM targeting impossible and motivating the shift to interest rate rules. This is a critique not of the IS-LM framework's logic but of one of its empirical components.

Despite these limitations, IS-LM provides genuinely indispensable conceptual infrastructure: the distinction between crowding out and monetary accommodation, the conditions for fiscal versus monetary policy effectiveness, the logic of the trilemma, and the qualitative analysis of demand shocks all emerge naturally and survive in more sophisticated models. The New Keynesian model, properly understood, is IS-LM with a microfounded IS curve, a microfounded Phillips curve in place of the static price level, and a Taylor rule in place of the LM curve.

---

## Chapter Summary

- The **IS curve** is downward-sloping: higher interest rates reduce investment, which reduces output through the multiplier. Its slope depends on investment interest sensitivity $|I'|$ and the MPC $C'$; it shifts by $\kappa_G \Delta G$ for a fiscal change.

- The **LM curve** is upward-sloping: higher income raises money demand, requiring a higher interest rate to maintain equilibrium with a fixed money supply. Special cases: the **liquidity trap** (LM horizontal, monetary policy powerless) and the **classical case** (LM vertical, fiscal multiplier zero).

- IS-LM equilibrium $Y^* = (h\bar{A} + b_r M/P)/(h + b_r k)$ implies a fiscal multiplier $h/(h+b_r k) < \kappa_G$ due to **crowding out** (fiscal expansion raises rates, partially offsetting private investment), and a monetary multiplier $b_r/(h+b_r k)$ equal to zero in the liquidity trap or with interest-inelastic investment.

- The **Mundell-Fleming** extension adds the BP curve for open-economy analysis. The **impossible trinity** states that a country cannot simultaneously fix the exchange rate, allow free capital flows, and maintain monetary independence — choosing two of the three characterizes all exchange rate regimes.

- Empirical evidence (Romer-Romer monetary shocks; Blanchard-Perotti fiscal shocks) broadly supports IS-LM's qualitative predictions. Money demand instability since the 1980s has displaced explicit LM targeting; modern central banks use interest rate rules (Taylor rule), making the IS curve alone the operative transmission mechanism.

---

*Next: Chapter 10 — The Phillips Curve*
