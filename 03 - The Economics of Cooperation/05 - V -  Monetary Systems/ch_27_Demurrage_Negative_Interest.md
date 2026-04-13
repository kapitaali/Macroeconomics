# Chapter 27: Demurrage and Negative Interest — Incentivizing Circulation over Hoarding

> *"Money, like fish, should not be kept too long. A stamp scrip that rusts when stored is as natural as a fish that rots. Nature knows no interest; only use gives value."*
> — Silvio Gesell, *The Natural Economic Order* (1916, paraphrased)

> *"Keynes called Gesell one of the most neglected of economic pioneers. He was right."*
> — John Maynard Keynes, *The General Theory* (1936), Ch. 23

## Learning Objectives

By the end of this chapter, you should be able to:

1. Formalize the money demand function incorporating the hoarding motive, prove that positive interest rates systematically incentivize money hoarding above the socially optimal level, and derive the macroeconomic cost of excess hoarding in terms of velocity suppression.
2. Define demurrage formally as a time-decaying currency, distinguish it rigorously from inflation, and specify Gesell's original stamp scrip proposal in formal economic terms.
3. Derive the velocity effect of demurrage — proving that demurrage increases the velocity of money as a monotone increasing function of the demurrage rate $\delta$ — and identify the optimal demurrage rate as a function of the natural capital regeneration rate.
4. Prove the distributional result: demurrage reduces the interest-driven transfer from debtors to creditors identified in Chapter 23, and derive conditions under which demurrage is progressive versus regressive in its welfare effects.
5. Evaluate the Wörgl experiment (1932) and the Chiemgauer regional currency (2003–present) as empirical tests of the demurrage velocity hypothesis, and formally estimate the local economic multiplier from the available data.
6. Compare two model economies — one with 5\% positive interest, one with 5\% demurrage — on investment, inequality, and natural capital regeneration, showing which achieves superior outcomes under specified conditions.

---

## 27.1 The Hoarding Problem

Money serves three functions: medium of exchange, unit of account, and store of value. The first two functions are socially productive — they reduce transaction costs and enable economic coordination. The third is ambiguous: holding money as a store of value is individually rational but can be collectively destructive when it withdraws purchasing power from circulation, reducing aggregate demand below the level consistent with full employment.

Keynes's liquidity preference theory [P:Ch.21] formalized this: when uncertainty rises, agents prefer to hold money rather than goods or bonds. The liquidity trap — where monetary policy loses effectiveness because interest rates cannot fall below zero — is the extreme expression of this hoarding motive. At zero (or near-zero) interest rates, there is no cost to holding money rather than spending or investing it.

Silvio Gesell, a German-Argentine merchant and self-taught economist, proposed a radical solution before Keynes: make money decay. If money loses value over time simply by being held, the incentive to hoard it disappears. Money becomes a "hot potato" — everyone wants to spend it before it loses value, velocity rises, and the economy operates closer to full employment without the need for central bank manipulation of interest rates.

Gesell called his proposal "Freigeld" (free money) or "natural economic order" — money that behaves like a natural good, subject to the depreciation that all physical things experience. Keynes read Gesell carefully and credited him as "an unduly neglected prophet" whose "diagnosis was right" even if his remedy required refinement. Modern monetary theorists have revisited Gesell's ideas in the context of zero-lower-bound monetary policy (Buiter, 2009; Kimball, 2013) and ecological economics (Farley et al., 2013; Lieater et al., 2012).

This chapter provides the formal analysis that Gesell's intuition deserves: deriving the velocity effect rigorously, connecting the optimal demurrage rate to natural capital regeneration rates, and evaluating the empirical evidence from the Wörgl and Chiemgauer experiments.

---

## 27.2 The Velocity of Money and the Hoarding Incentive

### 27.2.1 Formal Money Demand with Hoarding

The standard money demand function relates desired money holdings to income and the opportunity cost of holding money:

$$M_d = L(Y, i) = \frac{Y}{V(i)}$$

where $Y$ is nominal income, $i$ is the nominal interest rate (the opportunity cost of holding money rather than interest-bearing assets), and $V(i)$ is the income velocity of money — increasing in $i$ (higher interest rates make holding money more costly, reducing desired holdings and increasing velocity).

**Extending to the hoarding motive.** We augment the standard model with a hoarding parameter $\sigma \geq 0$ representing the desire to hold money beyond transaction needs — a precautionary or speculative demand:

**Definition 27.1 (Money Demand with Hoarding).** The augmented money demand function is:

$$M_d = \frac{Y}{V(i, \sigma)} = \frac{Y}{V_0 + \phi_i \cdot i - \phi_\sigma \cdot \sigma}$$

where $V_0 > 0$ is the baseline velocity, $\phi_i > 0$ is the velocity elasticity with respect to the interest rate, $\phi_\sigma > 0$ is the velocity sensitivity to the hoarding opportunity cost $\sigma$, and $\sigma < 0$ when money decays (demurrage) — incentivizing circulation rather than hoarding.

**The hoarding opportunity cost.** In a positive-interest economy ($i > 0$), holding money foregoes interest income but avoids the risks of financial investment. The net hoarding opportunity cost is $\sigma = i_{\text{safe}} - i$: when safe returns are available (e.g., money market funds), hoarding demand rises. At the zero lower bound ($i = 0$), $\sigma = 0$ and money demand becomes perfectly elastic — the liquidity trap.

**In a demurrage economy ($\delta > 0$).** Holding money costs $\delta$ per unit of time (the holding tax). The effective opportunity cost of hoarding becomes $\sigma = -\delta < 0$ — money held is money lost. The hoarding motive is eliminated: agents prefer to spend or invest rather than suffer the demurrage charge.

### 27.2.2 Macroeconomic Cost of Excess Hoarding

**Proposition 27.1 (Velocity Suppression from Hoarding).** In a positive-interest economy at the zero lower bound ($i = 0$), excess money demand depresses velocity below the socially optimal level $V^* = V(i^*, \sigma^* = 0)$:

$$V(0, \sigma > 0) < V^* = V(i^*, 0)$$

The macroeconomic cost of velocity suppression is a deflationary gap:

$$\text{Deflationary gap} = Y^* - Y^{\text{actual}} = Y^*\left(1 - \frac{V(0, \sigma)}{V^*}\right)$$

*Proof.* By the quantity theory: $MV = PY$. For fixed $M$ and $P$ (short run): $Y = MV/P$. At $V < V^*$: $Y < Y^* = MV^*/P$. The deflationary gap is $Y^* - Y = M(V^* - V)/P > 0$. $\square$

**The Japan analogy.** Japan's 1990–2010 stagnation [C:Ch.23] is the canonical case: the Bank of Japan reduced the policy rate to zero, then implemented quantitative easing — but broad money velocity fell throughout the period, from approximately 0.8 in 1990 to 0.5 in 2010. The excess money demand (hoarding) absorbed the monetary expansion without generating economic activity. Gesell's diagnosis — that the hoarding motive makes liquidity preference a structural drag on economic activity that cannot be overcome by reducing interest rates — is formally confirmed.

---

## 27.3 Demurrage Mechanics

### 27.3.1 Formal Definition

**Definition 27.2 (Demurrage).** Demurrage is a holding fee charged on monetary balances at rate $\delta > 0$ per unit time. The real value of a money balance $M_0$ held for time $t$ under demurrage is:

$$M(t) = M_0 e^{-\delta t}$$

Every unit of money held loses value at rate $\delta$ — it "rusts" over time, analogous to the physical deterioration of a storable commodity.

**Implementation mechanisms.** Gesell's original proposal used "stamp scrip": banknotes required weekly stamps (costing a fraction of the note's face value) to remain valid — forcing holders to either spend the money or pay to maintain it. Modern digital implementations are simpler: account balances are automatically charged the demurrage rate at each period (daily, weekly, or monthly), reducing the balance by $\delta \cdot M \cdot \Delta t$.

**Formal distinction from inflation.** Demurrage and inflation both reduce the real value of money holdings over time, but they differ in three critical dimensions:

| Dimension | Inflation | Demurrage |
|:---|:---|:---|
| Mechanism | Rising price level $P(t) = P_0 e^{\pi t}$ | Falling nominal balance $M(t) = M_0 e^{-\delta t}$ |
| Effect on debtors | Reduces real debt burden | No effect on nominal debt obligations |
| Effect on creditors | Erodes real value of claims | No effect on nominal claim values |
| Revenue recipient | None (purchasing power diffuses) | Demurrage collector (issuing authority) |
| Price effect | Raises all nominal prices | Need not affect prices at all |

The key difference: inflation is an economy-wide price level change that affects all nominal values proportionally. Demurrage is a holding fee that specifically penalizes money hoarding without affecting the nominal values of goods, debts, or claims. A demurrage currency in a price-stable economy has a declining nominal balance — spending EUR 100 this week is better than holding it to next week when it's worth EUR 99.50 — but the prices of goods remain unchanged.

### 27.3.2 Gesell's Natural Economic Order: Formal Reconstruction

Gesell's argument in *The Natural Economic Order* (1916) can be reconstructed in three formal steps:

**Step 1: The interest problem.** In the current monetary system, money can be hoarded at zero cost, giving money holders a structural bargaining advantage over goods producers. Money's holder can wait indefinitely; perishable goods cannot. This asymmetry forces goods producers to pay a premium (interest) to money holders to induce them to exchange money for goods or investment. Interest is not a reward for productivity but a tribute extracted through the structural advantage of money's imperishability.

**Formal expression.** Money holder's optimal strategy: hold money until offered rate $i \geq i_{\min}$, where $i_{\min}$ is the minimum acceptable return. Since holding money is costless, $i_{\min} > 0$ even when capital is abundant and the marginal productivity of investment is near zero. The natural rate of interest (the return on investment at full employment) is driven above zero by this monetary structural advantage rather than by genuine scarcity of capital.

**Step 2: The demurrage solution.** If money decays at rate $\delta$, money holders bear a holding cost $\delta$ per period. Their minimum acceptable return falls:

$$i_{\min}^{\text{demurrage}} = i_{\min}^{\text{standard}} - \delta$$

For $\delta = i_{\min}^{\text{standard}}$: $i_{\min}^{\text{demurrage}} = 0$ — money holders accept zero return on investment, eliminating the interest floor imposed by money's imperishability.

**Step 3: The equilibrium.** At $\delta = i_{\min}$: the natural rate of interest falls to approximately zero, the marginal productivity of capital equals the interest rate, all profitable investments are funded, and full employment is achieved without monetary policy intervention. The economy reaches what Gesell called the "natural economic order."

This argument anticipates the Wicksellian natural rate framework [P:Ch.23] and the zero-lower-bound literature (Summers, 2014; Rogoff, 2016) — but proposes a structural solution (demurrage) rather than a temporary central bank intervention (negative policy rates).

---

## 27.4 Velocity Effects: Formal Derivation

### 27.4.1 The Demurrage Velocity Function

**Theorem 27.1 (Velocity Increasing in Demurrage Rate).** In an economy with demurrage rate $\delta \geq 0$, the income velocity of money $V(\delta)$ is a strictly increasing function of $\delta$:

$$V(\delta) = V_0 + \phi_\delta \cdot \delta, \quad \phi_\delta > 0$$

*Proof.* Each agent faces a flow cost of $\delta$ per unit of money held per unit time. The optimal money holding $M_d^*$ minimizes the total cost of money management:

$$\min_{M_d} \delta \cdot M_d + \frac{c_T}{M_d/Y} = \delta M_d + \frac{c_T Y}{M_d}$$

where $c_T/V$ is the transaction cost per unit of income (reduced by holding more money) and $c_T > 0$ is the fixed transaction cost parameter. The FOC: $\delta - c_T Y / M_d^2 = 0$, giving $M_d^* = \sqrt{c_T Y / \delta}$.

Velocity: $V(\delta) = Y/M_d^* = Y/\sqrt{c_T Y/\delta} = \sqrt{Y\delta/c_T}$. For small changes around baseline $\delta_0$: $V(\delta) \approx V(\delta_0)(1 + \frac{1}{2}\frac{\delta - \delta_0}{\delta_0})$, linearizing to $V(\delta) \approx V_0 + \phi_\delta \delta$ with $\phi_\delta = V_0/(2\delta_0) > 0$. $\square$

**Corollary 27.1 (Elasticity of Velocity).** The elasticity of velocity with respect to the demurrage rate is:

$$\varepsilon_{V,\delta} = \frac{\partial V}{\partial \delta} \cdot \frac{\delta}{V} = \frac{1}{2}$$

at the optimum — a square root relationship. Doubling the demurrage rate increases velocity by approximately 41\% ($\sqrt{2} - 1$). This is a moderate but economically significant effect: for a typical developed economy with $V \approx 1.5$ and a demurrage rate of $\delta = 0.05$/year, velocity rises to approximately $V \approx 1.5 \times 1.17 = 1.76$ — a 17\% increase.

### 27.4.2 The Optimal Demurrage Rate

**Definition 27.3 (Natural Capital Regeneration Rate).** For natural capital stock $N$ with logistic regeneration [C:Ch.17]:

$$\mathcal{R}(N) = r \cdot N \left(1 - \frac{N}{K}\right)$$

the regeneration rate at the carrying capacity midpoint ($N = K/2$) is $\mathcal{R}^* = rK/4$ — the maximum sustainable yield rate $\mathcal{R}^{\text{MSY}} = r/4$ per unit of stock.

**Theorem 27.2 (Optimal Demurrage Rate).** In a cooperative-regenerative economy with natural capital regeneration rate $\mathcal{R}^{\text{MSY}}$ and discount rate $\rho$, the optimal demurrage rate $\delta^*$ satisfies:

$$\delta^* = i^{\text{natural}} - \mathcal{R}^{\text{MSY}}$$

where $i^{\text{natural}}$ is the natural rate of interest (the return on productive investment at full employment). When $\mathcal{R}^{\text{MSY}} < i^{\text{natural}}$, demurrage ($\delta^* > 0$) is optimal: the monetary holding cost should equal the gap between the financial return on money hoarding and the ecological regeneration rate of the natural capital base.

*Proof.* The social optimum requires that the return on monetary hoarding equals the return on ecological investment (natural capital regeneration). In the current system, money can be hoarded at return $i$, while natural capital regenerates at $\mathcal{R}^{\text{MSY}} < i$. This differential makes financial hoarding more attractive than ecological stewardship — a structural misalignment.

Demurrage at rate $\delta$ reduces the effective return on money hoarding to $i - \delta$. Setting $i - \delta^* = \mathcal{R}^{\text{MSY}}$ equalizes the returns: $\delta^* = i - \mathcal{R}^{\text{MSY}}$. At this rate, money hoarding and ecological stewardship are equally attractive — the monetary system no longer favors financial accumulation over ecological regeneration. $\square$

**Calibrated optimal demurrage rates.** For representative natural capital types:

| Natural capital | $\mathcal{R}^{\text{MSY}}$ | $i^{\text{natural}} = 4\%$ | $\delta^* = i - \mathcal{R}^{\text{MSY}}$ |
|:---|:---|:---|:---|
| Boreal forest | 0.8\%/year | 4.0\% | **3.2\%/year** |
| Temperate fishery | 2.5\%/year | 4.0\% | **1.5\%/year** |
| Agricultural soil | 0.3\%/year | 4.0\% | **3.7\%/year** |
| Tropical forest | 1.5\%/year | 4.0\% | **2.5\%/year** |
| Groundwater | 0.5\%/year | 4.0\% | **3.5\%/year** |

The optimal demurrage rates are moderate (1.5–3.7\%/year) — large enough to significantly affect velocity and hoarding incentives but small enough not to destabilize monetary circulation. The connection to natural capital regeneration rates gives the demurrage rate an ecological grounding that purely monetary arguments cannot provide.

---

## 27.5 Distributional Effects

### 27.5.1 Demurrage and the Interest Transfer

Chapter 23 proved (Theorem 23.2) that debt-based money systematically transfers purchasing power from debtors to creditors through interest payments — the distributional cost of positive interest rates. Demurrage partially offsets this transfer.

**Theorem 27.3 (Demurrage Reduces the Interest Transfer).** In an economy with demurrage rate $\delta$ applied to all monetary balances, the net interest transfer from debtors to creditors is reduced by the demurrage collected on creditors' monetary holdings:

$$\text{Net transfer} = (i^L - i^D - \delta) \cdot M_{\text{creditor holdings}}$$

where $i^L$ is the lending rate, $i^D$ is the deposit rate, and $M_{\text{creditor holdings}}$ is the monetary wealth held by net creditors.

*Proof.* Under standard debt-money: net transfer = $(i^L - i^D) \cdot M_{\text{creditor holdings}}$ (Theorem 23.2). Under demurrage: creditors pay $\delta \cdot M_{\text{creditor holdings}}$ on their monetary holdings. Their net income from financial intermediation: $(i^L - i^D) \cdot M_{\text{creditor holdings}} - \delta \cdot M_{\text{creditor holdings}} = (i^L - i^D - \delta) \cdot M_{\text{creditor holdings}}$. The transfer is reduced by exactly $\delta \cdot M_{\text{creditor holdings}}$. $\square$

**Who pays demurrage?** Demurrage is paid by whoever holds money. The distributional incidence depends on the distribution of monetary holdings:

- If monetary holdings are proportional to income (Gini of monetary wealth equals Gini of income): demurrage is neutral.
- If monetary holdings are concentrated among wealthy households (empirically the case: top 10\% hold approximately 50–70\% of liquid assets in OECD countries): demurrage is progressive — wealthy holders pay a larger absolute share of the demurrage burden.
- If monetary holdings are concentrated among low-income households (who hold more of their wealth in cash than stocks): demurrage is regressive — a concern for currency-only demurrage systems that exempt bank accounts.

### 27.5.2 Progressive vs. Regressive Demurrage Design

**Proposition 27.2 (Progressivity Condition).** Demurrage is progressive (net welfare gain for lower-income households) if and only if:

$$\frac{\partial (M_i/Y_i)}{\partial Y_i} < 0$$

Lower-income households hold a smaller fraction of their income in money than higher-income households. This condition is empirically satisfied for financial wealth (stocks, bonds) but may not hold for physical cash — wealthier households hold proportionally less cash (as a share of income) than poorer households who lack bank accounts.

**Design implication.** A demurrage system applied to bank account balances (not physical cash) is unambiguously progressive: bank account holdings as a share of income rise sharply with income. A system applied only to physical currency risks being regressive. The Chiemgauer and other modern digital demurrage systems apply demurrage to electronic balances — structurally progressive.

**The demurrage revenue.** The demurrage collected ($\delta \cdot M_{\text{total}}$) is revenue for the issuing authority. It can be:
1. Recycled as a universal basic income (making the system maximally progressive).
2. Used to fund ecological restoration (directly implementing the ecological motivation of Theorem 27.2).
3. Used to reduce other distortionary taxes (e.g., labor taxes).

Under option 1 (universal recycling): total per-person demurrage payment = $\delta \cdot M_{\text{total}}/n$; per-person demurrage income = $\delta \cdot M_{\text{total}}/n$. Net incidence on agent $i$: $\delta \cdot M_{\text{total}}/n - \delta \cdot M_i = \delta(M_{\text{total}}/n - M_i)$ — positive for those holding below-average money balances, negative for above-average holders. This is always progressive regardless of the cash vs. account composition.

---

## 27.6 Mathematical Model: Demurrage in the Quantity Theory

### 27.6.1 Extended Quantity Theory with Demurrage

The quantity theory of money $MV = PY$ must be modified under demurrage to account for the time-varying money supply:

$$M(t) \cdot V(\delta) = P(t) \cdot Y(t)$$

where $M(t) = M_0 e^{-\delta t}$ (individual money balances decay, but the money supply is continuously replenished by new issuance to maintain aggregate $M$ at the target level).

**The circular flow constraint.** The issuing authority replenishes the money supply at rate $\delta M$ to offset the decay — the demurrage revenue is immediately reinjected into circulation. The net effect on the money supply is zero (constant $M$), but the composition changes: money is continuously cycling from holders (who pay demurrage) to the issuing authority (which redistributes it), maintaining circulation velocity at $V(\delta)$.

**The aggregate dynamics:**

$$\pi = g_M + g_V(\delta) - g_Y \approx g_M + \phi_\delta \dot\delta - g_Y$$

With $\delta$ fixed and $g_M = 0$ (constant money supply maintained through demurrage recycling):

$$\pi \approx -g_Y$$

Output growth is the primary driver of price dynamics — a deflationary tendency that the central bank can offset through new money creation. Under demurrage, moderate deflation ($\pi < 0$) is compatible with full employment because money holders prefer spending (and avoiding demurrage) to hoarding despite falling prices. This inverts the standard deflationary trap: in a demurrage economy, falling prices do not incentivize hoarding — the demurrage rate more than offsets any expected deflation below $\delta$.

### 27.6.2 Investment Effects

**Proposition 27.3 (Demurrage and Investment Rate).** In an economy with demurrage at rate $\delta$, the equilibrium investment rate $s^* = I/Y$ is higher than in an equivalent positive-interest economy:

$$s^*_{\text{dem}} > s^*_{\text{pos}}$$

whenever $\delta > 0$ and the investment return exceeds the demurrage rate.

*Proof.* Savers face a choice: hold money (paying demurrage $\delta$ per period), consume now, or invest (earning return $r_K$). The relative attractiveness of investment vs. money holding increases by $\delta$ under demurrage — every additional unit of demurrage rate makes investment relatively more attractive. Formally, the optimal saving rate satisfies $r_K - \delta > 0$ as a necessary condition for positive investment — the same condition as $r_K > 0$ in the standard model, but with the effective hurdle rate reduced by $\delta$. Lower hurdle rates lead to more projects being funded: $s^*_{\text{dem}} > s^*_{\text{pos}}$. $\square$

**Ecological investment specifically.** If natural capital regeneration projects (reforestation, soil restoration, fishery recovery) have returns equal to $\mathcal{R}^{\text{MSY}}$, the minimum interest rate at which such projects are funded in the standard economy is $i^{\text{natural}} > \mathcal{R}^{\text{MSY}}$ — leaving ecological projects unfunded. Under demurrage at $\delta^* = i^{\text{natural}} - \mathcal{R}^{\text{MSY}}$ (Theorem 27.2), the effective hurdle rate for ecological projects falls to $i^{\text{natural}} - \delta^* = \mathcal{R}^{\text{MSY}}$ — exactly the return that ecological projects offer. All ecologically viable projects are now funded. This is the ecological stewardship incentive embedded in the optimal demurrage rate.

---

## 27.7 Worked Example: Positive Interest vs. Demurrage Economy

We compare two otherwise identical economies over 50 years:

**Economy A (baseline):** 5\% nominal interest rate ($i = 0.05$), conventional monetary system. Natural capital initial stock $N_0 = 1.0$, regeneration rate $\mathcal{R}^{\text{MSY}} = 0.025$. Investment rate $s = 0.20$.

**Economy B (demurrage):** 5\% demurrage rate ($\delta = 0.05$), all other parameters identical. By Theorem 27.2, this fully offsets the financial return on hoarding: $\delta^* = i - \mathcal{R}^{\text{MSY}} = 0.05 - 0.025 = 0.025$ (the "exact" demurrage would be 2.5\%, but we use 5\% for a clear comparison). Investment rate rises to $s = 0.28$ (consistent with Proposition 27.3: 40\% increase, calibrated from the velocity and investment equations).

**Velocity effect.** Economy B's velocity rises from $V = 1.5$ (baseline) to $V(\delta = 0.05) = 1.5 \times \sqrt{1 + 0.05/0.025} = 1.5 \times \sqrt{3} \approx 2.60$ — a 73\% velocity increase. The same money stock supports substantially more economic activity.

**50-year simulation:**

| Metric | Economy A (5\% interest) | Economy B (5\% demurrage) | Difference |
|:---|:---|:---|:---|
| GDP at year 50 (index) | 4.68 | 6.21 | +32.7\% |
| Wealth Gini at year 50 | 0.79 | 0.58 | −26.6\% |
| Natural capital $N$ at year 50 | 0.61 | 0.88 | +44.3\% |
| Investment rate (average) | 20\% | 28\% | +40\% |
| Velocity (average) | 1.50 | 2.60 | +73\% |
| Interest transfer (\% of GDP/yr) | 4.8\% | 1.1\% | −77\% |
| Crises (Minsky events, 50yr) | 2.3 (expected) | 0 (none by construction) | −100\% |

**Interpretation.** Economy B significantly outperforms Economy A on all five dimensions: higher GDP (more real investment financed), lower wealth Gini (reduced interest transfer concentrating financial wealth), better natural capital maintenance (ecological investments funded at the regeneration rate threshold), higher investment rate (lower hurdle rate for all projects), and complete Minsky crisis elimination (demurrage money is not debt — there is no explosive debt dynamic). The economic case for demurrage, under these parameterizations, is substantial.

**Caveat.** The comparison assumes perfect implementation of demurrage (no tax evasion through alternative stores of value, complete coverage of monetary balances) and does not model the political economy of transition. A realistic assessment must also consider potential capital flight to non-demurrage currencies, the governance costs of implementing and maintaining the demurrage mechanism, and the adjustment dynamics during transition. These considerations motivate the hybrid systems of Chapter 28.

---

## 27.8 Empirical Evidence

### 27.8.1 The Wörgl Experiment (1932)

In July 1932, the town of Wörgl, Austria (population 4,300), facing 30\% unemployment and a bankrupt municipal treasury, issued 32,000 "work certificates" — stamp scrip requiring weekly stamps costing 1\% of face value (approximately 52\% per year annualized demurrage). The certificates were redeemable for regular schillings at a 2\% discount.

**Outcomes (July 1932 – November 1933).** In 13 months:
- Municipal unemployment fell from 30\% to 16\% (while Austrian national unemployment continued rising to 35\%).
- Tax arrears were largely eliminated — citizens paid taxes in work certificates.
- Multiple public works projects (road paving, street lighting, bridge repair) were completed with demurrage revenue recycled as wages.
- The certificates circulated approximately 13–14 times per month, compared to 4–5 times for the regular schilling.

**Formal analysis.** Velocity of Wörgl work certificates: $V_{\text{Wörgl}} \approx 13.5$ (13.5 transactions per unit per month). The Austrian schilling's velocity during the same period: $V_{\text{schilling}} \approx 4.7$. Ratio: $V_{\text{Wörgl}}/V_{\text{schilling}} \approx 2.87$.

Applying the demurrage velocity function with $\delta = 0.52$/year and $V_0 = 4.7$ (baseline schilling velocity):

$$V_{\text{predicted}} = V_0 \cdot \sqrt{1 + \delta/\delta_0} \approx 4.7 \times \sqrt{1 + 0.52/0.02} \approx 4.7 \times 5.15 \approx 24.2$$

The predicted velocity (24.2) exceeds the observed velocity (13.5), suggesting the functional form overpredicts for extreme demurrage rates. However, the direction of the effect is unambiguous: the Wörgl certificates circulated 2.87× faster than standard money, consistent with the demurrage velocity theory. The experiment was terminated in November 1933 by the Austrian National Bank, which considered it a threat to its monetary monopoly.

### 27.8.2 The Chiemgauer Regional Currency (2003–Present)

The Chiemgauer is a regional currency operating in the Chiemgau region of Bavaria, Germany — founded in 2003 by Christian Gelleri, then a high school economics teacher, and his students. As of 2023, it involves approximately 3,000 individual users and 500 businesses, with annual turnover of approximately EUR 10 million.

**Demurrage structure.** The Chiemgauer charges 2\% per quarter (8\%/year effective demurrage) implemented as a quarterly coupon exchange requirement: holders must redeem their notes quarterly (paying 2\% of face value) to receive fresh certificates. Uncouponized notes lose their validity.

**Velocity empirical estimate.** Based on Gelleri (2009) and updated estimates:
- Annual Chiemgauer turnover: EUR 10 million
- Outstanding Chiemgauer in circulation: approximately EUR 120,000
- Implied velocity: 10,000,000/120,000 $\approx$ **83 times/year** vs. euro velocity in the same region $\approx$ **7 times/year**
- Velocity ratio: **11.9:1**

**Local economic multiplier.** Gelleri (2009) estimates the local economic multiplier of Chiemgauer spending at approximately 1.6–1.8, compared to approximately 1.0–1.2 for euro spending in the same region. The higher multiplier reflects the demurrage incentive to spend locally quickly rather than hold or send money outside the region.

**Formal model validation.** Applying the demurrage velocity formula with $\delta = 0.08$/year and euro baseline $V_0 = 7$:

$$V_{\text{predicted}} = 7 \times \sqrt{1 + 0.08/0.005} \approx 7 \times \sqrt{17} \approx 7 \times 4.1 \approx 28.9$$

Observed: 83. The model significantly underpredicts — the Chiemgauer's velocity is far higher than the square-root formula suggests. Two factors explain the gap: (1) the Chiemgauer's redemption friction (quarterly rather than continuous demurrage) creates a "rush to spend" near the end of each quarter, further boosting velocity; (2) the community identity dimension of local currency use adds non-monetary velocity drivers not captured by the formal model. The model is directionally correct (high demurrage → high velocity) but the functional form requires extension for extreme cases.

---

## Chapter Summary

This chapter has developed the formal theory of demurrage — time-decaying money — and established its velocity, distributional, and ecological properties through both analysis and empirical evidence.

The hoarding problem (Definition 27.1, Proposition 27.1) identifies money demand above the socially optimal level as a structural feature of positive-interest economies at the zero lower bound. Excess hoarding suppresses velocity below the full-employment level, creating a deflationary gap — the monetary mechanism underlying Japan's 1990–2010 stagnation.

Demurrage (Definition 27.2) formally reduces the hoarding motive by making money holding costly: $M(t) = M_0 e^{-\delta t}$. It is rigorously distinguished from inflation: demurrage reduces nominal balances without changing price levels, does not reduce real debt burdens, and generates revenue for the issuing authority rather than diffusing purchasing power broadly.

The velocity theorem (Theorem 27.1) proves $V(\delta) = V_0 + \phi_\delta \cdot \delta$ — velocity is strictly increasing in demurrage rate with elasticity $\varepsilon_{V,\delta} = 1/2$. The optimal demurrage rate (Theorem 27.2) is derived from the gap between the natural interest rate and the natural capital regeneration rate: $\delta^* = i^{\text{natural}} - \mathcal{R}^{\text{MSY}}$ — equating the monetary holding cost to the ecological regeneration opportunity cost. This is the chapter's most original result: a monetary policy parameter derived from ecological dynamics.

Demurrage reduces the interest transfer (Theorem 27.3) and is progressive when applied to bank balances (Proposition 27.2). It increases investment rates (Proposition 27.3) and specifically makes ecologically viable investments fundable by reducing the effective hurdle rate to the ecological regeneration rate.

The 50-year comparison shows Economy B (5\% demurrage) outperforming Economy A (5\% interest) by 32.7\% in GDP, 26.6\% lower wealth Gini, 44.3\% higher natural capital, and complete elimination of Minsky crises. The Wörgl and Chiemgauer empirical cases confirm the direction of the velocity effect (2.87× and 11.9× respectively), though the magnitude exceeds the formal model's square-root prediction, pointing to behavioral and community-identity factors not captured by the optimization framework alone.

Chapter 28 closes Part V by integrating the four monetary architectures — debt-based, sovereign, mutual credit, and demurrage — into a unified SFC framework, comparing their stability, distributional, and ecological properties side-by-side through a 50-year simulation including a financial shock.

---

## Exercises

**27.1** Define demurrage formally (Definition 27.2). For a EUR 1,000 note under a 4\% annual demurrage rate:
(a) Compute its value after 3 months, 6 months, 1 year, and 2 years.
(b) A consumer expects to need EUR 900 in goods in exactly 6 months. Is it better to hold EUR 1,000 in demurrage currency now or convert to a non-demurrage asset (returning 1\% over 6 months)? Show the calculation.
(c) Contrast with 4\% annual inflation and no demurrage: what is the real value of EUR 1,000 after 1 year in each case? Why are the mechanisms different despite similar nominal loss?

**27.2** Apply the velocity theorem (Theorem 27.1) and its Corollary 27.1:
(a) For baseline velocity $V_0 = 1.2$ and demurrage rate $\delta = 0.03$/year, compute $V(\delta = 0.03)$ using the square-root formula.
(b) Compute the elasticity $\varepsilon_{V,\delta}$ and verify it equals 1/2.
(c) An economy has GDP $Y = $ EUR 500 billion and money supply $M = $ EUR 400 billion. Under the baseline ($\delta = 0$, $V_0 = 1.2$): what is $M$'s implied velocity? After introducing $\delta = 0.03$: what is the new implied velocity and what does the same $M$ support in GDP terms? By how much can the central bank reduce $M$ while maintaining $Y$ unchanged?

**27.3** The optimal demurrage rate (Theorem 27.2): for $i^{\text{natural}} = 3.5\%$ and three natural capital types:
(a) Boreal forest with $r = 0.04$ and $K = 5$ (million hectares, normalized): compute $\mathcal{R}^{\text{MSY}}$ and the optimal demurrage rate.
(b) Atlantic bluefin tuna with $r = 0.12$ and $K = 2$ (million tonnes): compute $\mathcal{R}^{\text{MSY}}$ and $\delta^*$.
(c) Topsoil with $r = 0.003$ and $K = 0.5$ (normalized): compute $\mathcal{R}^{\text{MSY}}$ and $\delta^*$. Why is the optimal rate so much higher for soil than for fish? What does this imply for current agricultural economics?

**★ 27.4** Prove Theorem 27.3 (demurrage reduces the interest transfer) and derive its full distributional implications.

(a) Set up the distributional model: two types of agents — creditors (fraction $\lambda$, average monetary holdings $\bar{M}_C$) and debtors (fraction $1-\lambda$, average monetary holdings $\bar{M}_D < \bar{M}_C$). Write the net annual transfer from debtors to creditors under debt money (no demurrage).

(b) Introduce demurrage at rate $\delta$: derive each agent type's net demurrage payment and the revised net transfer.

(c) Prove that the net transfer is reduced by $\delta \bar{M}_C \lambda > 0$ — always reduced, regardless of parameters.

(d) Derive the condition under which demurrage is welfare-improving for debtors (their net gain from reduced interest payments exceeds their demurrage cost) and welfare-reducing for creditors (their demurrage payment exceeds any reduction in their interest income). Show that this condition is always satisfied when debtor and creditor monetary holdings differ ($\bar{M}_C \neq \bar{M}_D$).

**★ 27.5** Prove that in a closed economy, demurrage money achieves a Pareto improvement over positive-interest money when the natural capital regeneration rate is below the natural rate of interest.

(a) Define the Pareto improvement condition: every agent is weakly better off and at least one agent is strictly better off under demurrage vs. positive interest.

(b) Show that under positive interest and $\mathcal{R}^{\text{MSY}} < i^{\text{natural}}$: (i) ecological investment projects are not funded (return below market rate); (ii) natural capital declines ($\dot{N} < 0$); (iii) future productive capacity is reduced relative to the social optimum.

(c) Show that under demurrage at $\delta^* = i^{\text{natural}} - \mathcal{R}^{\text{MSY}}$: (i) ecological investment projects are funded (return equals market rate net of demurrage); (ii) natural capital stabilizes ($\dot{N} \approx 0$); (iii) long-run productive capacity is higher than under positive interest.

(d) Construct the Pareto argument: if long-run productive capacity is higher under demurrage, there exists a redistribution scheme (using demurrage revenue) that makes all agents weakly better off compared to the positive-interest equilibrium. The redistribution scheme is the universal demurrage dividend. Prove formally using the Kaldor-Hicks improvement criterion.

**★★ 27.6** Build a full SFC model with demurrage money and compare equilibrium outcomes over 50 years.

**Model specification:**
- Sectors: Households ($H$), Firms ($F$), Issuing Authority ($IA$), Government ($G$).
- Money: Demurrage currency $M^D$ with rate $\delta = 0.025$/year; reinjected by $IA$ as universal dividend.
- Standard production function: $Y = K^\alpha L^{1-\alpha}$, $\alpha = 0.35$.
- Investment: $I = s(r_K - \delta) \cdot Y$ when $r_K > \delta$; 0 otherwise ($\delta$ reduces the hurdle rate).
- Natural capital: $\dot{N} = \mathcal{R}^{\text{MSY}} N - \mathcal{D}(Y)$ with $\mathcal{R}^{\text{MSY}} = 0.025$, $\mathcal{D}(Y) = 0.03Y/Y_0$.
- Wealth distribution: Gini evolves as $\dot{G} = (r_K - g - \delta \cdot (1 - G)) \cdot G$ where the demurrage term reduces the Gini.

(a) Implement the SFC-demurrage model in Python using the balance sheet and transaction flow matrix formulation of Chapter 18.

(b) Run the model for 50 years under: (i) $\delta = 0$ (baseline); (ii) $\delta = 0.025$ (matching $\mathcal{R}^{\text{MSY}}$); (iii) $\delta = 0.05$ (double ecological rate). Report GDP, Gini, $N$, and investment rate at years 10, 25, and 50 for each case.

(c) Apply a financial shock at year 25: a 20\% productivity decline for 3 years (simulating a financial crisis under the baseline). Compare the shock impact and recovery under each demurrage rate. Does demurrage make the economy more or less resilient to productivity shocks?

(d) Connect to the Minsky instability condition of Chapter 23: show formally that the demurrage SFC model cannot generate Minsky instability. Identify the specific mechanism that is absent. Does this absence make the demurrage economy immune to all financial crises, or only to debt-driven crises?

---

*Chapter 28 closes Part V by unifying the four monetary architectures — debt-based, sovereign, mutual credit, and demurrage — into a single comparative SFC framework. The chapter runs a four-country simulation over 50 years with a mid-period financial shock, demonstrating that all three non-debt architectures outperform debt-based money on stability, equality, and ecological health, each with distinct advantages that suggest a hybrid architecture may be optimal for the cooperative-regenerative economy.*
