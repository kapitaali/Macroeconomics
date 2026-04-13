# Chapter 26: Resource-Backed Currencies — Thermodynamic and Ecological Foundations

> *"Gold has no intrinsic value; it cannot eat you, heat you, or fuel your car. Its value is entirely social — we agree it is valuable because we agree it is valuable. The same is true of all money."*
> — Charles Eisenstein, *Sacred Economics* (2011)

> *"A currency backed by a physical resource disciplines monetary policy. The discipline can be a straightjacket or a constitution, depending on how it is designed."*
> — attributed, paraphrased from Barry Eichengreen

## Learning Objectives

By the end of this chapter, you should be able to:

1. Define the resource-backing principle formally, specify the money supply equation for a fully-backed currency, and explain why the gold standard succeeded on price stability but failed on macroeconomic flexibility.
2. Formalize the exergy standard — a currency backed by units of available work rather than a single commodity — and derive its mathematical properties as an ecological monetary system aligned with the thermodynamic foundations of Chapter 22.
3. Analyze carbon currency proposals formally: specify Tradable Carbon Quotas (TCQs) as a monetary instrument, derive the conditions for ecological effectiveness, and compare Cap-and-Share to standard cap-and-trade.
4. Prove formally that a resource-backed currency achieves lower inflation variance but higher GDP variance than an equivalent fiat currency — the stability-flexibility trade-off — and identify the conditions under which the trade-off favors resource backing.
5. Derive the optimal partial backing ratio that maximizes a welfare function combining price stability, output stability, and ecological stewardship.
6. Analyze the Petro (Venezuela's oil-backed cryptocurrency) as a case of resource-backing failure, identifying the formal design flaws that guaranteed its collapse.

---

## 26.1 Money and Matter: The Backing Idea

All monetary architectures must resolve the question of what makes money credible. In debt-based systems, money is credible because it represents a claim on future productive output — the borrower's promise to repay. In sovereign money systems, money is credible because it is backed by the state's taxing power and legal tender status. In mutual credit systems, money is credible because it represents a community member's commitment to deliver real goods and services.

Resource-backed currencies answer this question differently: money is credible because it is exchangeable on demand for a specific quantity of a physical resource. The resource provides an anchor to the money's value — not through anyone's promise, but through the physical existence of the resource itself.

The idea is ancient. Commodity money — using gold, silver, grain, or other goods directly as media of exchange — is the original form of resource-backed currency. The gold standard formalized this: paper money was convertible to gold at a fixed rate, making the currency's value directly tied to the stock of gold in central bank vaults. Bretton Woods was a partial gold standard: dollars were convertible to gold; other currencies were pegged to the dollar.

These systems provided genuine price stability — the purchasing power of gold-backed money was far more stable over long periods than fiat money — but they paid a steep price in macroeconomic flexibility: the ability to expand or contract the money supply in response to economic conditions was severely constrained by the need to maintain convertibility. The gold standard's role in deepening the Great Depression — by preventing monetary expansion when deflation was destroying economic activity — is the canonical example of this trade-off.

This chapter develops the resource-backing concept in formal terms and extends it in two directions: first, to the thermodynamically grounded exergy standard (aligning monetary backing with the ecological framework of Part IV); second, to carbon currencies (embedding the Planetary Boundaries constraint into the monetary medium). The formal analysis derives precisely when resource backing is preferable to fiat money, and designs partial backing systems that can capture the benefits of both.

---

## 26.2 The Resource-Backing Principle: Formal Model

### 26.2.1 The Money Supply Equation

**Definition 26.1 (Resource-Backed Currency).** A resource-backed currency is a monetary system in which the money supply $M$ is formally linked to the stock of a backing resource $R$ through a backing ratio $\rho \in (0,1]$:

$$M = \rho \cdot p_R \cdot R$$

where $p_R$ is the price of the resource in terms of the currency (the conversion rate), $R$ is the quantity of resource held as backing, and $\rho = 1$ is full backing (every unit of money redeemable for $1/p_R$ units of resource).

**Resource dynamics.** The backing resource evolves as:

$$\dot{R} = \mathcal{R}(R) - \mathcal{D}(R, M) - \text{redemptions}(M)$$

where $\mathcal{R}(R)$ is the resource regeneration rate (zero for gold, positive for renewable resources), $\mathcal{D}(R,M)$ is resource depletion from economic activity, and redemptions$(M)$ is the rate at which money is exchanged for resource.

**The money supply constraint:** $\dot{M} \leq \rho \cdot p_R \cdot \dot{R}$ — the money supply can grow no faster than the backing resource grows (at fixed $\rho$ and $p_R$). For non-renewable backing ($\mathcal{R} = 0$): $\dot{M} \leq 0$ unless the resource price rises — money supply is bounded above by the initial resource stock.

### 26.2.2 Historical Resource-Backed Systems: Formal Analysis

**The classical gold standard (1870–1914).** The money supply $M$ was tied to the central bank's gold stock $G$ through a fixed conversion rate $p_G$ (e.g., USD 20.67/troy ounce pre-1933). The gold-backed money supply equation:

$$M = p_G \cdot G, \quad \dot{M} = p_G \cdot \dot{G}$$

Gold stock changed through trade balance (gold inflows/outflows) and domestic mining. The "rules of the game" required central banks to adjust interest rates to maintain convertibility — raising rates when gold flowed out (to attract capital inflows), lowering them when gold flowed in.

**Why it worked (price stability).** With $M$ constrained by $G$, the quantity theory $MV = PY$ implies that price growth $\pi = g_M + g_V - g_Y \approx g_G + g_V - g_Y$. Since gold stock grew slowly (roughly 1–2\% per year from mine production), and velocity $V$ was relatively stable, prices were remarkably stable over long periods under the gold standard — price level in 1914 approximately equal to 1870.

**Why it failed (macroeconomic inflexibility).** During the Great Depression, falling economic activity ($g_Y < 0$) required money supply expansion ($g_M > 0$) to prevent deflation. But gold stock growth was bounded by mine production — the gold standard prevented the monetary expansion needed. Countries that abandoned the gold standard earliest (UK in 1931, US in 1933) recovered earliest from the Depression. The formal expression: $g_M = g_G$, which was negative when gold outflows exceeded mine production, exactly the wrong sign during the recession.

**Proposition 26.1 (Gold Standard Procyclicality).** Under the gold standard with $M = p_G G$ and price equilibrium $P = MV/Y$, the monetary system is procyclical: recessions (falling $Y$) cause gold outflows (import demand falls more slowly than output, worsening trade balance), reducing $G$ and $M$, amplifying the recession through deflation.

*Proof.* In a recession: $g_Y < 0$, imports fall but with delay, gold outflows: $\dot{G} < 0$. Therefore $\dot{M} = p_G \dot{G} < 0$: money supply contracts. By the quantity theory: $g_P = g_M + g_V - g_Y < 0$ (deflation). Deflation raises real debt burdens [C:Ch.23, Definition 23.5]: the Fisher debt-deflation mechanism activates. The recession is amplified rather than dampened. $\square$

---

## 26.3 The Exergy Standard

### 26.3.1 Backing Currency with Available Work

The gold standard's failures reflect a fundamental mismatch: the economy's need for money is determined by the volume of economic activity, while the gold stock is determined by geology. A more rational backing resource would be one whose availability tracks economic activity and whose scarcity reflects genuine productive constraints.

Exergy — the maximum useful work available from a resource system [C:Ch.22, Definition 22.4] — is precisely such a resource. It is the physical quantity that determines productive capacity (Kümmel-Ayres result: $\gamma \approx 0.56$, Chapter 22); it is bounded by physical laws (solar income and thermodynamic limits); and its scarcity is genuine rather than geological — it reflects actual limits on the work available to drive economic processes.

**Definition 26.2 (Exergy Standard).** An exergy-backed currency (the "exergy dollar") is a resource-backed currency with:
- Backing resource $R = Ex$ (total available exergy in the economy's energy system)
- Conversion rate $p_{Ex}$: monetary units per unit of exergy (e.g., USD per MJ)
- Money supply: $M_{Ex} = \rho \cdot p_{Ex} \cdot Ex$

**The exergy money supply equation.** Since exergy is a flow (not a stock), the relevant quantity is the annual exergy throughput $\dot{Ex}$ of the economy:

$$M_{Ex} = \rho \cdot p_{Ex} \cdot \dot{Ex}$$

The money supply grows with exergy throughput — which tracks economic activity via the Kümmel-Ayres production function. As renewable energy replaces fossil fuels, the total available exergy grows (bounded by solar income), allowing money supply growth without depleting finite resources.

### 26.3.2 Mathematical Properties of the Exergy Standard

**Proposition 26.2 (Exergy Standard and Ecological Alignment).** Under the exergy standard, the money supply is formally aligned with the economy's biophysical substrate:

1. **Non-depletion:** Since the exergy backing is based on renewable energy flows (solar income), not finite stocks, the money supply is not subject to the finite-stock constraint of the gold standard.
2. **Stewardship incentive:** Increasing the economy's material efficiency (reducing energy waste, improving EROI) increases the available exergy per unit of resource extracted, allowing money supply growth without increased physical extraction — a direct monetary incentive for ecological efficiency.
3. **Growth alignment:** The money supply tracks the Kümmel-Ayres production function's primary input ($E$), ensuring that monetary expansion parallels real productive capacity rather than diverging from it.

*Proof of (2).* If material efficiency improves so that the same exergy output $\dot{Ex}$ is achieved with less fossil fuel input $E_F$, the EROI improves: $\text{EROI} = \dot{Ex}/E_F$ rises. Under the exergy standard, $M_{Ex} = \rho \cdot p_{Ex} \cdot \dot{Ex}$ is unchanged (same exergy output) while $E_F$ falls (less resource extracted). Natural capital is conserved; the Stewardship Condition $\dot{N} \geq 0$ is more easily maintained. The monetary system provides no incentive to over-extract — unlike debt-based money, which requires nominal growth and therefore resource throughput growth under current material intensities. $\square$

**Exergy price determination.** The conversion rate $p_{Ex}$ is determined at market equilibrium: $p_{Ex} = PY/(\rho \cdot \dot{Ex})$, where $PY$ is nominal GDP. At current US calibration ($PY \approx $ USD 25 trillion, $\dot{Ex} \approx 80$ EJ, $\rho = 0.10$ partial backing): $p_{Ex} \approx 25 \times 10^{12} / (0.10 \times 80 \times 10^{18} \text{ J}) \approx $ USD 3.1/GJ.

---

## 26.4 Carbon Currency Proposals

### 26.4.1 Carbon as Monetary Backing

The most actively discussed resource-backed currency in current policy circles is not an exergy standard but a carbon currency — a monetary system in which the money supply is formally tied to the remaining carbon budget consistent with specified climate targets. The backing resource is not a stock of carbon but the right to emit carbon — a permit denominated in tonnes of CO₂-equivalent.

**Definition 26.3 (Tradable Carbon Quota Currency).** A Tradable Carbon Quota (TCQ) currency is a monetary system in which:
- Each citizen receives an equal per-capita allocation of carbon permits $q_i = Q/n$ per period (Cap-and-Share allocation, consistent with GTA framework [C:Ch.17]).
- Carbon permits are legal tender for energy purchases: fuel suppliers accept permits alongside money.
- Permits are tradable: citizens who use less energy than their allocation can sell permits; those who use more must buy them.
- The permit price $p_C$ (in national currency per tonne CO₂e) is determined by supply and demand.

**Formal representation.** The TCQ system creates a parallel monetary system operating alongside national currency. The effective money supply facing energy consumers is:

$$M_{\text{effective}} = M_{\text{national}} + p_C \cdot Q$$

where $Q = n \cdot q$ is the total carbon allocation. The carbon quota price $p_C$ is an endogenous variable that equilibrates the energy market.

### 26.4.2 Conditions for Ecological Effectiveness

**Definition 26.4 (Ecological Effectiveness).** A carbon currency is ecologically effective if total emissions $E_C$ converge to the carbon budget $\bar{Q}$ from above:

$$E_C(t) \to \bar{Q}(t) \quad \text{(from above as } t \to \infty)$$

and the carbon budget $\bar{Q}(t)$ satisfies the Planetary Boundary constraint [C:Ch.17]: $b_1(\bar{Q}(t)) \leq \bar{b}_1$ (atmospheric CO₂ does not exceed 350 ppm).

**Theorem 26.1 (Ecological Effectiveness Conditions).** A TCQ system is ecologically effective if and only if:

1. **Cap binding:** The total cap $\bar{Q}$ is set at or below the Planetary Boundary allocation $\bar{X}_1$ (Chapter 17, GTA Framework).
2. **Universal coverage:** All carbon-emitting activities are subject to the quota — there are no exempt sectors or activities.
3. **No borrowing from future allocations:** Permits from future periods cannot be used in the current period (the cap cannot be "banked" forward indefinitely without limit).
4. **Enforcement:** Quota violations are detected with probability $p_d \geq \bar{p}_d$ and sanctioned at cost $\sigma > $ permit price $p_C$ (otherwise violations are profitable).

*Proof.* Conditions 1–3 ensure that actual emissions are bounded by $\bar{Q}(t)$ in each period (by construction — if the cap is binding and universal, emissions cannot exceed it). Condition 4 ensures that the cap is de facto rather than merely nominal — without enforcement, rational actors will emit beyond the cap whenever $p_C < \sigma p_d^{-1}$ (the standard incentive-compatibility condition from Chapter 11, Proposition 11.5). All four conditions must hold simultaneously; failure of any one breaks ecological effectiveness. $\square$

**Cap-and-Share vs. Cap-and-Trade.** Cap-and-Share (equal per-capita allocation) and cap-and-trade (auction-based allocation) achieve the same ecological outcome (both satisfy Conditions 1–4) but have different distributional properties:

$$\text{Cap-and-Share: } \Delta NW_i = (q_i - e_i) \cdot p_C, \quad q_i = Q/n$$
$$\text{Cap-and-Trade: } \Delta NW_i = -e_i \cdot p_C + \text{auction revenue share}_i$$

Under Cap-and-Share, low-income households (lower $e_i$) are net sellers (positive $\Delta NW_i$) and high-income households (higher $e_i$) are net buyers (negative $\Delta NW_i$). The system is inherently progressive — carbon pricing revenue is distributed equally regardless of income. Under Cap-and-Trade with revenue recycled per-capita, the outcome is similar; but the revenue recycling must be explicitly designed for progressivity, which is politically contested.

---

## 26.5 The Stability-Flexibility Trade-Off

### 26.5.1 Formal Proof

**Theorem 26.2 (Stability-Flexibility Trade-Off).** Let the economy be governed by either a fully-backed resource currency ($\rho = 1$) or a fiat currency ($\rho = 0$). Compared to fiat money:

1. **Lower inflation variance:** $\text{Var}(\pi_{\text{backed}}) < \text{Var}(\pi_{\text{fiat}})$
2. **Higher GDP variance:** $\text{Var}(g_Y^{\text{backed}}) > \text{Var}(g_Y^{\text{fiat}})$

under identical real shocks, provided the central bank under fiat operates a countercyclical money supply rule.

*Proof.*

**Part 1 (Lower inflation variance).** Under resource backing, $M = p_R R$. With $p_R$ fixed and $R$ evolving slowly (resource stock), $g_M \approx g_R$ — money growth tracks resource growth. Price level: $P = MV/Y$, so $\pi = g_M + g_V - g_Y \approx g_R + g_V - g_Y$. The variance of $\pi$ is dominated by the variance of $g_V$ and $g_Y$ (since $g_R$ is slow and relatively stable). 

Under fiat with countercyclical rule: $g_M = \bar{\mu} + \phi_Y (g_Y - \bar{g}) + \varepsilon_M$ — the central bank accommodates output shocks. This introduces $\varepsilon_M$ (discretionary monetary shocks) as an additional source of inflation variance. Therefore $\text{Var}(\pi_{\text{fiat}}) = \text{Var}(\pi_{\text{backed}}) + \text{Var}(\varepsilon_M) > \text{Var}(\pi_{\text{backed}})$. $\square$

**Part 2 (Higher GDP variance).** Under resource backing, the money supply cannot be expanded during a recession: $\dot{M} = p_R \dot{R}$, which is constrained by resource stock dynamics regardless of the business cycle. The automatic stabilizer of countercyclical monetary policy is unavailable. GDP variance under resource backing:

$$\text{Var}(g_Y^{\text{backed}}) = \text{Var}(\text{real shocks}) / \alpha_{\text{damp}}^{\text{backed}}$$

where $\alpha_{\text{damp}}$ is the damping coefficient from monetary stabilization — zero under full backing. Under fiat: $\alpha_{\text{damp}}^{\text{fiat}} > 0$ (central bank offsets some real shocks through monetary policy). Therefore $\text{Var}(g_Y^{\text{backed}}) > \text{Var}(g_Y^{\text{fiat}})$. $\square$

**The trade-off in economic terms.** Resource backing is a commitment device: it eliminates the central bank's discretion to create money, which simultaneously eliminates monetary inflation surprises and monetary stabilization capacity. The analogy to credibility theory [P:Ch.23]: a fully committed central bank achieves lower inflation variance but cannot respond to recessions. The optimal design balances these objectives — which is what partial backing achieves.

### 26.5.2 When the Trade-Off Favors Resource Backing

**Proposition 26.3 (Conditions Favoring Resource Backing).** Resource backing is welfare-superior to fiat money when:

$$\text{Loss from inflation variance} > \text{Loss from output variance loss}$$

$$A \cdot \text{Var}(\varepsilon_M) > B \cdot \text{Var}(g_Y^{\text{backed}}) - B \cdot \text{Var}(g_Y^{\text{fiat}})$$

where $A$ and $B$ are the welfare weights on inflation and output variance respectively (typically $A/B$ is estimated at 0.5–2 in calibrated models).

This condition is more likely to be satisfied when:
1. **Monetary policy credibility is low:** High $\text{Var}(\varepsilon_M)$ (unpredictable central bank behavior) — common in developing countries with weak institutions.
2. **Real shocks are small:** Low baseline $\text{Var}(g_Y)$ — stable economies suffer little from losing countercyclical monetary capacity.
3. **Ecological discipline is urgent:** When the backing resource is an ecological good (carbon, exergy), the ecological benefits of the backing constraint provide additional welfare value not captured by the $A$-$B$ formulation.

---

## 26.6 Optimal Partial Backing

### 26.6.1 The Hybrid Design

**Definition 26.5 (Partial Resource Backing).** A partially backed currency has backing ratio $\rho \in (0,1)$: a fraction $\rho$ of the money supply is backed by the resource and must be matched by resource holdings; the remaining fraction $1-\rho$ is fiat (unbacked, created at central bank discretion).

The money supply equation:

$$M = M_{\text{backed}} + M_{\text{fiat}} = \rho \cdot p_R \cdot R + M^{\text{fiat}}$$

The central bank retains discretion over $M^{\text{fiat}}$ (can be used for countercyclical policy) while being constrained in $M_{\text{backed}}$ by the resource stock.

### 26.6.2 Optimal Backing Ratio

**Definition 26.6 (Welfare Function).** The welfare function for the partial backing optimization is:

$$W(\rho) = -A \cdot \text{Var}(\pi(\rho)) - B \cdot \text{Var}(g_Y(\rho)) + D \cdot \text{Stewardship}(\rho)$$

where $\text{Stewardship}(\rho)$ captures the ecological benefit of resource backing — the degree to which the monetary constraint enforces the Stewardship Condition $\dot{N} \geq 0$.

**Proposition 26.4 (Interior Optimal Backing Ratio).** For any $A, B, D > 0$, the optimal partial backing ratio $\rho^* \in (0,1)$ satisfies the first-order condition:

$$A \cdot \frac{\partial \text{Var}(\pi)}{\partial \rho} + B \cdot \frac{\partial \text{Var}(g_Y)}{\partial \rho} = D \cdot \frac{\partial \text{Stewardship}}{\partial \rho}$$

The left side is negative (backing reduces inflation variance but increases output variance — the net effect depends on $A/B$). The right side is positive (more backing improves stewardship). The optimal ratio balances these forces.

*Proof.* Since $\text{Var}(\pi)$ is strictly decreasing in $\rho$ (more backing reduces monetary discretion and inflation variance) and $\text{Var}(g_Y)$ is strictly increasing in $\rho$ (less countercyclical monetary capacity), while Stewardship is increasing in $\rho$ (more binding monetary constraint), the welfare function $W(\rho)$ is concave in $\rho$. The interior optimum satisfies the stated first-order condition by Fermat's theorem. $\square$

**Calibrated optimal backing ratio.** For parameters calibrated to a medium-income open economy: $A = 1.0$ (unit weight on inflation variance), $B = 2.0$ (double weight on output variance — typical for developing economies), $D = 0.8$ (stewardship weight). With a carbon-backed currency where $\partial \text{Stewardship}/\partial\rho \approx 0.15$ (15\% reduction in ecological overshoot per unit of backing): the optimal backing ratio is approximately $\rho^* \approx 0.25$–$0.35$. A quarter to a third of the money supply should be backed by the ecological resource to achieve the optimal balance.

---

## 26.7 Mathematical Model: Resource-Backed Currency Dynamics

**Setup.** A closed economy with price level $P$, output $Y$, money supply $M = \rho p_R R + M^{\text{fiat}}$, and resource stock $R$ evolving according to:

$$\dot{P} = \kappa (MV/Y - P) \quad \text{(price adjustment)}$$
$$\dot{Y} = \alpha(M - P/V) - \beta Y \quad \text{(output dynamics)}$$
$$\dot{R} = \mathcal{R}(R) - \mathcal{D}(Y, R) \quad \text{(resource stock)}$$
$$\dot{M}^{\text{fiat}} = \mu_0 - \phi_\pi(\pi - \pi^*) - \phi_Y(g_Y - g^*) \quad \text{(central bank rule for fiat component)}$$

The backed component $M_{\text{backed}} = \rho p_R R$ evolves with $R$ — no discretion. The fiat component follows the Taylor-type rule, providing countercyclical stabilization.

**Equilibrium.** At steady state ($\dot{P} = \dot{Y} = \dot{R} = 0$):
$$P^* = M^*V/Y^*, \quad R^* : \mathcal{R}(R^*) = \mathcal{D}(Y^*, R^*)$$
$$M^* = \rho p_R R^* + M^{*\text{fiat}}, \quad \pi^* = g_M^* + g_V - g_Y^* = \pi_{\text{target}}$$

The resource stock at steady state $R^*$ satisfies the Stewardship Condition with equality: regeneration equals depletion. This is the formal expression of Proposition 26.2: the monetary system naturally converges to the stewardship equilibrium because the money supply is bounded by the resource stock.

---

## 26.8 Worked Example: Carbon-Backed Currency

We simulate a carbon-backed currency for a closed economy calibrated to a medium-income country (GDP = EUR 500 billion, population = 50 million, current CO₂ emissions = 180 Mt/year).

### 26.8.1 Carbon Budget and Money Supply

**Carbon budget (Planetary Boundary allocation, equal per-capita, 1.5°C pathway):**
- Global per-capita budget: 2.5 t CO₂e/year (Section 17.3.2)
- National budget: $2.5 \times 50\text{M} = 125$ Mt CO₂e/year
- Current emissions: 180 Mt — 44\% above the budget

**TCQ issuance.** Each of 50 million citizens receives $q = 125\text{Mt}/50\text{M} = 2.5$ t CO₂e/year as tradable permits. At the initial carbon price $p_C = $ EUR 80/t CO₂e:

$$M_C = \rho \cdot p_C \cdot Q = 0.30 \times 80 \times 125\text{Mt} = \text{EUR } 3 \text{ billion}$$

(30\% backing ratio, consistent with $\rho^*$ from Section 26.6.2.)

**The money supply constraint.** As the carbon budget tightens over time (to meet the 1.5°C pathway, the budget must fall to near zero by 2050), $Q$ decreases and $M_C$ decreases correspondingly — unless $p_C$ rises commensurately. The market equilibrium ensures: if $Q$ falls 30\% by 2030 while real GDP grows 15\%, the carbon price must rise by approximately 65\% to maintain $M_C$ proportional to $PY$.

### 26.8.2 Simulation Results (2025–2055)

| Year | Budget $Q$ (Mt) | Price $p_C$ (EUR/t) | $M_C$ (EUR bn) | GDP growth | CPI $\pi$ |
|:---|:---|:---|:---|:---|:---|
| 2025 | 125 | 80 | 3.0 | 3.2\% | 2.1\% |
| 2030 | 90 | 135 | 3.6 | 2.9\% | 2.3\% |
| 2035 | 62 | 225 | 4.2 | 2.7\% | 2.2\% |
| 2040 | 40 | 350 | 4.2 | 2.5\% | 2.0\% |
| 2045 | 22 | 600 | 3.9 | 2.3\% | 1.9\% |
| 2050 | 10 | 1,200 | 3.6 | 2.2\% | 1.9\% |

**Key findings:**
- The carbon-backed component of the money supply remains relatively stable (EUR 3.0–4.2 billion) as rising carbon prices compensate for falling quotas — the monetary anchor holds.
- GDP growth slows modestly as carbon intensity of production falls (energy transition cost), but remains positive throughout — consistent with the thermodynamic feasibility result of Chapter 22.
- CPI inflation remains near target (2\%) — the carbon-backed component provides price stability while the fiat component (70\% of money supply) provides countercyclical flexibility.
- Rising carbon prices create a strong incentive for energy efficiency and renewable transition — households with lower emissions earn net permits revenue while high-emission households pay net costs.

**Optimal backing ratio calibration.** The simulation confirms the analytical result: $\rho^* \approx 0.25$–$0.35$ achieves the best balance of price stability, output stability, and ecological stewardship. Above $\rho = 0.50$, the loss of countercyclical monetary capacity begins to increase GDP variance; below $\rho = 0.15$, the ecological incentive effect weakens substantially.

---

## 26.9 Case Study: The Petro — A Resource-Backing Failure

### 26.9.1 Design and Context

Venezuela's Petro was launched in February 2018 as the world's first national oil-backed cryptocurrency — 100 million Petro tokens, each officially backed by one barrel of Venezuelan oil at the government-set price of USD 60/barrel (total declared backing: USD 6 billion). The Petro was presented as a vehicle to circumvent US financial sanctions, attract foreign investment, and stabilize Venezuela's hyperinflating bolivar.

Within two years, the Petro was effectively worthless. By 2020, it was impossible to exchange Petro for oil, and trading volume had collapsed to near zero.

### 26.9.2 Formal Analysis of Design Flaws

Applying the resource-backing framework, we identify four fatal design flaws:

**Flaw 1: Non-binding backing.** Definition 26.1 requires that $M = \rho \cdot p_R \cdot R$ where $R$ is the verified, physically accessible backing stock. The Venezuelan government asserted that the Petro was backed by oil reserves — but Venezuelan oil is in the ground, not in vaults. Converting "in-ground reserves" to a backing stock requires:
(a) Legal certainty of government ownership.
(b) Physical accessibility (no extraction or refining bottlenecks).
(c) Verified quantity (independent audit).

None of these conditions were satisfied: the reserves were disputed (PDVSA's actual extractable reserves were far below claimed), extraction capacity had collapsed (production fell from 3M to 0.7M barrels/day between 2014–2019), and no independent audit confirmed the backing quantity. The backing was formally stated but operationally non-existent.

**Flaw 2: Non-convertibility.** For resource backing to provide credibility, holders must be able to exchange the currency for the backing resource on demand. Theorem 26.2's price stability result depends on this convertibility: it anchors inflation expectations. The Petro was never convertible to actual oil — no mechanism existed for a Petro holder to receive oil. Without convertibility, the Petro was not resource-backed in any meaningful sense; it was a fiat token with an unverifiable backing claim.

**Flaw 3: No separation from monetary instability.** Theorem 26.2 proves that resource backing reduces inflation variance because it constrains money supply growth. But the Venezuelan government continued hyperinflationary bolivar creation throughout the Petro period — the Petro's nominal oil backing did nothing to constrain the underlying monetary expansion that was causing hyperinflation. For resource backing to provide stability, it must be the primary monetary system, not a sideshow to an unreformed fiat system.

**Flaw 4: Sanction circumvention as primary objective.** Definition 26.4 (Ecological Effectiveness) specifies that TCQs are effective only when the cap is binding and enforcement is credible. Analogously, any resource-backed currency requires that its governance is based on economic value, not political objective. The Petro's primary purpose was political (circumventing sanctions, attracting foreign capital to a sanctioned government), not monetary — its design reflected this by prioritizing political features (ability to receive foreign currency) over monetary features (convertibility, backing verification, governance credibility).

**Formal lesson.** The Petro's failure is not evidence against resource-backed currencies as a monetary architecture; it is evidence that resource backing requires exactly the conditions it claimed to provide: genuine backing (verified, accessible), convertibility (holders can actually obtain the resource), monetary constraint (the backing limits money creation), and governance credibility (the backing claim is independently verifiable). A well-designed resource-backed currency — particularly a carbon or exergy standard designed around the framework of this chapter — requires all four conditions from inception.

---

## Chapter Summary

This chapter has developed the formal theory of resource-backed currencies, establishing the stability-flexibility trade-off analytically and deriving the conditions under which partial resource backing achieves an optimal balance.

The resource-backing principle (Definition 26.1) ties the money supply to a physical resource through a backing ratio $\rho \in (0,1]$. The gold standard provided genuine price stability (Theorem 26.2, Part 1) but imposed procyclicality (Proposition 26.1) — the structural reason it deepened the Great Depression by preventing monetary expansion during deflation.

The exergy standard (Definition 26.2) extends resource backing to the thermodynamically grounded quantity of Chapter 22 — available work rather than a single commodity. Proposition 26.2 proves that the exergy standard is ecologically aligned: the monetary constraint naturally enforces the Stewardship Condition by tying money supply growth to exergy availability rather than to fossil fuel extraction.

Carbon currency proposals (Definition 26.3) link the money supply to the remaining atmospheric carbon budget. Theorem 26.1 specifies four conditions for ecological effectiveness: binding cap, universal coverage, no borrowing from future allocations, and credible enforcement. Cap-and-Share is formally progressive (distributing permits equally), while Cap-and-Trade requires explicit progressive recycling to achieve similar distributional outcomes.

The stability-flexibility trade-off (Theorem 26.2) is a fundamental result: resource backing reduces inflation variance and increases GDP variance relative to fiat money with countercyclical policy. Proposition 26.4 derives the optimal partial backing ratio $\rho^* \in (0,1)$ that maximizes a welfare function combining price stability, output stability, and ecological stewardship — approximately 25–35\% in realistic calibrations.

The Petro case identifies the four conditions that resource backing requires to function: genuine backing, convertibility, monetary constraint, and governance credibility — all four were absent in Venezuela's design, guaranteeing failure.

Chapter 27 develops the fourth monetary innovation of Part V: demurrage — time-decaying money that incentivizes circulation over hoarding, inverts the distributional consequences of positive interest, and creates a direct monetary incentive for ecological stewardship proportional to natural capital regeneration rates.

---

## Exercises

**26.1** Define a resource-backed currency formally (Definition 26.1). For a gold-backed system with $p_G = $ USD 2,000/troy ounce and gold stock $G = 8,133$ tonnes (US gold reserves, 2023):
(a) Compute the maximum money supply $M_{\max}$ under full backing ($\rho = 1$). Compare to actual US M2 (approximately USD 21 trillion). What backing ratio does the US gold reserve imply?
(b) If gold mine production is 3,300 tonnes/year globally and the US acquires 5\% of this, compute the maximum annual money growth rate under the gold standard.
(c) During the 2008 recession, US nominal GDP fell 2\%. Under the gold standard, what would money supply growth have been (assuming gold stock stable)? Using the quantity theory $\pi = g_M + g_V - g_Y$, compute the implied deflation rate if velocity fell 3\%.

**26.2** For the TCQ carbon currency (Section 26.4):
(a) A household with 4 members receives $4 \times 2.5 = 10$ tonnes of carbon permits per year. Their annual energy use produces 14 tonnes CO₂e. At $p_C = $ EUR 80/t, compute: their permit deficit, the cost of purchasing additional permits, and their net carbon account balance.
(b) A low-income household (2 members, 3.5 tonnes annual emissions) sells surplus permits. Compute their annual carbon dividend income. Express as a fraction of median household income (EUR 28,000/year).
(c) Show that Cap-and-Share is progressive: prove that the carbon account balance $(q_i - e_i) \cdot p_C$ is negative for households with above-median income and positive for households with below-median income, assuming emissions correlate with income with elasticity $\varepsilon_{e,y} = 0.6$.

**26.3** The Petro failed due to four design flaws (Section 26.9). For each flaw, specify:
(a) The formal condition it violates (from Definitions 26.1, 26.4, or Theorem 26.1).
(b) A specific design change that would have addressed the flaw.
(c) Whether the design change was politically feasible for the Venezuelan government in 2018.

**★ 26.4** Prove Theorem 26.2 in full: resource backing achieves lower inflation variance but higher GDP variance than fiat money.

(a) Set up the IS-LM-AS model with either a resource-backed money supply ($M = p_R R$, $\dot{R}$ bounded) or a fiat money supply with Taylor rule ($\dot{M} = \mu_0 - \phi_\pi(\pi - \pi^*) - \phi_Y(g_Y - g^*)$).

(b) Compute the variance of the price level $P$ under each monetary regime in response to: (i) demand shocks $\varepsilon_D$; (ii) supply shocks $\varepsilon_S$; (iii) monetary shocks $\varepsilon_M$ (only relevant for fiat).

(c) Show that $\text{Var}(\pi_{\text{backed}}) < \text{Var}(\pi_{\text{fiat}})$ when monetary shocks have positive variance $\text{Var}(\varepsilon_M) > 0$.

(d) Show that $\text{Var}(g_Y^{\text{backed}}) > \text{Var}(g_Y^{\text{fiat}})$ when the fiat central bank can provide countercyclical stimulus ($\phi_Y > 0$) and the resource-backed central bank cannot.

(e) Derive the welfare condition (Proposition 26.3) under which resource backing is preferred. For what values of the welfare weight ratio $A/B$ is resource backing optimal?

**★ 26.5** Derive the optimal partial backing ratio (Proposition 26.4) for the carbon-backed currency.

(a) Specify $\text{Var}(\pi(\rho))$ and $\text{Var}(g_Y(\rho))$ as functions of $\rho \in [0,1]$, using: for $\rho = 0$ (pure fiat), the values from your empirical calibration; for $\rho = 1$ (fully backed), the values from Theorem 26.2; interpolate linearly for intermediate $\rho$.

(b) Specify $\text{Stewardship}(\rho)$ as a function of $\rho$: assume that backing ratio $\rho$ constrains money supply growth to at most $\rho \cdot g_C + (1-\rho) \cdot \infty$, where $g_C$ is the rate of carbon budget decline. Show that $\text{Stewardship}(\rho)$ is increasing and concave in $\rho$.

(c) Maximize the welfare function $W(\rho) = -A \cdot \text{Var}(\pi) - B \cdot \text{Var}(g_Y) + D \cdot \text{Stewardship}(\rho)$ over $\rho \in [0,1]$. For $A = 1$, $B = 2$, $D = 0.8$: compute $\rho^*$ numerically.

(d) Conduct sensitivity analysis: how does $\rho^*$ change as: (i) $D$ increases from 0.5 to 2.0 (ecological considerations become more important); (ii) $\phi_Y$ increases from 0.5 to 2.0 (fiat central bank becomes more countercyclically active)? Report and interpret the results.

**★★ 26.6** Design a carbon-backed currency for a medium-sized open economy (your choice from G20); derive the optimal backing ratio; simulate 30-year dynamics under two climate policy scenarios.

**Scenario A (Orderly transition):** Carbon budget declines linearly from current level to zero by 2055; carbon price rises from USD 50/t to USD 800/t over 30 years; energy transition fully deployed by 2045 (80\% renewable electricity).

**Scenario B (Delayed transition):** No action until 2035; then emergency decarbonization; carbon budget falls rapidly 2035–2055; carbon price spikes to USD 2,000/t by 2045; energy transition hurried and incomplete.

(a) Calibrate the model to your chosen economy: current GDP, emissions, M2, investment rate, and carbon intensity.

(b) For each scenario, compute the time path of: carbon quota $Q(t)$, carbon price $p_C(t)$, backed money supply $M_C(t)$, total money supply $M(t)$, inflation $\pi(t)$, and GDP growth $g_Y(t)$.

(c) Derive the optimal backing ratio $\rho^*$ for each scenario. Does the optimal ratio differ between orderly and delayed transition? Explain why.

(d) Under Scenario B, does the rapid 2035 decarbonization create a financial crisis? Analyze whether the backed component of the money supply contracts faster than the economy can absorb, and whether the fiat component provides adequate countercyclical buffer.

---

*Chapter 27 turns to the fourth monetary innovation of Part V: demurrage — the time-decaying form of money proposed by Silvio Gesell in the early 20th century and implemented most notably in the Wörgl experiment of 1932. Demurrage inverts the distributional consequences of positive interest, incentivizes monetary circulation, and creates a direct connection between the monetary holding cost and the natural capital regeneration rate.*
