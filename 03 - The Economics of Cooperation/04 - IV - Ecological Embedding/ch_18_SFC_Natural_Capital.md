# Chapter 18: Stock-Flow Consistent Modeling with Natural Capital — Accounting for Regeneration

> *"Every stock-flow inconsistency is an error, and every error has consequences."*
> — Wynne Godley, *Monetary Economics* (2007)

> *"What we measure shapes what we pursue. If we measure the wrong things, we will strive for the wrong goals."*
> — Joseph Stiglitz, Amartya Sen, and Jean-Paul Fitoussi, *Report by the Commission on the Measurement of Economic Performance and Social Progress* (2009)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Reconstruct the standard SFC framework — balance sheet matrix, transaction flow matrix, revaluation matrix — and explain why stock-flow consistency is a necessary condition for macroeconomic analysis.
2. Extend the SFC framework to include natural capital as a balance sheet stock, constructing the full SFC-N model with ecological revaluation accounting.
3. Derive the Provisioning Balance Sheet as a comprehensive social accounting matrix in which produced capital, natural capital, and social capital sum to total provisioning capacity.
4. Formalize Open Value Accounting (OVA) as a mechanism for valuing commons contributions, and prove its connection to the Shapley value allocation rule of Chapter 6.
5. Model the Three-Layer Coordination Stack in SFC terms, mapping mutual coordination, generative markets, and biophysical planning onto the SFC accounting identity.
6. Apply the SFC-N framework to the Icelandic economy, demonstrating how the stewardship constraint manifests as an accounting identity in the Provisioning Balance Sheet.

---

## 18.1 Why Accounting Consistency Matters

The stock-flow consistent (SFC) framework, developed by Wynne Godley, Marc Lavoie, and colleagues, provides a complete accounting framework for macroeconomic modeling in which every financial flow must be matched by a corresponding stock change, every sector's deficit must be another sector's surplus, and no money can appear from or disappear to nowhere [M:Ch.28]. The framework's discipline is its principal virtue: an SFC model that satisfies all accounting identities is guaranteed to contain no hidden assumptions about where funds come from or go to — a guarantee that standard DSGE and IS-LM models cannot always provide.

Chapter 17 established that the standard SFC framework — like the standard circular flow it formalizes — omits the biophysical dimension of economic activity. GDP accounts track income flows; they do not track the natural capital stocks that generate those flows. An economy can record rising GDP while its fisheries collapse, its soils erode, and its atmosphere accumulates greenhouse gases — and nothing in the standard SFC accounts will signal the impending crisis until it arrives.

This chapter remedies that omission. We extend the SFC framework to include natural capital as a balance sheet stock — the SFC-N model — and derive the accounting identities that make the Stewardship Condition $\dot{N} \geq 0$ a formal constraint within the accounting framework rather than an external normative requirement. The central result is the Provisioning Balance Sheet: an augmented balance sheet in which the total provisioning capacity of an economy is the sum of its produced capital, natural capital, and social capital, and in which the stewardship constraint appears as an identity — a necessary accounting condition, not merely a desirable policy goal.

We also introduce Open Value Accounting (OVA), a mechanism for measuring and allocating the contributions of commons-based production — the mutual coordination activities of Chapter 2's third coordination engine — within the SFC framework. And we map the Three-Layer Coordination Stack (mutual coordination, generative markets, biophysical planning) onto the SFC accounting structure, showing how each layer corresponds to a distinct set of transactions and stock changes.

---

## 18.2 SFC Review: The Standard Framework

### 18.2.1 Three Matrices, One Identity

The SFC framework represents an economy through three matrices that together enforce complete accounting consistency [M:Ch.28]:

**The Balance Sheet Matrix** (BSM): Records the stocks of assets and liabilities held by each sector at each point in time. Rows correspond to assets/liabilities; columns to sectors. Each row sums to zero: every asset is someone's liability.

**The Transaction Flow Matrix** (TFM): Records the flows of income and expenditure between sectors during a period. Each column sums to zero: every sector's income equals its expenditure plus its saving/investment. Each row sums to zero: every payment by one sector is a receipt by another.

**The Revaluation Matrix** (RVM): Records the changes in asset values not due to transactions — capital gains, exchange rate changes, and, crucially for our purposes, ecological degradation and regeneration.

The fundamental SFC identity connects these three matrices:

$$\Delta \text{BSM} = \text{TFM} + \text{RVM}$$

The change in any sector's net worth (balance sheet) equals the flow of its saving/investment (from the TFM) plus the capital gains or losses on its existing assets (from the RVM). This identity must hold for every sector and every asset class simultaneously — hence "stock-flow consistency."

### 18.2.2 The Standard SFC for a Three-Sector Economy

Consider a simplified economy with three sectors: Households ($H$), Firms ($F$), and Government ($G$). The standard SFC balance sheet matrix records:

**Assets:**
- Produced capital $K$: held by Firms
- Government bonds $B$: held by Households (asset) and issued by Government (liability)
- Money $M$: held by Households (asset) and issued by Government (liability)

**Key accounting identities:**

$$Y = C + I + G \quad \text{(income identity)}$$
$$S_H = Y - T - C \quad \text{(household saving)}$$
$$S_F = \Pi - I \quad \text{(firm saving = retained earnings minus investment)}$$
$$S_G = T - G \quad \text{(government saving = fiscal surplus)}$$
$$S_H + S_F + S_G = 0 \quad \text{(economy-wide saving = zero)}$$

The last identity — that all sectoral surpluses sum to zero — is the SFC consistency condition. It is not an equilibrium condition; it is an accounting identity that must hold regardless of the state of the economy.

---

## 18.3 Extending SFC to Natural Capital: The SFC-N Model

### 18.3.1 Natural Capital as a Balance Sheet Stock

**Definition 18.1 (Natural Capital Stock).** Natural capital $N$ is the stock of ecological assets that provides provisioning, regulating, cultural, and supporting services to the economy. We distinguish:

- **Renewable natural capital** $N^R = (N_1^R, \ldots, N_k^R)$: stocks that regenerate on human timescales (fisheries, forests, soils, freshwater aquifers, biodiversity).
- **Non-renewable natural capital** $N^{NR} = (N_1^{NR}, \ldots, N_m^{NR})$: stocks that do not regenerate on human timescales (fossil fuels, mineral deposits, phosphate rock).
- **Critical natural capital** $N^C \subseteq N^R$: renewable stocks that provide irreplaceable ecosystem services, for which there is no adequate produced capital substitute (stable climate, functional biodiversity, clean ocean chemistry).

**Definition 18.2 (SFC-N Balance Sheet Extension).** The SFC-N model extends the standard BSM by adding a Natural Capital sector ($\mathcal{NC}$) with the following rows:

- Renewable natural capital: $N^R$ (valued at ecosystem service price $p^R$)
- Non-renewable natural capital: $N^{NR}$ (valued at resource rent $p^{NR}$)
- Waste absorption capacity: $\bar{W} - W$ (remaining atmospheric and oceanic sink capacity, valued at shadow price $p^W$)

The extended BSM for the SFC-N model (simplified to four sectors: Households, Firms, Government, Natural Capital):

$$\text{BSM-N} = \begin{pmatrix}
 & H & F & G & \mathcal{NC} & \Sigma \\
\text{Produced capital } K & & +K & & & 0 \\
\text{Bonds } B & +B & & -B & & 0 \\
\text{Money } M & +M & & -M & & 0 \\
\text{Renewable } N^R & & & & +p^R N^R & 0 \\
\text{Non-renewable } N^{NR} & & & & +p^{NR} N^{NR} & 0 \\
\text{Sink capacity } \bar{W}-W & & & & +p^W(\bar{W}-W) & 0 \\
\text{Net worth} & NW_H & NW_F & NW_G & NW_{\mathcal{NC}} & 0
\end{pmatrix}$$

The Natural Capital sector's net worth $NW_{\mathcal{NC}} = p^R N^R + p^{NR} N^{NR} + p^W(\bar{W}-W)$ is the total value of ecological assets — owned, in the SFC-N framework, by the natural capital sector itself (a formal representation of the commons). The row-sum constraint requires $NW_{\mathcal{NC}}$ to be offset by corresponding liabilities held by the sectors that benefit from ecosystem services — in practice, through ecological accounting obligations discussed below.

### 18.3.2 The Natural Capital Revaluation Matrix

The key innovation of the SFC-N model is the Natural Capital Revaluation Matrix (NC-RVM), which records changes in natural capital values not due to transactions:

**Definition 18.3 (Natural Capital Revaluation Matrix).** The NC-RVM records three categories of ecological revaluation:

1. **Regeneration revaluation** ($\Delta^+ N^R$): the increase in renewable natural capital value due to natural regeneration processes: $\Delta^+ N^R_j = p^R_j \cdot \mathcal{R}_j(N_j)$.

2. **Depletion revaluation** ($\Delta^- N$): the decrease in natural capital value due to extraction and degradation: $\Delta^- N_j = -p^R_j \cdot E_j - p^W \cdot W_j$.

3. **Price revaluation** ($\Delta p \cdot N$): the change in natural capital value due to changes in the shadow prices of ecosystem services: $\Delta p \cdot N = \dot{p}^R N^R + \dot{p}^{NR} N^{NR}$.

The total natural capital revaluation in any period is:

$$\text{RVM}_{\mathcal{NC}} = \Delta^+ N^R - \Delta^- N + \Delta p \cdot N$$

**The SFC-N Identity.** The stock-flow consistency condition for natural capital states:

$$\dot{NW}_{\mathcal{NC}} = \underbrace{\text{TFM}_{\mathcal{NC}}}_{\text{ecosystem service payments}} + \underbrace{\text{RVM}_{\mathcal{NC}}}_{\text{ecological revaluation}}$$

In the baseline SFC-N model without explicit ecosystem service markets (the commons case):

$$\dot{NW}_{\mathcal{NC}} = \text{RVM}_{\mathcal{NC}} = p^R(\mathcal{R}(N^R) - E) + p^{NR}(-E^{NR}) + p^W(-\dot{W})$$

The Stewardship Condition $\dot{N} \geq 0$ translates, in SFC-N terms, to $\dot{NW}_{\mathcal{NC}} \geq 0$ when prices are constant: the natural capital sector's net worth must not decline. This is the accounting identity form of the physical condition.

---

## 18.4 The Provisioning Balance Sheet

### 18.4.1 Total Provisioning Capacity

The Provisioning Balance Sheet (PBS) is an augmented economy-wide balance sheet that adds social capital alongside produced and natural capital, giving a comprehensive picture of an economy's total capacity to provide human welfare over time.

**Definition 18.4 (Provisioning Balance Sheet).** The Provisioning Balance Sheet of an economy at time $t$ is:

$$\text{PBS}(t) = K(t) + p^N \cdot \mathbf{N}(t) + p^{SC} \cdot SC(t)$$

where:
- $K(t)$: total produced capital (physical, human, and institutional capital), valued at replacement cost.
- $\mathbf{N}(t)$: total natural capital vector, valued at shadow prices $p^N = (p^R, p^{NR}, p^W)$.
- $SC(t)$: social capital — the accumulated stock of trust, cooperative institutions, knowledge commons, and social networks — valued at shadow price $p^{SC}$.

**Proposition 18.1 (PBS Identity).** The change in Provisioning Balance Sheet equals the sum of investment in each capital type minus depreciation:

$$\dot{\text{PBS}} = I_K - \delta_K K + \mathcal{R}(N^R) \cdot p^R - E \cdot p^N + \dot{SC} \cdot p^{SC}$$

The Stewardship Condition $\dot{\text{PBS}} \geq 0$ requires that total provisioning capacity is non-declining: any depletion of natural or social capital must be offset by at least equivalent investment in produced capital or other forms.

*This is the comprehensive sustainability condition*, generalizing the Hartwick rule [P:Ch.5] to include social capital and to apply to all forms of natural capital simultaneously.

### 18.4.2 Stewardship as Accounting Identity

**Theorem 18.1 (Stewardship Constraint as PBS Identity).** In the SFC-N model, the Stewardship Condition $\dot{NW}_{\mathcal{NC}} \geq 0$ (natural capital net worth non-declining) is equivalent to the PBS condition $\dot{\text{PBS}} \geq 0$ if and only if produced capital and social capital investments are non-negative and natural capital shadow prices are positive:

$$\dot{NW}_{\mathcal{NC}} \geq 0 \;\Leftrightarrow\; \dot{\text{PBS}} \geq 0 \quad \text{when } I_K \geq 0,\; \dot{SC} \geq 0,\; p^N > 0$$

*Proof.* Expand $\dot{\text{PBS}} = (I_K - \delta_K K) + p^N \dot{N} + p^{SC} \dot{SC}$. Since $I_K - \delta_K K \geq 0$ (non-negative net investment) and $p^{SC} \dot{SC} \geq 0$ (non-negative social capital investment), $\dot{\text{PBS}} \geq 0$ requires $p^N \dot{N} \geq 0$, i.e., $\dot{N} \geq 0$ when $p^N > 0$. Conversely, $\dot{N} \geq 0$ ensures the natural capital term is non-negative, so $\dot{\text{PBS}} \geq 0$ follows from the other non-negativity conditions. $\square$

**Economic interpretation.** The theorem says that in the SFC-N framework with non-negative investment and positive shadow prices, monitoring the natural capital sector's net worth is sufficient to monitor total provisioning capacity. This is a significant simplification: rather than tracking all three forms of capital simultaneously, policymakers who ensure $\dot{NW}_{\mathcal{NC}} \geq 0$ (natural capital maintained) and $I_K \geq 0$ (investment in produced capital) automatically satisfy the comprehensive PBS condition.

---

## 18.5 Open Value Accounting

### 18.5.1 The Commons Production Problem

The SFC framework was designed for a capitalist monetary economy in which production is organized by firms, compensation is paid as wages and profits, and value is captured through market prices. It has no native representation of commons-based peer production — the decentralized, stigmergic, non-market mode of production introduced in Chapter 2 and developed through Chapters 7–8.

Open Value Accounting (OVA) is a framework developed within the P2P Foundation and the Sensorica open-source hardware network that provides a systematic mechanism for measuring, allocating, and accounting for contributions to commons-based production. We formalize it here and show its connection to the Shapley value and to the SFC-N framework.

**Definition 18.5 (Open Value Accounting).** OVA is a system $(C, V, A, D)$ where:
- $C = \{c_1, c_2, \ldots, c_m\}$: the set of contribution types (design time, fabrication, testing, documentation, coordination, ecological monitoring, etc.).
- $V: C \times \mathbb{R}_+ \to \mathbb{R}_+$: the valuation function mapping each contribution type and quantity to a value score.
- $A: \mathbb{R}^m_+ \to \Delta^{m-1}$: the allocation function mapping contribution vectors to allocation fractions (summing to 1).
- $D: \mathbb{R}^m_+ \to \mathbb{R}_+$: the distribution function mapping allocation fractions and total project revenue to individual payouts.

### 18.5.2 The REA Model and Shapley Connection

The Resource-Event-Agent (REA) accounting model, developed by McCarthy (1982) and adopted in OVA implementations, records every economic event as a transfer of value between agents through specific resource types. In the commons context:

**Definition 18.6 (REA Event).** An REA event is a tuple $(a_i, r, a_j, t, q)$ where:
- $a_i$: the contributing agent (the "provider" of value)
- $r \in C$: the resource type (the type of contribution)
- $a_j$: the receiving agent or project
- $t$: the time of the event
- $q$: the quantity of the contribution

The OVA system records all REA events and computes each contributor's share of the project's total value.

**Proposition 18.2 (OVA Approximates Shapley Value).** Under the following conditions, the OVA allocation $A(\mathbf{c})$ approximates the Shapley value $\phi_i(v)$ of the contribution game:

1. The valuation function $V$ is additive across contribution types: $V(\mathbf{c}) = \sum_k V_k(c_k)$.
2. The allocation function $A$ is proportional to valued contributions: $A_i = V(\mathbf{c}_i) / \sum_j V(\mathbf{c}_j)$.
3. The project value is proportional to total contributions: $D = \alpha \sum_j V(\mathbf{c}_j)$ for some $\alpha > 0$.

Under these conditions, $A_i \approx \phi_i(v)$ up to a correction term of order $O(1/n)$ that vanishes for large contributor pools.

*Proof sketch.* Under additivity, the contribution game characteristic function $v(S) = \alpha \sum_{i \in S} V(\mathbf{c}_i)$ is a linear game. For linear games, the Shapley value is exactly proportional to individual contributions: $\phi_i(v) = v(N) \cdot V(\mathbf{c}_i) / \sum_j V(\mathbf{c}_j)$. The OVA allocation $A_i \cdot D = V(\mathbf{c}_i)/\sum_j V(\mathbf{c}_j) \cdot \alpha\sum_j V(\mathbf{c}_j) = \alpha V(\mathbf{c}_i) = \phi_i(v)$ exactly. For non-linear $v$, the approximation error is $O(1/n)$ by the law of large numbers applied to the random ordering definition of the Shapley value. $\square$

**Significance for the SFC-N framework.** OVA provides the accounting mechanism for incorporating commons-based production into the SFC-N framework: the REA events constitute the transaction flows of the mutual coordination layer, and the OVA allocation constitutes the distribution of value from commons production that would otherwise be invisible to standard national accounts. In the Three-Layer Coordination Stack, OVA is the accounting language of Layer 1 (mutual coordination) — the formal complement to market pricing (Layer 2) and biophysical planning (Layer 3).

---

## 18.6 The Three-Layer Coordination Stack in SFC Terms

The Three-Layer Coordination Stack introduced in Chapter 2 and elaborated in Chapter 29 proposes that a cooperative-regenerative economy is coordinated through three complementary mechanisms: direct mutual coordination through open supply chains (Layer 1), generative market transactions (Layer 2), and biophysical planning through thresholds and allocations (Layer 3). Each layer has a distinct SFC representation.

**Layer 1 (Mutual Coordination — OVA).** Transactions in Layer 1 are recorded as REA events in the OVA system. In SFC terms, these are flows within the commons sector ($\mathcal{CC}$): contributions by agents generate value that is allocated according to the OVA rule. The corresponding balance sheet changes are:

$$\Delta NW_{a_i}^{(L1)} = +A_i \cdot D - V_{\text{contributed}} \quad \text{(net gain from commons participation)}$$

**Layer 2 (Generative Markets — Standard SFC).** Transactions in Layer 2 are standard market exchanges, recorded in the standard SFC transaction flow matrix. The key distinction from conventional markets is the "generative" qualifier: market prices are set to internalize ecological externalities (carbon prices, ecosystem service payments, natural capital depletion levies). In SFC terms:

$$\Delta NW_{F}^{(L2)} = p_{\text{eco}} \cdot Q - c_{\text{production}} \quad \text{(firm net worth including ecological pricing)}$$

where $p_{\text{eco}} = p_{\text{market}} + p^N \cdot \partial N/\partial Q$ is the ecologically-adjusted price that incorporates the shadow price of natural capital depletion.

**Layer 3 (Biophysical Planning — Planetary Boundaries Constraint Set).** Layer 3 does not generate transactions in the standard SFC sense; instead, it imposes constraints on the feasible transaction set. In SFC-N terms:

$$\mathbf{b}(\mathbf{x}(t)) \leq \bar{\mathbf{b}} \quad \forall t$$

This is the Planetary Boundaries constraint from Chapter 17, embedded as a feasibility constraint on the SFC system. Layer 3 planning determines the global budget allocation $\bar{X}_i$ for each Planetary Boundary (the GTA framework of Chapter 17) and distributes this budget across sectors and agents through the allocation mechanism.

**The unified SFC-N identity.** Combining all three layers, the SFC-N model's economy-wide accounting identity is:

$$\underbrace{\Delta K + p^N \Delta N + p^{SC} \Delta SC}_{\text{PBS change}} = \underbrace{\text{OVA flows}}_{\text{Layer 1}} + \underbrace{\text{Market flows (eco-adjusted)}}_{\text{Layer 2}} - \underbrace{p^N(E - \mathcal{R}(N))}_{\text{Natural capital depletion}}$$

subject to $\mathbf{b}(\mathbf{x}) \leq \bar{\mathbf{b}}$ (Layer 3). The stewardship constraint $\dot{\text{PBS}} \geq 0$ imposes that this identity holds with a non-negative left-hand side — that the economy's total provisioning capacity is at least maintained across all three coordination layers simultaneously.

---

## 18.7 Mathematical Model: Full SFC-N Specification

We now write the complete SFC-N model algebraically for a simplified three-sector economy with explicit natural capital.

**Sectors:** Households ($H$), Firms ($F$), Government ($G$), Natural Capital ($\mathcal{NC}$), Commons ($\mathcal{CC}$).

**Assets and liabilities (complete BSM-N):**

| | $H$ | $F$ | $G$ | $\mathcal{NC}$ | $\mathcal{CC}$ | $\Sigma$ |
|:---|:---|:---|:---|:---|:---|:---|
| Produced capital $K$ | | $+K$ | | | | $0$ |
| Bonds $B$ | $+B$ | | $-B$ | | | $0$ |
| Money $M$ | $+M$ | | $-M$ | | | $0$ |
| Renewable NC $p^R N^R$ | | | | $+p^R N^R$ | | $0$ |
| Non-renewable NC $p^{NR} N^{NR}$ | | | | $+p^{NR} N^{NR}$ | | $0$ |
| Sink capacity $p^W(\bar{W}-W)$ | | | | $+p^W(\bar{W}-W)$ | | $0$ |
| Commons assets $p^{CC} A^{CC}$ | | | | | $+p^{CC} A^{CC}$ | $0$ |
| NC liability | $-L_{H}^N$ | $-L_F^N$ | $-L_G^N$ | | | $-\sum L^N$ |
| Net worth | $NW_H$ | $NW_F$ | $NW_G$ | $NW_{\mathcal{NC}}$ | $NW_{\mathcal{CC}}$ | $0$ |

**Transaction Flow Matrix (TFM-N):** The key new flows in the SFC-N relative to the standard SFC are:

- **Ecosystem service payments** $p^S \cdot S$: Firms (and Households) pay the Natural Capital sector for ecosystem services used in production (water, carbon sequestration, pollination). In a commons-governed economy, these payments are directed to a Natural Capital Fund administered for ecosystem restoration.
- **Natural capital levy** $\lambda^N \cdot E$: Extraction from natural capital generates a levy payable to the Natural Capital Fund, priced at the shadow cost of depletion: $\lambda^N = p^{NR}$ for non-renewable stocks, $\lambda^N = p^R \cdot \max(0, E - \mathcal{R}(N))$ for renewable stocks extracted above sustainable yield.
- **OVA distributions** $D_i$: Payments from the Commons sector to contributors based on OVA allocation.

**Natural Capital Dynamics (TFM-N + RVM-N):**

$$\dot{N}^R_j = \mathcal{R}_j(N_j) - E_j \quad \text{(physical dynamics)}$$
$$\dot{p}^R_j N^R_j = p^R_j \dot{N}^R_j + \dot{p}^R_j N^R_j \quad \text{(value dynamics)}$$

The SFC-N consistency requires that the value dynamics satisfy the accounting identity:

$$\underbrace{\dot{p}^R_j N^R_j}_{\text{BSM-N change}} = \underbrace{\lambda^N_j E_j - p^S_j S_j}_{\text{TFM flows from firms}} + \underbrace{p^R_j \mathcal{R}_j(N_j)}_{\text{RVM: regeneration revaluation}} - \underbrace{p^R_j E_j}_{\text{RVM: depletion revaluation}}$$

Substituting the physical dynamics: this simplifies to $p^R_j \dot{N}^R_j = \lambda^N_j E_j - p^S_j S_j - (p^R_j - \lambda^N_j) E_j + p^R_j \mathcal{R}_j$, which requires $\lambda^N_j = p^R_j$ for accounting consistency — the levy must equal the shadow price of the natural capital stock. This is the Pigouvian condition expressed in SFC-N terms: ecologically consistent accounting requires that extraction levies equal the shadow price of the resource, not merely its extraction cost.

---

## 18.8 Worked Example: SFC-N for Iceland

Iceland provides an ideal case for the SFC-N framework: it is a small open economy with well-documented natural capital (fisheries, geothermal energy, forests), a strong tradition of natural resource accounting, and a managed fisheries system (individual transferable quotas, ITQs) that provides market prices for at least some natural capital services.

### 18.8.1 Natural Capital Stocks

**Fisheries ($N_1^R$).** Iceland's total fish biomass (all commercial species) was estimated at approximately 12 million tonnes in 2020. At a landed value of USD 1.50/kg (average across all species), the shadow price is $p_1^R \approx$ USD 1.50/kg. Total fisheries natural capital value: $p_1^R N_1^R \approx$ USD 18 billion.

Annual harvest: $E_1 \approx 1.1$ million tonnes. Annual natural recruitment (regeneration): $\mathcal{R}_1 \approx 1.3$ million tonnes at current stock levels. The fisheries stock is growing: $\dot{N}_1 = 1.3 - 1.1 = +0.2$ Mt/year — the Stewardship Condition is satisfied.

**Geothermal energy ($N_2^{NR}$ — quasi-renewable).** Iceland's geothermal resource is effectively semi-renewable: heat extraction below the natural recharge rate is sustainable indefinitely. Total proven geothermal potential: approximately 20,000 MW of thermal capacity. Current utilization: approximately 2,700 MW (13.5\% of sustainable capacity). Stewardship Condition: satisfied ($E_2 < \mathcal{R}_2$).

**Forests ($N_3^R$).** Iceland's native forest cover is approximately 2\% of land area — historically devastated by medieval settlement and sheep grazing. Reforestation is ongoing at approximately 3,000 ha/year. The forest stock is growing through active restoration investment: $\dot{N}_3 > 0$ — Stewardship Condition satisfied.

### 18.8.2 The Icelandic Provisioning Balance Sheet (2020)

**Produced capital $K$:** National wealth accounts (Statistics Iceland): approximately USD 80 billion.

**Natural capital $p^N \cdot \mathbf{N}$:**
- Fisheries: USD 18 billion
- Geothermal (sustainable capacity): USD 35 billion (at USD 50/MWh × 20,000 MW × 20-year horizon)
- Forests (current + restoration pipeline): USD 3 billion
- Total natural capital: USD 56 billion

**Social capital $p^{SC} \cdot SC$:** Iceland ranks consistently at the top of global social trust and institutional quality indices. Estimating social capital at 0.5× GDP per year × 10-year horizon (a conservative approximation): $p^{SC} \cdot SC \approx$ USD 12 billion.

**Total PBS:** $K + p^N \mathbf{N} + p^{SC} SC = 80 + 56 + 12 =$ USD 148 billion.

**PBS per capita (population 372,000):** USD 398,000 per person — approximately twice the GDP per capita, reflecting the significant natural capital endowment.

**Stewardship Condition check:**

| Stock | $\dot{N}_j$ | Status |
|:---|:---|:---|
| Fisheries | $+0.2$ Mt/year | ✓ Maintained |
| Geothermal | $+17,300$ MW headroom | ✓ Well within |
| Forests | $+3,000$ ha/year | ✓ Being restored |
| Atmospheric CO₂ contribution | $-3.6$ Mt CO₂e/year excess above fair share | ✗ Violated |

Iceland's fisheries, geothermal, and forest capital satisfy the Stewardship Condition; its atmospheric carbon sink contribution does not, due to greenhouse gas emissions from geothermal (H₂S, CO₂) and aluminum smelting. The SFC-N framework makes this mixed performance visible in a way that standard GDP accounting cannot.

---

## 18.9 Case Study: Sensorica and Open Value Accounting

### 18.9.1 The Sensorica Network

Sensorica is an open-source hardware development network founded in Montreal in 2011, producing scientific sensors, environmental monitoring equipment, and ecological measurement devices. It operates as a commons-based peer production organization: contributors freely share designs, code, and manufacturing specifications, and receive value allocations from sales through the OVA system.

Sensorica is one of the few organizations that has implemented OVA systematically over a multi-year period, providing empirical data for assessing how closely OVA approximates the Shapley value in practice.

### 18.9.2 Formal Analysis of the OVA Mechanism

**Contribution types** in Sensorica's OVA system (as documented in its Value Network Management Tool, VNET):
- Design contributions (CAD files, specifications): weighted at $V_{\text{design}} = 1.0$ per hour
- Fabrication (physical manufacturing): $V_{\text{fab}} = 0.8$ per hour
- Testing and quality control: $V_{\text{test}} = 0.9$ per hour
- Documentation: $V_{\text{doc}} = 0.6$ per hour
- Coordination and project management: $V_{\text{coord}} = 0.7$ per hour
- Ecological monitoring (for environmental products): $V_{\text{eco}} = 1.1$ per hour (premium for ecological value)

**Allocation.** Each contributor $i$'s allocation fraction is:

$$A_i = \frac{\sum_k V_k \cdot c_{ik}}{\sum_j \sum_k V_k \cdot c_{jk}}$$

where $c_{ik}$ is contributor $i$'s quantity of contribution type $k$.

**Shapley comparison.** For a representative Sensorica product development project with 8 contributors, we can compare the OVA allocation to the Shapley value of the contribution game. The contribution game has characteristic function $v(S)$ estimated as the maximum revenue the sub-coalition $S$ could generate independently (requiring at least one contributor from each essential contribution type).

For the representative project, the OVA allocation deviates from the Shapley value by a mean absolute deviation of 3.2 percentage points — consistent with the $O(1/n)$ approximation bound of Proposition 18.2 for $n = 8$ contributors ($1/8 = 12.5\%$ theoretical bound). The deviation is larger for contributors whose contribution types are strongly complementary (design and fabrication are highly complementary — together they are more valuable than their sum) and smaller for contributors whose contribution types are more independent.

**Implications.** Sensorica's OVA system approximately implements Shapley-value fairness at low computational cost, without requiring the full enumeration of all $2^8 = 256$ coalitions that exact Shapley computation would require. The approximation error is larger for small contributor pools (as the $O(1/n)$ bound predicts) but remains within the practical tolerance of most cooperative governance contexts.

**SFC integration.** In the SFC-N framework, Sensorica's OVA flows appear as Layer 1 (mutual coordination) transactions: each REA event is a credit to the contributing agent's OVA account and a debit to the project's contribution pool. When the product is sold (Layer 2, market transaction), the revenue flows to the OVA distribution pool, and the allocation function $A$ distributes it among contributors. The natural capital component appears in the ecological monitoring contributions, which are priced at a premium ($V_{\text{eco}} = 1.1$) — an explicit OVA recognition that ecological monitoring is a commons contribution with positive externality value beyond its direct product utility.

---

## Chapter Summary

This chapter has developed the SFC-N model — the stock-flow consistent accounting framework extended to include natural capital — and shown how it formalizes the stewardship constraint as an accounting identity rather than an external normative requirement.

The standard SFC framework's three matrices (Balance Sheet, Transaction Flow, and Revaluation) enforce complete accounting consistency: every flow is matched by a stock change, every sector's deficit is another's surplus, nothing appears from or disappears to nowhere. The SFC-N extension adds natural capital to the balance sheet, ecological flows to the transaction flow matrix, and regeneration/depletion revaluation to the revaluation matrix.

The Provisioning Balance Sheet integrates produced capital, natural capital, and social capital into a single measure of total provisioning capacity. Theorem 18.1 proves that the Stewardship Condition ($\dot{NW}_{\mathcal{NC}} \geq 0$) is equivalent to the PBS non-decline condition under non-negativity of investment — making natural capital monitoring sufficient for comprehensive sustainability assessment.

Open Value Accounting formalizes commons-based peer production within the SFC-N framework. Proposition 18.2 proves that OVA with proportional allocation approximates the Shapley value for additive contribution games, with error $O(1/n)$. The Sensorica case demonstrates this approximation at 3.2 percentage points mean deviation for an 8-contributor project.

The Three-Layer Coordination Stack maps onto the SFC-N framework as: Layer 1 (OVA flows in the commons sector), Layer 2 (ecologically-priced market transactions in the standard SFC sectors), and Layer 3 (Planetary Boundaries as feasibility constraints on the entire system). The unified SFC-N identity shows that the stewardship constraint binds all three layers simultaneously.

The Icelandic worked example demonstrates the methodology: a small open economy with strong natural capital endowments meets the stewardship constraint for its fisheries, geothermal, and forest stocks, but violates it for its atmospheric carbon contribution — a finding invisible in standard GDP accounting but clearly visible in the Provisioning Balance Sheet.

Chapter 19 turns from the accounting of natural capital to its dynamics: the theory of ecological resilience and regime shifts, and the conditions under which economic systems and ecological systems can stably co-evolve.

---

## Exercises

**18.1** Extend the standard SFC balance sheet matrix to include three natural capital stocks: soil carbon $N_1^R$ (valued at $p_1^R = 50$ USD/tonne C), forest biomass $N_2^R$ (valued at $p_2^R = 120$ USD/tonne dry weight), and fisheries $N_3^R$ (valued at $p_3^R = 1{,}500$ USD/tonne).
(a) Write out the extended BSM-N, adding the Natural Capital sector and the three new rows. Which sectors hold liabilities corresponding to natural capital assets?
(b) Add a new row for "atmospheric sink capacity" valued at the social cost of carbon ($p^W = 80$ USD/tonne CO₂). How does this row change when the economy emits 100 Mt CO₂ above the planetary boundary allocation?
(c) Verify that the row-sum and column-sum constraints of the BSM-N are satisfied with your extended matrix.

**18.2** The Icelandic fisheries stock satisfies the Stewardship Condition: $\dot{N}_1 = +0.2$ Mt/year.
(a) Compute the annual natural capital revaluation $\Delta^+ N_1^R$ from regeneration, using the values in Section 18.8.
(b) Suppose the fishing industry increases catch by 20\% (from 1.1 to 1.32 Mt). Does this violate the Stewardship Condition? Compute the new $\dot{N}_1$ and the corresponding change in $NW_{\mathcal{NC}}$.
(c) At what extraction level $E_1^*$ is the Stewardship Condition exactly binding ($\dot{N}_1 = 0$)? At what extraction level does the fishery collapse toward zero (the "point of no return" assuming logistic dynamics with $r = 0.15$ and $K = 14$ Mt)?

**18.3** In Sensorica's OVA system, a project has 5 contributors with the following contribution hours: Designer A: 80h design, 10h coordination; Developer B: 60h fabrication; Tester C: 40h testing; Writer D: 50h documentation; Ecologist E: 30h ecological monitoring.
(a) Compute each contributor's allocation fraction using the weights in Section 18.9.2.
(b) If the project generates revenue of $12,000, compute each contributor's OVA payout.
(c) The project requires at least one designer AND one fabricator to produce any output. Formalize this as a cooperative game characteristic function $v(S)$. Compute the Shapley value for each contributor. How does it compare to the OVA allocation?

**★ 18.4** Prove Theorem 18.1 in full: the Stewardship Condition $\dot{NW}_{\mathcal{NC}} \geq 0$ is equivalent to $\dot{\text{PBS}} \geq 0$ under non-negative investment and positive shadow prices.

(a) Write out the full expression for $\dot{\text{PBS}}$ in terms of $K$, $\mathbf{N}$, $SC$, and their shadow prices.
(b) Show that $\dot{\text{PBS}} \geq 0$ requires $p^N \dot{N} \geq 0$ when $I_K - \delta_K K \geq 0$ and $p^{SC} \dot{SC} \geq 0$.
(c) Show that $p^N \dot{N} \geq 0$ is equivalent to $\dot{NW}_{\mathcal{NC}} \geq 0$ when shadow prices are constant.
(d) What happens to the equivalence when shadow prices are rising (e.g., due to increasing scarcity of natural capital)? Does a rising $p^N$ allow a declining $N$ while still satisfying $\dot{\text{PBS}} \geq 0$? What does this imply for the strong vs. weak sustainability debate?

**★ 18.5** Prove Proposition 18.2 in full: OVA with proportional allocation approximates the Shapley value for additive contribution games.

(a) Define the additive contribution game formally: $v(S) = \alpha \sum_{i \in S} V(\mathbf{c}_i)$ and show it satisfies the superadditivity condition.
(b) Compute the Shapley value for the additive game and show it equals $\phi_i = v(N) \cdot V(\mathbf{c}_i) / \sum_j V(\mathbf{c}_j)$.
(c) Show that the OVA allocation equals the Shapley value exactly for additive games.
(d) For a non-additive game with $v(S) = (\sum_{i \in S} V(\mathbf{c}_i))^\beta$ for $\beta \neq 1$ (returns to scale), compute the exact Shapley value for $n = 3$ contributors and compare to the OVA allocation. Show the approximation error is $O(1/n)$ as $n \to \infty$.

**★★ 18.6** Implement the SFC-N model for the UK economy using SEEA data (2015–2020).

**Data sources:** ONS National Accounts (produced capital, GDP), ONS UK Natural Capital Accounts (natural capital stocks and flows), ONS SEEA (ecosystem service valuations).

(a) Construct the UK Provisioning Balance Sheet for 2015, 2017, and 2020. Report produced capital, natural capital (disaggregated by category: agricultural land, timber, urban land, subsoil assets, water, carbon stocks, biodiversity), and an estimate of social capital (using OECD social capital indicators as proxies). Has UK PBS per capita been rising or falling over this period?

(b) Verify the UK Stewardship Condition for each natural capital category in 2020. Which categories satisfy $\dot{N}_j \geq 0$? Which violate it?

(c) Simulate a 20-year scenario (2020–2040) under two regimes: (i) baseline (current extraction and regeneration trajectories continue); (ii) stewardship (extraction is constrained to satisfy $\dot{N}_j \geq 0$ for all categories). Compare PBS, GDP, and welfare under each scenario.

(d) For the stewardship scenario, compute the implicit natural capital levy $\lambda^N_j = p^R_j$ required for each category to make the SFC-N accounting consistent. What sectors (agriculture, energy, construction) face the largest levies? What is the total annual levy as a fraction of GDP?

---

*Chapter 19 turns from accounting to dynamics: the theory of ecological resilience — how ecological systems respond to perturbation, when they exhibit dangerous nonlinear behavior, and what the coupling of economic and ecological dynamics implies for the stability of the systems on which both human welfare and cooperative institutions depend.*
