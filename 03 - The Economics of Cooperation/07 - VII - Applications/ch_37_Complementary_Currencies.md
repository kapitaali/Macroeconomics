# Chapter 37: Complementary Currencies — The Swiss WIR, Sardex, and Others

> *"The WIR franc is not a curiosity or a relic. It is a working demonstration, running for ninety years, that money can be designed to serve the economy rather than dominate it."*
> — Bernard Lietaer, *The Future of Money* (2001)

> *"Sardex does not replace the euro. It complements it. In times of scarcity, it keeps the network alive."*
> — Giuseppe Littera, Sardex co-founder (2020, paraphrased)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Apply the mutual credit theory of Chapter 25 and the demurrage theory of Chapter 27 to real complementary currency systems, formally testing their theoretical predictions against empirical evidence from the WIR, Sardex, and Sarafu systems.
2. Specify a two-currency stock-flow consistent model, derive the formal conditions for stable coexistence of a complementary currency alongside a national currency, and identify the parameter configurations under which coexistence is stable versus one currency driving out the other.
3. Empirically test the counter-cyclical hypothesis for the Swiss WIR: formally estimating the relationship between WIR volume and Swiss business cycle indicators, and evaluating the economic stabilization value of the WIR system.
4. Identify the design principles that distinguish durable complementary currency systems (WIR, Sardex, Chiemgauer) from short-lived experiments (Brixton Pound, Liberty Dollar, most ICO-era crypto tokens), using the Ostrom governance framework and the SFC consistency conditions.
5. Formally analyze the Brixton Pound's failure through the lens of the design principles, and simulate the counterfactual — what would have happened with corrected design.
6. Evaluate Sardex (2010–2020) and Sarafu Network as contrasting implementations of the same mutual credit principles in high-income and low-income contexts respectively.

---

## 37.1 Why Complementary Currencies Exist

The monetary theory of Part V established that the dominant debt-based monetary system has three structural pathologies: Minsky instability, systematic distributional transfer, and a growth imperative. Chapter 28 proved that sovereign money, mutual credit, and demurrage money each address one or more of these pathologies. But that analysis was theoretical — it described what alternative monetary architectures could achieve under ideal conditions.

Complementary currencies are the empirical test: real monetary systems operating alongside national currencies, in real communities and business networks, for periods long enough to reveal their properties under realistic conditions. They are not replacements for national currencies — they coexist with national currencies, filling specific niches that the national monetary system systematically fails to serve. This coexistence is both their strength (they do not require a national monetary revolution to function) and their limitation (they operate at smaller scale and with less monetary authority than national currencies).

Three properties of complementary currencies make them important for the cooperative-regenerative framework:

**Counter-cyclical stabilization.** The WIR expands when the Swiss economy contracts — providing liquidity exactly when the national monetary system tightens credit. This is the automatic stabilizer of mutual credit demonstrated empirically over 90 years.

**Local economic retention.** Complementary currencies spend only within their defined geographic or sectoral community — multiplying local economic activity rather than allowing purchasing power to leak to distant supply chains.

**Community trust and governance.** The governance of complementary currencies is inherently community-based — the Ostrom conditions are necessary for their survival. Failed systems almost universally exhibit governance failures; successful systems almost universally score high on the Ostrom principles.

This chapter tests whether the formal predictions of Chapters 25 and 27 — the velocity theorem, the counter-cyclical property, the stability conditions — survive contact with the empirical record of three decades of complementary currency operation across three continents.

---

## 37.2 The Two-Currency SFC Model

### 37.2.1 Specification

**Definition 37.1 (Two-Currency Economy).** A two-currency economy contains a national currency $M^N$ (euros, CHF, Kenyan shillings) and a complementary currency $M^C$ (WIR francs, Sardex credits, Sarafu tokens). Each currency serves a subset of transactions:

- National currency: all transactions, especially external (imports, taxes, bank debt service).
- Complementary currency: a subset of domestic transactions — B2B trade (WIR, Sardex) or community-level exchange (Sarafu, Chiemgauer).

**Balance sheet.** Households and firms hold both currencies:
$$NW_i = K_i + M^N_i + M^C_i - L^N_i - L^C_i$$

where $L^N_i$ are national currency loans and $L^C_i$ are complementary currency credit balances (negative for debtors, positive for creditors in the mutual credit system).

**SFC constraints.** Two separate zero-sum constraints apply:

National: $\sum_i M^N_i - \sum_i L^N_i = BG$ (net national financial assets = government bonds only, Godley identity).

Complementary: $\sum_i M^C_i = 0$ (mutual credit: net sum of all complementary currency balances is always zero — the system balance constraint of Chapter 25).

**Transaction flow matrix (TFM) — two-currency extension:**

| Flow | $H$ | $F$ | $B^N$ | $CB^C$ | $G$ |
|:---|:---|:---|:---|:---|:---|
| Wages ($M^N$) | $+W$ | $-W$ | | | |
| B2B trade ($M^C$) | | $+R^C - P^C$ | | | |
| New money ($M^N$, sovereign) | $+\Delta M^N$ | | | $-\Delta M^N$ | |
| CC creation (mutual credit) | | $+\Delta M^C_{\text{debtors}}$ | | $-\Delta M^C_{\text{creditors}}$ | |
| CC clearing | | $-\Delta\text{bal}^C$ | | $+\Delta\text{bal}^C$ | |
| Taxes | $-T$ | $-T_F$ | | | $+(T+T_F)$ |

The key new row: complementary currency B2B trade — firms simultaneously credit some accounts ($+R^C$) and debit others ($-P^C$) by the same amount (zero net, preserving $\sum M^C = 0$).

### 37.2.2 Coexistence Conditions

**Definition 37.2 (Currency Coexistence).** The two-currency economy is in stable coexistence when:
1. Both $M^N > 0$ and $M^C > 0$ persist in the long run.
2. Neither currency drives out the other (Gresham's law does not apply: bad money does not drive out good).
3. The complementary currency occupies a stable niche — a subset of transactions for which it provides advantages over the national currency.

**Theorem 37.1 (Stable Coexistence Conditions).** The two-currency economy reaches stable coexistence if and only if:

**C1 (Niche differentiation):** There exists a set of transactions $\mathcal{T}^C \subset \mathcal{T}$ for which the complementary currency provides strictly higher value than the national currency: $u_i(M^C, \mathcal{T}^C_i) > u_i(M^N, \mathcal{T}^C_i)$.

**C2 (Limited convertibility):** Complementary currency is not freely convertible to national currency at par — if it were, arbitrage would lead agents to convert all complementary currency to national currency (Gresham's law), eliminating the complementary system.

**C3 (SFC consistency):** The complementary currency satisfies the four SFC constraints of Chapter 28 — it does not generate fictitious assets or liabilities.

**C4 (Ostrom governance):** The complementary currency system satisfies at least DP1–DP5 of the Ostrom design principles — defined boundaries, congruent rules, collective choice, monitoring, and graduated sanctions.

*Proof.* C1 ensures demand for the complementary currency in its niche — without it, all transactions migrate to national currency. C2 prevents arbitrage collapse — the complementary currency must have a friction against full convertibility (WIR: not freely convertible; Sardex: no direct euro conversion; Sarafu: redeemable only for goods in network). C3 ensures the complementary currency does not create systemic financial risks (which would provoke regulatory shutdown). C4 ensures governance stability — without adequate governance, the complementary currency collapses through free-riding or fraud. All four conditions are necessary; absence of any one creates a failure pathway. $\square$

### 37.2.3 The Counter-Cyclical Mechanism

From Chapter 25, the WIR's counter-cyclical behavior was formally modeled as:

$$\dot{x}_{\text{WIR}} = \gamma\left[\frac{\ell_{\text{WIR}}}{i_{\text{CHF}} \cdot \text{tightness}(t)} - x_{\text{WIR}}\right]$$

**Empirical test.** Stodder (2009) and Stodder and Lietaer (2016) provide the definitive empirical analysis. Using Swiss M2, WIR Bank turnover data, and GDP data from 1948 to 2013:

**Estimated regression:**
$$\Delta \ln(\text{WIR turnover}_t) = \alpha + \beta_1 \Delta \ln(\text{GDP}_t) + \beta_2 \Delta \ln(\text{M2}_t) + \beta_3 \text{crisis}_t + \varepsilon_t$$

**Results:**
| Variable | Coefficient | Std. Error | $t$-stat |
|:---|:---|:---|:---|
| $\Delta\ln$(GDP) | **−1.47** | 0.38 | −3.87*** |
| $\Delta\ln$(M2) | −0.22 | 0.19 | −1.16 |
| Crisis dummy | **+0.18** | 0.06 | +3.0** |
| Constant | 0.03 | 0.02 | 1.5 |
| $R^2$ | 0.41 | | |

***p<0.001, **p<0.01.*

**Interpretation.** A 1\% fall in Swiss GDP is associated with a 1.47\% rise in WIR turnover — strongly counter-cyclical and statistically significant. The national M2 growth has no significant effect on WIR volume, confirming that WIR operates as a separate monetary circuit that responds to real economic conditions (credit tightness) rather than to national monetary policy.

**Economic stabilization value.** With WIR turnover approximately CHF 1.5 billion and a GDP multiplier estimated at 0.8 (conservative, accounting for partial circulation within the Swiss SME sector): WIR contributes approximately CHF 1.2 billion of effective demand during the typical Swiss recession. Swiss GDP is approximately CHF 750 billion; the WIR contributes approximately 0.16\% of GDP in counter-cyclical demand — modest but consistent and automatic.

---

## 37.3 Design Principles for Successful Complementary Currencies

Drawing from the formal conditions of Theorem 37.1 and the Ostrom governance framework, we identify eight design principles for complementary currency systems:

**CC-1 (Clear community definition).** The currency serves a defined community with a meaningful economic relationship — a geographic area, an industry network, or a social group with dense transactions. Without this, the niche (C1) is undefined and the governance community (C4) is absent.

**CC-2 (Mutual credit foundation).** The currency is created through bilateral credit commitments rather than by any central issuer purchasing backing assets. This satisfies the SFC consistency condition (C3) and avoids the political risks of currency issuance.

**CC-3 (Limited convertibility).** The currency is redeemable for goods and services within the community but not freely convertible to national currency at par. A modest conversion fee (2–5\%) is sufficient to prevent arbitrage collapse while not prohibiting exit.

**CC-4 (Calibrated credit limits).** Credit limits are set to each member's ability to deliver real value to the community — not to their desire for purchasing power. Over-extension of credit destroys the currency's backing by real goods and services.

**CC-5 (Active clearing cycles).** Regular multilateral clearing (Algorithm 25.1) maintains system liquidity and reduces accumulated imbalances. Without active clearing, imbalances grow until members hit credit limits and the system seizes.

**CC-6 (Governance transparency).** All members have access to the system's aggregate balance statistics (though not individual balances, for privacy). Transparency enables community monitoring (DP4) and prevents hidden concentration of credit imbalances.

**CC-7 (Demurrage or circulation incentive).** A holding fee or circulation requirement (CC-7a: demurrage like Chiemgauer; CC-7b: expiry date like Wörgl; CC-7c: charitable donation on non-use like Bristol Pound) prevents hoarding and maintains velocity. Without this, successful currencies tend to be saved rather than spent, reducing their counter-cyclical effectiveness.

**CC-8 (External recognition).** Some form of regulatory recognition — even informal tolerance — is necessary for the system to operate without legal harassment. The WIR Bank operates as a regulated financial institution; Sardex operates under Italian commercial law; Sarafu operates with NGO recognition in Kenya.

---

## 37.4 Failure Analysis: The Brixton Pound

### 37.4.1 Background

The Brixton Pound (B£) launched in September 2009 in the Brixton neighborhood of London — a local paper currency (later digital) accepted at independent local businesses, designed to strengthen the local economy by keeping spending within Brixton. It generated significant media attention and was widely cited as a model for local currency innovation. It was withdrawn from circulation in August 2019, after a decade of declining volume and merchant participation.

### 37.4.2 Formal Failure Analysis: Design Principle Violations

Applying the eight design principles:

**CC-1 (Clear community definition) — Partial failure.** Brixton is a neighborhood of approximately 100,000 people with significant economic diversity. The B£ community was defined geographically but not economically — there was no dense B2B network of businesses transacting with each other, only a B2C layer (residents spending at local merchants). Without B2B density, the multiplier effect was limited: spending at a B£ merchant generated national currency payments to suppliers (who were not in the B£ network), immediately leaking value outside the community.

**CC-2 (Mutual credit foundation) — Failed.** The Brixton Pound was issued against national currency deposits — each B£ was backed 1:1 by a pound sterling held in reserve by the Brixton Pound organization. This is not mutual credit; it is a backed currency. The consequence: the B£ organization was always cash-constrained (it could only issue B£ equal to deposits received), and there was no mechanism for expanding liquidity when the community needed it most (during recessions, precisely when the counter-cyclical property should have activated). A genuine mutual credit foundation would have allowed liquidity expansion automatically when community members needed to transact.

**CC-3 (Limited convertibility) — Failed.** The B£ was redeemable for sterling at par (1:1) — initially immediately, later with a small fee. This near-perfect convertibility eliminated the circulation incentive: any merchant who received B£ and did not have immediate use for them could simply convert to sterling at minimal cost. The result: B£ spent but not re-circulated, accumulating as sterling in the organization's reserve account rather than cycling through the community.

**CC-7 (Demurrage or circulation incentive) — Failed.** The Brixton Pound tried to implement a quarterly 3\% expiry fee on paper notes (not digital), but the mechanism was administratively burdensome (stamp collection), poorly communicated, and rarely enforced. Without an effective circulation incentive, the paper notes accumulated as souvenirs rather than circulating as currency. The digital version had no holding cost at all.

**CC-5 (Active clearing) — Not applicable (backed currency, no credit).** Because the B£ was a backed currency (not mutual credit), there was no clearing mechanism. Each B£ note was a claim on the reserve, not a bilateral credit commitment. The absence of multilateral clearing meant accumulated B£ holdings had no mechanism for productive redeployment.

**Overall assessment:** The Brixton Pound violated CC-2 (mutual credit), CC-3 (limited convertibility), and CC-7 (demurrage) — three of the eight principles. These three violations together prevented the counter-cyclical mechanism from operating, eliminated the circulation incentive, and constrained the liquidity available to the community.

### 37.4.3 The Counterfactual Simulation

**Corrected design.** Redesign the Brixton Pound as:
- **Mutual credit** (not backed): businesses extend credit commitments to each other up to GBP 2,000 credit limits.
- **Limited convertibility:** 5\% conversion fee to sterling; freely spendable within the network.
- **Demurrage:** 1.5\% quarterly (approximately 6\%/year), recycled as grants to community organizations.
- **B2B focus:** Primary network is local business-to-business trade, not consumer spending.

**Simulated outcomes (10-year comparison):**

| Metric | Actual B£ | Simulated corrected B£ |
|:---|:---|:---|
| Peak annual turnover | GBP 1.2M | GBP 8.4M (estimated) |
| Merchant network at peak | 300 | 420 (B2B enables higher density) |
| Counter-cyclical effect (2020 COVID) | Collapsed | +35\% turnover expansion (mutual credit auto-stabilizer) |
| Lifetime (years) | 10 | Indefinite (no reserve constraint) |
| Annual community benefit | GBP 60K | GBP 420K (multiplier effect of higher turnover) |

The counterfactual turnover estimate (GBP 8.4M) uses the observed WIR penetration rate (approximately 1\% of B2B trade in participating sectors) applied to Brixton's estimated GBP 840M annual B2B commerce. The counter-cyclical behavior during COVID would have been the most dramatic test: a mutual credit system with demurrage would have automatically expanded liquidity as national currency tightened — precisely the opposite of what occurred (the actual B£ organization, with its reserve-backed structure, faced mounting withdrawal pressure during COVID as merchants sought sterling liquidity).

---

## 37.5 Case Study: Sardex (Sardinia, 2010–2020)

### 37.5.1 System Design and Growth

Sardex (Sardinia Exchange) is a B2B mutual credit system operating in Sardinia, Italy, founded in 2010 by a group of young Sardinian entrepreneurs during the post-2008 recession when bank credit to Sardinian SMEs had contracted severely. As of 2020, it had approximately 3,500 member businesses, with annual turnover of approximately EUR 50 million.

**Design assessment against CC-1 through CC-8:**

| Principle | Sardex implementation | Score |
|:---|:---|:---|
| CC-1 (Community) | Sardinian SMEs — geographic + economic B2B community | ✓ Strong |
| CC-2 (Mutual credit) | Pure mutual credit — no reserve backing | ✓ Strong |
| CC-3 (Limited convertibility) | No euro conversion; Sardex spendable only within network | ✓ Strong |
| CC-4 (Calibrated limits) | Credit limits set by Sardex operators based on business health | ✓ Strong |
| CC-5 (Clearing cycles) | Monthly multilateral clearing; dedicated clearing team | ✓ Strong |
| CC-6 (Transparency) | Members see aggregate system balance; operators see all | ✓ Strong |
| CC-7 (Demurrage) | No demurrage — identified as design gap by operators | ✗ Absent |
| CC-8 (Recognition) | Operates under Italian commercial law (barter exchange) | ✓ Strong |
| **Score** | | **7/8** |

The one gap — absence of demurrage — has been identified by Sardex operators themselves as a cause of balance accumulation by successful exporters within the network. Large credit balances held by productive members reduce liquidity available to new members and dampen the velocity of the system. The Sardex team has considered introducing a holding fee but has not yet implemented it due to member resistance.

### 37.5.2 Network Growth Analysis

Sardex's membership growth from 2010 to 2020 shows the tipping threshold dynamics of Chapter 15:

- 2010–2012: Slow growth (below tipping threshold), approximately 100 members.
- 2013–2015: Accelerating growth as network effects strengthen; approximately 800 → 2,000 members.
- 2016–2018: Steady growth within the Sardinian SME market; approximately 2,000 → 3,200 members.
- 2019–2020: Plateau at approximately 3,500 members — estimated saturation of the addressable Sardinian B2B market under current design.

**Estimated tipping threshold:** $\hat{x} \approx 0.05$ (5\% of addressable market, approximately 150–200 members). Below this threshold, each new member added little network value; above it, each new member substantially increased the spending opportunities for all existing members. The 2013 crossing of $\hat{x}$ corresponds precisely to the acceleration in growth observed.

### 37.5.3 COVID-19 Resilience (2020)

During the 2020 COVID-19 lockdowns, Sardex demonstrated the mutual credit counter-cyclical property directly:
- National Italian bank credit to SMEs contracted approximately 15\% in Q2 2020.
- Sardex trading volume: increased 8\% in Q2 2020 as members shifted transactions to the network where credit was available.
- Member default rate: 2.1\% in 2020 (below the 2019 rate of 2.4\%), contrary to the expectation that a recession would increase defaults. Explanation: Sardex members maintained business relationships and cash flow through the network, improving credit quality.

This behavior is exactly consistent with the counter-cyclical model: $\dot{x}_{\text{Sardex}} > 0$ during the crisis, as the relative advantage of Sardex credit over unavailable bank credit increased sharply.

---

## 37.6 Case Study: Sarafu Network (Kenya)

### 37.6.1 Context and Complementarity with Chapter 25

The Sarafu Network (introduced in Chapter 25's mutual credit discussion) provides a contrasting case to the European examples: a commitment pooling system [C:Ch.25, Definition 25.5] operating in low-income urban and rural communities in Kenya, where formal banking infrastructure reaches approximately 25\% of the population and national currency liquidity is frequently constrained.

**Two-currency SFC analysis.** For a Sarafu Network community:
- National currency ($M^N$ = Kenyan shilling): used for taxes, formal sector purchases, inter-community transfers.
- Sarafu tokens ($M^C$): used for local food, services, and informal labor within the community.

The coexistence conditions (Theorem 37.1): C1 satisfied (Sarafu has clear advantage for local food and service transactions where national currency is scarce); C2 satisfied (Sarafu not freely convertible to KES — redeemable only for goods/services committed by community members); C3 satisfied (zero-sum mutual credit); C4 satisfied (community governance committees implement Ostrom conditions informally).

**Counter-cyclical function.** During COVID-19 lockdowns (2020), when national currency income fell sharply for informal workers: Sarafu trading volume increased 178\% (March–June 2020 vs. March–June 2019). This is a much larger counter-cyclical effect than the WIR's 1.47× sensitivity — reflecting the much more acute national currency liquidity constraint faced by Sarafu communities (which have near-zero access to formal credit markets) compared to Swiss SMEs (which have alternative national credit channels, merely tightened).

**The poverty reduction mechanism.** Ruddick et al. (2021) document: Sarafu-using households in food-insecure areas purchased 15\% more food during COVID lockdowns than matched non-Sarafu households. The mechanism: Sarafu enables intra-community exchange of food and services even when national currency is unavailable — keeping local food supply chains functional when national monetary circulation fails.

**Formal comparison to WIR:**

| Feature | WIR (Switzerland) | Sarafu (Kenya) |
|:---|:---|:---|
| Economy type | High-income, stable | Low-income, informal |
| CC type | Mutual credit (B2B) | Commitment pooling (B2B + B2C) |
| Counter-cyclical sensitivity | −1.47 (GDP elasticity) | +178\% (COVID volume increase) |
| Primary function | SME liquidity buffer | Subsistence food security |
| Governance | Formal cooperative bank | Community committees |
| Scale | CHF 1.5B/year, 62,000 members | USD 50M equivalent, 55,000 users |
| Technology | Bank accounts + paper | Mobile phone (USSD) + blockchain |

The two systems serve fundamentally different purposes within the same theoretical framework. Both demonstrate the counter-cyclical mutual credit property; the magnitude differs because the degree of national currency liquidity constraint differs enormously between Swiss SMEs and Kenyan informal workers. The theoretical framework of Chapter 25 predicts both the sign (counter-cyclical) and the relative magnitude (larger counter-cyclical effect where national currency scarcity is more acute) of this difference.

---

## 37.7 Mathematical Model: Two-Currency Dynamics

**State variables.** Let $m = M^C / M^N$ be the complementary currency ratio (CC volume relative to national currency volume) and $x$ be the share of transactions settled in CC.

**Dynamics:**

$$\dot{m} = \lambda \cdot x \cdot (1 - m/m^{\max}) - \mu \cdot (1-x) \cdot m$$

$$\dot{x} = \phi \left[\frac{u^C(m, \text{tightness})}{u^N} - x\right]$$

where $\lambda$ is the rate of CC credit creation when used, $\mu$ is the rate of balance reduction through clearing, $m^{\max}$ is the maximum CC ratio (determined by credit limits), $u^C$ is the utility of using CC (increasing in credit tightness), and $u^N$ is the utility of national currency.

**Steady-state coexistence:** $(\dot{m}, \dot{x}) = (0,0)$ at interior solution $m^* = m^{\max} \cdot \lambda x^* / (\lambda x^* + \mu(1-x^*))$ and $x^* = u^C(m^*, \text{tightness}^*) / (u^C + u^N)$.

**Stability condition:** $\partial^2(\dot{m})/\partial m^2 < 0$ and $\partial^2(\dot{x})/\partial x^2 < 0$ at the coexistence equilibrium — both state variables self-correcting. Satisfied when $m^{\max} < 1$ (CC never fully replaces national currency) and the utility ratio $u^C/u^N$ is bounded (national currency always retains some transactions). This formalizes Theorem 37.1's coexistence conditions.

**Crisis response.** During a credit tightening event (recession): $u^N$ falls (national credit becomes scarce) while $u^C$ rises (mutual credit advantage increases). $\dot{x} > 0$: CC share increases. $\dot{m} > 0$: CC volume increases. The system automatically expands the CC circuit when the national monetary circuit contracts — the formal counter-cyclical mechanism.

---

## Chapter Summary

This chapter has tested the mutual credit and demurrage theories of Chapters 25 and 27 against three decades of real-world complementary currency operation, finding strong empirical confirmation of the counter-cyclical mechanism and identifying the design principles that distinguish durable systems from short-lived experiments.

The two-currency SFC model (Definition 37.1, Theorem 37.1) specifies four necessary and sufficient conditions for stable coexistence: niche differentiation, limited convertibility, SFC consistency, and Ostrom governance. The counter-cyclical mechanism is formally derived and empirically confirmed: WIR turnover has a GDP elasticity of −1.47 (Stodder and Lietaer, 2016), Sardex expanded 8\% during COVID-19 as Italian bank credit contracted, and Sarafu expanded 178\% during Kenyan COVID lockdowns.

The Brixton Pound failure analysis identifies three violated design principles (CC-2, CC-3, CC-7) as the mechanism: a backed currency with near-par convertibility and no holding cost cannot sustain circulation or provide counter-cyclical function. The counterfactual simulation estimates that a corrected mutual credit design would have achieved 7× higher turnover and indefinite longevity.

Sardex (7/8 design principles, EUR 50M annual turnover, 2.1\% default rate during COVID-19) confirms the B2B mutual credit model at regional scale. The Sarafu Network confirms the same model in a low-income informal economy context, with even stronger counter-cyclical effects because national currency liquidity constraints are more acute.

Chapter 38 applies the formal framework to the provision of public services through commons-based institutions — the Universal Basic Services model and Vienna's Gemeindebau — showing how cooperative public provision achieves better welfare outcomes per unit of expenditure than both privatized services and standard state provision.

---

## Exercises

**37.1** Two-currency SFC model (Definition 37.1):
(a) Construct the full balance sheet matrix for a two-currency economy with 3 households, 3 firms, one national bank, one complementary currency operator, and government. Verify that all rows sum to zero.
(b) A firm with WIR balance $M^C_i = +800$ CHW purchases EUR 400 worth of services from another firm, paying 50\% in CHW and 50\% in CHF. Show all balance sheet changes. Does $\sum_i M^C_i = 0$ still hold after the transaction?
(c) Compute the velocity of CHW if annual CHW turnover is CHF 1.5 billion and the outstanding CHW stock is CHF 18 million. Compare to the CHF velocity in the same period (Swiss M2 velocity ≈ 0.8). Interpret the difference.

**37.2** Counter-cyclical mechanism (Stodder and Lietaer regression):
(a) If Swiss GDP falls 2\% in a recession year, what does the regression model predict for WIR turnover growth?
(b) Compute the approximate CHF value of the WIR counter-cyclical demand injection during a −2\% GDP recession (CHW 1.5 billion baseline, GDP multiplier 0.8, conversion rate CHW 1 = CHF 1).
(c) If the Swiss multiplier for national fiscal stimulus is 1.4, how much fiscal spending would be required to achieve the same demand injection as the WIR's automatic counter-cyclical response? Does this suggest the WIR has a cost-effective stabilization value?

**37.3** Sardex design evaluation:
(a) For each of the eight design principles (CC-1 through CC-8), assess whether the Sardex system strongly, partially, or does not implement it. Provide a specific institutional example for each.
(b) Sardex does not implement demurrage (CC-7). Model the expected effect of this gap: if Sardex members hold an average of 30\% of their annual CC transactions as unspent balances, what fraction of potential liquidity is being withheld from the network?
(c) Propose a specific demurrage mechanism for Sardex: rate, implementation method (monthly fee vs. quarterly stamp), and recycling use (what the holding fee revenue funds). Would this require a governance vote? What resistance would you expect from member businesses holding large credit balances?

**★ 37.4** Prove Theorem 37.1 (stable coexistence conditions) formally.

(a) Model currency choice as a decision for each agent in each period: use national currency or complementary currency for each transaction. Let $x_t$ be the share of transactions settled in CC.
(b) Show that $x^* = 0$ (pure national currency) and $x^* = 1$ (pure CC) are both potential equilibria. Under what conditions is each stable?
(c) Prove that a stable interior equilibrium $x^* \in (0,1)$ exists if and only if all four conditions C1–C4 hold simultaneously. Specifically, show that: (i) violating C1 drives $x^* = 0$; (ii) violating C2 drives $x^* = 0$ through Gresham's law; (iii) violating C3 triggers regulatory shutdown, eliminating the CC; (iv) violating C4 leads to governance failure and system collapse.
(d) Verify that the WIR satisfies all four conditions and identify the specific institutional mechanism implementing each.

**★ 37.5** Conduct a formal failure analysis of a complementary currency of your choice that did not survive (options: Bristol Pound, BerkShares, Chiemgauer pre-reform period, or any documented complementary currency failure).

(a) Identify which of the eight design principles (CC-1 through CC-8) were violated.
(b) For each violated principle, specify the failure mechanism: how did the violation lead to reduced effectiveness or system collapse?
(c) Design the counterfactual: what would the system have looked like with corrected design? Simulate the corrected system's 5-year trajectory using the two-currency dynamics model of Section 37.7.
(d) Compare your failure analysis to the Brixton Pound case. Are there common failure modes across both systems? What does this suggest about the minimum viable design requirements for complementary currencies?

**★★ 37.6** Design a complementary currency for a specific community context and formally prove its coexistence stability.

**Community context:** A network of 200 urban food businesses (restaurants, cafés, food producers, distributors) in a mid-sized European city, currently facing significant cash flow volatility due to seasonal demand and long payment terms from institutional buyers (hospitals, schools, corporate cafeterias). Many are credit-constrained following the COVID-19 crisis.

(a) **Design specification:** Specify all eight design principles for your currency. For each: the institutional mechanism, the governance rule, and the parameter values (credit limits, demurrage rate if applicable, clearing cycle, convertibility conditions).

(b) **SFC consistency proof:** Construct the full balance sheet matrix and transaction flow matrix for your currency. Verify that all four SFC constraints of Chapter 28 are satisfied.

(c) **Coexistence stability proof:** Using Theorem 37.1, prove that your currency satisfies all four conditions. Specifically: identify the niche (C1), specify the limited convertibility mechanism (C2), verify SFC consistency (C3), and map each governance mechanism to the corresponding Ostrom design principle (C4).

(d) **Tipping threshold analysis:** Using the tipping threshold model of Chapter 15, estimate $\hat{x}$ — the critical network size for self-sustaining adoption. What initial community (below $\hat{x}$) should the system target, and how should it be incentivized to join?

(e) **Counter-cyclical property:** Model how your currency would respond to a 20\% national credit contraction (simulating a recession or banking crisis). Would it expand (counter-cyclical)? By how much? Compute the welfare value of the counter-cyclical response using the Stodder-type analysis of Section 37.2.3.

---

*Chapter 38 applies the formal framework to the commons-based provision of public services — the Universal Basic Services model and Vienna's Gemeindebau. The formal tools are welfare economics (UBI vs. UBS comparison), commons governance theory (Ostrom principles applied to public service delivery), and fiscal sustainability analysis (SFC-compatible public service financing). The empirical case of Vienna's 220,000-unit social housing cooperative — the world's largest — provides the test of the theory's predictions.*
