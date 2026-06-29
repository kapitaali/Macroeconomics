# Chapter 44: Capstone Project — Designing a Local Cooperative-Regenerative Economy

> *"An economist is someone who, when they see something working in practice, wonders if it can work in theory. A cooperative-regenerative economist reverses this: they see something that works in theory and ask how to make it work in practice."*
> — attributed, original

> *"The purpose of a design is not to be correct. It is to be good enough to act on, and honest enough to learn from."*
> — Herbert Simon, *The Sciences of the Artificial* (1969, paraphrased)

## Preamble: The Nature of This Chapter

This chapter is a structured design project, not a lecture. It does not present new theory — Parts II–VI have done that. It does not present new empirical evidence — Part VII has done that. It does not develop transition strategies — Part VIII has done that. Instead, it asks you to do what the entire trilogy has been building toward: take the complete analytical toolkit and apply it to the design of a real economy.

The project specification requires five components, each drawing on multiple parts of the book. Component 1 (institutional architecture) draws on Parts II, III, and V. Component 2 (ecological embedding) draws on Part IV. Component 3 (cooperative game analysis) draws on Part II. Component 4 (simulation) draws on Parts II–V. Component 5 (transition pathway) draws on Part VIII. The integration of all five components is the capstone: a coherent, internally consistent design for a local cooperative-regenerative economy that can be defended analytically and is actionable in practice.

The worked partial solution (Section 44.2) demonstrates the methodology for a specific locality — a medium-sized UK market town — completing Components 1 and 2 in full and sketching the other three. The worked solution is not meant to be copied; it is meant to show the quality of reasoning and the level of formal specificity that a complete design requires, so that when you design your own locality, you know what the standard looks like.

---

## 44.1 Project Specification

### The Locality

Choose a real or hypothetical locality with population between 50,000 and 500,000. The locality should be:

- **Specific:** A named place (or a clearly specified hypothetical place with realistic parameters) rather than a generic "average city."
- **Data-accessible:** Sufficient public data exists to calibrate the models (national census, local authority accounts, regional input-output tables, ecological surveys).
- **Contextually grounded:** The design should reflect the locality's specific economic structure, natural capital endowment, governance capacity, and social history.

If you choose a real locality, acknowledge what you do not know and what assumptions you are making. If you choose a hypothetical locality, specify its parameters with sufficient precision that the design is analytically tractable.

---

### Component 1: Institutional Architecture

**1a. Governance structure.** Using the Three-Layer Coordination Stack (Chapter 29, Definition 29.2) and the polycentric governance framework (Chapter 14), specify the governance structure of your locality's cooperative-regenerative economy:

- **Layer 1 (Direct mutual coordination):** Which activities are coordinated through commons governance, cooperative enterprise, and peer-to-peer networks? Specify the governance structure for at least three Layer 1 institutions (e.g., a housing cooperative, a food commons, a mutual credit network), including their Ostrom design principle compliance (score each on DP1–DP8).
- **Layer 2 (Generative markets):** Which markets operate at Layer 2, and how are ecological prices (natural capital levies, carbon pricing) incorporated? Specify the price adjustments relative to current market prices for at least three key goods/services.
- **Layer 3 (Biophysical planning):** Which constraints from the Planetary Boundaries and the local GTA allocation [C:Ch.17] are binding for your locality? How are they governed?

**1b. Enterprise mix.** Design the cooperative enterprise sector for your locality:

- Current enterprise mix: proportion of employment and output in conventional firms, cooperatives, public sector, and informal economy.
- Target enterprise mix at year 20 of the transition.
- Specific cooperative enterprise designs for three sectors (using the Chapter 34 sector design templates: manufacturing, professional services, platform technology, or others appropriate to your locality).
- The Furubotn-Pejovich problem: how is under-investment in cooperative capital addressed?

**1c. Monetary system.** Design the hybrid monetary system for your locality (Chapter 28, Proposition 28.1):

- What fraction of the money supply is sovereign money, mutual credit, and demurrage?
- Specify the mutual credit network: membership criteria, credit limits, clearing cycles, governance structure (Ostrom score).
- Specify the demurrage token: rate, coverage (which transactions), revenue use.
- Verify SFC consistency: show that the transaction flow matrix balances.

---

### Component 2: Ecological Embedding

**2a. Material and energy flow analysis.** Using the MFA framework of Chapter 17:

- Identify the five largest material flows through your locality (inputs, outputs, internal stocks).
- Compute the local DMC (Domestic Material Consumption) per capita and compare to the GTA-based Planetary Boundary allocation.
- Identify which Planetary Boundaries are currently violated by your locality's material and energy footprint.

**2b. Natural capital accounting.** Using the SFC-N framework of Chapter 18:

- Identify the five most important natural capital stocks in your locality (e.g., agricultural soil, urban green space, local watercourses, atmospheric carbon, biodiversity).
- For each stock: current level, regeneration rate $r_j$, depletion rate $\mathcal{D}_j$, critical threshold $N_j^{\text{critical}}$, and shadow price $p^N_j$.
- Construct the Provisioning Balance Sheet (PBS) for your locality. Is $\dot{\text{PBS}} \geq 0$ (is the locality maintaining its natural capital stocks)?

**2c. Regeneration condition verification.** For each natural capital stock:

- Is the Stewardship Condition $\dot{N}_j \geq 0$ currently met?
- If not: what is the investment needed to restore $\dot{N}_j \geq 0$? At what shadow price?
- How does the cooperative institutional design of Component 1 contribute to meeting the Stewardship Condition (through cooperative landscape governance, regenerative agriculture PES, etc.)?

---

### Component 3: Cooperative Game Analysis

**3a. Characteristic function.** Construct the cooperative game for your locality's economy:

- Sectors as players: identify $n$ major economic sectors (agriculture, manufacturing, retail, services, public sector, cooperative enterprise sector, etc.).
- Estimate $v(S)$ for all singleton coalitions (competitive equilibrium allocation) and the grand coalition $v(\mathcal{M})$ (cooperative-regenerative equilibrium allocation).
- Compute the cooperative surplus fraction: $\sigma_v = (v(\mathcal{M}) - \sum_i v(\{i\})) / v(\mathcal{M})$.

**3b. Core stability analysis.** Verify that the cooperative allocation is core-stable (Chapter 6, Theorem 6.1, Bondareva-Shapley):

- Is the game balanced? Check the Bondareva-Shapley condition for at least three coalition configurations.
- If the core is empty: identify which design modifications (institutional reforms, governance changes) would restore core stability.

**3c. Shapley value allocation.** Compute the Shapley value for each sector:

- Using the formula from Chapter 3 (exact for small $n \leq 6$; Monte Carlo approximation for larger $n$).
- Compare the Shapley allocation to the competitive equilibrium (factor income shares). Which sectors gain and which lose?
- Compute the Gini coefficient of sector-level income under the competitive and Shapley allocations. What is the inequality reduction from cooperative institutions?

---

### Component 4: Simulation

**4a. Agent-based model sketch.** Using the BDI framework of Chapter 10:

- Specify the agent population: at least three agent types (households, cooperative enterprises, conventional firms), with heterogeneous beliefs, desires, and intentions.
- Specify the interaction rules: which agents interact in which layers of the coordination stack, and what determines outcomes of each interaction.
- Specify the ecological feedback: how do agents' production and consumption decisions affect the natural capital stocks of Component 2?
- Run the ABM for 20 periods (or specify the pseudocode and explain what a 20-period simulation would show). Report the key outcome variables: employment, income distribution, natural capital stocks, cooperative sector share.

**4b. SFC macro model.** Using the SFC-N framework of Chapters 18 and 28:

- Specify the sector-level TFM (Transaction Flow Matrix) for your locality.
- Verify that all rows sum to zero (SFC consistency).
- Run the SFC model for 20 years under: (i) baseline (current institutions) and (ii) cooperative-regenerative design (Component 1). Report the key macroeconomic variables.

**4c. Sensitivity analysis.** Identify the three model parameters with greatest uncertainty and conduct sensitivity analysis:

- For each parameter, compute the change in the key outcome variable (IPI at year 20) per unit change in the parameter.
- Identify the "robust" findings — those that hold across the full plausible range of parameter values — and the "fragile" findings — those that reverse under alternative parameter assumptions.

---

### Component 5: Transition Pathway

**5a. Policy package.** Using the eleven instruments of Chapter 41 and the predistribution instruments of Chapter 41b:

- Select the six most relevant instruments for your locality, including at least one from each of: Chapter 41's price/subsidy instruments; Chapter 41's institutional reforms; and Chapter 41b's predistributive instruments (employment guarantee, working time reduction, or care OVA).
- For each: specify the rate/level, the implementation timeline, the revenue generated or cost incurred, and the expected effect on $\theta$, $v_{\text{CRE}}$, and $\phi_{\text{landscape}}$.
- Verify fiscal balance: total revenue ≥ total cost across all instruments.
- Compute the updated Transition Tipping Point $\hat{x}^{\text{CRE}}$ under the policy package. Is it below the current cooperative sector share?

**5b. Democratic governance framework.** Using Chapter 41c's institutional tools:

- Design a local citizens' assembly process for determining the welfare weight vector $\mathbf{w}$ for your locality's MPD planning. Specify: stratification criteria, evidence base, justificatory constraint, supermajority threshold.
- Specify the SIW accounts for your locality: which dimensions of $\mathcal{W}$ (current wellbeing), $\mathcal{S}$ (sustainability), and $\mathcal{I}$ (inclusion) are measurable with available data?
- Design a Participatory Social Audit for your locality's cooperative enterprise zone, specifying all five PSA components (monitors, verified accounts, audit cycle, sanctioning authority, escalation mechanism).

**5c. Experimentation strategy.** Using the MAB framework of Chapter 42:

- Identify three institutional design variants (arms) that should be tested in parallel before commitment to the full design.
- Specify the outcome metrics for each arm and the evaluation timeline.
- Using Thompson Sampling logic, specify when you would commit to the best arm and begin full replication.

**5d. Stakeholder engagement.** Using the political economy model of Chapter 41 and the coalition analysis:

- Identify the five most important stakeholder groups in your locality and their expected position on the cooperative-regenerative transition (support, neutral, oppose).
- Design the engagement process: how are opponents identified and addressed, how are supporters organized, and what early wins build the coalition?
- Apply the Proposition 41.1 (Winning Coalition) analysis: does your coalition constitute a democratic majority?

---

## 44.2 Worked Partial Solution: Darlington, County Durham (UK)

We demonstrate the methodology for Darlington — a market town of approximately 107,000 in County Durham, northeast England, with a historical industrial economy (railway engineering, manufacturing) that has diversified into service sectors following deindustrialization. It has a small existing cooperative sector (principally the Darlington Building Society, a mutual mortgage provider) and sits within one of England's highest-deprivation regions by the IMD (Index of Multiple Deprivation).

Darlington makes an interesting design case for three reasons: its cooperative history (it was near the Rochdale pioneers' era of cooperative organizing in northern England); its post-industrial transition (typical of many northern English towns facing deindustrialization); and its ecological situation (the River Tees catchment provides significant natural capital but faces agricultural runoff and historical industrial contamination challenges).

*Note: The worked solution completes Components 1–2 in full and provides the framework for Components 3–5. Students completing the project for a different locality should use this as a methodological template, not a design to be copied.*

---

### Component 1: Institutional Architecture (Darlington — Full)

#### 1a. Governance Structure

**Layer 1 — Direct Mutual Coordination:**

*Institution 1a-i: Darlington Housing Cooperative.*

Design: 800-unit community land trust (CLT) housing cooperative, developed on brownfield sites along the River Tees corridor. Governance: Cosmo-Local structure — building committees (10–15 members each), area assemblies (50–80 members), general assembly (all 800 member-households).

Ostrom assessment:
| DP | Implementation | Score |
|:---|:---|:---|
| DP1 | Membership defined by income limit (£40K individual/£60K household) and Darlington residency | 2/2 |
| DP2 | Rules differentiated by property type (family, elder, young professional), location (riverside, town centre) | 2/2 |
| DP3 | Building committees govern maintenance; area assemblies govern shared services; general assembly governs major decisions (2/3 majority) | 2/2 |
| DP4 | Annual building surveys; energy monitoring through smart meters; quarterly financial reporting to members | 2/2 |
| DP5 | Warning → mediation → lease review for anti-community behaviour; transparent and consistently applied | 1.5/2 |
| DP6 | In-house mediation officer; Darlington Council ombudsman as backstop | 1.5/2 |
| DP7 | Community Land Trust legal structure (UK CLT Act 2008); registered social landlord status | 2/2 |
| DP8 | Building committees → Area assemblies → CLT board → Darlington Community Wealth Alliance (see below) | 2/2 |
| **Total** | | **15/16** |

*Institution 1a-ii: Darlington Mutual Credit Network (DMCN).*

Design: B2B mutual credit for Darlington's 1,800 SMEs. Credit limits: 5% of annual turnover (minimum £2,000, maximum £40,000). Monthly multilateral clearing using Algorithm 25.1. Governance: cooperative association of member businesses; five-member elected board. Ostrom score: designed to achieve 14/16 (DP6 conflict resolution is the typical weakness — DMCN will implement an online arbitration panel from year 1 to address this).

*Institution 1a-iii: Tees Valley Food Commons.*

Design: A regional food commons governing the agricultural land within 20km of Darlington for cooperative food production, community market gardens, and regenerative landscape management. Membership: 85 farms, 12 market gardens, 6 food processing cooperatives, 3 community food hubs. Governance: landscape cooperative structure from Chapter 36 (Proposition 36.3). Annual GTA allocation for the Tees Valley agricultural landscape: calibrated from Chapter 17's equal per-capita allocation.

**Layer 2 — Generative Markets:**

Three key price adjustments for Darlington's economy:

- **Agricultural land:** Current market price ≈ £10,000/ha (productive arable). Shadow price including soil carbon, water regulation, biodiversity: £10,000 + £3,600 (ecosystem services) = £13,600/ha effective value. LVT at 1.5% of unimproved land value captures £150/ha/year, redirected to Tees Valley Food Commons PES payments. This reduces the incentive gap between conventional and regenerative farming.

- **Commercial property:** Darlington town centre has a 28% retail vacancy rate (2023). LVT on commercial land (assessed at 1.5% of unimproved value) reduces speculative holding and incentivizes productive use. Estimated yield: £2.1M/year from commercial land LVT.

- **Energy:** Carbon levy at £80/tonne CO₂e applied to all commercial energy use within Darlington (approximately 180,000 tonnes CO₂e/year from commercial sector): £14.4M/year. Recycled 60% as Cap-and-Share dividend (£85/resident/year), 40% to ecological restoration fund.

**Layer 3 — Biophysical Planning:**

Binding constraints for Darlington:
- Atmospheric carbon: Darlington's share of UK carbon budget (1.5°C pathway) = 107,000/67M × 0.55 Gt CO₂e = approximately 0.88 Mt CO₂e/year. Current emissions: approximately 0.72 Mt CO₂e/year. Currently within budget, but transport sector (0.31 Mt) is above proportional allocation and requires focused intervention.
- Land system change: The Tees catchment requires 15% woodland cover for flood attenuation (currently 8%). Deficit: 7% = approximately 2,100 ha of new woodland required. Governed through Layer 1's Tees Valley Food Commons and funded through Layer 2's natural capital levy.
- Freshwater: River Tees nitrogen loading from agriculture is above the safe operating space. Target: 25% reduction in agricultural nitrogen runoff within 10 years.

#### 1b. Enterprise Mix

**Current (2024):**
- Conventional private firms: 68% of employment
- Public sector: 24% of employment
- Cooperatives and social enterprises: 5% of employment
- Informal and self-employed: 3%

**Target (Year 20):**
- Conventional private firms: 42% of employment
- Public sector: 22% of employment
- Cooperatives and social enterprises: 30% of employment
- Informal and self-employed: 6%

**Three sector cooperative designs:**

*Manufacturing (Darlington Engineering Cooperative):* 150-member precision engineering cooperative, building on Darlington's historical engineering skills. Wage structure: 1:5 ratio (lowest to highest). Capital: member accounts (20% of income retained, 4% interest) plus Darlington Community Finance cooperative credit line. Governance: three-tier Cosmo-Local structure. Key challenge: Furubotn-Pejovich problem addressed through 20-year asset depreciation smoothing in member accounts — members who leave receive the full present value of their capital account contribution as a structured annuity, eliminating the horizon problem.

*Professional services (North East Law Cooperative):* 40-member legal and financial advisory cooperative serving cooperatives, social enterprises, and community organizations. OVA allocation with contribution dimensions: billable hours (weight 0.7), client origination (0.4), community governance (0.5), mentorship (0.3). Maximum OVA multiplier: 3.5× base (£45,000/year).

*Platform (Darlington Local Platform Cooperative):* Multi-stakeholder cooperative platform connecting local service providers (tradespeople, childminders, tutors, delivery workers) with local customers. Governance: 60% worker members / 25% customer members / 15% steward members (technology maintenance). Smart contract governance using the template of Chapter 35. Cold-start solution: anchor institution (Darlington Council) commits to using the platform for all local service procurement from year 1.

#### 1c. Monetary System

**Hybrid design (60/25/15):**

- **Sovereign money (60%):** Supplied by National Bank of the UK (in the sovereign money transition scenario) or, short-term, by the Darlington Community Wealth Alliance accessing sovereign money through a municipal bond backed by LVT and carbon levy revenue.

- **Mutual credit (25%):** The DMCN provides approximately £28M of effective liquidity to Darlington's SME sector (25% × £112M estimated annual B2B commerce). Credit limits: £2,000–£40,000 per member. Monthly clearing.

- **Demurrage tokens (15%):** "Darlo" tokens issued at 2% quarterly demurrage. Accepted by all cooperative enterprises, public sector facilities, and participating conventional businesses (targeting 40% of city businesses by year 3). Revenue from demurrage (2% × £6.7M circulation = £134,000/year) directed to: 50% Tees Valley ecological restoration, 30% community social fund, 20% Darlington Community Wealth Alliance operating costs.

**SFC consistency (partial verification):**
- Total monetary stock: £[sovereign] + £[DMCN credits] + £[Darlo tokens] = total ≡ productive capacity + investment claims.
- DMCN row: $\sum_i b_i^{\text{DMCN}} = 0$ by construction (zero-sum mutual credit).
- Darlo row: $\sum_i M_i^{\text{Darlo}} = M^{\text{Darlo}}_{\text{outstanding}}$ (IA liability); demurrage revenue flows from holders to IA; IA reissues as ecological restoration grants. TFM balance: confirmed.

---

### Component 2: Ecological Embedding (Darlington — Full)

#### 2a. Material and Energy Flow Analysis

**Five largest material flows (Darlington, estimated from BEIS, Defra, and ONS regional data):**

| Material flow | Annual quantity | GTA-based allocation | Current status |
|:---|:---|:---|:---|
| Carbon emissions (all sectors) | 720,000 t CO₂e | 880,000 t CO₂e | Within budget |
| Food imports | 42,000 t | Primarily local (target) | 78% imported — overshoot |
| Construction materials | 185,000 t | 120,000 t | +54% above allocation |
| Water extraction (Tees) | 18.5 Mm³ | 22.0 Mm³ | Within allocation |
| Nitrogen (agricultural) | 3,200 t to Tees | 2,000 t (safe threshold) | +60% above safe level |

**Two Planetary Boundaries currently violated:** Land system change (food import intensity, construction overshoot) and biogeochemical flows (nitrogen to Tees).

#### 2b. Natural Capital Accounting

**Provisioning Balance Sheet (Darlington, £M):**

| Natural capital stock | Current value (£M) | Annual regeneration | Annual depletion | Shadow price | $\dot{N}$ |
|:---|:---|:---|:---|:---|:---|
| Agricultural soil (Tees Valley) | 285 | +8.5 | −12.3 | £2,400/t C | **−£3.8M/yr** |
| Urban green space | 45 | +1.8 | −0.9 | £800/ha/yr | +£0.9M/yr |
| River Tees water quality | 62 | +3.1 | −5.8 | £150/ML | **−£2.7M/yr** |
| Local biodiversity | 38 | +1.2 | −2.9 | £1,200/ha/yr | **−£1.7M/yr** |
| Carbon sequestration capacity | 118 | +4.2 | −7.1 | £80/t CO₂e | **−£2.9M/yr** |
| **PBS total** | **548** | **+18.8** | **−29.0** | | **−£10.2M/yr** |

**Stewardship Condition violated:** $\dot{\text{PBS}} = -£10.2\text{M/yr} < 0$. Darlington is liquidating its natural capital stock at £10.2M/year — the formal measure of its ecological unsustainability.

#### 2c. Regeneration Condition Restoration

**Investment needed to achieve $\dot{N}_j \geq 0$ for each stock:**

- Agricultural soil: Regenerative agriculture transition for 4,200 ha of Tees Valley farmland; PES payment needed = £3.8M/yr (Chapter 36, externality gap formula).
- River Tees: Agricultural nitrogen reduction (25% of current load = 800 t/yr); investment in buffer strips, wetland restoration = £2.7M/yr.
- Biodiversity: Landscape connectivity improvement (new wildlife corridors, hedgerow restoration) = £1.7M/yr.
- Carbon: 2,100 ha of new woodland; cost = £4,200/ha establishment + £400/ha/yr maintenance = £2.9M/yr annualized.

**Total investment for Stewardship Condition: £11.1M/year.** Funding sources: natural capital levy (£2.1M from commercial land LVT), carbon levy ecological restoration portion (£5.8M), cooperative landscape governance PES (£3.2M from Tees Valley Food Commons collective application). **Total: £11.1M — exactly meets the requirement.**

The cooperative institutional design of Component 1 (Tees Valley Food Commons as the landscape cooperative governance body, cooperative landscape governance achieving first-best social optimum per Proposition 36.3) is not merely complementary to the ecological embedding — it is the mechanism by which the £3.2M PES is channeled to the farmers who can most efficiently deliver it.

---

### Components 3–5: Framework (To Be Completed)

The following provides the analytical framework for Components 3–5. A full solution requires the data collection, calculation, and judgment specific to the chosen locality.

#### Component 3 Framework: Cooperative Game Analysis

For Darlington, the relevant sectors are: Agriculture (A), Manufacturing (M), Services (S), Public sector (P), Cooperative enterprises (C), and Housing/land (H) — a six-player game. Compute $v(\{A\}) + v(\{M\}) + \ldots + v(\{H\})$ from Darlington's regional I-O table (ONS NUTS3 data for County Durham). Compute $v(\mathcal{M})$ using the cooperative surplus formula calibrated to Chapter 29's Danish result — starting estimate: $\sigma_v \approx 0.15$ (15% cooperative surplus). Run exact Shapley calculation for $n = 6$ (feasible: $6! = 720$ permutations). Assess Bondareva-Shapley condition for at least five coalition configurations.

**Expected finding:** The cooperative sector (C) receives a Shapley premium above its competitive factor income share, because its marginal contribution to the grand coalition (through supply chain integration with local firms, governance public goods, and mutual credit network provision) exceeds its direct productive output. The housing/land sector (H) receives below its competitive factor share under Shapley — because its competitive income reflects speculative land rent rather than productive contribution to the local economy.

#### Component 4 Framework: Simulation

**ABM specification sketch:** Three agent types — Household (H: 46,000 households; heterogeneous income, preference for cooperative employment, housing, and food commons); Enterprise (E: 2,800 firms; heterogeneous sector, size, governance type); Natural Capital (NC: five stocks as defined above). Interaction rules: firms hire from household pool; households purchase from firms; Layer 1 institutions govern commons access; ecological stocks evolve according to Section 2b dynamics.

**SFC model specification:** Seven sectors (Households, Cooperative Enterprises, Conventional Firms, DMCN Mutual Credit Operator, Darlo Token Issuer, Government, Natural Capital). TFM has 18 flow rows (wages, investment, taxes, mutual credit creation/clearing, demurrage collection/reinjection, natural capital levies, PES payments). Full TFM in the companion computational model.

#### Component 5 Framework: Transition Pathway

**Policy package (Darlington-specific six instruments):** LVT on commercial land (raising £2.1M/yr), carbon levy at £80/t (£14.4M/yr), cooperative procurement mandate (20% of Darlington Council's £95M annual procurement = £19M from cooperatives), CEDF (£1.5M initial capitalization from Council reserves), regenerative PES through Tees Valley Food Commons (£3.2M/yr from natural capital levy), and DMCN formation support (£200K one-time setup, self-sustaining by year 2).

**Transition tipping point:** With the six instruments deployed, estimated $\Delta\theta = -0.14$, $\Delta v_{\text{CRE}} = +0.22$, $\Delta\phi = +0.12$. Updated $\hat{x}^{\text{CRE}} = (0.65 - 0.14 - 0.12)/(0.35 + 0.22) = 0.39/0.57 = 0.68$. Current cooperative share: 5%. Year 5 projected: 14%. Year 10 projected: 22%. **The tipping threshold (68%) is not crossed within 10 years** — but year 20 projection: 30%, still below 68%. This indicates that Darlington, as a single locality, cannot achieve system-level transition in isolation; national policy instruments are required to reduce $\theta$ further. The Darlington experiment contributes institutional learning and political evidence; it requires national partners to achieve systemic transition.

**The honest conclusion.** This is the most important finding of the Darlington exercise: local experimentation at the scale of a 107,000-person market town is necessary but not sufficient for systemic transition. The tipping threshold remains above the achievable cooperative share within 20 years under local policy alone. This is not a reason to abandon the local experiment — it is the clearest argument for why local experimentation must be coordinated with national policy reform (Chapter 41) and municipalist network learning (Chapter 42). Darlington becomes one node in the network that collectively pushes the national $\hat{x}^{\text{CRE}}$ below the tipping threshold.

---

## 44.3 Assessment Criteria

A complete capstone project should satisfy the following criteria:

**Analytical completeness.** All five components are addressed with the level of formal detail demonstrated in the Darlington worked solution. No component is reduced to a vague description — each has specific numbers, specific governance structures, specific calculations.

**Internal consistency.** The five components are mutually consistent: the monetary system (1c) is funded by the natural capital accounting (2b) revenues; the cooperative game analysis (3) reflects the enterprise mix (1b); the SFC model (4b) incorporates the TFM entries implied by the monetary system (1c). Inconsistencies are acknowledged and explained.

**Honest uncertainty.** The sensitivity analysis (4c) identifies the most important uncertainties, and the design is robust to plausible parameter variations where possible. Claims are qualified appropriately — the design is "good enough to act on" (as Simon's epigraph requires), not presented as optimal or certain.

**Contextual grounding.** The design is specific to the chosen locality — it could not be copy-pasted to a different place without substantial modification. The locality's specific economic structure, ecological endowment, and social history shape the design choices.

**The honest finding.** Like the Darlington analysis, a good capstone project should be willing to report what the analysis reveals even when the finding is not what was hoped. If the tipping threshold is not reachable locally, say so and explain why national coordination is needed. If the Stewardship Condition cannot be met within the design's funding capacity, say so and identify the gap. Honest analysis that reveals limitations is more valuable than optimistic analysis that conceals them.

---

## Chapter Summary

This chapter has provided the complete specification and a partial worked solution for the capstone project — the synthesis exercise that applies the full analytical toolkit of the trilogy to the concrete design of a local cooperative-regenerative economy.

The five components — institutional architecture, ecological embedding, cooperative game analysis, simulation, and transition pathway — require drawing on every major part of the book. Component 1 applies the Three-Layer Coordination Stack (Chapter 29), cooperative enterprise design (Chapter 34), polycentric governance (Chapter 14), and the hybrid monetary system (Chapter 28) — and now also the Employment Guarantee and care OVA from Chapter 41b, which are required elements of any cooperative enterprise zone's labour market design. Component 2 applies the MFA framework (Chapter 17), SFC-N natural capital accounting (Chapter 18), and the regeneration condition (Chapter 20). Component 3 applies cooperative game theory (Chapters 3, 6) and Shapley value analysis — extended to care labour following Proposition 41b.1. Component 4 applies agent-based modeling (Chapter 10) and SFC macro modeling (Chapter 28). Component 5 applies the policy instruments (Chapter 41), the predistribution instruments (Chapter 41b — employment guarantee, working time reduction, care OVA), the democratic planning and SIW governance framework (Chapter 41c — citizens' assembly process, social audit design, Democratic Accountability Loop), and the MAB experimentation methodology (Chapter 42).

The Darlington worked solution demonstrates the methodology at full detail for Components 1 and 2, and at framework level for Components 3–5. Its most important finding — that local experimentation is necessary but not sufficient for systemic transition, requiring national policy coordination — is the honest analytical conclusion that the project's formal framework produces. It is not a counsel of despair but a specification of what is needed beyond the local.

Chapter 45 provides the synthesis conclusion: assembling the book's argument, situating it in the history of economic thought, and extending the invitation to the reader to become a contributor to the research program that this book has opened rather than closed.

---

## Project Submission Guidance

**Length:** 8,000–12,000 words (approximately 25–35 pages) for an individual project; 12,000–20,000 words for a team project.

**Format:** The project should be organized by the five components, with clear section headers. Equations should be numbered and referenced. All data sources should be cited. Models should be reproducible — pseudocode or working code should be provided for all simulations.

**Collaboration:** The project is designed to be completed by teams of 2–4 students, reflecting the interdisciplinary nature of cooperative-regenerative economics. Individual components may be divided by student, but the integration (internal consistency across components) is a joint responsibility.

**Presentation:** Projects should be presented in a 20-minute seminar format, with 10 minutes for questions. The seminar should focus on: the most important design choices and their justification, the key uncertainties and how the design handles them, and the honest finding — what the analysis reveals about the limits of local-scale transformation.

---

*Chapter 45 provides the synthesis conclusion: the argument assembled, the place in history, the role of rigor, and the invitation to the reader to join the ongoing project of building an economics adequate to the challenges of the 21st century.*
