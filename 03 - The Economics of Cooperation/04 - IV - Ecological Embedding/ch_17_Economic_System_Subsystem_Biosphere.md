# Part IV: Core Theory III — Regeneration and Ecological Embedding

The economic theories developed in Parts II and III treat the economy as a social system: agents bargain, cooperate, form institutions, manage information, and govern their shared resources. This framing is necessary but not sufficient. The economy is also a physical system — a subsystem of the biosphere that extracts materials and energy from natural stocks, transforms them through production processes, and returns waste to natural sinks. No cooperative institution, however well-designed, can sustain itself if the biophysical substrate on which it depends is being depleted faster than it regenerates.

Part IV provides the formal tools for analyzing this physical dimension. Six chapters develop, in sequence, the material and energy accounting framework (Chapter 17), the stock-flow consistent model extended to natural capital (Chapter 18), the theory of ecological resilience and regime shifts (Chapter 19), ecological network analysis as a measure of economic-ecological health (Chapter 20), the optimization model of circular economy design (Chapter 21), and the thermodynamic foundations that bound all economic activity (Chapter 22). The central result, assembled across these chapters, is a formal proof of the Regeneration Condition: the set of biophysical constraints that any indefinitely sustainable economy must satisfy.

Part IV changes the register of the book. The earlier parts drew primarily from economics, game theory, and network science. Part IV draws additionally from ecology, thermodynamics, and Earth system science — each borrowed with care, translated explicitly into economic terms, and deployed with full acknowledgment of the limits of the analogy. The reader who has worked through Books 1 and 2 will recognize the analytical posture: we follow the evidence, state the assumptions explicitly, and resist the temptation to claim more than the models can deliver.

The stewardship constraint $\dot{N} \geq 0$ — introduced in Chapter 2 and formalized in the Stewardship Objective Function — is the thread connecting all of Part IV. Each chapter adds a dimension to its content: Chapter 17 identifies what $N$ includes; Chapter 18 shows how to account for it; Chapter 19 shows when it exhibits dangerous dynamics; Chapter 20 shows how to measure its health; Chapter 21 shows how to minimize its depletion; and Chapter 22 shows the thermodynamic limits that constrain it absolutely.

---

# Chapter 17: The Economic System as a Subsystem of the Biosphere — Material and Energy Flows

> *"The economy is a wholly owned subsidiary of the environment, not the other way around."*
> — Herman Daly, *Beyond Growth* (1996)

> *"We do not inherit the earth from our ancestors; we borrow it from our children."*
> — attributed to Antoine de Saint-Exupéry and others; provenance uncertain, truth uncontested

## Learning Objectives

By the end of this chapter, you should be able to:

1. Construct the augmented circular flow of income that incorporates material extraction, energy flows, and waste generation, identifying the biophysical stocks and flows that the standard circular flow omits.
2. Apply material flow analysis (MFA) to compute the Total Material Requirement (TMR) and Domestic Material Consumption (DMC) of a national economy, and verify the material balance principle.
3. Formalize the nine Planetary Boundaries as a system of inequality constraints on economic activity, identify which have been crossed, and derive the economic implications of boundary violations.
4. Specify the Global Thresholds and Allocations framework formally as a two-stage decision problem and solve the allocation stage as a linear program for a specified fairness criterion.
5. Derive the natural capital dynamics equation $\dot{N} = \mathcal{R}(N) - \mathcal{D}(C)$ and prove that the Stewardship Condition $\dot{N} \geq 0$ is a necessary (though not sufficient) condition for indefinitely sustainable economic production.
6. Apply the full MFA framework to the Finnish economy and assess it against the Planetary Boundaries.

---

## 17.1 The Standard Circular Flow and What It Omits

The standard circular flow of income, as presented in introductory macroeconomics [P:Ch.3], describes an economy as a closed system of transactions between households and firms: households supply labor and capital to firms, receive income in return, and spend that income on goods and services produced by firms. The circular flow is a useful representation of the monetary economy — the flows of income, expenditure, and value added that national accounts measure.

It is a dangerously incomplete representation of the physical economy.

Two categories of physical flow are entirely absent from the standard circular flow. First, the flow of materials and energy from the natural environment into the production process: firms do not create the atoms they use; they extract them from natural stocks — mineral deposits, agricultural soils, forests, fisheries, aquifers, and the atmosphere. Second, the flow of waste and emissions from the production process back to the natural environment: firms and households do not make matter disappear; they transform it into lower-quality forms — heat, greenhouse gases, chemical pollutants, degraded soils — that return to natural sinks.

These omissions are not merely an accounting inconvenience. They are the source of the central pathology of modern economic systems: the economy can appear to grow — more output, more income, more welfare by any monetary measure — while simultaneously depleting the biophysical substrate on which all future production depends. An economy that liquidates its natural capital while recording rising GDP is analogous to a household that finances current consumption by drawing down savings while recording rising income: the flow statement looks healthy while the balance sheet deteriorates. Chapter 18 develops the accounting framework that captures this balance sheet dimension; this chapter develops the physical measurement tools — material flow analysis and the Planetary Boundaries framework — that tell us what the natural capital stocks are and how fast they are being depleted.

---

## 17.2 The Augmented Circular Flow

**Definition 17.1 (Augmented Circular Flow).** The augmented circular flow extends the standard model by adding three biophysical flows:

1. **Resource extraction flow** $E(t)$: the rate at which the economy extracts materials and energy from natural stocks, measured in mass (tonnes per year) or energy units (petajoules per year).

2. **Waste generation flow** $W(t)$: the rate at which the economy returns materials to natural sinks, in degraded forms — greenhouse gas emissions, wastewater, solid waste, chemical pollutants, heat — measured in equivalent mass units.

3. **Ecosystem service flow** $S(t)$: the rate at which natural systems provide services directly to human welfare without passing through the market — climate regulation, water purification, pollination, nutrient cycling, flood control, aesthetic value — measured in monetary equivalents or biophysical service units.

**The Material Balance Principle.** Conservation of mass in a closed system states that all matter extracted must either remain in the economy as productive stocks or return to the environment as waste:

$$E(t) = \Delta K_{\text{material}}(t) + W(t)$$

where $\Delta K_{\text{material}}$ is the net addition to the economy's material capital stock (buildings, infrastructure, equipment, consumer durables) and $W(t)$ is the waste flow. This is the first law of thermodynamics applied to materials: matter is neither created nor destroyed. At steady state ($\Delta K_{\text{material}} = 0$): $E = W$ — extraction equals waste generation exactly.

**Proposition 17.1 (Waste Cannot Be Eliminated).** In any physical economy at steady state, total waste generation $W = E > 0$ as long as any positive material extraction occurs. Waste can be reduced by reducing extraction, by improving material efficiency (reducing extraction for a given level of output), or by cycling materials within the economy — but it cannot be eliminated without reducing material extraction to zero.

*Proof.* From the material balance principle at steady state: $W = E - \Delta K_{\text{material}} = E - 0 = E$. Since $E > 0$ for any economy with positive material production, $W > 0$. $\square$

The corollary is that "zero waste" strategies, while valuable for reducing $W$, cannot achieve their stated goal without first achieving zero extraction — which in practice means radical reductions in material throughput and/or perfect recycling (zero dissipative losses), both of which face thermodynamic limits analyzed in Chapter 22.

---

## 17.3 Material Flow Analysis

### 17.3.1 The MFA Framework

Material Flow Analysis (MFA) is a systematic accounting framework that tracks the physical flows of materials through an economy, consistent with the material balance principle. We present the formal framework following Eurostat (2001) and the System of Environmental-Economic Accounting (SEEA).

**Definition 17.2 (Total Material Requirement).** The Total Material Requirement (TMR) of an economy is:

$$\text{TMR} = \text{DE} + \text{RDE} + \text{ID} + \text{TRDE}$$

where:
- $\text{DE}$ = Domestic Extraction (materials physically extracted within the territory)
- $\text{RDE}$ = Rucksack of Domestic Extraction (hidden material flows associated with domestic extraction — overburden from mining, processing wastes, etc.)
- $\text{ID}$ = Imports (mass of imported goods)
- $\text{TRDE}$ = Total Rucksack of Imports (hidden flows associated with imports — materials used in production in exporting countries)

**Definition 17.3 (Domestic Material Consumption).** The Domestic Material Consumption (DMC) is:

$$\text{DMC} = \text{DE} + \text{ID} - \text{EX}$$

where $\text{EX}$ = Exports (mass of exported goods). DMC measures the materials actually consumed within the territory, excluding exports and the hidden flows associated with both imports and domestic extraction.

**Definition 17.4 (Material Productivity).** Material productivity (resource efficiency) is:

$$\text{MP} = \frac{\text{GDP}}{\text{DMC}}$$

measured in USD per tonne of material consumed. Higher material productivity indicates more economic value generated per unit of material input — the key indicator of decoupling between economic activity and material throughput.

### 17.3.2 The Material Balance at National Scale

For a national economy, the material balance principle takes the form:

$$\text{DE} + \text{ID} = \text{DMC} + \text{EX} + \text{TW}$$

where $\text{TW}$ = Total Waste (including emissions, water releases, and solid waste returned to the environment). All terms are measured in mass units (typically million tonnes per year).

This identity is the physical analogue of the national income accounting identity [P:Ch.3]: it says that everything that enters the economy (domestic extraction plus imports) must either remain in the economy (DMC net of exports) or return to the environment as waste. Violating this identity implies a measurement error; confirming it validates the accounting.

### 17.3.3 Extraction by Material Category

**Definition 17.5 (Material Categories).** Domestic extraction is disaggregated into four categories:

1. **Biomass** ($E_B$): agricultural crops, grasses and fodder, wood, fish and other aquatic biomass, hunted and gathered biomass.
2. **Fossil fuels** ($E_F$): coal, crude oil, natural gas, peat, and oil shales.
3. **Metal ores** ($E_M$): iron ores, non-ferrous metal ores, precious metal ores, uranium and thorium ores.
4. **Non-metallic minerals** ($E_N$): limestone, gypsum, gravel, sand, clay, salt, and other construction minerals.

$$\text{DE} = E_B + E_F + E_M + E_N$$

Each category has different ecological implications: biomass extraction is renewable (if not exceeding regeneration rates), fossil fuel extraction is non-renewable and generates greenhouse gas emissions, metal ore extraction leaves permanent mine waste, and non-metallic mineral extraction causes landscape transformation. The category breakdown is therefore essential for assessing the ecological sustainability of the material flow, not merely its volume.

---

## 17.4 The Planetary Boundaries

### 17.4.1 Nine Boundaries, Nine Constraints

The Planetary Boundaries framework (Rockström et al., 2009; updated in Steffen et al., 2015 and Richardson et al., 2023) identifies nine Earth system processes that, when destabilized by human activity beyond critical thresholds, risk triggering non-linear, potentially irreversible changes in Earth system functioning.

**Definition 17.6 (Planetary Boundaries).** The nine Planetary Boundaries are formal inequality constraints on economic activity:

$$b_i(\mathbf{x}) \leq \bar{b}_i, \quad i = 1, 2, \ldots, 9$$

where $b_i(\mathbf{x})$ is the control variable for boundary $i$ (a measurable biophysical quantity) as a function of economic activity $\mathbf{x}$, and $\bar{b}_i$ is the boundary value above which the risk of dangerous destabilization becomes unacceptably high.

| $i$ | Earth system process | Control variable $b_i$ | Boundary $\bar{b}_i$ | 2023 status |
|:---|:---|:---|:---|:---|
| 1 | Climate change | Atmospheric CO₂ (ppm) | 350 ppm | **Crossed** (421 ppm) |
| 2 | Biosphere integrity | Extinction rate (species/Myr) | 10 E/MSY | **Crossed** (~100–1000 E/MSY) |
| 3 | Land-system change | Forest cover (% of original) | 75% | **Crossed** (~60%) |
| 4 | Freshwater change | Green/blue water flow alteration | Various | **Crossed** (green water) |
| 5 | Biogeochemical flows | P and N flows (Tg/yr) | P: 11, N: 62 | **Crossed** (P: 22, N: ~150) |
| 6 | Ocean acidification | Carbonate ion saturation ($\Omega$) | 80% of pre-industrial | Approaching |
| 7 | Atmospheric aerosols | Aerosol optical depth | Regional thresholds | Regionally crossed |
| 8 | Stratospheric ozone | O₃ concentration (DU) | 276 DU | Recovering (Montreal Protocol) |
| 9 | Novel entities | Functional diversity loss | Not yet quantified | Unknown/likely crossed |

As of the 2023 assessment by Richardson et al., six of the nine boundaries have been crossed — a number that has increased from three at the time of the original 2009 assessment.

### 17.4.2 Economic Implications of Boundary Violations

Each crossed boundary has direct economic implications.

**Climate change (boundary 1).** The DICE model [P:Ch.37] estimates the social cost of carbon at USD 50–200 per tonne of CO₂ under standard discount rates, rising to USD 500+ under near-zero discount rates. At 421 ppm CO₂ (71 ppm above the boundary), the implied liability of current atmospheric CO₂ above the safe boundary is approximately USD 50–200 trillion globally — a stock liability that does not appear on any national balance sheet.

**Biosphere integrity (boundary 2).** Ecosystem services from biodiversity — pollination, pest control, water purification, disease regulation — have been estimated at USD 125–145 trillion per year (Costanza et al., 2014), approximately 1.5× global GDP. Extinction rates 10–100× above background rates imply rapid erosion of this service base, with non-linear risks as keystone species are lost and ecosystem functions collapse.

**Biogeochemical flows (boundary 5).** Nitrogen and phosphorus loading in excess of planetary boundaries creates dead zones in aquatic ecosystems, acidifies freshwater, and contributes to greenhouse gas emissions (nitrous oxide, $N_2O$, is 298× more potent than CO₂ over 100 years). The economic cost of hypoxic dead zones alone has been estimated at several billion USD per year in lost fisheries and water treatment costs.

**Formal representation.** The Planetary Boundaries form a constraint set $\mathcal{B} \subset \mathbb{R}^9$ in the space of biophysical control variables. Any sustainable economic trajectory must satisfy:

$$\mathbf{b}(\mathbf{x}(t)) \in \mathcal{B} \quad \forall t \geq 0$$

where $\mathbf{b}: \mathbb{R}^n \to \mathbb{R}^9$ maps economic activity vector $\mathbf{x}$ to biophysical indicators. Economic growth that pushes any $b_i$ above $\bar{b}_i$ is formally unsustainable — the Stewardship Condition is violated for the corresponding natural capital stock.

---

## 17.5 Global Thresholds and Allocations

### 17.5.1 The Two-Stage Framework

The Planetary Boundaries framework identifies what cannot be exceeded globally; it does not specify who has the right to use the remaining biophysical capacity. The Global Thresholds and Allocations (GTA) framework addresses this second question.

**Definition 17.7 (Global Thresholds and Allocations Framework).** The GTA framework is a two-stage decision problem:

**Stage 1 (Thresholds — scientific):** Determine the maximum globally permissible level of activity $\bar{X}_i$ for each Planetary Boundary $i$, consistent with $b_i(\bar{X}_i) \leq \bar{b}_i$. This is a scientific question, answered by Earth system scientists, not by economic optimization.

**Stage 2 (Allocations — normative):** Distribute the permissible activity $\bar{X}_i$ across nations $j = 1, \ldots, n$, sectors $s = 1, \ldots, m$, and generations $g = 1, 2, \ldots$ according to a fairness criterion:

$$\sum_j \sum_s x_{ij,s} \leq \bar{X}_i \quad \forall i \quad \text{(planetary constraint)}$$
$$x_{ij,s} \geq 0 \quad \forall i, j, s$$

with an objective function $f(\{x_{ij,s}\})$ expressing the chosen fairness criterion.

**Fairness criteria.** Three standard fairness criteria generate different allocations:

1. **Equal per-capita:** $x_{ij,s} = \bar{X}_i \cdot \frac{P_j}{P_{\text{global}}} \cdot \frac{1}{m}$ — each person receives an equal share of the global budget, regardless of current emissions or economic status.

2. **Grandfathering:** $x_{ij,s} \propto$ historical emissions — countries that have historically emitted more retain more of the remaining budget.

3. **Capability-adjusted:** $x_{ij,s}$ accounts for both per-capita equity and economic capacity — richer countries reduce faster, poorer countries retain larger per-capita budgets.

**Proposition 17.2 (No Allocation is Fair to All Criteria).** There exists no allocation $\{x_{ij,s}\}$ satisfying all three fairness criteria simultaneously, for any configuration of national populations and emissions histories in which historical emitters differ from per-capita leaders.

*Proof.* Consider a two-country case (North, South) with North having 20\% of global population but 60\% of historical emissions. Equal per-capita gives North 20\% of the global budget; grandfathering gives North 60\%; capability-adjusted gives North somewhere between, weighted by income. These three allocations cannot simultaneously be equal — they represent genuinely conflicting normative principles. By induction, the same holds for $n > 2$ countries. $\square$

The allocation problem is therefore inherently political, not merely technical. The science determines $\bar{X}_i$ (Stage 1); the distribution of $\bar{X}_i$ across nations and generations requires collective political decisions that cannot be derived from the science alone. This is the formal basis for the claim, in the UNFCCC and Paris Agreement process, that equity and ambition must be decided together.

### 17.5.2 The Allocation LP

For a given fairness criterion, the allocation problem can be formulated as a linear program. We present the equal per-capita formulation with a budget flexibility term:

$$\min_{\{x_{ij}\}} \sum_j c_j \cdot \left|x_{ij} - \frac{P_j}{P} \cdot \bar{X}_i\right|$$

$$\text{subject to:} \quad \sum_j x_{ij} \leq \bar{X}_i, \quad x_{ij} \geq x_{ij}^{\min}$$

where $c_j$ is the cost weight for country $j$ deviating from the equal per-capita allocation (reflecting political feasibility constraints) and $x_{ij}^{\min}$ is the minimum subsistence allocation for country $j$. The $L^1$ objective (absolute deviation from per-capita fairness) makes this a linear program in the absolute value sense, solvable by standard LP methods [M:Ch.1].

---

## 17.6 Natural Capital Dynamics

### 17.6.1 The Dynamics Equation

Natural capital $N(t)$ encompasses the stocks of ecological assets that provide the material and energetic inputs to economic production and the waste absorption capacity that receives its outputs. We model $N$ as a vector $\mathbf{N} = (N_1, N_2, \ldots, N_k)$ of $k$ distinct natural capital stocks (soils, forests, fisheries, groundwater, atmospheric capacity, biodiversity, etc.).

**Definition 17.8 (Natural Capital Dynamics).** Each natural capital stock $N_j$ evolves according to:

$$\dot{N}_j = \mathcal{R}_j(N_j, \mathbf{N}_{-j}) - \mathcal{D}_j(C, E_j)$$

where:
- $\mathcal{R}_j(N_j, \mathbf{N}_{-j})$: the regeneration function — the rate at which stock $j$ recovers through natural processes, potentially depending on other stocks $\mathbf{N}_{-j}$ (ecological interdependencies).
- $\mathcal{D}_j(C, E_j)$: the depletion function — the rate at which stock $j$ is drawn down by consumption $C$ and direct extraction $E_j$.

For the logistic regeneration model (appropriate for renewable biological stocks):

$$\mathcal{R}_j(N_j) = r_j N_j\left(1 - \frac{N_j}{K_j}\right)$$

where $r_j$ is the intrinsic growth rate and $K_j$ is the carrying capacity of stock $j$.

For fossil fuels and non-renewable minerals: $\mathcal{R}_j = 0$ — no natural regeneration on human timescales.

**Definition 17.9 (Stewardship Condition).** The Stewardship Condition for natural capital stock $j$ is:

$$\dot{N}_j \geq 0 \quad \Leftrightarrow \quad \mathcal{R}_j(N_j, \mathbf{N}_{-j}) \geq \mathcal{D}_j(C, E_j)$$

The aggregate Stewardship Condition requires:

$$\sum_j \mu_j \dot{N}_j \geq 0$$

where $\mu_j > 0$ are shadow prices of natural capital stocks — the condition that the value-weighted natural capital stock is non-declining (weak sustainability).

**Theorem 17.1 (Stewardship Condition as Necessary Condition for Sustainable Production).** If the Stewardship Condition is violated for any essential natural capital stock $j$ — i.e., if $N_j(t) \to 0$ as $t \to \infty$ — then economic production $Y(t) \to 0$ as $t \to \infty$, provided $j$ is an essential input to the production function.

*Proof.* Let $Y = F(K, N_j, L, \ldots)$ be the production function with $\partial F/\partial N_j > 0$ (natural capital is productive) and $F(K, 0, L, \ldots) = 0$ (natural capital is essential — production is impossible without it). If $\dot{N}_j < 0$ for all $t$ sufficiently large, $N_j(t) \to 0$ in finite or infinite time. Since $F$ is continuous and $F(K, 0, L) = 0$, $Y(t) = F(K, N_j(t), L, \ldots) \to 0$. $\square$

The proof is simple but the implication is profound: an economy that depletes its essential natural capital stocks is not merely reducing future welfare — it is setting itself on a trajectory toward zero production. This is not a claim about climate preferences or environmental values; it is a claim about physical production possibility.

---

## 17.7 Mathematical Model: The Formal MFA System

We now develop the complete formal MFA model, integrating material balance, Planetary Boundaries, and natural capital dynamics into a unified system.

**State variables.** Let:
- $\mathbf{N}(t) = (N_1, \ldots, N_k)$: natural capital stocks (tonnes or ecological units)
- $\mathbf{K}(t)$: produced capital stocks (monetary units)
- $\mathbf{E}(t) = (E_1, \ldots, E_k)$: extraction flows from each stock (tonnes/year)
- $\mathbf{W}(t)$: waste and emission flows (tonnes/year)
- $\mathbf{b}(t) = (b_1, \ldots, b_9)$: Planetary Boundary control variables

**Dynamics:**

$$\dot{N}_j = \mathcal{R}_j(N_j) - E_j(t) \quad \forall j \tag{17.1}$$

$$\dot{K} = F(\mathbf{K}, \mathbf{N}, L) - C - \delta_K K \quad \tag{17.2}$$

$$\mathbf{b}(t) = \mathcal{B}(\mathbf{E}(t), \mathbf{W}(t), \mathbf{N}(t)) \quad \tag{17.3}$$

$$\mathbf{b}(t) \leq \bar{\mathbf{b}} \quad \text{(Planetary Boundaries constraint)} \tag{17.4}$$

**The material balance:**

$$\sum_j E_j(t) = \Delta K_{\text{material}}(t) + \sum_j W_j(t) \quad \tag{17.5}$$

This system couples the economic dynamics (17.2) to the ecological dynamics (17.1) through extraction flows $\mathbf{E}$, subject to the Planetary Boundaries constraints (17.4) and the material balance (17.5). The system is the formal representation of the augmented circular flow of Section 17.2.

**The Stewardship Program.** The socially optimal extraction path solves the Stewardship Objective Function [C:Ch.2, Definition 2.7] subject to the Planetary Boundaries constraints:

$$\max_{\{\mathbf{E}(t), C(t)\}} \int_0^\infty e^{-\rho t} U(C(t), S(t))\,dt$$

$$\text{subject to: } (17.1)-(17.5), \quad \dot{N}_j \geq 0 \; \forall j, \quad \mathbf{b}(t) \leq \bar{\mathbf{b}}$$

This is a constrained optimal control problem [Appendix A] with state variables $(\mathbf{N}, \mathbf{K})$, control variables $(\mathbf{E}, C)$, and two types of constraints: the stock maintenance constraints ($\dot{N}_j \geq 0$) and the Planetary Boundaries constraints ($\mathbf{b}(t) \leq \bar{\mathbf{b}}$). We develop the full solution in Chapters 18–20; the present chapter establishes the problem.

---

## 17.8 Worked Example: MFA for the Finnish Economy

Finland provides a useful worked example because its national statistical office (Statistics Finland) publishes detailed material flow accounts consistent with Eurostat standards, and because Finland's combination of strong natural capital (forests, freshwater, minerals) and high economic development creates interesting tensions between growth and stewardship.

### 17.8.1 Material Flows: 2020 Data

We use Statistics Finland's MFA data for 2020 (the most recent complete year at time of writing). All flows in million tonnes (Mt) per year.

**Domestic Extraction by category:**
- Biomass: $E_B = 88$ Mt (dominated by wood, food crops, grasses)
- Fossil fuels: $E_F = 2$ Mt (Finland has minimal domestic fossil fuel reserves)
- Metal ores: $E_M = 45$ Mt (significant mining sector: nickel, copper, chromite)
- Non-metallic minerals: $E_N = 90$ Mt (construction materials dominate)

$$\text{DE} = 88 + 2 + 45 + 90 = 225 \text{ Mt}$$

**Imports:** ID = 62 Mt (Finland imports oil, manufactured goods, food)
**Exports:** EX = 78 Mt (paper, pulp, metals, machinery)

**DMC:**
$$\text{DMC} = 225 + 62 - 78 = 209 \text{ Mt}$$

**Per capita DMC:** Finland's 2020 population was 5.53 million. $\text{DMC per capita} = 209/5.53 = 37.8$ tonnes/person/year.

For comparison: EU average DMC per capita ≈ 14 tonnes/person/year; global average ≈ 12 tonnes/person/year. Finland's per capita DMC is approximately 3× the EU average, driven primarily by the large mining and forestry sectors.

**Material productivity:**
$$\text{MP} = \frac{\text{GDP}}{\text{DMC}} = \frac{\text{EUR } 236\text{bn}}{209\text{ Mt}} = 1{,}130 \text{ EUR/tonne}$$

EU average material productivity ≈ 2,400 EUR/tonne. Finland's material productivity is below the EU average, reflecting the material intensity of its natural-resource-based economy.

### 17.8.2 Planetary Boundary Assessment

We assess Finland against each Planetary Boundary for which national-level data is available.

**Climate change (boundary 1).** Finland's total greenhouse gas emissions in 2020: 48.4 Mt CO₂-equivalent. Per capita: 8.8 t CO₂e/person/year. A fair-share equal per-capita budget consistent with the 1.5°C pathway (IPCC AR6): approximately 2.5 t CO₂e/person/year by 2030. Finland currently emits 3.5× its fair-share budget. Status: **Crossed**.

**Biogeochemical flows (boundary 5).** Finland's reactive nitrogen surplus in agricultural soils: approximately 8 kg N/ha/year (below the European average of ~25 kg/ha/year but above the planetary boundary for nitrogen flow per capita). Phosphorus: Finland has implemented effective phosphorus capture in wastewater treatment; agricultural P surplus ≈ 2 kg P/ha/year, close to the boundary. Status: **Approaching/borderline**.

**Land-system change (boundary 3).** Finland maintains approximately 75\% of its territory under forest cover — meeting the Planetary Boundary threshold locally. However, net deforestation in Finland's managed forests (harvest minus growth) has been negative in recent years (more harvested than growing), raising questions about whether the stock is truly maintained. Status: **Nominally meeting, but declining**.

**Freshwater change (boundary 4).** Finland has abundant freshwater resources; annual freshwater withdrawal is approximately 15\% of renewable freshwater availability. Status: **Within boundary**.

### 17.8.3 Throughput Reduction Required for Compliance

To estimate the reduction in economic throughput required for Finland to comply with its fair-share Planetary Boundary allocations, we solve the inverse problem: given current technology (material productivity MP), what level of GDP is consistent with the Planetary Boundaries?

**For climate change:** Fair-share CO₂ budget = 2.5 t/person/year × 5.53 million = 13.8 Mt CO₂e/year. At Finland's current GDP-CO₂ intensity of 0.205 kg CO₂e/EUR: maximum GDP consistent with climate boundary = 13.8 × 10⁶ t / (0.205 kg/EUR × 10⁻³ t/kg) = EUR 67 billion — approximately **72\% reduction** from 2020 GDP (EUR 236 billion) at current intensity.

This is not, however, the only feasible path. Decarbonization raises the denominator: if Finland achieves the EU's projected 2030 carbon intensity of 0.040 kg CO₂e/EUR, the compliant GDP level rises to EUR 345 billion — above current GDP. The throughput reduction required depends critically on whether decarbonization proceeds on schedule.

**For material flows:** At the EU average material productivity of 2,400 EUR/tonne, Finland's current GDP of EUR 236 billion would require DMC of 98 Mt — a 53\% reduction from actual DMC of 209 Mt. Much of this reduction would come from shifting from primary extraction to recycled inputs in the mining and construction sectors.

**Summary:** Finland's material and carbon metabolism significantly exceeds its planetary fair-share allocation. Meeting the Planetary Boundaries requires either large reductions in throughput or large improvements in material and energy productivity — ideally both. The magnitude of the challenge is substantial but not unprecedented: Finland's per-capita emissions fell by approximately 40\% between 1990 and 2020 through energy system transformation, demonstrating that significant decoupling is achievable.

---

## 17.9 Case Study: The Doughnut Economics Framework

### 17.9.1 The Framework

Kate Raworth's Doughnut Economics (2017) provides a framework that combines the Planetary Boundaries ceiling with a social foundation floor — a set of minimum social conditions (food, water, health, education, income, political voice, social equity, gender equality, housing, networks, energy, democratic governance) that every person must be able to access. The resulting "doughnut" shape defines a safe and just space for humanity between the social floor (the inner ring) and the Planetary Boundaries ceiling (the outer ring).

**Formal specification.** Define the Doughnut as the feasible set:

$$\mathcal{D} = \{\mathbf{x} \in \mathbb{R}^n : \mathbf{b}(\mathbf{x}) \leq \bar{\mathbf{b}} \text{ and } \mathbf{s}(\mathbf{x}) \geq \underline{\mathbf{s}}\}$$

where $\mathbf{s}(\mathbf{x}) = (s_1, \ldots, s_{12})$ is the vector of social indicators (food, water, health, etc.) and $\underline{\mathbf{s}}$ is the social floor — the minimum acceptable level of each indicator.

An economy at activity vector $\mathbf{x}$ is:
- In the "safe and just space" if $\mathbf{x} \in \mathcal{D}$.
- "Overshooting" the ecological ceiling if $b_i(\mathbf{x}) > \bar{b}_i$ for any $i$ (too much environmental pressure).
- "Falling short" of the social foundation if $s_j(\mathbf{x}) < \underline{s}_j$ for any $j$ (too little social provision).

Most countries today fall into one of three categories: (1) high-income countries that meet the social foundation but overshoot the ecological ceiling (most OECD); (2) low-income countries that meet the ecological ceiling but fall short of the social foundation (most sub-Saharan Africa); (3) middle-income countries that partially meet both or partially violate both. No country currently meets both conditions simultaneously in all dimensions.

### 17.9.2 Empirical Assessment

A comprehensive empirical assessment of the Doughnut for 150 countries (O'Neill et al., 2018, updated for 2020 data) finds:

- **No country** achieves all 12 social indicators above the floor while meeting all 7 ecological indicators below the ceiling.
- **High-income countries** (OECD) score well on social indicators (median social shortfall: 0.2/12 dimensions) but consistently overshoot ecological ceilings (median ecological overshoot: 5.1/7 dimensions).
- **Low-income countries** (LICs) rarely overshoot ecological ceilings (median: 0.8/7) but consistently fall short on social indicators (median: 7.3/12).
- The **cross-over point** — countries that come closest to meeting both simultaneously — tends to occur at a level of biophysical throughput roughly 1.5–2× above the fair-share planetary allocation. This is the quantitative expression of the development dilemma: current technologies do not allow universal achievement of the social floor within the ecological ceiling.

**The formal implication.** The Doughnut framework formalizes what the Stewardship Objective Function introduced in Chapter 2 implies: welfare is not a function of consumption alone but of provisioning services $S(t)$, which depend on both produced capital and natural capital. Achieving the social floor requires adequate provisioning services; respecting the ecological ceiling requires maintaining natural capital. The Doughnut's safe and just space $\mathcal{D}$ is the feasible set of the Stewardship Program — the set of economic activity vectors that simultaneously maintain natural capital (ceiling compliance) and generate sufficient provisioning services (floor compliance).

The technical challenge of the 21st century is to design economic systems that can operate within $\mathcal{D}$: that achieve universal social provision without ecological overshoot. The remaining chapters of Part IV and Parts V–VIII develop the theory and practice of such systems.

---

## Chapter Summary

This chapter has established the biophysical foundations of the economics of cooperation.

The standard circular flow of income omits two categories of physical flow that are essential for sustainability analysis: resource extraction from natural stocks and waste generation to natural sinks. The material balance principle — conservation of mass — implies that these flows are inseparable: everything extracted must either accumulate in the economy or return to the environment as waste. "Zero waste" is achievable only by achieving zero extraction.

Material Flow Analysis provides the formal accounting framework for tracking physical flows through an economy. Total Material Requirement (TMR) captures the full material footprint including hidden flows; Domestic Material Consumption (DMC) captures the materials actually consumed within a territory; material productivity (GDP/DMC) measures the efficiency of material use.

The Planetary Boundaries framework formalizes nine Earth system constraints as inequality conditions $b_i(\mathbf{x}) \leq \bar{b}_i$. Six of nine boundaries have been crossed as of the 2023 assessment, implying that the current global economic trajectory violates the Stewardship Condition for multiple essential natural capital stocks. The Global Thresholds and Allocations framework separates the scientific question (what are the thresholds?) from the normative question (how should the remaining biophysical capacity be allocated?), with three fairness criteria generating three different allocation rules — none universally satisfactory.

The natural capital dynamics equation $\dot{N}_j = \mathcal{R}_j(N_j) - \mathcal{D}_j(C, E_j)$ formalizes the Stewardship Condition: for sustainable production, regeneration must at least match depletion for every essential natural capital stock. Theorem 17.1 proves this is a necessary condition: economies that deplete essential natural capital set themselves on a trajectory toward zero production.

The Finnish MFA demonstrates the methodology: a high-income country with strong natural capital endowments nonetheless faces substantial overshots of its fair-share Planetary Boundary allocations, requiring either deep decarbonization or reduced material throughput to return to the safe-and-just space of the Doughnut framework.

Chapter 18 develops the accounting framework needed to track natural capital alongside produced capital in a stock-flow consistent model — providing the balance sheet perspective that the flow analysis of this chapter requires for completeness.

---

## Exercises

**17.1** State the material balance principle formally (Proposition 17.1). Apply it to the following economy:
- Domestic extraction: 150 Mt biomass, 30 Mt fossil fuels, 20 Mt metal ores, 80 Mt non-metallic minerals.
- Imports: 45 Mt; Exports: 60 Mt.
- Addition to material capital stocks ($\Delta K_{\text{material}}$): 35 Mt.

(a) Compute DMC and TMR.
(b) Compute total waste generation $W$ from the material balance.
(c) If 60\% of waste consists of CO₂ emissions from fossil fuel combustion, compute the CO₂ emissions from this economy and compare to the Planetary Boundary threshold.

**17.2** Six of the nine Planetary Boundaries have been crossed as of 2023. For the three boundaries nearest to being crossed (ocean acidification, atmospheric aerosols, novel entities):
(a) Describe the Earth system process governed by each boundary and the primary economic activities that push it.
(b) For ocean acidification: the control variable is seawater carbonate ion saturation. At current atmospheric CO₂ of 421 ppm, the saturation level is approximately 90\% of pre-industrial. The boundary is 80\%. How much additional CO₂ can the atmosphere absorb before the boundary is crossed? Express in Gt CO₂.
(c) Why is the novel entities boundary (boundary 9) particularly difficult to quantify? What measurement challenges prevent its formal specification as an inequality constraint?

**17.3** Apply the equal per-capita allocation rule from Section 17.5.2 to the OECD's 2020 CO₂ emissions data.
(a) The Paris Agreement 1.5°C budget remaining from 2020 is approximately 400 Gt CO₂. Distributed equally across the 2020 global population of 7.8 billion, compute the per-capita budget.
(b) For three representative countries — USA ($\approx$4.5 t CO₂/person/year in 2020), EU average ($\approx$7.3 t/person/year), India ($\approx$1.9 t/person/year) — compute each country's budget under the equal per-capita rule and compare to current emissions. Who must reduce most?
(c) Formulate the allocation LP of Section 17.5.2 for these three countries with a capability weight $c_j$ inversely proportional to per-capita income. How does the capability-adjusted allocation differ from equal per-capita?

**★ 17.4** Using the natural capital dynamics model (Definition 17.8) with logistic regeneration, derive the condition for minimum sustainable extraction.

(a) For stock $j$ with $\mathcal{R}_j(N_j) = r_j N_j(1-N_j/K_j)$ and constant extraction $E_j$, find the steady-state stock $N_j^*$ as a function of $E_j$.
(b) Show that the maximum sustainable yield (MSY) is $E_j^{\text{MSY}} = r_j K_j / 4$ at $N_j^* = K_j/2$.
(c) Derive the condition on consumption $C$ and extraction efficiency $\eta$ such that $E_j(C, \eta) \leq E_j^{\text{MSY}}$ — the formal statement of the Stewardship Condition in terms of economic variables.
(d) How does the condition change if ecological interdependencies are present? Let $\mathcal{R}_j(N_j, N_l) = r_j N_j(1-N_j/K_j) + \gamma_{jl} N_l$ where $\gamma_{jl}$ is the positive externality from stock $l$ on stock $j$'s regeneration. Derive the new MSY and show it is higher than the single-stock MSY when $\gamma_{jl} > 0$.

**★ 17.5** Prove Theorem 17.1 (Stewardship Condition as necessary condition) under more general conditions.

(a) Prove the theorem for a production function of the form $Y = F(K, N_j, L)$ satisfying the Inada conditions: $\partial F/\partial N_j \to \infty$ as $N_j \to 0$ (natural capital becomes infinitely productive in scarcity). Show that not only does $Y \to 0$ as $N_j \to 0$, but the rate of convergence is faster than linear.
(b) Now consider weak sustainability: suppose produced capital $K$ can substitute for natural capital $N_j$ (Solow-Hartwick rule [P:Ch.5]). Under what conditions does the production function $Y = F(K, N_j, L)$ allow $Y$ to remain positive as $N_j \to 0$? What elasticity of substitution $\sigma$ is required?
(c) The coral reef ecosystem provides a concrete example of a natural capital stock for which the Inada condition approximately holds: tourism revenue, fisheries, and coastal protection all decline sharply as coral cover falls toward zero. Calibrate a simple version of the model to the Maldives case and estimate the critical reef cover threshold below which economic production collapses.

**★★ 17.6** Using publicly available data from the World Bank SEEA (System of Environmental-Economic Accounting) and Eurostat's Material Flow Accounts, conduct a full Doughnut assessment for two countries: one high-income (your choice from OECD) and one upper-middle-income (your choice from non-OECD).

(a) Compute TMR, DMC, and material productivity for each country for the most recent available year.
(b) Assess each country against all nine Planetary Boundaries for which national-level data is available. For each boundary, report the control variable value, the boundary threshold, and the status (within/crossed/approaching).
(c) Assess each country against the 12 social foundation indicators (using UN SDG data as proxies). For each indicator, report the value, the minimum threshold, and the status (above/below floor).
(d) For each country, identify the minimum set of interventions (in both material throughput and social provision) that would bring them into the Doughnut's safe and just space. Estimate the magnitude of each intervention.
(e) Is there a trajectory from each country's current position to the safe and just space that does not require reducing GDP? What assumptions about decoupling (material and carbon intensity improvement) are required?

---

*Chapter 18 extends the stock-flow consistent accounting framework of Book 2 [M:Ch.28] to incorporate natural capital as a balance sheet stock, developing the Provisioning Balance Sheet that integrates produced capital, natural capital, and social capital into a unified accounting framework — the formal foundation for the stewardship approach to macroeconomic measurement.*
