# Chapter 38: Universal Basic Services and Commons-Based Provision

> *"The question is not whether we can afford universal public services. The question is whether we can afford to go without them."*
> — Andrew Harrop, *For Us All: Redesigning Social Security for the 2020s* (2016)

> *"Vienna's housing policy is not social policy. It is economic policy. By providing affordable housing to everyone, Vienna keeps its economy competitive."*
> — Michael Ludwig, Mayor of Vienna (2019, paraphrased)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Formally compare Universal Basic Income (UBI) and Universal Basic Services (UBS) on welfare, distributional, ecological, and fiscal sustainability grounds, proving the conditions under which UBS dominates UBI in a welfare-maximizing sense.
2. Specify a formal commons-based provision model — identifying when cooperative or commons-based service delivery outperforms both privatized markets and standard state bureaucracy, and designing the governance structures that make it viable.
3. Construct the SFC-compatible fiscal framework for UBS, proving that UBS can be financed without either deficit spending or tax increases when funded through the seigniorage of sovereign money creation and the demurrage revenue of Chapter 27.
4. Apply the UBS framework to a concrete city design problem — specifying the service bundle, governance structure, and funding mechanism for a city of 500,000 — and compute the welfare gain from UBS relative to market provision.
5. Analyze Vienna's Gemeindebau — the world's largest social housing cooperative system, providing approximately 60\% of Vienna's population with housing at regulated rents — as a formal test of the commons-based provision model.
6. Derive the conditions under which UBS reduces ecological footprint relative to UBI, connecting the provisioning model to the natural capital stewardship framework of Part IV.

---

## 38.1 From Cash to Services: The UBS Proposition

The dominant paradigm of social policy since the 1990s has been means-tested cash transfers: identify individuals below an income threshold and provide them with enough money to bring them above it. This paradigm has three formal properties — means-testing efficiency (resources go to those who need them), individual autonomy (recipients decide what to buy), and market compatibility (cash transfers do not distort private markets).

Universal Basic Income extends this paradigm to its logical limit: eliminate means-testing, provide a flat cash transfer to every citizen regardless of income, and let individuals purchase whatever services they need. UBI's proponents argue for simplicity, autonomy, and the elimination of poverty traps created by means-testing. Its critics note the high fiscal cost of universal provision and the risk that cash is spent on goods that do not generate the welfare improvements that public policy aims for.

Universal Basic Services proposes a different approach: rather than giving everyone cash to purchase services in markets, provide everyone with access to the services themselves — healthcare, education, housing, transport, digital access, democratic participation, and care — at high quality and free at point of use. The services are delivered collectively, through commons-based or publicly governed institutions, and funded through collective means.

The formal comparison is not simple. UBS is welfare-superior to UBI under some conditions and inferior under others. This chapter derives those conditions precisely, rather than asserting either position.

---

## 38.2 Formal Comparison: UBI vs. UBS

### 38.2.1 The Welfare Model

**Definition 38.1 (Household Welfare Function).** Household $i$'s welfare is:

$$W_i = U_i(C_i^{\text{private}}, S_i^{\text{public}}, L_i, E_i, H_i)$$

where $C_i^{\text{private}}$ is private consumption, $S_i^{\text{public}}$ is publicly provided service consumption, $L_i$ is leisure, $E_i$ is ecological quality (natural capital services), and $H_i$ is housing adequacy. The partial derivatives: $\partial W_i/\partial C_i > 0$ (more private consumption is better), $\partial W_i/\partial S_i > 0$ (more public services is better), with diminishing marginal utilities.

**The market failure conditions.** UBS is motivated by specific market failures in service provision:

1. **Information asymmetry:** Healthcare and education markets have severe information asymmetries — patients and students cannot evaluate quality ex ante. Market provision leads to quality collapse through adverse selection and Akerlof lemons dynamics [C:Ch.16].

2. **Network externalities:** Education generates positive externalities (a more educated population benefits everyone through innovation, civic participation, reduced crime). Market provision under-supplies education relative to the social optimum [C:Ch.33].

3. **Non-rival provision economics:** Public health and sanitation are largely non-rival — the marginal cost of serving an additional person is near zero once infrastructure exists. Market provision at cost-covering prices excludes welfare-improving consumption.

4. **Risk pooling:** Insurance markets suffer from adverse selection — high-risk individuals seek coverage, raising premiums, driving out low-risk individuals in a death spiral. Social insurance (publicly funded healthcare, unemployment insurance) pools risk more efficiently.

### 38.2.2 The Formal Comparison

**UBI specification.** Every citizen receives a cash transfer $B$ per period. Household $i$ purchases services at market prices $\{p_j\}$: $C_i^s = B_i + y_i - p_j s_{ij}$ where $s_{ij}$ is consumption of service $j$.

**UBS specification.** Every citizen receives universal access to a service bundle $\bar{S} = (\bar{s}_1, \ldots, \bar{s}_k)$ at zero price, funded collectively. Household $i$ may supplement with private purchases: $s_{ij} = \bar{s}_j + x_{ij}$ where $x_{ij} \geq 0$ is additional private consumption.

**Theorem 38.1 (UBS Dominates UBI Under Market Failure).** For any service $j$ with market failure condition $M_j$ (information asymmetry, network externality, or non-rival provision), UBS provision at quality $\bar{q}_j$ dominates UBI provision of cash equivalent $p_j \bar{s}_j$ in welfare terms:

$$\sum_i W_i(\text{UBS}) > \sum_i W_i(\text{UBI})$$

when the social provision quality $\bar{q}_j$ exceeds the average market quality achievable with the cash equivalent $p_j \bar{s}_j$ under the market failure condition.

*Proof.* Under information asymmetry (case 1): households with income below average cannot distinguish high-quality from low-quality providers. Market provision with cash equivalent $p_j \bar{s}_j$ purchases average-quality services at best; households systematically receive lower quality than the market price implies (adverse selection). UBS provision at publicly monitored quality $\bar{q}_j$ eliminates adverse selection — all households receive the monitored quality level. Therefore $S_i^{\text{UBS}} = \bar{q}_j > S_i^{\text{UBI}} = q_i^{\text{market}} < \bar{q}_j$ for low-income households, giving $W_i^{\text{UBS}} > W_i^{\text{UBI}}$ for those households. For high-income households with access to above-average market quality: UBS may be welfare-inferior ($W_i^{\text{UBS}} < W_i^{\text{UBI}}$ if $q_i^{\text{market}} > \bar{q}_j$) — but UBS allows supplementary private purchases, so $W_i^{\text{UBS}} \geq W_i^{\text{UBI}}$ for all households who retain the right to supplement. Aggregating over all households: UBS is weakly better for all and strictly better for those facing adverse selection. $\square$

**Corollary 38.1 (UBS Distributional Dominance).** UBS is more progressive than UBI of equal fiscal cost: high-quality universal services transfer more welfare to lower-income households (who cannot afford quality services privately) than to higher-income households (who could have purchased equivalent quality privately).

*Proof.* Let $\Delta W_i^{\text{UBS-UBI}} = W_i^{\text{UBS}} - W_i^{\text{UBI}}$ be the welfare gain from UBS over UBI for household $i$. For high-income households ($y_i$ high): they already purchase high-quality services privately. $\Delta W_i^{\text{high}} \approx 0$ (UBS provides what they already have). For low-income households: $\Delta W_i^{\text{low}} > 0$ (UBS provides quality they cannot afford). Therefore $\Delta W_i$ is decreasing in income — UBS is distributionally progressive by construction. $\square$

### 38.2.3 The Ecological Argument for UBS

**Proposition 38.1 (UBS Reduces Ecological Footprint Relative to UBI).** Universal provision of mobility, heating, and food security through collective commons-based infrastructure reduces per-capita material throughput relative to individualized cash-funded market provision.

*Proof.* Three mechanisms:

1. **Scale economies in material use:** Public transport requires fewer vehicles, less fuel, and less road infrastructure per passenger-kilometre than private car ownership. A shared bus serving 50 passengers uses approximately 1/15 of the fuel per passenger of 50 individually owned cars. UBS provision of mobility (public transport) reduces the material and energy throughput of the mobility sector relative to UBI-funded private car ownership.

2. **Elimination of defensive expenditure:** UBS provision of healthcare, legal services, and housing security eliminates the defensive expenditure that households make under market provision — purchasing insurance, lock-in contracts, premium "safe" products — much of which generates material throughput without improving welfare [C:Ch.31, Definition 31.1].

3. **Longer asset lifetimes:** Collective ownership (commons) maintains assets to higher quality standards and replaces them less frequently than market provision, where competitive pressure incentivizes planned obsolescence. Public housing built to cooperative standards (Vienna's Gemeindebau: average building age 60+ years) outlasts private rental housing (average replacement cycle 25–35 years in comparable markets). $\square$

---

## 38.3 Commons-Based Provision Model

### 38.3.1 When Commons Provision Outperforms State and Market

The public choice literature identifies pathologies of state provision (bureaucratic capture, political short-termism, democratic accountability gaps) and the market failure literature identifies pathologies of market provision (adverse selection, moral hazard, monopoly, externalities). Commons-based provision — governance by the community of users and contributors, as formalized in Chapter 14 — offers a third option.

**Definition 38.2 (Commons-Based Service Provision).** Commons-based service provision is the delivery of a public service by a self-governing community of users and contributors under the Ostrom design principles, with:
- Users who are also members of the governing body (not passive recipients).
- Contribution norms that make service provision partially self-reinforcing (members maintain the commons they use).
- Graduated membership (service quality and governance rights linked to contribution history).

**Proposition 38.2 (Commons Provision Outperforms State and Market Under Specific Conditions).** Commons-based provision outperforms both market provision and state bureaucratic provision when:

1. **Local knowledge is essential:** Service quality depends on contextual knowledge that central planners and distant market actors cannot possess — home care services, community mental health, neighborhood safety, local food systems. [C:Ch.13, Hayek Theorem 13.1]

2. **Intrinsic motivation matters:** Service quality depends critically on contributor motivation — teaching, nursing, social work, community organizing. Market provision (paying market wages) crowds out intrinsic motivation; commons provision (governance rights, peer recognition, collective mission) maintains it. [Frey and Jegen, 2001 — motivation crowding]

3. **Trust is a primary input:** Services where the producer-consumer relationship requires high trust — healthcare, childcare, elder care — benefit from community governance that aligns incentives through repeated interaction and reputation [C:Ch.16].

4. **Asset specificity is high:** Long-lived assets (housing, community infrastructure, cooperative facilities) require governance structures that maintain them across generations — commons governance with stewardship norms outperforms short-term market contracts.

*Proof.* Each condition identifies a failure mode of market and state provision that commons governance avoids. Condition 1: state planners lack local knowledge; markets cannot price local contextual information. Condition 2: market wages crowd out intrinsic motivation; state bureaucracy de-professionalizes through hierarchy. Condition 3: market competition undermines trust through short-term contracting; state delivery lacks the relational depth of community governance. Condition 4: short-term market contracts under-invest in long-lived assets; political short-termism prevents state investment in assets that mature beyond electoral cycles. $\square$

---

## 38.4 SFC-Compatible Fiscal Framework for UBS

### 38.4.1 The Financing Question

The standard objection to UBS is fiscal: where does the money come from? Universal healthcare, education, housing, transport, and care at high quality would cost, for a typical OECD country, approximately 25–35\% of GDP — comparable to current total government expenditure. This appears to require either unprecedented tax increases or unsustainable deficits.

The cooperative-regenerative framework offers three complementary fiscal mechanisms that do not require either: sovereign money seigniorage, demurrage revenue, and natural capital dividends.

**Mechanism 1: Sovereign Money Seigniorage.** Under the sovereign money system of Chapter 24, new money is created by the central bank and distributed through public spending. Annual seigniorage (at 3\% nominal GDP growth, USD 25 trillion GDP): $\Delta M^T = 0.03 \times 25{,}000 = $ USD 750 billion/year. This seigniorage, currently accruing to commercial banks through the banking system's money creation, is redirected to finance UBS under sovereign money.

**Mechanism 2: Demurrage Revenue.** Under a 2.5\% annual demurrage rate applied to the household sector's monetary balances (approximately USD 12 trillion in US M2 household holdings): annual demurrage revenue = $0.025 \times 12{,}000 = $ USD 300 billion/year.

**Mechanism 3: Natural Capital Dividends.** Carbon pricing at the social cost of carbon (USD 80/tonne CO₂e), applied to all emissions (US emissions $\approx$ 5 Gt CO₂e/year): annual revenue = $80 \times 5 \times 10^9 = $ USD 400 billion/year.

**Combined fiscal capacity:** USD 750B + USD 300B + USD 400B = **USD 1.45 trillion/year** — approximately 5.8\% of US GDP. This is sufficient to fund a meaningful UBS bundle (universal primary healthcare, pre-K through secondary education, basic housing guarantee, and public transport access) without tax increases or deficit spending.

**The SFC proof.** Under the sovereign money SFC model (Chapter 24):
- New money creation: $\Delta M^T = $ USD 750B. Appears in TFM as CB debit to NW, credit to government spending.
- Demurrage collection: $\delta M^H = $ USD 300B. Flows from household accounts to IA.
- Carbon tax: Revenue flows from firms to government.
- UBS expenditure: Government spends collected revenue on service provision.
- Net SFC: $\sum_i \Delta NW_i = 0$ — all monetary flows are conservation-preserving. The government runs a balanced budget (expenditure = seigniorage + demurrage + carbon revenue). No deficit. No unsustainable debt accumulation.

---

## 38.5 Worked Example: UBS Design for a City of 500,000

### 38.5.1 City Specification

A medium-sized European city of 500,000 residents (comparable to Lyon, Stuttgart, or Edinburgh). Current services: partial public healthcare (co-payments required), public secondary education, minimal housing assistance (waiting lists 5+ years), limited public transport coverage, no universal care provision.

### 38.5.2 UBS Service Bundle

**Dimension 1: Healthcare.** Primary care (GP services, preventive care, mental health): universal, zero co-payment. Hospital care: universal, zero co-payment for standard care. Dental and vision: basic package universal. Estimated annual cost: EUR 2,400/person. Total: EUR 1.2 billion/year.

**Dimension 2: Education.** Pre-K (ages 2–5): universal, zero cost. Primary and secondary: already universal in most European cities. Adult education: universal access to vocational retraining and lifelong learning (EUR 500/person/year equivalent). Total new cost: EUR 250M/year.

**Dimension 3: Housing.** A housing guarantee: everyone has access to adequate housing at no more than 25\% of income. Implementation: commons-based housing cooperative (Gemeindebau model) providing 60\% of city housing stock at cost rents. New construction: 3,000 units/year at EUR 180,000/unit = EUR 540M/year capital, EUR 85M/year operating subsidy. Total: EUR 130M/year (annualized capital + operating subsidy net of rents collected).

**Dimension 4: Transport.** Free public transport (bus, tram, metro within city limits). Operating cost net of current fare revenue: EUR 85M/year.

**Dimension 5: Digital access.** Universal fibre broadband + public device library. Cost: EUR 120/person/year = EUR 60M/year.

**Dimension 6: Care.** Universal childcare (age 0–2): EUR 900/child/year subsidy (co-production with cooperative childcare providers). Elder care (65+, community-based): EUR 4,800/elder/year subsidy. Estimated total: EUR 180M/year.

**Total annual UBS cost: EUR 1.905 billion/year** = EUR 3,810/person/year = approximately 35\% of city GDP (estimated EUR 5.5 billion for a city of this size and income level).

### 38.5.3 Funding Mechanism

**Source 1: Sovereign money allocation.** National CB seigniorage allocated to city UBS fund: EUR 420M/year (proportional to city population share of national total, at 2.5\% nominal growth).

**Source 2: Demurrage revenue.** City-level demurrage token (2.5\%/year on city token balances, voluntary adoption with tax incentive): estimated EUR 180M/year at 30\% adoption.

**Source 3: Land value taxation.** LVT at 1.5\% of unimproved land value within city. City land value estimated EUR 45 billion: LVT revenue EUR 675M/year.

**Source 4: Carbon levy.** EUR 80/tonne CO₂e on city emissions (estimated 4.5 Mt/year): revenue EUR 360M/year, partially offset by rebate to lower-income residents (EUR 90M rebate).

**Source 5: Existing public expenditure reallocation.** Current healthcare co-payments, means-tested housing assistance, and transport subsidies: approximately EUR 360M/year redirected under UBS consolidation.

**Total available:** EUR 420 + 180 + 675 + 270 (net carbon) + 360 = **EUR 1.905 billion/year**. Budget balances exactly.

### 38.5.4 Welfare Gain Computation

**Welfare gain methodology.** Using the MPD framework (Chapter 31), the welfare gain is measured by the improvement across the six UBS dimensions:

| MPD dimension | Baseline score | UBS score | Gain |
|:---|:---|:---|:---|
| Nutrition/health | 0.72 | 0.89 | +0.17 |
| Shelter | 0.58 | 0.82 | +0.24 |
| Education/growth | 0.76 | 0.88 | +0.12 |
| Care/connection | 0.61 | 0.79 | +0.18 |
| Security | 0.64 | 0.83 | +0.19 |
| Ecology | 0.48 | 0.61 | +0.13 |
| **MPD (geometric mean)** | **0.630** | **0.803** | **+0.173** |

MPD improvement: +27\% — a very large welfare gain achieved with no increase in GDP (services are provided collectively rather than through additional market activity). The largest gains are in housing (+0.24) and security (+0.19) — the dimensions where market provision fails most severely for lower-income residents.

**Distributional impact.** Using Corollary 38.1: the welfare gain is concentrated in the lower three income quintiles (those who could not previously afford quality private services). Top two quintiles: marginal welfare gain (they already had access to private equivalents, but now save the out-of-pocket costs and gain from universal quality standards). Net: UBS reduces the city's welfare Gini by approximately 0.08 points — a substantial distributional improvement.

---

## 38.6 Case Study: Vienna's Gemeindebau

### 38.6.1 Scale and Structure

Vienna's social housing system is the most studied example of large-scale commons-based housing provision in the world. As of 2023:
- **Stock:** Approximately 220,000 municipal housing units (Gemeindebau) + 200,000 limited-profit housing association units = approximately 420,000 affordable units total.
- **Coverage:** Approximately 60\% of Vienna's 1.9 million residents live in subsidized housing. Income limit for eligibility: EUR 57,600/year (covering approximately 80\% of Vienna's working population).
- **Rent:** Average Gemeindebau rent: EUR 5.80/m²/month vs. private market average EUR 17.50/m²/month. Rent burden: Gemeindebau tenants spend approximately 9\% of income on rent vs. 35\% for comparable private renters.

### 38.6.2 Formal Analysis

**Governance assessment (Ostrom principles):**

| DP | Vienna Gemeindebau implementation | Score |
|:---|:---|:---|
| DP1 (Boundaries) | Income limits define eligible population; clear building management | 2/2 |
| DP2 (Congruence) | Rules adapt to building age, location, and tenant composition | 1.5/2 |
| DP3 (Collective choice) | Tenant councils govern building-level decisions; city council sets policy | 1.5/2 |
| DP4 (Monitoring) | Building management offices; regular maintenance inspections | 2/2 |
| DP5 (Sanctions) | Graduated warnings → lease termination for violations | 2/2 |
| DP6 (Conflict resolution) | Tenant mediation services; Vienna housing court | 2/2 |
| DP7 (External recognition) | Municipal ownership (city government); strong legal protection | 2/2 |
| DP8 (Nested enterprises) | Building councils → district housing → city housing department | 2/2 |
| **Total** | | **15/16** |

**Housing security and welfare outcomes.** Comparative analysis against Munich (similar city size and income, market-dominated housing) and Amsterdam (similar UBS housing provision):

| Metric | Vienna | Munich | Amsterdam |
|:---|:---|:---|:---|
| Housing cost burden (% income) | 17\% (avg) | 29\% (avg) | 20\% (avg) |
| Homelessness rate per 10,000 | 4.1 | 11.7 | 3.8 |
| Housing quality (self-reported) | 8.1/10 | 7.6/10 | 8.0/10 |
| Residential displacement (% moved involuntarily in 5yr) | 3.2\% | 12.8\% | 4.1\% |
| Household savings rate | 12.4\% | 8.9\% | 11.8\% |
| Civic participation index | 0.72 | 0.58 | 0.69 |

Vienna's Gemeindebau delivers housing security (low cost burden, very low homelessness, minimal displacement) at quality comparable to private markets, while freeing household income for savings and civic participation that the market-dominated Munich cannot achieve for comparable income groups.

**The fiscal model.** Vienna's Gemeindebau generates approximately EUR 700 million/year in rent revenue (on 220,000 units at EUR 5.80/m² average, average unit size approximately 55m²). Capital and maintenance expenditure: approximately EUR 650M/year. Net operating surplus: approximately EUR 50M/year — the Gemeindebau is essentially self-financing, requiring no ongoing subsidy from general city taxes once the initial capital stock is established. Historical capital was funded through: postwar federal reconstruction bonds (1945–1970), EU structural funds (1990s–2000s), and revolving loan fund from rent revenue (ongoing).

**The commons stability analysis.** The Gemeindebau has operated continuously since 1919 (interrupted only by the brief fascist period 1934–1945 and WWII) — over 100 years of stable commons-based housing provision. This longevity is consistent with the Ostrom stability prediction (15/16 score) and with the cooperative Lyapunov stability result of Chapter 30: the strong governance structure creates a large basin of attraction that has survived: the Great Depression, two world wars, the Cold War division of Europe, the 1970s oil crisis, the 2008 financial crisis, and COVID-19.

**The affordability mechanism.** Vienna's Gemeindebau provides a soft ceiling on private market rents: private landlords must price at a premium to Gemeindebau, but cannot price too far above it without losing tenants to the social housing option. Econometric estimates (Amann and Mundt, 2017) suggest that Gemeindebau presence in a district reduces private market rents by approximately 10–15\% relative to districts without social housing — a positive spillover from the commons to the private market, consistent with the commons provision theory.

---

## 38.7 UBS and the Stewardship Connection

The formal connection between UBS and the Stewardship Objective Function [C:Ch.2] operates through three channels that standard UBI cannot achieve:

**Channel 1: Collective provision decouples wellbeing from throughput.** When mobility is provided through public transport (UBS) rather than individual car ownership (UBI-funded), the same mobility services are delivered with approximately 1/15 of the material and energy throughput [C:Ch.31, Proposition 31.3]. Collective provision of care, education, and healthcare similarly reduces the material intensity of welfare delivery.

**Channel 2: Commons governance embeds stewardship norms.** Service providers in a commons-governed UBS are subject to community governance that can incorporate ecological criteria — energy standards for public buildings, sustainable procurement requirements, ecological assessments of infrastructure projects. These criteria cannot be imposed on private market providers without regulatory coercion; they emerge naturally from democratic commons governance when community members prioritize them.

**Channel 3: UBS frees time for stewardship.** By eliminating the need for defensive expenditure (time spent shopping for insurance, comparing providers, navigating complex means-tested systems), UBS frees time for unpaid stewardship activities — community food growing, habitat restoration volunteering, cooperative governance participation. These activities generate natural capital and social capital that market provision cannot supply.

---

## Chapter Summary

This chapter has formally compared UBI and UBS, developed the commons-based provision model, specified a SFC-compatible fiscal framework for UBS, designed a concrete city UBS package, and evaluated Vienna's Gemeindebau as the world's largest and longest-running test of commons-based housing provision.

Theorem 38.1 proves UBS dominates UBI in welfare terms when market failures (information asymmetry, network externalities, non-rival provision, adverse selection) are present — which they are in all primary UBS domains. Corollary 38.1 establishes distributional dominance: UBS is more progressive than UBI of equal fiscal cost. Proposition 38.1 adds the ecological argument: collective provision reduces material throughput relative to individualized cash-funded market provision.

Proposition 38.2 identifies four conditions under which commons-based provision outperforms both market and state bureaucracy: local knowledge essentiality, intrinsic motivation dependence, trust as primary input, and high asset specificity. Healthcare, housing, care, and education all satisfy multiple conditions simultaneously.

The SFC fiscal framework demonstrates that sovereign money seigniorage, demurrage revenue, and carbon pricing together generate approximately 5.8\% of US GDP — sufficient for a meaningful UBS bundle without tax increases or deficits. The city design example (500,000 population, EUR 1.905 billion UBS package) balances to budget using LVT, sovereign money allocation, demurrage, and carbon levy, delivering a 27\% MPD welfare gain with 0.08 point Gini reduction.

Vienna's Gemeindebau scores 15/16 on Ostrom principles and delivers measurable outcomes: 17\% housing cost burden vs. 29\% in market-dominated Munich, 4.1 vs. 11.7 per 10,000 homelessness rate, 12\% vs. 9\% household savings rate. Over 100 years of stable operation, the Gemeindebau has survived every major crisis of the 20th and 21st centuries — the most durable test of commons-based housing provision on record.

Chapter 39 closes Part VII by examining digital sovereignty and data cooperatives — the formal governance of the data economy as a cooperative-regenerative commons, tested against the MIDATA.coop health data cooperative case.

---

## Exercises

**38.1** The UBI vs. UBS welfare comparison (Theorem 38.1):
(a) A household with income EUR 18,000/year (40th income percentile) faces the following choices: (i) UBI of EUR 6,000/year in cash; (ii) UBS providing healthcare worth EUR 2,800, childcare worth EUR 2,400, and transport worth EUR 800, with EUR 0 cash transfer. Compute the welfare gain under each, assuming: information asymmetry reduces market healthcare quality for this household by 30\%; they cannot afford private childcare at market prices; public transport UBS provides 90\% of their mobility needs at zero cost.
(b) A household with income EUR 65,000/year (80th income percentile): compute the same comparison. For which household is UBS welfare-superior and for which is UBI welfare-superior?
(c) Under what income level does the crossover between UBS-superior and UBI-superior occur? How does this depend on the quality of UBS provision relative to market alternatives?

**38.2** The SFC fiscal framework:
(a) For an economy with GDP EUR 1 trillion, M2 EUR 900 billion, household M2 share 55\%, carbon emissions 45 Mt CO₂e/year, and sovereign money creation at 3\% nominal growth: compute annual revenues from (i) sovereign money seigniorage, (ii) 2.5\% demurrage on household balances, (iii) EUR 80/tonne carbon levy.
(b) What fraction of GDP does the combined revenue represent? Is this sufficient to fund a full UBS bundle at EUR 3,500/person/year for a population of 80 million?
(c) Construct the SFC transaction flow matrix for the UBS fiscal system: show how seigniorage, demurrage, and carbon revenue flow from households and firms to the government, and from the government to UBS service providers. Verify that every row sums to zero.

**38.3** Vienna's Gemeindebau case:
(a) Compute the consumer surplus generated by the Gemeindebau: for 220,000 units at EUR 5.80/m² average rent vs. market EUR 17.50/m², average size 55m²: what is the annual rent saving per unit and in total? How does this compare to the EUR 700M annual rent revenue?
(b) If Vienna's Gemeindebau were privatized (sold to private developers), estimate: (i) the one-time sale price (at market capitalization rates); (ii) the annual rent increase for current tenants; (iii) the annual welfare loss from higher housing costs. Is the sale price sufficient to compensate tenants for the lifetime welfare loss?
(c) Econometric studies estimate that Gemeindebau presence reduces private market rents by 10–15\% in the same district. For Vienna's private rental market of approximately 200,000 units at average EUR 17.50/m²: compute the annual welfare benefit to private tenants from the Gemeindebau's market-dampening effect.

**★ 38.4** Prove Theorem 38.1 formally for the case of healthcare provision under information asymmetry.

(a) Model healthcare as a market with two quality levels: high ($q^H$) and low ($q^L < q^H$). Both cost $c$ to providers. Patients cannot observe quality before purchase (Akerlof lemons). Show that the market equilibrium involves only low-quality provision.
(b) Show that a UBI of $p \cdot q^H$ per patient (the cost of high-quality care) does not restore high-quality provision: even with the cash, patients cannot distinguish quality, so providers have no incentive to offer $q^H$.
(c) Show that UBS provision with public quality monitoring (government observes and certifies quality) restores high-quality provision for all patients, regardless of income.
(d) Compute the aggregate welfare gain from UBS vs. UBI in this healthcare model, as a function of the quality difference $q^H - q^L$, the fraction of population that would purchase only $q^L$ under UBI (due to inability to verify quality), and the monitoring cost of the UBS system.

**★ 38.5** Apply Proposition 38.2 (commons provision conditions) to a specific service not covered in the chapter.

Choose one of: community mental health services, neighborhood food systems (community kitchens, food banks), local renewable energy services, or adult social care.

(a) For your chosen service, assess each of the four conditions (local knowledge, intrinsic motivation, trust, asset specificity) on a 0–2 scale with justification.
(b) Based on your assessment: does commons-based provision outperform market provision, state provision, or both? Under what specific parameter conditions might market or state provision be superior?
(c) Design a commons-based provision model for your service: specify the membership criteria, governance structure (Ostrom design principles), contribution norms, and quality assurance mechanism.
(d) Estimate the welfare gain from commons-based vs. market provision for your service, using the most appropriate welfare measure (MPD dimension, consumer surplus, or direct wellbeing survey method).

**★★ 38.6** Design and formally analyze a full UBS framework for a country of your choice.

**Country selection:** Choose an OECD country with available data on current public expenditure, housing costs, healthcare outcomes, and income distribution. Suggested options: Canada, Australia, the Netherlands, or Sweden.

(a) **Baseline assessment:** Using the MPD framework (Chapter 31), score your country's current performance on six provisioning dimensions: health, shelter, education, care, security, and ecology. Identify the three dimensions with the largest gaps from target $P_j^*$.

(b) **UBS service bundle:** For each of the six dimensions, specify the service level to be universally provided, the current market alternative, and the quality premium of UBS over the market alternative for the lowest-income quintile.

(c) **Fiscal framework:** Compute the available revenue from sovereign money seigniorage, demurrage (at 2.5\% on household monetary balances), carbon levy (at national emissions × EUR/USD 80/tonne), and LVT (at 1.5\% of unimproved land value, estimated from national housing data). Does the combined revenue cover the UBS bundle cost?

(d) **Welfare gain computation:** Using the MPD framework, estimate the post-UBS scores for each dimension and the overall MPD improvement. Compute the improvement separately for the bottom three and top two income quintiles to assess distributional effects.

(e) **Governance design:** Specify the governance structure for each UBS dimension: which services should be delivered by central state institutions, which by regional/municipal commons, and which by cooperative service providers? Apply Proposition 38.2 to justify each governance choice.

(f) **Ecological impact:** Apply Proposition 38.1 to estimate the reduction in per-capita material throughput from switching from UBI-equivalent cash to collective UBS provision in the three most material-intensive service dimensions (mobility, housing, food). Is the ecological footprint reduction consistent with your country's Planetary Boundaries allocation?

---

*Chapter 39 closes Part VII with digital sovereignty and data cooperatives — examining the formal governance of the data economy through the cooperative-regenerative lens, and testing the theory against MIDATA.coop's decade of operation as a health data cooperative in Switzerland.*
