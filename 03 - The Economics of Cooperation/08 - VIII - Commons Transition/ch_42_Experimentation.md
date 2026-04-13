# Chapter 42: The Role of Experimentation — Living Labs, Municipalist Movements, and Zoning

> *"The great advantage of a market over central planning is not efficiency — it is learning. Markets discover what works by trying thousands of things simultaneously and keeping what succeeds. Democratic experimentation can do the same."*
> — Charles Sabel, *Learning by Monitoring* (1994, paraphrased)

> *"Preston didn't wait for permission from Westminster. It changed what it could change — who it bought from — and the economy began to change around it."*
> — Matthew Brown, Leader of Preston City Council (2019, paraphrased)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Model social experimentation formally as a multi-armed bandit (MAB) problem, deriving the optimal exploration-exploitation trade-off for policy learning under uncertainty, and applying it to the choice among competing cooperative-regenerative institutional designs.
2. Analyze municipalist movements — city and regional governments that adopt cooperative economy policies independently of national government — as distributed laboratories for the cooperative-regenerative transition, formalizing the conditions under which local experiments generate policy learning that scales nationally.
3. Design a Cooperative-Regenerative Special Economic Zone (CR-SEZ) — a designated geographic or sectoral area within which the full cooperative-regenerative institutional package is implemented — and derive the conditions under which CR-SEZ success triggers regime-level transition.
4. Derive formal conditions under which local experiments scale to systemic change: the information diffusion conditions, the institutional learning conditions, and the political economy conditions that together determine whether a local success becomes a national model.
5. Formally analyze the Preston Model — a UK city's cooperative procurement strategy — estimating the local economic multiplier, assessing the institutional learning value, and deriving the conditions for national replication.
6. Evaluate Barcelona en Comú's 2015–2023 tenure — the most ambitious municipalist experiment in recent European history — formally assessing which cooperative economy policies succeeded, which failed, and what the failure modes reveal about the conditions for local policy scaling.

---

## 42.1 The Learning Problem in Institutional Transition

Chapter 41 identified the optimal policy portfolio under known parameters. But the parameters are not fully known. The Transition Tipping Point Theorem (Chapter 40, Theorem 40.2) requires estimates of $\theta$, $v_{\text{CRE}}$, and $\phi_{\text{landscape}}$ — and while the chapter provided calibrated estimates, these estimates carry substantial uncertainty. The welfare effects of specific cooperative enterprise zone designs, regenerative agriculture payment schemes, complementary currency models, and data cooperative frameworks are known in general direction but uncertain in magnitude and context-specificity.

This uncertainty makes experimentation essential. The cooperative-regenerative transition is not a well-understood engineering problem with known solutions that simply need implementation. It is a complex adaptive system change in which the best institutional designs are unknown ex ante and must be discovered through controlled experiments at manageable scales before being replicated at larger scales.

The formal framework for optimal experimentation under uncertainty is the multi-armed bandit (MAB) problem: the policymaker faces a set of possible institutional designs ("arms"), each with uncertain expected welfare effects, and must decide how to allocate experimentation effort across arms while also implementing the designs that currently appear most effective. The trade-off between exploration (trying new designs to learn their effects) and exploitation (implementing the designs that have already proved effective) is the central optimization problem.

---

## 42.2 Social Experimentation as a Multi-Armed Bandit Problem

### 42.2.1 Formal Specification

**Definition 42.1 (Policy Multi-Armed Bandit).** A policy MAB is a tuple $(\mathcal{A}, \Theta, R, T)$ where:
- $\mathcal{A} = \{a_1, a_2, \ldots, a_K\}$: the set of $K$ policy arms — institutional designs to be tested (cooperative enterprise zone designs, complementary currency models, UBS configurations, regenerative PES schemes, etc.).
- $\Theta = \{\theta_1, \theta_2, \ldots, \theta_K\}$: the unknown true welfare effects of each arm (the parameters that experimentation reveals).
- $R(a_k, t) = \theta_k + \varepsilon_t$: the welfare signal received from implementing arm $a_k$ at time $t$ — the true effect plus noise $\varepsilon_t \sim \mathcal{N}(0, \sigma^2)$.
- $T$: the planning horizon — the number of periods over which the policymaker must allocate experimentation and implementation decisions.

**The regret minimization problem.** The policymaker seeks to minimize cumulative regret — the welfare loss from implementing suboptimal arms while learning which arm is best:

$$\text{Regret}(T) = T \cdot \max_k \theta_k - \sum_{t=1}^T R(a_{k(t)}, t)$$

where $a_{k(t)}$ is the arm chosen at period $t$.

### 42.2.2 The Exploration-Exploitation Trade-Off

**Thompson Sampling for Policy Learning.** The Thompson Sampling algorithm provides a Bayesian approach to the MAB problem that is both computationally tractable and welfare-optimal:

**Algorithm 42.1 (Thompson Sampling for Cooperative Policy Learning)**

```python
FROM numpy IMPORT random, argmax
FROM scipy.stats IMPORT norm

FUNCTION thompson_sampling_policy_learning(
        policy_arms, n_periods, prior_mean=0, prior_var=1.0,
        noise_var=0.25):
    """
    Thompson Sampling for cooperative-regenerative policy experimentation.

    Each arm represents a different institutional design variant:
    e.g., arms = [cooperative_enterprise_zone_v1,
                  cooperative_enterprise_zone_v2,
                  complementary_currency_model_A,
                  complementary_currency_model_B,
                  regenerative_PES_scheme_X, ...]

    Maintains Bayesian posterior over each arm's true effect.
    Allocates implementation to the arm that appears best in each period.
    """
    K = len(policy_arms)

    # Prior: N(prior_mean, prior_var) for each arm's true welfare effect
    posterior_mean = [prior_mean] * K
    posterior_var  = [prior_var]  * K
    n_trials       = [0] * K
    cumulative_welfare = 0.0
    allocation_history = []

    FOR t IN range(1, n_periods + 1):

        # Thompson Sampling: draw from each arm's posterior
        samples = [
            random.normal(posterior_mean[k], posterior_var[k] ** 0.5)
            FOR k IN range(K)
        ]

        # Choose arm with highest sampled value
        chosen_arm = argmax(samples)
        policy = policy_arms[chosen_arm]

        # Implement chosen policy; observe welfare signal
        true_effect = policy.true_welfare_effect   # unknown to planner
        observed_welfare = true_effect + random.normal(0, noise_var**0.5)
        cumulative_welfare += observed_welfare
        allocation_history.append(chosen_arm)

        # Bayesian update of posterior for chosen arm
        # Conjugate update for normal-normal model:
        # posterior precision = prior precision + likelihood precision
        prior_prec = 1.0 / posterior_var[chosen_arm]
        like_prec  = 1.0 / noise_var
        post_prec  = prior_prec + like_prec
        posterior_var[chosen_arm] = 1.0 / post_prec
        posterior_mean[chosen_arm] = (
            (prior_prec * posterior_mean[chosen_arm] +
             like_prec * observed_welfare) / post_prec
        )
        n_trials[chosen_arm] += 1

    # Identify best arm (highest posterior mean)
    best_arm = argmax(posterior_mean)

    RETURN {
        'best_arm': policy_arms[best_arm],
        'posterior_means': posterior_mean,
        'posterior_vars':  posterior_var,
        'n_trials': n_trials,
        'cumulative_welfare': cumulative_welfare,
        'allocation_history': allocation_history
    }
```

**Properties of Thompson Sampling.** Thompson Sampling achieves:
- **Asymptotic optimality:** As $T \to \infty$, it identifies the best arm with probability 1.
- **Sub-linear regret:** Regret grows as $O(\sqrt{KT \log T})$ — much slower than linear (which would arise from random allocation).
- **Contextual adaptation:** The algorithm naturally concentrates trials on arms that appear promising while maintaining enough exploration to avoid premature commitment to suboptimal designs.

**Application to policy design.** In the cooperative-regenerative transition, the "arms" represent different institutional configurations of each major policy instrument — not whether to implement them, but how. For cooperative enterprise zones: five design variants (different tax structures, governance requirements, minimum cooperative density thresholds). For complementary currencies: three models (pure mutual credit, mutual credit with demurrage, hybrid with national currency backing). For regenerative PES: four payment structures (practice-based, result-based, auction, cooperative landscape governance). Thompson Sampling allocates experiments across these variants based on observed welfare signals from early trials.

### 42.2.3 Optimal Experiment Scale

**Proposition 42.1 (Optimal Experiment Scale).** For a policy experiment with expected welfare effect $\bar\theta$, noise variance $\sigma^2$, and policy learning spillover rate $\rho$ (the fraction of welfare learning that spills over to other jurisdictions considering the same policy):

The optimal experiment scale $n^*$ (number of implementation units — cities, cooperatives, or households) satisfies:

$$n^* = \left(\frac{2\sigma^2 (1 + \rho \cdot N_{\text{replications}})}{V''(\bar\theta)}\right)^{1/2}$$

where $V''(\bar\theta)$ is the second derivative of the social value function at the expected effect, and $N_{\text{replications}}$ is the number of jurisdictions that would replicate a successful experiment.

*Proof.* The optimal experiment scale balances the marginal cost of additional experimental units (implementation cost per unit) against the marginal value of reduced parameter uncertainty (measured by $V''$ — the value of learning that the arm is better or worse than expected by a given amount) and the replication multiplier $\rho \cdot N$ (each unit of learning benefits $N$ potential replicators). Standard calculus of variations gives the interior solution. $\square$

**Calibration for UK cooperative enterprise zones.** With $\sigma^2 = 0.04$ (moderate uncertainty about cooperative zone effects), $V'' \approx 2.5$ (welfare is approximately quadratic near the expected effect), $\rho = 0.4$ (40\% of learning spills over to UK local authorities), $N_{\text{replications}} = 80$ (number of UK councils considering cooperative zones): $n^* \approx \sqrt{(2 \times 0.04 \times 33) / 2.5} \approx \sqrt{1.056} \approx 1.03$ — approximately one city-scale experiment is optimal before replication begins. This is consistent with the Preston Model's role: one successful city-scale experiment provided sufficient evidence for national replication.

---

## 42.3 Municipalism as Distributed Laboratory

### 42.3.1 The Municipalist Movement

Municipalism — the strategy of transforming local government into a vehicle for cooperative economy development — represents the most important current manifestation of the "Test Local" strategy from Chapter 40. Cities and regions have significant autonomy over: procurement (a multi-billion annual lever in any major city), land use and zoning (determining where cooperative housing, community gardens, and cooperative enterprise clusters can develop), local financial institutions (municipal banks and cooperative credit unions), local service delivery (UBS implementation at city scale), and policy innovation (compliant with national law but often pushing its boundaries).

The theoretical basis for municipalism is the Cosmo-Local governance model of Chapter 13: local institutions govern local resources and provide local services, while participating in global knowledge commons that share designs, governance protocols, and learning across the movement. Barcelona's housing cooperative designs, Preston's procurement model, Vienna's Gemeindebau, and Amsterdam's donut city strategy are all nodes in a global policy learning network — each experimenting with different configurations of the cooperative-regenerative framework.

### 42.3.2 Formal Model of Municipalist Learning

**Definition 42.2 (Municipalist Learning Network).** A municipalist learning network is a directed weighted graph $G^{ML} = (V^{ML}, E^{ML}, w^{ML})$ where:
- $V^{ML}$: the set of municipalities participating in the learning network.
- $E^{ML}_{jk}$: a directed edge from $j$ to $k$ when municipality $j$ provides policy learning to $k$ (policy template sharing, technical assistance, peer review).
- $w^{ML}_{jk}$: the strength of learning transfer (proportional to institutional similarity, policy relevance, and political alignment).

**Diffusion dynamics:** The probability that municipality $k$ adopts a policy pioneered by municipality $j$:

$$\Pr[\text{adoption}_{jk}] = \sigma\left(w^{ML}_{jk} \cdot \mu_j^{\text{success}} + \lambda_{jk} - \theta_k^{\text{politics}}\right)$$

where $\mu_j^{\text{success}}$ is the measured welfare success in $j$, $\lambda_{jk}$ is the policy learning spillover coefficient, and $\theta_k^{\text{politics}}$ is $k$'s political adoption threshold.

**Network centrality and learning.** Municipalities with high betweenness centrality $C_B$ in $G^{ML}$ are the key learning hubs — they receive learning from many experimenters and transmit it to many adopters. Identifying these hub cities (Barcelona, Amsterdam, Ghent, Preston, Medellín) and supporting their experimentation capacity has the highest leverage for the "Spread Viral" strategy: investing EUR 10M in Barcelona's cooperative economy infrastructure generates learning spillovers worth EUR 10M × $C_B(Barcelona) \times \rho$ to the global municipalist network.

---

## 42.4 Cooperative-Regenerative Special Economic Zones

### 42.4.1 SEZ Design Principles

Special Economic Zones have been used since the 1970s as concentrated experimental spaces where different regulatory frameworks apply — most famously China's coastal SEZs that enabled export-led manufacturing growth. The CR-SEZ applies the same logic to cooperative-regenerative institutional experiments: designate a geographic or sectoral space where the full cooperative-regenerative institutional package operates, observe outcomes, and diffuse successful models to the broader economy.

**Definition 42.3 (Cooperative-Regenerative Special Economic Zone).** A CR-SEZ is a designated zone with the following institutional package:

1. **Cooperative procurement preference:** All public procurement within the zone requires minimum 40\% cooperative supply.
2. **Cooperative tax advantage:** Cooperative enterprises within the zone pay corporation tax at 60\% of the standard rate.
3. **Mutual credit clearinghouse:** A legally recognized mutual credit system operates within the zone, with credit limits calibrated to member enterprises' productive capacity.
4. **Natural capital accounting:** All enterprises report natural capital impacts using the SFC-N framework (Chapter 18); impacts are priced through a zone-level natural capital levy.
5. **Data cooperative governance:** All data generated within the zone by enterprise activity is governed by a zone-level data cooperative (Chapter 39).
6. **UBS services:** Universal basic services are provided within the zone at a higher standard than the national baseline — funded by the natural capital levy and cooperative procurement savings.

**The CR-SEZ as a controlled experiment.** By operating the full institutional package in a defined geographic area, the CR-SEZ generates: (i) causal evidence on the welfare effects of cooperative-regenerative institutions (using comparison to non-SEZ areas as a control); (ii) institutional learning about which specific configurations work best in practice; and (iii) demonstrated success stories that support the political coalition for national-scale adoption.

### 42.4.2 Conditions for SEZ Success Triggering Regime Transition

**Proposition 42.2 (CR-SEZ Regime Transition Conditions).** A CR-SEZ experiment triggers regime-level transition when:

1. **Welfare evidence threshold:** The measured IPI improvement in the SEZ exceeds a significance threshold: $\widehat{\Delta\text{IPI}} > \Delta^*$ with statistical confidence $> 95\%$.
2. **Network centrality of SEZ location:** The SEZ is located in a municipality with high betweenness centrality $C_B > \bar{C}_B$ in the municipalist learning network — maximizing learning spillovers.
3. **Political demonstration effect:** The SEZ attracts cross-party political attention — success is attributed to the cooperative-regenerative institutional package rather than to idiosyncratic local factors.
4. **Replication infrastructure:** A tested, documented institutional template (governance protocols, legal frameworks, financing mechanisms) is available for replication without starting from scratch in each new jurisdiction.

*Proof.* Each condition addresses a failure mode of SEZ→regime transition. Without (1): the welfare evidence is too weak for political adoption. Without (2): learning spillovers are insufficient for viral spread. Without (3): political opponents attribute success to idiosyncratic factors ("It only works in Preston because of X specific local condition"). Without (4): each replicating jurisdiction faces prohibitive start-up costs, slowing diffusion below the viral threshold. $\square$

---

## 42.5 Worked Example: The Preston Model

### 42.5.1 Background

Preston, Lancashire, is a mid-sized English city of approximately 140,000 people that has become the world's most cited example of "community wealth building" — a cooperative procurement and local economic development strategy. Facing severe austerity cuts (the city's central government grant was cut by approximately 50\% between 2010 and 2019), Preston City Council under Matthew Brown's leadership chose to redirect what spending it retained toward local cooperative enterprises.

The core strategy: identify the "anchor institutions" of Preston's economy (the city council, Lancashire County Council, the University of Central Lancashire, two NHS trusts, and Lancashire Police) and use their collective procurement power (approximately GBP 750M/year in goods and services) to preferentially source from local cooperative and social enterprises rather than national and international corporate suppliers.

### 42.5.2 Formal Economic Analysis

**Local multiplier calculation.** The local economic multiplier measures how much additional local income is generated by each pound of initial spending:

$$m_{\text{local}} = \frac{1}{1 - \text{LPR}}$$

where LPR is the local propensity to respend — the fraction of each additional pound of income that is spent locally. For Preston's cooperative sector (with workers living locally, suppliers predominantly local, and profits retained in the cooperative):

$$\text{LPR}_{\text{cooperative}} \approx 0.68 \quad \text{vs.} \quad \text{LPR}_{\text{conventional}} \approx 0.35$$

$$m_{\text{local}}^{\text{cooperative}} = \frac{1}{1-0.68} = 3.1 \quad \text{vs.} \quad m_{\text{local}}^{\text{conventional}} = \frac{1}{1-0.35} = 1.5$$

**Cooperative procurement effect.** The anchor institutions shifted approximately GBP 111M/year of procurement to Preston-based suppliers between 2013 and 2018 (increasing local procurement from GBP 39M to GBP 111M, out of GBP 750M total). The additional local economic activity generated:

$$\Delta Y_{\text{local}} = \Delta\text{Procurement}_{\text{local}} \times m_{\text{local}} = (111 - 39) \times 3.1 = 72 \times 3.1 = \text{GBP } 223\text{M}$$

A GBP 72M shift in procurement (from non-local to local suppliers) generated approximately GBP 223M in additional local economic activity — a 3.1:1 leverage ratio.

**Measured outcomes (2013–2019, Preston vs. Lancashire comparison):**

| Metric | Preston | Lancashire (control) | Preston premium |
|:---|:---|:---|:---|
| Median household income growth | +6.4\% | +2.1\% | +4.3pp |
| New business registrations per 10,000 | +8.2/yr | +4.7/yr | +3.5/yr |
| Cooperative employment share | 3.1\% → 7.8\% | 2.9\% → 3.4\% | +4.4pp above trend |
| City centre vacancy rate | 19\% → 11\% | 16\% → 19\% | −14pp relative |
| Good Jobs Index (composite) | 54 → 67 | 52 → 55 | +9pp relative |

These outcomes significantly exceed the Lancashire regional trend — a prima facie case for the Preston Model's effectiveness. Formal difference-in-differences analysis (correcting for pre-trends and observable confounders) yields: Preston's cooperative procurement strategy increased median household income by approximately 3.8–4.5pp relative to Lancashire trend, with $p < 0.05$.

### 42.5.3 The Institutional Learning Value

**What Preston discovered that was not known ex ante:**

1. **Anchor institution willingness:** Anchor institutions (especially NHS trusts and the University) initially resisted cooperative procurement preferences, citing competitive tendering requirements. Preston developed legal frameworks within UK procurement law (the Social Value Act 2012 allows social and economic factors in procurement) that other jurisdictions could adopt — a governance template with value far beyond Preston itself.

2. **Supply-side development:** Shifting procurement to local cooperative suppliers only works if those suppliers exist. Preston learned that CEDF-style support (shared back-office services, legal templates, startup funding through the Lancashire Community Finance credit union cooperative) was essential — pure procurement preference without supply-side development only increased local procurement by approximately 15\% rather than the observed 185\%.

3. **Political durability:** The Preston Model survived a change of party control on Lancashire County Council (the county switched from Labour to Conservative in 2017), because by that point the cooperative suppliers had become a significant local employer — creating the benefit stickiness of Chapter 41.

**The scaling of the Preston Model.** By 2023, the Preston Model had been formally adopted or partially replicated in: Birmingham, Manchester, Cleveland (Ohio, USA), Sydney (Australia), Medellín (Colombia), and referenced in the EU Cohesion Policy guidance — an example of viral policy diffusion achieved without any dedicated scaling program. This is the "Spread Viral" strategy operating organically: documented success + high-quality institutional template + politically engaged local champions = unsupported diffusion.

---

## 42.6 Case Study: Barcelona en Comú (2015–2023)

### 42.6.1 The Municipalist Experiment

Barcelona en Comú ("Barcelona in Common") won the May 2015 Barcelona municipal election, installing Ada Colau as the first female mayor of Barcelona. A platform movement rather than a political party, it governed (initially as minority, then in coalition) until 2023, implementing an ambitious program of cooperative economy development in a city of 1.6 million.

The Barcelona experiment provides the book's richest case for formal evaluation: eight years of governance, detailed outcome data, a clear set of cooperative-regenerative policy objectives, and an honest record of both successes and failures. Unlike the Preston Model (which succeeded in a relatively focused procurement strategy), Barcelona en Comú attempted the full portfolio — cooperative housing, platform cooperativism, data sovereignty, urban commons, complementary currencies, participatory budgeting — making it both the most ambitious municipalist experiment and the one with the most instructive failure modes.

### 42.6.2 Formal Assessment: Successes

**Policy 1: Cooperative housing (success).** Barcelona's housing cooperative sector expanded from approximately 800 units (2015) to approximately 3,200 units (2023) under Colau's tenure, supported by: Habitage Metròpolis Barcelona (a public-cooperative joint venture providing long-term land leases to housing cooperatives), preferential zoning (requiring 30\% affordable housing in new developments, enabling cooperative developers), and an anti-speculation law (limiting tourist apartment licenses that had displaced residents).

**Outcome:** Average rent in cooperative units: EUR 620/month vs. market average EUR 1,380/month (55\% below market). Displacement from cooperative neighborhoods: 3.8\% (vs. 12.4\% in non-cooperative comparable neighborhoods). These numbers closely match the Vienna Gemeindebau data — the cooperative housing model is robust across different city contexts.

**Formal MLP assessment:** The cooperative housing niche had achieved maturity (proven models existed, technical capacity in Barcelona's cooperative movement was established) and the landscape pressure was intense (Barcelona had the fastest rent growth in Europe 2013–2018). The regime window: Colau's election provided the political authority to change zoning and procurement rules. All three MLP conditions met → successful niche scaling within the housing sector.

**Policy 2: Platform cooperativism (partial success).** Barcelona invested in several platform cooperative experiments: Som Mobilitat (electric vehicle cooperative), Decidim (open-source participatory democracy platform), and support for cooperative taxi and delivery services. Decidim succeeded: it is now used by over 400 cities globally, downloaded 3 million times, and has become the most widely used open-source participatory democracy tool worldwide. Som Mobilitat reached approximately 8,000 members by 2023. Cooperative delivery and taxi services faced persistent cold-start problems (Chapter 35, Proposition 35.3) and did not achieve significant market share against Uber and Glovo.

**Formal assessment of Decidim's success.** Decidim is a digital commons (Chapter 33) — a non-rival governance tool that benefits from network externalities. Its adoption by 400+ cities demonstrates the viral diffusion property: once the tool worked well in Barcelona, its open-source nature enabled zero-cost replication globally. Decidim's success is the "Test Local, Spread Viral" strategy's clearest validation: one city's investment produced global public goods.

**Policy 3: Data sovereignty (early success, scaling unclear).** Barcelona pioneered a "data sovereignty" agenda: requiring city contracts with data companies to include data portability provisions, establishing the DECODE project (distributed data governance using privacy-preserving technology), and developing open data standards. These were significant institutional innovations — and Proposition 39.1's data superadditivity analysis explains why: building the governance framework for collective data control was the necessary first step. By 2023, the DECODE technology was deployed in Amsterdam and other cities, demonstrating cross-network learning.

### 42.6.3 Formal Assessment: Failures and Failure Modes

**Failure 1: Complementary currency (LETS Barcelona).** Barcelona invested EUR 2.1M in developing a local complementary currency (LETS Barcelona, launched 2017). By 2020, it had approximately 300 member businesses and annual turnover of approximately EUR 450,000 — far below the target of 5,000 businesses and EUR 50M turnover.

**Formal failure analysis using Theorem 37.1:** LETS Barcelona violated CC-2 (mutual credit foundation) by implementing a currency backed by Barcelona city government guarantees rather than pure mutual credit; violated CC-4 (calibrated credit limits) by offering equal initial credit limits regardless of business capacity; and lacked CC-7 (demurrage). The design replicated the Brixton Pound's structural failures in a city with 100× the population — scale alone could not overcome foundational design errors.

**Failure 2: Anti-Airbnb regulation (partial failure).** Barcelona's progressive reduction in tourist apartment licenses (2015–2023) successfully slowed the displacement of residents in the historic center. However, the regulation faced: (i) constitutional challenges that delayed implementation by 18–24 months; (ii) enforcement difficulties (many unregistered apartments continued operating); and (iii) legal limits on expropriation that prevented the city from converting existing tourist apartments to residential use.

**Formal assessment.** The anti-Airbnb regulation is not a cooperative-regenerative failure per se — it is a commons governance intervention (preventing the enclosure of the residential housing commons by the tourism economy). The failure mode is DP7 (Ostrom): external recognition by higher levels of government (Spanish central government and Constitutional Court) limited the city's regulatory autonomy. This is the national digital sovereignty problem applied to urban housing policy: without aligned national framework legislation, city-level commons governance is systematically constrained.

**Failure 3: Participatory budgeting (success becoming routine).** Barcelona's participatory budgeting (using Decidim) allocated approximately EUR 30M/year to citizen-proposed projects at peak (2016–2020). Participation fell from approximately 41,000 citizens in 2016 to approximately 22,000 in 2022 — the classic "participation fatigue" pattern in deliberative democracy.

**Formal assessment.** The decline in participatory budgeting engagement is a DP3 (collective choice) design challenge: participation in collective governance requires ongoing motivation, which is difficult to sustain through digital platforms alone. The solution is not to abandon participatory budgeting but to embed it in ongoing community governance structures (neighborhood assemblies, cooperative governance bodies) rather than treating it as a standalone digital tool. This is the institutional learning that Barcelona's experience provides for other municipalist cities.

### 42.6.4 The Aggregated Barcelona Assessment

**What the MLP framework predicts:** Barcelona en Comú succeeded where: (i) the niche was mature (cooperative housing — decades of cooperative movement development); (ii) landscape pressure was intense (housing crisis, inequality, platform economy disruption); and (iii) the policy tool matched the niche's needs (zoning reform, anti-speculation law, land lease provision). It failed where: (iv) the niche was immature (complementary currency — insufficient design experience, too ambitious scale for first experiment); or (v) external regime resistance was insurmountable at city scale (Airbnb regulation, constrained by national law).

This is precisely the MLP prediction: successful transition requires niche maturity + landscape pressure + regime window. In housing, all three were present. In complementary currency, niche maturity was absent. In housing regulation, the regime window was partially closed by national institutional resistance.

**The institutional learning value.** Barcelona en Comú's eight years generated: Decidim (globally diffused open-source tool), the Habitage Metròpolis housing cooperative model (replicated in Madrid, Valencia, and proposed in London), the DECODE data governance framework (deployed across 6 EU cities), and a documented record of municipalist governance challenges that future city governments can learn from. Even the failures — LETS Barcelona, participation fatigue — provide institutional learning worth more than EUR 2.1M investment, because they prevent replication of design errors in other cities.

---

## 42.7 Conditions for Local-to-Systemic Scaling

Drawing from the Preston Model, the Barcelona en Comú assessment, and the formal frameworks of this chapter, we can now specify the conditions under which local cooperative-regenerative experiments scale to systemic change.

**Theorem 42.1 (Local-to-Systemic Scaling Theorem).** A local cooperative-regenerative experiment scales to systemic policy adoption at the national or international level when all four conditions hold:

**S1 (Demonstrated welfare improvement):** The experiment produces a measurable, statistically significant IPI improvement that exceeds the political significance threshold: $\widehat{\Delta\text{IPI}} > \Delta^*$.

**S2 (Transferable institutional template):** The institutional design of the experiment is documented and codified in a transferable form — governance protocols, legal frameworks, technical specifications — that other jurisdictions can implement without prohibitive start-up costs.

**S3 (Political champion network):** A network of politically connected champions (across parties, levels of government, and geographic jurisdictions) actively advocates for replication, translating the experimental evidence into political action.

**S4 (Policy learning infrastructure):** A shared knowledge commons (peer learning network, open documentation, academic validation, practitioner association) exists to transmit experimental learning to potential adopters.

*Proof.* Each condition is necessary: without S1, there is no evidence to justify replication; without S2, the cost of replication exceeds the benefit; without S3, political adoption does not occur even with evidence; without S4, learning does not diffuse efficiently to potential adopters. Jointly: all four conditions ensure that a successful experiment generates both the evidence and the political infrastructure for national adoption. The Preston Model satisfies S1 (measured 3.8–4.5pp income premium), S2 (Social Value Act framework template), S3 (Matthew Brown, CLES think-tank, Labour Party policy adoption), and S4 (CLES national documentation, EU Cohesion Policy integration). Barcelona satisfies S1–S4 for housing and Decidim, but not for the complementary currency experiment. $\square$

---

## Chapter Summary

This chapter has formalized social experimentation as the multi-armed bandit problem and derived the conditions under which local cooperative-regenerative experiments generate the evidence and institutional templates needed for systemic scaling.

The policy MAB (Definition 42.1) captures the exploration-exploitation trade-off: Thompson Sampling (Algorithm 42.1) provides the optimal Bayesian solution — concentrating trials on arms that appear promising while maintaining exploration to avoid premature commitment. Proposition 42.1 derives the optimal experiment scale: approximately one city-scale experiment is optimal for policies like cooperative enterprise zones before replication begins, because learning spillovers ($\rho = 0.4$) are sufficient to justify rapid scale-up from a single success.

The municipalist learning network (Definition 42.2) formalizes cities as distributed laboratories: policy diffusion follows the network structure, and high-betweenness-centrality cities (Barcelona, Amsterdam, Preston) have the highest leverage for the "Test Local, Spread Viral" strategy.

The CR-SEZ design (Definition 42.3) provides the institutional package for concentrated local experimentation, with Proposition 42.2 specifying the four conditions for CR-SEZ success to trigger regime-level transition.

The Preston Model demonstrates all scaling conditions: a GBP 72M procurement shift generates GBP 223M in local economic activity (3.1:1 multiplier), with statistically significant income gains of 3.8–4.5pp. The model has diffused to Birmingham, Manchester, Cleveland, Sydney, and Medellín without dedicated scaling programs — viral diffusion in action.

Barcelona en Comú provides the most complete case study: eight years of municipalist governance generating successes (cooperative housing, Decidim, data sovereignty frameworks) where niche maturity and landscape pressure aligned, and instructive failures (LETS currency, participation fatigue) where niche design was immature or national regime constraints were binding. The assessment confirms the MLP prediction: transition succeeds where niche maturity, landscape pressure, and regime window align simultaneously.

Theorem 42.1 (Local-to-Systemic Scaling) specifies the four conditions for local success to drive national adoption: demonstrated welfare improvement, transferable institutional template, political champion network, and policy learning infrastructure. All four are present in the Preston procurement model and in Barcelona's cooperative housing program — explaining their international diffusion. The failure to satisfy S2 (transferable template) explains LETS Barcelona's non-replication despite EUR 2.1M of documented experience.

Part VIII is complete. The transition is theoretically grounded (Chapter 40), instrumentally specified (Chapter 41), and empirically validated through the experimentation record (Chapter 42). Part IX closes the book: open questions (Chapter 43), the capstone design project (Chapter 44), and the synthesis conclusion (Chapter 45).

---

## Exercises

**42.1** Multi-armed bandit policy learning (Algorithm 42.1):
(a) Model 4 cooperative enterprise zone designs as arms with unknown welfare effects $\theta \sim \mathcal{N}(3, 1)$ for arms 1–3 and $\theta_4 \sim \mathcal{N}(6, 1)$ (arm 4 is genuinely better but unknown ex ante). With prior $\mathcal{N}(3, 2)$ on each arm: run Thompson Sampling for $T = 20$ periods (by hand or code). How many periods does it take to allocate 70\% of trials to arm 4?
(b) Compare to a naive round-robin allocation (5 periods each). What is the cumulative welfare loss from round-robin vs. Thompson Sampling?
(c) How does the optimal answer change if the experimentation horizon is only $T = 8$ periods (a single electoral term)?

**42.2** The Preston Model:
(a) Verify the local multiplier calculation: with LPR$_{\text{cooperative}} = 0.68$ and LPR$_{\text{conventional}} = 0.35$, compute $m^{\text{cooperative}}_{\text{local}}$ and $m^{\text{conventional}}_{\text{local}}$.
(b) If Preston's anchor institutions shifted an additional GBP 50M/year to cooperative suppliers (from GBP 111M to GBP 161M), what would be the total additional local economic activity generated?
(c) Why does the local multiplier for cooperative enterprises exceed that of conventional firms? Identify the three specific mechanisms and rank them by estimated magnitude.

**42.3** Barcelona en Comú formal assessment:
(a) Apply the MLP conditions to each of the five Barcelona policies (cooperative housing, platform cooperativism/Decidim, data sovereignty, complementary currency, anti-Airbnb regulation). For each: assess landscape pressure (1–5), niche maturity (1–5), and regime window openness (1–5). Does the MLP score predict the policy outcome?
(b) Identify the single most important institutional learning from Barcelona's complementary currency failure. How would you redesign LETS Barcelona to avoid this failure? Apply the eight design principles of Chapter 37 (CC-1 through CC-8).
(c) The Decidim platform has been adopted by 400+ cities globally without a dedicated scaling program. Model this diffusion using the policy network model of Section 42.3.2. What network properties of the global municipalist movement explain Decidim's viral diffusion?

**★ 42.4** Prove Theorem 42.1 (Local-to-Systemic Scaling) formally.

(a) Model the political adoption decision of a national government as a threshold decision rule: adopt policy $i$ if $\Pr[\text{IPI}(i) > 0 | \text{evidence}] > p^*$. Show that S1 (demonstrated welfare improvement) is necessary for $\Pr > p^*$ to be achievable with finite experimentation evidence.
(b) Model the replication decision of a second city: adopt if the expected replication cost $C^{\text{rep}} < E[\text{benefit}]$. Show that S2 (transferable template) is necessary for $C^{\text{rep}}$ to be sufficiently low.
(c) Model the political action decision of the champion network as a public goods game. Show that S3 (political champion network) is necessary for collective political action to occur — without a network, individual champions free-ride on others' advocacy.
(d) Model the information diffusion dynamics. Show that S4 (policy learning infrastructure) is necessary for diffusion to exceed the viral threshold ($R_0 > 1$) — without a shared knowledge commons, each potential adopter must independently discover the institutional template, slowing diffusion below the viral rate.

**★ 42.5** Design a CR-SEZ for a context of your choice and formally analyze its scaling potential.

Choose one of: a rural agricultural district in transition from conventional to regenerative farming; an urban district in a post-industrial city seeking cooperative enterprise development; a port or logistics zone transitioning to circular economy principles; or a digital industry cluster developing cooperative platform alternatives.

(a) **CR-SEZ design:** For your chosen context, specify which elements of Definition 42.3's six-part institutional package are most relevant and how they should be adapted. What additional elements are needed for your specific context?

(b) **MLP assessment:** Assess landscape pressure, niche maturity, and regime window for your context. Based on the MLP conditions, what is the probability of successful scaling within 5 years? Within 10 years?

(c) **Scaling potential analysis:** Apply Proposition 42.2 and Theorem 42.1. Does your proposed CR-SEZ satisfy all four scaling conditions? If not, which conditions are missing and how would you address them?

(d) **Thompson Sampling design:** Specify 3–5 institutional variants of your CR-SEZ design as "arms" in a policy MAB. What outcome metrics would you use as the welfare signal $R(a_k, t)$? Design the experimentation protocol (number of sites, measurement frequency, evaluation timeline).

**★★ 42.6** Conduct a full evaluation of a municipalist policy experiment in a city of your choice.

**Selection criteria:** Choose a city that: (i) has a documented cooperative economy or commons policy agenda; (ii) has outcome data spanning at least 4 years; and (iii) has comparison cities or regions available for difference-in-differences analysis.

**Suggested options:** Ghent (Belgium) — cooperative procurement and food commons; Amsterdam (Netherlands) — donut city strategy; Cleveland (Ohio, USA) — Evergreen Cooperatives; Medellín (Colombia) — innovation districts and cooperative enterprise; Marinaleda (Spain) — cooperative village economy.

(a) **Policy documentation:** Catalog the cooperative-regenerative policies adopted by your chosen city, organized by the eleven instruments of Chapter 41. Which instruments were deployed? What design choices were made?

(b) **MLP assessment:** For each major policy: assess landscape pressure, niche maturity, and regime window at the time of adoption. Does the MLP framework predict the outcomes observed?

(c) **Formal outcome evaluation:** Using difference-in-differences or synthetic control methods (with a matched comparison city or region), estimate the causal effect of your chosen city's cooperative policies on: median household income, cooperative employment share, housing cost burden, carbon footprint, and an available welfare measure.

(d) **Scaling assessment:** Apply Theorem 42.1 to assess the scaling potential of each policy. Which policies have been replicated? Which have not? Does the theorem's four conditions predict the actual replication pattern?

(e) **Design recommendations:** Based on your evaluation, specify three changes to the city's cooperative economy policies that would: (i) improve welfare outcomes in the city; (ii) increase the scaling potential for national diffusion; and (iii) be politically feasible within the city's current coalition.

---

*Part IX opens with Chapter 43. The formal framework is complete; the transition theory is specified; the empirical cases confirm the model's core predictions while identifying its limits. The final Part turns from analysis to synthesis: the open questions that define the research frontier, the capstone project that applies the full toolkit, and the conclusion that situates the cooperative-regenerative economics framework within the longer arc of economic thought.*
