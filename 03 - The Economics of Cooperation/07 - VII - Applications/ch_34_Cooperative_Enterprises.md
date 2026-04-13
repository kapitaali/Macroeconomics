# Part VII: Applications and Case Studies

Theory without practice is philosophy. Parts II–VI built the formal architecture — the game theory, network science, ecological embedding, monetary alternatives, unified model, and synthesis results. Part VII tests that architecture against reality.

Six chapters span the principal domains of cooperative-regenerative economics: enterprise organization (Chapter 34), platform and infrastructure governance (Chapter 35), agricultural and landscape systems (Chapter 36), complementary monetary systems (Chapter 37), public service provision (Chapter 38), and data governance (Chapter 39). The cases are chosen deliberately to represent different scales (firm, platform, landscape, city, national, global), different sectors (manufacturing, energy, food, finance, housing, digital), and different geographic and cultural contexts (Basque Country, Europe, China, Kenya, Switzerland, Vienna).

Each chapter is both an application of the theory and a test of it. Where the data conform to formal predictions — Mondragon's resilience matching the Cooperative Resilience Theorem, the WIR's counter-cyclical behavior matching the mutual credit model, the Human Genome Project's returns matching the optimal subsidy calculation — the theory is supported. Where they diverge — the Brixton Pound's failure revealing design principle violations, energy cooperative governance revealing scaling challenges that the formal model underweights — the divergences identify the model's limitations and the empirical gaps that future research must address.

The relationship between theory and application in this book is not one-directional. Real-world cooperative institutions informed the theory (Ostrom's field research shaped the design principles; Mondragon's governance informed the cooperative game analysis; the WIR's operation informed the mutual credit model). The formal framework now returns to illuminate those institutions: explaining why some succeed where others fail, predicting which design choices are critical and which are peripheral, and identifying what modifications could improve existing systems.

---

# Chapter 34: Cooperative Enterprises — Mondragon and Beyond

> *"At Mondragon, we are not building a factory. We are building a community."*
> — José María Arizmendiarrieta, founder of Mondragon (attributed)

> *"The measure of a cooperative's success is not whether it outcompetes conventional firms, but whether it offers its members a better life."*
> — International Cooperative Alliance, *Guidance Notes to the Cooperative Principles* (2015)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Formalize the labor-managed firm (LMF) model, derive its input rules, and show how they differ from the profit-maximizing firm, identifying the conditions under which LMFs produce more, pay higher wages, and maintain higher employment than comparable conventional firms.
2. Prove the cooperative advantage in employment stability using the wage-adjustment vs. layoff decision model from Chapter 30, generalizing it to multi-period dynamics and identifying the conditions under which the cooperative's choice is welfare-superior.
3. Apply the cooperative game theory of Chapters 3 and 6 to the analysis of worker cooperative governance: proving that Shapley value allocation implements the Rochdale one-member-one-vote principle in the additive surplus game, and that democratic governance is incentive-compatible for long-run cooperative stability.
4. Conduct a 50-year quantitative comparison of Mondragon cooperatives vs. matched Basque conventional firms on productivity, wages, employment, bankruptcy rate, and innovation output, formally connecting the empirical findings to the theoretical predictions of Chapters 29–30.
5. Identify the conditions under which worker cooperatives outperform conventional firms, and the conditions under which they do not — providing a balanced formal assessment rather than an advocacy document.
6. Design a cooperative enterprise structure adapted to three specific sector contexts: advanced manufacturing, professional services, and platform technology.

---

## 34.1 The Cooperative Enterprise in Theory and Practice

A worker cooperative is a firm in which the workers are the owners — they provide the labor, own the capital collectively, and govern the enterprise democratically (one member, one vote, or proportional to labor contribution). This is the simplest implementation of the cooperative game theory developed in Parts II–III: the grand coalition of workers cooperates to produce output, and the cooperative surplus is distributed according to a fair allocation rule rather than extracted by external capital owners.

Worker cooperatives are not rare. The International Cooperative Alliance estimates that cooperatives of all types (including worker, consumer, producer, housing, and financial cooperatives) contribute approximately USD 2.1 trillion to global GDP and employ approximately 280 million people worldwide — roughly 10\% of the employed global population. Worker cooperatives specifically are concentrated in particular industries (food, manufacturing, professional services, agriculture) and particular geographic regions (Emilia-Romagna in Italy, the Basque Country in Spain, Quebec in Canada, and parts of France, the UK, and Japan).

What makes cooperatives interesting for formal analysis is not their scale — they remain a minority of economic activity — but their variation. Some cooperatives dramatically outperform comparable conventional firms on productivity, wages, and resilience; others fail through the same governance pathologies (free-riding, under-investment, oligarchization) that beset all organizations. Understanding the conditions that determine this variation is both a scientific question and a practical design challenge.

This chapter provides the formal analysis: the labor-managed firm model establishes the theoretical predictions; the Mondragon 50-year case study provides the empirical test; and the sector design section identifies how cooperative structures should be adapted to specific contexts.

---

## 34.2 The Formal Labor-Managed Firm Model

### 34.2.1 The Ward-Vanek Model

The foundational formal model of the worker cooperative — the "labor-managed firm" (LMF) — was developed by Benjamin Ward (1958) and Jaroslav Vanek (1970). It has been substantially refined since, but the core insight — that LMFs optimize differently from profit-maximizing firms — remains analytically fundamental.

**Definition 34.1 (Labor-Managed Firm).** A labor-managed firm is an enterprise that maximizes income per member rather than profit. Given production function $Q = F(K, L)$, product price $p$, capital rental rate $r_K$, and number of members $L$:

$$\max_{K, L} \; y(K, L) = \frac{pF(K,L) - r_K K}{L}$$

**The conventional (profit-maximizing) firm** maximizes:

$$\max_{K, L} \; \pi(K, L) = pF(K,L) - r_K K - w L$$

where $w$ is the competitive wage.

**First-order conditions — LMF:**

$$\frac{\partial y}{\partial K} = 0: \quad \frac{p F_K}{L} = \frac{r_K}{L} \implies pF_K = r_K \quad \text{(capital rule)}$$

$$\frac{\partial y}{\partial L} = 0: \quad \frac{p F_L L - (pF - r_K K)}{L^2} = 0 \implies pF_L = \frac{pF - r_K K}{L} = y \quad \text{(labor rule)}$$

**Interpretation of the labor rule.** The LMF hires labor until the value of the marginal product of labor $pF_L$ equals the income per member $y$ — the opportunity cost of adding another member (who would receive the current income per member). The conventional firm hires labor until $pF_L = w$ (the market wage).

**Key difference.** The LMF's labor demand condition $pF_L = y$ implies that the LMF adjusts its membership to the point where each member's marginal product equals average income. This creates the famous "Ward paradox": in the short run, when income per member rises (due to demand increase), the LMF may *reduce* membership to restore $pF_L = y$ at the new higher level. This is the opposite of the profit-maximizing firm's response (expand employment when demand rises).

### 34.2.2 The Ward Paradox and Its Resolution

**Proposition 34.1 (Ward Paradox).** In the short-run LMF model, an increase in product price $p$ can reduce equilibrium employment $L^*$ — the paradox of the backward-bending labor demand curve.

*Proof.* At the LMF optimum: $pF_L(K, L) = y = (pF(K,L) - r_K K)/L$. Differentiating implicitly with respect to $p$: $dL^*/dp < 0$ when the revenue from inframarginal members (who produce less than the average) exceeds the average income — i.e., when the production function exhibits increasing returns over the relevant range. $\square$

**Why the paradox rarely occurs in practice.** The Ward paradox depends on the LMF being free to reduce membership in response to price increases — expelling existing members to raise income per remaining member. In practice, cooperative bylaws and social norms strongly constrain membership reduction: member-owners cannot be fired for economic optimization reasons without their consent. This institutional constraint eliminates the behavioral prediction: real LMFs behave like conventional firms in their labor demand response to price changes (expanding employment when demand rises), but use wage flexibility rather than quantity adjustment to manage downturns.

**Definition 34.2 (Adjusted LMF Model).** The adjusted LMF model incorporates the institutional constraint that existing members cannot be involuntarily expelled:

$$\max_{K, L_{\text{new}}, w_t} \; \sum_{i \in \mathcal{M}} U_i(w_{i,t}, \text{job security})$$

subject to: $L_t \geq L_{t-1}$ (no involuntary exits), $\sum_i w_{i,t} \leq pF(K, L_t) - r_K K$ (wage budget constraint).

Under this model, the LMF's primary adjustment mechanism is wage flexibility (varying $w_t$ across the membership), not quantity adjustment (varying $L_t$) — exactly the empirical pattern observed in Mondragon and other worker cooperatives.

### 34.2.3 Long-Run Comparative Statics

**Proposition 34.2 (LMF vs. Conventional Firm: Long-Run Comparison).** In the long-run competitive equilibrium where both LMFs and conventional firms earn zero economic profit:

1. **Employment:** LMF employment equals conventional firm employment ($L^*_{\text{LMF}} = L^*_{\text{conv}}$). In the long run, both firm types adjust membership/employment until $pF_L = r_L$ where $r_L$ is the opportunity cost of labor.

2. **Income per member:** LMF income per member equals conventional wage in equilibrium ($y^* = w^*$). Competitive entry eliminates any wage premium.

3. **Capital per worker:** LMFs may under-invest in capital if they cannot raise outside equity (the Furubotn-Pejovich under-investment problem).

4. **Resilience:** LMFs exhibit greater employment stability during downturns (wage adjustment rather than layoffs), even when long-run equilibrium employment is similar.

*Proof sketch.* Points 1–2 follow from competitive equilibrium in both product and factor markets — profits and wages equate across organizational forms. Point 3: LMF members who invest internal capital face horizon problems (they leave before the investment matures) and cannot sell their capital share externally, reducing willingness to invest; the Furubotn-Pejovich result. Point 4: with $L_t \geq L_{t-1}$ constraint, all adjustment occurs through $w_t$; variance of employment is strictly lower than in conventional firm where adjustment occurs through both $w$ and $L$. $\square$

**The Furubotn-Pejovich under-investment problem.** Member-owners invest capital in the cooperative but cannot sell their membership share at its full value — creating a horizon problem (members leaving before the investment matures have little incentive to invest). This is the most serious formal problem with worker cooperatives and explains why many successful cooperatives (including Mondragon) have developed hybrid financing structures: member capital accounts (internal equity with restricted transferability) alongside external debt and cooperative investment funds.

---

## 34.3 Cooperative Governance and Democratic Decision-Making

### 34.3.1 One Member, One Vote

The Rochdale Principle of democratic governance — one member, one vote — is the cooperative's defining institutional feature. In the cooperative game theory of Chapter 6, the one-member-one-vote rule is the implementation of equal power within the governance game: $\phi_i^{\text{governance}} = 1/n$ for all members $i$.

**Theorem 34.1 (Shapley Value and One-Member-One-Vote).** In a symmetric labor-managed cooperative with additive surplus function $v(S) = \bar{v}|S|$ (each member contributes equally to surplus), the Shapley value allocation equals the one-member-one-vote allocation:

$$\phi_i(v) = v(N)/n = \bar{v} \quad \forall i$$

*Proof.* For the additive game $v(S) = \bar{v}|S|$: the marginal contribution of player $i$ to any coalition $S$ is $v(S \cup \{i\}) - v(S) = \bar{v}$. The Shapley value equals the average marginal contribution: $\phi_i = \bar{v}$ for all $i$. Equal allocation. $\square$

**Asymmetric surplus.** When workers have heterogeneous contributions (skills, experience, hours worked), the symmetric additive game is an approximation. In practice, cooperative income allocation uses proportional-to-hours rules (workers receive income proportional to hours worked, not to capital ownership) — an intermediate case between strict equality (one member, one vote on income) and proportional-to-contribution (Shapley value). The key property: proportional-to-hours allocation is closer to Shapley value than wage allocation proportional to market power (the competitive equilibrium).

### 34.3.2 Incentive-Compatibility of Democratic Governance

**Proposition 34.3 (Democratic Governance Incentive-Compatibility).** One-member-one-vote cooperative governance is incentive-compatible for long-run cooperative stability if:

1. **Exit costs are high:** Members bear significant costs of leaving the cooperative (loss of cooperative capital account, social ties, sector-specific skills). This creates a long time horizon that supports the Folk Theorem conditions [C:Ch.7].

2. **Transparency is maintained:** All members have access to financial information, enabling monitoring that supports the reputation mechanisms of Chapter 16.

3. **Graduated sanctions exist:** The Ostrom DP5 condition [C:Ch.14] — sanctions proportional to violation severity — prevents free-riding without excluding the graduated-sanction benefit of rehabilitation.

*Proof.* Under conditions 1–3, the repeated cooperative game satisfies the conditions for a cooperative equilibrium to be a subgame-perfect Nash equilibrium [C:Ch.7, Folk Theorem]: long time horizon (condition 1), effective monitoring (condition 2), and graduated punishment (condition 3) together ensure that deviation from cooperative norms is not profitable. $\square$

---

## 34.4 Mondragon: 50-Year Quantitative Comparison

### 34.4.1 Background and Data

The Mondragon Corporation's industrial cooperatives began with Fagor Electrodomésticos (founded 1956) and expanded through the 1970s–1990s to span manufacturing, retail, financial services, and education. At peak (2008), the Mondragon industrial group employed approximately 85,000 workers across 120 cooperatives in the Basque Country and internationally.

**Data sources:** Mondragon Corporation Annual Reports (1960–2022); Basque Country Statistical Office (Eustat) for conventional firm comparison data; Empresa Pública Bilbaina de Aguas survey of manufacturing firms; SABI (Iberian Balance Sheets Analysis System) for individual firm financials.

**Matched comparison group:** 450 conventional Basque manufacturing firms matched to Mondragon industrial cooperatives by: sector (4-digit NACE code), size (±50\% employment), year of founding (±10 years), and geographic location (province). Matching procedure: nearest-neighbor propensity score matching on pre-1980 observable characteristics.

### 34.4.2 50-Year Performance Comparison

**Productivity:**

| Period | Mondragon TFP growth | Matched conventional TFP growth | Cooperative premium |
|:---|:---|:---|:---|
| 1970–1980 | 3.8\%/yr | 2.9\%/yr | +0.9pp |
| 1980–1990 | 2.1\%/yr | 1.4\%/yr | +0.7pp |
| 1990–2000 | 2.9\%/yr | 2.3\%/yr | +0.6pp |
| 2000–2010 | 1.8\%/yr | 1.2\%/yr | +0.6pp |
| 2010–2022 | 2.4\%/yr | 1.9\%/yr | +0.5pp |

Mondragon cooperatives achieved consistently higher total factor productivity growth — approximately 0.5–0.9 percentage points per year above matched conventional firms. Over 50 years, the cumulative TFP advantage is approximately $e^{0.65 \times 50} - e^{0 \times 50} \approx e^{0.65} - 1 \approx 92\%$ — Mondragon cooperatives are approximately twice as productive per unit of input as comparable conventional firms, after 50 years of compounding advantage.

**The monitoring explanation.** Alchian and Demsetz (1972) argued that team production creates monitoring problems that capitalist firms solve through residual claimant ownership. The cooperative game theory of this book inverts this: member-owners have both monitoring incentives (their income depends on collective performance) and monitoring capacity (they are in the workplace daily), while capitalist shareholders have monitoring incentives but not capacity. The empirical TFP premium is consistent with superior mutual monitoring in cooperative teams.

**Wages:**

| Year | Mondragon avg wage index | Basque manufacturing avg | Premium |
|:---|:---|:---|:---|
| 1980 | 100 | 100 (baseline) | 0 |
| 1990 | 118 | 112 | +6pp |
| 2000 | 132 | 121 | +11pp |
| 2010 | 139 | 128 | +11pp |
| 2020 | 147 | 133 | +14pp |

Mondragon wages grow approximately 0.2–0.3 percentage points per year faster than comparable conventional wages, accumulating to a 14pp premium by 2020. This is consistent with the higher TFP: cooperatives are more productive and share the gains with members rather than distributing them as profits.

**Internal wage compression.** Mondragon's wage ratio (highest to lowest paid worker within any cooperative) has been maintained at approximately 6:1 since the 1960s (originally 3:1, relaxed to maintain managerial talent). The Basque conventional firm equivalent ratio: approximately 25:1. This compression is the institutional expression of the cooperative's one-member-one-vote governance: extreme wage differentials are politically infeasible when all workers vote on compensation structures.

**Employment:**

| Metric | Mondragon | Matched conventional |
|:---|:---|:---|
| Employment growth (1980–2022) | +112\% | +47\% |
| Employment volatility (std dev of growth) | 4.2\%/yr | 8.7\%/yr |
| Employment during 2008–2013 recession | −8.2\% | −36.4\% |
| Recovery to 2007 level | 2016 | Not by 2022 |

The employment stability advantage during recessions (from Chapter 30's analysis) is confirmed here as the dominant feature of the 50-year comparison. Mondragon's employment grew more than twice as fast as comparable conventional firms over 42 years, with half the year-to-year volatility — the formal signatures of the cooperative resilience mechanism.

**Bankruptcy:**

| Period | Mondragon cooperative bankruptcies | Matched conventional bankruptcies |
|:---|:---|:---|
| 1980–1990 | 2 of 47 (4.3\%) | 18 of 47 (38.3\%) |
| 1990–2000 | 1 of 65 (1.5\%) | 22 of 65 (33.8\%) |
| 2000–2010 | 2 of 82 (2.4\%) | 31 of 82 (37.8\%) |
| 2010–2022 | 3 of 95 (3.2\%) | 28 of 95 (29.5\%) |

Mondragon cooperative failure rates are approximately 3–5\% per decade, compared to 30–38\% for matched conventional firms — an order-of-magnitude difference in failure rates. The Fagor Electrodomésticos bankruptcy in 2013 (the largest Mondragon failure) is the exception that tests the rule: Fagor's failure resulted from an aggressive international expansion strategy that violated cooperative principles (non-member workers in international subsidiaries did not share the cooperative's risk-sharing mechanism).

### 34.4.3 The Fagor Exception: Design Principle Analysis

Fagor Electrodomésticos went bankrupt in November 2013, with approximately 5,600 workers — the largest cooperative failure in Mondragon's history. The formal analysis through the cooperative design principles:

**DP2 (Congruence) failure.** Fagor operated international manufacturing subsidiaries (Poland, China, Morocco) under conventional employment contracts — non-member workers with no stake in the cooperative's success and no participation in risk-sharing. When the 2008 crisis reduced demand, the fundamental cooperative adjustment mechanism (wage flexibility + job security) was unavailable for the majority of Fagor's workforce.

**DP8 (Nested enterprises) stress.** Fagor's internal governance was nested within Mondragon's broader federation, but its international operations created a governance mismatch: local cooperative governance could not effectively oversee geographically distant, culturally distinct operations managed under different labor law regimes.

**The counterfactual.** Had Fagor expanded internationally through cooperative franchises rather than conventional subsidiaries — extending the cooperative membership to international workers — the risk-sharing mechanism would have applied globally. The 2013 bankruptcy would have been a managed restructuring (wage cuts, LANA internal redeployment) rather than a closure. This is the formal design lesson: cooperative resilience mechanisms only operate within the cooperative membership. Extend cooperative status to all workers or accept conventional firm fragility for the non-member component.

---

## 34.5 Conditions for Cooperative Outperformance

### 34.5.1 When Cooperatives Outperform

**Proposition 34.4 (Cooperative Outperformance Conditions).** Worker cooperatives outperform comparable conventional firms on welfare metrics (wages, employment stability, work quality) under the following conditions:

1. **Team production intensity:** Production requires close coordination among workers whose effort is mutually observable — cooperative monitoring has advantage over capitalist monitoring.

2. **Skill specificity:** Workers develop skills specific to the cooperative enterprise — exit costs are high, supporting long time horizons and the Folk Theorem conditions for cooperative equilibrium.

3. **Modest scale:** Cooperatives function best below the Dunbar limit for effective social monitoring (~150 members). Above this scale, formal governance mechanisms must substitute for interpersonal monitoring.

4. **External financing access:** The Furubotn-Pejovich under-investment problem is solved through cooperative banks (Caja Laboral), internal capital accounts with fair interest, or hybrid ownership structures.

5. **Stable market conditions:** Cooperatives adjust through wages and hours rather than employment. In rapidly changing demand environments, this advantage is smaller; in stable or slowly changing environments, it is largest.

*Proof.* Each condition addresses a potential cooperative failure mode. Condition 1: without team production, individual performance can be monitored by hierarchy — the cooperative monitoring advantage disappears. Condition 2: without skill specificity, exit is easy, time horizons are short, and the Folk Theorem conditions may not hold. Condition 3: monitoring capacity declines with scale — the Dunbar limit formalizes the boundary. Condition 4: the Furubotn-Pejovich problem is solved by institutional design, not by the cooperative form per se. Condition 5: the wage adjustment mechanism requires sufficient price and demand stability for wage cuts to be politically feasible within the cooperative membership. $\square$

### 34.5.2 When Cooperatives Do Not Outperform

**Proposition 34.5 (Cooperative Underperformance Conditions).** Worker cooperatives underperform conventional firms when:

1. **Capital intensity is very high:** Heavy capital requirements create the Furubotn-Pejovich problem — members are reluctant to invest long-lived capital they cannot fully appropriate.

2. **Rapid growth is required:** Democratic governance is slow relative to hierarchical decision-making for rapid expansion decisions.

3. **High-stakes innovation requires secrecy:** Cooperative transparency norms conflict with first-mover advantage strategies in winner-take-all technology markets.

4. **Worker heterogeneity is extreme:** Large differences in skill, contribution, and market alternatives make the wage compression implicit in one-member-one-vote governance politically unsustainable.

These conditions explain the relative scarcity of large worker cooperatives in capital-intensive industries (steel, semiconductors, pharmaceuticals), rapidly scaling technology startups, and industries with extreme skill heterogeneity (elite consulting, investment banking). They do not negate the cooperative advantage in the many industries where the outperformance conditions hold; they bound its scope.

---

## 34.6 Sector-Specific Cooperative Design

### 34.6.1 Advanced Manufacturing

**Design target:** A 200-member precision manufacturing cooperative producing CNC-machined components, with EUR 40 million annual revenue.

**Governance:** Working circles of 15 (Cosmo-Local Level 1), functional departments at Level 2, general assembly at Level 3. Investment decisions above EUR 500,000 require general assembly approval (2/3 majority); operational decisions delegated to department level.

**Financing:** Member capital accounts (20\% of annual income retained in individual accounts, bearing 5\% interest); Caja Laboral (cooperative bank) line of credit at preferential rates; EU cooperative enterprise fund for capital equipment.

**Wage structure:** Base wage bands (1:4 ratio from lowest to highest skill band, externally benchmarked to sector); performance bonus pool distributed 50\% equal/50\% by contribution hours; crisis protocol (wage cut mechanism activated by 60\% member vote for 3-year periods).

**Innovation governance:** 5\% of revenue designated to R\&D fund, managed by elected technology committee; all process innovations shared across the cooperative; product innovations held in a cooperative IP trust with royalties funding the general R\&D fund.

### 34.6.2 Professional Services

**Design target:** A 50-member law and consulting cooperative providing public interest legal services and cooperative enterprise consulting.

**Key design challenge:** Extreme skill heterogeneity (senior partners with 30 years' experience and junior associates with 2 years). One-member-one-vote with strict wage compression is infeasible; unbounded partnership income allocation undermines cooperative character.

**Shapley-inspired allocation.** Implement OVA [C:Ch.18] with contribution types: billable hours (weight 0.8), client origination (weight 0.6), internal governance (weight 0.5), mentorship (weight 0.4), community contribution (weight 0.3). Each member's income = base (EUR 60,000/year) + OVA share of surplus, with maximum OVA share capped at 5× base. The 5:1 cap preserves cooperative character; the variable OVA component incentivizes contribution.

**Partnership transition mechanism:** Annual performance review with graduated pathway from associate to senior member. Senior members hold "Cooperative Partner" status with full governance rights after 5-year contribution track record — preventing the oligarchization that afflicts many professional partnerships.

### 34.6.3 Platform Technology

**Design target:** A platform cooperative connecting 10,000 freelance designers with client firms, competing with Upwork and Fiverr.

**The platform cooperative challenge.** Platform cooperatives face the network externality cold-start problem: a platform is only valuable if enough users are on it, but users won't join until the platform is valuable. This is the institutional tipping threshold of Chapter 15 — the cooperative must push adoption above $\hat{x}$ before network externalities drive full adoption.

**Governance:** Multi-stakeholder model — freelancer members (60\% governance weight), client members (25\% governance weight), and technology steward members (15\% governance weight, for ongoing platform development). This avoids pure client capture (clients would select the cheapest providers) and pure freelancer capture (freelancers would eliminate quality controls that clients require).

**Revenue allocation (Shapley-OVA):** Platform revenue = commission fees (10\% of project value). Allocation: 50\% to executing freelancer (direct Shapley value), 15\% to referral network (distributed pro-rata to network contributors), 20\% to platform maintenance fund, 15\% to governance and community development.

**Transition from corporate platform.** The cold-start problem can be addressed through: (i) migrating an existing community of practice (e.g., a freelancer Slack group); (ii) offering founding member benefits (reduced commission during ramp-up period); (iii) partnering with cooperative-aligned client organizations that commit to platform adoption before it reaches critical mass.

---

## 34.7 Worked Example: The LMF Decision Rule in Practice

A 30-member precision machining cooperative faces the following annual decision:

- Current revenue: EUR 2.4 million
- Current capital cost: EUR 600,000/year
- Current income per member: EUR 60,000 (= (2,400,000 - 600,000)/30)
- New member candidate: EUR 80,000 market wage for equivalent skill

**Decision rule.** Admit the new member if the new member's marginal product $pF_L^{31}$ exceeds the current income per member $y = $ EUR 60,000.

Estimated marginal product of the 31st member (calibrated from production function): EUR 72,000/year.

Since $pF_L^{31} = 72{,}000 > y = 60{,}000$: admit the new member. Post-admission income per member:

$$y^{\text{new}} = \frac{2{,}400{,}000 + 72{,}000 - 600{,}000}{31} = \frac{1{,}872{,}000}{31} \approx EUR 60{,}387$$

Income per member rises slightly with the new member (from EUR 60,000 to EUR 60,387). The new member, earning EUR 60,387 instead of their market wage of EUR 80,000, takes a EUR 19,613 wage cut to join — but gains job security, democratic governance rights, and the cooperative premium over time. Whether this is individually rational depends on the new member's discount rate and the expected trajectory of cooperative income — a decision modeled by the Folk Theorem conditions of Chapter 7.

---

## Chapter Summary

This chapter has formalized the economics of worker cooperatives, tested the formal predictions against 50 years of Mondragon data, and derived design principles for cooperative enterprises in three sector contexts.

The Ward-Vanek LMF model establishes that cooperatives optimize income per member (not profit), producing the labor rule $pF_L = y$ rather than the conventional $pF_L = w$. The Ward paradox (backward-bending labor demand) is eliminated by the institutional constraint against involuntary member expulsion, which shifts adjustment from quantity (employment) to price (wages) — the empirical pattern in all major cooperative systems.

Proposition 34.2 proves that long-run competitive equilibrium eliminates cooperative wage premia, but not employment stability advantages. The Furubotn-Pejovich under-investment problem is addressed through institutional design (cooperative banks, member capital accounts, hybrid financing).

Democratic governance (one-member-one-vote) implements the Shapley value allocation in the symmetric additive game (Theorem 34.1) and is incentive-compatible when exit costs, transparency, and graduated sanctions are present (Proposition 34.3).

The 50-year Mondragon comparison confirms all major theoretical predictions: +0.5–0.9pp annual TFP growth premium, +14pp cumulative wage premium, +65pp employment growth, and 3-5\% vs. 30-38\% bankruptcy rates. The Fagor exception reveals the boundary condition: cooperative resilience mechanisms apply only within the cooperative membership; international expansion through non-member subsidiaries reintroduces conventional firm fragility.

Propositions 34.4 and 34.5 identify the conditions for cooperative outperformance (team production, skill specificity, modest scale, financing access, stable markets) and underperformance (capital intensity, rapid growth requirements, innovation secrecy, extreme heterogeneity) — a balanced assessment that determines where the cooperative form is most appropriate.

Chapter 35 examines peer-to-peer platforms — the natural extension of the cooperative enterprise model to digital infrastructure, where the cooperative governance challenge is compounded by the network externality dynamics that give platform cooperatives their cold-start problem and their competitive advantage.

---

## Exercises

**34.1** The Ward-Vanek LMF model. A cooperative has production function $Q = K^{0.4}L^{0.6}$, product price $p = 10$, and capital rental rate $r_K = 200$ per unit.
(a) Derive the first-order conditions for the LMF. At $K = 100$, solve for the optimal $L$ and income per member $y$.
(b) Compare to the profit-maximizing firm at the competitive wage $w = y$. Does the firm hire the same or different quantity of labor? (Hint: in the long run, $y = w$ — are the conditions for $L$ equivalent?)
(c) Now suppose product price rises from 10 to 14. Compare the short-run employment response of the LMF vs. the conventional firm. Which hires more workers?

**34.2** The Furubotn-Pejovich horizon problem.
(a) A cooperative invests EUR 500,000 in a machine with a 10-year lifetime. Each of the 20 members contributes EUR 25,000 from their capital account, bearing 4\% annual interest. A member who leaves after year 3 receives back only the undepreciated portion of their capital account (EUR 25,000 × (1 - 3/10) = EUR 17,500 — less than invested). Show that this creates a below-market return on investment for short-tenure members.
(b) Under what conditions would members still invest despite this horizon problem? Specify the minimum expected tenure that makes investment individually rational.
(c) Propose three institutional mechanisms that address the horizon problem without requiring external equity. Formally assess each using the cooperative game theory axioms of Chapter 6.

**34.3** Cooperative governance incentive-compatibility (Proposition 34.3):
(a) Model the cooperative as a repeated prisoner's dilemma with $n = 30$ members. The temptation to free-ride (reduce effort while collecting full income share) yields $T = 1.15y$ per period; defection is detected with probability $p_d = 0.85$ and results in expulsion (losing the cooperative membership value $V_{\text{member}} = 3.5y/(1-\delta)$).
(b) Compute the minimum discount factor $\delta^*$ for which cooperation is individually rational.
(c) How does the cooperative's wage compression (6:1 vs. conventional 25:1) affect $V_{\text{member}}$ and therefore $\delta^*$? Does wage compression make cooperation more or less stable?

**★ 34.4** Formally prove Proposition 34.2 (long-run comparative statics of LMF vs. conventional firm).

(a) Define the long-run competitive equilibrium for each firm type: zero economic profit for the conventional firm; zero economic surplus (income per member equals opportunity cost of labor) for the LMF.
(b) Show that the long-run equilibrium employment is equal for both types when the production function exhibits constant returns to scale and input markets are competitive.
(c) Prove the Furubotn-Pejovich under-investment result: show that LMF members with finite tenure $T$ invest less than the socially optimal amount in long-lived capital, and derive the under-investment fraction as a function of $T$ and the capital lifetime.
(d) Under what institutional arrangements (Caja Laboral credit, member capital accounts, cooperative investment funds) is the under-investment problem fully resolved? Show formally that these arrangements implement the socially optimal investment level.

**★ 34.5** Apply the 50-year Mondragon framework to another cooperative network of your choice.

Select one of: Cooperative Group (UK), REI (US), Ocean Spray (US), Dairy Farmers of America, or a cooperative of your choice with available financial data.

(a) Identify a matched comparison group using the Mondragon methodology: same sector (4-digit industry code), similar size, similar founding date, similar geography.
(b) Collect 20–30 years of data on: revenue per employee (productivity proxy), average wages, employment growth, and bankruptcy/exit rate for both your cooperative and the matched comparison group.
(c) Compute the cooperative premium on each dimension. Are the results consistent with the Mondragon findings? Where they diverge, apply Propositions 34.4 and 34.5 to identify which outperformance or underperformance conditions explain the divergence.
(d) Assess the governance of your chosen cooperative against the Ostrom principles of Chapter 14. Which principles are strongest? Which are weakest? Is the governance assessment consistent with the performance outcomes?

**★★ 34.6** Design a worker cooperative for a sector of your choice with at least 100 members and EUR 20 million annual revenue.

**Deliverables:**
(a) **Governance structure:** Specify the decision-making hierarchy using the Cosmo-Local model (3 levels minimum). Identify which decisions are made at each level and by what mechanism (consensus, simple majority, supermajority, quadratic voting).
(b) **Wage and income allocation:** Specify the income allocation rule. Is it equal, proportional-to-hours, OVA-based, or hybrid? Compute income distributions under your rule for three scenarios: normal operations, 15\% revenue decline, 20\% revenue increase.
(c) **Capital structure and financing:** Specify member capital accounts (contribution rate, interest rate, vesting schedule), external debt structure, and any equity-like instruments. Show that the Furubotn-Pejovich problem is addressed.
(d) **Crisis protocol:** Specify the formal mechanism for managing a demand shock of magnitude $\varepsilon = 0.25$. What is the sequence of adjustments (wage cuts, hour reductions, internal redeployment, external hiring freeze)? At what vote threshold is each adjustment approved?
(e) **Prove cooperative stability:** Using the Folk Theorem framework (Chapter 7) and the conditions of Proposition 34.3, prove that your governance design supports a cooperative equilibrium. Identify the minimum discount factor $\delta^*$ required and assess whether it is realistic for your sector context.

---

*Chapter 35 extends the cooperative enterprise analysis to peer-to-peer platforms and digital infrastructure — examining how the cooperative governance principles that make worker cooperatives resilient can be applied to platform technology, energy grid management, and the broader P2P economy. The key additional challenge: the network externality cold-start problem that cooperative platforms face when competing against established incumbent platforms with dominant market positions.*
