# Chapter 33: The Role of Digital Commons — Knowledge, Software, and Data as Shared Resources

> *"Information wants to be free. Information also wants to be expensive. That tension will not go away."*
> — Stewart Brand, *The Media Lab* (1987)

> *"Software is eating the world. Open-source software is eating software."*
> — Marc Andreessen (2011), paraphrased; and the open-source community

## Learning Objectives

By the end of this chapter, you should be able to:

1. Formally specify digital commons goods as non-rival and partially excludable, prove that their optimal pricing under standard welfare theory is zero (marginal cost pricing), and derive the dynamic efficiency implications of under-provision through intellectual property restrictions.
2. Estimate the social value of digital commons using willingness-to-accept and consumer surplus methods, demonstrating that GDP systematically under-counts the value of non-market digital goods by an order of magnitude.
3. Apply the Ostrom design principles to digital commons (Wikipedia, Linux, OpenStreetMap), scoring each principle and formally analyzing the stability of digital commons governance against the predictions of Chapter 14.
4. Formalize data as a partially rival common-pool resource, specify the governance framework for data trusts and data cooperatives using Shapley value contribution allocation, and analyze GDPR as an institutional experiment in data commons governance.
5. Derive the formal conditions under which AI development benefits versus undermines the digital commons, and specify governance structures for open-source AI development consistent with the cooperative-regenerative framework.
6. Estimate the social value of the Python ecosystem, compute the optimal public subsidy implied by welfare theory, and compare to actual public and private investment levels.

---

## 33.1 The Information Paradox

Information has always been awkward for economics. The classical theory of value — whether labor value or marginal utility — was built around goods that are scarce: a bushel of wheat consumed by one person cannot be consumed by another; an hour of a craftsman's labor spent on one task is unavailable for another. Scarcity is the condition under which markets allocate efficiently. Remove scarcity, and the market mechanism becomes not merely imperfect but structurally inappropriate.

Knowledge goods — information, software, scientific discoveries, artistic works, and data — are non-rival: one person's use of a mathematical theorem, a software library, or a genome sequence does not diminish its availability to others. In the digital age, they are also replicable at essentially zero marginal cost: copying a software program costs a few milliseconds of computation; sharing a scientific paper costs a few bytes of bandwidth. The scarcity conditions that justify private property and market allocation simply do not hold.

Yet the dominant institutional response to non-rival knowledge goods has been to manufacture artificial scarcity through intellectual property law — patents, copyrights, trade secrets — and to treat information as if it were a rivalrous commodity. The welfare consequences of this choice are the subject of this chapter. We prove formally that IP law is welfare-reducing under standard assumptions, estimate the magnitude of value that the digital commons already generates outside the IP framework, and design governance institutions that sustain innovation without requiring artificial scarcity.

This chapter closes Part VI by establishing the information infrastructure of the cooperative-regenerative economy. The Three-Layer Coordination Stack requires massive information exchange — supply chain commitments, ecological monitoring data, governance protocols. The Planetary Ledger requires distributed data verification. The open value accounting of Chapter 18 requires transparent contribution records. All of this depends on digital commons: freely shared protocols, open-source software, open data standards, and the governance norms that sustain them. The digital commons are not peripheral to the cooperative-regenerative economy; they are its informational substrate.

---

## 33.2 Digital Commons as Non-Rival Goods

### 33.2.1 Formal Properties

**Definition 33.1 (Non-Rival Good).** A good $g$ is non-rival if its consumption by agent $i$ does not reduce its availability for agent $j$:

$$\frac{\partial u_j}{\partial x_j^g} \bigg|_{x_i^g > 0} = \frac{\partial u_j}{\partial x_j^g} \bigg|_{x_i^g = 0} \quad \forall i \neq j$$

Agent $j$'s marginal utility from good $g$ is independent of how much agent $i$ consumes. The supply of $g$ is not depleted by consumption.

**Definition 33.2 (Digital Commons Good).** A digital commons good is a non-rival good with:
- **Zero marginal reproduction cost:** $mc_{\text{copy}} \approx 0$ — creating additional copies is essentially free.
- **Partial excludability:** Access can be technically restricted (through DRM, paywalls, proprietary licenses) but such restriction is not a natural property of the good — it must be artificially imposed.
- **Network externalities:** The value of the good increases with the number of users: $u_i(x_i^g; n) > u_i(x_i^g; n-1)$ — more users make the standard more valuable, the community larger, the software better tested.

**Examples:** mathematical theorems, scientific publications, software libraries, programming languages, machine learning models, genome sequences, linguistic corpora, collaborative encyclopedias, GPS signal, open protocols (HTTP, TCP/IP), geographic data (OpenStreetMap).

### 33.2.2 Optimal Pricing of Non-Rival Goods

**Theorem 33.1 (Zero Optimal Price for Non-Rival Goods).** For a non-rival good $g$ with zero marginal reproduction cost ($mc = 0$) and positive consumption externalities (value increases with number of users), the welfare-maximizing price is zero: $p^* = 0$.

*Proof.* 

**Step 1 (Allocative efficiency).** The standard welfare theorem: set price equal to marginal cost ($p = mc$) for allocative efficiency (maximize consumer surplus without dead-weight loss). For $mc = 0$: $p^* = mc = 0$.

**Step 2 (Network externality reinforcement).** With positive consumption externalities, each additional user creates value not only for themselves but for existing users. The social marginal benefit of an additional user ($MB_{\text{social}} = MB_{\text{private}} + \text{externality value}$) exceeds the private marginal benefit. For $p > 0$: some potential users whose private $MB < p$ but social $MB > 0$ are excluded, creating a double deadweight loss (lost private surplus plus lost externality value). For $p = 0$: all users with non-negative private MB access the good, maximizing both private and social value.

**Step 3 (Formally).** Total welfare $W = \int_0^{n^*(p)} [MB_{\text{social}}(x) - mc] dx - F$ where $F$ is the fixed cost of production. For $mc = 0$: $W = \int_0^{n^*(0)} MB_{\text{social}}(x) dx - F$, maximized when $n^*(p) = n^*(0)$ — all users with positive social marginal benefit are served. Any $p > 0$ reduces $n^*$ below this optimum. $\square$

**The dynamic efficiency problem.** Theorem 33.1 establishes that zero pricing is allocatively efficient — it maximizes the use of existing knowledge goods. But it creates a dynamic efficiency problem: if price is zero, how are producers compensated for the fixed cost $F$ of creating the good in the first place? Intellectual property solves this by granting producers a temporary monopoly (patent or copyright), allowing them to charge $p > mc = 0$ to recover $F$. This creates a welfare trade-off: static inefficiency (dead-weight loss from $p > 0$) is traded for dynamic efficiency (innovation incentives).

**Proposition 33.1 (IP Law is Welfare-Reducing When Commons Alternatives Exist).** Intellectual property law generates a net welfare loss compared to public funding plus open access when:

$$\text{DWL}_{\text{IP}} > F^{\text{public}} - F^{\text{IP}} + \text{Innovation premium}_{\text{IP}}$$

where $\text{DWL}_{\text{IP}}$ is the dead-weight loss from IP pricing, $F^{\text{public}} - F^{\text{IP}}$ is the difference in production costs between public and private R\&D, and $\text{Innovation premium}_{\text{IP}}$ is any additional innovation induced by IP incentives vs. public funding.

Empirical evidence (Boldrin and Levine, 2008; Williams, 2013): the DWL from pharmaceutical patents is estimated at USD 40–90 billion/year in the US alone; the innovation premium is contested but likely much smaller; and public R\&D produces comparable innovation per dollar to private IP-protected R\&D in most fields. The condition of Proposition 33.1 is satisfied in most knowledge-intensive industries — IP law reduces net welfare.

---

## 33.3 The Hidden Value of the Digital Commons

### 33.3.1 Why GDP Misses Digital Commons Value

GDP counts market transactions. Non-market digital commons goods — open-source software, Wikipedia, open scientific data, open protocols — are produced and consumed without market transactions. They generate enormous consumer surplus (the difference between what users would pay and what they actually pay — zero) that is entirely invisible to GDP. This is the information-age amplification of the GDP-welfare gap identified in Chapter 31.

**Definition 33.3 (Consumer Surplus of Digital Commons).** The consumer surplus of a digital commons good $g$ available at price $p = 0$ is:

$$CS^g = \int_0^{\bar{n}^g} WTP_i^g \, di$$

where $WTP_i^g$ is agent $i$'s willingness to pay for access to $g$ and $\bar{n}^g$ is the total number of users.

### 33.3.2 Estimates of Digital Commons Value

**Wikipedia.** Brynjolfsson, Eggers, and Gannamaneni (2019) surveyed willingness to accept (WTA) for losing access to various digital goods. Wikipedia: median WTA = USD 3,600/year. With 1.7 billion unique monthly users globally: total consumer surplus $\approx$ USD 3,600 × 1.7 × 10⁹ / 12 months ≈ **USD 510 billion/year** — approximately 0.6\% of global GDP, generated at a cost of approximately USD 100 million/year to the Wikimedia Foundation (a leverage ratio of 5,100:1). GDP counts the USD 100 million in donations; it misses the USD 510 billion in value.

**Open-source software.** Hoffmann, Nagle, and Zhou (2024, Harvard Business School) estimated the value of open-source software to the global economy using the "cost approach" (how much would it cost to replace OSS with proprietary equivalents?). Result: the demand-side value of OSS to businesses globally is approximately **USD 8.8 trillion/year** — approximately 8.5\% of global GDP. The supply side cost to produce this OSS (largely volunteer and employer-sponsored developer time) is approximately USD 4.15 billion/year — another 2,100:1 leverage ratio of social value to production cost.

**Open scientific data.** The Human Genome Project (case study, Section 33.8) generated approximately USD 965 billion in economic value from a USD 3.8 billion investment (250:1 leverage) according to a 2013 Battelle Memorial Institute study. Open-access mandates for publicly funded research (Wellcome Trust, NIH Public Access Policy) have been estimated to generate 3-5× additional citations and broader commercial applications compared to subscription-gated equivalents.

**Aggregate estimate.** Brynjolfsson and Collis (2019) estimate that digital commons goods (freely available internet services, open-source software, Wikipedia, social media) generate consumer surplus of approximately **USD 30–50 trillion/year** globally — roughly 35–58\% of global GDP — that is entirely invisible to conventional GDP accounting. The cooperative-regenerative economy's information layer is already generating value comparable to the market economy; it simply is not counted.

---

## 33.4 Governance of Digital Commons

### 33.4.1 Ostrom Principles Applied

Chapter 14 proved that Ostrom's eight design principles are sufficient for commons governance stability. We apply them to three major digital commons: Wikipedia, the Linux kernel, and OpenStreetMap.

**Wikipedia (English language, ~6.7 million articles, ~270,000 active editors):**

| DP | Implementation | Score |
|:---|:---|:---|
| DP1 (Boundaries) | Registered editors; notability criteria define resource scope | 1.5/2 |
| DP2 (Congruence) | Policies vary by article type (biography, science, politics) | 2/2 |
| DP3 (Collective choice) | RfC, ArbCom, community votes on policy | 2/2 |
| DP4 (Monitoring) | Automated bots + volunteer patrollers | 2/2 |
| DP5 (Graduated sanctions) | Warning → talk page notice → block → indefinite ban | 2/2 |
| DP6 (Conflict resolution) | Talk pages → Mediation → ArbCom, but slow and opaque | 1/2 |
| DP7 (External recognition) | 501(c)(3) legal status; Wikimedia Foundation | 2/2 |
| DP8 (Nested enterprises) | Wikipedia → Wikimedia Foundation → Wikimedia movement | 1.5/2 |
| **Total** | | **14/16** |

Wikipedia's score of 14/16 predicts high stability — consistent with its sustained operation since 2001. Its weakness at DP6 (conflict resolution: ArbCom decisions are perceived as slow and opaque by many editors) is a genuine institutional vulnerability correlated with the well-documented decline in editor retention and the resulting systemic bias toward topics covered by remaining editors (predominantly Western, male, and English-speaking).

**Linux kernel (~30 million lines, ~4,200 active contributors/year):**

Previously scored 14/16 in Chapter 11. The Linux kernel's BDFL (Linus Torvalds) structure partially violates DP3 (collective choice) but provides structural coherence that the Wikipedia model cannot achieve — centralized technical direction without centralized editorial control.

**OpenStreetMap (~10 million registered contributors, 500 million+ map features):**

Previously scored 10/16 in Chapter 11. OSM's weakest area remains DP5 and DP6 — vandalism sanctions are weak and conflict resolution is ad hoc. This is partly structural: geographic data can be degraded by subtle inaccuracies that are difficult to detect without local knowledge, creating monitoring (DP4) challenges that text-based commons do not face in the same way.

### 33.4.2 The Non-Rival Contribution Game

**Definition 33.4 (Digital Commons Contribution Game).** The contribution game for a digital commons good $G$ is a public goods game with characteristic function:

$$v(S) = B\left(\sum_{i \in S} c_i\right) - \sum_{i \in S} \text{cost}(c_i)$$

where $c_i$ is agent $i$'s contribution (code commits, edits, data submissions), $B(\cdot)$ is the collective benefit function (concave — diminishing returns to additional contributions), and cost$(c_i)$ is the individual cost of contribution.

**Key difference from rivalrous commons.** In a rivalrous commons, contributions are subtracted (the more others contribute, the less remains for me). In a digital commons, contributions are additive — more contributors make the good better, not scarce. The free-rider problem still exists (I can use the good without contributing), but the mechanism is different: free-riding is a welfare problem (under-provision), not a sustainability problem (the good does not deplete).

**Proposition 33.2 (Digital Commons Under-Provision).** The Nash equilibrium of the digital commons contribution game achieves contribution level $\sum c_i^{\text{NE}} < \sum c_i^*$ (social optimum), by the standard public goods under-provision result. The under-provision gap is:

$$\sum c_i^* - \sum c_i^{\text{NE}} = \frac{n-1}{n} \cdot \left[\sum c_i^*\right] \cdot \frac{B'(\sum c^*)}{B'(\sum c^*) + \text{cost}''(c^*)}$$

*Proof.* Standard Bergstrom-Blume-Varian public goods result: the Nash equilibrium provides $1/n$ of the socially optimal contribution in the symmetric case, because each agent internalizes only $1/n$ of the social benefit of their contribution (the rest accrues to the other $n-1$ contributors). $\square$

**Corollary 33.1 (Shapley Allocation as Contribution Incentive).** The Shapley value of the digital commons contribution game allocates to each contributor their average marginal contribution to the collective good. A digital commons platform that implements Shapley value allocation — distributing revenue or recognition in proportion to average marginal contribution — internalizes the externality gap and achieves closer to socially optimal contribution levels.

This is the theoretical foundation for OVA (Open Value Accounting, Chapter 18) applied to digital commons: by measuring and rewarding marginal contributions through a Shapley-approximating mechanism, cooperative digital platforms can overcome the public goods under-provision that pure volunteer systems face.

---

## 33.5 Data as a Commons

### 33.5.1 Data: Between Non-Rival and Rival

Data occupies an interesting intermediate position between pure public goods and private goods.

**Data is non-rival in processing:** The same dataset can be analyzed by multiple researchers simultaneously without any reduction in data quality. One hospital's access to a patient dataset does not reduce another's access to the same dataset.

**Data is partially rival in privacy and market value.** Personal data, once disclosed to third parties, cannot be "un-disclosed" — the information asymmetry is irreversibly changed. In markets, exclusive access to proprietary data confers competitive advantage: if Amazon shares its sales data with suppliers, suppliers gain insights that reduce Amazon's information advantage. Data is a commons in processing but a private good in strategic disclosure.

**Definition 33.5 (Data Commons).** A data commons is a dataset governed under the Fifth Magisterium of the Commons [C:Ch.14], with:
- **Bounded membership:** Defined contributors and authorized users.
- **Collective governance:** Decisions about data access, use, and sharing governed by the contributing community.
- **Stewardship obligation:** Data quality, privacy, and long-run integrity maintained for future users.
- **Non-commodification norm:** Raw personal data not sold to third parties; insights derived from the data may be commercialized under revenue-sharing agreements.

**The data cooperative.** A data cooperative is an institutional form that implements the data commons: members contribute their personal data (health records, energy consumption, mobility patterns, financial transactions), the cooperative governs access and use, and members share in the value generated from the data pool. The data trust is a related form in which a trustee holds data on behalf of contributors under fiduciary obligations.

### 33.5.2 GDPR as Institutional Experiment

The European Union's General Data Protection Regulation (GDPR, 2018) represents the most ambitious attempt to govern personal data as a commons rather than as private property of the platforms that collect it. Its formal content:

- **Article 15-22:** Data subject rights — access, rectification, erasure, portability, and objection. These rights convert data from a platform asset to a contributor asset.
- **Article 25:** Data protection by design — privacy must be embedded in systems from inception, not added as an afterthought. This is the data governance equivalent of the DfD (Design for Disassembly) condition of Chapter 21.
- **Article 35:** Data Protection Impact Assessments for high-risk processing — the data governance equivalent of ecological impact assessment.

**Formal assessment against Ostrom principles.** The GDPR partially implements the commons governance structure:

| Ostrom DP | GDPR Implementation | Compliance |
|:---|:---|:---|
| DP1 (Defined boundaries) | Jurisdiction: EU residents' data | Partial (extraterritorial gaps) |
| DP2 (Congruence) | Sector-specific rules (health, financial) | Moderate |
| DP3 (Collective choice) | Data subject rights (individual) but no collective governance | **Absent** |
| DP4 (Monitoring) | DPO requirements; supervisory authorities | Moderate |
| DP5 (Graduated sanctions) | Tiered fines (2\% / 4\% of global turnover) | Strong |
| DP6 (Conflict resolution) | National supervisory authority complaints | Moderate |
| DP7 (External recognition) | EU law; Schrems II limits US transfers | Strong (EU only) |
| DP8 (Nested enterprises) | National DPAs → EDPB | Moderate |

GDPR's most significant gap: **DP3 (collective choice)**. GDPR gives individual data subjects rights but creates no collective governance mechanism. Data contributors cannot collectively decide how their pooled data is used — only individually consent or object to uses of their own data. The data cooperative fills this gap: it provides the collective governance layer that GDPR's individual rights framework cannot provide.

---

## 33.6 AI and the Digital Commons

### 33.6.1 The Dual Relationship

AI development has a dual relationship with the digital commons: it depends on the commons (training data, open-source frameworks, scientific publications) and it threatens the commons (by extracting value without reciprocating, concentrating AI capabilities in private hands, and generating synthetic content that degrades the commons).

**AI as commons beneficiary.** Large language models (GPT-4, Claude, Llama) are trained on datasets including Wikipedia (≈USD 510 billion in social value), Common Crawl (open web data), GitHub (open-source code), arXiv (open scientific papers), and books1/books2 (copyrighted works under fair use claims). The AI training datasets are themselves partially composed of commons goods — goods whose open accessibility was a deliberate design choice by commons communities.

**Formal statement.** The value of AI system $A$ trained on digital commons dataset $\mathcal{D}$ is:

$$V(A) = f\left(\mathcal{D}, \theta\right)$$

where $\theta$ is the model architecture and training investment. The commons contribution $\mathcal{D}$ is an input to AI value creation. The Shapley value of the commons contribution $\mathcal{D}$ to the AI's value:

$$\phi_{\mathcal{D}}(V) = \mathbb{E}_{S \subseteq \{A, \mathcal{D}\}}[V(S \cup \{\mathcal{D}\}) - V(S)]$$

This is the additional value the AI system generates because of the commons training data — which rightfully belongs to the commons community under Shapley allocation logic.

**AI as commons threat.** Three mechanisms through which AI can undermine digital commons:

1. **Content extraction without reciprocation:** AI systems train on commons content but do not contribute back — they are free riders at industrial scale. The Wikipedia article that a contributor spent hours writing generates training signal worth far more (in AI capability terms) than the contributor receives in return.

2. **Synthetic content degradation:** AI systems generate text, code, and data that flood back into commons pools — potentially degrading quality, spreading misinformation, or homogenizing content in ways that reduce the diversity value of the commons.

3. **Concentration of capability:** If AI development is concentrated in a few private firms with proprietary models, the AI capability commons (the collective ability of society to use AI) becomes a private good — accessible only to those who can pay, or on terms set by private firms with misaligned incentives.

### 33.6.2 Governance Conditions for Beneficial AI

**Theorem 33.2 (Commons-Compatible AI Conditions).** AI development benefits rather than undermines the digital commons if and only if:

1. **Reciprocity:** AI systems that extract value from commons training data contribute back to the commons — either through data donation, model open-sourcing, or financial contribution to commons maintenance ($\phi_{\mathcal{D}}$ returned to the commons community).

2. **Open-source models:** The model architecture, weights, and training methodology are published under open licenses — making AI capability itself a digital commons good rather than a proprietary asset.

3. **Content provenance:** AI-generated content is labeled as such and governed by quality standards that prevent degradation of commons quality — analogous to the quality standards that Wikipedia's DP2 (congruence) and DP4 (monitoring) enforce for human contributions.

4. **Distributed governance:** AI development governance involves multi-stakeholder democratic processes [C:Ch.13, Cosmo-Local model] rather than unilateral decisions by private AI labs — analogous to IETF's rough consensus model for internet protocol governance.

*Proof.* Conditions 1–4 together ensure: (1) the Shapley value is returned to the commons (no net extraction from the commons); (2) AI capability is itself a commons good (no privatization of collective intelligence); (3) commons quality is maintained (no degradation from synthetic content); (4) commons governance is maintained (no capture by private AI interests). Any of the four conditions failing creates a pathway for AI to undermine the commons. $\square$

---

## 33.7 Mathematical Model: The Public Goods Contribution Game

**Formal specification.** A digital commons with $n$ potential contributors, each choosing contribution level $c_i \geq 0$ at private cost $\text{cost}(c_i) = \frac{1}{2}\kappa c_i^2$. The collective benefit:

$$B\left(\sum_i c_i\right) = \beta \ln\left(1 + \sum_i c_i\right)$$

(logarithmic — diminishing returns to aggregate contribution). Each agent receives equal share of the benefit: $u_i = \beta \ln\left(1 + \sum_j c_j\right) - \frac{1}{2}\kappa c_i^2$.

**Nash equilibrium.** FOC for agent $i$: $\beta/(1 + \sum_j c_j) = \kappa c_i$. In symmetric equilibrium: $c^{\text{NE}} = \beta/(\kappa(1 + nc^{\text{NE}}))$. Solving: $c^{\text{NE}}(1 + nc^{\text{NE}}) = \beta/\kappa$, giving $nc^{(NE)2} + c^{\text{NE}} - \beta/\kappa = 0$. For large $n$: $c^{\text{NE}} \approx \sqrt{\beta/(\kappa n)}$ — each individual contributes approximately $1/\sqrt{n}$ as the community grows.

**Social optimum.** Maximize $\sum_i u_i = n\beta\ln(1+\sum c_i) - \frac{\kappa}{2}\sum c_i^2$. FOC: $n\beta/(1+\sum c_i^*) = \kappa c_i^*$. In symmetric optimum: $c^* = \sqrt{n\beta/\kappa}$ — each agent contributes $\sqrt{n}$ times their Nash contribution.

**Under-provision gap:**

$$\frac{c^*}{c^{\text{NE}}} = n \quad \text{(factor of } n\text{ under-provision)}$$

The social optimum requires $n$ times the Nash equilibrium contribution. For a community of $n = 10,000$ contributors: each individual contributes 10,000× less than the social optimum in the Nash equilibrium. This is the formal quantification of open-source under-provision — the problem that platform economics, Shapley allocation, and OVA mechanisms aim to address.

**The Shapley allocation fix.** Under Shapley value allocation of the collective benefit: $\phi_i = B(\sum c_j)/n + [\phi_i^{\text{surplus}}]$. With each contributor receiving their full average marginal contribution, the effective return to contribution rises from $\beta/n$ (proportional equal sharing) toward $\beta$ (full social benefit). This closes the under-provision gap: contributors internalize a larger fraction of the social benefit of their contribution, approaching the social optimum.

---

## 33.8 Worked Example: The Social Value of the Python Ecosystem

### 33.8.1 The Python Ecosystem

Python is the world's most widely used programming language (as of 2024 TIOBE index), with an ecosystem of approximately 500,000 open-source packages hosted on PyPI (the Python Package Index). The Python language itself (CPython implementation) is open-source, maintained by the Python Software Foundation with approximately 1,500 active contributors. Key packages — NumPy, pandas, SciPy, scikit-learn, TensorFlow, PyTorch — are used in scientific research, data analysis, machine learning, and web development across every major industry.

### 33.8.2 Value Estimation

**Replacement cost approach.** If the Python ecosystem did not exist, organizations would need to either: (i) build proprietary equivalents; or (ii) use inferior alternatives. Nagle et al. (2024, extending the methodology of Hoffmann et al.):

- Python language (CPython): replacement cost approximately USD 2.5 billion (estimated from comparable proprietary language development costs).
- Top 20 packages by usage: NumPy, pandas, requests, etc. — combined replacement cost approximately USD 15.2 billion.
- Full PyPI ecosystem (500,000 packages): extrapolating from top-20 sample, approximately USD 48 billion total replacement cost.

**Consumer surplus approach.** Brynjolfsson et al. WTA methodology: surveyed Python users' willingness to accept compensation for losing access. Median WTA: USD 1,800/year (professional users); USD 400/year (casual users). Distribution: approximately 25 million professional Python users globally, 75 million casual users.

$$CS^{\text{Python}} = 25\text{M} \times \$1{,}800 + 75\text{M} \times \$400 = \$45\text{B} + \$30\text{B} = \mathbf{\$75 \text{ billion/year}}$$

**Annual production cost.** Python Software Foundation budget: approximately USD 6 million/year. Volunteer developer time (estimated): 250,000 contributor-hours/year at USD 150/hour = USD 37.5 million. Corporate sponsorships: approximately USD 30 million/year. Total production cost: approximately **USD 73 million/year**.

**Leverage ratio:** $75{,}000 \text{M} / 73\text{M} \approx 1{,}027$:1 — the Python ecosystem generates approximately USD 1,000 in social value for every USD 1 invested in its production.

### 33.8.3 The Optimal Public Subsidy

**Proposition 33.3 (Optimal Public Subsidy for Digital Commons).** The optimal public subsidy $S^*$ for a digital commons with social value $CS^g$ and private production $F^{\text{private}}$ satisfies:

$$S^* = CS^g - F^{\text{private}}$$

when the commons provides a pure public good and private investment is zero in the absence of subsidy (no excludability, so private returns = 0). For partial excludability (some private return possible through consulting, support, training):

$$S^* = CS^g - F^{\text{private}} - \pi^{\text{private}}$$

where $\pi^{\text{private}}$ is the private return from excludable complementary services.

For Python: $CS^g = $ USD 75 billion/year, $F^{\text{private}} = $ USD 73 million, $\pi^{\text{private}} = $ USD 500 million (Python consulting, training, commercial distributions):

$$S^* = 75{,}000 - 73 - 500 = \mathbf{\$74.4 \text{ billion/year}}$$

**Actual public investment in Python:** PSF government grants $\approx$ USD 2 million/year; NSF funding for scientific Python packages $\approx$ USD 50 million/year; EU Horizon funding for open-source $\approx$ USD 20 million/year. Total public investment $\approx$ **USD 72 million/year** — approximately **0.1\% of the optimal subsidy level**.

The Python ecosystem is massively undersubsidized by public investment — an inevitable consequence of GDP accounting that misses the USD 75 billion in annual consumer surplus it generates. If public R\&D investment were allocated to maximize social welfare rather than measured GDP, Python and similar digital commons goods would receive orders-of-magnitude more public funding.

---

## 33.9 Case Study: The Human Genome Project as Digital Commons

### 33.9.1 The Open-Access Decision

In 1998, the Human Genome Project (HGP) — a publicly funded international consortium — faced a critical decision: whether to patent and restrict access to the genome sequences being produced, or to release them into the public domain through the "Bermuda Principles" (release all data publicly within 24 hours of production).

Simultaneously, Celera Genomics (a private company) entered the race to sequence the human genome, with the explicit intention of patenting portions and selling database access. The competition between HGP (open commons) and Celera (private property) became the most consequential intellectual property experiment in the history of science.

### 33.9.2 Economic Consequences of the Open-Access Choice

**The Bermuda Principles (HGP).** All sequence data entered the public domain immediately. Researchers worldwide could access the genome freely, accelerating research across all fields without licensing costs, patent barriers, or access restrictions.

**Economic value generated.** A 2013 study by the Battelle Memorial Institute estimated the economic impact of HGP:
- Total public investment in HGP: USD 3.8 billion (1988–2003)
- Economic output attributable to HGP: USD 796 billion (2003–2013)
- ROI: 209:1 in 10 years
- Jobs created: 310,000
- Tax revenue generated: USD 244 billion

The Williams (2013, Harvard) natural experiment: comparing gene segments that Celera did patent (restricting access) to those that HGP released (open access). Celera's patented segments received 30\% fewer subsequent scientific citations and 20\% less follow-on innovation than comparable HGP segments. The Williams study provides the cleanest causal evidence in the literature that IP restrictions reduce innovation spillovers — exactly the mechanism Proposition 33.1 predicts.

**The counterfactual.** If HGP had followed Celera's proprietary model: estimated 30\% fewer follow-on innovations, delayed cancer treatments by 3–8 years (based on citation delays in Williams 2013), and reduced economic value generated by approximately USD 240 billion over 10 years.

**Formal validation.** Applying the optimal subsidy formula (Proposition 33.3) to the human genome: $CS^g \approx $ USD 796 billion, $F^{\text{private}} = $ USD 3.8 billion, $\pi^{\text{private}} \approx 0$ (sequence data non-excludable under Bermuda Principles). $S^* = 796 - 3.8 - 0 \approx $ USD 792 billion — the correct public subsidy was the full consumer surplus, appropriately achieved through public funding plus open access. Actual public investment of USD 3.8 billion was USD 792 billion less than optimal, but the open access decision captured the full social value by preventing private appropriation.

**The institutional lesson.** The HGP case demonstrates that the optimal governance structure for knowledge goods is: public funding for production (solving the dynamic efficiency problem) plus open access for distribution (solving the allocative efficiency problem). Intellectual property law does neither: it creates private incentives for production (solving dynamic efficiency imperfectly) while imposing allocative inefficiency (restricting distribution). When the public sector can fund production directly, public funding plus open access dominates IP law on both static and dynamic efficiency grounds.

---

## Chapter Summary

This chapter has established the formal economics of digital commons, proving that non-rival goods require fundamentally different governance than rival goods, estimating the enormous hidden value of existing digital commons, and designing governance structures that sustain innovation and equity simultaneously.

The zero optimal price theorem (Theorem 33.1) proves that allocative efficiency requires zero pricing for non-rival goods with zero marginal reproduction costs. IP law is welfare-reducing when commons alternatives exist (Proposition 33.1) — confirmed empirically by the Williams (2013) human genome natural experiment. The social value of digital commons — approximately USD 30–50 trillion/year globally — is systematically invisible to GDP, amplifying the GDP-welfare gap analyzed in Chapter 31.

The Ostrom analysis of Wikipedia (14/16), Linux (14/16), and OpenStreetMap (10/16) confirms that digital commons stability is governed by the same principles as physical commons — with DP6 (conflict resolution) and DP3 (collective choice) as the most common failure points. The non-rival contribution game (Proposition 33.2) quantifies the under-provision gap: $n$-fold under-contribution relative to social optimum in symmetric Nash equilibrium. Shapley value allocation (Corollary 33.1) partially corrects this by internalizing the contribution externality.

Data occupies a middle position — non-rival in processing but partially rival in privacy and market value. The data cooperative fills GDPR's DP3 gap by providing collective governance over data use. The human genome project and the GDPR analysis together illustrate the institutional design space: public funding plus open access dominates IP law; individual data rights without collective governance is insufficient.

AI development is a dual-edged phenomenon for digital commons: it extracts massive value from commons training data while potentially concentrating AI capabilities in private hands and degrading commons quality with synthetic content. Theorem 33.2 specifies four conditions for commons-compatible AI: reciprocity, open-source models, content provenance, and distributed governance. The Python ecosystem worked example quantifies the under-investment problem: USD 72 million in actual public investment vs. USD 74 billion optimal — a 1,000:1 gap that reflects the systematic under-counting of digital commons value in GDP-based resource allocation.

Part VI is now complete. Five chapters have assembled the unified cooperative-regenerative framework: Chapter 29 proved existence of the CRE and the Cooperative Stewardship Theorem; Chapter 30 demonstrated structural resilience; Chapter 31 proved post-growth prosperity; Chapter 32 derived the inequality-reducing mechanisms of cooperative institutions; and this chapter established the information infrastructure of the cooperative-regenerative economy. Part VII grounds all of this theory in empirical case studies across six sectors and geographic contexts.

---

## Exercises

**33.1** Prove that a knowledge good with zero marginal reproduction cost should be priced at zero for allocative efficiency (Theorem 33.1).
(a) Set up the welfare maximization problem for a non-rival good with $mc = 0$ and $n$ potential users with heterogeneous valuations $v_i$.
(b) Show that any price $p > 0$ excludes users with $v_i < p$, creating dead-weight loss $\int_0^p [n \cdot \bar{v}(x)]dx > 0$.
(c) Prove that at $p = 0$, all users with $v_i > 0$ access the good and total welfare is maximized.
(d) **Dynamic efficiency problem:** If the producer charges $p = 0$, how are the fixed production costs $F$ recovered? Identify three alternative funding mechanisms other than IP law and assess each on both static and dynamic efficiency grounds.

**33.2** Apply Ostrom's eight design principles to Wikipedia (Section 33.4.1). For each principle:
(a) Identify the specific Wikipedia institution or policy that implements it.
(b) Assess the strength of implementation on a 0–2 scale with justification.
(c) For the principle with the lowest score: propose a specific institutional reform that would raise it by 1 point. Would the reform be politically feasible given Wikipedia's existing community governance norms?

**33.3** The Python ecosystem social value:
(a) The Python Software Foundation has a budget of USD 6 million/year. Using Proposition 33.3, compute the funding gap between actual and optimal public investment.
(b) If governments allocated 1\% of the optimal subsidy to Python maintenance, what governance changes would be necessary to ensure the funds are spent efficiently? Use the Cosmo-Local model [C:Ch.13] to specify the governance structure.
(c) Several companies (Google, Microsoft, Meta) sponsor Python development significantly. Compute the Shapley value of each company's contribution to Python's social value, using their documented contribution hours and financial support as inputs. Is the current sponsorship level consistent with their Shapley value share?

**★ 33.4** Prove formally that the social optimum for digital commons under-provision requires $n$-times the Nash equilibrium contribution.

(a) Set up the public goods game with $n$ agents, logarithmic benefit function $B(\sum c_i) = \beta\ln(1+\sum c_i)$, and quadratic costs $\text{cost}(c_i) = \kappa c_i^2/2$.
(b) Derive the symmetric Nash equilibrium contribution $c^{\text{NE}}$.
(c) Derive the symmetric social optimum contribution $c^*$.
(d) Compute the ratio $c^*/c^{\text{NE}}$ and show it equals $n$.
(e) How does this ratio change if agents are heterogeneous (different $\kappa_i$ and different private benefit shares $b_i < 1$)? Under what conditions does heterogeneity improve upon the symmetric Nash equilibrium?

**★ 33.5** Prove Theorem 33.2 (conditions for commons-compatible AI) and evaluate current leading AI systems against it.

(a) Formalize "reciprocity" (Condition 1): define the Shapley value $\phi_{\mathcal{D}}$ of the training data commons in the AI value creation game, and specify the reciprocity condition as $\phi_{\mathcal{D}} \leq \text{Value returned to commons}$.
(b) For a leading AI system of your choice (GPT-4, Claude, Llama, Gemini): estimate $\phi_{\mathcal{D}}$ using available information about training data composition. What fraction of the AI system's value derives from digital commons training data vs. proprietary investment?
(c) Evaluate each of the four conditions (reciprocity, open-source, content provenance, distributed governance) for your chosen AI system. Score each condition 0–2 (as with the Ostrom principles).
(d) Based on your assessment: is your chosen AI system commons-compatible? What specific changes would bring it into compliance with all four conditions?

**★★ 33.6** Design a data cooperative for a 10,000-patient medical dataset; specify the full governance structure, contribution valuation method, and revenue distribution mechanism.

**Dataset description:** 10,000 patients contributed electronic health records (EHR), genomic sequences, wearable device data, and consent to participate in research studies. The cooperative's goal: govern data sharing with pharmaceutical companies, academic researchers, and health systems in a way that benefits patients collectively.

**(a) Governance structure (Ostrom principles):**

For each of the eight Ostrom design principles, specify the concrete institutional mechanism that implements it:
- DP1: Who are the members? What data is in scope?
- DP2: Are rules uniform or adapted to subgroups (e.g., genetic data vs. EHR)?
- DP3: How do members collectively decide which research access requests to approve?
- DP4: How is data use monitored? Who audits access logs?
- DP5: What sanctions apply to unauthorized data use?
- DP6: How do members dispute data use decisions?
- DP7: What legal structure protects the cooperative's governance autonomy?
- DP8: How does the cooperative nest within national health data governance?

**(b) Contribution valuation (Shapley value):**

Model the dataset's value as a cooperative game: the value of any subset $S$ of patients' data is $v(S) = f(|S|, \text{richness}(S))$ where richness includes diversity of conditions, genomic variation, and longitudinal completeness. Specify:
- The characteristic function $v(S)$.
- The Shapley value $\phi_i(v)$ for each of the five most common patient profiles.
- The OVA implementation: how are Shapley contributions calculated and updated as new data is contributed?

**(c) Revenue distribution mechanism:**

The cooperative licenses data access to three pharmaceutical companies, generating EUR 8 million/year. Specify:
- The allocation between: (i) patient dividends (in proportion to Shapley value); (ii) cooperative maintenance and governance costs; (iii) reinvestment in data infrastructure.
- The minimum patient dividend (EUR/patient/year) if 70\% of revenue is distributed proportionally to Shapley value.
- The governance mechanism for deciding how to use the 30\% retained for cooperative purposes.

**(d) Prove incentive-compatibility:**

Show that under the specified governance and revenue distribution, it is individually rational for each patient to: (i) contribute their data; (ii) maintain data quality; (iii) participate in governance. Use the Folk Theorem framework [C:Ch.7] and the Ostrom incentive conditions [C:Ch.14].

---

*Part VII opens with Chapter 34. The formal framework of Parts I–VI is now complete: the game theory, network science, ecological embedding, monetary alternatives, unified model, stability analysis, post-growth economics, inequality dynamics, and digital commons governance all assembled. Part VII tests this framework against reality in six real-world applications — cooperative enterprises, P2P platforms, regenerative agriculture, complementary currencies, universal basic services, and data cooperatives — identifying where the theory is confirmed and where it requires extension.*
