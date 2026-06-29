# Chapter 41c: Democratic Economic Governance
## Planning Beyond GDP, Citizens' Deliberation, and Post-Growth Accountability

> *"GDP measures everything except that which makes life worth living."*
> — Robert F. Kennedy, University of Kansas (1968)

> *"The 'encasement' of economic activity in technocratic institutions and legal rules has resulted in
> the subordination of popular, democratic priorities to market imperatives that protect
> economic decision-making from public demands."*
> — *Roadmap for Eradicating Poverty Beyond Growth*, UN Special Rapporteur (2026)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Distinguish formally between a technocratic governance equilibrium — in which economic
   decision-making is insulated from democratic deliberation through delegation to independent
   institutions — and a democratic governance equilibrium, in which economic objectives are
   defined, monitored, and revised through participatory processes, and prove that the latter
   achieves strictly higher IPI under a condition of systematic objective misspecification.
2. Model national democratic planning as constrained multi-objective optimisation over the
   Multidimensional Provisioning Dashboard (MPD), connecting it to the policy optimisation
   problem of Chapter 41 and identifying the formal role of participatory deliberation in
   determining the welfare weight vector.
3. Specify the Sustainable and Inclusive Wellbeing (SIW) accounting framework as a formal
   extension of the SFC-N system of Chapter 18, integrating economic flows, ecological stocks,
   and distributional outcomes into a single internally consistent triple of accounts, and prove
   that SIW and GDP diverge in direction under the conditions of extractive growth.
4. Model citizens' assemblies as a Bayesian preference aggregation mechanism operating under
   a justificatory deliberation constraint, and derive the conditions under which the assembly
   process converges to a Condorcet winner in the MPD weight space.
5. Formalise participatory social auditing as an extension of Ostrom's DP4 monitoring
   principle to multi-level economies, prove that it reduces resource misallocation relative to
   internal audit, and demonstrate how the Democratic Accountability Loop connects planning,
   implementation, auditing, and sanction into a stable iterative governance cycle.

---

## 41c.1 The Governance Gap

The cooperative-regenerative framework, as it stands at the end of Chapter 41, is formally
complete in one sense: we have proved that cooperative institutions produce higher welfare
(Theorem 29.2), that the transition is analytically tractable (Theorem 40.2), and that a specific
policy portfolio accelerates it cost-effectively (Chapter 41). But there is a prior question that
none of this resolves — a question so familiar it risks appearing obvious: who decides?

Who decides what the IPI target is? Who decides which ten dimensions of the Multidimensional
Provisioning Dashboard matter, and in what proportions? Who decides when a Planetary Boundary
violation is severe enough to trigger the Layer 3 biophysical planning mechanisms of Chapter 29?
Who holds cooperative enterprises accountable when their OVA allocations diverge from
Shapley-consistent norms? Who monitors the natural capital accounts that Chapter 18 formalised,
and acts on them when stewardship conditions are violated?

These are not technical questions with uniquely correct answers derivable from economic
analysis. They are inherently political questions — questions about whose values prevail when
values conflict, whose knowledge is recognised when knowledge is contested, and whose interests
are protected when interests diverge. The history of economic governance is, in significant part,
the history of delegating these questions to institutions whose claim to authority rests on technical
expertise rather than democratic mandate. Central banks set interest rates; fiscal councils assess
budget trajectories; environmental agencies adjudicate pollution permits; credit rating agencies
determine the cost of sovereign borrowing. Each of these institutions makes deeply political
decisions — about whose welfare is prioritised, over what time horizon, against what ecological
constraints — while presenting them as technical conclusions that democratic deliberation cannot
improve.

The problem is not that technical expertise is worthless. It is that the dominant governance
architecture has fused the appropriate claim of expertise — that technical analysis can inform
decisions — with the inappropriate claim that technical analysis can replace them. GDP growth
has been the organising logic of this architecture for eighty years. It is not merely an economic
indicator. It is a political choice, embedded in accounting conventions, model structures, and
institutional mandates, that systematically privileges output over distribution, the present over
the future, and the measurable over the real.

This chapter provides the cooperative-regenerative framework's answer to the governance
question. It does not replace technical analysis with populism; it specifies the institutional
architecture through which democratic deliberation and technical expertise can combine
productively — each doing what the other cannot. Democratic deliberation defines objectives
and weights values; technical analysis models the consequences of different policy choices for
those democratically determined objectives. The Sustainable and Inclusive Wellbeing framework
provides the accounting system. The citizens' assembly provides the deliberative mechanism.
Participatory social auditing provides the accountability mechanism. Together, they constitute
what we call the Democratic Accountability Loop — the steering architecture without which the
cooperative-regenerative economy remains a framework without a helm.

---

## 41c.2 Technocracy and Democracy: A Formal Comparison

The choice between technocratic and democratic governance of economic objectives is not a
rhetorical debate about political philosophy. It has precise formal content, and the formal
analysis reveals welfare consequences that the political philosophy literature has rarely made
explicit.

**Definition 41c.1 (Technocratic Governance Equilibrium).** A Technocratic Governance
Equilibrium (TGE) is an economic governance arrangement $\mathcal{G}^T = (\mathcal{O},
\mathcal{I}, \mathcal{A}, \mathcal{D})$ where:

- $\mathcal{O}$: a fixed set of policy objectives, constitutionally or legislatively mandated
  and deliberately insulated from routine democratic revision (e.g., price stability, fiscal
  balance, GDP growth);
- $\mathcal{I}$: independent institutions with delegated authority to pursue $\mathcal{O}$,
  whose staff are appointed, not elected, and whose operational independence is legally protected;
- $\mathcal{A}$: a narrow accountability mechanism — typically parliamentary ratification of
  appointments and periodic legislative review — that does not include continuous democratic
  oversight of day-to-day decisions;
- $\mathcal{D}$: a deliberative space restricted to credentialed technical experts, excluding
  lay deliberation from high-stakes economic decisions.

**Definition 41c.2 (Democratic Governance Equilibrium).** A Democratic Governance
Equilibrium (DGE) is an economic governance arrangement $\mathcal{G}^D = (\mathcal{O}^D,
\mathcal{I}^D, \mathcal{A}^D, \mathcal{D}^D)$ where:

- $\mathcal{O}^D$: objectives defined through inclusive democratic deliberation —
  including the SIW framework's three dimensions (current wellbeing, sustainability,
  inclusion) — and subject to revision through the same process at each planning cycle;
- $\mathcal{I}^D$: institutions accountable to democratic bodies through continuous
  reporting requirements, not merely periodic appointment review;
- $\mathcal{A}^D$: accountability that is continuous and participatory — social auditing,
  citizens' assemblies with binding follow-through, community monitoring;
- $\mathcal{D}^D$: a deliberative space that includes citizens' assemblies, participatory
  budgeting processes, and trade union social dialogue as co-equal partners with technical
  expert bodies.

These are not idealisations. The TGE describes the governance architecture that has prevailed in
OECD countries since the 1990s: independent central banks, fiscal councils, environmental
agencies, and international financial institutions with mandates set by treaty rather than
democratic revision. The DGE describes a genuine institutional alternative, with existing
instantiations in New Zealand's Wellbeing Budget, Wales's Well-being of Future Generations
Act (Chapter 41), Ireland's Citizens' Assembly, and France's Convention Citoyenne pour le
Climat. Neither is pure; all real governance systems mix elements of both. But the formal
distinction is analytically useful because it identifies the welfare consequences of moving along
the spectrum from TGE toward DGE.

**Theorem 41c.1 (DGE Welfare Advantage).** Let $W^{\text{true}}$ be the true social welfare
function — the function that correctly captures the values of the population, including
distributional equity, ecological sustainability, and long-run flourishing. Let $W^{\text{TGE}}$
be the welfare function actually maximised by the TGE, which we assume to be a weighted sum
of GDP, price stability, and fiscal balance. Define the objective misspecification error as:

$$\delta_{\mathcal{O}} = \sup_{\mathbf{u} \in U} \bigl[ W^{\text{true}}(\mathbf{u}) - W^{\text{TGE}}(\mathbf{u}^*_{\text{TGE}}) \bigr] > 0$$

where $\mathbf{u}^*_{\text{TGE}}$ is the TGE's optimal policy. Then, provided the DGE's
deliberative mechanism reduces $\delta_{\mathcal{O}}$ toward zero as the number of deliberative
cycles increases:

$$\text{IPI}(\mathcal{G}^D) > \text{IPI}(\mathcal{G}^T)$$

for all planning horizons $T$ beyond a threshold $\bar{T}$ determined by the speed of the DGE's
deliberative convergence.

*Proof.* $\text{IPI}(\mathcal{G}^T) = \sum_t \beta^t W^{\text{true}}(\mathbf{u}^*_{\text{TGE}}(t))$.
Since the TGE maximises $W^{\text{TGE}} \neq W^{\text{true}}$, there exists at each period $t$
a policy $\tilde{\mathbf{u}}(t)$ with $W^{\text{true}}(\tilde{\mathbf{u}}(t)) > W^{\text{true}}(\mathbf{u}^*_{\text{TGE}}(t))$ — specifically, the policy that directly maximises $W^{\text{true}}$
subject to the ecological and fiscal constraints. The DGE selects this policy through its
deliberative update of $\mathcal{O}^D$. The IPI gap $\text{IPI}(\mathcal{G}^D) - \text{IPI}(\mathcal{G}^T) = \sum_t \beta^t [W^{\text{true}}(\mathbf{u}^*_{\text{DGE}}) - W^{\text{true}}(\mathbf{u}^*_{\text{TGE}})] > 0$
for $T > \bar{T}$, since the DGE's welfare gain compounds as the objective misspecification
error is progressively corrected. $\square$

**The objective misspecification error has three components.** In the current TGE:

*(i) Distributional omission:* GDP does not measure how income is distributed; the TGE
optimising GDP growth can achieve its objective while inequality rises without bound.
Chapter 32 showed that $r > g$ is sufficient for this — the TGE's own growth imperative
produces the distributional dynamics that violate the true welfare function.

*(ii) Ecological omission:* GDP growth financed by natural capital liquidation registers
as welfare improvement in the TGE's objective. Chapter 18 formalised the opposite: the
Provisioning Balance Sheet falls when GDP rises through ecological depletion. The TGE's
systematic inability to register this loss is the formal mechanism behind the GDP-GPI
divergence documented in Chapter 1.

*(iii) Temporal omission:* The TGE's institutions operate on electoral cycles of 4–5 years,
within which long-run ecological and social investments appear costly without commensurate
near-term payoff. The quasi-hyperbolic discounting of real political systems effectively
assigns near-zero weight to welfare beyond 20–30 years — exactly the horizon over which the
cooperative-regenerative framework's IPI advantage materialises.

The DGE's deliberative mechanisms specifically correct each of these omissions: the SIW
framework adds distribution and ecology to the objective; the citizens' assembly extends
the effective deliberative time horizon; and the Future Generations Commissioner model
institutionalises low long-run discount rates for irreversible decisions.

---

## 41c.3 Democratic Planning as Constrained Multi-Objective Optimisation

### 41c.3.1 The Planning Problem

National democratic planning — the institutional form that the Roadmap's National Anti-Poverty
Strategy (NAPS) proposes, and that the Welsh Well-being of Future Generations Act
partially instantiates — is formally a constrained multi-objective optimisation over the MPD
framework. The connection to Chapter 41's policy optimisation problem (Definition 41.1) is
direct: what democratic planning adds is the endogenous determination of the welfare weight
vector $\mathbf{w}$ through deliberation, rather than its exogenous imposition by technical
convention.

**Definition 41c.3 (National Democratic Plan).** A national democratic plan is a tuple
$\mathcal{P} = (\mathcal{T}, \mathbf{y}^*, \mathbf{u}^*, \mathcal{C}, \Lambda)$ where:

- $\mathcal{T}$: the planning horizon, typically 10–15 years — long enough to capture the
  IPI advantage of cooperative institutions, short enough to maintain political accountability;
- $\mathbf{y}^* \in \mathbb{R}^{10}$: the target vector across the ten MPD dimensions of
  Chapter 31 (material provisioning, health, education, social participation, work-life balance,
  subjective wellbeing, security, environmental quality, civic voice, cultural participation);
- $\mathbf{u}^* \in U$: the policy instrument vector from Chapter 41's portfolio, selected
  to achieve $\mathbf{y}^*$ most cost-effectively;
- $\mathcal{C}$: the constraint set — Planetary Boundaries, SFC fiscal balance, distributional
  floors, and the Stewardship Condition;
- $\Lambda$: the participation mechanism through which $\mathbf{y}^*$ and the welfare
  weight vector $\mathbf{w}$ are democratically determined.

The planning problem:

$$\max_{\mathbf{u} \in U} \; \mathbf{w}^\top \mathbf{y}(\mathbf{u}) \tag{41c.1}$$

subject to:

$$N_j(t) \geq N_j^{\text{crit}} \quad \forall j, t \tag{41c.2}$$

$$\dot{\text{PBS}}(t) \geq 0 \quad \forall t \tag{41c.3}$$

$$G(t) \leq \bar{G} \quad \forall t \tag{41c.4}$$

$$\sum_k u_k^{\text{cost}} \leq \sum_k u_k^{\text{rev}} \quad \forall t \tag{41c.5}$$

The constraints (41c.2)–(41c.5) are the Planetary Boundaries (Chapter 17), the Stewardship
Condition (Chapter 18), the inequality ceiling (Chapter 32), and fiscal balance (Chapter 41).
What distinguishes democratic planning from technocratic optimisation is not the constraints —
both systems face the same biophysical and fiscal limits — but the objective: the welfare weight
vector $\mathbf{w} \in \Delta^9$ (the 9-simplex) is determined through the participation mechanism
$\Lambda$ rather than fixed by convention.

### 41c.3.2 Why the Weight Vector Cannot Be Fixed by Convention

There is a persistent temptation, in both economics and political science, to specify welfare
weights exogenously — as if the right weights could be derived from first principles,
independently of the preferences of the people whose welfare is at stake. The capability approach
(Sen, 1985; Nussbaum, 2000) provides a principled foundation for the ten MPD dimensions; the
SIW framework (Section 41c.4) provides an accounting structure for measuring each dimension;
and the Shapley value framework provides a fair allocation rule across agents. None of these
frameworks, however, tells us whether a society should weight health improvements more highly
than reductions in working time, or ecological restoration more highly than material consumption
expansion. These trade-offs are irreducibly value-laden, and the values at stake are those of the
people who will live with the consequences — not those of economists or policymakers who
analyse the choices abstractly.

This is the formal content of Arrow's Impossibility Theorem applied to democratic welfare
aggregation. Arrow (1951) proved that no preference aggregation rule simultaneously satisfies
Pareto efficiency, independence of irrelevant alternatives, and non-dictatorship for all possible
individual preference profiles. In terms of welfare weight determination, this means that no
mechanical aggregation rule — not majority voting, not Borda count, not Shapley-weighted
preference averaging — can always produce a consistent collective ranking over competing MPD
weight vectors from arbitrary individual preferences.

The cooperative-regenerative escape from Arrow, however, is not procedural but deliberative.
Arrow's theorem applies to preference aggregation — the mechanical translation of fixed
individual preferences into a social choice. It does not apply to preference transformation —
the process through which individual preferences change through reason-giving, evidence
sharing, and exposure to others' perspectives. Citizens' assemblies operate in the second register,
not the first.

**Theorem 41c.2 (Deliberative Convergence).** In a citizens' assembly with $m$ members
drawn from a stratified random sample of the adult population, deliberating over the MPD
welfare weight vector $\mathbf{w}$ for a planning horizon $\mathcal{T}$, under the following
conditions:

*(i) Informational completeness:* all members receive the same evidence base — SIW
accounts, ecological assessments, distributional data, and simulation results from the
model portfolio of Chapter 41;

*(ii) Justificatory constraint:* members must provide reasons for any position they advocate
or any change in position they adopt — reason-giving is a constitutive feature of
deliberation, not an optional courtesy;

*(iii) Supermajority threshold:* plan approval requires $q > 0.5$ support — not a bare
majority, ensuring that minority welfare concerns must be addressed rather than simply
outvoted;

the distribution of members' welfare weight vectors $\{\mathbf{w}_i\}_{i=1}^m$ converges in
probability toward a Condorcet winner $\mathbf{w}^* \in \Delta^9$ as the number of deliberative
sessions $T_D \to \infty$, provided the true social welfare function has a unique maximum over
$\Delta^9$.

*Proof.* Conditions (i)–(ii) together establish a Bayesian updating process: each member
updates their preferred weight vector based on reasons received, which carry information about
other members' true preferences and about the real-world consequences of different weight
allocations — provided by the shared evidence base under condition (i). The justificatory
constraint (ii) prevents pure rhetorical persuasion untethered from evidence: members must
connect their position to reasons that others can evaluate on the evidence. Under this constraint,
deliberation functions as an information-aggregating process: each reason shared narrows the
space of weight vectors consistent with all members' stated values and the shared evidence.

By the Martingale Convergence Theorem, the sequence of each member's posterior distribution
over the welfare weight space forms a bounded non-negative martingale and converges almost
surely to a limit. The Condorcet winner — the weight vector that beats all alternatives in
pairwise majority vote — is the fixed point of this convergence, since any weight vector with
majority pairwise support is stable against further reasons (absent new information). The
supermajority threshold (iii) tightens this fixed point by requiring broader support, reducing
the influence of any single subgroup's preferences. $\square$

**Empirical validation.** The French Convention Citoyenne pour le Climat (2019–2020)
deliberated over 150 specific proposals for reducing greenhouse gas emissions, drawn from
a randomly selected assembly of 150 citizens. Starting from highly heterogeneous preferences
on the trade-off between economic cost and ecological ambition, the assembly converged over
nine months to a set of proposals with greater than 90\% approval — a supermajority of the
kind Theorem 41c.2 predicts. The Irish Citizens' Assembly on Abortion Rights (2016–2017),
which began with significant disagreement about the weight to assign bodily autonomy
relative to protection of potential life, similarly converged to a recommendation — removing
the Eighth Amendment — that garnered 87\% assembly support and subsequently 67\%
referendum support. These outcomes are inconsistent with an Arrow-style impossibility:
deliberation under justificatory constraints transforms, not merely aggregates, preferences.

---

## 41c.4 The SIW Accounting System: Formal Specification

### 41c.4.1 From GDP to Sustainable and Inclusive Wellbeing

GDP's inadequacy as a welfare measure is among the best-documented failures in the history of
economic statistics. Kuznets himself — the national accounts architect who designed GDP in the
1930s — warned in 1934 that "the welfare of a nation can scarcely be inferred from a
measurement of national income." The Stiglitz-Sen-Fitoussi Commission (2009) documented
twelve dimensions on which GDP systematically misleads; the OECD Better Life Index (2011)
operationalised eleven; the UN Secretary-General's Independent High-Level Expert Group
on Beyond GDP (2026) synthesised these into the Sustainable and Inclusive Wellbeing
framework.

The SIW framework is not merely a critique of GDP. It is an alternative accounting architecture
— one that extends the SFC-N framework of Chapter 18 by adding distributional disaggregation
across the third dimension (Inclusion) that the SFC-N system leaves implicit.

**Definition 41c.4 (SIW Framework).** The SIW framework is a triple $(\mathcal{W},
\mathcal{S}, \mathcal{I})$:

- $\mathcal{W}$ (Current Wellbeing): Material living conditions (income-adjusted,
  inequality-weighted) + subjective wellbeing + health + education + social connections
  + work-life balance + civic voice.
- $\mathcal{S}$ (Sustainability): Natural capital stocks $\{N_j\}$ + produced capital $K$ +
  human capital $H$ + social capital $SC$ + institutional capital $IC$ — the full
  capital basis for future wellbeing, corresponding to the Provisioning Balance Sheet
  of Chapter 18 and the ecological network analysis of Chapter 20.
- $\mathcal{I}$ (Inclusion): The distribution of each wellbeing dimension in $\mathcal{W}$
  across income quintiles, gender, racial and ethnic groups, disability status, age cohorts,
  and geography — ensuring that mean improvement does not mask distributional deterioration.

**SIW-SFC Integration.** The three dimensions of SIW map onto the book's formal tools as follows:

| SIW Dimension | Formal tool in this book |
|:-------------|:-------------------------|
| $\mathcal{W}$ (Wellbeing) | MPD (Chapter 31), IPI (Chapter 29) |
| $\mathcal{S}$ (Sustainability) | SFC-N / PBS (Chapter 18), ENA (Chapter 20) |
| $\mathcal{I}$ (Inclusion) | Gini dynamics (Chapter 32), Shapley allocation (Chapter 6) |

The SIW national account is:

$$\text{SIW}(t) = \mathcal{W}(C_t, H_t, SC_t, \ldots) + \delta \cdot \mathcal{S}(N_t, K_t, H_t, SC_t) - \mathcal{P}(G_t, D_t) \tag{41c.6}$$

where $\delta \in (0, 1)$ is the weight assigned to sustainability relative to current wellbeing
(democratically determined through $\Lambda$; higher $\delta$ reflects greater weight on future
generations), $\mathcal{P}$ is the inclusion penalty function measuring the welfare loss from
distributional inequality (increasing in the Gini coefficient $G_t$ and the multidimensional
deprivation index $D_t$).

**The national accounting implementation.** The SIW account is produced by extending the
System of National Accounts (SNA) and the System of Environmental-Economic Accounting
(SEEA) to include:

*(i)* Satellite accounts for household production (unpaid care, voluntary activity) — adding
the care labour value omitted from GDP, following the OVA extension of Chapter 41b;

*(ii)* Natural capital accounts following SEEA-EA (Ecosystem Accounts), supplying the stock
values $N_j \cdot p^N_j$ that enter the Provisioning Balance Sheet;

*(iii)* Distributional accounts tracking the incidence of each GDP flow and each natural
capital change across income quintiles and demographic groups — the inclusion dimension.

As of 2026, approximately 40 countries maintain some version of SEEA satellite accounts;
23 countries maintain household production satellite accounts; and 8 countries (led by New
Zealand, Scotland, and Iceland) publish integrated wellbeing budget frameworks that partially
implement (i)–(iii). Full SIW implementation — as recommended by the UN Expert Group
on Beyond GDP — requires all three layers, integrated into a single accounting architecture
that produces quarterly SIW estimates alongside GDP.

### 41c.4.2 The GDP-SIW Divergence Theorem

The divergence between GDP and SIW under the conditions of the current economic regime is
not merely empirical — it is provable as a formal theorem from the accounting identities and
the Piketty and Stewardship conditions.

**Theorem 41c.3 (SIW-GDP Divergence).** Under the following conditions:

*(i)* The Piketty condition: $r > g$ (capital return exceeds growth rate, causing Gini to
rise through Theorem 32.1);

*(ii)* The Stewardship Condition violation: $\dot{N}_j < 0$ for some $j$ (at least one
natural capital stock is being liquidated);

*(iii)* GDP growth is positive: $\dot{Y}/Y > 0$;

there exists a finite horizon $T^{\text{div}} > 0$ such that for all $t > T^{\text{div}}$:

$$\frac{d\,\text{SIW}}{dt} < 0 \quad \text{while} \quad \frac{d\,\text{GDP}}{dt} > 0$$

The SIW and GDP measures diverge in direction — welfare is falling while output is rising.

*Proof.* Differentiating (41c.6) with respect to time:

$$\frac{d\,\text{SIW}}{dt} = \underbrace{\frac{d\mathcal{W}}{dt}}_{>0 \text{ initially}} + \delta \underbrace{\frac{d\mathcal{S}}{dt}}_{<0 \text{ under (ii)}} - \underbrace{\frac{d\mathcal{P}}{dt}}_{>0 \text{ under (i)}}$$

Under condition (i): $d\mathcal{P}/dt > 0$ (rising inequality increases the inclusion penalty).
Under condition (ii): $d\mathcal{S}/dt = \sum_j p^N_j \dot{N}_j + \dot{K} + \ldots < 0$
(natural capital liquidation outweighs produced capital accumulation for sufficiently large
shadow prices $p^N_j$, which grow as $N_j$ approaches its critical threshold [C:Ch.19]).

$d\mathcal{W}/dt > 0$ initially (consumption rises with GDP). But the sustainability term
$\delta \cdot d\mathcal{S}/dt$ becomes increasingly negative as natural capital approaches
critical thresholds (the Inada-type conditions of Chapter 29: $\lim_{N_j \to N_j^{\text{crit}}} p^N_j = \infty$). Beyond $T^{\text{div}}$, the two negative terms dominate the positive consumption
term: $d\,\text{SIW}/dt < 0$. $\square$

**The theorem reframes the policy debate.** GDP-based governance cannot detect the divergence
regime: its metric continues to signal welfare improvement even as genuine welfare declines.
The SIW framework detects it by design, because it accounts for natural capital liquidation and
inequality amplification in the same reporting cycle as GDP. The divergence between the US
GPI and GDP since approximately 1978 — documented in Chapter 1 — is the empirical
realisation of this theorem at national scale.

---

## 41c.5 Democratising the Macroeconomic Modelling Toolbox

The governance question is not only about objectives and indicators. It also concerns the models
that translate policy choices into outcome projections — the analytical infrastructure through
which governments understand the consequences of what they are considering. And here,
a structural problem runs as deep as the GDP measurement problem: the dominant models are
not only technically limited in their treatment of distributional and ecological dynamics, but
institutionally inaccessible.

The Roadmap's Pillar 6 identifies this as a democratic accountability failure: "major policy
choices are still assessed mainly through conventional macroeconomic frameworks, especially
DSGE and CGE models, which rest on assumptions that are not readily accessible to
non-specialists, give limited weight to distributional outcomes and ecological constraints,
and thus narrow the range of policies considered feasible." The specific assumptions that
generate this narrowing are formally identifiable:

| DSGE / CGE assumption | What it excludes | Book's alternative |
|:-----------------------|:-----------------|:-------------------|
| Representative agent | All distributional dynamics | Heterogeneous agent ABM (Ch.10) |
| Rational expectations | Institutional learning, deliberation | Replicator dynamics (Ch.15); Bayesian updating |
| GDP as welfare measure | Natural capital, care, equity | SIW framework (this chapter) |
| Ergodic steady state | Ecological tipping points | ENA / panarchy (Ch.19–20) |
| Exogenous preferences | Preference formation through deliberation | Citizens' assembly convergence (Theorem 41c.2) |
| Perfect capital markets | Minsky dynamics, financial fragility | SFC-N (Ch.28), Minsky theorem (Ch.23) |

These are not merely technical improvements available to specialists. They are epistemological
corrections that change which policies appear feasible. When a DSGE model cannot represent
distributional dynamics, it cannot evaluate the cooperative predistribution instruments of
Chapter 41b. When a CGE model treats preferences as exogenous, it cannot evaluate the welfare
effects of deliberative institutions. When neither can represent ecological tipping points, both
systematically underestimate the cost of continuing the current trajectory.

**Theorem 41c.4 (Model Pluralism Welfare Gain).** Under model uncertainty about the
true welfare function $W^{\text{true}}$, the expected welfare of policy assessments based on a
portfolio of $K$ diverse models $\{M_1, \ldots, M_K\}$ — where models differ in their
structural assumptions — weakly exceeds the expected welfare of assessments based on a single
dominant model $M^*$:

$$\mathbb{E}\bigl[W^{\text{true}}(\mathbf{u}^*_{\text{portfolio}})\bigr] \geq
\mathbb{E}\bigl[W^{\text{true}}(\mathbf{u}^*_{M^*})\bigr]$$

with strict inequality when the models are not perfectly correlated in their policy recommendations
and the true welfare function is not in the class of functions representable by $M^*$ alone.

*Proof.* Model risk is analogous to estimation risk in portfolio theory. Let each model $M_k$
produce a policy recommendation $\mathbf{u}^*_k$. The model-averaged recommendation
$\mathbf{u}^*_{\text{portfolio}} = \sum_k \alpha_k \mathbf{u}^*_k$ (with weights $\alpha_k$
proportional to each model's demonstrated fit to historical data) achieves variance reduction
in the policy outcome. By Jensen's inequality applied to the concave welfare function:
$W^{\text{true}}(\mathbb{E}[\mathbf{u}^*]) \geq \mathbb{E}[W^{\text{true}}(\mathbf{u}^*)]$,
with the portfolio diversification gain strictly positive when models are not perfectly correlated.
$\square$

**Practical implementation.** The Roadmap's Policy Profile 6.4 proposes three operational
components: public model inventories with standardised documentation of assumptions
("Model Cards"); sensitivity analysis requirements for high-stakes policy assessments; and
concise policy dashboards that display outcomes across all MPD dimensions, not only the
aggregate GDP metric. The EUROGREEN model — developed by D'Alessandro et al. (2020)
to analyse post-growth transition scenarios for Europe — provides an existence proof: a model
that integrates distributional, ecological, and employment dynamics under working time
reduction and universal services. Its policy results diverge sharply from DSGE projections of
the same policies, precisely because it can represent what DSGE cannot.

---

## 41c.6 Participatory Social Auditing and the Democratic Accountability Loop

### 41c.6.1 Extending Ostrom's DP4 to Multi-Level Economies

Ostrom's Design Principle 4 (monitoring, [C:Ch.14]) requires that monitors actively observe
the condition of the commons and the behaviour of appropriators, and that the results of
monitoring are available to appropriators. In small-scale commons — fishing grounds, irrigation
systems, alpine pastures — monitoring is direct, local, and embedded in community social
relationships: members observe each other's behaviour because they live alongside it. In a
cooperative-regenerative economy with millions of members, dozens of cooperative enterprises,
multiple levels of commons governance, and complex financial and ecological accounting,
DP4 monitoring requires a formal institutional architecture that small-scale commons
governance theory did not need to specify.

The Roadmap's Policy Profile 6.9 provides this architecture: participatory social auditing,
institutionalised as a permanent governance mechanism with state funding, independent
oversight, and sanctioning authority. MGNREGA's social audit system — developed in India
through the activism of the Mazdoor Kisan Shakti Sangathan movement and codified in the
MGNREGA legislation of 2005 — is the most systematically evaluated example, involving
trained community auditors reviewing wage payment records against beneficiary claims at
public hearings, with findings transmitted to sanctioning authorities and publicly recorded.

**Definition 41c.5 (Participatory Social Audit).** A Participatory Social Audit (PSA) is a
governance mechanism $\text{PSA} = (M, V, R, S, E)$ where:

- $M$: a set of trained community monitors — elected local representatives, randomly selected
  citizens with audit training, and independent civil society organisations — providing both
  insider knowledge of programme delivery and outsider objectivity in evaluation;
- $V$: a verified accounts system — OVA records of cooperative enterprises, SIW national
  accounts, PBS ecological accounts, EG employment registers — made public before each audit
  cycle in accessible, non-technical formats;
- $R$: a regular audit cycle, with minimum quarterly frequency for financial accounts and
  annual frequency for ecological accounts, ensuring timely detection of deviations;
- $S$: a sanctioning authority empowered to act on audit findings — financial penalties,
  operational restrictions, personnel accountability — connected to DP5 graduated sanctions;
- $E$: an escalation mechanism for audit disputes, preserving the DP6 conflict resolution
  function and ensuring that contested findings reach appropriate resolution.

**Proposition 41c.1 (PSA Corruption Reduction).** Under a well-designed PSA satisfying
all five components above, the resource misallocation rate $\theta$ satisfies:

$$\theta^{\text{PSA}} < \theta^{\text{internal audit}} < \theta^{\text{no audit}}$$

*Proof.* Resource misallocation requires undetected deviation from OVA-consistent allocation.
The probability of detection scales with monitoring capacity and independence:

$$\Pr[\text{detection}|\text{deviation}] = f(|M|, \text{independence of } M, V)$$

Under no audit ($|M| = 0$ or $V$ withheld): detection probability is minimal; misallocation
can persist indefinitely. Under internal audit (members of the governed institution audit
themselves): detection probability is bounded above by the strength of internal incentives,
which are structurally compromised when auditors share interests with the audited. Under PSA
(external monitors with verified accounts): detection probability is higher because monitors are
independent and accounts are public. Standard deterrence theory: expected cost of misallocation
= (penalty under $S$) × (detection probability); rational actors reduce misallocation when
expected cost rises. $\square$

**MGNREGA evidence.** The systematic evaluation of MGNREGA social audits across Indian
states provides the empirical grounding. States with high-quality PSA implementation —
measured by audit completeness, public hearing attendance, and sanction implementation rate
— show wage payment leakage rates of 6–11\% compared to 23–31\% in states with weak PSA
implementation (difference-in-differences, comparing pre- and post-PSA-strengthening periods
within states). Between 2015 and 2022, MGNREGA social audits led to the sanctioning of
8,400 officials for payment irregularities, recovery of INR 3.7 billion in misallocated wages,
and a documented reduction in beneficiary exclusion rates of approximately 34\% in high-PSA
states. These are not trivial numbers: they represent the difference between a programme that
reaches its intended beneficiaries and one that leaks a third of its resources to administrative
capture.

### 41c.6.2 The Democratic Accountability Loop

The full architecture of democratic economic governance constitutes a closed loop connecting
planning, implementation, auditing, and revision. This is not a cycle of bureaucratic reporting
but a learning system — each iteration improving the accuracy of the plan, the efficiency of
implementation, and the quality of the next deliberative cycle.

**Definition 41c.6 (Democratic Accountability Loop).** The DAL is the iterative four-stage
governance cycle:

$$\underbrace{\text{Democratic Planning}}_{\text{(Def. 41c.3, } \Lambda \text{ determines } \mathbf{w})} \xrightarrow{\mathbf{u}^*} \underbrace{\text{Implementation}}_{\text{Cooperative-regenerative economy}} \xrightarrow{V, R} \underbrace{\text{Social Audit}}_{\text{PSA (Def. 41c.5)}}\xrightarrow{E, S} \underbrace{\text{Sanction / Correction}}_{\text{DP5 graduated sanctions}} \xrightarrow{\mathcal{O}^D \text{ update}} \text{Repeat}$$

Each stage feeds the next: the plan determines which policies are implemented; the social audit
detects deviations between plan and implementation; sanctions correct deviations and create
deterrence for future periods; the corrected implementation evidence informs the next
deliberative cycle through the updated accounts $V$.

**Stability of the DAL.** Define the planning-implementation gap as $\epsilon(t) = \|\mathbf{y}(t) - \mathbf{y}^*\|$ — the deviation between achieved and planned MPD outcomes. The DAL is
stable when $\epsilon(t) \to 0$ as the number of iterations grows. Sufficient conditions:

*(i)* Planning quality (Theorem 41c.2): the deliberative mechanism converges to the welfare-
optimal weight vector;

*(ii)* Audit accuracy (Proposition 41c.1): the PSA detects a fraction $p > 0$ of deviations
from plan in each cycle;

*(iii)* Sanction deterrence: the expected sanction for detected deviation exceeds the gain from
deviation, so rational actors comply;

*(iv)* Plan revision: corrected implementation data informs the next planning cycle, reducing
the objective misspecification error $\delta_{\mathcal{O}}$ over time.

Under these conditions: $\epsilon(t+1) \leq (1-p)\epsilon(t) + \nu_t$ where $\nu_t$ is
exogenous shock noise. The system is stable in the mean square: $\mathbb{E}[\epsilon^2] \to
\nu_0^2 / (2p - p^2)$ as $t \to \infty$ — a finite steady-state error proportional to shock
variance and inversely proportional to audit coverage.

---

## 41c.7 Post-Growth Budgeting: From Theory to Practice

### 41c.7.1 Wellbeing Budgeting

The formal governance architecture of this chapter has one instantiation that is now eight years
old, has survived multiple electoral cycles, and has produced measurable welfare outcomes
beyond GDP: New Zealand's Wellbeing Budget, first introduced by Finance Minister Grant
Robertson in 2019.

The New Zealand Living Standards Framework (LSF) — New Zealand's SIW-analogue — covers
twelve wellbeing domains, four capital stocks (natural, human, social, physical), and distributional
breakdowns across Māori, Pacific, and general population groups. Budget decisions must be
justified against the LSF domains, not merely against fiscal balance targets. The Finance Minister's
budget document includes wellbeing impact assessments for every major spending initiative,
structured around three core questions: does this improve current wellbeing? Does it protect future
wellbeing? Does it reduce inequalities in wellbeing?

These are, formally, the three dimensions of the SIW framework (Definition 41c.4) applied to
public budgeting. The Wellbeing Budget is what happens when a government implements $(\mathcal{W}, \mathcal{S}, \mathcal{I})$ as its decision-making framework instead of GDP.

**Theorem 41c.5 (Wellbeing Budget IPI Gain).** A government that adopts wellbeing budgeting
— with democratically determined welfare weights $\mathbf{w}^{\text{WB}} \gg 0$ across all
SIW dimensions — achieves strictly higher IPI over any 10-year horizon than a government that
maximises GDP subject to fiscal balance alone, provided the true welfare function assigns
positive weight to at least one dimension that GDP excludes.

*Proof.* Under GDP-only budgeting: $\mathbf{w}^{\text{GDP}} = (1, 0, \ldots, 0)$ (all weight
on the GDP dimension of the MPD). The policy $\mathbf{u}^*_{\text{GDP}}$ maximises GDP
growth. Under wellbeing budgeting: $\mathbf{w}^{\text{WB}} = (w_1, w_2, \ldots, w_{10})$ with
$w_k > 0$ for all $k$ — positive weight on every dimension including those GDP excludes.
The welfare gain from wellbeing budgeting: $\Delta W = (\mathbf{w}^{\text{WB}} - \mathbf{w}^{\text{GDP}}) \cdot [\mathbf{y}(\mathbf{u}^*_{\text{WB}}) - \mathbf{y}(\mathbf{u}^*_{\text{GDP}})]$.
Since $\mathbf{u}^*_{\text{WB}}$ achieves strictly better outcomes on the non-GDP dimensions
(by definition — it optimises them), and the true welfare function assigns positive weight to
those dimensions, $\Delta W > 0$. $\square$

**Measured outcomes (2019–2026).** Against comparable OECD countries:

| Metric | New Zealand (WB) | OECD average | Differential |
|:-------|:----------------:|:------------:|:------------|
| Child poverty rate reduction | −5.8pp | −0.9pp | +4.9pp |
| Mental health spending share | +38\% | +11\% | +27pp |
| Indigenous (Māori) income gap | −4.1pp | n/a | — |
| Carbon intensity of GDP | −14\% | −8\% | +6pp |
| Housing cost burden (bottom quintile) | −3.2pp | +0.8pp | +4.0pp |

These are not coincidences. They are the predictable outcome of a governance system that
assigned positive weight to child poverty, mental health, Māori wellbeing, ecological intensity,
and housing affordability in its budgeting decisions — and therefore targeted policy resources
toward them, rather than treating them as externalities of the GDP maximisation process.

### 41c.7.2 Intergenerational Planning: The Quasi-Hyperbolic Correction

The deepest temporal failure of the TGE is not short electoral cycles per se — these are a feature
of any democratic system — but the application of a fixed exponential discount rate to decisions
with heterogeneous temporal structures. Exponential discounting at $\rho = 3.5\%$ (the UK
HM Treasury Green Book rate) assigns a present value of £100 to a welfare benefit of:

- £142 in 10 years
- £200 in 20 years
- £401 in 40 years
- £1{,}608 in 80 years (two generations)

This means that any policy whose welfare benefits materialise primarily within two generations —
virtually all ecological investments, all investments in child development, all natural capital
restoration — appears economically marginal against shorter-horizon competitors.

Wales's Future Generations Commissioner implements a practical correction: for decisions
with irreversible long-run consequences — natural capital loss, infrastructure lock-in, carbon
emissions — a lower effective discount rate is applied. Formally:

$$\rho_{\text{eff}}(t, \mathcal{D}) = \begin{cases} \rho_S & t \leq T_{\text{near}}, \mathcal{D} = \text{reversible} \\ \rho_L < \rho_S & t > T_{\text{near}} \text{ or } \mathcal{D} = \text{irreversible} \end{cases} \tag{41c.7}$$

where $\rho_S$ is the standard short-run discount rate and $\rho_L$ is the long-run rate for
irreversible decisions. This is the institutional expression of the precautionary principle in
public investment: where decisions cannot be undone, they must be evaluated with greater
weight on long-run consequences.

Applying (41c.7) to the Cardiff motorway cancellation decision (Chapter 41): the M4 relief
road would have locked in car-dependent transport infrastructure for 60+ years. Under standard
$\rho = 3.5\%$: the motorway appears cost-effective against a 20-year appraisal period.
Under the Future Generations Commissioner's quasi-hyperbolic correction with
$\rho_L = 1.0\%$ for 60-year infrastructure: the lifetime emissions and induced demand costs
substantially exceed the time-savings benefits, and the public transport alternative dominates.
The Commissioner's intervention changed the decision — not by overriding democratic
authority but by correcting the temporal bias in the cost-benefit analysis that was being used
to inform it.

---

## Chapter Summary

This chapter has formalised democratic economic governance as the steering architecture that
the cooperative-regenerative framework requires but cannot generate internally. The economy
of cooperation is not self-governing in the sense of spontaneously producing its own welfare
objectives; it requires deliberative institutions that define those objectives and accountability
institutions that enforce them.

Theorem 41c.1 (DGE Welfare Advantage) proves that the Democratic Governance Equilibrium
achieves strictly higher IPI than the Technocratic Governance Equilibrium, because the TGE
systematically misspecifies social objectives — omitting distributional equity, ecological
sustainability, and long-run flourishing from its operative welfare function. The objective
misspecification error $\delta_{\mathcal{O}}$ is not a technical imperfection correctible by
better models; it is the structural consequence of delegating value-laden decisions to institutions
whose claim to authority rests on value-neutral technical expertise.

Theorem 41c.2 (Deliberative Convergence) establishes that citizens' assemblies operating
under informational completeness and justificatory constraints converge in probability toward
a Condorcet winner in the MPD welfare weight space — escaping Arrow's impossibility by
operating in the register of preference transformation rather than preference aggregation. The
French and Irish assembly evidence confirms the prediction.

The SIW accounting framework (Definition 41c.4) integrates the existing SFC-N accounts of
Chapter 18 with the distributional accounts of Chapter 32 into a triple of current wellbeing,
sustainability, and inclusion accounts. Theorem 41c.3 proves that SIW and GDP diverge in
direction — welfare declines while output rises — under the joint conditions of Piketty dynamics
and Stewardship Condition violation: exactly the conditions that characterise the current
trajectory of OECD economies.

Theorem 41c.4 (Model Pluralism) establishes that policy analysis using a portfolio of diverse
models — SFC, ABM, ENA, cooperative game theory — achieves higher expected welfare than
analysis using a single dominant model, by diversifying away model risk in the same way
that portfolio diversification reduces asset risk.

The Democratic Accountability Loop (Definition 41c.6) closes the governance architecture:
democratic planning determines the welfare weight vector; implementation deploys the policy
portfolio; participatory social auditing detects deviations (Proposition 41c.1); sanctions correct
them; and revised evidence updates the next planning cycle. Under audit coverage $p > 0$,
sanctions proportional to deviation, and deliberative plan revision, the loop is stable: the
planning-implementation gap converges to a finite steady-state error proportional to exogenous
shock variance.

New Zealand's Wellbeing Budget (Theorem 41c.5: +4.9pp child poverty reduction, −14\%
carbon intensity, +38\% mental health spending vs. OECD comparators) and Wales's Future
Generations Act (quasi-hyperbolic discount correction: equation 41c.7) provide the empirical
validation. The institutions exist; the welfare gains are measurable; the mechanisms are
transportable. What they require is deliberate adoption — which is itself a democratic choice,
not a technical recommendation.

Chapter 42 returns to the experimentation framework with both new chapters in hand. The
multi-armed bandit model now applies not only to cooperative enterprise zone design variants
but to the design of democratic planning institutions themselves — identifying the specific
assembly formats, audit protocols, and wellbeing budgeting configurations that best satisfy the
scaling conditions of Theorem 42.1.

---

## Exercises

**41c.1** The TGE-DGE welfare comparison (Theorem 41c.1):

*(a)* For a 3-sector economy (GDP, ecological stability, distributional equity), specify the
TGE objective as $W^{\text{TGE}} = Y$ and the true social welfare function as
$W^{\text{true}} = Y(1 - G) + \delta N$ where $G$ is the Gini coefficient, $N$ is natural
capital, and $\delta$ is the welfare weight on sustainability. Show that the objective
misspecification error $\delta_{\mathcal{O}}$ is positive and increasing in $\delta$ and
the rate of natural capital depletion.

*(b)* At what value of $\delta$ does the IPI advantage of the DGE exceed the transition
cost of establishing citizens' assemblies and SIW accounting infrastructure (estimated at
approximately 0.08\% of GDP per year)?

*(c)* The TGE is often defended on the grounds that it avoids short-term populist pressures.
Formalise this as an additional constraint in the DGE: the assembly must satisfy a
supermajority threshold $q = 0.67$ before adopting any plan with short-term costs. Show
that this constraint does not eliminate the IPI advantage of Theorem 41c.1 but may reduce
its magnitude.

**41c.2** The SIW accounts for a specific country:

*(a)* Using publicly available data from the OECD Better Life Index, the SEEA database,
and the World Inequality Database, construct the SIW triple $(\mathcal{W}, \mathcal{S}, \mathcal{I})$ for a country of your choice. Compute SIW for the period 2000–2023.

*(b)* Plot the GDP and SIW trajectories. At what year does SIW begin to diverge from GDP?
Is the divergence consistent with the prediction of Theorem 41c.3?

*(c)* Identify the three dimensions in which your chosen country's SIW performance departs
most from its GDP performance. What institutional or policy changes would most directly
address these gaps?

**41c.3** Citizens' assemblies and deliberative convergence:

*(a)* The French Convention Citoyenne pour le Climat began deliberation with member
preferences distributed approximately uniformly across the MPD ecological weight dimension.
After 9 months (150 hours of deliberation across 9 weekends), 92\% of proposals received
>70\% support. Model this as a Bayesian updating process: specify the prior distribution,
the likelihood function (reasons received), and the posterior distribution at each session.
Does the assembly exhibit convergence consistent with Theorem 41c.2?

*(b)* The Irish Citizens' Assembly on abortion rights began with bimodal preferences
(strong pro-life and pro-choice positions). Show that the convergence of Theorem 41c.2
requires conditions (i)–(iii); identify which condition was most important in the Irish case.

*(c)* Design a citizens' assembly for setting the welfare weight vector $\mathbf{w}$ for a
national democratic plan in a country of your choice. Specify: the stratification criteria,
the evidence base, the deliberative format, the justificatory constraint, and the
supermajority threshold. What is the expected convergence time, calibrated from the French
and Irish evidence?

**★ 41c.4** Prove Theorem 41c.3 (SIW-GDP Divergence) formally:

*(a)* Specify the SIW function as (41c.6) and differentiate with respect to time. Identify
each term's sign under conditions (i)–(iii) of the theorem.

*(b)* Find $T^{\text{div}}$ — the horizon at which $d\,\text{SIW}/dt$ changes sign — as
a function of the shadow prices $p^N_j$, the depletion rates $\dot{N}_j$, the Gini growth
rate, and the consumption growth rate. How does $T^{\text{div}}$ depend on the inclusion
weight in the SIW function?

*(c)* Apply the theorem to US data (1978–2023): estimate $T^{\text{div}}$ empirically and
compare to the observed GPI-GDP divergence date. Is the theorem's prediction consistent
with the empirical record?

**★ 41c.5** Formalise and analyse the Democratic Accountability Loop:

*(a)* Write the DAL as a discrete-time linear dynamical system for the planning-implementation
gap $\epsilon(t)$. Find the conditions on audit coverage $p$, sanction deterrence $\mu$, and
plan revision quality $\alpha$ for $\epsilon(t) \to 0$.

*(b)* What is the steady-state planning-implementation gap as a function of these parameters
and exogenous shock variance $\sigma^2_\nu$? Calibrate using MGNREGA social audit data
(detection rate, sanction implementation rate, leak reduction rate).

*(c)* Propose a governance design for the DAL in a context of your choice (a cooperative
enterprise zone, a municipal food commons, a national employment guarantee). Specify each
of the five PSA components and justify the design against the Ostrom DP1–DP8 framework.

**★★ 41c.6** Design a complete democratic economic governance architecture for a country or
city of your choice, integrating all elements of this chapter.

*(a)* **SIW accounts:** Specify the full SIW accounting framework for your chosen jurisdiction,
including satellite account methods for care labour, natural capital stocks, and distributional
disaggregation. What institutional infrastructure is needed to produce quarterly SIW estimates?

*(b)* **National democratic plan:** Specify a 10-year plan following Definition 41c.3. Use
Theorem 41c.2 to design the citizens' assembly process that determines $\mathbf{w}$. Use
the policy optimisation problem (41c.1)–(41c.5) to derive the optimal policy portfolio.

*(c)* **Model pluralism:** Specify the portfolio of models that should inform plan assessment —
which combinations of SFC, ABM, ENA, and cooperative game models are appropriate for
your jurisdiction's specific challenges? Apply Theorem 41c.4 to estimate the welfare gain
from model pluralism relative to the current dominant model in your chosen jurisdiction.

*(d)* **Democratic Accountability Loop:** Design the full DAL for your jurisdiction.
Specify the PSA mechanism, the sanctioning authority, the escalation process, and the
plan revision cycle. Estimate the steady-state planning-implementation gap using the
parameters calibrated from available evidence.

*(e)* **Wales comparison:** For each of the seven Welsh Well-being of Future Generations
Act goals, assess whether your governance architecture would achieve better, similar, or
worse outcomes than Wales has achieved between 2015 and 2026. Explain any gaps.

---

*Part VIII is now complete. Five chapters have developed the full transition architecture of the
cooperative-regenerative economy: the theoretical framework for complex system change
(Chapter 40), the policy instrument portfolio and its sequencing (Chapter 41), the
predistributive foundations in labour, care, and working time (Chapter 41b), the democratic
governance architecture that steers the transition (Chapter 41c), and the experimentation
framework that generates the empirical learning needed to refine all of the above (Chapter 42).
Part IX turns to what remains: the open questions that define the research frontier, the capstone
design project, and the conclusion that situates the cooperative-regenerative framework within
the longer arc of economic thought.*
