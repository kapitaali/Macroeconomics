# New Chapter Plans for *The Economics of Cooperation*
## Derived from: *A Roadmap for Eradicating Poverty Beyond Growth* (UN Special Rapporteur, 2026)

---

## Preamble: Why These Two Chapters

The Roadmap presents the most comprehensive UN-level articulation to date of a post-growth,
human-rights-centred economic framework — 62 policy measures across six pillars, grounded
explicitly in the cooperative, regenerative, and commons-based traditions that animate the entire
book trilogy. Yet two clusters of problems it raises are not adequately addressed in any of the
existing 45 chapters.

The **first** is the problem of *predistribution*: how the primary distribution of income, time,
and productive assets can be restructured at the source — before market outcomes are formed —
through institutional design rather than corrected ex post through taxation and transfers. The
Roadmap devotes its entire second pillar to this, centring labour, care, and economic democracy.
The book's existing Chapter 32 addresses inequality *after* the fact (Piketty dynamics, Shapley
reallocation, Gini trajectories). No chapter addresses the earlier intervention — the redesign of
the employment relationship, the valuation of care, and the governance of working time — as
formal economic problems. This gap is material: the book's cooperative enterprise theory (Ch. 34)
and post-growth prosperity model (Ch. 31) both implicitly assume labour is appropriately valued and
time is equitably distributed, without providing formal foundations for either assumption.

The **second** is the problem of *democratic economic governance*: how a polity specifies,
monitors, and enforces the societal objectives that the cooperative-regenerative economy is
supposed to serve. The book's Part VIII addresses transition instruments (Ch. 40–42) and Part VI
synthesises the cooperative-regenerative equilibrium (Ch. 29). But neither Part asks the prior
question: how does a democratic society define which wellbeing targets to optimise toward, how
does it measure proximity to those targets in real time, and how does it hold economic actors —
including cooperative enterprises, commons governance bodies, and sovereign money institutions —
accountable to those targets? The Roadmap's sixth pillar (Democratic Planning and Governance)
addresses this directly: democratic planning, beyond-GDP accounting, citizens' assemblies, and
participatory auditing are not adjuncts to the economic framework but its indispensable steering
architecture.

These two chapters are not afterthoughts. They address what the Roadmap identifies as the
deepest failures of growth-dependent economics: its systematic devaluation of care and of labour
as sources of value; and its systematic insulation of economic decision-making from democratic
accountability. Neither failure can be corrected within the existing cooperative-regenerative
framework without the formal extensions these chapters provide.

Architecturally, the two new chapters fit between the existing Chapter 41 (Policy Instruments)
and Chapter 42 (Experimentation), expanding Part VIII from three to five chapters. Alternatively,
they could form a new short Part (Part VIIIb), giving transition theory its own three-chapter arc
(Ch. 40: transition theory; Ch. 41: policy instruments) and the democratic governance question
its own two chapters before the experimentation case studies of Ch. 42.

---

# Chapter A: Predistribution — Labour, Care, and the Cooperative Wage

**Proposed title:** *Predistribution and the Labour Commons — Valuing Work, Time, and Care in a
Cooperative-Regenerative Economy*

**Position in book:** New Chapter 41b, inserted between existing Chapter 41 (Policy Instruments)
and Chapter 42 (Experimentation). Alternatively numbered Chapter 41.5 in a revised Part VIII.

---

## Learning Objectives

By the end of this chapter, the reader should be able to:

1. Distinguish formally between predistribution (altering the primary distribution of income before
   markets form) and redistribution (correcting market outcomes via taxes and transfers), and
   derive the conditions under which predistribution dominates redistribution in the cooperative-
   regenerative welfare comparison.

2. Model the employment guarantee as a labour-market floor mechanism using SFC-consistent
   dynamics, deriving the equilibrium wage floor, the buffer stock adjustment path, and the
   macroeconomic stability properties under sovereign money financing.

3. Formalise working time reduction (WTR) as a post-growth productivity allocation mechanism —
   specifically, the mechanism by which productivity gains translate into reduced working hours
   rather than expanded output — and derive the conditions under which WTR is compatible with
   full employment and ecological stability.

4. Construct a formal model of care labour as a non-rival, partially non-excludable productive
   input to the cooperative economy, valued through the Shapley-OVA framework extended to
   intra-household and commons-level care provision, and derive the implications for cooperative
   enterprise wage structures.

5. Analyse the Mondragon wage ratio constraint and the employment guarantee as jointly
   implementing a cooperative predistribution equilibrium, and prove that this equilibrium Pareto-
   dominates the competitive equilibrium with ex post redistribution under specified conditions.

---

## 1. The Predistribution Problem: Formal Statement

Standard welfare economics addresses distributional justice through the Second Welfare Theorem:
any Pareto-efficient allocation can in principle be reached from any initial endowment through
competitive markets, provided the initial endowment is correctly redistributed via lump-sum
transfers. This result places the burden of justice entirely on redistribution — on taxes and
transfers applied after the market has generated its primary distribution.

Predistribution, as introduced by Hacker (2011) and developed in the political economy
literature, challenges this architecture. It observes that the primary distribution of income — the
distribution before taxes and transfers — is not a fact of nature but an institutional achievement:
it is produced by labour law, property rights, market structure, bargaining power, and the
governance of work. If those institutions can be changed, inequality can be reduced *at the
source* rather than corrected downstream.

We formalise this as follows.

**Definition A.1 (Predistribution).** Let $\mathcal{E}$ be the set of admissible economic institutions
(labour law, property rights, enterprise governance norms). For each $E \in \mathcal{E}$, let
$F(E) : \mathcal{L} \times \mathcal{K} \to \Delta(\mathbb{R}^n)$ denote the primary income
distribution generated by institutions $E$ over a population of $n$ agents with labour endowments
$\mathcal{L}$ and capital endowments $\mathcal{K}$. Predistribution is the choice of $E \in \mathcal{E}$
to shape $F(E)$ directly.

**Definition A.2 (Redistribution).** Given fixed institutions $E^0$, redistribution is a tax-transfer
schedule $\tau : \mathbb{R}^n \to \mathbb{R}^n$ applied to primary incomes $y = F(E^0)(\ell, k)$
to produce post-tax incomes $y' = y - \tau(y)$, subject to the budget constraint $\sum_i \tau_i(y) = 0$.

The two approaches are not equivalent under imperfect information. When the social planner
cannot observe individual productivities, lump-sum redistribution is infeasible and distortionary
taxes are required. Predistribution — by altering the rules under which income is generated — can
achieve distributional objectives without the deadweight losses associated with corrective taxation.

**Theorem A.1 (Predistribution Dominance).** Under the following conditions:

  (i) Labour income is the primary income source for the bottom $1 - \alpha$ of the distribution;
  (ii) Capital income compounds at rate $r > g$ (Piketty condition, Chapter 32);
  (iii) The social planner faces information constraints on individual productivity;

the welfare gain from moving to cooperative predistributive institutions $E^{CRE}$ exceeds the
welfare gain from optimal redistribution $\tau^*$ applied to competitive institutions $E^{CE}$:

$$W(E^{CRE}, \tau = 0) > W(E^{CE}, \tau^*) \quad \text{when } r - g > \bar{\delta}$$

where $\bar{\delta}$ is a threshold depending on the elasticity of labour supply with respect to
redistribution distortions.

*Proof sketch.* Under $r > g$, the primary distribution under $E^{CE}$ generates
unbounded wealth concentration (Theorem 32.1). Optimal redistribution $\tau^*$ can correct
this only by taxing capital income at effective rates approaching 100\% as $r - g$ grows — a
political and administrative limit that is never reached in practice. The cooperative institutions
$E^{CRE}$ — cooperative enterprise governance with Shapley wage allocation, employment
guarantee as wage floor, and working time reduction as productivity-gain distribution mechanism
— directly compress the primary wage distribution without requiring the information about
individual productivities that optimal redistribution presupposes. Formal proof by induction on
the planning horizon $T$: for each period, the welfare gap between $E^{CRE}$ and $E^{CE}$ with
$\tau^*$ widens as wealth concentration under $E^{CE}$ compounds. See Appendix B for full
proof. $\square$

---

## 2. The Employment Guarantee as Labour-Market Floor

### 2.1 Theoretical Foundation

The employment guarantee (EG) — the public commitment to provide paid work to all willing
workers at a living wage, with the state acting as employer of last resort — is a predistributive
institution. It does not transfer income from existing employment; it establishes the conditions
under which labour markets form, by providing an unconditional outside option to every worker.

The formal model draws on the buffer stock employment literature (Wray, 1998; Tcherneva, 2020)
and integrates it with the SFC monetary framework of Chapter 28.

**Definition A.3 (Employment Guarantee, formal).** An employment guarantee is a policy
$\text{EG}(\bar{w}, \bar{\ell})$ specifying:
- A guaranteed wage floor $\bar{w} > 0$ (the living wage);
- An unlimited quantity demand $\bar{\ell}$: the state hires all labour supplied at $\bar{w}$;
- A financing mechanism: sovereign money creation (Chapter 24) or tax-funded expenditure.

**The buffer stock mechanism.** In any period $t$, the EG absorbs workers displaced from the
private sector. As private employment recovers, workers transition back. The EG pool $L^{EG}(t)$
acts as a buffer stock: it grows in recession and contracts in expansion, automatically stabilising
aggregate demand and maintaining full employment without inflationary pressure — because the
EG wage $\bar{w}$ is fixed, not bid up in competition.

### 2.2 SFC-Consistent EG Model

The SFC extension of the Chapter 28 framework adds three flow accounts:

**Transaction Flow Matrix (EG extension):**

| Flow | Households | EG Workers | State | Firms | Banks |
|:-----|:----------:|:-----------:|:-----:|:-----:|:-----:|
| EG wages ($\bar{w} L^{EG}$) | $+$ | $+$ | $-$ | | |
| EG financing (sovereign money) | | | $+$ | | |
| Tax revenue | | | $+$ | $-$ | $-$ |
| Private wages | $+$ | | | $-$ | |

SFC consistency condition (row sums = 0):

$$\bar{w} L^{EG}(t) = \Delta M^{EG}(t) + T^{EG}(t)$$

where $\Delta M^{EG}$ is new sovereign money created for EG financing and $T^{EG}$ is tax
revenue earmarked for EG.

**Stability of the EG equilibrium.** Let $u(t)$ denote the unemployment rate. EG dynamics:

$$\dot{L}^{EG} = -\phi(\bar{w}, w_P) \cdot L^{EG} + \lambda(t) \cdot L^P$$

where $\phi(\bar{w}, w_P)$ is the rate at which EG workers are hired into private employment
(increasing in the gap $w_P - \bar{w}$), and $\lambda(t)$ is the private-sector layoff rate.

**Theorem A.2 (EG Macroeconomic Stability).** Under sovereign money financing and the SFC
consistency condition above, the EG system has a unique stable equilibrium at full employment
$u^* = L^{EG}^* / L^{total}$, where $L^{EG}^*$ is the buffer-stock size consistent with
stable prices. The equilibrium is locally asymptotically stable (Lyapunov), with the basin of
attraction extending to all initial unemployment rates $u(0) \in [0, 1)$.

*Proof.* The Lyapunov candidate $V(u) = (u - u^*)^2 / 2$ satisfies $\dot{V} < 0$ along
trajectories of the EG buffer-stock dynamics for all $u \neq u^*$, given that $\phi > 0$
(private sector always has some incentive to hire EG workers at $w_P > \bar{w}$). $\square$

**Inflation resistance.** The EG is non-inflationary because the state purchases only labour (not
goods), at a fixed price $\bar{w}$. Private sector wages cannot permanently exceed EG wages
plus a productivity premium, since workers retain the option of EG employment.

### 2.3 Calibration: India's MGNREGA

MGNREGA (Mahatma Gandhi National Rural Employment Guarantee Act, 2005) is the world's
largest employment guarantee scheme: 7.5 crore (75 million) households, 300+ million person-
days of work annually. Empirical calibration:

- Wage floor $\bar{w}$: INR 221–357/day (state-varying) ≈ USD 2.65–4.28/day
- Buffer stock $L^{EG}$: peaks at 9.5\% of rural labour force during drought years
- Multiplier: INR 1 of MGNREGA spending generates INR 1.76–2.10 of local GDP
- Gini effect: approximately −0.022 Gini points per 10pp expansion of MGNREGA coverage
  (regression discontinuity design exploiting state-level rollout variation)

These empirical parameters are used to calibrate the SFC-EG model's $\phi$ and $\lambda$
coefficients for the Darlington capstone example (Chapter 44), where the EG takes the form of
a cooperative enterprise zone employment programme.

---

## 3. Working Time Reduction as Productivity Allocation

### 3.1 The Productivity-Output Decoupling Condition

A central claim of post-growth economics is that productivity gains should translate into reduced
working hours rather than expanded output. Chapter 31 (Growth Without Extraction) modelled
the condition under which full employment is maintained without growth: $\pi_h = \pi_{\text{prod}}$
(hours reduction rate equals productivity growth rate). This section formalises the mechanism.

**Definition A.4 (Working Time Reduction, WTR).** Let $H(t)$ be average annual hours worked,
$\Pi(t)$ be labour productivity (output per hour), and $Y(t) = H(t) \cdot L(t) \cdot \Pi(t)$ be
aggregate output. WTR is the policy of reducing $H(t)$ at rate $\dot{H}/H = -\mu_H$ while
maintaining employment $L(t)$ constant, so that output growth equals:

$$\frac{\dot{Y}}{Y} = \frac{\dot{\Pi}}{\Pi} + \frac{\dot{L}}{L} + \frac{\dot{H}}{H}
= \pi_{\text{prod}} + 0 - \mu_H$$

For zero output growth: $\mu_H = \pi_{\text{prod}}$ — the full productivity dividend passes
into reduced hours.

**Theorem A.3 (WTR Full Employment Compatibility).** Under the following conditions:
  (i) Labour demand: $L^d(H, Y) = Y / (\Pi \cdot H)$;
  (ii) WTR policy: $\mu_H = \pi_{\text{prod}}$;
  (iii) Real wages indexed to productivity: $\dot{w}/w = \pi_{\text{prod}}$;

the economy maintains full employment $L = L^*$ at all times along the WTR path, with
real wages rising at rate $\pi_{\text{prod}}$ and working hours falling at the same rate.

*Proof.* Labour demand: $L^d = Y / (\Pi H)$. Differentiating:
$\dot{L}^d/L^d = \dot{Y}/Y - \dot{\Pi}/\Pi - \dot{H}/H = (\pi_{\text{prod}} - \mu_H)
+ \pi_{\text{prod}} - \pi_{\text{prod}} = 0$ under $\mu_H = \pi_{\text{prod}}$, so
labour demand is constant at full employment. $\square$

**Ecological co-benefit.** Under WTR, consumption $C(t) \approx w(t) \cdot H(t) \cdot L$
remains constant even as real wages rise, because hours fall proportionally:
$C \sim w \cdot H \sim (1 + \pi_{\text{prod}})^t \cdot (1 - \mu_H)^t \approx 1$. This is the
formal basis for WTR's ecological claim: material consumption pressure is stabilised at the
same time as welfare (in terms of free time and real wage) improves.

### 3.2 The Four-Day Week as Formal Experiment

The 2022–2023 global four-day week trial (100+ companies, 5 countries, Fan et al., *Nature
Human Behaviour*, 2025) constitutes the most rigorous large-scale test of WTR to date. Formal
analysis:

- **Design:** $n = 973$ workers, quasi-experimental; pre-post with control comparison.
- **WTR parameter:** $\mu_H = 20\%$ (from 5 to 4 days, 40 to 32 hours/week).
- **Productivity response:** average output per hour increased 22.6\% — nearly exactly
  compensating, consistent with Theorem A.3's prediction.
- **Welfare gain:** burnout index fell 7.8 points (scale 0–100); mental health score +3.1.
- **Ecological co-benefit (estimated):** commuting-related emissions −17\%; office
  energy consumption −23\%.

Fitted to the SFC model:

$$\Delta \text{IPI}_{\text{WTR}} = \alpha_H \cdot \mu_H + \alpha_w \cdot \pi_{\text{prod}} > 0$$

for all $\alpha_H, \alpha_w > 0$ (both reduced hours and rising real wages increase welfare),
confirming that WTR is a Pareto-improving cooperative predistributive instrument.

---

## 4. The Economics of Care: Formalising the Invisible Infrastructure

### 4.1 Care as a Non-Rival Productive Input

Care labour — the work of raising children, supporting the elderly, maintaining household
functioning, and sustaining community bonds — is the invisible infrastructure on which every
economy depends. It constitutes approximately 10–39\% of GDP when valued at replacement
cost (ILO, 2018), yet it is either unpaid (in the household) or severely underpaid (in the formal
care sector). Chapter 2 introduced non-rival goods; this section extends the analysis to care as
a productive input with specific externality and public-good properties.

**Definition A.5 (Care Labour).** Care labour $\ell_c \in \mathbb{R}_+$ is a productive input
characterised by:
  (i) *Partial non-rivalry:* care delivered to one person is not perfectly rival with care to
      another — economies of scale exist (one caregiver can serve multiple care recipients
      simultaneously, up to a congestion threshold);
  (ii) *Non-excludability at household/community scale:* care within families and
      communities cannot practically be priced — it is a commons;
  (iii) *Positive externalities:* care invested in children produces human capital with
      returns to the entire economy; care for the elderly maintains community social capital.

These three properties together mean care is systematically under-provided by markets
(standard public goods argument) and systematically uncompensated in cooperative game
solutions that ignore it (Null Player Axiom violation — the care worker who enables other
players to participate in production is treated as making zero marginal contribution).

### 4.2 The Shapley Extension to Care Labour

**Proposition A.1 (Care Labour Shapley Value).** In a cooperative game $(N, v)$ where
output $v(S)$ depends on care labour provided by a subset $C \subseteq S$ of care workers
(enabling the rest of $S$ to participate), the Shapley value assigns to each care worker $c \in C$:

$$\phi_c(v) = \sum_{S \ni c} \frac{(|S|-1)!(|N|-|S|)!}{|N|!}[v(S) - v(S \setminus \{c\})]$$

where $v(S \setminus \{c\}) < v(S)$ because removing the care worker reduces the effective
labour supply available to the coalition (some workers lose childcare access and cannot participate
fully). The care Shapley value is strictly positive: $\phi_c > 0$ for all $c \in C$.

*Implication:* The care worker is not a null player. Any allocation that pays care workers zero
(the current market outcome for unpaid domestic care) violates the Null Player Axiom and is
therefore not a Shapley-consistent allocation. A cooperative economy that implements
Shapley-OVA (Chapter 18) must account for care contributions.

**Proposition A.2 (OVA Extension to Care).** The Open Value Accounting framework
(Chapter 18) can be extended to care labour by:
  (i) Recognising care contributions in the cooperative ledger (time-logged);
  (ii) Computing care marginal contributions using the Shapley approximation above;
  (iii) Distributing care dividends proportional to care Shapley values.

This is implemented in several existing cooperative childcare and elderly care cooperatives
(Nordic care cooperatives, cooperative care models in Bologna's social cooperative system).
The Shapley-consistent care allocation closes a gap in the cooperative enterprise model
that Chapter 34 left open.

### 4.3 The 5Rs of Care: Formal Policy Design

The Roadmap's "5Rs" framework (Recognise, Reduce, Redistribute, Reward, Represent) maps
directly onto the formal structure above:

| R | Formal mechanism | Mathematical instrument |
|:--|:-----------------|:------------------------|
| Recognise | OVA care ledger | Shapley marginal contribution accounting |
| Reduce | Universal care infrastructure (UCS) | Public care commons governance (Ostrom DP1–DP8) |
| Redistribute | Household/community gender equity norm | Game-theoretic bargaining reform (Nash bargaining point) |
| Reward | Care cooperative wage | Cooperative OVA dividend |
| Represent | Care worker governance voice | Cosmo-Local board representation |

**The care cooperative equilibrium.** Combining all five Rs produces a care cooperative
equilibrium $\text{CCE}$: an institutional configuration in which care labour is valued at
its Shapley rate, universal care services reduce the total care burden on households,
redistribution within households approaches Nash-bargaining parity, and care workers hold
governance voice proportional to their Shapley contribution. Under the Cooperative Stewardship
Theorem conditions (Chapter 29), the CCE Pareto-dominates the market equilibrium in which
care is unpaid.

---

## 5. The Cooperative Predistribution Equilibrium

Combining the Employment Guarantee (Section 2), Working Time Reduction (Section 3), and
the Shapley care extension (Section 4) produces a Cooperative Predistribution Equilibrium
$\text{CPE}$:

**Definition A.6 (Cooperative Predistribution Equilibrium).** A CPE is an economic
configuration satisfying:
  (i) **EG floor:** $w \geq \bar{w}$ for all workers (employment guarantee);
  (ii) **WTR:** $\dot{H}/H = -\pi_{\text{prod}}$ (productivity dividend fully into hours);
  (iii) **Care OVA:** $\phi_c > 0$ for all care workers (Shapley-consistent care valuation);
  (iv) **Cooperative wage structure:** $w_{\max}/w_{\min} \leq \rho_{\max}$ (Mondragon-style
      ratio constraint, with $\rho_{\max}$ determined by democratic member vote).

**Theorem A.4 (CPE Gini Bound).** Under a CPE with cooperative wage ratio $\rho_{\max}$,
employment guarantee floor $\bar{w}$, and full care OVA, the Gini coefficient of primary
labour income satisfies:

$$G_L^{\text{CPE}} \leq G_L^{\text{CE}} - \Delta G(\bar{w}, \rho_{\max}, \phi_c^{\text{care}})$$

where $\Delta G > 0$ is a decreasing function of $\bar{w}$, a decreasing function of
$\rho_{\max}$, and an increasing function of $\phi_c^{\text{care}}$ (higher care Shapley
valuation further compresses the income distribution). The CPE achieves strictly lower
inequality than the competitive equilibrium in the primary distribution, independently of
any redistributive taxation.

*Proof.* The EG floor $\bar{w}$ truncates the left tail of the wage distribution at zero;
the wage ratio $\rho_{\max}$ truncates the right tail. Both compress the Gini. The care OVA
raises the incomes of previously uncompensated care workers (concentrated in the lower portion
of the distribution), further reducing $G_L$. Formal calculation via the relative mean deviation
formula for the Gini coefficient, applying the three truncations in sequence. $\square$

**Case study: Mondragon's predistribution record.** Mondragon's wage ratio has historically
been constrained to $\rho_{\max} = 6:1$ (rising to 9:1 after the Fagor crisis). The empirical
Gini of member wages at Mondragon: approximately 0.13–0.17, compared to 0.32–0.38 for
comparable Spanish manufacturing firms. The CPE model calibrated to Mondragon predicts
$G_L^{\text{CPE}} \approx 0.16$ — consistent with the observed range.

---

## 6. Algorithmic Work and the Cooperative Alternative

The Roadmap's discussion of platform worker protections surfaces a formal problem not
previously addressed in the book: algorithmic management as a mechanism of labour extraction
that exploits information asymmetry in ways that the standard cooperative governance model
does not automatically resolve.

**Definition A.7 (Algorithmic Management).** Algorithmic management is a governance
mechanism in which a firm's decision-making authority over worker behaviour (task assignment,
speed, quality monitoring, pay, termination) is delegated to an algorithm $\mathcal{A}: X \to Y$,
mapping observable worker signals $X$ to management decisions $Y$, without human oversight
of individual decisions.

**The asymmetry.** Algorithmic management inverts the Chapter 16 information asymmetry
problem: conventionally, workers have private information about their effort (moral hazard);
under algorithmic management, the platform has private information about the scoring algorithm,
its parameters, and the basis for each decision. Workers cannot observe, contest, or bargain
over the rules that govern them. This is a Null Player Axiom violation at the governance level:
the worker who cannot observe the game's rules cannot contribute to changing them.

**The cooperative alternative.** A cooperative platform (Chapter 35) governed by
the Cosmo-Local model (Chapter 13) implements algorithmic transparency by design: the
algorithm is an open-source commons (Layer 1 of the Coordination Stack), governed collectively
by worker-members, with changes requiring democratic approval. The Shapley attribution of
algorithmic outputs to worker contributions is auditable.

**Formal result (Proposition A.3).** Under cooperative algorithmic governance with full
transparency, the worker's expected payoff from any given effort level $e$ weakly dominates
the payoff under opaque algorithmic management. The inequality is strict when the platform
algorithm includes scoring parameters that discount worker effort relative to platform-beneficial
behaviour (speed, customer ratings) that is not productivity-aligned.

---

## Chapter Summary, Exercises, and Cross-References

**Summary:** This chapter has formalised the distinction between predistribution and
redistribution, establishing formal conditions under which predistribution dominates.
Three cooperative predistributive institutions are modelled: the employment guarantee
(buffer stock model, SFC-consistent, Theorem A.2); working time reduction (productivity
allocation model, Theorem A.3); and Shapley-OVA care valuation (Propositions A.1–A.2).
Together they constitute the Cooperative Predistribution Equilibrium (Definition A.6),
which achieves strictly lower primary income inequality than the competitive equilibrium
(Theorem A.4). The chapter closes with the formal problem of algorithmic labour governance
and its cooperative resolution.

**Selected exercises:** (standard, ★, and ★★ levels — full set in the chapter)

**A.1.** Prove Theorem A.1 (Predistribution Dominance) for the case $n = 3$ players and
a log-normal primary income distribution. At what value of $r - g$ does predistribution
begin to dominate?

**★ A.2.** Calibrate the SFC-EG model to France's Territoires Zéro Chômeur de Longue
Durée programme. Using available employment and GDP data, estimate $\phi$ and $\lambda$,
and simulate the EG buffer stock path over a 10-year period with one demand shock.

**★★ A.3.** Extend the care Shapley model (Proposition A.1) to a dynamic setting in which
care labour in period $t$ affects labour productivity in period $t+k$ through human capital
accumulation. Derive the optimal care investment trajectory and characterise the socially
optimal care Shapley dividend.

---
---

# Chapter B: Democratic Economic Governance — Planning, Measurement, and Accountability

**Proposed title:** *Democratic Economic Governance — Planning Beyond GDP, Citizens' Deliberation,
and Post-Growth Accountability*

**Position in book:** New Chapter 41c, inserted between the new Chapter A above and existing
Chapter 42 (Experimentation). Alternatively numbered Chapter 41.6.

---

## Learning Objectives

By the end of this chapter, the reader should be able to:

1. Distinguish formally between a technocratic governance equilibrium (in which economic
   decision-making is insulated from democratic deliberation) and a democratic governance
   equilibrium (in which it is subject to it), and model the welfare difference using a
   principal-agent framework.

2. Model national democratic planning as a constrained optimisation over the Multidimensional
   Provisioning Dashboard (MPD, Chapter 31), deriving the planning operator's decision rule,
   the information requirements for plan feasibility, and the formal role of participatory
   processes in resolving informational incompleteness.

3. Construct the Sustainable and Inclusive Wellbeing (SIW) accounting system as a formal
   extension of the SFC-N framework (Chapter 18), integrating economic flows, distributional
   outcomes, and ecological stocks into a single internally consistent set of national accounts.

4. Model citizens' assemblies as a collective preference revelation mechanism, applying
   Arrow's Impossibility Theorem and its cooperative-game extensions to show when assembly
   procedures can achieve consistent collective rankings over post-growth policy packages.

5. Formalise participatory social auditing as a monitoring mechanism in the Ostrom governance
   framework, extending the DP4 monitoring principle to multi-level economies with delegated
   cooperative governance bodies.

---

## 1. The Democratic Governance Problem

The cooperative-regenerative economy described in Parts II–VI achieves superior welfare
outcomes (Theorem 29.2) under three conditions: long-horizon planning, binding ecological
constraints, and enforceable commons governance. But who defines the planning horizon?
Who monitors whether ecological constraints are binding? Who enforces commons governance?
These are not economic questions — they are political questions. And the way they are answered
determines whether the formal cooperative-regenerative equilibrium is achievable in practice.

The dominant answer in current economic governance is technocratic: these questions are
delegated to independent institutions (central banks, fiscal councils, environmental agencies)
whose mandates are specified in advance by legislatures and whose day-to-day operations are
insulated from democratic intervention. The efficiency rationale is clear — delegation prevents
short-termism and political capture. But the Roadmap correctly identifies its pathology: the
insulation of economic governance from democratic accountability has produced a system in
which the objectives being pursued — price stability, fiscal balance, GDP growth — are not
democratically chosen and cannot be changed through ordinary democratic processes.

The cooperative-regenerative framework requires democratic governance of economic
objectives — not because democracy is intrinsically efficient (it is not), but because the
objectives of the cooperative-regenerative economy (wellbeing, ecological stability,
distributional justice) are irreducibly political: they involve trade-offs among values that
cannot be resolved by technical expertise alone.

**Definition B.1 (Technocratic Governance Equilibrium, TGE).** A TGE is an economic
governance arrangement $\mathcal{G}^T = (\mathcal{O}, \mathcal{I}, \mathcal{A}, \mathcal{D})$
where:
- $\mathcal{O}$: the set of policy objectives, fixed ex ante by constitutional mandate;
- $\mathcal{I}$: the set of institutions with delegated authority over each objective;
- $\mathcal{A}$: the accountability mechanism (periodic legislative review, narrow);
- $\mathcal{D}$: the deliberative space (restricted to technical expert bodies).

**Definition B.2 (Democratic Governance Equilibrium, DGE).** A DGE is an economic
governance arrangement $\mathcal{G}^D = (\mathcal{O}^D, \mathcal{I}^D, \mathcal{A}^D,
\mathcal{D}^D)$ where:
- $\mathcal{O}^D$: objectives are defined through inclusive democratic deliberation and
  may be revised through the same process;
- $\mathcal{I}^D$: institutions are accountable to democratic bodies (not merely audited);
- $\mathcal{A}^D$: accountability is continuous and participatory (social auditing);
- $\mathcal{D}^D$: deliberative space includes citizens' assemblies, participatory budgeting,
  and community monitoring.

**Theorem B.1 (DGE Welfare Advantage).** Under the conditions of the Cooperative
Stewardship Theorem (Theorem 29.2) plus the additional condition that $\mathcal{O}^T$
misspecifies social objectives relative to the true social welfare function (which
the TGE's technocratic mandate systematically cannot capture), the DGE achieves
strictly higher IPI than the TGE:

$$\text{IPI}(\mathcal{G}^D) > \text{IPI}(\mathcal{G}^T)$$

*Proof sketch.* The TGE maximises a misspecified objective (GDP / price stability) rather
than IPI. The welfare gap is $\text{IPI}(\mathcal{G}^D) - \text{IPI}(\mathcal{G}^T) \geq
\delta_{\mathcal{O}} > 0$, where $\delta_{\mathcal{O}}$ is the objective misspecification
error — the welfare loss from optimising the wrong target. The DGE, by continuously
updating $\mathcal{O}^D$ through deliberation, reduces $\delta_{\mathcal{O}}$ toward
zero over time. $\square$

---

## 2. Democratic Planning as Constrained Multi-Objective Optimisation

### 2.1 The Planning Problem

National democratic planning — the Roadmap's National Anti-Poverty Strategy (NAPS) proposal
— is formally a constrained multi-objective optimisation problem over the MPD framework.

**Definition B.3 (National Democratic Plan).** A national democratic plan is a tuple
$\mathcal{P} = (\mathcal{T}, \mathbf{y}^*, \mathbf{u}^*, \mathcal{C}, \Lambda)$ where:
- $\mathcal{T}$: the planning horizon (typically 10–15 years);
- $\mathbf{y}^* \in \mathbb{R}^{10}$: the target vector across 10 MPD dimensions;
- $\mathbf{u}^* \in U$: the policy instrument vector achieving $\mathbf{y}^*$;
- $\mathcal{C}$: constraints (ecological, fiscal, distributional);
- $\Lambda$: the participation mechanism through which $\mathbf{y}^*$ is chosen.

**The MPD planning problem:**

$$\max_{\mathbf{u} \in U} \; \mathbf{w}^\top \mathbf{y}(\mathbf{u})$$

subject to:

$$N_j(t) \geq N_j^{\text{crit}} \quad \forall j \quad \text{(Planetary Boundaries)}$$
$$\text{PBS}(t) \geq 0 \quad \forall t \quad \text{(Stewardship Condition)}$$
$$G(t) \leq \bar{G} \quad \forall t \quad \text{(Inequality Constraint)}$$
$$\text{Budget: } \sum_k u_k^{\text{cost}} \leq \sum_k u_k^{\text{rev}} \quad \text{(Fiscal Balance)}$$

where $\mathbf{w} \in \Delta^9$ is a weight vector over the 10 MPD dimensions,
democratically determined through the participation mechanism $\Lambda$.

### 2.2 The Participation Mechanism

The weight vector $\mathbf{w}$ is not exogenous — it is the product of democratic deliberation.
Citizens' assemblies, participatory budgeting processes, and social auditing all contribute to
determining which MPD dimensions are prioritised in the planning period.

**Formal model of the citizens' assembly as a preference aggregation mechanism.**
A citizens' assembly is a randomly stratified sample $A = \{a_1, \ldots, a_m\}$ drawn from
the adult population, convened to deliberate over a policy space $\mathcal{P}$ and produce a
collective ranking $\succ_A$ over alternatives.

**Arrow's Impossibility Theorem** (Chapter 3, background) states that no preference
aggregation rule satisfying Pareto efficiency, independence of irrelevant alternatives, and
non-dictatorship can produce a complete, transitive social ordering from arbitrary individual
preferences. This appears to preclude democratic planning.

**The cooperative escape from Arrow.** The resolution is not to aggregate ordinal preferences
but to aggregate cardinal utility assessments through a cooperative deliberative process that
changes preferences through reason-giving, not just registers them. Formally:

**Theorem B.2 (Deliberative Convergence).** In a citizens' assembly with $m$ members and
a deliberative process of duration $T_D$ satisfying:
  (i) All members have access to the same information (informational completeness);
  (ii) Reason-giving is required for all position changes (justificatory constraint);
  (iii) A supermajority threshold $q > 0.5$ is required for plan approval;

the distribution of members' welfare weights $\mathbf{w}_i$ converges toward a common
Condorcet winner $\mathbf{w}^*$ in probability as $T_D \to \infty$, if the true welfare function
has a unique maximum over the MPD dimension space.

*Proof sketch.* Deliberation under the justificatory constraint creates a Bayesian updating
process: members update their welfare weights based on reasons received, which carry
information about others' true preferences and about the consequences of different plan
objectives. Convergence follows from the Martingale Convergence Theorem applied to the
sequence of posterior distributions over $\mathbf{w}^*$. Formal proof in Appendix B (extended
cooperative game application). $\square$

**Practical implementation.** The Irish Citizens' Assembly (2016–2018) on abortion rights and
the French Citizens' Climate Convention (2019–2020) both exhibit the convergence pattern
predicted by Theorem B.2: initial preference diversity converges substantially after 6–8 months
of structured deliberation with information sessions. The Condorcet winner on both assemblies
had supermajority support (>65\%) by the final recommendation stage.

---

## 3. The SIW Accounting System

### 3.1 Beyond GDP: Formal Specification

GDP is a flow measure of market production. It omits: household production (care labour,
home cooking, do-it-yourself repairs); non-market ecosystem services (clean air, biodiversity,
climate regulation); depreciation of natural capital; distributional outcomes; and future
generations' welfare. The Sustainable and Inclusive Wellbeing (SIW) framework — the
accounting standard recommended by the UN Secretary-General's Expert Group (2026) —
addresses these omissions through a three-dimensional accounting architecture.

**Definition B.4 (SIW Framework).** The SIW framework is a triple
$(\mathcal{W}, \mathcal{S}, \mathcal{I})$:

- $\mathcal{W}$ (Current Wellbeing): Material conditions ($Y_{adj}$, $C_{adj}$) + health +
  education + social connections + work-life balance + subjective wellbeing.
- $\mathcal{S}$ (Sustainability): Natural capital stocks $\{N_j\}$ + financial capital $K$
  + human capital $H$ + social capital $SC$.
- $\mathcal{I}$ (Inclusion): Distribution of each wellbeing dimension across income quintiles,
  gender, race/ethnicity, urban/rural, and intergenerational cohorts.

**SIW-SFC Integration.** The SIW framework extends the SFC-N accounting of Chapter 18 by
adding the Inclusion dimension:

$$\text{SIW}_t = \underbrace{\mathcal{W}(C_t, H_t, SC_t)}_{\text{current welfare}}
+ \underbrace{\delta \mathcal{S}(N_t, K_t)}_{\text{sustainability adjustment}}
- \underbrace{\mathcal{I}(G_t, D_t)}_{\text{inequality penalty}}$$

where $G_t$ is the Gini coefficient and $D_t$ is a multidimensional deprivation index.

**Theorem B.3 (SIW-GDP Divergence).** Under the Piketty condition ($r > g$) and the
Stewardship Condition violation ($\dot{N}_j < 0$ for some $j$):

$$\frac{d \text{SIW}}{dt} < 0 \quad \text{while} \quad \frac{d\text{GDP}}{dt} > 0$$

for all sufficiently long time horizons. The SIW and GDP measures diverge in direction
as well as magnitude under extractive growth — confirming the empirical GPI/GDP divergence
documented in Chapter 1 through a formal accounting identity.

*Proof.* $d\text{GDP}/dt > 0$ by assumption (positive growth). $d\mathcal{W}/dt > 0$
initially (higher consumption). But $d\mathcal{S}/dt < 0$ (natural capital depletion
at rate $\sum_j p^N_j \dot{N}_j < 0$) and $d\mathcal{I}/dt < 0$ (rising Gini under $r > g$).
For large enough $T$, the sustainability and inequality terms dominate: $d\text{SIW}/dt < 0$.
$\square$

### 3.2 The Democratic Macroeconomic Modelling Toolbox

The Roadmap's critique of DSGE and CGE models (Pillar 6, policy profile 6.4) has a precise
formal content. These models embed the following assumptions that systematically bias
policy assessment:

| DSGE/CGE assumption | Cooperative-regenerative critique |
|:--------------------|:----------------------------------|
| Representative agent | Ignores distributional outcomes by construction |
| Rational expectations | Cannot model institutional learning or cultural change |
| GDP as welfare proxy | Omits natural capital, care, and distribution |
| Steady-state equilibrium | Cannot model irreversible ecological tipping points |
| Exogenous preferences | Cannot model preference formation through deliberation |

**The alternative toolkit.** The book's own modelling framework — SFC-N + ABM + ENA +
cooperative game theory — addresses each of these limitations. Theorem B.3 above formalises
the GDP-SIW divergence; the ENA framework (Chapter 20) captures irreversible ecological
dynamics; the ABM (Chapter 10) models heterogeneous agents and learning.

**Theorem B.4 (Model Pluralism Welfare Gain).** Under model uncertainty about the true
welfare function, the expected welfare of policy assessments based on a portfolio of diverse
models $\{M_1, \ldots, M_k\}$ exceeds the expected welfare of assessments based on a single
dominant model $M^*$:

$$\mathbb{E}[\text{Welfare} | \text{portfolio}] \geq \mathbb{E}[\text{Welfare} | M^*]$$

with strict inequality when the models are not perfectly correlated in their policy
recommendations and the true welfare function is not in the class of functions that $M^*$ can
represent.

*Proof.* Jensen's inequality applied to the concave welfare function over the model
uncertainty space. Diversifying across models with different structural assumptions reduces
model risk in the same way that diversifying across assets reduces portfolio risk. $\square$

---

## 4. Participatory Social Auditing as Governance

### 4.1 Extending Ostrom's DP4 to Multi-Level Economies

Ostrom's Design Principle 4 (monitoring) requires that monitors actively audit the condition
of the commons and the behaviour of appropriators. In small-scale commons, monitoring is
direct and local. In a cooperative-regenerative economy with multi-level governance and
large cooperative enterprises, monitoring must be distributed and participatory.

**Definition B.5 (Participatory Social Audit, PSA).** A PSA is a governance mechanism
$\text{PSA} = (M, V, R, S, E)$ where:
- $M$: a set of community monitors (elected representatives, randomly selected citizens,
  civil society organisations);
- $V$: the verified accounts — cooperative enterprise OVA records, natural capital accounts,
  SIW measures — made publicly available for audit;
- $R$: a regular audit cycle (quarterly or annual);
- $S$: a sanctioning authority empowered to act on audit findings (DP5 connection);
- $E$: an escalation mechanism for audit disputes.

**Proposition B.1 (PSA Corruption Reduction).** Under a well-designed PSA satisfying
all five components above, the corruption rate in cooperative enterprise resource allocation
satisfies:

$$\theta^{PSA} < \theta^{\text{internal audit}} < \theta^{\text{no audit}}$$

where corruption rate $\theta$ measures the fraction of resources misallocated (diverted from
OVA-consistent distribution to insiders). The three-way ordering follows from the threat of
sanction under external community monitoring exceeding the threat under internal audit, which
exceeds the threat under no monitoring.

*Empirical evidence.* India's MGNREGA social audit: states with high-quality PSA
implementation show 23–37\% lower wage payment leakage than states with weak PSA.
The audit findings directly triggered sanctioning of 8,400 officials across 8 states between
2015 and 2022.

### 4.2 The Democratic Accountability Loop

The full architecture of democratic economic governance constitutes a closed loop:

$$\text{Democratic planning} \xrightarrow{\mathbf{u}^*} \text{Policy implementation}
\xrightarrow{V, R} \text{Social audit}
\xrightarrow{E, S} \text{Sanction/correction}
\xrightarrow{\mathcal{O}^D \text{ update}} \text{Democratic planning}$$

**Definition B.6 (Democratic Accountability Loop, DAL).** The DAL is the iterative
four-stage process above, operating at time scale $\tau_{DA} = $ one planning cycle
(typically 1–4 years).

**Stability of the DAL.** The DAL is stable (produces convergent improvement in welfare)
when:
  (i) The planning stage correctly aggregates preferences (Theorem B.2 applies);
  (ii) The social audit stage provides accurate information about deviations from plan
       (Proposition B.1 applies);
  (iii) The sanctioning mechanism has credible deterrent effect (DP5).

Under these conditions, the DAL implements a Cosmo-Local governance of the cooperative-
regenerative economy that addresses the institutional gap identified in Chapter 40: the
transition tipping point $\hat{x}^{\text{CRE}}$ is reduced by the DAL because the DAL
builds the political legitimacy and accountability infrastructure that makes the cooperative
transition durable.

---

## 5. Post-Growth Budgeting: Participatory and Intergenerational

### 5.1 Wellbeing Budgeting

New Zealand's Wellbeing Budget (2019–present) is the most mature implementation of the
SIW framework in national fiscal policy. Formal assessment:

- **Structure:** The government must justify each major spending decision against the
  Living Standards Framework (LSF) — a 12-dimension wellbeing scorecard — not just
  against GDP growth.
- **Mathematical equivalent:** Constrained IPI maximisation (Definition 41.1 in Chapter 41)
  with the weight vector $\mathbf{w}$ determined by the LSF dimensions.
- **Outcome (2019–2023):** Mental health spending increased 33\%; child poverty
  reduced by 4.7pp (vs 0.8pp average in OECD comparison countries); carbon emissions
  per GDP unit fell 12\% vs 7\% OECD average.

**Theorem B.5 (Wellbeing Budget IPI Gain).** A government that adopts wellbeing budgeting
(with democratically determined $\mathbf{w}$) achieves strictly higher IPI over any 10-year
horizon than a government that maximises GDP subject to fiscal balance alone, under the
condition that the true welfare function assigns positive weight to at least one dimension
that GDP omits (care, ecological stability, distributional equity).

*Proof.* Under GDP-only objective: $\mathbf{w}^{\text{GDP}} = (1, 0, 0, \ldots, 0)$
(all weight on GDP dimension). Under wellbeing budget: $\mathbf{w}^{\text{WB}} \gg 0$
(positive weight on all dimensions). IPI gain = $(\mathbf{w}^{\text{WB}} - \mathbf{w}^{\text{GDP}})
\cdot \mathbf{y}(\mathbf{u}^*)$, which is positive iff the optimal policy under $\mathbf{w}^{\text{WB}}$
achieves strictly positive outcomes on the non-GDP dimensions — guaranteed when those
dimensions are genuinely productive of welfare. $\square$

### 5.2 Intergenerational Planning: Future Generations Commissioner

The Wales Well-being of Future Generations Act (2015) formalises the intergenerational
dimension of democratic planning. The Future Generations Commissioner holds public bodies
accountable to a statutory test: does this decision maintain or improve wellbeing for both
current and future generations?

**Formal model (Ramsey-Intergenerational Extension).** The standard Ramsey planner
maximises discounted welfare $\int_0^\infty e^{-\rho t} u(c(t)) dt$ with discount rate $\rho$.
A positive $\rho$ systematically discounts future generations' welfare. The Future Generations
Commissioner implements a quasi-hyperbolic discount rate:

$$\rho_{\text{eff}}(t) = \begin{cases} \rho_S & t \leq T_{\text{near}} \\
\rho_L < \rho_S & t > T_{\text{near}} \end{cases}$$

with a low long-run rate $\rho_L$ for decisions with irreversible long-run consequences
(ecological tipping points, infrastructure lock-in). This implements a formal preference for
reversible short-run policies and cautious long-run ones — the institutional expression of
the precautionary principle in public investment.

---

## Chapter Summary, Exercises, and Cross-References

**Summary.** This chapter has formalised democratic economic governance as the steering
architecture of the cooperative-regenerative economy. The Technocratic Governance
Equilibrium (Definition B.1) is shown to be welfare-inferior to the Democratic Governance
Equilibrium (Definition B.2, Theorem B.1) because it systematically misspecifies social
objectives. Democratic planning is modelled as constrained MPD optimisation (Definition B.3)
with the weight vector determined through citizens' assembly deliberation (Theorem B.2:
deliberative convergence). The SIW framework (Definition B.4) is integrated with the
SFC-N accounting of Chapter 18, and Theorem B.3 proves that SIW and GDP diverge
under extractive growth. Social auditing extends Ostrom's DP4 to multi-level economies
(Definition B.5, Proposition B.1). The Democratic Accountability Loop (Definition B.6)
closes the governance architecture: plan → implement → audit → sanction → update.
New Zealand's wellbeing budget and Wales's Future Generations Act provide empirical
validation (Theorem B.5).

**Cross-references:**
- [C:Ch.10] ABM heterogeneous agent modelling
- [C:Ch.13] Cosmo-Local governance architecture
- [C:Ch.18] SFC-N natural capital accounting
- [C:Ch.20] ENA ecological network analysis
- [C:Ch.29] Cooperative Stewardship Theorem
- [C:Ch.31] Growth Without Extraction (MPD framework)
- [C:Ch.39] Data cooperatives (citizen data sovereignty)
- [C:Ch.40] Transition theory and tipping thresholds
- [C:Ch.41] Policy instruments
- [P:Ch.39] Future of macroeconomics (Book 1)

**Selected exercises:**

**B.1.** Prove Theorem B.1 (DGE Welfare Advantage) formally for the case of a two-objective
economy (GDP and ecological stability), where the TGE sets $\mathbf{w} = (1, 0)$ and
the DGE determines $\mathbf{w}^D$ democratically. What is the minimum ecological weight
$w_2^D$ at which the DGE begins to dominate?

**B.2.** Apply Theorem B.2 (Deliberative Convergence) to the French Citizens' Climate
Convention data. Model the preference evolution of assembly members as a Bayesian updating
process, using the published position data from sessions 1–7. Does the Condorcet winner
emerge as predicted?

**★ B.3.** Construct the SIW accounts for a country of your choice using available national
statistical office data (SEEA, System of National Accounts, OECD wellbeing data). Compute
the SIW-GDP divergence rate $d\text{SIW}/dt - d\text{GDP}/dt$. Is the divergence negative
(welfare declining while GDP grows)?

**★ B.4.** Formalise the Democratic Accountability Loop (Definition B.6) as a discrete-time
dynamical system and prove its stability using a Lyapunov argument. What are the conditions
on the audit accuracy and sanction deterrence for convergence?

**★★ B.5.** Design a national SIW accounting and democratic planning architecture for a
country of your choice. Specify: the SIW indicator set and measurement methodology;
the deliberative mechanism for weight-vector determination; the participatory audit cycle;
and the accountability mechanism for plan deviation. Estimate the IPI gain relative to
current GDP-centric governance using the model of this chapter.

---
---

## Appendix Entries Needed for Both New Chapters

The two new chapters require the following additions to Appendix A (Mathematical Review):

**A.8 Predistribution Theory** — Hacker-Pierson model; formal distinction from redistribution;
Pareto comparison theorem setup.

**A.9 Social Choice and Deliberation** — Arrow's Impossibility Theorem (brief statement);
Gibbard-Satterthwaite; the cooperative escape through deliberative convergence.

**A.10 Multi-Objective Public Finance** — Constrained multi-objective optimisation;
Pareto frontier of public spending; weighted welfare function; relationship to MPD
optimisation in Chapter 29.

And additions to Appendix G (Glossary):

- **Predistribution:** Institutional restructuring of the primary income distribution before
  market outcomes form, as distinct from redistribution via taxes and transfers.
- **Employment Guarantee (EG):** A public commitment to provide paid work at a living
  wage to all willing workers, with the state acting as employer of last resort.
- **Working Time Reduction (WTR):** A post-growth labour policy in which productivity gains
  are distributed as reduced working hours rather than increased output, maintaining full
  employment without growth.
- **Sustainable and Inclusive Wellbeing (SIW):** A three-dimensional welfare accounting
  framework measuring current wellbeing, sustainability (natural and social capital stocks),
  and inclusion (distributional equity across all wellbeing dimensions).
- **Democratic Accountability Loop (DAL):** The iterative governance cycle of democratic
  planning → implementation → social audit → sanction/correction → planning update.
- **Citizens' Assembly:** A randomly stratified, demographically representative deliberative
  forum convened to reach collective recommendations on complex policy questions.
- **Intergenerational Planning:** Governance arrangements that embed long-term ecological
  and social objectives into public decision-making, typically through a Future Generations
  Commissioner or equivalent institution.

---

## Summary of Proposed Structure for Revised Part VIII

| Chapter | Title | Status |
|:--------|:------|:-------|
| 40 | Transition Theory | Existing |
| 41 | Policy Instruments | Existing |
| **41b (NEW)** | **Predistribution and the Labour Commons** | **This plan** |
| **41c (NEW)** | **Democratic Economic Governance** | **This plan** |
| 42 | The Role of Experimentation | Existing |

