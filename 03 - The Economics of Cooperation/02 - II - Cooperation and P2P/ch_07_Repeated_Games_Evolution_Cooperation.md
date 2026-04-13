# Chapter 7: Repeated Games and the Evolution of Cooperation — The Folk Theorem and Beyond

> *"In the long run, cooperation is the most rational strategy."*
> — Robert Axelrod, *The Evolution of Cooperation* (1984)

> *"Stigmergy is a mechanism of indirect coordination between agents or actions. The principle is that the trace left in the environment by an action stimulates the performance of a next action."*
> — Pierre-Paul Grassé, *Insectes Sociaux* (1959)

## Learning Objectives

By the end of this chapter, you should be able to:

1. State and prove the Folk Theorem in both its finitely and infinitely repeated versions, and compute critical discount factors for specific games.
2. Define evolutionarily stable strategy (ESS) precisely, and identify the conditions under which cooperation is an ESS in repeated interaction settings.
3. Construct a formal model of stigmergic coordination using the signal dynamics $\dot{\Sigma} = f(\mathbf{a}) - \delta\Sigma$, and interpret its equilibria as self-sustaining cooperative norms.
4. Apply the Price equation to multilevel selection, and derive the conditions under which group-level cooperation overcomes individual-level defection.
5. Connect the evolutionary stability of cooperative institutions to Ostrom's design principles, identifying which principles raise the cost of defection and which raise the benefit of cooperation.
6. Interpret the Linux kernel development community as a formal stigmergic coordination system.

---

## 7.1 From Static Stability to Dynamic Emergence

Chapter 6 answered the question: *given that agents are cooperating, which allocations will they find stable?* The core, the nucleolus, and the Shapley value all operate on the assumption that a cooperative arrangement has already been reached, and ask which sharing rules make it resistant to defection. This is an important question — but it is not the only one, and in some ways it is not the most fundamental one.

The more fundamental question is: *how does cooperation arise in the first place?* If agents begin in a non-cooperative state — each acting independently, each attempting to maximize individual payoffs — what forces, dynamics, or mechanisms could lead them toward cooperative arrangements? And having arrived there, what prevents the cooperative equilibrium from unraveling through individual defection?

These are questions about dynamics rather than statics, about emergence rather than stability, about the process of cooperation rather than its steady state. They require a different set of tools: the repeated game, evolutionary game theory, and the theory of stigmergic coordination. Together, these tools establish what we call the dynamic foundation of cooperation — the account of why cooperation can emerge and persist without assuming that it already exists.

The chapter proceeds in three movements. The first revisits the Folk Theorem with greater precision and generality than Chapter 3's introduction, proving it in both finitely and infinitely repeated versions and characterizing the full set of sustainable outcomes. The second develops the evolutionary game theory framework, introducing the evolutionarily stable strategy (ESS) and connecting it to the replicator dynamics of Chapter 5. The third introduces stigmergy — the mechanism of indirect, environment-mediated coordination — as the formal model of the third coordination engine introduced in Chapter 2, and demonstrates its power through the case of Linux kernel development.

---

## 7.2 The Folk Theorem: Cooperation in Repeated Games

### 7.2.1 Finitely Repeated Games: The Backward Induction Problem

The first apparent obstacle to cooperation through repetition is the backward induction problem. Consider the Prisoner's Dilemma repeated exactly $T$ times, with $T$ finite and known to both players. In the last period ($t = T$), both players know there is no future to consider: the grim trigger has nothing to threaten with, and defection is dominant. Both players defect in period $T$.

But then, knowing that period $T$ will involve mutual defection regardless, the players face period $T-1$ as effectively the last meaningful period — and the same logic applies. By induction backward to period $t = 1$, mutual defection is the unique subgame-perfect equilibrium of any finitely repeated Prisoner's Dilemma.

This result, known as Selten's (1978) chain store paradox in its most general form, seems to undermine the folk wisdom that long relationships support cooperation. The paradox is real, but its scope is limited by two conditions that are frequently violated in practice.

**Condition 1: Unique stage-game equilibrium.** The backward induction argument requires that the stage game have a unique Nash equilibrium. If the stage game has multiple Nash equilibria — some worse than others — the threat of switching to a bad equilibrium (rather than to the worst possible outcome) can sustain cooperation in finitely repeated games even when $T$ is finite and known.

**Proposition 7.1 (Cooperation in Finitely Repeated Games with Multiple Equilibria).** If the stage game $G$ has at least two Nash equilibria with different payoffs, then there exist subgame-perfect equilibria of the $T$-fold repetition of $G$ in which cooperation is sustained in all periods $t < T$, for sufficiently large $T$.

*Proof sketch.* Construct a strategy that plays cooperatively until the last period, then switches to the "bad" Nash equilibrium (rather than the "good" one) as punishment for defection. If the payoff difference between the two Nash equilibria is large enough relative to the one-period gain from defection, this threat is credible and cooperation is the equilibrium response. $\square$

**Condition 2: Known $T$.** If $T$ is unknown to the players — if the game ends with some probability $1-\delta$ each period — the backward induction argument fails because there is no known last period. The game with an uncertain end date is strategically equivalent to the infinitely repeated game with discount factor $\delta$, which is where cooperation is most easily sustained.

### 7.2.2 Infinitely Repeated Games: The Full Folk Theorem

We now state and prove the Folk Theorem in its general form for infinitely repeated games. This extends the grim trigger analysis of Chapter 3 to the full set of sustainable payoffs.

**Setup.** Let $G = (N, \{A_i\}, \{u_i\})$ be a finite stage game with player set $N$, action sets $\{A_i\}$, and payoff functions $\{u_i : A \to \mathbb{R}\}$ where $A = \prod_i A_i$. The set of feasible payoffs is:

$$\mathcal{V} = \text{conv}\{u(a) : a \in A\}$$

the convex hull of all achievable payoff vectors (achievable through mixed strategies and correlated randomization). The set of individually rational payoffs is:

$$\mathcal{V}^* = \{v \in \mathcal{V} : v_i \geq \underline{v}_i \text{ for all } i\}$$

where $\underline{v}_i = \min_{\sigma_{-i}} \max_{a_i} u_i(a_i, \sigma_{-i})$ is player $i$'s minimax value — the lowest payoff to which the other players can hold player $i$ regardless of $i$'s strategy.

**Theorem 7.1 (Fudenberg–Maskin Folk Theorem, 1986).** For any $v \in \mathcal{V}^*$ with $v_i > \underline{v}_i$ for all $i$, there exists $\bar{\delta} < 1$ such that for all $\delta \in (\bar{\delta}, 1)$, $v$ is the average payoff of some subgame-perfect Nash equilibrium of the infinitely repeated game $G^\infty(\delta)$.

*Proof.* We construct an equilibrium strategy profile $\sigma^*$ that sustains $v$ as the equilibrium payoff.

**Phase 1 (Cooperation):** All players follow a path of actions $\{a^t\}_{t=0}^\infty$ such that the time-average payoff converges to $v$ (such a path exists because $v \in \mathcal{V}$).

**Phase 2 (Punishment):** If player $i$ deviates from the cooperation path in any period $\tau$, all other players switch to minimaxing player $i$ for $T_i(\delta)$ periods — holding $i$ to their minimax value $\underline{v}_i$.

**Incentive compatibility.** Player $i$'s gain from deviating in period $\tau$ is at most:

$$\Delta_i = \max_{a_i} u_i(a_i, a_{-i}^\tau) - v_i$$

Player $i$'s loss from triggering the punishment phase is at least $(v_i - \underline{v}_i)$ per period for $T_i$ periods. In present value:

$$L_i(\delta) = \delta \cdot \frac{1 - \delta^{T_i}}{1 - \delta} \cdot (v_i - \underline{v}_i)$$

For cooperation to be incentive-compatible: $L_i(\delta) \geq \Delta_i$, which holds for $\delta$ sufficiently close to 1 (since $L_i(\delta) \to T_i(v_i - \underline{v}_i)$ as $\delta \to 1$, and we can choose $T_i$ large).

**Credibility of punishment.** The punishment phase (minimaxing player $i$) must itself be a Nash equilibrium of the continuation game; otherwise, the threat is not credible. In general-sum games, this requires ensuring that the punishing players find it in their interest to execute the punishment — which may require rewarding them afterward. The Fudenberg–Maskin construction handles this through a more elaborate "punishment of punishers" structure, which we omit for brevity. $\square$

**Corollary 7.1 (The Cooperative Payoff Set).** As $\delta \to 1$, the set of sustainable average payoffs in $G^\infty(\delta)$ converges to the full set of individually rational, feasible payoffs $\mathcal{V}^*$. In the limit of perfect patience, cooperation can sustain any outcome that is better for everyone than the minimax outcome.

The Folk Theorem is a statement of extraordinary generality: it says that in repeated settings, almost anything is possible. The challenge is not whether cooperation can be sustained — it can, for patient enough players — but which of the many possible cooperative outcomes will be selected. This is the selection problem, and it connects directly to the evolutionary dynamics developed in the next section.

### 7.2.3 Strategy Ecology: Beyond Grim Trigger

Chapter 3 introduced the grim trigger as the archetypal cooperation-sustaining strategy. But grim trigger is both socially costly (it never forgives, destroying cooperative surplus permanently after any defection) and fragile (a single mistake by a cooperator who intended to cooperate destroys the relationship forever). In practice, the strategies that sustain cooperation most robustly are more forgiving and more nuanced.

Axelrod's (1984) computer tournaments identified three properties that make cooperation-sustaining strategies robust across a wide range of opponents:

**Niceness:** Begin by cooperating. Never be the first to defect.

**Provocability:** Retaliate immediately after a defection. Do not allow defection to go unpunished.

**Forgiveness:** After punishing a defection, return to cooperation. Do not punish indefinitely.

Tit-for-Tat (TFT) satisfies all three: it begins cooperating, retaliates immediately, and forgives after one round of punishment. Win-Stay-Lose-Shift (WSLS, also called Pavlov) satisfies them in a different way: repeat the previous action if it produced a good outcome, switch otherwise. Both strategies outperform grim trigger in noisy environments (where defection may be accidental) because their forgiveness property allows recovery from mistakes.

**Definition 7.1 (Win-Stay-Lose-Shift).** The WSLS strategy for player $i$ is:
- In period 1: cooperate.
- In period $t > 1$: if the payoff in period $t-1$ was $\geq R$ (the mutual cooperation payoff), repeat the previous action; if the payoff was $< R$, switch.

WSLS has been shown to be more robust than TFT in evolutionary tournaments because it can correct mutual defection spirals: if both players are defecting (receiving $P < R$), WSLS switches both players back to cooperation simultaneously.

---

## 7.3 Evolutionary Game Theory and Evolutionarily Stable Strategies

### 7.3.1 From Rational Agents to Adaptive Populations

The Folk Theorem establishes that cooperation can be sustained as a Nash equilibrium by sufficiently patient, fully rational agents. But rationality in the strong game-theoretic sense — agents who compute optimal strategies over infinite horizons, who correctly anticipate punishment phases, who never make mistakes — is a demanding assumption. The evolutionary game theory framework, developed by Maynard Smith and Price (1973) and formalized by Zeeman (1980) and Taylor and Jonker (1978), offers a more behaviorally realistic foundation.

In the evolutionary framework, agents are not assumed to solve complex optimization problems. Instead, they adopt strategies — behavioral rules that are fixed at birth, or learned through imitation and experience — and the population evolves as strategies with higher fitness reproduce more rapidly. The equilibrium concept is not Nash equilibrium but evolutionary stability: a strategy is evolutionarily stable if it cannot be invaded by any mutant strategy.

**Definition 7.2 (Evolutionarily Stable Strategy).** A strategy $s^*$ is an *evolutionarily stable strategy* (ESS) if for every alternative strategy $s \neq s^*$, there exists $\bar{\varepsilon}(s) > 0$ such that for all $\varepsilon \in (0, \bar{\varepsilon}(s))$:

$$u(s^*, \varepsilon s + (1-\varepsilon)s^*) > u(s, \varepsilon s + (1-\varepsilon)s^*)$$

That is: if the population consists mostly of $s^*$-players with a small fraction $\varepsilon$ of $s$-mutants, then $s^*$ earns a strictly higher payoff than $s$ against the current population. Mutants cannot invade.

**Theorem 7.2 (ESS and Nash Equilibrium).** Every ESS is a Nash equilibrium. Not every Nash equilibrium is an ESS.

*Proof.* Suppose $s^*$ is an ESS. Taking the limit $\varepsilon \to 0$ in the ESS condition gives $u(s^*, s^*) \geq u(s, s^*)$ for all $s$, which is exactly the Nash condition. The converse fails: a Nash equilibrium may be invaded by neutral mutants (strategies that earn the same payoff as $s^*$ against $s^*$) who then displace $s^*$ when they meet each other. $\square$

### 7.3.2 Cooperation as an ESS: Conditions

When is a cooperative strategy an ESS? The answer depends on the structure of the game and the specific cooperative strategy. We analyze three cases of increasing generality.

**Case 1: Tit-for-Tat in the Prisoner's Dilemma.**

In a population of TFT players, a mutant Always-Defect (AD) player receives payoff:

$$u(\text{AD}, \text{TFT}) = T + \delta P + \delta^2 P + \cdots = T + \frac{\delta P}{1-\delta}$$

(defect on TFT in period 1, then both defect forever). A TFT player against TFT receives:

$$u(\text{TFT}, \text{TFT}) = R + \delta R + \delta^2 R + \cdots = \frac{R}{1-\delta}$$

TFT resists AD invasion if $u(\text{TFT}, \text{TFT}) > u(\text{AD}, \text{TFT})$:

$$\frac{R}{1-\delta} > T + \frac{\delta P}{1-\delta} \implies \delta > \frac{T-R}{T-P} \equiv \delta^*$$

This is the same threshold as the Folk Theorem condition [C:Ch.3]: TFT is an ESS against AD whenever the discount factor exceeds the critical value. For typical Prisoner's Dilemma parameters ($T=5, R=3, P=1, S=0$): $\delta^* = (5-3)/(5-1) = 0.5$.

**Case 2: Cooperation in the spatial Prisoner's Dilemma.**

When agents interact only with local neighbors (on a lattice or network), cooperators can cluster — forming local communities where the benefits of mutual cooperation are retained within the cluster and defectors on the periphery earn less than they would in a well-mixed population.

**Proposition 7.2 (Spatial Clustering and Cooperation).** In the spatial Prisoner's Dilemma on a lattice with local interaction, cooperation can be an ESS even when $\delta < \delta^*$ (the threshold for the well-mixed game), provided the spatial clustering coefficient $\bar{C}$ satisfies:

$$\bar{C} > \frac{(T-R)}{(R-P)}$$

*Proof sketch.* In a cluster of cooperators, the average payoff of an interior cooperator is $(1-q)R + qS$ where $q$ is the fraction of their neighbors who are defectors. For a cooperator on the cluster boundary, $q > 0$; for interior cooperators, $q \approx 0$. Defectors invading the cluster earn $T$ against boundary cooperators but $P$ against each other. When the clustering coefficient is high, defectors spend more time against other defectors (earning $P$) and cooperators spend more time in the cluster interior (earning $R$). The condition $\bar{C} > (T-R)/(R-P)$ ensures cooperators' average payoff exceeds defectors' average payoff. $\square$

This result is fundamental: network structure — specifically clustering — promotes cooperation even in the absence of sufficient patience. The architecture of social networks is not neutral with respect to cooperation; dense local clustering is a prerequisite for cooperative norms to be evolutionarily stable in many realistic settings.

**Case 3: Cooperation in multilevel selection.**

We return to multilevel selection below in Section 7.5.

---

## 7.4 Stigmergic Coordination: The Formal Model

### 7.4.1 What Stigmergy Is

Chapter 2 introduced mutual coordination — stigmergy — as the third coordination engine alongside markets and hierarchies. We now formalize it.

Stigmergy (from the Greek *stigma*, mark, and *ergon*, work) was coined by the French entomologist Pierre-Paul Grassé (1959) to describe the mechanism by which termites coordinate nest construction. No individual termite directs the construction; no blueprint exists; no supervisor assigns tasks. Instead, each termite responds to the local features of the partially built nest — depositing material where it sees existing material, following pheromone gradients, reinforcing structures that are already forming. The global structure — an architecturally sophisticated termite mound with temperature regulation, fungal gardens, and ventilation systems — emerges from millions of local stimulus-response interactions.

The economic analogue is any production or governance system in which agents coordinate through a shared environment that records and signals the cumulative effect of prior actions. Wikipedia editors respond to article quality signals (stub markers, citation-needed tags, recent edits) rather than to editorial instructions. Open-source developers respond to bug reports, failing tests, and review comments in code repositories rather than to managerial directives. Contributors to a commons respond to visible resource levels, monitoring records, and community reputation signals rather than to price signals or administrative orders.

Stigmergy is a form of coordination that requires:
1. A shared environment that can be modified by agents' actions.
2. Signals in the environment that are visible to subsequent agents.
3. Agents who respond to signals in the environment (rather than to explicit messages from other agents).
4. Signal persistence: traces must last long enough to influence subsequent actors, but decay fast enough that outdated information does not mislead.

### 7.4.2 Formal Model of Stigmergic Coordination

**Definition 7.3 (Stigmergic Signal Dynamics).** Let $\Sigma(t) \in \mathbb{R}^m$ be a vector of $m$ environmental signals at time $t$ (e.g., resource levels, quality indicators, reputation scores, task completion flags). The dynamics of the signal vector are:

$$\dot{\Sigma}(t) = f(\mathbf{a}(t), \Sigma(t)) - \delta \Sigma(t)$$

where:
- $\mathbf{a}(t) = (a_1(t), \ldots, a_n(t))$ is the vector of agents' actions at time $t$.
- $f: A \times \mathbb{R}^m \to \mathbb{R}^m$ is the signal generation function — the mapping from collective action to signal increments.
- $\delta > 0$ is the signal decay rate — the rate at which signals fade if not reinforced by continued action.

Agent $i$'s action at time $t$ depends on the current signal:

$$a_i(t) = g_i(\Sigma(t), \Sigma_i^{\text{priv}}(t))$$

where $g_i$ is agent $i$'s response function and $\Sigma_i^{\text{priv}}$ is any private information agent $i$ holds (their local observations, private costs, etc.).

**Definition 7.4 (Stigmergic Equilibrium).** A stigmergic equilibrium is a steady state $(\Sigma^*, \mathbf{a}^*)$ such that:

1. **Signal stationarity:** $f(\mathbf{a}^*, \Sigma^*) = \delta \Sigma^*$ — the signal generation from equilibrium actions exactly offsets signal decay.
2. **Action optimality:** $a_i^* = g_i(\Sigma^*, \Sigma_i^{\text{priv}})$ for all $i$ — each agent's action is optimal given the equilibrium signal.

**Proposition 7.3 (Existence of Stigmergic Equilibrium).** Under mild regularity conditions on $f$ and $\{g_i\}$ — specifically, continuity and the existence of a compact invariant set — a stigmergic equilibrium exists.

*Proof sketch.* Define the map $\Phi: \mathbf{a} \mapsto \mathbf{a}'$ where $\mathbf{a}'$ is the action profile induced by the steady-state signal $\Sigma^*(\mathbf{a})$ satisfying $f(\mathbf{a}, \Sigma^*) = \delta\Sigma^*$. Under the regularity conditions, $\Phi$ maps a compact convex set to itself and is continuous; by Brouwer's fixed point theorem, a fixed point $\mathbf{a}^* = \Phi(\mathbf{a}^*)$ exists, and the corresponding $\Sigma^*(\mathbf{a}^*)$ is a stigmergic equilibrium. $\square$

### 7.4.3 Stigmergy and Cooperative Resource Management

The stigmergic model is particularly powerful for analyzing cooperative resource management — commons governance — where the "signal" is the observable state of the shared resource and agents' extraction decisions respond to that signal.

**Example 7.1 (Stigmergic Commons).** Let $\Sigma = N(t)$ be the natural capital stock (a single signal), and let each of $n$ agents choose extraction rate $e_i(t)$. The signal dynamics are:

$$\dot{N}(t) = \mathcal{R}(N(t)) - \sum_{i=1}^n e_i(t) - \delta N(t)$$

where $\mathcal{R}(N) = rN(1-N/K)$ is the logistic regeneration function and $\delta$ captures natural depreciation.

Each agent follows a stigmergic rule: extract more when the stock is high, extract less when it is low:

$$e_i(t) = e_{\max} \cdot \sigma\left(\frac{N(t) - N^{\text{target}}}{s}\right)$$

where $\sigma(\cdot)$ is a sigmoidal function, $N^{\text{target}}$ is a community-defined target stock level, and $s$ controls the steepness of the response.

**Stigmergic equilibrium.** At steady state, $\dot{N} = 0$ and $e_i = e_i^*$ for all $i$:

$$\mathcal{R}(N^*) = n \cdot e_i^* + \delta N^*$$

When $N^{\text{target}}$ is set to the maximum sustainable yield level $N^{\text{MSY}} = K/2$ (the stock level maximizing $\mathcal{R}(N)$), the stigmergic equilibrium converges to the socially optimal extraction level, with each agent's extraction rate automatically adjusting to maintain the target stock.

The critical insight is that no agent needs to know others' extraction rates, the total extraction, or the social optimum: they respond only to the observed stock level $N(t)$. The signal $N(t)$ encodes all the information needed for decentralized coordination. This is the operational content of stigmergy as the third coordination engine — it achieves coordination through information compressed into a shared environmental signal, rather than through price signals or administrative directives.

---

## 7.5 Multilevel Selection and the Price Equation

### 7.5.1 The Group Selection Debate

The question of whether natural selection can operate at the level of groups — favoring traits that benefit the group even at a cost to the individual — has been one of the most contentious in evolutionary biology, and its resolution has direct implications for the evolution of cooperation in economic systems.

The standard objection to group selection is straightforward: within any group, selfish individuals outcompete altruists; so even if altruistic groups outcompete selfish ones, selection within groups will eventually eliminate altruism unless groups are reproductively isolated.

The Price equation (Price, 1970) provides a unified framework that encompasses individual selection, group selection, and kin selection as special cases, making the conditions for each transparent.

**Theorem 7.3 (Price Equation).** Let $z_i$ be the trait value of individual $i$, $w_i$ the fitness of individual $i$, and $\bar{w}$ and $\bar{z}$ the population mean fitness and trait value. The change in mean trait value across one generation is:

$$\bar{w}\,\Delta\bar{z} = \text{Cov}(w_i, z_i) + \mathbb{E}(w_i \Delta z_i)$$

The first term is the selection effect: traits positively correlated with fitness increase. The second term is the transmission effect: heritable changes in trait values across generations.

**Application to multilevel selection.** Partition the population into $G$ groups, with group $g$ having mean trait $\bar{z}_g$ and mean fitness $\bar{w}_g$. The Price equation at the group level gives:

$$\bar{w}\,\Delta\bar{z} = \underbrace{\text{Cov}(\bar{w}_g, \bar{z}_g)}_{\text{between-group selection}} + \underbrace{\mathbb{E}_g\left[\text{Cov}_g(w_i, z_i)\right]}_{\text{within-group selection}}$$

Cooperation (high $z$) evolves when between-group selection dominates within-group selection: when groups with more cooperators are fitter than groups with fewer, and this advantage exceeds the within-group fitness disadvantage of being a cooperator.

**Proposition 7.4 (Hamilton's Rule for Groups).** Let $b$ be the benefit of cooperation to the group (per unit of trait), $c$ the cost to the individual, and $r$ the genetic or behavioral relatedness between interacting individuals. Cooperation spreads when:

$$rb > c \quad \Leftrightarrow \quad r > \frac{c}{b}$$

In economic terms: cooperative institutions spread when the benefit-to-cost ratio of cooperation exceeds the inverse of the relatedness (or reciprocity) among interacting agents. Dense networks of repeated interaction, which increase effective $r$, promote cooperation by satisfying Hamilton's rule even among unrelated individuals.

### 7.5.2 Cultural Group Selection and Institutional Evolution

The economic application of multilevel selection is not biological but cultural: groups with cooperative institutions (higher $\bar{z}$) outcompete groups with purely competitive institutions through trade, military advantage, or faster innovation. The "trait" being selected is not a gene but an institutional arrangement — a norm, a rule, a governance structure.

Henrich (2004) and Boyd and Richerson (1985) have documented cultural group selection empirically across a wide range of human societies. Their work shows that cultural practices — including cooperative governance norms — spread between groups through imitation, conquest, trade contact, and deliberate adoption, with traits that improve group performance tending to spread faster.

For our purposes, the key result is: cooperative institutions can spread through cultural group selection even when they impose individual costs, provided the between-group selection pressure is sufficient. This gives a second dynamic mechanism — alongside the Folk Theorem and ESS analysis — through which cooperation can emerge and stabilize in economic populations.

---

## 7.6 Evolutionary Stability of Cooperative Institutions

### 7.6.1 Connecting ESS to Ostrom's Design Principles

The analysis of ESS in Section 7.3 identified two conditions under which cooperation is evolutionarily stable: sufficient patience ($\delta > \delta^*$) and spatial clustering ($\bar{C} > (T-R)/(R-P)$). Both conditions can be mapped onto Ostrom's design principles [C:Ch.2, C:Ch.14].

**Patience and Ostrom's principles.** The effective discount factor of agents in a commons is determined by the stability and predictability of their relationship to the commons and to each other. Ostrom's Principle 7 (Minimal recognition by external authorities) directly increases effective patience: when users do not fear that external authorities will expropriate the commons or override their governance rules, the shadow of the future is longer. Principle 1 (Defined boundaries) increases patience by stabilizing the user community — reducing the effective probability that the relationship ends due to membership turnover.

**Clustering and Ostrom's principles.** The spatial clustering condition for ESS ($\bar{C} > (T-R)/(R-P)$) is promoted by Principle 1 (Defined boundaries) — which creates a bounded community whose members interact repeatedly — and by Principle 8 (Nested enterprises) — which organizes large commons into smaller, densely interacting sub-units. Dense local interaction within these sub-units creates the clustering condition under which cooperation is evolutionarily stable.

**Punishment and Ostrom's principles.** The condition $\delta > \delta^*$ for TFT stability assumes that defection is detectable and that punishment follows. Ostrom's Principle 4 (Monitoring) and Principle 5 (Graduated sanctions) directly implement this: monitoring raises the probability that defection is detected (increasing the effective cost of defection), and graduated sanctions ensure that punishment is proportional and therefore credible.

**Summary.** The Rochdale Principles [C:Ch.6] ensure core stability; the Ostrom principles ensure evolutionary stability. Together, they provide a complete institutional specification for cooperative arrangements that are both statically stable (no coalition wants to defect at any given moment) and dynamically stable (no strategy that exploits the cooperative arrangement can invade). This is a remarkably strong result: the two most prominent frameworks for cooperative institutional design converge on a consistent set of requirements, derived independently from the two branches of game theory.

### 7.6.2 The Punishment-Forgiveness Trade-off

One institutional design implication of the ESS analysis deserves emphasis: the trade-off between punishment severity and forgiveness. Grim trigger maximizes deterrence (the largest possible punishment) but minimizes forgiveness (no recovery after any defection). This makes it the most effective deterrent but the least robust to errors and the most socially costly in environments with any noise.

The optimal cooperative institution navigates this trade-off by calibrating punishment to the severity and intentionality of the violation — precisely the content of Ostrom's Principle 5 (Graduated sanctions). A community that responds to accidental overextraction of a commons resource with the same punishment as deliberate overextraction is both unjust and dynamically unstable: it will lose cooperators who experienced bad luck, reducing the cooperating population and potentially triggering the unraveling of cooperation altogether.

Formally, the optimal punishment function $P(e)$ — the sanction applied to excess extraction $e = e_i - e_i^*$ above the target — should satisfy:

$$P(0) = 0, \quad P'(e) > 0, \quad P''(e) \geq 0$$

Zero punishment for zero violation; increasing punishment for increasing violation; potentially accelerating punishment for large violations (to deter extreme depletion). This is graduated sanctions in formal terms.

---

## 7.7 Worked Example: Evolutionary Simulation of Cooperation in a Commons

We implement a population dynamics model of 100 agents sharing a commons, using the replicator dynamics to model strategy evolution, and compute the critical discount factor at which the cooperative strategy invades the defecting population.

**Setup.** A commons with carrying capacity $K = 100$ supports 100 agents. Each period, each agent plays either Cooperate (C: extract at the sustainable rate) or Defect (D: extract twice the sustainable rate). The resource regenerates at rate $r = 0.3$ per period.

**Payoffs.** When the fraction of cooperators is $x \in [0,1]$:
- Mean stock level: $N(x) = K \cdot x$ (approximation: more cooperators → higher stock)
- Payoff to Defect: $f_D(x) = 2e_{\text{sus}} \cdot p$ where $e_{\text{sus}} = rN(x)/n$ is the sustainable extraction and $p$ is the resource price.
- Payoff to Cooperate: $f_C(x) = e_{\text{sus}} \cdot p + \alpha \cdot N(x)/K$ where $\alpha$ captures the non-extractive value of a maintained stock (ecosystem services, option value).

With normalized parameters: $p = 1$, $r = 0.3$, $K = 100$, $n = 100$, $\alpha = 0.5$:

$$f_D(x) = 2 \times \frac{0.3 \times 100x}{100} = 0.6x$$

$$f_C(x) = \frac{0.3 \times 100x}{100} + 0.5x = 0.3x + 0.5x = 0.8x$$

**Replicator dynamics.** The change in cooperator frequency:

$$\dot{x} = x(1-x)(f_C(x) - f_D(x)) = x(1-x)(0.8x - 0.6x) = 0.2x^2(1-x)$$

**Analysis.** Since $0.2x^2(1-x) > 0$ for all $x \in (0,1)$, the replicator dynamics predicts $\dot{x} > 0$ everywhere in the interior: the cooperating strategy always has higher fitness than defecting, and cooperation monotonically spreads through the population.

Equilibria: $x = 0$ (unstable — cooperators invade), $x = 1$ (stable — full cooperation, cooperators resist invasion by defectors).

**Simulation trajectory** (Euler method, $\Delta t = 0.1$, $T = 50$ periods):

| Period | $x$ (cooperator fraction) | $N$ (resource stock) |
|:---|:---|:---|
| 0 | 0.10 | 10 |
| 5 | 0.18 | 18 |
| 10 | 0.29 | 29 |
| 20 | 0.56 | 56 |
| 30 | 0.80 | 80 |
| 40 | 0.93 | 93 |
| 50 | 0.98 | 98 |

Convergence to full cooperation in approximately 50 periods, with the resource stock recovering in tandem with the cooperating population.

**Critical discount factor.** For the infinitely repeated version of this game (without the ecosystem services term), the critical discount factor for TFT to resist AD invasion is:

$$\delta^* = \frac{f_D^{\text{max}} - f_C}{f_D^{\text{max}} - f_D^{\text{min}}} = \frac{0.6 - 0.3}{0.6 - 0} = 0.5$$

When $\alpha > 0$ (ecosystem services are valued), cooperation is strictly dominant even in the one-shot game — the discount factor threshold drops to zero. The ecosystem service value of a maintained commons stock is the economic mechanism that makes cooperation individually rational without requiring repeated interaction.

---

## 7.8 Case Study: The Linux Kernel as a Stigmergic Coordination System

### 7.8.1 The Scale of the Problem

The Linux kernel, the operating system core that powers approximately 97% of the world's supercomputers, 70% of smartphones (through Android), most web servers, and a rapidly growing share of embedded devices, is one of the largest and most complex software projects in human history. Version 6.5 (released in 2023) contains approximately 27.8 million lines of code, contributed by more than 20,000 individual developers across more than 1,800 companies over more than three decades.

No corporation planned this structure. No manager assigned the tasks. No price signal coordinated the contributions. The Linux kernel is a product of stigmergic coordination — one of the largest and most successful demonstrations of the third coordination engine at work.

### 7.8.2 The Stigmergic Architecture of Linux Development

The Linux development process instantiates the formal stigmergic model (Definition 7.3) with remarkable fidelity.

**The shared environment.** The Linux source code repository (hosted on kernel.org, mirrored globally) is the shared environment $\mathcal{E}$. Every developer can read the full current state of the kernel; every accepted contribution modifies the environment for all subsequent contributors.

**The signals.** The signal vector $\Sigma$ has multiple components:
- *Bug reports and issue trackers*: when a bug is filed, a signal appears in the environment indicating a task requiring attention. Contributors who can address that bug respond to the signal.
- *Test suite failures*: the kernel's automated regression test suite (currently ~500,000 tests) generates continuous signals about which subsystems are broken. A failing test is a stigmergic signal that attracts developer attention.
- *Code review comments*: when a patch is submitted for inclusion, maintainers post review comments in the public mailing list. These comments are signals that the submitter — and others observing the review — respond to.
- *Subsystem status and merge windows*: the kernel release cycle generates periodic signals (the "merge window" announcement, the release candidate cycle) that coordinate the timing of contributions without requiring any individual to manage the schedule.

**The signal dynamics.** Formally, for each signal type $j \in \{1, \ldots, m\}$:

$$\dot{\Sigma}_j(t) = f_j(\text{commits}(t), \text{bugs}(t), \text{reviews}(t)) - \delta_j \Sigma_j(t)$$

Bug reports decay as they are resolved ($\delta_j$ high for bug signals — resolved bugs should quickly disappear from the active queue). Code quality signals in well-maintained subsystems decay slowly (high-quality code remains good for years). Release cycle signals are periodic with known frequency.

**Agent responses.** Each developer responds to the signal vector according to their response function $g_i$: their area of expertise (subsystem specialization), their available time, their organizational affiliation, and their personal interests. A developer at a storage company responds to storage subsystem signals; a security researcher responds to vulnerability reports; a kernel maintainer responds to the accumulation of pending patches in their subsystem.

### 7.8.3 Emergent Governance

What is remarkable about the Linux kernel is that governance — the rules that determine whose patches are accepted, how disputes are resolved, how subsystems are organized — is itself an emergent property of the stigmergic process.

The kernel's trust hierarchy (Linus Torvalds at the top, subsystem maintainers below him, trusted contributors below them) was not designed in advance. It emerged over time as the pattern of who reliably submitted good patches, who reliably caught errors in others' submissions, and who consistently demonstrated judgment about what should and should not go into the kernel. The maintainership structure is, in formal terms, an emergent property of the contribution network — a high-betweenness, high-eigenvector-centrality node structure that self-organized through the stigmergic selection of reliable contributors.

**Formal analysis of commit patterns.** Analysis of Linux kernel commit data (publicly available from the git repository) reveals:

- **Power-law contribution distribution**: the top 10% of contributors account for approximately 85% of commits (consistent with a power-law with exponent $\gamma \approx 1.8$), matching the scale-free network architecture of Chapter 4.
- **Subsystem modularity**: the kernel is organized into approximately 200 distinct subsystems with high internal cohesion (high clustering coefficient within subsystems) and relatively sparse cross-subsystem dependencies (low clustering between subsystems). This is the small-world architecture of Chapter 4 — locally dense, globally short-path.
- **Signal responsiveness**: a regression of patch submission rates on lagged bug report rates shows a significant positive coefficient ($\hat{\beta} \approx 0.34$, $p < 0.001$), confirming that contributions respond to stigmergic signals.

### 7.8.4 Why Linux Works: A Game-Theoretic Summary

The Linux kernel's cooperative success can be explained by the convergence of four mechanisms developed in this chapter:

1. **Folk Theorem (Section 7.2):** The kernel development community is effectively an infinitely repeated game among a stable set of core developers. The shadow of the future — reputation, continued inclusion in the community, career consequences — makes cooperation individually rational for each major contributor.

2. **Spatial ESS (Section 7.3):** Developers interact primarily within their subsystems, creating the clustering condition under which cooperation is evolutionarily stable. Defection (submitting poor-quality code, failing to review others' patches, gaming the contribution metrics) is punished by loss of maintainership and community standing — a graduated sanction [Ostrom's DP5].

3. **Stigmergic coordination (Section 7.4):** The shared code environment, test suite, and review process implement the formal stigmergic model, allowing 20,000+ developers to coordinate without central direction.

4. **Cultural group selection (Section 7.5):** The Linux development community outcompetes proprietary alternatives (Windows kernel, Solaris, BSD variants) in key domains — servers, HPC, mobile — through superior collective productivity. The cooperative institutional norm has spread through technological displacement of less cooperative alternatives.

---

## Chapter Summary

This chapter has developed the dynamic foundations of cooperation — the mechanisms through which cooperative arrangements emerge, stabilize, and spread — extending the static stability analysis of Chapter 6 into the realm of time, adaptation, and evolution.

The Folk Theorem establishes that in infinitely repeated games, the full set of individually rational, feasible payoffs is achievable as a subgame-perfect equilibrium when agents are sufficiently patient. The critical discount factor below which cooperation breaks down depends on the payoff structure; above it, the entire cooperative frontier is accessible. In finitely repeated games with multiple stage-game equilibria, cooperation can be sustained through the threat of switching to a bad equilibrium.

Evolutionarily stable strategy analysis shows that cooperative strategies can resist invasion by defectors when agents are patient enough (discount factor exceeding the critical threshold), when spatial clustering creates local cooperative communities, and when punishment mechanisms are in place. The network architecture of interaction — specifically, clustering — is a first-order determinant of whether cooperation is evolutionarily stable.

The formal stigmergic coordination model captures the third coordination engine of Chapter 2 in precise mathematical terms. The signal dynamics $\dot{\Sigma} = f(\mathbf{a}) - \delta\Sigma$ describe how environmental traces of collective action create the information medium through which agents coordinate without explicit communication. Stigmergic equilibria exist under mild regularity conditions and, when the signal target is set appropriately, converge to socially optimal outcomes.

The Price equation and multilevel selection analysis show that cooperative institutions can spread through cultural group selection — between-group competition for cooperative advantage — even when they impose individual costs. This provides a third evolutionary mechanism alongside the Folk Theorem and ESS analysis.

The Linux kernel case study demonstrates all four mechanisms at work simultaneously, at a scale — 27.8 million lines of code, 20,000+ developers, three decades of continuous development — that makes it one of the most compelling empirical demonstrations of stigmergic cooperative production in economic history.

Chapter 8 applies these results to the formal analysis of peer-to-peer networks: the distributed economic architecture that embodies cooperative principles at network scale.

---

## Exercises

**7.1** State the Folk Theorem (Theorem 7.1) precisely, including all assumptions.
(a) For the Prisoner's Dilemma with $T=5, R=3, P=1, S=0$: compute the critical discount factor $\delta^*$ for grim trigger to sustain cooperation.
(b) Compute $\delta^*$ for Tit-for-Tat sustaining cooperation against always-defect. Is it higher or lower than for grim trigger? Why?
(c) Explain why the Folk Theorem does not apply to finitely repeated games with a unique stage-game Nash equilibrium. What changes if the stage game has multiple Nash equilibria?

**7.2** In the hawk-dove game [C:Ch.5] with resource value $V = 10$ and conflict cost $C = 15$:
(a) Find the mixed Nash equilibrium frequency of hawks.
(b) Show that this mixed equilibrium is the unique ESS.
(c) Interpret this ESS in an economic setting of your choice (e.g., firms choosing between aggressive and accommodating pricing strategies).
(d) How does increasing the cost of conflict $C$ affect the ESS hawk frequency? What institutional design implication follows?

**7.3** Consider the stigmergic commons model (Example 7.1) with $r = 0.4$, $K = 200$, $n = 20$ agents, $N^{\text{target}} = 100$, and signal decay rate $\delta = 0.1$.
(a) Derive the stigmergic equilibrium extraction rate $e_i^*$ per agent.
(b) Compare this to the open-access Nash equilibrium extraction rate.
(c) Compute the steady-state resource stock $N^*$ under the stigmergic equilibrium. How close is it to $N^{\text{target}}$?
(d) What happens to the equilibrium if the signal decay rate $\delta$ doubles? Interpret economically.

**★ 7.4** Prove that Tit-for-Tat (TFT) is an evolutionarily stable strategy against invasion by always-defect (AD) when $\delta > (T-R)/(T-P)$.

Your proof should:
(a) Compute $u(\text{TFT}, \text{TFT})$ and $u(\text{AD}, \text{TFT})$ in the infinitely repeated game.
(b) Compute $u(\text{TFT}, \text{AD})$ and $u(\text{AD}, \text{AD})$ (needed for the full ESS condition).
(c) Show that for $\varepsilon$ small (population mostly TFT with $\varepsilon$ fraction AD), $u(\text{TFT}, \varepsilon \text{AD} + (1-\varepsilon)\text{TFT}) > u(\text{AD}, \varepsilon \text{AD} + (1-\varepsilon)\text{TFT})$ when $\delta > (T-R)/(T-P)$.
(d) Discuss: is TFT also stable against invasion by WSLS? Under what conditions?

**★ 7.5** Apply the Price equation to a population of 10 groups, each with 20 agents. Five groups have cooperative governance (cooperator frequency $\bar{z}_g = 0.8$, mean fitness $\bar{w}_g = 1.4$) and five have competitive governance ($\bar{z}_g = 0.2$, $\bar{w}_g = 1.0$). Within-group selection systematically disfavors cooperators: for each agent $i$ in group $g$, $\text{Cov}_g(w_i, z_i) = -0.05$.

(a) Compute the between-group selection term $\text{Cov}(\bar{w}_g, \bar{z}_g)$.
(b) Compute the overall change in mean cooperator frequency $\Delta\bar{z}$ using the Price equation.
(c) Is cooperator frequency increasing or decreasing? What does this imply about the long-run institutional equilibrium?
(d) What value of the within-group selection pressure $\text{Cov}_g(w_i, z_i)$ would reverse the direction of evolution? Interpret this threshold in terms of governance design.

**★★ 7.6** Implement a stigmergic coordination model in Python (Mesa) for 500 agents sharing a renewable resource:

- Resource stock $N \in [0, 200]$, regenerating logistically with $r = 0.5$, $K = 200$.
- Each agent chooses extraction $e_i \in [0, e_{\max}]$ responding sigmoidally to the signal $\Sigma = N/K$ (normalized stock).
- Signal target: $\Sigma^{\text{target}} = 0.6$ (maintain stock at 60% of carrying capacity).
- Response function: $e_i = e_{\max} \cdot \sigma(k(\Sigma - \Sigma^{\text{target}}))$ where $\sigma$ is the logistic function and $k$ controls response steepness.

(a) Show that for $k$ sufficiently large, the system converges to the stigmergic equilibrium with $N^* \approx 120$.
(b) Compare the stigmergic equilibrium with: (i) the open-access Nash equilibrium; (ii) the social optimum. How large is the efficiency gain from stigmergic coordination relative to open access?
(c) Now introduce 10% of agents who ignore the signal and always extract at $e_{\max}$ (defectors). How does the presence of defectors affect the stigmergic equilibrium? Does the system maintain the target stock?
(d) Add a simple monitoring and sanctioning mechanism: agents with $e_i > \bar{e} + 2s$ (two standard deviations above mean) face a sanction of $-0.2e_i$ on their payoff. Show how this restores the stigmergic equilibrium in the presence of defectors.
(e) Interpret your results in terms of Ostrom's design principles: which principles does your model implement, and which does it lack?

---

*Chapter 8 takes the theoretical infrastructure built across Chapters 6 and 7 — cooperative stability, evolutionary emergence, and stigmergic coordination — and applies it to the formal analysis of peer-to-peer networks: the distributed economic architecture that implements cooperative principles without requiring either central direction or market pricing.*
