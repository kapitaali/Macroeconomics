# Chapter 15: Emergent Institutions — How Rules Arise from Interactions

> *"Institutions are the rules of the game in a society, or, more formally, are the humanly devised constraints that shape human interaction."*
> — Douglass North, *Institutions, Institutional Change and Economic Performance* (1990)

> *"An institution is a regularity in social behavior that is agreed to by all members of a society, specifies behaviors in specific recurrent situations, and is either self-policed or policed by some external authority."*
> — Andrew Schotter, *The Economic Theory of Social Institutions* (1981)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Define an institution formally as an equilibrium of a coordination game, characterize the conditions for institutional stability, and distinguish designed from emergent institutions.
2. Apply replicator dynamics to populations of agents choosing among institutional variants, prove the basin-of-attraction dominance result, and identify when efficient institutions fail to replace inefficient ones.
3. Construct a formal model of stigmergic institutional formation — institutions as self-sustaining signal-behavior patterns — and identify the conditions for institutional crystallization.
4. Model institutional entrepreneurship formally as a cost-minimization problem, and derive the conditions under which crises reduce the cost of institutional transition.
5. Specify the formal conditions for lock-in to an inefficient institutional equilibrium and characterize the minimum perturbation required to escape.
6. Reconstruct the lex mercatoria as a stigmergic institution and explain what made it stable without state enforcement.

---

## 15.1 From Governance Design to Institutional Emergence

Chapters 13 and 14 treated governance as something designed: an authority structure chosen to optimize information aggregation, externality internalization, or cooperative stability. The Cosmo-Local model and the Ostrom principles are prescriptions — frameworks for constructing governance institutions that perform well by specified criteria.

But most institutions are not designed. They emerge from the accumulated interactions of agents who are pursuing local objectives, adapting to each other's behavior, and leaving traces in a shared environment that shape subsequent interactions. Property rights were not invented by economists; they evolved from patterns of occupancy, defense, and community recognition that crystallized into enforceable norms. Money was not decreed by a central authority; it emerged from the repeated selection of commodities that happened to be useful as media of exchange. Contract law was not legislated from scratch; it grew from the informal enforcement mechanisms of merchant communities. Even the most deliberately designed institutions — constitutions, regulatory frameworks, international treaties — are embedded in, and constrained by, a vast ecology of informal institutions that no one designed.

This chapter develops the formal theory of institutional emergence: how institutions arise from interaction, why some persist and others dissolve, and when emergent institutions are efficient or trapped in costly equilibria. The tools are the evolutionary game theory of Chapter 7, the stigmergic coordination model of the same chapter, and the complexity economics of Chapter 5. The result is a framework that explains not just what institutions should look like (Chapters 13–14) but what they do look like and why.

---

## 15.2 What Is an Institution? Three Formal Definitions

The word "institution" is used in at least three analytically distinct senses in the economics literature, each capturing something real and each having different implications for the theory of institutional emergence. We present all three and show how they relate.

### 15.2.1 North: Institutions as Rules

**Definition 15.1 (North, 1990).** An institution is a set of formal rules $\mathcal{F}$ and informal constraints $\mathcal{I}$ that structure human interaction, together with an enforcement mechanism $\mathcal{E}$:

$$\mathcal{N}\text{-institution} = (\mathcal{F}, \mathcal{I}, \mathcal{E})$$

Formal rules include laws, regulations, contracts, and property rights — rules that are codified and backed by state authority. Informal constraints include norms of behavior, conventions, and self-imposed codes of conduct — rules that are uncodified and backed by social sanction. The enforcement mechanism specifies how violations of either type of rule are detected and punished.

North's definition emphasizes the rule-setting function of institutions — their role in reducing uncertainty by providing a stable structure of expectations. Its weakness, from a dynamic perspective, is that it treats institutions as given rather than as outcomes of a process — it describes what institutions are but not how they come to be.

### 15.2.2 Ostrom: Institutions as Shared Rules-in-Use

**Definition 15.2 (Ostrom, 2005).** An institution is a set of rules-in-use: rules that are actually followed by agents in their interactions, as distinguished from rules-on-paper that are formally codified but behaviorally ineffective:

$$\mathcal{O}\text{-institution} = \{\text{rules-in-use} : \text{agents actually follow them}\}$$

Ostrom's definition focuses on behavioral effectiveness: a rule is part of an institution only if it actually governs behavior, not merely if it is written down. This distinction matters enormously in practice — most legal systems contain rules that are formally in force but behaviorally dead letters, while many behavioral regularities are enforced with great effectiveness by social norms that have no legal status.

### 15.2.3 Aoki: Institutions as Self-Sustaining Equilibria

**Definition 15.3 (Aoki, 2001).** An institution is a self-sustaining equilibrium of a game played repeatedly among agents:

$$\mathcal{A}\text{-institution} = \{(s_1^*, s_2^*, \ldots, s_n^*) : \text{a Nash equilibrium of the repeated interaction game}\}$$

More precisely, an institution is a system of shared beliefs about how the game is played, together with strategies that are mutually consistent with those beliefs. An institution exists when agents' beliefs about each other's behavior are correct (equilibrium beliefs) and when those beliefs are self-confirming — agents have no incentive to deviate from the behavior that sustains them.

Aoki's definition is the one most directly amenable to game-theoretic analysis, and it is the one we formalize most extensively in this chapter. It captures the central insight: an institution is stable not because it is enforced from outside but because it is in everyone's interest to maintain it given that everyone else is maintaining it. Institutions are, in this sense, self-enforcing coordination devices.

**The relationship between the three definitions.** North's rules provide the formal scaffolding within which Aoki's equilibria operate; Ostrom's rules-in-use identify which of the formally specified rules are actually part of the behavioral equilibrium. A complete account of institutions requires all three: what rules exist (North), which are actually followed (Ostrom), and why the followed rules are stable (Aoki).

---

## 15.3 Evolutionary Models of Institutional Emergence

### 15.3.1 Institutions as Coordination Game Equilibria

The simplest model of institutional emergence treats an institution as an equilibrium of a coordination game — a game in which multiple equilibria exist and agents have a common interest in coordinating on the same one.

**Definition 15.4 (Coordination Game).** A coordination game is a symmetric $n$-player game with action set $\mathcal{I} = \{I_1, I_2, \ldots, I_m\}$ (the available institutions) and payoff function:

$$u_i(I_k, I_{-i}) = \begin{cases} v(I_k) & \text{if } I_j = I_k \text{ for all } j \neq i \text{ (coordination)} \\ 0 & \text{if any } I_j \neq I_k \text{ (miscoordination)} \end{cases}$$

where $v(I_k) > 0$ is the coordination value of institution $I_k$. Every pure strategy profile in which all agents play the same $I_k$ is a Nash equilibrium. The game has $m$ pure Nash equilibria.

**Proposition 15.1 (Institution as Nash Equilibrium).** Every pure-strategy Nash equilibrium of a coordination game corresponds to an institution in the Aoki sense: it is a self-sustaining behavioral regularity in which no agent benefits from unilateral deviation.

This mapping between institutions and Nash equilibria has two immediate implications. First, the number of possible institutions equals the number of Nash equilibria of the underlying game — typically large, which explains the diversity of institutions observed across societies facing similar fundamental coordination problems. Second, institutional change is the problem of transitioning between equilibria — not merely adjusting an equilibrium parameter, but moving from one basin of attraction to another.

### 15.3.2 Replicator Dynamics and Institutional Selection

When do efficient institutions replace inefficient ones? The replicator dynamics of Chapter 5 and Chapter 7 provide the formal answer.

**Definition 15.5 (Institutional Replicator Dynamics).** Let $x_k(t)$ be the fraction of agents following institution $I_k$ at time $t$. The institutional replicator equation is:

$$\dot{x}_k = x_k \left[f_k(\mathbf{x}) - \bar{f}(\mathbf{x})\right]$$

where $f_k(\mathbf{x}) = \sum_l x_l u(I_k, I_l)$ is the fitness of institution $I_k$ given the current institutional mix $\mathbf{x} = (x_1, \ldots, x_m)$, and $\bar{f}(\mathbf{x}) = \sum_k x_k f_k(\mathbf{x})$ is the population mean fitness.

For the pure coordination game above:

$$f_k(\mathbf{x}) = x_k \cdot v(I_k) \quad (\text{fitness only when meeting a coordinated partner})$$

$$\dot{x}_k = x_k \left[x_k v(I_k) - \sum_l x_l^2 v(I_l)\right]$$

**Theorem 15.1 (Basin-of-Attraction Dominance).** In a population of $n$ agents playing the institutional coordination game under replicator dynamics, the institution $I^*$ with the largest basin of attraction — the largest set of initial conditions $\mathbf{x}(0)$ from which the dynamics converge to the all-$I^*$ equilibrium — dominates in the long run.

*Proof.* The replicator dynamics has $m$ stable fixed points corresponding to the $m$ pure Nash equilibria (each all-$I_k$ population) and $m-1$ unstable mixed fixed points (one between each pair of adjacent pure equilibria). The boundary between basins of attraction passes through these unstable fixed points. Any random perturbation starting within the basin of $I^*$ converges to $I^*$; any perturbation starting outside it converges to a different equilibrium. The long-run invariant distribution under random perturbations therefore concentrates mass on the equilibrium with the largest basin. $\square$

**Corollary 15.1 (Efficiency Does Not Guarantee Selection).** The institution with the largest coordination value $v(I^*)$ is not necessarily the one with the largest basin of attraction. An institution with lower coordination value but an initial frequency advantage — perhaps because it was adopted first by a large community — may have a larger basin and will therefore dominate in the long run despite being less efficient.

This is the formal statement of institutional lock-in: the evolutionary process selects for historical incumbency as much as for efficiency. QWERTY keyboard layouts, VHS tapes, internal combustion engines, and the US customary system of measurement are all institutions with smaller basins of attraction than their more efficient alternatives, but which achieved sufficient initial penetration to lock in. The implication for cooperative institutional design is direct: introducing a more efficient institutional arrangement is insufficient — it must achieve sufficient initial adoption to push the system past the basin boundary.

---

## 15.4 Stigmergic Institutional Formation

### 15.4.1 Institutions as Signal-Behavior Patterns

Chapter 7 introduced stigmergy as the coordination mechanism through which agents modify a shared environment in ways that structure subsequent behavior. We now apply this framework to institutional formation: institutions emerge when behavioral patterns leave environmental traces that reinforce themselves.

**Definition 15.6 (Stigmergic Institution).** A stigmergic institution is a stable fixed point $(\Sigma^*, \mathbf{a}^*)$ of the stigmergic signal dynamics (Definition 7.3 [C:Ch.7]):

$$\dot{\Sigma} = f(\mathbf{a}) - \delta\Sigma, \quad a_i = g_i(\Sigma)$$

such that:
(i) The signal $\Sigma^*$ represents a behavioral pattern (an institutional regularity);
(ii) Agents' best responses $g_i(\Sigma^*)$ reproduce the signal: $f(\mathbf{a}^*) = \delta\Sigma^*$;
(iii) The fixed point is stable: small perturbations to $\Sigma$ or $\mathbf{a}$ return to $(\Sigma^*, \mathbf{a}^*)$.

The key feature of a stigmergic institution is that it is self-maintaining through the mutual reinforcement of signals and behaviors: agents respond to the signal by taking the institutionally appropriate action, which generates the signal, which induces the action. The institution persists as long as the signal-response loop is stable.

**Example 15.1 (Property Rights as Stigmergic Institution).** Consider a community sharing a common-pool resource. Initially, the resource is ungoverned (open access). Individual agents begin marking the patches they regularly use — physically (fences, markers) or socially (verbal claims, witnessed occupation). Each marker is a stigmergic signal: it tells subsequent agents "this patch is claimed." Agents who observe the signal respond by avoiding the claimed patch (to avoid conflict). Their avoidance behavior reinforces the original claimant's exclusive use, which reinforces the signal, completing the loop.

The property rights institution crystallizes when the signal-response pattern reaches a stable fixed point: all patches are marked, all agents respect the marks, and the community has collectively evolved an institution that no one explicitly designed.

**Definition 15.7 (Institutional Crystallization).** A stigmergic institution crystallizes when the fraction of agents following the institutionally appropriate response to the signal exceeds the critical threshold $\bar{x}$ above which the signal is self-sustaining:

$$x(t) \geq \bar{x} \implies \Sigma(t) \to \Sigma^* \text{ (institution crystallizes)}$$
$$x(t) < \bar{x} \implies \Sigma(t) \to 0 \text{ (institution dissolves)}$$

The crystallization threshold $\bar{x}$ is the institutional analogue of the tipping point of Chapter 5: above it, the institutional signal is self-reinforcing; below it, it decays.

### 15.4.2 The Crystallization Threshold

**Proposition 15.2 (Crystallization Threshold).** For the linear stigmergic dynamics $\dot{\Sigma} = \alpha x \Sigma - \delta\Sigma$ (where $\alpha$ is the signal-reinforcement rate and $x$ is the fraction of agents following the institution), the crystallization threshold is:

$$\bar{x} = \frac{\delta}{\alpha}$$

Above $\bar{x}$, the signal dynamics have a stable fixed point at $\Sigma^* = (\alpha x - \delta)/\delta \cdot \Sigma_0 > 0$. Below $\bar{x}$, the only stable fixed point is $\Sigma = 0$ (no institution).

*Proof.* At steady state $\dot{\Sigma} = 0$: $\alpha x \Sigma^* = \delta\Sigma^*$, giving $x = \delta/\alpha \equiv \bar{x}$ as the bifurcation point. For $x > \bar{x}$: $\alpha x > \delta$ and any positive $\Sigma$ grows toward the non-zero fixed point. For $x < \bar{x}$: $\alpha x < \delta$ and $\Sigma$ decays to zero. $\square$

**Economic interpretation.** The crystallization threshold $\bar{x} = \delta/\alpha$ depends on the ratio of signal decay rate to signal reinforcement rate. Institutions that generate strong behavioral signals (high $\alpha$) — visible, salient, easy to read and respond to — have low crystallization thresholds and emerge easily. Institutions that generate weak signals (low $\alpha$) or decay rapidly (high $\delta$) require a larger fraction of initial adopters before they can become self-sustaining.

This explains why some institutions crystallize rapidly from small initial conditions (language conventions, handshake greetings, traffic rules) while others require sustained collective effort to establish (new property right systems, novel governance forms, unfamiliar monetary arrangements): the former have high $\alpha$ (highly visible, easy to imitate signals) while the latter have low $\alpha$ or high $\delta$ (complex, subtle, easily forgotten signals).

---

## 15.5 The Role of Entrepreneurs and Crises

### 15.5.1 Institutional Entrepreneurship

Institutional change does not always wait for evolutionary selection or spontaneous stigmergic crystallization. Some agents — institutional entrepreneurs — deliberately invest resources in shifting the institutional equilibrium, changing the payoff structure of the coordination game, or constructing new stigmergic signals to bootstrap a new institution past the crystallization threshold.

**Definition 15.8 (Institutional Entrepreneur).** An institutional entrepreneur is an agent $i^*$ who chooses to invest resources $k_{i^*} \geq 0$ in institutional change, with the objective:

$$\max_{k_{i^*}} B(x(k_{i^*})) - C(k_{i^*})$$

where $B(x)$ is the entrepreneur's benefit from moving the institutional frequency from $x_0$ to $x(k_{i^*})$ (the new frequency achievable with investment $k$), and $C(k)$ is the cost of the investment.

The entrepreneur's problem is to find the minimum investment $k^*$ such that $x(k^*) \geq \bar{x}$ — sufficient to push the institution past the crystallization threshold. Below this threshold, the investment is wasted (the institution dissolves when the entrepreneur stops investing). Above it, the institution becomes self-sustaining and the entrepreneur can stop investing.

**Proposition 15.3 (Minimum Viable Investment).** The minimum investment to crystallize a new institution with initial frequency $x_0 < \bar{x}$ and crystallization threshold $\bar{x}$ is:

$$k^* = C^{-1}\left(\frac{\bar{x} - x_0}{\eta}\right)$$

where $\eta > 0$ is the marginal effectiveness of investment in raising institutional frequency and $C^{-1}$ is the inverse cost function. The investment is socially optimal if $B(1) > k^*$ — if the total benefit of the new institution exceeds the minimum investment cost.

This result has a direct implication for cooperative institutional design: cooperatives, municipalities, and activist communities that want to promote new institutional arrangements should focus their resources on pushing the new institution past its crystallization threshold, not on persuading the entire population. Once the threshold is crossed, the institution becomes self-sustaining and further investment is unnecessary.

### 15.5.2 Crises as Institutional Opportunities

The coordination game model implies that switching from one institutional equilibrium to another requires perturbing the system out of the incumbent equilibrium's basin of attraction. Under normal conditions, this requires either large investment (institutional entrepreneurship) or rare large shocks. Crises — economic collapses, wars, pandemics, ecological disasters — provide precisely the large perturbations that can dislodge incumbent institutions.

**Definition 15.9 (Crisis-Induced Institutional Window).** A crisis opens an institutional window of duration $T_{window}$ during which:
(i) The payoff to the incumbent institution is temporarily reduced: $v(I_{\text{incumbent}})$ falls by $\Delta v > 0$.
(ii) The cost of adopting an alternative institution is temporarily reduced: $C(I_{\text{new}})$ falls by $\Delta C > 0$.
(iii) The basin of attraction of the incumbent institution is temporarily contracted.

**Proposition 15.4 (Crisis and the Institutional Window).** A crisis opens an institutional window if and only if:

$$\Delta v + \Delta C > v(I_{\text{incumbent}}) - v(I_{\text{new}})$$

— the combined reduction in incumbent payoff and adoption cost exceeds the initial institutional advantage of the incumbent.

*Proof.* Under normal conditions, the incumbent is selected because $v(I_{\text{incumbent}}) > v(I_{\text{new}})$ and $C(I_{\text{new}}) > 0$. A crisis reduces the incumbent payoff by $\Delta v$ and adoption cost by $\Delta C$. The window opens when the net advantage of the incumbent is reversed: $(v(I_{\text{incumbent}}) - \Delta v) - (v(I_{\text{new}}) - 0) - C(I_{\text{new}}) + \Delta C < 0$, which simplifies to the condition given. $\square$

Historical examples: the Great Depression opened an institutional window for New Deal financial regulation; the 2008 financial crisis opened a window for unconventional monetary policy; the COVID-19 pandemic opened windows for universal healthcare advocacy, remote work norms, and accelerated vaccine development governance. In each case, the crisis reduced the perceived payoff to existing institutions and reduced the adoption costs of alternatives — sometimes temporarily, sometimes permanently changing the institutional landscape.

---

## 15.6 Inefficient Institutional Equilibria and Lock-In

### 15.6.1 The Lock-In Condition

Not all emergent institutions are efficient. The basin-of-attraction dominance theorem (Theorem 15.1) established that historical incumbency can override efficiency in institutional selection. The formal condition for lock-in to an inefficient institution is:

**Definition 15.10 (Lock-In).** An institution $I^{\text{bad}}$ is locked-in relative to a superior alternative $I^*$ if:

1. $v(I^*) > v(I^{\text{bad}})$ (the alternative is more efficient);
2. The basin of attraction of $I^{\text{bad}}$ under replicator dynamics contains the current institutional frequency: $\mathbf{x}(0) \in \mathcal{B}(I^{\text{bad}})$;
3. The cost of exiting $\mathcal{B}(I^{\text{bad}})$ — moving to $I^*$'s basin of attraction — exceeds the private benefit of doing so for any single agent.

**Theorem 15.2 (Lock-In Persistence).** An institution $I^{\text{bad}}$ remains locked-in until one of the following conditions is met:
(i) A crisis reduces $v(I^{\text{bad}})$ sufficiently to contract $\mathcal{B}(I^{\text{bad}})$ beyond the current $\mathbf{x}(0)$ (Proposition 15.4).
(ii) An institutional entrepreneur invests sufficient resources to push $\mathbf{x}$ out of $\mathcal{B}(I^{\text{bad}})$ (Proposition 15.3).
(iii) The crystallization threshold $\bar{x}$ of $I^*$ falls below the current frequency of $I^*$ adherents through gradual cultural change (Proposition 15.2, with $\delta$ slowly declining).

*Proof.* By Theorem 15.1, replicator dynamics starting inside $\mathcal{B}(I^{\text{bad}})$ converge to $I^{\text{bad}}$. The only exits are basin boundary crossings, which require either the basin to contract toward the current state (conditions (i) and (iii)) or the state to move across the basin boundary (condition (ii)). $\square$

### 15.6.2 Escaping Bad Equilibria: Coordination Mechanisms

When a society is locked into an inefficient institution, individual agents face a collective action problem: transitioning to the better institution requires many agents to change simultaneously, but each agent benefits from the transition only if others change first. This is the institutional version of the chicken-and-egg problem.

Three coordination mechanisms can resolve this:

**Mechanism 1 (Salience and focal points).** If one alternative institution is cognitively salient — it stands out from alternatives in a way that makes it an obvious coordination point — agents may coordinate on it without explicit communication (Schelling, 1960). The formal condition: institution $I^*$ is a focal point if $\Pr[\text{agent } i \text{ plays } I^* \mid \text{agent } j \text{ plays } I^*] > \bar{x}$ for all $i, j$ in the relevant community.

**Mechanism 2 (Commitment devices).** If agents can credibly commit to adopting the new institution before others do, they can bootstrap the coordination. A commitment device is a costly action that makes defection from the new institution more expensive than maintaining it. Formal legal commitments, public declarations, and reputational investments serve this function.

**Mechanism 3 (Sequenced adoption).** If the community can identify a subgroup whose adoption of the new institution alone pushes the frequency above the crystallization threshold $\bar{x}$, sequencing adoption through this critical mass can achieve system-wide transition at minimum cost. The critical mass is the minimum coalition $S^*$ such that $|S^*|/n \geq \bar{x}$.

---

## 15.7 Mathematical Model: Property Rights Emergence in a Commons

We now model the emergence of property rights in a commons as a coordination game and trace the replicator dynamics, showing how different initial conditions lead to either the property rights institution or the open-access institution.

**Setup.** A commons of $n$ agents shares a resource $R$ with carrying capacity $K$. Each agent can adopt one of two institutional strategies:

- **Open Access ($I_0$):** Extract without regard to others' claims. Payoff: $v_0(x_0) = A - bx_0 n \bar{e}$ (the Cournot payoff declining in total extraction by $x_0 n$ open-access users).

- **Property Rights ($I_1$):** Respect a system of territorial claims; extract only from one's designated patch. Payoff: $v_1(x_1) = A - bx_1 n e^* - (1-x_1)n\bar{e}$ where $e^* < \bar{e}$ is the sustainable extraction level under property rights.

**Payoff structure.** The payoff from property rights exceeds open access when resource depletion under open access is severe:

$$v_1 > v_0 \iff A - bx_1 n e^* - (1-x_1)n\bar{e} > A - bx_0 n\bar{e}$$

For $x_1 = x_0 = x$ (equal frequencies transitioning simultaneously):

$$v_1 > v_0 \iff b(\bar{e} - e^*)n(1-x) > 0$$

which always holds when $\bar{e} > e^*$ (sustainable extraction is lower than open-access extraction) — property rights always dominate in a homogeneous population. The problem is the transition: at $x = 0$ (all open-access), a single agent adopting property rights earns only a small coordination benefit while bearing the full cost of respecting others' claims (which they don't reciprocate).

**Replicator dynamics.** Let $x = x_1$ (the fraction using property rights). The fitness difference is:

$$f_1(x) - f_0(x) = bn(\bar{e} - e^*)(x - \hat{x})$$

where $\hat{x}$ is the unstable mixed equilibrium (the basin boundary):

$$\hat{x} = \frac{(1-x)\bar{e}n \cdot b - (b n e^*) \cdot x}{bn(\bar{e}-e^*)} \approx \frac{\bar{e}}{\bar{e} + e^*}$$

For $\bar{e} = 2e^*$ (open-access extraction is double sustainable): $\hat{x} = 2/3$. The property rights institution requires two-thirds of the population to adopt before it becomes self-sustaining under replicator dynamics.

**Phase portrait.** The replicator dynamics:
$$\dot{x} = x(1-x) \cdot bn(\bar{e}-e^*)(x - \hat{x})$$

has three fixed points: $x=0$ (stable, all open-access), $x=\hat{x}$ (unstable, the tipping threshold), and $x=1$ (stable, all property rights). The dynamics converge to:
- $I_0$ (open access) if $x(0) < \hat{x}$
- $I_1$ (property rights) if $x(0) > \hat{x}$

**When does the "wrong" institution win?** If the community begins with fewer than two-thirds of members willing to respect property claims — perhaps because the commons was previously ungoverned and most agents are habituated to open access — the replicator dynamics converge to the open-access equilibrium even though property rights are more efficient. This is institutional lock-in in its simplest form: the historically prior institution (open access) prevents the adoption of the superior alternative (property rights) through the dynamics of coordination failure.

**The institutional entrepreneur's role.** An entrepreneur who can credibly signal property rights adoption to a critical mass — through physical markers, community organizing, or legal registration — can push $x(0)$ above $\hat{x}$ and tip the system toward the property rights equilibrium. The minimum investment is proportional to $\hat{x} - x(0)$, the distance from the current frequency to the tipping threshold.

---

## 15.8 Case Study: The Lex Mercatoria as a Stigmergic Institution

### 15.8.1 Medieval Commerce Without State Enforcement

Between approximately the 11th and 16th centuries, a pan-European system of merchant law — the lex mercatoria or "law merchant" — governed commercial transactions across political boundaries in which no common state authority existed and no single legal system applied. Merchants from Genoa, Venice, Bruges, Hamburg, and London transacting at the great fairs of Champagne, the North Sea ports, and the Baltic markets needed mechanisms to enforce contracts, resolve disputes, and certify the quality of goods — all without recourse to any state.

How did they manage? The conventional answer — that merchants developed their own courts, their own norms, and their own enforcement mechanisms — is correct but incomplete. The formal question is: what made the lex mercatoria a stable institution without the coercive authority that standard contract theory assumes is necessary?

### 15.8.2 The Lex Mercatoria as Stigmergic Equilibrium

**The signals.** Medieval merchants generated and responded to several classes of stigmergic signals:

*Merchant reputation signals:* Letters of credit, bills of exchange, and commercial correspondence created a paper trail of reputation. A merchant who defaulted left a documentary trace — in the merchant's own ledger (which third parties might inspect), in the records of merchant courts, and in the communal memory of the merchant guilds and trading communities. This trace persisted long after the original transaction and was readable by subsequent parties.

*Quality certification signals:* Guild marks, standardized measures (the Flemish ell, the pound troy), and the physical characteristics of certified goods created signals of product quality that agents could read without direct knowledge of the seller's history. A bolt of cloth bearing the seal of the Cloth Merchants' Guild of Florence carried information about its quality that required no bilateral relationship to interpret.

*Community membership signals:* Membership in a merchant guild, registration in a town's commercial records, and visible participation in the commercial fairs were signals of accountability — they indicated that a merchant was embedded in the social network within which reputation mattered and sanctions were enforceable.

**The response functions.** Merchants responded to these signals by: (i) extending credit to agents with strong reputation signals and withholding it from those with weak signals; (ii) accepting guild-certified goods at standard prices and discounting uncertified goods; and (iii) trading preferentially with other guild members, whose institutional embeddedness made default costly.

**Formal reconstruction.** The lex mercatoria is the stigmergic equilibrium of the medieval merchant game:

$$\Sigma^* = \{\text{reputation records, quality marks, guild membership certificates}\}$$
$$\mathbf{a}^* = \{\text{extend credit to good reputation, reject bad; accept certified quality, discount uncertified}\}$$

satisfying the fixed-point conditions: $f(\mathbf{a}^*) = \delta\Sigma^*$ (merchants' actions generate the signals at exactly the rate at which they decay) and $a_i = g_i(\Sigma^*)$ (each merchant's best response to the observed signals is the institutionally appropriate action).

### 15.8.3 Stability Without State Enforcement

What made the lex mercatoria stable in the absence of state coercive authority? Four features of the signal-response architecture:

**1. Low signal decay rate.** Reputation information was recorded in ledgers and court records that persisted for decades. $\delta$ was low — signals from past transactions continued to influence future ones for long periods, extending the shadow of the future and making defection costly.

**2. High signal visibility.** Physical guild marks, standardized weights and measures, and the social visibility of fair participation made signals immediately readable without investigation. $\alpha$ (the signal reinforcement rate) was high — compliant behavior generated clear, visible signals that attracted future business.

**3. Community closure.** The merchant community was small enough (a few thousand major merchants in any regional network) and interconnected enough that reputation information spread rapidly and reliably. This is the small-world property of Chapter 4 applied to reputation networks: high clustering within merchant communities (strong local reputation) combined with short path lengths to other communities (rapid spread of reputation information across communities).

**4. Graduated sanctions.** The sanction for default was not sudden expulsion but gradual exclusion: reduced credit, higher collateral requirements, exclusion from preferred trading partnerships. This is DP5 (graduated sanctions) of the Ostrom framework operating without formal legal authority — sanctions administered by the commercial community itself.

The lex mercatoria collapsed not because the stigmergic equilibrium became unstable but because the state absorbed and formalized it: as national commercial codes developed in the 16th and 17th centuries, the informal merchant law was replaced by state-enforced contract law with higher $\delta$ (formal records decay more slowly — legal records are permanent) but lower $\alpha$ (legal proceedings are slow and cumbersome, reducing the feedback speed of the signal loop). The institution was absorbed into the state magisterium precisely because state enforcement could replicate some but not all of the lex mercatoria's properties.

---

## Chapter Summary

This chapter has developed the formal theory of institutional emergence, connecting the evolutionary game theory of Chapter 7, the stigmergic coordination model, and the complexity economics of Chapter 5 into a unified account of how institutions arise from interaction.

Institutions are defined in three complementary ways — as rules (North), as rules-in-use (Ostrom), and as self-sustaining equilibria (Aoki) — each capturing a different aspect of institutional reality. The Aoki definition is the most amenable to game-theoretic analysis: an institution is a Nash equilibrium of the repeated interaction game, stable because it is in everyone's interest to maintain it given that others are doing so.

The institutional replicator dynamics selects institutions by their basin of attraction rather than their efficiency: the historically prior institution may dominate even when a more efficient alternative exists (Theorem 15.1, Corollary 15.1). This is the formal foundation of institutional lock-in.

Stigmergic institutional formation models institutions as stable signal-behavior patterns. The crystallization threshold $\bar{x} = \delta/\alpha$ is the minimum adoption frequency above which the institution becomes self-sustaining — below it, the institutional signal decays and the institution dissolves. High-signal-strength institutions (high $\alpha$) and low-decay institutions (low $\delta$) crystallize more easily from smaller initial frequencies.

Crises reduce incumbent payoffs and adoption costs simultaneously, opening institutional windows when the combined reduction exceeds the incumbent's advantage (Proposition 15.4). Institutional entrepreneurs who can push frequency past the crystallization threshold enable system-wide institutional change at minimum cost.

The lex mercatoria exemplifies all these dynamics: a stigmergic institution that crystallized from the repeated commercial interactions of medieval merchants, sustained by low signal decay, high signal visibility, community closure, and graduated sanctions — stable for five centuries without state enforcement, then absorbed into the state magisterium when formal contract law could replicate its coordination function at lower operational cost.

Chapter 16 closes Part III with the formal economics of information asymmetry in networks — the theory of how trust, reputation, and reciprocity overcome the problems of adverse selection and moral hazard in economic relationships where information is unequally distributed.

---

## Exercises

**15.1** Define an institution formally as a Nash equilibrium (Definition 15.3 / Aoki). Give two economic examples, specifying for each: the players, the strategy space $\mathcal{I}$, and the payoff structure that makes the institution self-sustaining.

**15.2** In the property rights emergence model (Section 15.7):
(a) Compute the tipping threshold $\hat{x}$ for $\bar{e} = 3e^*$ (open-access extraction is triple sustainable). How does this compare to the $\bar{e} = 2e^*$ case?
(b) Under what conditions does the tipping threshold approach 0 (property rights always emerge)? Approach 1 (property rights never emerge)?
(c) An institutional entrepreneur can move $x(0)$ from 0.40 to 0.55 by investing in community organizing. If $\hat{x} = 0.50$, is this investment sufficient? What is the minimum investment to guarantee emergence regardless of $\hat{x}$ within the range $[0.45, 0.60]$?

**15.3** The lex mercatoria was stable for approximately 500 years without state enforcement. Analyze it as a stigmergic institution:
(a) Identify the signal vector $\Sigma^*$, the response functions $g_i$, and the signal generation function $f(\mathbf{a})$. Specify the signal decay rates $\delta_j$ for each signal type.
(b) Estimate the crystallization threshold $\bar{x} = \delta/\alpha$ for the reputation signal component. What does this imply about the minimum size of the initial merchant community needed to sustain the institution?
(c) When national commercial codes began replacing the lex mercatoria in the 16th century, did this represent: (i) a crisis-induced institutional window; (ii) an institutional entrepreneur investment; or (iii) a gradual reduction in the crystallization threshold? Justify your answer using the formal framework.

**★ 15.4** Prove Theorem 15.1: in a population of $n$ agents playing the institutional coordination game under replicator dynamics, the institution with the larger basin of attraction dominates in the long run.

(a) For a two-institution coordination game with coordination values $v_1 > v_2 > 0$ but initial frequencies $x_1(0) < x_2(0)$, write out the replicator equation and find all fixed points.
(b) Identify the unstable mixed equilibrium $\hat{x}$ (the basin boundary) and show it depends on $v_1$ and $v_2$.
(c) Prove that if $x_1(0) < \hat{x}$, then $x_1(t) \to 0$ even though $v_1 > v_2$. Interpret this as institutional lock-in.
(d) For a population undergoing random drift (each period, a fraction $\varepsilon$ of agents randomly switch institutions), compute the stationary distribution over the two institutional equilibria. Show that for small $\varepsilon$, the stationary distribution concentrates on the institution with the larger basin of attraction rather than the higher coordination value.

**★ 15.5** Model the 2008 financial crisis as an institutional window.

(a) Identify the incumbent institution $I^{\text{bad}}$ (the pre-2008 financial regulation regime), the alternative $I^*$ (alternative candidates: narrow banking, Tobin tax, higher capital requirements, public banking), and the coordination values $v(I^{\text{bad}})$ and $v(I^*)$.
(b) Using Proposition 15.4, compute the change in incumbent payoff $\Delta v$ and adoption cost reduction $\Delta C$ implied by the crisis. Did the crisis open an institutional window for each of the four alternatives?
(c) The post-2008 regulatory response was primarily Dodd-Frank (higher capital requirements) rather than more radical alternatives. Using the basin-of-attraction framework, explain why the more moderate alternative was selected even if a more radical alternative was more efficient.
(d) What would have been required — in terms of institutional entrepreneur investment or crisis severity — for the window to have produced a more fundamental institutional change (e.g., full narrow banking)?

**★★ 15.6** Implement an ABM of institutional emergence for 500 agents on a network.

**Model specification:**
- 500 agents, each choosing between two institutional strategies $I_0$ (open access) and $I_1$ (property rights).
- Network: Watts-Strogatz small-world with $n=500$, $\bar{k}=8$, $\beta=0.15$.
- Payoffs: $f_0(x) = A - bxn\bar{e}$, $f_1(x) = A - bxne^* - b(1-x)n\bar{e}$ with $A=10$, $b=0.01$, $\bar{e}=0.3$, $e^*=0.15$.
- Strategy update: Fermi learning rule (agents copy successful neighbors with probability $1/(1 + e^{-\beta_l(w_j-w_i)})$, $\beta_l=0.1$).
- Initial conditions: vary $x_1(0)$ from 0.1 to 0.9 in steps of 0.1.

(a) Run 100 periods for each initial condition across 20 replications. At each initial condition, record whether the simulation converges to $I_0$, $I_1$, or remains mixed.

(b) Estimate the empirical tipping threshold $\hat{x}^{\text{ABM}}$ and compare to the analytical threshold $\hat{x}$ derived in Section 15.7. Are they consistent? If not, explain the discrepancy in terms of network structure.

(c) Run the same experiment on a Barabási-Albert scale-free network with $m=4$. How does the tipping threshold change? Connect to the clustering condition for ESS (Proposition 7.2 [C:Ch.7]) and the hub-vulnerability result of Chapter 12.

(d) Introduce an institutional entrepreneur: at period 0, a single high-degree hub agent (the most connected agent in the network) commits to $I_1$ regardless of payoffs. Does this shift the empirical tipping threshold? By how much? Interpret this in terms of Proposition 15.3 (minimum viable investment).

(e) Add a crisis at period 50: the open-access payoff drops temporarily by 30\% for 10 periods. How does this affect institutional trajectories starting from $x_1(0) = 0.35$ (below the analytical threshold)? Does the crisis open an institutional window (Proposition 15.4), and if so, which institution emerges after the window closes?

---

*Chapter 16 closes Part III with the formal economics of trust, reputation, and reciprocity in networks — the mechanisms through which information asymmetry is overcome in economic relationships that cannot rely on formal contract enforcement. The theoretical tools connect directly back to Chapter 7's evolutionary stability analysis and forward to the design of cooperative institutions in Parts V and VI.*
