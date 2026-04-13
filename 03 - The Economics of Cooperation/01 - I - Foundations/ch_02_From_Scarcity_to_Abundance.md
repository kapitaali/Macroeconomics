# Chapter 2: From Scarcity to Abundance — Rethinking Resources and Value

> *"The ideas of economists and political philosophers, both when they are right and when they are wrong, are more powerful than is commonly understood. Indeed the world is ruled by little else."*
> — John Maynard Keynes, *The General Theory* (1936)

> *"An economist is someone who, when they see something working in practice, wonders whether it could also work in theory."*
> — attributed, after Ronald Coase

## Learning Objectives

By the end of this chapter, you should be able to:

1. Distinguish between rival and non-rival goods, and between excludable and non-excludable goods, and explain the economic consequences of each classification.
2. Explain formally why knowledge, software, and other digital goods violate the scarcity postulate and why standard pricing theory breaks down in their presence.
3. State the Provisioning Question and explain why it reframes the central problem of economics from allocation to maintenance.
4. Define the Stewardship Objective Function and derive its first-order conditions.
5. Describe the Three Coordination Engines — markets, hierarchies, and mutual coordination — and explain the information-processing logic of each.
6. State Ostrom's eight design principles for commons governance and articulate the conditions under which commons management outperforms both private and state alternatives.

---

## 2.1 The Scarcity Postulate and Its Domain

Lionel Robbins, in his 1932 *Essay on the Nature and Significance of Economic Science*, offered the definition of economics that became canonical for the twentieth century: economics is "the science which studies human behaviour as a relationship between ends and scarce means which have alternative uses." For Robbins, scarcity was not merely an empirical observation about the current state of technology and resources; it was the constitutive condition of economic life, the fact that makes economics necessary as a discipline.

This definition has served the field well. It provides a clear object of study, a natural optimization framework, and a rationale for the price mechanism as a device for allocating scarce goods to their highest-valued uses. It also carries a hidden assumption that has become increasingly consequential as the structure of the economy has shifted: it treats scarcity as a feature of the goods themselves rather than as a property of the institutional arrangements under which they are produced and distributed.

Consider three goods: a kilogram of wheat, a river's worth of salmon, and a piece of software. All three satisfy real human needs. The wheat is unambiguously scarce: if I eat it, you cannot. The salmon are scarcer than they once were, but whether they remain scarce in the future depends critically on how many boats fish the river, how the fishery is governed, and whether the river's ecosystem is maintained — that is, it depends on institutions, not on the physical properties of fish. The software is not scarce in any meaningful physical sense: if I use a copy of an operating system, your copy is unaffected. Its "scarcity" is entirely the product of intellectual property law — an institutional choice, not a natural constraint.

The wheat, the salmon, and the software thus inhabit quite different economic worlds, and the analytical tools appropriate to each differ correspondingly. A single framework built on the scarcity postulate will handle wheat well, handle salmon contingently (depending on governance), and systematically mishandle software — producing pricing recommendations that are welfare-reducing, and a policy orientation toward artificial restriction of what could, at zero marginal cost, be universally available.

This chapter develops the conceptual vocabulary for distinguishing these cases, introduces a positive theory of the commons as a governance mechanism suited to common-pool resources like salmon, and argues that the appropriate response to non-rival goods like software is neither private pricing nor state provision but a third mode of social organization: the commons-based peer production that has, since the 1990s, built much of the world's digital infrastructure.

### 2.1.1 The Goods Taxonomy

The classical taxonomy of goods, due in its formal version to Paul Samuelson (1954) and later extended by Ostrom and colleagues, rests on two independent dimensions:

**Definition 2.1 (Rivalry).** A good is *rival* if consumption of one unit by one agent reduces the quantity available to others. It is *non-rival* if consumption by one agent does not diminish availability for others.

Formally, good $g$ is non-rival if, for any two agents $i$ and $j$ and any consumption levels $x_i, x_j \geq 0$:

$$\frac{\partial^2 u_i}{\partial x_i \partial x_j} = 0 \quad \text{and} \quad \frac{\partial u_j}{\partial x_i} = 0$$

That is, the consumption of agent $i$ has no effect on either the utility or the marginal utility of the good for agent $j$. This is the condition that the good's benefits are genuinely non-subtractive.

**Definition 2.2 (Excludability).** A good is *excludable* if it is technically and legally feasible to prevent any agent from consuming it. It is *non-excludable* if preventing consumption is infeasible or prohibitively costly.

The combination of these two properties generates a fourfold taxonomy that is worth understanding precisely, because the appropriate institutional response differs in each quadrant:

| | **Excludable** | **Non-excludable** |
|:---|:---|:---|
| **Rival** | Private goods *(bread, land, cars)* | Common-pool resources *(fish stocks, groundwater, atmosphere)* |
| **Non-rival** | Club goods *(toll roads, cable TV, patents)* | Pure public goods *(national defense, open science, climate stability)* |

Two features of this table deserve immediate emphasis. First, the categories are not fixed by nature; they are partly determined by technology and institutions. Encryption technology can transform a non-excludable digital good into a club good. Overfishing governance can transform a common-pool resource into something closer to a private good (through individual transferable quotas) or a commons (through collective management). The relevant question is always not merely "what kind of good is this?" but "what institutional arrangement is appropriate given its physical properties?"

Second, the standard market mechanism is designed for private goods — rival and excludable — and its performance degrades as goods deviate from this archetype. Markets underprovide public goods (as we saw in Chapter 1), handle common-pool resources badly without additional governance, and, as we argue below, generate welfare-reducing outcomes for non-rival goods whenever they are made artificially excludable through intellectual property rights.

### 2.1.2 Non-Rival Goods and the Digital Economy

The digital revolution has made non-rival goods the fastest-growing sector of the global economy. Software, scientific knowledge, recorded music, film, databases, trained machine-learning models, and the accumulated text of human culture: all of these are non-rival in the strict sense of Definition 2.1. Once produced, they can be reproduced and transmitted at near-zero marginal cost.

This creates a fundamental tension with standard pricing theory. The efficiency condition for any good is that price equals marginal cost:

$$p^* = MC$$

For a non-rival good with zero marginal reproduction cost, this condition requires $p^* = 0$. But if the good was costly to produce — which most software, scientific research, and creative works are — then pricing at zero marginal cost means the producer cannot recover fixed production costs and will not produce the good in the first place. We face a genuine dilemma: static efficiency (price equal to marginal cost) is incompatible with dynamic efficiency (sufficient incentive to produce) for any non-rival good with positive fixed costs [P:Ch.36].

The standard policy response — intellectual property (IP) law — resolves this dilemma by creating artificial scarcity. Patents, copyrights, and trade secrets make otherwise non-excludable goods excludable, allowing producers to charge above marginal cost and recover fixed costs. The welfare loss from this artificial scarcity is the deadweight loss of monopoly pricing: the consumer surplus destroyed by the gap between price and marginal cost. In the worked example of Section 2.6, we show that for a good with truly zero marginal cost, this deadweight loss equals the entire consumer surplus — a result that motivates the analysis of open-source and commons-based alternatives.

The digital economy does not abolish scarcity; it relocates it. Physical infrastructure (data centers, network cables, energy) remains rival. Human attention remains rival. The time of skilled software engineers remains rival. What the digital economy does is separate, for the first time in economic history, the non-rival information component of production from the rival physical and human components. This separation creates new possibilities — a piece of software can simultaneously serve millions of users at no additional cost — and new problems — the institutions designed around the scarcity of physical goods are poorly adapted to the abundance of digital goods.

---

## 2.2 The Provisioning Question

The discussion of scarcity in Section 2.1 is largely a discussion of *allocation*: how should scarce goods be distributed among competing uses? This is the central question of standard economics, and it is an important one. But it is not the most fundamental question for an economy that must sustain human life and welfare over time.

The most fundamental question is one of *maintenance*: which assets are required for long-term human provisioning — and are we maintaining them?

We call this the Provisioning Question, and its implications reach far beyond the standard framework. Allocation asks how to distribute a given stock of goods today. The Provisioning Question asks whether the stock itself — the productive capacity of the economy, including its natural capital, its institutions, its knowledge base, and its social fabric — is being maintained or depleted.

The shift from allocation to maintenance has three immediate consequences for economic analysis.

**First**, the appropriate unit of welfare analysis changes. If the question is allocation, the natural welfare metric is utility from current consumption: $U(C_t)$. If the question is provisioning, the natural metric includes the flow of services from all relevant capital stocks:

$$U(C_t, S_t)$$

where $S_t$ represents provisioning services — the capacity to sustain human welfare that the current asset base provides. These are not the same thing. A community can sustain high current consumption $C_t$ while degrading $S_t$ through soil depletion, aquifer overextraction, or deforestation — exactly as a household can sustain high current expenditure by liquidating savings. The household in this situation is not prosperous; it is insolvent.

**Second**, the appropriate planning horizon changes. Allocation problems can often be analyzed period by period, with future periods discounted at the rate of time preference. The Provisioning Question is inherently long-run: we are asking whether the economy will remain capable of sustaining human life not for the next few years but for the indefinite future. This requires a treatment of natural capital that does not discount its value to zero in finite time.

**Third**, the appropriate governance framework changes. Allocation is typically governed through ownership and price: private owners respond to price signals by allocating resources to their highest-valued uses. Maintenance is typically governed through stewardship: the commitment of custodians to preserve an asset for future use even when current liquidation might be more profitable. The economics of stewardship is not the economics of allocation; it requires a different theory of institutions, a different theory of value, and — as we will see — a different coordination mechanism.

The Provisioning Question is, in this sense, not a supplement to the standard economic framework but a reorientation of it. It shifts the central problem from "how do we allocate scarce goods?" to "how do we maintain the productive capacity that generates the goods we need?" This reorientation is the pivot on which the rest of this book turns.

---

## 2.3 Social and Ecological Value: A Multi-Dimensional Framework

Standard economics recognizes one dimension of value: exchange value, expressed as a market price. A good is worth what someone will pay for it in a voluntary transaction. This operationalization has practical virtues — it is observable, it aggregates individual preferences, and it provides a common metric for comparison — but it also has systematic blind spots.

Exchange value captures only the value that can be made excludable and transactable. It misses:

- **Use value** that is not traded: the value of a home garden, care work, community volunteering, or self-produced food.
- **Social value** that is non-rival and non-excludable: the value of a culture, a language, a body of shared knowledge, a functional democracy.
- **Ecological value** that has no market: the value of a functioning watershed, a stable climate, a biodiversity-rich ecosystem, or a healthy soil microbiome.

We propose a multi-dimensional value framework in which the total value $V$ of a good or service decomposes as:

$$V = \alpha V_{\text{exchange}} + \beta V_{\text{use}} + \gamma V_{\text{social}} + \delta V_{\text{ecological}}$$

where $\alpha, \beta, \gamma, \delta \geq 0$ are weights reflecting the relative importance of each dimension, and $\alpha + \beta + \gamma + \delta = 1$ under a normalization convention.

This is not a formula to be computed with precision; it is a conceptual scaffold. Its purpose is to make visible the dimensions of value that standard accounting systematically omits, and to motivate the institutional arrangements — commons governance, natural capital accounting, open-source production — that can sustain those dimensions without necessarily routing them through markets.

The ecological dimension deserves particular emphasis. The natural systems that sustain human life — the atmosphere, the hydrological cycle, soil fertility, pollinator populations, ocean chemistry — generate a flow of ecosystem services estimated in the trillions of dollars annually by ecological economists (Costanza et al., 1997; 2014). These services do not appear in GDP, are not owned by any agent with an incentive to maintain them, and are not protected by any market mechanism. Their degradation is therefore costless to any individual economic actor, which is precisely why it proceeds so rapidly.

The multi-dimensional value framework does not solve the measurement problem — commensuration across dimensions remains genuinely difficult — but it correctly diagnoses why standard accounting fails to protect natural and social capital: it does not measure what it does not price, and it does not price what has no owner. The institutional implication is that maintaining natural and social capital requires governance mechanisms beyond the market.

---

## 2.4 The Three Coordination Engines

Chapter 1 introduced Glushkov's Second Information Barrier: the claim that neither markets nor hierarchies can coordinate a planetary economy. In this section, we develop the positive counterpart: the claim that there is a third coordination engine — mutual coordination, or stigmergy — that complements markets and hierarchies and is specifically suited to the coordination challenges they cannot handle.

We define each engine precisely.

**Definition 2.3 (Market Coordination).** Market coordination is a decentralized mechanism in which agents signal their preferences and constraints through the prices they are willing to pay or accept. The price vector $p \in \mathbb{R}^L_+$ aggregates dispersed information about relative scarcity and relative desire, and coordinates production and consumption decisions without requiring any central authority.

The informational efficiency of market coordination rests on Hayek's (1945) insight: the knowledge required for efficient allocation is dispersed among millions of agents, each holding fragments of local information about their own circumstances. No central authority can collect and process this knowledge; but the price system can aggregate it through voluntary exchange. The market price of a commodity encodes, in a single number, the supply and demand conditions of all participants — a remarkable compression of distributed information.

The information processing capacity of market coordination can be characterized formally. With $L$ goods and $I$ agents, competitive equilibrium requires $L$ prices (not $I \times L$ bilateral exchange rates, which barter would require). This is the first information barrier that Glushkov identified money as solving: the reduction from $O(I^2)$ to $O(L)$ information requirements.

The limit of market coordination is equally well-defined: it fails when the information relevant to coordination cannot be expressed as a price. Non-rival goods have no efficient price (zero, which cannot sustain production). Ecosystem services have no price (they are not owned and not traded). Future generations' welfare has no price (the unborn cannot participate in markets). And the dynamic interactions of a complex adaptive system — ecological tipping points, systemic financial risk, the feedback between inequality and growth — cannot be represented in any static price vector, no matter how complete.

**Definition 2.4 (Hierarchical Coordination).** Hierarchical coordination is a centralized mechanism in which a principal (a manager, a regulator, a state) issues directives to agents, who act according to instructions rather than autonomous choice. The information processing occurs at the center, and the coordination is achieved through command rather than price.

Hierarchical coordination performs well when the relevant information can be centralized, when the principal's objectives can be specified precisely, when compliance can be monitored, and when the coordination problem is stable enough that centrally designed rules remain valid. Its pathologies — the Hayekian knowledge problem, principal-agent distortions, rent-seeking, bureaucratic rigidity — arise when these conditions fail.

In the context of the Provisioning Question, hierarchies face a specific problem: the information required to maintain complex ecological systems is distributed across millions of local agents (farmers, fishers, foresters, communities) who understand their local conditions in ways that no central authority can replicate. At the same time, the constraints imposed by planetary boundaries are genuinely global and require coordination beyond any local authority. Hierarchies work at neither scale well: too centralized for local knowledge, too fragmented for global coordination.

**Definition 2.5 (Mutual Coordination).** Mutual coordination is a distributed mechanism in which agents coordinate their actions by reading and responding to signals in a shared environment — signals generated by the prior actions of other agents — without requiring either price formation or central direction. The signal medium is the environment itself, modified by collective action.

The biological analogue is stigmergy: the phenomenon by which social insects (ants, termites, bees) coordinate complex collective behavior without central control, through the indirect communication of environmental signals — pheromone trails, structural features of the nest, the positioning of stored food. Each agent responds to local signals; the aggregate behavior that emerges is globally coordinated and often highly efficient.

The economic analogue is commons-based peer production (Benkler, 2006): the coordination of productive activity through shared platforms, open protocols, reputation systems, and community norms, without either market prices or managerial direction. The Linux kernel — built by tens of thousands of contributors, coordinating through shared code repositories, bug trackers, mailing lists, and community conventions — is the most studied example. Wikipedia, OpenStreetMap, scientific preprint servers, and the entire ecosystem of open-source software represent the same mode of production at different scales.

The information processing logic of mutual coordination differs from both markets and hierarchies. Rather than aggregating information into a single price vector or channeling it up a command structure, mutual coordination propagates information laterally through a shared environment, allowing agents to contribute and adapt in parallel. This makes it particularly suited to:

- Non-rival goods, where the zero marginal cost of reproduction makes price signals degenerate.
- Complex adaptive problems, where the required response cannot be specified in advance.
- Long-horizon stewardship, where the relevant information is distributed across local custodians who cannot be effectively monitored by a central authority.
- Cross-boundary commons, where neither market ownership nor state jurisdiction is well-defined.

We formalize mutual coordination as a stigmergic signaling system in Chapter 7. Here, we note only the comparative information processing properties of the three engines:

| Coordination engine | Signal medium | Information aggregation | Suited to |
|:---|:---|:---|:---|
| Market | Price | $O(L)$ prices | Rival, excludable goods; local, short-horizon |
| Hierarchy | Directive | $O(\text{org depth})$ layers | Specifiable tasks; monitorable agents |
| Mutual (stigmergy) | Shared environment | Distributed, parallel | Non-rival goods; adaptive problems; stewardship |

None of the three engines dominates the others across all domains. The argument of this book is not that mutual coordination should replace markets and hierarchies, but that it is a necessary complement to them — one that the existing framework systematically ignores because it has no formal representation in standard economic theory.

---

## 2.5 The Commons: Theory and Evidence

The concept of the commons has suffered, in mainstream economics, from one of the most consequential misreadings in the history of thought. Garrett Hardin's 1968 essay "The Tragedy of the Commons" is widely interpreted as demonstrating that shared resources are inevitably depleted — that the rational self-interest of individual users destroys any commons in the absence of either private property or state control. This reading has justified decades of policy recommendations for privatization of natural resources and regulatory centralization.

Hardin's argument was wrong, not because his logic was invalid, but because his model was mislabeled. What he described was not a commons but an open-access regime — a resource to which everyone has unrestricted access and which no one is responsible for maintaining. The distinction is fundamental.

A genuine commons is a shared resource governed by a defined community of users according to collectively developed rules. Open access is an ungoverned resource that anyone can exploit without restriction. Hardin's analysis correctly shows that open access leads to depletion; it says nothing about commons governance. Elinor Ostrom demonstrated this empirically across hundreds of case studies spanning millennia and continents: from Swiss alpine meadows to Japanese inshore fisheries to irrigation systems in Spain, commons institutions have successfully managed shared resources for generations without either privatization or state control [C:Ch.14].

### 2.5.1 The Commons Game: A Formal Analysis

We formalize the commons problem as an $n$-player game to make the distinction between open access and commons governance precise.

**Setup.** Consider $n$ agents sharing a common-pool resource with stock $R$. Each agent $i$ chooses an extraction level $e_i \geq 0$. The resource regenerates at rate $f(R)$ and is depleted by total extraction $E = \sum_i e_i$. The net resource dynamics are:

$$\dot{R} = f(R) - E, \quad f(R) = rR\left(1 - \frac{R}{K}\right)$$

where $r$ is the intrinsic growth rate and $K$ is the carrying capacity. Agent $i$'s payoff from a single period of extraction is $\pi_i(e_i, E) = p \cdot e_i - c(e_i)$, where $p$ is the price of the extracted resource and $c$ is the cost of extraction, increasing in $e_i$.

**Proposition 2.1 (Open Access Tragedy).** Under open access (no governance), the Nash equilibrium of the single-shot extraction game satisfies:

$$p = c'(e^{NA}_i) \quad \forall i$$

and total extraction $E^{NA} = \sum_i e^{NA}_i > E^{SO}$, where $E^{SO}$ is the socially optimal extraction level satisfying:

$$p = c'(e^{SO}_i) + \frac{p \cdot e^{SO}_i}{R - E^{SO}/r} \quad \forall i$$

The second term in the social optimum condition is the shadow price of the resource stock — the opportunity cost of extraction in terms of foregone future regeneration. Individual agents in the Nash equilibrium do not internalize this cost, leading to overextraction relative to the social optimum.

*Proof sketch.* Each agent maximizes $\pi_i(e_i, E_{-i})$ taking others' extraction $E_{-i}$ as given. The first-order condition $p = c'(e_i)$ ignores the effect of $e_i$ on the future resource stock (the externality). The social optimum maximizes $\sum_i \pi_i$ jointly, which requires internalizing the stock externality. The difference is precisely the shadow price term. $\square$

**This is Hardin's tragedy — but it is the tragedy of open access, not of the commons.** The game just described has no governance: agents act independently, with no communication, no rules, and no enforcement. This is the opposite of a commons, which is precisely defined by the presence of governance.

Now consider the repeated game. Suppose the same $n$ agents interact over an infinite horizon with common discount factor $\delta \in (0,1)$. The socially optimal extraction path — which maintains the resource stock indefinitely — is a feasible outcome of the repeated game. Is it an equilibrium?

**Proposition 2.2 (Cooperation in the Repeated Commons Game).** If agents use grim-trigger strategies — cooperate (extract at the socially optimal level) unless any agent deviates, then switch to Nash forever — and if:

$$\delta \geq \frac{\pi^{D} - \pi^{C}}{\pi^{D} - \pi^{N}}$$

where $\pi^C$ is the per-period payoff under cooperation, $\pi^D$ is the payoff from deviating optimally given others cooperate, and $\pi^N$ is the Nash equilibrium payoff, then the cooperative outcome is a subgame-perfect equilibrium.

*Proof.* At any history in which all agents have cooperated, agent $i$'s payoff from continued cooperation is $\pi^C / (1-\delta)$. The payoff from deviating (exploiting others' cooperation for one period, then receiving the Nash payoff forever) is $\pi^D + \delta \pi^N / (1-\delta)$. Cooperation is preferred if and only if $\pi^C / (1-\delta) \geq \pi^D + \delta \pi^N / (1-\delta)$, which rearranges to the condition stated. $\square$

The condition $\delta \geq (\pi^D - \pi^C)/(\pi^D - \pi^N)$ is the formal statement of a familiar intuition: cooperation is sustainable when agents are sufficiently patient — when the future matters enough relative to the present gain from defection. For a stable community with ongoing relationships, the discount factor is typically high. The one-shot character of Hardin's game is therefore not a description of how commons actually work; it is a description of how they would work if governance were entirely absent.

### 2.5.2 Ostrom's Design Principles

Elinor Ostrom's empirical research identified eight design principles common to long-lived, successful commons institutions. We state these as formal conditions on the governance game:

**Definition 2.6 (Ostrom Design Principles, formalized).** A commons governance system satisfies the Ostrom conditions if:

1. **DP1 (Defined boundaries):** The set of legitimate users $\mathcal{U} \subset \{1, \ldots, n\}$ and the boundaries of the resource system $R$ are clearly defined.

2. **DP2 (Congruence):** The rules governing extraction ($e_i \leq \bar{e}_i$) are adapted to local conditions; the cost-benefit distribution of governance is proportional to local circumstances.

3. **DP3 (Collective choice):** Users can participate in modifying the rules; the governance game has an endogenous rule-setting stage.

4. **DP4 (Monitoring):** Compliance with extraction rules is monitored, either by users themselves or by agents accountable to users.

5. **DP5 (Graduated sanctions):** Violations trigger sanctions proportional to their severity and context; the sanction game is not binary.

6. **DP6 (Conflict resolution):** There are accessible, low-cost mechanisms for resolving disputes among users and between users and rule-setters.

7. **DP7 (Minimal recognition):** External authorities (states, corporations) recognize the legitimacy of the commons governance system and do not override it.

8. **DP8 (Nested enterprises):** For large-scale commons, governance is organized in nested layers, with smaller units embedded in larger ones.

A governance system satisfying DP1–DP8 changes the structure of the extraction game in three ways: it restricts the action space (DP1, DP2), it increases the probability that defection is detected (DP4), and it makes the cost of defection credible and calibrated (DP5, DP6). Together, these changes raise $\delta^*$ — the minimum discount factor required for cooperation — by raising $\pi^N$ (the Nash payoff, which is less attractive when sanctions apply) and lowering $\pi^D$ (the deviation payoff, which is less attractive when detection probability is high).

The formal result is that commons governance, when it satisfies the Ostrom conditions, expands the set of parameters under which cooperation is the equilibrium — sometimes to the point where cooperation is dominant strategy regardless of the discount factor. We develop the full formal analysis in Chapter 14.

---

## 2.6 Mathematical Model: The Stewardship Objective Function

We now introduce the formal optimization problem that underpins the provisioning framework. This will serve as the master objective function against which we measure the performance of cooperative, regenerative economies throughout the book.

**Setup.** A representative social planner chooses a consumption path $\{C_t\}_{t \geq 0}$ to maximize intertemporal welfare, subject to the requirement that both produced and natural capital be maintained.

**Definition 2.7 (Stewardship Objective Function).** The Stewardship Objective Function is:

$$\max_{\{C_t, I_{K,t}, I_{N,t}\}} \int_0^\infty e^{-\rho t} U(C_t, S_t)\, dt$$

subject to:

$$\dot{K}_t = I_{K,t} - \delta_K K_t \geq 0 \quad \text{(produced capital constraint)}$$
$$\dot{N}_t = \mathcal{R}(N_t) - \mathcal{D}(C_t, E_t) + I_{N,t} \geq 0 \quad \text{(natural capital constraint)}$$
$$S_t = \mathcal{S}(K_t, N_t) \quad \text{(provisioning services)}$$
$$C_t + I_{K,t} + I_{N,t} = Y_t = F(K_t, N_t, L) \quad \text{(resource constraint)}$$

where $\rho > 0$ is the discount rate, $K_t$ is produced capital, $N_t$ is natural capital, $S_t$ is the flow of provisioning services, $I_{K,t}$ and $I_{N,t}$ are investments in produced and natural capital respectively, $\delta_K$ is the depreciation rate of produced capital, $\mathcal{R}(N_t)$ is the natural regeneration function, $\mathcal{D}(C_t, E_t)$ is the depletion of natural capital from consumption and extraction $E_t$, and $F$ is the production function.

The key distinction from the standard Ramsey-Cass-Koopmans problem [P:Ch.5] is the natural capital constraint $\dot{N}_t \geq 0$: the Stewardship Objective Function requires that natural capital not decline. This is the Stewardship Constraint, and it is not merely a welfare-improving addition to the standard problem; it is a necessary condition for indefinite provisioning.

**Deriving the First-Order Conditions.** Form the current-value Hamiltonian:

$$\mathcal{H} = U(C_t, S_t) + \mu_K \left[I_{K,t} - \delta_K K_t\right] + \mu_N \left[\mathcal{R}(N_t) - \mathcal{D}(C_t, E_t) + I_{N,t}\right]$$

where $\mu_K$ and $\mu_N$ are the co-state variables (shadow prices) of produced and natural capital respectively. The optimality conditions are:

$$\frac{\partial U}{\partial C_t} = \mu_K + \mu_N \frac{\partial \mathcal{D}}{\partial C_t} \quad \text{(consumption)}$$

$$\dot{\mu}_K = \rho \mu_K - \mu_K \left[F_K - \delta_K\right] - \mu_N \frac{\partial \mathcal{S}}{\partial K_t} \frac{\partial U}{\partial S_t} \quad \text{(produced capital)}$$

$$\dot{\mu}_N = \rho \mu_N - \mu_N \mathcal{R}'(N_t) - \frac{\partial U}{\partial S_t} \frac{\partial \mathcal{S}}{\partial N_t} \quad \text{(natural capital)}$$

The consumption condition is the modified Keynes-Ramsey rule: the marginal utility of consumption equals the shadow price of produced capital $\mu_K$ plus the cost of natural capital depletion induced by that consumption, valued at $\mu_N$. The standard Keynes-Ramsey rule has only the first term; the second term is the stewardship correction.

The natural capital condition shows that the shadow price $\mu_N$ evolves according to the natural regeneration rate $\mathcal{R}'(N_t)$ and the direct utility from provisioning services $\partial U / \partial S_t \cdot \partial \mathcal{S} / \partial N_t$. When $\mu_N > 0$ — when natural capital is genuinely scarce relative to the stewardship constraint — the optimal path involves active investment in natural capital maintenance, not merely passive extraction.

**The Stewardship Constraint as a Side Condition.** When the constraint $\dot{N}_t \geq 0$ binds, the shadow price $\mu_N > 0$ and the consumption path must be modified to accommodate natural capital maintenance. In the standard problem without the constraint, $\mu_N = 0$ and natural capital depletion is optimal whenever it increases present consumption. The Stewardship Objective Function formally rules this out: it treats the maintenance of natural capital not as an optional welfare improvement but as a binding constraint on the planning problem.

This distinction has operational significance. A planner solving the standard Ramsey problem will deplete exhaustible natural capital whenever the return to depletion (in terms of current consumption) exceeds the discount rate. A planner solving the Stewardship Objective will not deplete non-renewable natural capital below its critical threshold $N_{\text{critical}}$, regardless of the discount rate. The two problems have different solutions in the empirically relevant case where $\rho > 0$ and natural capital is non-renewable or slowly renewable.

---

## 2.7 Worked Example: The Welfare Economics of a Non-Rival Good

We demonstrate formally that monopoly pricing of a non-rival digital good with zero marginal cost generates a deadweight loss equal to the entire consumer surplus, and that open-source provision eliminates this loss.

**Setup.** Consider a digital good (a piece of software, a dataset, a research paper) with the following properties:
- Fixed production cost: $F > 0$ (developer time, infrastructure)
- Marginal reproduction cost: $MC = 0$
- Inverse demand: $P(Q) = a - bQ$, where $a, b > 0$ and $Q$ is the number of users

**Under monopoly pricing.** A profit-maximizing monopolist sets marginal revenue equal to marginal cost:

$$MR = a - 2bQ^m = MC = 0 \implies Q^m = \frac{a}{2b}, \quad P^m = \frac{a}{2}$$

Monopoly profit: $\pi^m = P^m Q^m - F = \frac{a^2}{4b} - F$

For the monopoly to be viable, we require $\pi^m \geq 0$, i.e., $F \leq \frac{a^2}{4b}$.

Consumer surplus under monopoly:

$$CS^m = \frac{1}{2}(a - P^m) Q^m = \frac{1}{2} \cdot \frac{a}{2} \cdot \frac{a}{2b} = \frac{a^2}{8b}$$

**Under efficient (zero-price) provision.** The socially efficient allocation has $P^* = MC = 0$, so:

$$Q^* = \frac{a}{b}, \quad P^* = 0$$

Total consumer surplus under efficient provision:

$$CS^* = \frac{1}{2} a \cdot \frac{a}{b} = \frac{a^2}{2b}$$

**The deadweight loss.** The deadweight loss from monopoly pricing of a non-rival good is the difference between efficient consumer surplus and the sum of monopoly consumer surplus and monopoly profit:

$$DWL = CS^* - (CS^m + \pi^m) = \frac{a^2}{2b} - \left(\frac{a^2}{8b} + \frac{a^2}{4b} - F\right) = \frac{a^2}{8b} + F$$

The first term $a^2/8b$ is the standard deadweight loss triangle. The second term $F$ is more striking: under monopoly, the fixed cost $F$ is a private benefit to the producer that comes entirely at the expense of consumers — it is transferred, not created.

But now consider: what is the total potential consumer surplus if the good were provided at zero price? It is $CS^* = a^2/(2b)$. Under monopoly, consumers receive only $CS^m = a^2/(8b)$ — exactly one quarter of the potential value. The monopolist captures $\pi^m + F = a^2/4b$ and a deadweight loss of $a^2/8b$ is destroyed entirely.

**Numerical illustration.** Let $a = 100$, $b = 1$, $F = 200$:

| Allocation | Price | Users ($Q$) | CS | Profit | Total welfare |
|:---|:---|:---|:---|:---|:---|
| Monopoly | 50 | 50 | 1,250 | 2,300 | 3,550 |
| Efficient ($P=0$) | 0 | 100 | 5,000 | −200 | 4,800 |
| Open-source ($P=0$, $F$ via grants/community) | 0 | 100 | 5,000 | 0 | 5,000 |

The efficient outcome requires a mechanism to cover the fixed cost $F$ other than user pricing. Open-source production — funded through voluntary contribution, corporate sponsorship, foundation grants, or public investment — achieves this: users receive the full $CS^*$ and the fixed cost is covered through mechanisms that do not restrict access. The welfare gain from open-source provision relative to monopoly, in this example, is $5{,}000 - 3{,}550 = 1{,}450$ — a 41% improvement in total welfare.

This example is stylized, but its logic is robust. The key insight is that for any non-rival good with $MC = 0$, monopoly pricing destroys value relative to open provision by precisely the amount of consumer surplus foregone — and this amount grows without bound as the number of potential users increases. The larger the potential audience, the greater the case for open provision.

---

## 2.8 Case Study: Wikipedia as an Empirical Test of Non-Rival Good Theory

Wikipedia represents one of the most significant natural experiments in the economics of non-rival goods. In 2001, the encyclopedia market was dominated by two models: the traditional commercial print encyclopedia (Encyclopædia Britannica) and the digitally distributed commercial encyclopedia (Microsoft Encarta, launched 1993). Both were club goods — non-rival in consumption but made excludable through price and access control. Wikipedia launched as a pure public good: freely readable, freely editable, sustained entirely through volunteer labor and donation funding.

Standard theory predicted that Wikipedia should fail. The free-rider problem implies that contributors — who provide labor at private cost to produce a public good — will be undersupplied. Hayek's knowledge problem implies that decentralized, uncoordinated editing will produce inconsistency and error. And the competitive economics of information goods implies that a better-funded commercial competitor (Encarta, backed by Microsoft; Britannica, backed by institutional subscriptions) should be able to produce higher-quality content.

None of these predictions held.

**Scale.** By 2024, English Wikipedia contained approximately 6.8 million articles, compared to Britannica's 120,000 (its full database including online) and Encarta's approximately 60,000 at peak (Encarta was discontinued in 2009). Wikipedia receives approximately 15 billion page views per month globally across all language editions. This is not a niche product: it is among the five most visited websites in the world, sustained without advertising revenue.

**Quality.** A 2005 peer review in *Nature* compared Wikipedia and Britannica on 42 science articles, finding that Wikipedia averaged 3.86 errors per article versus Britannica's 2.92 — a difference not statistically significant at conventional levels and far smaller than most observers expected. Subsequent studies have found that Wikipedia quality on well-trafficked articles has improved substantially and now approaches professional reference quality in many domains (Rector, 2008; Fichman and Kow, 2011).

**Contribution dynamics.** As of 2024, English Wikipedia has approximately 45,000 active editors (making at least five edits per month) out of a total registered user base in the tens of millions. The distribution of contributions is highly skewed — a small fraction of editors contribute a large fraction of total edits, consistent with the power-law distributions typical of peer production systems [C:Ch.4] — but the tail is long and the core is stable.

**The commercial counterfactual.** Microsoft discontinued Encarta in April 2009, explicitly citing Wikipedia competition as a major factor. Encyclopædia Britannica ended print publication in 2012. The commercial model for encyclopedia production was not disrupted by a better-funded competitor; it was displaced by a non-rival good produced through mutual coordination.

**Interpreting the evidence.** Wikipedia does not refute the free-rider logic; it shows that the free-rider problem can be solved by institutional design that creates intrinsic motivation (the satisfaction of contributing to a shared knowledge project), social recognition (reputation among the Wikipedia community), and effective stigmergic coordination (the wiki platform as a shared environment that enables distributed, parallel, asynchronous editing).

The wiki software functions as the stigmergic medium: each edit leaves a trace (in edit history, on talk pages, in article quality ratings) that signals to subsequent editors what work remains, what disputes are unresolved, and what sources need attention. This is mutual coordination in its clearest empirical expression: agents reading environmental signals left by prior actors and contributing in response, without any market price or managerial directive.

Wikipedia thus confirms, across twenty years and hundreds of millions of edits, the core claim of this chapter: non-rival goods are best produced and governed through mutual coordination rather than market pricing or hierarchical direction. The argument is not only theoretical; it has been empirically tested at extraordinary scale.

---

## 2.9 Abundance as a Design Property

We close the chapter by returning to its opening provocation: scarcity is partly a feature of institutions rather than a fact of nature. The converse of this observation is equally important: abundance can be a design property of economic systems that are organized to create and sustain it.

Regenerative ecological systems exhibit this property. A forest managed for long-run productivity — harvested below its regenerative capacity, with soil, hydrology, and biodiversity maintained — produces timber, carbon sequestration, watershed services, and biodiversity indefinitely. The forest does not deplete; it regenerates. The key is that the rate of extraction does not exceed the rate of regeneration: $\mathcal{D}(E) \leq \mathcal{R}(N)$, using the notation of the Stewardship Objective Function [Definition 2.7].

This regenerative logic applies beyond forests. Soil managed with cover crops, reduced tillage, and organic matter addition regenerates fertility over time; soil managed with continuous monoculture and synthetic inputs depletes. Fisheries managed below maximum sustainable yield regenerate; fisheries managed beyond it collapse. The difference between a depleting and a regenerating system is not the physical nature of the resource; it is the governance regime under which extraction occurs.

The same logic extends to social and intellectual capital. A community that invests in education, trust, shared infrastructure, and knowledge commons experiences increasing returns to cooperation: each additional contribution makes the commons more valuable, which attracts more contribution. Open-source software exhibits this dynamic: the Linux kernel, as it grows more sophisticated and more widely deployed, attracts more developer contributions, which make it more sophisticated, in a self-reinforcing cycle that has made it the dominant operating system for servers, mobile devices, and supercomputers.

The economics of abundance is therefore not the economics of cornucopia — it does not claim that resources are infinite or that scarcity has been abolished. It claims, more precisely, that the scarcity or abundance of many economically important goods and systems is determined by institutional design, and that institutions organized around stewardship, mutual coordination, and regeneration can sustain and expand the productive capacity available to human communities, rather than depleting it.

This is the promise of the economics of cooperation — and the remaining chapters of this book are dedicated to making that promise rigorous.

---

## Chapter Summary

This chapter has moved from critique to construction. Having established in Chapter 1 where the standard framework fails, we have begun to build the alternative.

The taxonomy of goods by rivalry and excludability shows that the standard market mechanism is well-suited only to private goods — rival and excludable — and performs poorly across large and growing domains of the economy: common-pool resources, public goods, and non-rival digital goods. For non-rival goods, monopoly pricing imposes a deadweight loss equal to the entire foregone consumer surplus relative to open provision, motivating commons-based alternatives.

The Provisioning Question reframes the central problem of economics from allocation (how do we distribute scarce goods today?) to maintenance (are we sustaining the productive capacity required for indefinite human welfare?). The Stewardship Objective Function formalizes this reframing: it maximizes intertemporal utility while imposing binding constraints on both produced and natural capital maintenance. Its first-order conditions reveal that natural capital has a positive shadow price whenever the stewardship constraint binds — a result with immediate policy implications for natural capital accounting and ecosystem service valuation.

The Three Coordination Engines framework establishes that markets and hierarchies are both essential and insufficient: essential for the domains they handle well, insufficient for the coordination challenges of a planetary economy embedded in a biosphere with long-run dynamics. Mutual coordination — the stigmergic, commons-based mode of production exemplified by Wikipedia, open-source software, and long-lived commons institutions — provides the third engine that complements the other two.

The commons is not the absence of governance; it is a specific form of governance whose conditions — formalized in Ostrom's eight design principles — can be stated precisely and whose performance can be assessed empirically. The repeated game analysis shows that cooperation in the commons is a subgame-perfect equilibrium whenever agents are sufficiently patient and governance is sufficiently robust.

Chapter 3 develops the game-theoretic foundations in full, examining the conditions under which cooperation is not merely possible but optimal.

---

## Exercises

**2.1** Classify each of the following goods into the four quadrants (private good, club good, common-pool resource, public good) and justify your classification: (a) a fish in the open ocean; (b) access to a toll highway; (c) a mathematical theorem; (d) a city park; (e) a brand name; (f) atmospheric CO₂ absorption capacity.

**2.2** Explain why GDP fails to measure welfare adequately for an economy that is simultaneously increasing consumption and depleting natural capital. What additional information would be needed to assess whether such an economy is on a sustainable path? (Reference [P:Ch.3] for the GDP accounting framework.)

**2.3** A digital streaming service produces a film at a fixed cost of $F = \$10$ million and distributes it at zero marginal cost. If the inverse demand for streams is $P(Q) = 20 - 0.001Q$ (with $Q$ in thousands of streams and $P$ in dollars), compute: (a) the monopoly price, quantity, consumer surplus, and profit; (b) the efficient price and quantity; (c) the deadweight loss from monopoly pricing.

**★ 2.4** Prove that Hardin's tragedy of the commons requires a one-shot interaction (or equivalently, an infinitely short time horizon). Specifically: (a) show that the open-access Nash equilibrium involves overextraction relative to the social optimum; (b) show that with infinite repetition and discount factor $\delta \geq (\pi^D - \pi^C)/(\pi^D - \pi^N)$, the socially optimal extraction level is sustainable as a subgame-perfect equilibrium under grim-trigger strategies; (c) discuss what values of $\delta$ are plausible for a stable fishing community and for a mobile population of newcomers, and what this implies for governance design.

**★ 2.5** Consider the Stewardship Objective Function [Definition 2.7] with the specific functional forms $U(C, S) = \ln C + \alpha \ln S$, $S = N^\beta$, $\mathcal{R}(N) = rN(1 - N/K)$, $\mathcal{D}(C, E) = \gamma C$, and $F(K, N, L) = AK^{1-\nu} N^\nu$. (a) Write out the current-value Hamiltonian. (b) Derive the first-order conditions. (c) Show that the steady state requires $\dot{N} = 0$, and find the steady-state natural capital stock $N^*$ as a function of the model parameters. (d) Compare $N^*$ to the natural capital stock that would prevail in the standard Ramsey problem (in which the natural capital constraint is absent).

**★★ 2.6** The Three Coordination Engines framework assigns different information processing functions to markets, hierarchies, and mutual coordination. (a) Formalize the information processing capacity of a price system with $L$ goods and $I$ agents: show that competitive equilibrium requires $L$ prices, not $I(I-1)/2$ bilateral exchange rates. (b) Formalize the information distortion in a hierarchy of depth $d$ and branching factor $k$: show that under a simple error model, the total information loss grows as $O(d \cdot \varepsilon)$ where $\varepsilon$ is the per-layer error rate. (c) Propose a formal model of stigmergic information processing in which each agent observes a local signal $\sigma_i$ from the shared environment and updates their action $a_i$ accordingly. Under what conditions does the aggregate action profile $\{a_i\}$ converge to the socially optimal profile? What are the analogues of the conditions for price-clearing in markets and compliance in hierarchies?

---

*Chapter 3 deepens the game-theoretic foundation. Having introduced the commons as a governance institution and cooperation as an equilibrium of the repeated game, we now develop cooperative game theory in full — the branch of game theory that takes binding agreements as its starting point and asks how the gains from cooperation should be shared.*
