# Chapter 3: The Game-Theoretic Primacy of Cooperation — Foundations in Cooperative Game Theory

> *"The best for the group comes when everyone in the group does what's best for himself and the group."*
> — John Nash, as portrayed in *A Beautiful Mind* (2001) — and promptly corrected by his own mathematics

> *"There is a large cooperative surplus which can be reaped only if people coordinate their behaviour. The central question in social theory is: how is this coordination achieved?"*
> — Michael Taylor, *The Possibility of Cooperation* (1987)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Distinguish cooperative game theory from non-cooperative game theory and identify the economic settings in which each is the appropriate model.
2. Define a cooperative game through its characteristic function, and identify the properties of superadditivity and convexity with their economic interpretations.
3. Define the core of a cooperative game and state the Bondareva–Shapley theorem on core existence.
4. Derive the Shapley value from its four axioms and compute it for small games.
5. Explain the Folk Theorem and its implication that cooperation can emerge as a rational equilibrium in repeated interactions.
6. Apply these concepts to supply chain bargaining, labor markets, and international environmental agreements.

---

## 3.1 Two Branches of Game Theory

Game theory is the formal study of strategic interaction — situations in which the outcome for each agent depends on the choices of others. Since its axiomatic foundations were laid by von Neumann and Morgenstern in 1944 and extended by Nash in the early 1950s, game theory has become one of the most productive frameworks in economics, providing rigorous foundations for the analysis of bargaining, competition, mechanism design, and social norms.

The field divides, at its deepest level, into two branches that differ not in mathematical technique but in their fundamental conception of what agents can do.

**Non-cooperative game theory** treats each agent as acting independently. Agents may communicate, but they cannot make binding commitments: any agreement reached must be self-enforcing — that is, it must be in each agent's independent interest to comply, given what others are doing. The central solution concept is the Nash equilibrium: a profile of strategies from which no agent has a unilateral incentive to deviate. Non-cooperative game theory is the appropriate framework when institutions for enforcing agreements are absent, unreliable, or prohibitively costly.

**Cooperative game theory** treats groups of agents as the fundamental unit of analysis. It assumes that binding agreements are possible — that agents can form coalitions and commit credibly to joint strategies. The central questions shift: not "what will each agent do independently?" but "which coalitions will form, and how will the gains from coalition be distributed?" Cooperative game theory is the appropriate framework when institutions for commitment exist — when courts can enforce contracts, when reputational mechanisms sustain cooperation, when social norms make defection costly, or when the interaction is repeated in ways that make long-run relationships valuable.

The distinction matters for economic analysis in a concrete way. Much of the theoretical case against cooperation — the prisoner's dilemma, the tragedy of the commons, the free-rider problem — is derived from non-cooperative game theory applied to one-shot interactions. It shows, correctly, that independent rational agents without commitment mechanisms will fail to cooperate. But this result does not generalize to the world of repeated interaction, existing institutions, and enforceable agreements that characterizes most real economic relationships. Cooperative game theory is the appropriate tool for that world, and it yields systematically different — and, we will argue, more realistic — predictions.

This chapter develops the foundations of cooperative game theory that will be used throughout the book. The treatment is rigorous but self-contained; no prior knowledge of game theory beyond [P:Ch.2] is assumed.

---

## 3.2 Cooperative Games and the Characteristic Function

The primitive object of cooperative game theory is the characteristic function, which assigns to each possible coalition of agents the maximum total payoff that coalition can guarantee for its members, regardless of what non-members do.

**Definition 3.1 (Cooperative Game).** A cooperative game with transferable utility (TU game) is a pair $(N, v)$ where:
- $N = \{1, 2, \ldots, n\}$ is the set of players.
- $v: 2^N \to \mathbb{R}$ is the characteristic function, assigning to each coalition $S \subseteq N$ the value $v(S)$ that $S$ can guarantee its members, with the normalization $v(\emptyset) = 0$.

The assumption of transferable utility — that the total payoff to a coalition can be freely distributed among its members in any proportion — is a substantive one. It holds when payoffs are monetary and players have quasi-linear utility (utility linear in money), so that interpersonal comparisons of utility are unambiguous. We relax this assumption in Chapter 6, where we consider non-transferable utility games. For the foundational analysis of this chapter, TU games provide sufficient generality.

Two properties of the characteristic function are of central importance.

**Definition 3.2 (Superadditivity).** A game $(N, v)$ is superadditive if for all disjoint coalitions $S, T \subseteq N$ with $S \cap T = \emptyset$:

$$v(S \cup T) \geq v(S) + v(T)$$

Superadditivity says that coalitions can always do at least as well together as they can separately. It is the formal expression of the gains from cooperation: the value of joint action is at least as large as the sum of values from independent action. Superadditivity is the minimal condition for cooperation to be individually rational — it ensures that joining the grand coalition is always weakly better than acting alone.

Superadditivity holds naturally in many economic settings. A supply chain where upstream and downstream firms coordinate production is more valuable than the same firms operating independently. A group of workers who pool their skills for a complex project can produce more than the sum of their individual outputs. A coalition of countries that coordinates carbon abatement achieves environmental targets at lower total cost than countries acting alone.

**Definition 3.3 (Convexity).** A game $(N, v)$ is convex if for all $S \subseteq T \subseteq N$ and all $i \notin T$:

$$v(S \cup \{i\}) - v(S) \leq v(T \cup \{i\}) - v(T)$$

Convexity says that the marginal contribution of a player is non-decreasing in the size of the coalition she joins. This is the cooperative analogue of increasing returns: adding player $i$ to a larger coalition creates more value than adding her to a smaller one. Convexity is a stronger condition than superadditivity; it implies superadditivity but not conversely.

Convexity holds when coalition formation generates positive externalities — when the presence of more members makes each additional member more productive. Knowledge-sharing cooperatives exhibit this property: each new member brings information that increases the value of all existing members' knowledge. Network industries exhibit it too: each new user makes the network more valuable to all existing users.

**Example 3.1 (A Simple Three-Player Game).** Three firms — a raw materials supplier ($i=1$), a manufacturer ($i=2$), and a distributor ($i=3$) — can collaborate to produce and sell a product. Their individual and coalition values (in units of $10^3$ EUR) are:

| Coalition $S$ | $v(S)$ | Interpretation |
|:---|:---|:---|
| $\{1\}$ | 10 | Supplier sells raw materials on spot market |
| $\{2\}$ | 15 | Manufacturer buys inputs at market price, sells output |
| $\{3\}$ | 8 | Distributor sources and resells products |
| $\{1,2\}$ | 40 | Supplier and manufacturer integrate vertically |
| $\{1,3\}$ | 20 | Supplier and distributor form a trading partnership |
| $\{2,3\}$ | 35 | Manufacturer and distributor integrate forward |
| $\{1,2,3\}$ | 70 | Full supply chain cooperation |

We verify superadditivity: $v(\{1,2,3\}) = 70 \geq v(\{1\}) + v(\{2\}) + v(\{3\}) = 33$. The gains from full cooperation are $70 - 33 = 37$ units — more than the sum of all individual values.

We check convexity by computing marginal contributions. Player 1's marginal contribution to $\{2\}$ is $v(\{1,2\}) - v(\{2\}) = 40 - 15 = 25$. Player 1's marginal contribution to $\{2,3\}$ is $v(\{1,2,3\}) - v(\{2,3\}) = 70 - 35 = 35 > 25$. The marginal contribution increases with coalition size, consistent with convexity (a full check for all players and all nested pairs of coalitions would confirm this).

---

## 3.3 The Core: Stability and Its Conditions

The most fundamental solution concept in cooperative game theory is the core: the set of allocations of the grand coalition's value such that no subset of players would prefer to break away and form their own coalition.

**Definition 3.4 (The Core).** The core of a game $(N, v)$ is the set of payoff vectors $x = (x_1, \ldots, x_n) \in \mathbb{R}^n$ such that:

1. **Efficiency:** $\sum_{i \in N} x_i = v(N)$ — the grand coalition's value is fully distributed.
2. **Coalitional rationality:** For all $S \subseteq N$, $\sum_{i \in S} x_i \geq v(S)$ — no coalition can improve on its members' payoffs by defecting.

The core is the cooperative analogue of Nash equilibrium: it is the set of allocations from which no group of players has a collective incentive to deviate. A payoff vector is in the core if and only if it is stable against all possible defections.

The core may be empty — there may exist no allocation of $v(N)$ that satisfies all coalitional rationality constraints simultaneously. When the core is empty, the grand coalition is unstable: any proposed allocation will be blocked by some coalition, and cooperation is difficult to sustain without side payments or institutional arrangements outside the model.

**Theorem 3.1 (Bondareva–Shapley Theorem, 1963).** The core of a TU game $(N, v)$ is non-empty if and only if the game is balanced.

A game is balanced if for every balanced collection of coalitions $\mathcal{B} = \{S_1, \ldots, S_k\}$ — a collection such that for each player $i$, $\sum_{S \ni i} \lambda_S = 1$ for some weights $\lambda_S > 0$ — the following holds:

$$\sum_{S \in \mathcal{B}} \lambda_S v(S) \leq v(N)$$

*Proof sketch.* The core existence problem can be written as a linear program: find $x$ satisfying $\sum_{i \in N} x_i = v(N)$ and $\sum_{i \in S} x_i \geq v(S)$ for all $S \subseteq N$. By LP duality, the primal is feasible (core is non-empty) if and only if the dual has a bounded optimum — which is precisely the balancedness condition. $\square$

The Bondareva–Shapley theorem has a direct economic interpretation: the core is non-empty if and only if no weighted combination of sub-coalition values exceeds the grand coalition value. In other words, cooperation is stable whenever the grand coalition does at least as well as any convex combination of sub-coalitions — a condition that is closely related to, but weaker than, full convexity.

**Proposition 3.1 (Convex Games Have Non-Empty Cores).** Every convex game has a non-empty core.

*Proof.* It suffices to show that any convex game is balanced. For a convex game and any balanced collection $\mathcal{B}$ with weights $\{\lambda_S\}$, the balancedness condition $\sum_S \lambda_S v(S) \leq v(N)$ follows from the supermodularity property of convex games by induction on the number of players. $\square$

**Example 3.1 (continued).** Returning to the three-firm supply chain, we compute the core constraints:

For all individual players: $x_1 \geq 10$, $x_2 \geq 15$, $x_3 \geq 8$.

For all two-player coalitions: $x_1 + x_2 \geq 40$, $x_1 + x_3 \geq 20$, $x_2 + x_3 \geq 35$.

Efficiency: $x_1 + x_2 + x_3 = 70$.

Adding the three two-player constraints: $2(x_1 + x_2 + x_3) \geq 40 + 20 + 35 = 95$, which gives $x_1 + x_2 + x_3 \geq 47.5$. Since $v(N) = 70 \geq 47.5$, no contradiction arises, and the core is non-empty.

One allocation in the core: $x = (15, 37, 18)$. Verification: $15 \geq 10$ ✓; $37 \geq 15$ ✓; $18 \geq 8$ ✓; $15 + 37 = 52 \geq 40$ ✓; $15 + 18 = 33 \geq 20$ ✓; $37 + 18 = 55 \geq 35$ ✓; $15 + 37 + 18 = 70$ ✓. This allocation is stable: no subset of firms can do better by forming their own coalition.

---

## 3.4 The Shapley Value: Fairness in Cooperation

The core tells us which allocations are stable, but it does not tell us which stable allocation will be chosen. When the core is large, many allocations satisfy the stability constraints, and additional criteria are needed to select among them. When the core is empty, we need a solution concept that provides a unique recommendation even in the absence of stable allocations.

The Shapley value, introduced by Lloyd Shapley in 1953, is the uniquely determined allocation that satisfies four axioms of fairness. It is one of the most elegant results in cooperative game theory, and it will appear repeatedly throughout this book as a tool for fair allocation in cooperatives, supply chains, data economies, and public goods provision.

### 3.4.1 The Four Axioms

**Axiom 3.1 (Efficiency).** The Shapley value $\phi(v) = (\phi_1(v), \ldots, \phi_n(v))$ distributes the entire value of the grand coalition:

$$\sum_{i \in N} \phi_i(v) = v(N)$$

**Axiom 3.2 (Symmetry).** If players $i$ and $j$ are symmetric — meaning $v(S \cup \{i\}) = v(S \cup \{j\})$ for all $S \subseteq N \setminus \{i, j\}$ — then $\phi_i(v) = \phi_j(v)$.

**Axiom 3.3 (Null Player).** If player $i$ contributes nothing to any coalition — $v(S \cup \{i\}) = v(S)$ for all $S \subseteq N \setminus \{i\}$ — then $\phi_i(v) = 0$.

**Axiom 3.4 (Additivity).** For any two games $(N, v)$ and $(N, w)$: $\phi_i(v + w) = \phi_i(v) + \phi_i(w)$, where $(v + w)(S) = v(S) + w(S)$.

**Theorem 3.2 (Shapley, 1953).** There exists a unique value function $\phi: v \mapsto \phi(v)$ satisfying Axioms 3.1–3.4, given by:

$$\phi_i(v) = \sum_{S \subseteq N \setminus \{i\}} \frac{|S|!\,(n - |S| - 1)!}{n!} \left[v(S \cup \{i\}) - v(S)\right]$$

*Proof of uniqueness.* Every game $(N, v)$ can be written as a linear combination of unanimity games $u_T$, defined by $u_T(S) = 1$ if $T \subseteq S$ and $0$ otherwise, which form a basis for the space of all TU games. By Additivity (Axiom 3.4), it suffices to determine $\phi$ on unanimity games. By Symmetry (Axiom 3.2), all players in $T$ receive equal shares of $u_T$. By the Null Player axiom (Axiom 3.3), players not in $T$ receive zero. By Efficiency (Axiom 3.1), the equal shares must sum to 1. Thus $\phi_i(u_T) = 1/|T|$ for $i \in T$ and $0$ otherwise. Since $\phi$ is determined on the basis, it is unique. $\square$

### 3.4.2 The Marginal Contribution Interpretation

The Shapley value formula has an illuminating interpretation: $\phi_i(v)$ is player $i$'s expected marginal contribution to the coalition they join, when the order of coalition formation is chosen uniformly at random.

To see this: imagine that all $n!$ orderings of the players are equally likely. In a given ordering, player $i$ joins a coalition $S$ (the set of players who precede $i$ in the ordering) and contributes $v(S \cup \{i\}) - v(S)$. The Shapley value is the expected value of this marginal contribution across all orderings. The coefficient $|S|!(n-|S|-1)!/n!$ is precisely the probability that ordering places exactly the members of $S$ before player $i$ and the remaining $n - |S| - 1$ players after.

This interpretation makes the Shapley value's fairness properties transparent. It allocates to each player their average marginal contribution — not their contribution to the grand coalition alone, nor to the smallest coalition, but their expected contribution when cooperation builds up from scratch in a random order. It rewards players for what they bring, not for their bargaining power or their coalition position.

**Example 3.1 (continued).** Computing the Shapley value for the three-firm supply chain:

For player 1 (supplier), the marginal contributions across all $3! = 6$ orderings are:

| Ordering | Coalition $S$ when 1 joins | $v(S \cup \{1\}) - v(S)$ |
|:---|:---|:---|
| $(1,2,3)$ | $\emptyset$ | $v(\{1\}) - v(\emptyset) = 10$ |
| $(1,3,2)$ | $\emptyset$ | $10$ |
| $(2,1,3)$ | $\{2\}$ | $v(\{1,2\}) - v(\{2\}) = 25$ |
| $(2,3,1)$ | $\{2,3\}$ | $v(\{1,2,3\}) - v(\{2,3\}) = 35$ |
| $(3,1,2)$ | $\{3\}$ | $v(\{1,3\}) - v(\{3\}) = 12$ |
| $(3,2,1)$ | $\{2,3\}$ | $35$ |

$$\phi_1(v) = \frac{1}{6}(10 + 10 + 25 + 35 + 12 + 35) = \frac{127}{6} \approx 21.2$$

By analogous computation: $\phi_2(v) \approx 30.8$ and $\phi_3(v) \approx 18.0$.

Check: $21.2 + 30.8 + 18.0 = 70.0 = v(N)$ ✓.

Compare this to the core allocation $(15, 37, 18)$ computed earlier. The Shapley value allocates less to the manufacturer (30.8 vs. 37) and more to the supplier (21.2 vs. 15), reflecting the supplier's high marginal contribution when joining the manufacturer-distributor coalition. Both are in the core, but the Shapley value selects a specific point based on average marginal contribution rather than bargaining leverage.

---

## 3.5 The Folk Theorem and the Emergence of Cooperation

The cooperative game framework assumes that binding agreements are possible. But where do the institutions that make agreements binding come from? This question leads to a second major result in game theory: the Folk Theorem, which demonstrates that cooperative outcomes can emerge as equilibria of repeated non-cooperative games, even without external enforcement mechanisms. In this sense, cooperation is not merely assumed by cooperative game theory — it can be derived from first principles through repeated interaction.

### 3.5.1 Setup: The Repeated Prisoner's Dilemma

Consider the Prisoner's Dilemma, the canonical model of the tension between individual and collective rationality:

|  | Cooperate ($C$) | Defect ($D$) |
|:---|:---|:---|
| **Cooperate ($C$)** | $(R, R)$ | $(S, T)$ |
| **Defect ($D$)** | $(T, S)$ | $(P, P)$ |

with $T > R > P > S$ and $2R > T + S$ (so mutual cooperation is socially efficient). In the one-shot game, defection is the dominant strategy: regardless of what the other player does, defection yields a higher payoff ($T > R$ and $P > S$). The unique Nash equilibrium is mutual defection, yielding $(P, P)$, despite the fact that $(R, R)$ Pareto-dominates it.

Now suppose the game is repeated infinitely, with each player discounting future payoffs by factor $\delta \in (0,1)$ per period. The payoff from any strategy profile is the present discounted value of the infinite sequence of stage-game payoffs.

### 3.5.2 The Grim Trigger Strategy and Folk Theorem

**Definition 3.5 (Grim Trigger Strategy).** The grim trigger strategy for player $i$ is:
- Play $C$ in period 1.
- Play $C$ in period $t > 1$ if all players have played $C$ in all periods $1, \ldots, t-1$.
- Play $D$ forever if any player has ever played $D$.

**Theorem 3.3 (Folk Theorem, Discount Factor Version).** In the infinitely repeated Prisoner's Dilemma, mutual cooperation sustained by grim trigger strategies is a subgame-perfect Nash equilibrium if and only if:

$$\delta \geq \delta^* \equiv \frac{T - R}{T - P}$$

*Proof.* If both players use grim trigger, player $i$'s payoff from continued cooperation is:
$$V^C = R + \delta R + \delta^2 R + \cdots = \frac{R}{1 - \delta}$$

If player $i$ deviates (plays $D$ while the other plays $C$), they receive $T$ in the deviation period, then $P$ forever (since the grim trigger is activated):
$$V^D = T + \delta P + \delta^2 P + \cdots = T + \frac{\delta P}{1 - \delta}$$

Cooperation is preferred if $V^C \geq V^D$:
$$\frac{R}{1 - \delta} \geq T + \frac{\delta P}{1 - \delta}$$
$$R \geq T(1 - \delta) + \delta P$$
$$R - P \geq (T - P)(1 - \delta)$$
$$\delta \geq 1 - \frac{R - P}{T - P} = \frac{T - R}{T - P} \equiv \delta^*$$

Subgame perfection requires that the threat to revert to permanent defection is credible — which it is, since permanent defection is itself a Nash equilibrium of the stage game. $\square$

**Interpretation.** The critical discount factor $\delta^* = (T-R)/(T-P)$ is the ratio of the short-run temptation gain $(T-R)$ to the maximum possible punishment loss relative to cooperation $(T-P)$. Cooperation is sustained when agents are sufficiently patient — when the long-run loss from triggering the punishment outweighs the short-run gain from defection.

**The General Folk Theorem.** The result extends far beyond the Prisoner's Dilemma. In the general version (Fudenberg and Maskin, 1986), any feasible payoff vector that Pareto-dominates the Nash equilibrium payoffs of the stage game can be sustained as a subgame-perfect equilibrium of the infinitely repeated game, provided players are sufficiently patient. The set of sustainable payoffs expands as $\delta \to 1$, converging to the full set of individually rational, feasible payoffs.

This is a profound result. It means that virtually any cooperative outcome — any distribution of payoffs that rational players would individually prefer to the non-cooperative baseline — can, in principle, be sustained as an equilibrium of repeated interaction. The institution of cooperation is not a deus ex machina imposed from outside the model; it emerges from the rational calculation of forward-looking agents who value their long-run relationships.

### 3.5.3 Beyond Grim Trigger: Robustness and Proportionality

The grim trigger strategy has an important limitation: it is unforgiving. Any mistake — a misread signal, a noisy observation, an inadvertent defection — triggers permanent punishment, which is both welfare-destroying and potentially disproportionate. In practice, the most robust cooperative strategies are more lenient.

Robert Axelrod's famous computer tournaments (1984) identified Tit-for-Tat — cooperate on the first move, then mirror the other player's previous action — as the most successful strategy across a wide range of opponents. Tit-for-Tat is cooperative (starts friendly), retaliatory (punishes immediately), forgiving (returns to cooperation after one punishment cycle), and transparent (simple enough that opponents can quickly learn it). These properties make it robust to invasion by defectors while avoiding the welfare destruction of permanent punishment.

Formally, the set of discount factors under which Tit-for-Tat sustains cooperation in the Prisoner's Dilemma is:

$$\delta \geq \frac{T - R}{T - S}$$

This threshold is weakly lower than the grim trigger threshold $\delta^* = (T-R)/(T-P)$ when $P \geq S$ (which holds by assumption), meaning Tit-for-Tat requires weakly less patience than grim trigger to sustain cooperation. The more forgiving strategy is also the more efficient one.

We develop the full evolutionary analysis of cooperative strategies in Chapter 7, including the stigmergic mechanisms through which cooperative norms propagate across populations. Here, the Folk Theorem provides the game-theoretic foundation: cooperation is not a departure from rational self-interest but a rational strategy for agents who interact repeatedly and care about the future.

---

## 3.6 Applications: Supply Chains, Labor Markets, and Environmental Agreements

The abstract framework of cooperative game theory connects to concrete economic problems through a set of recurring applications. We survey three, each of which reappears in later chapters.

### 3.6.1 Bargaining in Supply Chains

A supply chain is a cooperative game in which upstream and downstream firms must agree on how to share the surplus from vertical coordination. The characteristic function $v(S)$ captures the value that any subset of supply chain participants can generate through their joint activities. The Shapley value provides a principled baseline for how this surplus should be distributed — one that reflects each firm's average marginal contribution rather than its bargaining power.

In practice, supply chain contracts rarely implement the Shapley value exactly; but the theory provides a benchmark against which actual contracts can be assessed. When actual contracts deviate significantly from the Shapley value, this typically reflects either unequal bargaining power (large buyers or suppliers extracting rents through market power) or incomplete information (firms not knowing others' outside options). Both deviations create inefficiencies that cooperative institutions — supply chain cooperatives, long-term relational contracts, open-book costing — can partially correct.

The supply chain game also illustrates why the core is often a more natural solution concept than the Shapley value in bilateral negotiations: when there are only two parties (buyer and seller), the core reduces to the set of individually rational, jointly efficient agreements — the classical bargaining set — and the Shapley value selects the Nash bargaining solution (the midpoint of the surplus, under equal bargaining weights).

### 3.6.2 Coalition Formation in Labor Markets

Worker cooperatives and labor unions are cooperative games in which workers pool their collective bargaining power to negotiate with employers. The characteristic function captures the value that any subset of workers can produce — either on their own (in a cooperative) or in negotiation with an employer.

The Shapley value has a natural interpretation in this context: it is the fair wage for each worker, reflecting their marginal contribution to the collective enterprise. In a well-governed cooperative where the Shapley value is approximately implemented, workers are paid according to what they contribute — a principle that aligns private and collective incentives far more closely than the standard wage contract, which pays workers their marginal product in competitive equilibrium but does not share the surplus from cooperation.

The stability of labor coalitions depends on core non-emptiness. A union is stable (in the cooperative game sense) if no subset of workers would be better off negotiating separately with the employer. This condition is more likely to hold when workers are complementary — when the value of any team exceeds the sum of individual values — and less likely when workers are substitutes — when the employer can replace any coalition with an equivalent one.

### 3.6.3 International Environmental Agreements

International environmental agreements (IEAs) are cooperative games among sovereign nations, with the characteristic function capturing each coalition's ability to reduce emissions and the economic cost of doing so. The central challenge is that no supranational enforcement mechanism exists: any agreement must be self-enforcing, placing IEAs at the boundary between cooperative and non-cooperative game theory.

The formal analysis of IEAs as cooperative games yields three important results. First, the core of a global climate agreement is typically empty: any allocation of abatement costs that is acceptable to the grand coalition can be blocked by a coalition of nations that prefer to free-ride on others' efforts. This is not a failure of diplomacy; it is a mathematical consequence of the structure of the problem. Second, the maximum stable coalition — the largest group of nations that can sustain a self-enforcing agreement — is typically small (Barrett, 1994), because the incentive to free-ride increases with the size of the coalition (since larger coalitions reduce emissions more, making free-riding more attractive). Third, side payments — transfers of value between nations, which correspond to the difference between the Shapley value and each nation's unilateral payoff — can dramatically expand the set of stable agreements, by compensating nations that bear disproportionate abatement costs.

We return to this analysis in Chapter 32, where we examine the relationship between cooperative institutions and the global commons, and in the worked example on the Kyoto Protocol in the exercises.

---

## 3.7 Worked Example: A Three-Firm Supply Chain

We now work through the three-firm supply chain game [Example 3.1] in full, computing the core, the Shapley value, and the Nash bargaining solution, and comparing their implications for supply chain governance.

**The game.** Recall the characteristic function:

$$v(\{1\}) = 10, \quad v(\{2\}) = 15, \quad v(\{3\}) = 8$$
$$v(\{1,2\}) = 40, \quad v(\{1,3\}) = 20, \quad v(\{2,3\}) = 35$$
$$v(\{1,2,3\}) = 70$$

**Step 1: Core.** The core requires:

$$x_1 \geq 10, \quad x_2 \geq 15, \quad x_3 \geq 8$$
$$x_1 + x_2 \geq 40, \quad x_1 + x_3 \geq 20, \quad x_2 + x_3 \geq 35$$
$$x_1 + x_2 + x_3 = 70$$

The binding constraints are $x_1 + x_2 \geq 40$ and $x_2 + x_3 \geq 35$. Combined with efficiency ($x_1 + x_2 + x_3 = 70$), these imply $x_3 \leq 30$ and $x_1 \leq 35$.

The core is non-empty and takes the form:

$$\mathcal{C}(v) = \{(x_1, x_2, x_3) : x_1 \in [10, 35],\; x_2 \in [15, 52],\; x_3 \in [8, 30],\; x_1+x_2+x_3=70, \text{ all constraints met}\}$$

This is a two-dimensional polytope (the intersection of a hyperplane with several half-spaces in $\mathbb{R}^3$). The core is large: many allocations are stable. Additional selection criteria are needed.

**Step 2: Shapley value.** We already computed $\phi(v) \approx (21.2, 30.8, 18.0)$.

Full computation via the formula:

$$\phi_1(v) = \frac{0!\cdot 2!}{3!}[v(\{1\})-0] + \frac{1!\cdot 1!}{3!}[v(\{1,2\})-v(\{2\})] + \frac{1!\cdot 1!}{3!}[v(\{1,3\})-v(\{3\})] + \frac{2!\cdot 0!}{3!}[v(\{1,2,3\})-v(\{2,3\})]$$

$$= \frac{1}{3}(10) + \frac{1}{6}(25) + \frac{1}{6}(12) + \frac{1}{3}(35) = \frac{10}{3} + \frac{25}{6} + \frac{12}{6} + \frac{35}{3} = \frac{20+25+12+70}{6} = \frac{127}{6} \approx 21.17$$

$$\phi_2(v) = \frac{1}{3}(15) + \frac{1}{6}(25) + \frac{1}{6}(27) + \frac{1}{3}(35) = \frac{15}{3} + \frac{25}{6} + \frac{27}{6} + \frac{35}{3} = \frac{30+25+27+70}{6} = \frac{152}{6} \approx 25.33$$

Wait — let us recompute $\phi_2$ carefully. The marginal contributions of player 2 are:

- To $\emptyset$: $v(\{2\}) - v(\emptyset) = 15$
- To $\{1\}$: $v(\{1,2\}) - v(\{1\}) = 30$
- To $\{3\}$: $v(\{2,3\}) - v(\{3\}) = 27$
- To $\{1,3\}$: $v(\{1,2,3\}) - v(\{1,3\}) = 50$

$$\phi_2(v) = \frac{1}{3}(15) + \frac{1}{6}(30) + \frac{1}{6}(27) + \frac{1}{3}(50) = 5 + 5 + 4.5 + \frac{50}{3} = 14.5 + 16.67 = 31.17$$

And $\phi_3(v) = 70 - 21.17 - 31.17 = 17.67$.

Rounding: $\phi(v) \approx (21.2, 31.2, 17.6)$, summing to $70$. ✓

The Shapley value is in the core: we verify $21.2 + 31.2 = 52.4 \geq 40$ ✓, $21.2 + 17.6 = 38.8 \geq 20$ ✓, $31.2 + 17.6 = 48.8 \geq 35$ ✓.

**Step 3: Nash bargaining solution.** The Nash bargaining solution maximizes the product of gains over the disagreement point $d = (v(\{1\}), v(\{2\}), v(\{3\})) = (10, 15, 8)$:

$$\max_{x \in \mathcal{C}(v)} (x_1 - 10)(x_2 - 15)(x_3 - 8) \quad \text{s.t.} \quad x_1 + x_2 + x_3 = 70$$

The feasible gains are $(x_1 - 10, x_2 - 15, x_3 - 8)$ with $\sum (x_i - d_i) = 70 - 33 = 37$. For a product-maximizing allocation with a fixed sum, the solution is equal gains:

$$x_1 - 10 = x_2 - 15 = x_3 - 8 = \frac{37}{3} \approx 12.33$$

Giving $x^{NBS} = (22.33, 27.33, 20.33)$.

**Summary and interpretation.**

| Solution | $x_1$ (Supplier) | $x_2$ (Manufacturer) | $x_3$ (Distributor) | Notes |
|:---|:---|:---|:---|:---|
| Shapley value | 21.2 | 31.2 | 17.6 | Rewards average marginal contribution |
| Nash bargaining | 22.3 | 27.3 | 20.3 | Rewards outside option strength equally |
| One core allocation | 15.0 | 37.0 | 18.0 | Manufacturer captures most of bilateral surplus |
| Another core allocation | 32.0 | 20.0 | 18.0 | Supplier captures most of bilateral surplus |

The range within the core shows the scope for bargaining power to determine outcomes within the set of stable agreements. The Shapley value and Nash bargaining solution select specific points within (or near) the core, using different fairness criteria. The choice between them is not merely technical; it reflects a normative commitment to how cooperation should be rewarded.

In a cooperative enterprise — a supply chain cooperative, a worker-owned firm, or a commons-based production system — the Shapley value provides the most defensible allocation: it rewards each party for what they actually contribute, averaged across all possible orderings of cooperation. We will return to this principle throughout the book, from the design of cooperative enterprises in Chapter 34 to the fair allocation of data rents in Chapter 39.

---

## 3.8 Case Study: OPEC as a Cooperative Game

The Organization of the Petroleum Exporting Countries (OPEC), founded in 1960, is one of the most studied and least stable cooperative coalitions in economic history. Its goal is to coordinate oil production among member states to maintain prices above the competitive level — a classic cartel problem. Its history of partial success and recurrent failure provides a rich test of cooperative game theory.

### 3.8.1 The OPEC Game

Model OPEC as a cooperative game among $n$ oil-producing nations. Each nation $i$ has production capacity $q_i^{\max}$ and extraction cost $c_i$ per barrel. The market inverse demand is $P(Q) = a - bQ$. The characteristic function $v(S)$ for coalition $S$ is the maximum profit that $S$ can earn by jointly choosing production levels, taking non-members' production as given at their Nash equilibrium levels.

The grand coalition OPEC is superadditive (all members together can restrict output more effectively than any subset), but the core exhibits a fundamental problem: the free-rider incentive grows with coalition size.

### 3.8.2 Why the Core Is Often Empty

Consider a simplified version with three producer types: a large producer (Saudi Arabia, with low extraction costs and high capacity), a medium producer (e.g., Iraq or Iran), and a high-cost small producer. The coalition's value depends critically on whether low-cost producers are willing to bear a disproportionate share of the output reduction — since high-cost producers benefit more from the price increase relative to their counterfactual profits.

The core of this game is typically empty for the following reason: the grand coalition requires large, low-cost producers to restrict output heavily (since they have the most capacity to cut), while high-cost producers benefit substantially from the price increase. But the low-cost producers can form their own sub-coalition — bilateral coordination between major producers — and do better than their Shapley value in the grand coalition. This blocking coalition destabilizes the grand coalition, shrinking the core.

Formally, for the OPEC game with $n > 3$ sufficiently heterogeneous producers, one can show:

**Proposition 3.2 (OPEC Core Instability).** When extraction costs are sufficiently heterogeneous across producers and output quotas must be assigned in proportion to capacity, the core of the OPEC cooperative game is empty.

The intuition: any quota allocation that keeps all members in the coalition is blocked by the coalition of low-cost producers, who can restrict output to a lesser degree while still achieving a higher price than under competition.

### 3.8.3 What Happens When the Core Is Empty

When the core is empty, stable cooperative outcomes are impossible without side payments — transfers between coalition members that compensate those who bear disproportionate costs. OPEC's historical record reflects this structure precisely.

The 1973 oil embargo succeeded because the political solidarity of Arab member states provided a side payment (political credibility, shared geopolitical goals) that supplemented the economic incentives. The price collapses of 1986 and 2014 occurred when side payments became insufficient: Saudi Arabia, facing a blocking coalition of smaller producers who continued to overproduce, responded by flooding the market — effectively defecting from the grand coalition and accepting the Nash equilibrium outcome.

The 2016 OPEC+ agreement, which expanded the coalition to include Russia and other non-OPEC producers, represents an attempt to change the structure of the game by expanding the coalition to include the most important free-riders. Whether this changes the core existence conditions depends on whether the enlarged coalition is sufficiently homogeneous in cost structure — an empirical question that the cooperative game framework makes precise.

### 3.8.4 Lessons for Cooperative Design

The OPEC case illustrates three principles that will recur throughout this book:

1. **Superadditivity does not guarantee stability.** Even when cooperation creates value, that value may not be distributable in a way that keeps all parties in the coalition. Core non-emptiness requires more than superadditivity; it requires the balancedness condition of Theorem 3.1.

2. **Side payments expand stability.** When the core is empty, side payments — transfers between coalition members — can restore stability by compensating parties who bear disproportionate costs. The design of fair side payment mechanisms is a central problem in cooperative economics, to which we return in the context of international climate agreements [Chapter 32] and complementary currency design [Chapter 25].

3. **Coalition structure matters as much as the grand coalition.** The relevant unit of analysis is not just the grand coalition but the entire lattice of possible sub-coalitions. Understanding which sub-coalitions are likely to form — and which are likely to block the grand coalition — is essential for the design of robust cooperative institutions.

---

## Chapter Summary

This chapter has developed the theoretical core of cooperative game theory and established the game-theoretic case for the primacy of cooperation.

The characteristic function formalizes what each coalition of agents can achieve by acting jointly. Superadditivity — the condition that coalitions can always do at least as well together as separately — is the formal expression of the gains from cooperation, and it is widely satisfied across economic settings. Convexity — the condition that marginal contributions are non-decreasing in coalition size — is the stronger condition of increasing returns to cooperation, which holds in network industries, knowledge economies, and ecological systems.

The core is the set of allocations from which no coalition would defect. The Bondareva–Shapley theorem provides a complete characterization of core non-emptiness: the core exists if and only if the game is balanced, a condition that is implied by (but weaker than) convexity. When the core is non-empty, cooperation is stable without side payments; when it is empty, stability requires institutional mechanisms beyond the payoff structure of the game.

The Shapley value is the unique allocation satisfying four axioms of fairness: efficiency, symmetry, the null player property, and additivity. It allocates to each player their average marginal contribution across all orderings of coalition formation. It provides a principled benchmark for the distribution of cooperative surplus — one that appears repeatedly in cooperative enterprise design, commons governance, and the economics of the digital economy.

The Folk Theorem establishes that cooperation can emerge as a rational equilibrium of repeated interaction, even without external enforcement. The critical discount factor below which cooperation breaks down depends on the stage-game payoffs; above this threshold, any efficient cooperative outcome is sustainable. The result anchors the game-theoretic case for cooperation in rational choice theory rather than social preferences.

Chapter 4 turns from the structure of games to the structure of the networks through which economic agents interact — the architecture that determines who can form coalitions with whom, and how information and value flow through the economy.

---

## Exercises

**3.1** Define the characteristic function of a three-player public goods game in which each player can contribute 0 or 1 unit of a public good at cost $c = 0.8$, and each player receives a benefit of $0.5g$ where $g$ is total contribution. Write out $v(S)$ for all $S \subseteq \{1,2,3\}$. Is the core non-empty? If so, find one allocation in the core.

**3.2** Explain in plain language why the Nash bargaining solution differs from the Shapley value as an allocation of cooperative surplus. In which economic settings would you prefer to use each as a normative benchmark, and why?

**3.3** In the Folk Theorem, the critical discount factor for sustaining cooperation under grim trigger is $\delta^* = (T-R)/(T-P)$. For the standard numerical Prisoner's Dilemma with $T=5, R=3, P=1, S=0$: (a) compute $\delta^*$; (b) show that Tit-for-Tat has a weakly lower threshold; (c) explain why a forgiving strategy can require less patience than an unforgiving one.

**★ 3.4** Prove that every convex game has a non-empty core. Your proof should proceed by showing that any convex game is balanced. (*Hint:* use the marginal contribution interpretation of the Shapley value and the fact that convexity implies the Shapley value lies in the core of a convex game — a result due to Shapley, 1971.)

**★ 3.5** Consider a four-firm supply chain with characteristic function values: $v(\{i\}) = 10$ for all $i$; $v(S) = 30|S|$ for $|S| = 2$; $v(S) = 80|S|$ for $|S| = 3$; $v(N) = 400$.
(a) Show that this game is convex.
(b) Compute the Shapley value.
(c) Show that the Shapley value lies in the core.
(d) Interpret the Shapley value allocation in terms of a cooperative enterprise governance principle.

**★★ 3.6** Formalize the Kyoto Protocol as a cooperative game. Let $N$ be the set of signatory nations, each with abatement cost function $c_i(a_i)$ (convex, increasing in abatement level $a_i$) and benefit function $b(A)$ where $A = \sum_i a_i$ is total global abatement (concave, increasing). The characteristic function is $v(S) = \max_{\{a_i\}_{i \in S}} \left[\sum_{i \in S} b_i(A^S + A^{N\setminus S}^*) - c_i(a_i)\right]$ where $A^{N\setminus S}^*$ is the Nash equilibrium abatement of non-members.
(a) Show that the game is superadditive.
(b) Under what conditions on $c_i(\cdot)$ and $b(\cdot)$ is the game convex?
(c) Show that when cost functions are heterogeneous, the core of the game without side payments may be empty.
(d) Compute the minimum side payment from high-benefit to high-cost nations required to make full participation individually rational. How does this relate to the Green Climate Fund?

---

*Chapter 4 moves from the structure of cooperation to the structure of the networks through which agents interact. The architecture of economic networks — who is connected to whom, how many connections each agent maintains, how information flows through the graph — turns out to be a first-order determinant of economic outcomes, not a secondary feature of institutional context.*
