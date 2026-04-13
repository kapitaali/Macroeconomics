# Part II: Core Theory I — Cooperation and Peer-to-Peer Dynamics

The five chapters of Part I assembled a foundation: a precise diagnosis of where the standard framework fails, a reorientation toward provisioning and stewardship, a first encounter with cooperative game theory and network structure, and an introduction to the complex adaptive dynamics that make equilibrium analysis insufficient on its own. What Part I did not do — deliberately — is prove that the alternative is better. It identified the cage; it did not yet demonstrate what lies outside it.

Part II builds the positive case. Over six chapters, we establish through formal proof, simulation, and case evidence that cooperation is not merely normatively attractive but mathematically superior to competition under conditions that are regularly met in real economies. We then develop the peer-to-peer framework — the organizational architecture that embodies cooperative principles at network scale — and extend it through the formal analysis of flat hierarchies, agent-based simulation, and the economics of distributed ledgers.

The mathematical level rises in this Part. Chapter 3 introduced the core and the Shapley value; Chapter 6 proves the Bondareva–Shapley theorem in full and derives the nucleolus. Chapter 5 introduced replicator dynamics; Chapter 7 formalizes stigmergic coordination and proves evolutionary stability conditions for cooperative norms. Chapter 4 introduced algebraic connectivity; Chapter 8 applies it to the design of resilient peer-to-peer production systems. The tools from Books 1 and 2 — optimization [M:Ch.1], linear algebra [M:Ch.2], differential equations [M:Ch.4], and dynamic programming [M:Ch.15] — are used throughout without re-derivation; readers who need a refresher should consult Appendix A.

The central result of Part II, proved progressively across the six chapters and synthesized in the unified model of Part VI, is the Cooperative Advantage Theorem: under the conditions of superadditivity, repeated interaction, and network reciprocity that characterize most real economic relationships, cooperative institutions achieve higher total welfare, greater resilience, and more equitable distribution than their competitive counterparts. This is not a claim about ideal conditions. It is a claim about the real economy — grounded in the mathematics of game theory, network science, and evolutionary dynamics, and confirmed by the empirical record of cooperative enterprises, peer-to-peer platforms, and commons-based institutions worldwide.

---

# Chapter 6: Cooperative Games and the Core — Mathematical Conditions for Stable Cooperation

> *"A theory is exactly as good as its ability to handle limiting cases."*
> — attributed to Paul Samuelson

> *"The question is not whether men will cooperate, but on what terms."*
> — Elinor Ostrom, *Governing the Commons* (1990)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Define balanced collections and state the Bondareva–Shapley theorem precisely; reproduce its proof.
2. Prove that convex cooperative games always have non-empty cores, and identify the economic conditions that generate convexity.
3. Define the nucleolus and compute it for small games using the lexicographic excess-minimization procedure.
4. Apply the core, Shapley value, and nucleolus to cost allocation problems in airport operations, river pollution abatement, and cooperative enterprises.
5. Analyze the Rochdale Principles as a cooperative game design and explain formally why they promote core stability.

---

## 6.1 The Core Revisited: Blocking Coalitions and Stability

Chapter 3 introduced the core as the set of allocations from which no coalition would defect. We gave an informal proof of the Bondareva–Shapley theorem and computed the core for a three-player supply chain. In this chapter we develop the theory in full, proving the theorem rigorously, extending the analysis to the nucleolus, and grounding everything in economic applications of sufficient complexity to be practically useful.

We begin by sharpening the concept of the core.

**Definition 6.1 (Excess and Blocking Coalition).** Given a cooperative game $(N, v)$ and a payoff vector $x \in \mathbb{R}^n$ with $\sum_{i \in N} x_i = v(N)$, the *excess* of coalition $S$ at $x$ is:

$$e(S, x) = v(S) - \sum_{i \in S} x_i$$

The excess measures how much better coalition $S$ could do by defecting from $x$ and claiming its stand-alone value. If $e(S, x) > 0$, coalition $S$ has an incentive to block the allocation $x$: its members can collectively do better. If $e(S, x) \leq 0$ for all $S \subseteq N$, then no coalition has an incentive to defect — $x$ is in the core.

**Definition 6.2 (Core, restated).** The core $\mathcal{C}(v)$ is:

$$\mathcal{C}(v) = \left\{x \in \mathbb{R}^n : \sum_{i \in N} x_i = v(N) \text{ and } e(S, x) \leq 0 \text{ for all } S \subseteq N\right\}$$

The core is a convex polytope (the intersection of finitely many half-spaces with a hyperplane). It may be empty, a single point, or a high-dimensional set. Understanding when it is non-empty — and, when it is, how to select among its elements — is the central problem of cooperative solution theory.

Three observations about the core deserve emphasis before we proceed to the existence theorem.

**First**, the core is a stability concept, not a fairness concept. An allocation in the core may be highly unequal: a player with high outside options may claim a large share of the total surplus, leaving others little more than their individual values. The core guarantees that no coalition can do better by defecting; it says nothing about whether the outcome is just. The Shapley value and nucleolus, developed later in this chapter, add fairness criteria to the stability requirement.

**Second**, the core is sensitive to the characteristic function. Small changes in $v(S)$ for key coalitions can shift the core significantly, empty it, or dramatically enlarge it. This sensitivity makes the design of cooperative institutions — choosing the rules that determine what each coalition can achieve — a powerful lever for governance, as we explore in the case study on the Rochdale Principles.

**Third**, core stability is not the same as dynamic stability. An allocation in the core is stable against one-shot defection; it may not be stable against sequential or multilateral renegotiation. Dynamic stability requires the additional apparatus of the repeated game [C:Ch.3] or the evolutionary dynamics [C:Ch.7]. In this chapter we focus on the static core; the dynamic stability of cooperative institutions is developed in Chapter 30.

---

## 6.2 Balanced Collections and the Bondareva–Shapley Theorem

The Bondareva–Shapley theorem provides the complete characterization of core non-emptiness. Its proof requires the concept of a balanced collection of coalitions.

**Definition 6.3 (Balanced Collection).** A collection of coalitions $\mathcal{B} = \{S_1, S_2, \ldots, S_k\} \subseteq 2^N \setminus \{\emptyset\}$ is *balanced* if there exist positive weights $\{\lambda_S\}_{S \in \mathcal{B}}$ such that for every player $i \in N$:

$$\sum_{S \in \mathcal{B}: i \in S} \lambda_S = 1$$

The weights $\lambda_S$ can be interpreted as the fraction of their time that the members of coalition $S$ devote to coalition $S$, with each player's total time summing to one. A balanced collection is a feasible time-sharing arrangement across coalitions.

**Examples.** For $N = \{1, 2, 3\}$:
- The grand coalition alone: $\mathcal{B} = \{N\}$ with $\lambda_N = 1$. Balanced trivially.
- All singleton coalitions: $\mathcal{B} = \{\{1\}, \{2\}, \{3\}\}$ with $\lambda_{\{i\}} = 1$. Balanced.
- The three two-player coalitions: $\mathcal{B} = \{\{1,2\}, \{1,3\}, \{2,3\}\}$ with $\lambda_S = 1/2$ for each. Check: player 1 is in $\{1,2\}$ and $\{1,3\}$, so $1/2 + 1/2 = 1$ ✓. Balanced.

**Definition 6.4 (Balanced Game).** A game $(N, v)$ is *balanced* if for every balanced collection $\mathcal{B}$ with weights $\{\lambda_S\}$:

$$\sum_{S \in \mathcal{B}} \lambda_S v(S) \leq v(N)$$

The balancedness condition says: the grand coalition produces at least as much as any time-sharing arrangement among sub-coalitions. It is the formal expression of the idea that full cooperation is weakly more valuable than any partial arrangement.

**Theorem 6.1 (Bondareva–Shapley Theorem, 1963).** The core of a TU game $(N, v)$ is non-empty if and only if the game is balanced.

*Proof.* We prove both directions.

**($\Rightarrow$) If the core is non-empty, then the game is balanced.**

Let $x \in \mathcal{C}(v)$ and let $\mathcal{B}$ be any balanced collection with weights $\{\lambda_S\}$. Since $x$ is in the core, $\sum_{i \in S} x_i \geq v(S)$ for all $S$. Multiplying by $\lambda_S$ and summing over $S \in \mathcal{B}$:

$$\sum_{S \in \mathcal{B}} \lambda_S v(S) \leq \sum_{S \in \mathcal{B}} \lambda_S \sum_{i \in S} x_i = \sum_{i \in N} x_i \underbrace{\sum_{S \in \mathcal{B}: i \in S} \lambda_S}_{= 1} = \sum_{i \in N} x_i = v(N)$$

where we used the definition of balancedness (weights sum to 1 for each player) and efficiency of $x$. This holds for all balanced collections, so the game is balanced. $\square$

**($\Leftarrow$) If the game is balanced, then the core is non-empty.**

Consider the linear program:

$$\text{(P)} \quad \min_{x} \; 0 \quad \text{subject to} \quad \sum_{i \in S} x_i \geq v(S) \; \forall S \subseteq N, \quad \sum_{i \in N} x_i = v(N)$$

The core is non-empty if and only if (P) is feasible. By LP duality, (P) is feasible if and only if its dual (D) has a bounded optimum.

The dual of (P) is:

$$\text{(D)} \quad \max_{\lambda} \sum_{S \subseteq N} \lambda_S v(S) \quad \text{subject to} \quad \sum_{S: i \in S} \lambda_S = 1 \; \forall i \in N, \quad \lambda_S \geq 0 \; \forall S \subseteq N$$

The dual constraint says that $\{\lambda_S\}$ with $\sum_{S: i \in S} \lambda_S = 1$ for all $i$ defines a balanced collection. The dual objective is $\sum_S \lambda_S v(S)$. If the game is balanced, then $\sum_S \lambda_S v(S) \leq v(N)$ for every feasible dual solution. Therefore the dual objective is bounded above by $v(N)$, the dual has an optimal solution, and by LP duality, the primal (P) is feasible — that is, the core is non-empty. $\square$

The Bondareva–Shapley theorem reduces core non-emptiness to a linear programming condition, which is both theoretically clean and computationally tractable. For any specific game, one can verify balancedness (and hence core non-emptiness) by solving a linear program with $2^n - 1$ constraints and $n$ variables — manageable for small $n$, requiring computational tools for large coalitional games.

---

## 6.3 Convex Games and Core Non-Emptiness

The balancedness condition of the Bondareva–Shapley theorem is necessary and sufficient but not always easy to verify. A sufficient condition that is more directly interpretable economically is convexity.

**Definition 6.5 (Convex Game, restated).** A game $(N, v)$ is *convex* (or *supermodular*) if for all $S, T \subseteq N$:

$$v(S \cup T) + v(S \cap T) \geq v(S) + v(T)$$

This is equivalent to the condition that marginal contributions are non-decreasing in coalition size [Definition 3.3], and to the condition that $v$ is a supermodular set function.

**Theorem 6.2 (Shapley, 1971 — Convex Games Have Non-Empty Cores).** Every convex game has a non-empty core. Moreover, the Shapley value of a convex game always lies in the core.

*Proof.* We show that the Shapley value $\phi(v)$ lies in the core of any convex game.

For any coalition $S \subseteq N$, we must show $\sum_{i \in S} \phi_i(v) \geq v(S)$.

Fix an arbitrary ordering $\sigma$ of $S = \{s_1, s_2, \ldots, s_{|S|}\}$. By convexity, the marginal contribution of $s_k$ to the coalition $\{s_1, \ldots, s_{k-1}\}$ is at least as large as the marginal contribution of $s_k$ to the subset $\{s_1, \ldots, s_{k-1}\} \cap T$ for any $T$:

$$v(\{s_1, \ldots, s_k\}) - v(\{s_1, \ldots, s_{k-1}\}) \geq v(S_k) - v(S_{k-1})$$

where $S_k = \{s_1, \ldots, s_k\}$. Summing over $k = 1, \ldots, |S|$:

$$\sum_{i \in S}\left[v(S_{\sigma(i)}) - v(S_{\sigma(i)} \setminus \{i\})\right] = v(S) - v(\emptyset) = v(S)$$

This telescoping sum gives the marginal contributions along the ordering $\sigma$. Since the Shapley value averages marginal contributions across all orderings of $N$, and convexity ensures that marginal contributions of members of $S$ are at least as large when we restrict to the ordering $\sigma$ within $S$ as in any other ordering, we obtain $\sum_{i \in S} \phi_i(v) \geq v(S)$. Since efficiency holds by definition ($\sum_{i \in N} \phi_i(v) = v(N)$), we have $\phi(v) \in \mathcal{C}(v)$.

Since $\mathcal{C}(v)$ is non-empty (it contains $\phi(v)$), the game has a non-empty core. $\square$

**Economic conditions generating convexity.** Convexity arises naturally in three important economic settings:

*Complementarities.* When agents' skills or resources are complementary — each addition makes others' contributions more valuable — the production function exhibits supermodularity, which translates to a convex characteristic function. Multi-skilled teams, knowledge-sharing cooperatives, and ecosystems with keystone species all exhibit this structure.

*Network externalities.* When the value of participation increases with the number of participants — platforms, communication networks, standards — the characteristic function is convex in membership. Each additional member increases the marginal value of the next.

*Increasing returns to scale.* When production exhibits increasing returns — fixed costs spread over more output, learning-by-doing, economies of scope — larger coalitions are disproportionately more productive, generating convexity.

In contrast, convexity fails when coalitions face capacity constraints, when there is redundancy among members (two firms with the same specialty add less than twice the value of one), or when coordination costs rise superlinearly with coalition size. Understanding which conditions generate convexity is essential for the design of cooperative institutions: convex settings support large, stable coalitions; non-convex settings require more carefully engineered sharing rules.

---

## 6.4 The Nucleolus: The Most Stable Point in the Core

When the core is non-empty and large — as it often is for convex games — additional selection criteria are needed. The Shapley value provides one such criterion (fairness as average marginal contribution). The nucleolus provides another, from a different angle: rather than asking what is fair, it asks what is most stable.

**Definition 6.6 (Excess Vector).** Given a game $(N,v)$ and an allocation $x$, the *excess vector* $\theta(x) \in \mathbb{R}^{2^n - 2}$ lists the excesses $e(S, x) = v(S) - \sum_{i \in S} x_i$ for all non-trivial coalitions $S \subset N$, arranged in non-increasing order:

$$\theta_1(x) \geq \theta_2(x) \geq \cdots \geq \theta_{2^n - 2}(x)$$

The largest excess $\theta_1(x)$ is the complaint of the most aggrieved coalition at allocation $x$ — how much better off its members would be if they defected and claimed their coalition value.

**Definition 6.7 (Nucleolus).** The *nucleolus* $\eta(v)$ is the unique allocation in the core (if the core is non-empty; otherwise in the pre-imputation set) that lexicographically minimizes the excess vector:

$$\eta(v) = \arg \min_x^{\text{lex}} \; \theta(x)$$

where $\min^{\text{lex}}$ means: first minimize $\theta_1(x)$; among all $x$ achieving the minimum $\theta_1$, minimize $\theta_2(x)$; and so on.

The nucleolus minimizes the maximum complaint first, then the second-largest complaint, and so on — a lexicographic minimax procedure. It is the allocation that makes the most dissatisfied coalition as satisfied as possible, then does the same for the next most dissatisfied coalition, iterating until the allocation is uniquely determined.

**Theorem 6.3 (Schmeidler, 1969 — Existence and Uniqueness).** For any game $(N, v)$ with a non-empty imputation set, the nucleolus exists and is unique.

*Proof sketch.* At each stage of the lexicographic minimization, the feasible set is a non-empty compact convex set (the intersection of the efficiency constraint with constraints on excess values achieved in previous stages). The minimum of a continuous function on a compact set is attained; convexity ensures that each minimum is achieved on a convex subset, and the process terminates in a unique point within finitely many stages. $\square$

**Properties of the nucleolus:**

1. **Core membership:** If $\mathcal{C}(v) \neq \emptyset$, then $\eta(v) \in \mathcal{C}(v)$.
2. **Uniqueness:** Unlike the core (which may be a polytope) and unlike the Shapley value (which is unique but not always in the core), the nucleolus is always both unique and in the core when the core is non-empty.
3. **Coalitional rationality:** It satisfies coalitional rationality by construction (since it minimizes the maximum excess, which is non-positive for core allocations).
4. **Covariance:** Like the Shapley value, the nucleolus is covariant under positive affine transformations of the characteristic function.

**Computing the nucleolus.** For small games, the nucleolus can be computed through a sequence of linear programs. At stage 1:

$$\min_{x, \varepsilon} \; \varepsilon \quad \text{subject to} \quad \sum_{i \in N} x_i = v(N), \quad v(S) - \sum_{i \in S} x_i \leq \varepsilon \; \forall S \subsetneq N$$

This finds the minimum achievable maximum excess $\varepsilon_1^*$. At stage 2, fix $\varepsilon = \varepsilon_1^*$ for all coalitions achieving this maximum, and minimize the next-largest excess over the remaining coalitions. Continue until the allocation is unique. For the nucleolus of games with up to $n = 6$ players, this sequential LP approach is computationally tractable; for larger games, specialized algorithms exist (Kohlberg, 1971).

**Economic interpretation.** The nucleolus is the allocation a mediator would propose if their objective were to minimize social conflict — to make it as hard as possible for any coalition to make a credible complaint about its share of the surplus. It is not necessarily the fairest allocation, nor the most efficient path to agreement, but it is the most resistant to destabilization. In settings where coalition formation is a genuine threat — international negotiations, labor-management bargaining, multi-stakeholder consortia — the nucleolus provides a principled focal point for negotiation.

---

## 6.5 Mathematical Model: Cost Allocation as a Cooperative Game

A recurring theme in applied cooperative game theory is cost allocation: when a group of agents jointly uses an infrastructure or service, how should the total cost be shared? Cost allocation games are the mirror image of profit-sharing games — the characteristic function represents costs rather than benefits, and the solution concepts minimize rather than maximize.

**Definition 6.8 (Cost Allocation Game).** A cost allocation game $(N, c)$ is a cooperative game in which $c(S)$ represents the minimum total cost of providing a service to the coalition $S$. A cost allocation $x$ is in the core if:

$$\sum_{i \in N} x_i = c(N) \quad \text{and} \quad \sum_{i \in S} x_i \leq c(S) \; \forall S \subseteq N$$

The core condition now requires that no coalition pays more than it would cost them to operate independently: $\sum_{i \in S} x_i \leq c(S)$.

**Proposition 6.1 (Concave Cost Games Have Non-Empty Cores).** If $c(S)$ is submodular (concave) — meaning $c(S \cup T) + c(S \cap T) \leq c(S) + c(T)$ for all $S, T$ — then the cost allocation game has a non-empty core.

*Proof.* The game $v(S) = c(N) - c(N \setminus S)$ transforms the cost game into an equivalent profit game. Submodularity of $c$ implies supermodularity of $v$ — convexity — and Theorem 6.2 then guarantees a non-empty core. $\square$

Submodularity of costs arises naturally when there are economies of scale: the marginal cost of serving an additional agent decreases as the coalition grows. Airport runway costs, electricity grid infrastructure, shared research facilities, and municipal water systems all exhibit this property.

---

## 6.6 Worked Example: Airport Cost Allocation

The airport game, introduced by Littlechild and Owen (1973), is among the most studied cost allocation problems in cooperative game theory. It provides a clean, empirically grounded illustration of the core, Shapley value, and nucleolus.

**Setup.** An airport is used by four airlines: $A$ (small regional carrier), $B$ (medium-haul carrier), $C$ (long-haul carrier), and $D$ (international heavy carrier). Each airline requires a runway of minimum length to land its aircraft. The costs of constructing runways of increasing length are:

| Runway length | Required by | Marginal cost ($\times 10^3$) | Cumulative cost ($\times 10^3$) |
|:---|:---|:---|:---|
| Short | A only | 100 | 100 |
| Medium | B and above | 80 | 180 |
| Long | C and above | 120 | 300 |
| Full | D only | 150 | 450 |

The key feature of the airport game is that a runway that serves a larger aircraft automatically serves all smaller aircraft. Therefore, the cost of providing service to any coalition $S$ equals the cost of the longest runway required by any member of $S$.

**Characteristic function.** Let $c(S)$ denote the cost for coalition $S$:

$$c(\{A\}) = 100, \quad c(\{B\}) = 180, \quad c(\{C\}) = 300, \quad c(\{D\}) = 450$$
$$c(\{A,B\}) = 180, \quad c(\{A,C\}) = 300, \quad c(\{A,D\}) = 450$$
$$c(\{B,C\}) = 300, \quad c(\{B,D\}) = 450, \quad c(\{C,D\}) = 450$$
$$c(\{A,B,C\}) = 300, \quad c(\{A,B,D\}) = 450, \quad c(\{A,C,D\}) = 450, \quad c(\{B,C,D\}) = 450$$
$$c(\{A,B,C,D\}) = 450$$

This is submodular (concave): $c(S \cup T) + c(S \cap T) \leq c(S) + c(T)$ because the cost is determined by the maximum type in the coalition, and the max function is submodular. By Proposition 6.1, the core is non-empty.

**Step 1: Shapley value computation.**

For the cost game, the Shapley value allocates to each player their average marginal contribution to cost:

$$\phi_i(c) = \sum_{S \subseteq N \setminus \{i\}} \frac{|S|!(n-|S|-1)!}{n!}\left[c(S \cup \{i\}) - c(S)\right]$$

Player $A$'s marginal cost contributions: $A$ only adds cost when joining a coalition that has no other member (i.e., $A$ is the first and only small carrier). In all coalitions already containing $B$, $C$, or $D$, adding $A$ contributes zero marginal cost (the runway already covers $A$).

Marginal contributions of $A$:
- To $\emptyset$: $c(\{A\}) - c(\emptyset) = 100$
- To $\{B\}$: $c(\{A,B\}) - c(\{B\}) = 180 - 180 = 0$
- To $\{C\}$: $c(\{A,C\}) - c(\{C\}) = 300 - 300 = 0$
- To $\{D\}$: $c(\{A,D\}) - c(\{D\}) = 450 - 450 = 0$
- To $\{B,C\}$: $0$; to $\{B,D\}$: $0$; to $\{C,D\}$: $0$; to $\{B,C,D\}$: $0$

Therefore $\phi_A(c) = \frac{3! \cdot 0!}{4!} \cdot 100 = \frac{1}{4} \cdot 100 = 25$.

By analogous reasoning (marginal cost arises only when the player introduces the need for a new runway tier):

$$\phi_A = 25, \quad \phi_B = \frac{1}{4} \cdot 80 + \frac{1}{12} \cdot 80 = 20 + 6.67 = 26.67$$

More carefully — $B$ introduces the medium runway tier when $B$ is the first member from $\{B,C,D\}$ to join. The probability of this (in a random ordering) is $1/4$ when $B$ enters before $C$, $D$, and $A$ is already present or not; the marginal cost is $80$ when $B$ is the first of $\{B,C,D\}$ and zero otherwise.

The fraction of orderings in which $B$ is first among $\{B, C, D\}$ is $1/3$ (since each of the three is equally likely to be first among them). Therefore:

$$\phi_B = \frac{1}{3} \cdot 80 = 26.67$$

Similarly, $C$ is first among $\{C, D\}$ with probability $1/2$:

$$\phi_C = \frac{1}{2} \cdot 120 = 60$$

And $D$ always introduces the full-length runway (the only member requiring it):

$$\phi_D = 150$$

Check: $\phi_A + \phi_B + \phi_C + \phi_D = 25 + 26.67 + 60 + 150 = 261.67 \neq 450$.

We have made an error — we need to account for the full cost structure more carefully. The Shapley value for cost games should satisfy $\sum_i \phi_i = c(N) = 450$. Let us recompute using the standard formula directly.

The total cost $c(N) = 450$ can be decomposed as: the cost of the short runway ($100$), which is used by all four airlines; the incremental cost of the medium runway ($80$), used by $B$, $C$, and $D$; the incremental cost of the long runway ($120$), used by $C$ and $D$; the incremental cost of the full runway ($150$), used by $D$ alone.

By the Shapley value axioms (specifically, the additivity axiom applied to this decomposition):

$$\phi_A = \frac{100}{4} = 25$$

$$\phi_B = \frac{100}{4} + \frac{80}{3} = 25 + 26.67 = 51.67$$

$$\phi_C = \frac{100}{4} + \frac{80}{3} + \frac{120}{2} = 25 + 26.67 + 60 = 111.67$$

$$\phi_D = \frac{100}{4} + \frac{80}{3} + \frac{120}{2} + \frac{150}{1} = 25 + 26.67 + 60 + 150 = 261.67$$

Check: $25 + 51.67 + 111.67 + 261.67 = 450$ ✓.

The decomposition reveals the economic logic of the Shapley value for airport games: each runway tier is shared equally among the airlines that use it. The full-length runway costs $150$ and is used only by $D$, so $D$ bears the entire $150$. The long runway extension costs $120$ and is used by $C$ and $D$, so each bears $60$. The medium extension costs $80$ and is shared among $B$, $C$, and $D$, giving $26.67$ each. The short runway costs $100$ and is shared equally among all four airlines.

**Step 2: Core.**

The core constraints require that no coalition pays more than its stand-alone cost. Key constraints:

$$\phi_D \leq c(\{D\}) = 450 \; \checkmark \quad (261.67 \leq 450)$$
$$\phi_C + \phi_D \leq c(\{C,D\}) = 450 \; \checkmark \quad (373.33 \leq 450)$$
$$\phi_B + \phi_C + \phi_D \leq c(\{B,C,D\}) = 450 \; \checkmark \quad (425 \leq 450)$$
$$\phi_A \leq c(\{A\}) = 100 \; \checkmark \quad (25 \leq 100)$$

The Shapley value is in the core for this game — consistent with the convexity result.

**Step 3: Nucleolus.**

The nucleolus minimizes the maximum excess. For the airport game, one can show (by the LP procedure) that the nucleolus coincides with the Shapley value for this specific decomposable structure. This is a known result: for airport games, the Shapley value equals the nucleolus.

**Step 4: Comparison and interpretation.**

| Airline | Stand-alone cost | Shapley / Nucleolus | Proportional to use | Equal share |
|:---|:---|:---|:---|:---|
| A (small) | 100 | **25** | 100 | 112.5 |
| B (medium) | 180 | **51.67** | 180 | 112.5 |
| C (long) | 300 | **111.67** | 300 | 112.5 |
| D (heavy) | 450 | **261.67** | 450 | 112.5 |
| **Total** | — | **450** | — | **450** |

The Shapley/nucleolus allocation is dramatically more equitable than proportional-to-use (which would simply charge each airline its stand-alone cost — the same as no cooperation). It is also more equitable than equal shares (which would undercharge $D$ and overcharge $A$). The cooperative solution rewards the small carrier for the cost savings it provides to others by sharing the short runway tier.

This result generalizes: cooperative cost allocation, whenever the characteristic function is submodular, produces sharing arrangements that reward smaller users more equitably than either equal-shares or proportional-to-cost allocations — a property of direct relevance to cooperative infrastructure projects, shared renewable energy facilities, and municipal services.

---

## 6.7 Case Study: The Rochdale Principles as Cooperative Game Design

The Rochdale Society of Equitable Pioneers, founded in 1844 in Rochdale, England by a group of working-class artisans, established the cooperative movement on a set of principles that have since been codified, refined, and adopted by the International Cooperative Alliance as the foundation of cooperative enterprise worldwide. The Rochdale Principles are often presented as ethical commitments; in this case study, we analyze them as cooperative game design — a set of rules that shape the characteristic function $v(S)$ and the core $\mathcal{C}(v)$ of the enterprise game.

**The seven principles:**
1. Voluntary and open membership
2. Democratic member control (one member, one vote)
3. Member economic participation (members contribute equitably to capital)
4. Autonomy and independence
5. Education, training, and information
6. Cooperation among cooperatives
7. Concern for community

We analyze each through the lens of core stability.

**Principle 1 (Open membership) and core stability.** A cooperative with closed membership may be in the core with respect to its current members, but it creates external blocking coalitions: excluded potential members who could form an alternative cooperative that captures some of the surplus. Open membership prevents the formation of competing external coalitions by ensuring that anyone who could credibly threaten the existing cooperative can instead join it.

Formally: if $v_{\text{open}}(S)$ is the characteristic function of an open cooperative (where $S$ can include any potential member) and $v_{\text{closed}}(S)$ is the function for a closed cooperative (where $S$ is restricted to current members), then $v_{\text{open}}(N) \geq v_{\text{closed}}(N)$ whenever new members add positive value. Open membership expands the core by enlarging the feasible set of grand coalition allocations.

**Principle 2 (Democratic control) and fairness.** One member, one vote is not merely an ethical commitment; it is an allocation rule that implements something close to the Shapley value. In a game where all members have equal voting rights, the characteristic function treats members symmetrically (the Shapley symmetry axiom), implying equal Shapley values — and hence equal shares of cooperative surplus — for members with equal contributions. Departures from one-member-one-vote (weighted voting, capital-proportional voting) break symmetry and shift the allocation toward the Shapley value of a less symmetric game, typically benefiting capital-rich members at the expense of labor-rich ones.

**Principle 3 (Member economic participation) and the no-free-rider condition.** Members contribute equitably to the cooperative's capital and share in its surplus in proportion to their participation. This is the formal expression of the null player axiom: a member who contributes nothing (who is a null player in the game) receives no share of the surplus. The principle prevents free-riding on the cooperative's capital by ensuring that those who benefit from the cooperative's assets are also those who invest in them.

**Principle 6 (Cooperation among cooperatives) and coalition formation.** The principle of inter-cooperative cooperation formally enlarges the characteristic function: $v(S \cup T) \geq v(S) + v(T)$ when $S$ and $T$ are cooperatives that cooperate with each other. This is precisely the superadditivity condition, and its fulfillment is not automatic — it requires deliberate institutional investment in inter-cooperative networks, shared services, and federated governance. The Mondragon cooperative federation [C:Ch.34] and the cooperative banking system (Caja Laboral) are empirical instantiations of Principle 6 at scale.

**The Rochdale game.** We can now characterize the Rochdale Principles collectively as a set of rules that:

1. Maximize the size of the feasible grand coalition (Principle 1).
2. Ensure symmetric treatment of members (Principle 2), moving the allocation toward the Shapley value.
3. Eliminate null players from the surplus distribution (Principle 3), satisfying the null player axiom.
4. Make the game superadditive and, when inter-cooperative networks are well-developed, convex (Principle 6).
5. Invest in the long-run productivity of members (Principle 5), increasing $v(N)$ over time and expanding the core.

A cooperative that implements all seven principles fully is — in the formal language of cooperative game theory — a mechanism for constructing a convex game with a Shapley-value allocation. The Rochdale Pioneers arrived at these design principles empirically, through trial and error across the 1840s cooperative movement. The cooperative game theory of Bondareva, Shapley, and Schmeidler, developed a century later, provides the formal foundation that explains why the principles work: they create the conditions — balanced coalitions, symmetric players, null-player exclusion, superadditivity — under which stable, fair cooperation is mathematically guaranteed.

---

## Chapter Summary

This chapter has developed the full theoretical apparatus of the core in cooperative game theory, extending the introduction of Chapter 3 to a complete and rigorous treatment.

The Bondareva–Shapley theorem establishes that a cooperative game has a non-empty core if and only if it is balanced — if no time-sharing arrangement among sub-coalitions can outperform the grand coalition. The proof proceeds by LP duality: core non-emptiness is equivalent to the boundedness of the dual program, which is the formal statement of balancedness. This result connects cooperative game theory to linear programming and provides a computationally tractable method for verifying core existence in specific applications.

Convex games — in which marginal contributions are non-decreasing in coalition size — always have non-empty cores, and their Shapley values always lie within the core. Convexity arises naturally in settings with complementarities, network externalities, and increasing returns — precisely the settings that characterize cooperative enterprises, knowledge economies, and ecological systems.

The nucleolus selects the unique allocation in the core that lexicographically minimizes the excess vector — the allocation that minimizes the complaint of the most aggrieved coalition, then the next most aggrieved, and so on. It is the most stable point in the core in the sense of minimizing the maximum incentive to defect.

The airport cost allocation problem demonstrates the three solution concepts in a concrete, calibrated setting and reveals the economic logic of cooperative cost sharing: costs are shared in proportion to the incremental burden each agent imposes on the infrastructure, distributed equitably across the agents who require each tier of service.

The Rochdale Principles, analyzed as cooperative game design, show that the canonical institutional rules of the cooperative movement — open membership, democratic control, economic participation, inter-cooperative cooperation — are formal implementations of the conditions that guarantee core non-emptiness, Shapley-value fairness, and long-run coalition stability.

Chapter 7 turns to the dynamic foundation of cooperation: the repeated game and its evolutionary extensions. Having established when cooperation is stable (the core), we now ask how it emerges — and how cooperative norms spread and persist in populations of adapting agents.

---

## Exercises

**6.1** Consider a three-player joint venture with characteristic function: $v(\{1\}) = 20$, $v(\{2\}) = 30$, $v(\{3\}) = 10$, $v(\{1,2\}) = 70$, $v(\{1,3\}) = 40$, $v(\{2,3\}) = 50$, $v(\{1,2,3\}) = 100$.
(a) Verify that this game is superadditive.
(b) Check whether it is convex (test the supermodularity condition for all pairs of coalitions).
(c) If the core is non-empty, find one allocation in the core and verify it satisfies all constraints.
(d) Compute the Shapley value. Is it in the core?

**6.2** Explain intuitively why a cooperative enterprise that implements the Rochdale Principle of "cooperation among cooperatives" (Principle 6) tends to make its cooperative game convex. What economic mechanism is responsible for this?

**6.3** In the airport cost allocation example, suppose a fifth airline $E$ joins, requiring a full-length runway identical to $D$. Recompute:
(a) The characteristic function $c(S)$ for all relevant coalitions.
(b) The Shapley value allocation.
(c) Does $E$'s presence reduce or increase the cost burden on $A$, $B$, and $C$? Why?

**★ 6.4** Prove Theorem 6.2: the Shapley value of a convex game lies in the core.

Your proof should:
(a) Fix an arbitrary coalition $S$ and an arbitrary ordering $\sigma$ of $S$.
(b) Show that under convexity, the marginal contribution of $s_k$ to $S_k = \{s_1, \ldots, s_k\}$ is at least the marginal contribution to $S_k \cap T$ for any coalition $T \supseteq S_k$.
(c) Use this to show $\sum_{i \in S} \phi_i(v) \geq v(S)$ for all $S$.
(d) Conclude that $\phi(v) \in \mathcal{C}(v)$.

**★ 6.5** A mutual insurance cooperative has 10 members. Member $i$ faces an annual loss $\ell_i$ with probability $p_i$, independent across members. The cost of self-insurance for member $i$ alone is $c_i = \ell_i \cdot p_i + k\sqrt{\ell_i \cdot p_i(1-p_i)}$, where $k > 0$ reflects a risk-loading factor. For any coalition $S$, the cooperative's insurance cost is the expected loss plus the same risk loading applied to the aggregate portfolio:

$$c(S) = \sum_{i \in S} \ell_i p_i + k\sqrt{\sum_{i \in S} \ell_i^2 p_i(1-p_i)}$$

(a) Show that $c(S)$ is submodular (concave) in $S$, and therefore the cost allocation game has a non-empty core.

(b) With $n = 10$ members and parameters $\ell_i = 1000$, $p_i = 0.05$, $k = 1.645$ (95th percentile of the normal distribution) for all $i$, compute $c(\{i\})$ for a single member and $c(N)$ for the grand coalition.

(c) Compute the Shapley value allocation for the symmetric case (all members identical). Interpret the result: how much does each member save by joining the cooperative?

(d) Now suppose member 10 has a higher risk: $p_{10} = 0.15$, $\ell_{10} = 1000$. Recompute the Shapley value. Does the high-risk member pay more or less than in the symmetric case? Is the allocation still in the core?

**★★ 6.6** The nucleolus of the airport game (Section 6.6) coincides with the Shapley value. This is a known result for airport games, but it does not hold in general.

(a) Construct a four-player game where the nucleolus and Shapley value differ. (Hint: try a game where one coalition has a much higher excess at the Shapley value than at the nucleolus.)

(b) Implement the sequential LP algorithm for computing the nucleolus for your game. Report the nucleolus and verify it is in the core.

(c) Interpret the difference between the nucleolus and the Shapley value for your game: which allocation would a risk-averse cooperative prefer, and why?

(d) Design a ten-firm circular supply chain cooperative in which each firm supplies to two downstream firms and receives inputs from two upstream firms. Specify the characteristic function, verify core non-emptiness, compute the Shapley value, and discuss whether the Shapley value or nucleolus is the more appropriate allocation rule for this governance structure.

---

*Chapter 7 extends the analysis from static stability to dynamic emergence. Having proved when cooperation is stable, we now ask how it arises — how cooperative norms spread through populations of agents who are not assumed to be rational in the full game-theoretic sense, but who adapt their strategies based on what works. The bridge between game theory and evolutionary dynamics turns out to run directly through the mechanism of stigmergy: the distributed, environment-mediated coordination that Chapter 2 identified as the third engine of economic organization.*
