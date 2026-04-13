# Chapter 25: Mutual Credit and Clearing Systems — Mathematical Models of Decentralized Liquidity

> *"Money is not a thing; it is a relationship. Mutual credit makes that relationship visible."*
> — Thomas Greco, *The End of Money and the Future of Civilization* (2009)

> *"The WIR franc is not inflationary because it is always backed by real goods and services in the trading network."*
> — James Stodder, *Reciprocal Exchange Networks: Implications for Macroeconomic Stability* (2009)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Specify the double-entry mechanics of mutual credit formally — showing how credit commitments between trading partners create liquidity without prior savings or central authority — and represent the system as a directed weighted graph with balances as node weights.
2. Formalize the clearing problem as a minimum-cost flow linear program, and prove that multilateral clearing achieves weakly greater liquidity efficiency than bilateral netting for any network configuration.
3. Prove the Optimal Clearing Theorem: that a mutual credit system with complete multilateral clearing achieves the same real allocation as a monetary economy with perfect credit markets, without requiring any central money creation.
4. Analyze the stability and resilience of mutual credit networks using Lyapunov methods, identifying the network topologies that sustain mutual credit under member failure and credit limit violations.
5. Formalize Commitment Pooling — the four-function model of community currency issuance — and analyze the Sarafu Network as a case of commitment pooling in practice.
6. Apply the framework to the Swiss WIR system, assessing its counter-cyclical stability properties through the formal model.

---

## 25.1 The Mutual Credit Insight

The preceding two chapters developed monetary alternatives that operate at the system level: sovereign money requires a central bank to create and manage all money; even the CBDC transition requires coordinated state action. Mutual credit is different in kind: it requires no central authority, no prior savings, and no pre-existing money. It is money created at the point of transaction, by the parties to the transaction, from nothing but their mutual commitment to future delivery.

The insight is ancient — older than coinage — but the formal theory is modern. When a baker promises to deliver bread in exchange for a carpenter's promise to deliver furniture, each has extended credit to the other. If these commitments are made on a common ledger, and if that ledger is sufficiently connected, the commitments can circulate as a medium of exchange: the baker's promise can be used to pay a third party, who uses it to pay a fourth, and so on. No gold, no government, no bank is required. Only a ledger and a community willing to honor commitments.

The formal theory of mutual credit connects three bodies of theory developed earlier in this book: cooperative game theory (Chapter 6 — mutual credit is a cooperative allocation mechanism), P2P network theory (Chapter 8 — mutual credit is decentralized coordination), and the clearing problem (a minimum-cost flow problem, formalizing the intuition that optimal clearing eliminates as much bilateral debt as possible through multilateral offset). This chapter unifies these connections into a coherent formal framework and applies it to the oldest and most successful mutual credit system currently operating: Switzerland's WIR franc, 90 years old and still counter-cyclically stabilizing the Swiss economy.

---

## 25.2 Mutual Credit Mechanics: Double-Entry Credit

### 25.2.1 The Accounting Structure

**Definition 25.1 (Mutual Credit System).** A mutual credit system is a tuple $(\mathcal{M}, \mathcal{L}, \mathcal{R}, b, \ell)$ where:
- $\mathcal{M} = \{1, 2, \ldots, n\}$: the set of member agents.
- $\mathcal{L}$: the shared ledger — a distributed record of all credit relationships.
- $\mathcal{R}$: the rule set governing credit limits, clearing, and membership (the Ostrom governance structure [C:Ch.14]).
- $b_i \in \mathbb{R}$: member $i$'s current balance (positive = creditor, negative = debtor).
- $\ell_i > 0$: member $i$'s credit limit — the maximum negative balance permitted.

**System constraint.** Mutual credit is by construction zero-sum: for every debit there is a corresponding credit. The system balance constraint is:

$$\sum_{i \in \mathcal{M}} b_i = 0 \quad \text{at all times}$$

This is the mutual credit equivalent of the SFC row-sum identity: every credit in the system has a corresponding debit. The system creates no net purchasing power — it redistributes existing commitments.

**Transaction mechanics.** When member $i$ purchases goods worth $v$ from member $j$:
- $i$'s balance: $b_i \to b_i - v$ (debited; $i$ owes more to the system)
- $j$'s balance: $b_j \to b_j + v$ (credited; $j$ is owed more by the system)
- System balance: $b_i + b_j$ unchanged (net change = $-v + v = 0$)

The transaction creates no new money — it transfers a claim. However, it does create liquidity: member $i$ was able to purchase without prior money, and member $j$ received a claim redeemable from any network member, not only from $i$. This is the distributed liquidity provision that distinguishes mutual credit from barter.

**Proposition 25.1 (Mutual Credit Creates Liquidity Without Money Creation).** In a mutual credit system with $n$ members and credit limits $\{\ell_i\}$, the total liquidity available for trade is:

$$\mathcal{L}_{\text{total}} = \sum_{i \in \mathcal{M}} \ell_i$$

without requiring any pre-existing monetary stock. A monetary economy achieving the same transaction volume requires a money supply of at least $M \geq \mathcal{L}_{\text{total}}/V$ where $V$ is monetary velocity.

*Proof.* In the mutual credit system, each member can purchase up to $\ell_i$ of goods from others without prior accumulation — credit is extended at the point of transaction. Total trade capacity is $\sum_i \ell_i$. In a monetary economy, the same trade volume requires the money stock to turn over $V$ times: $M \cdot V = PY \geq \mathcal{L}_{\text{total}}$, so $M \geq \mathcal{L}_{\text{total}}/V$. Mutual credit achieves the same liquidity with $M = 0$. $\square$

### 25.2.2 Comparison with Debt-Based and Sovereign Money

| Feature | Debt-based money | Sovereign money | Mutual credit |
|:---|:---|:---|:---|
| Money creator | Commercial banks | Central bank | Members at transaction |
| Prior savings required | No (banks create) | No (CB creates) | No |
| Growth imperative | Yes (Theorem 23.3) | No | No |
| Minsky instability | Structural | Eliminated | Absent by construction |
| Seigniorage beneficiary | Banking sector | Public | None (no seigniorage) |
| Governance | Regulatory | Democratic | Commons (Ostrom) |
| Scale | National/global | National | Local/regional |
| Systemic risk | High | Low | Low/distributed |

Mutual credit has no seigniorage — no party benefits from the creation of the credit medium itself. This is its most distinctive property: the monetary system is a pure coordination device, not a source of rent. It is also its principal limitation: mutual credit cannot fund large capital investments (which require accumulated purchasing power beyond current credit commitments), and its resilience depends on the density and trust-fulness of the trading network.

---

## 25.3 The Mutual Credit Network: Graph-Theoretic Model

### 25.3.1 The Directed Weighted Graph

**Definition 25.2 (Mutual Credit Network).** The mutual credit network at time $t$ is a directed weighted graph $G^{MC}(t) = (V, E, w, b)$ where:
- $V = \mathcal{M}$ (members are nodes).
- $E \subseteq V \times V$ (directed edges represent credit relationships or pending transactions).
- $w_{ij} \geq 0$: the amount of outstanding credit owed by member $i$ to member $j$.
- $b_i = \sum_j w_{ji} - \sum_j w_{ij}$: member $i$'s net balance (credits received minus credits owed).

**The balance vector** $\mathbf{b} = (b_1, b_2, \ldots, b_n)^T$ represents the current state of the mutual credit system. It satisfies $\mathbf{1}^T \mathbf{b} = 0$ (system balance constraint). The balance vector is related to the incidence matrix of the credit graph:

$$\mathbf{b} = B \mathbf{w}$$

where $B$ is the $n \times |E|$ signed incidence matrix (columns correspond to edges; entry $+1$ for the head, $-1$ for the tail) and $\mathbf{w}$ is the vector of edge weights.

**Example 25.1 (Three-Member Network).** Members: Baker (B), Carpenter (C), Farmer (F). Outstanding obligations: C owes B 100 units ($w_{CB} = 100$), B owes F 70 units ($w_{BF} = 70$), F owes C 50 units ($w_{FC} = 50$).

Balance vector:
- $b_B = w_{CB} - w_{BF} = 100 - 70 = +30$ (net creditor)
- $b_C = w_{FC} - w_{CB} = 50 - 100 = -50$ (net debtor)
- $b_F = w_{BF} - w_{FC} = 70 - 50 = +20$ (net creditor)

System check: $b_B + b_C + b_F = 30 - 50 + 20 = 0$. ✓

The circular structure ($C \to B \to F \to C$) permits clearing: if all obligations are settled simultaneously, every member's balance returns toward zero, releasing credit capacity for new transactions.

---

## 25.4 The Clearing Problem

### 25.4.1 Why Clearing Matters

Gross obligations in a mutual credit network exceed net obligations: in Example 25.1, gross obligations total $100 + 70 + 50 = 220$ units, while net obligations (what would remain after all circularity is cleared) total only $|b_B| + |b_C| + |b_F|/2 = (30 + 50 + 20)/2 = 50$ units. Clearing — the process of offsetting circular obligations against each other — reduces the effective debt burden by $220 - 2 \times 50 = 120$ units (77\%) without anyone paying anyone anything.

This clearing efficiency is the formal foundation of the mutual credit advantage: networks of trade obligations can be dramatically reduced through multilateral offset, releasing credit capacity that would otherwise be tied up in gross positions.

### 25.4.2 The Minimum-Cost Flow Formulation

**Definition 25.3 (Clearing Problem).** Given a mutual credit network $G^{MC} = (V, E, w, b)$, the clearing problem is: find a flow $f: E \to \mathbb{R}_{\geq 0}$ that reduces all balances toward zero at minimum total transaction cost:

$$\min_{f} \sum_{(i,j) \in E} c_{ij} f_{ij}$$

subject to:
$$\sum_j f_{ij} - \sum_j f_{ji} = -b_i \quad \forall i \in V \quad \text{(balance reduction)} \tag{BR}$$
$$0 \leq f_{ij} \leq w_{ij} \quad \forall (i,j) \in E \quad \text{(flow bounds)} \tag{FB}$$
$$\sum_j f_{ij} - \sum_j f_{ji} \geq -b_i - \varepsilon_i \quad \forall i \in V \quad \text{(partial clearing tolerance)} \tag{PC}$$

where $c_{ij} \geq 0$ is the unit transaction cost on edge $(i,j)$ and $\varepsilon_i \geq 0$ allows partial clearing when full clearing is infeasible (e.g., when some debtor lacks sufficient creditors in the network).

**This is a minimum-cost flow problem** [C:Ch.21 — same structure as circular economy optimization], solvable in polynomial time by network simplex or successive shortest paths algorithms.

**Theorem 25.1 (Clearing Feasibility).** The clearing problem has a feasible solution (achieving $f$ that reduces all balances to within $\varepsilon_i$ of zero) if and only if the mutual credit network is strongly connected — there exists a directed path from every member to every other member.

*Proof.* Strong connectivity ensures that for any debtor $i$ and creditor $j$, there exists a path of obligations from $i$ to $j$ along which flow can be routed. In a non-strongly-connected network, some debtors have no path to creditors and cannot be cleared multilaterally — bilateral settlement (or default) is required. $\square$

**Corollary 25.1 (Clearing Efficiency and Network Density).** The fraction of gross obligations that can be cleared multilaterally is increasing in the algebraic connectivity $\lambda_2(L_{G^{MC}})$ of the mutual credit network — denser, better-connected networks achieve greater clearing efficiency.

This connects to Chapter 12's result: networks with higher Fiedler value are more resilient and, in the mutual credit context, more efficient at clearing. Designing mutual credit networks with high algebraic connectivity is therefore simultaneously a resilience strategy and a clearing-efficiency strategy.

### 25.4.3 The Optimal Clearing Theorem

**Theorem 25.2 (Optimal Clearing Theorem).** A mutual credit system with complete multilateral clearing achieves the same real allocation of goods and services as a monetary economy with perfect credit markets (no transaction costs, no credit rationing, no information asymmetry), without requiring any central money creation.

*Proof.* 

**Step 1: Equivalence of real allocations.** In a monetary economy with perfect credit markets, any agent $i$ can purchase goods worth up to their credit limit $\ell_i$ by borrowing. The resulting real allocation is determined by agents' preferences and credit limits, not by prior money holdings. The same is true in mutual credit: any member $i$ can purchase goods worth up to $\ell_i$ from other members, with the debit creating a credit elsewhere in the network. The set of feasible real allocations is identical.

**Step 2: Multilateral clearing achieves feasibility.** After a sequence of transactions, the resulting balance vector $\mathbf{b}$ may have some members near their credit limits (preventing further purchases). Multilateral clearing reduces all balances toward zero, releasing credit capacity — enabling the next round of transactions that would be infeasible with uncleaned balances. Complete multilateral clearing (solving the clearing LP to optimality) achieves the maximum release of credit capacity.

**Step 3: No money required.** The clearing process operates on the ledger of credit commitments — it does not require the transfer of any money, only the simultaneous offset of mutual obligations. The real economy achieves its efficient allocation through credit commitments and their clearing, with money playing no necessary role. $\square$

**Theorem 25.3 (Multilateral Dominates Bilateral Clearing).** For any mutual credit network and any balance vector $\mathbf{b}$, multilateral clearing achieves weakly greater clearing efficiency than bilateral netting:

$$\text{Cleared}^{\text{multi}} \geq \text{Cleared}^{\text{bilateral}} \quad \text{(in terms of gross obligations reduced)}$$

with strict inequality whenever the network contains at least one cycle of length $\geq 3$.

*Proof.* Bilateral netting reduces each pair $(i,j)$'s obligations to their net bilateral position: it replaces $\{w_{ij}, w_{ji}\}$ with $\{\max(0, w_{ij} - w_{ji}), 0\}$ or $\{0, \max(0, w_{ji} - w_{ij})\}$. This eliminates circular obligations of length 2 but cannot eliminate cycles of length $\geq 3$. Multilateral clearing (solving the clearing LP) finds all cycles of any length and eliminates them — strictly more efficient whenever any cycle of length $\geq 3$ exists in the obligation graph. $\square$

---

## 25.5 Stability and Resilience of Mutual Credit Networks

### 25.5.1 Lyapunov Stability Analysis

**Definition 25.4 (Mutual Credit System State).** The state of the mutual credit system at time $t$ is the balance vector $\mathbf{b}(t) \in \mathbb{R}^n$ with $\mathbf{1}^T \mathbf{b} = 0$.

**Candidate Lyapunov function.** A natural measure of system stress is the sum of squared deviations from zero balance — the "imbalance energy":

$$V(\mathbf{b}) = \frac{1}{2} \|\mathbf{b}\|^2 = \frac{1}{2} \sum_i b_i^2$$

$V(\mathbf{b}) = 0$ iff all balances are zero (fully cleared). $V(\mathbf{b}) > 0$ iff some members carry non-zero balances.

**Dynamics.** Between clearing events, balances evolve as members transact:

$$\dot{\mathbf{b}} = A(t) \mathbf{f}(t)$$

where $A$ is the signed incidence matrix and $\mathbf{f}(t)$ is the transaction flow vector. During a clearing event at time $t_c$:

$$\mathbf{b}(t_c^+) = \mathbf{b}(t_c^-) - \mathbf{\delta b}_{\text{cleared}}$$

where $\mathbf{\delta b}_{\text{cleared}} \geq 0$ (component-wise: clearing reduces the magnitude of all balances).

**Proposition 25.2 (Lyapunov Stability under Clearing).** The mutual credit system with periodic optimal clearing is Lyapunov stable: $V(\mathbf{b}(t)) \leq V(\mathbf{b}(0))$ for all $t \geq 0$, with $V$ decreasing at each clearing event.

*Proof.* Between clearing events, $V$ may increase or decrease depending on transaction patterns. At each clearing event, the optimal clearing LP minimizes the post-clearing $\|\mathbf{b}\|^2$ subject to the clearing constraints — by definition, $V$ does not increase at a clearing event. The trajectory $\{V(t)\}$ is therefore bounded above by $V(\mathbf{b}(0))$ and decreasing at clearing events. $\square$

**Resilience to member failure.** Suppose member $k$ fails — cannot honor their credit obligations. The impact on the system:

$$\Delta V = \frac{1}{2} b_k^2 \quad \text{(removed from system, redistributed as loss)}$$

The loss is bounded by $b_k^2/2 \leq \ell_k^2/2$ — proportional to the square of the failed member's credit limit. This is a localized impact (unlike debt-based banking where one failure can trigger cascading failures through the financial network [C:Ch.12]). The mutual credit system absorbs member failure locally, without systemic contagion, as long as the credit limit $\ell_k$ is calibrated to the member's capacity to deliver real goods and services.

### 25.5.2 Network Design Principles for Resilient Mutual Credit

Drawing on the network resilience analysis of Chapter 12 and the Lyapunov stability result:

**Principle MC-1 (Calibrated Credit Limits).** Credit limits $\ell_i$ should be calibrated to each member's productive capacity — the value of goods and services they can plausibly deliver to the network. Over-extended credit limits create systemic stress ($\Delta V$ large); under-extended limits reduce liquidity provision.

**Principle MC-2 (High Clustering for Clearing).** High clustering coefficient $\bar{C}$ [C:Ch.4] increases the density of short clearing cycles, improving clearing efficiency and reducing imbalance energy. Small-world network topology [C:Ch.12] balances clustering (short cycles) with short path lengths (efficient clearing propagation).

**Principle MC-3 (Regular Clearing Cycles).** Clearing frequency should be calibrated to the transaction volume — more frequent clearing reduces accumulated imbalances and releases credit capacity faster. Daily clearing (as in interbank systems) is appropriate for high-volume networks; weekly or monthly clearing for lower-volume community systems.

**Principle MC-4 (Redundant Clearing Paths).** High algebraic connectivity $\lambda_2(L_{G^{MC}})$ ensures multiple clearing paths between any debtor-creditor pair, reducing dependence on any single clearing route and making the system robust to edge failures [C:Ch.4, Theorem 4.2].

---

## 25.6 Commitment Pooling

### 25.6.1 The Formal Model

Commitment Pooling (CP), developed in the context of the Sarafu Network and formalized by Ruddick et al. (2021), extends mutual credit to communities that lack formal accounting infrastructure or legal enforcement capacity. It is particularly suited to low-income and informal economy settings where trust is high within communities but institutional capacity is low.

**Definition 25.5 (Commitment Pool).** A commitment pool is a tuple $(\mathcal{M}, \mathcal{CP}, \Theta_{\text{cur}}, \Theta_{\text{val}}, \Theta_{\text{lim}}, \Theta_{\text{exc}})$ where:
- $\mathcal{M}$: member community.
- $\mathcal{CP}$: the pool of committed value — a shared accounting of members' productive commitments.
- $\Theta_{\text{cur}}$: **Curation function** — the governance process by which members are admitted and commitments are accepted as valid.
- $\Theta_{\text{val}}$: **Valuation function** — the process by which commitments are given a common unit of account, enabling exchange across different types of commitment (food, labor, care, craft).
- $\Theta_{\text{lim}}$: **Limitation function** — the process by which credit limits are set, balances are monitored, and over-extension is prevented.
- $\Theta_{\text{exc}}$: **Exchange function** — the mechanism by which committed value is exchanged between members.

**The four functions.** These four functions correspond directly to the Ostrom design principles [C:Ch.14]: $\Theta_{\text{cur}}$ implements DP1 (defined boundaries — who belongs); $\Theta_{\text{val}}$ implements DP2 (congruence — rules adapted to local conditions, including what counts as valuable); $\Theta_{\text{lim}}$ implements DP4 and DP5 (monitoring and graduated sanctions); and $\Theta_{\text{exc}}$ implements DP3 (collective choice — the community governs the exchange mechanism).

**Formal commitment structure.** Each member $i$ commits to deliver value $q_{ij}$ of type $s_{ij}$ (e.g., 3 kg of maize, 2 hours of child care, 1 haircut) to any requesting member $j$, up to their credit limit $\ell_i$. The valuation function converts these heterogeneous commitments to a common unit of account:

$$v(q_{ij}, s_{ij}) = \Theta_{\text{val}}(q_{ij}, s_{ij}) \in \mathbb{R}_+$$

The pool's aggregate value is:

$$\mathcal{CP} = \sum_i \sum_{j,s} v(q_{ij}, s_{ij}) \leq \sum_i \ell_i$$

---

## 25.7 Mathematical Model: Clearing LP and Stability

We develop the complete clearing LP and stability analysis for a 50-member B2B mutual credit network — the worked example of Section 25.8. Here we present the mathematical structure abstractly.

**Setup.** Network $G^{MC} = (V, E, w, b)$ with $n = 50$ members, $|E| = 320$ directed credit edges (average degree $\bar{k} = 6.4$, Watts-Strogatz small-world with $\beta = 0.1$). Balance vector $\mathbf{b}$ after one trading period (before clearing).

**Clearing LP formulation:**

$$\min_{\mathbf{f}} \;\; \mathbf{c}^T \mathbf{f}$$

$$\text{s.t.} \quad B \mathbf{f} = -\mathbf{b} \quad \text{(clear all balances)}$$
$$\mathbf{0} \leq \mathbf{f} \leq \mathbf{w} \quad \text{(flow within capacity)}$$

where $B \in \mathbb{R}^{n \times |E|}$ is the signed incidence matrix, $\mathbf{c} \in \mathbb{R}^{|E|}_+$ is the cost vector (zero for automated digital clearing), and $\mathbf{w} \in \mathbb{R}^{|E|}_+$ is the capacity vector (existing credit obligations).

**When $\mathbf{c} = \mathbf{0}$ (free clearing):** The LP minimizes zero cost subject to clearing all balances. Any feasible flow is optimal. The network simplex algorithm finds a feasible clearing in $O(n^2)$ time for sparse networks.

**When full clearing is infeasible** (not all balances can be cleared — network not strongly connected or some obligations exceed capacity): solve the relaxed LP with partial clearing tolerance:

$$\min_{\mathbf{f}, \boldsymbol{\varepsilon}} \;\; \mathbf{1}^T \boldsymbol{\varepsilon}$$

$$\text{s.t.} \quad B \mathbf{f} = -\mathbf{b} + \boldsymbol{\varepsilon}, \quad \mathbf{0} \leq \mathbf{f} \leq \mathbf{w}, \quad \boldsymbol{\varepsilon} \geq \mathbf{0}$$

This minimizes uncleared balances (sum of residuals $\varepsilon_i$), finding the maximum feasible clearing given the network structure.

**Algorithm 25.1: Optimal Multilateral Clearing**

```python
FROM networkx IMPORT DiGraph, min_cost_flow
FROM numpy IMPORT array

FUNCTION optimal_clearing(members, obligations, balances, costs=None):
    """
    Compute optimal multilateral clearing for a mutual credit network.

    Parameters:
    - members: list of n member IDs
    - obligations: dict {(i,j): amount} — gross bilateral obligations
    - balances: dict {i: balance} — current balance per member (sum = 0)
    - costs: dict {(i,j): cost} — clearing cost per unit (default 0)

    Returns:
    - clearing_flows: dict {(i,j): flow} — clearing flows per edge
    - residual_balances: dict {i: balance} — post-clearing balances
    - clearing_efficiency: float — fraction of gross obligations cleared
    """

    # Build directed graph for NetworkX min_cost_flow
    G = DiGraph()

    # Add nodes with demand = -balance (NetworkX convention:
    # positive demand = source, negative demand = sink)
    FOR i IN members:
        G.add_node(i, demand=-balances[i])
        # demand = -(balance_i): creditors are sinks (demand<0),
        # debtors are sources (demand>0)

    # Add edges with capacity = obligation amount, weight = cost
    total_gross = 0
    FOR (i, j), amount IN obligations.items():
        IF amount > 0:
            cost_ij = costs.get((i,j), 0) IF costs ELSE 0
            G.add_edge(i, j,
                      capacity=amount,   # max flow = obligation amount
                      weight=cost_ij)    # clearing cost per unit
            total_gross += amount

    # Solve min-cost flow (raises NetworkXUnfeasible if infeasible)
    TRY:
        flow_dict = min_cost_flow(G)

        # Extract clearing flows and compute residual balances
        clearing_flows = {}
        cleared_volume = 0
        FOR i IN G.nodes():
            FOR j, flow_data IN G[i].items():
                f = flow_dict[i].get(j, 0)
                IF f > 0:
                    clearing_flows[(i,j)] = f
                    cleared_volume += f

        # Compute residual balances (should all be ~0 if fully feasible)
        residual = {i: balances[i] for i in members}
        FOR (i,j), f IN clearing_flows.items():
            residual[i] += f    # flow out reduces debit
            residual[j] -= f    # flow in reduces credit

        efficiency = 1 - sum(abs(r) for r in residual.values()) / \
                         (2 * sum(abs(b) for b in balances.values()))

        RETURN {
            'flows': clearing_flows,
            'residual_balances': residual,
            'clearing_efficiency': efficiency,
            'gross_obligations_cleared': cleared_volume,
            'total_gross': total_gross
        }

    EXCEPT NetworkXUnfeasible:
        # Fall back to partial clearing
        RETURN partial_clearing(G, members, balances)
```

**Scalability.** NetworkX's min-cost flow implementation uses the network simplex algorithm with complexity $O(VE \log V \log(V C_{\max}))$ where $C_{\max}$ is the maximum cost. For a 50-node, 320-edge network with integer capacities, clearing completes in $< 10$ milliseconds. For a 10,000-node network (large regional B2B system): $< 5$ seconds. For national scale (1 million nodes): a distributed multi-agent clearing protocol is needed — see the ★★ exercise.

---

## 25.8 Worked Example: 50-Firm B2B Mutual Credit Network

We design a mutual credit system for a 50-firm business-to-business network in a regional economy (inspired by the structure of regional business barter exchanges operating in northern Italy and Switzerland).

### 25.8.1 System Design

**Member firms:** 50 SMEs across 5 sectors: manufacturing (12), services (15), agriculture (8), retail (10), construction (5). Average annual turnover: EUR 1.2 million. Suggested credit limit: 5\% of annual turnover = EUR 60,000 per firm.

**Credit limit calibration (Principle MC-1).** Credit limits are set in proportion to each firm's average monthly invoicing within the network — calibrated after a 3-month observation period. Initial limits: EUR 30,000 (conservative half-limit during probationary period).

**Network topology (Principle MC-2).** Watts-Strogatz small-world with $n = 50$, $\bar{k} = 8$, $\beta = 0.10$: clustering coefficient $\bar{C} \approx 0.48$ (high, supporting short clearing cycles), average path length $\bar{d} \approx 3.1$ (short, supporting clearing propagation), Fiedler value $\lambda_2 \approx 1.7$ (robust connectivity).

### 25.8.2 Optimal Clearing: A Trading Period

After one trading month, the 50 firms have generated the following balance distribution:
- Mean balance: 0 (system constraint)
- Standard deviation: EUR 12,400
- Maximum debit: EUR 28,600 (firm 23, construction sector — high input purchases, lower reciprocal sales this month)
- Maximum credit: EUR 31,200 (firm 41, agricultural supplier — high sales, lower purchases)

Gross bilateral obligations: EUR 847,000 (sum of all individual transaction amounts).
Net bilateral positions after bilateral netting: EUR 412,000 (51\% reduction from gross).
Net after optimal multilateral clearing (Algorithm 25.1): EUR 124,000 (85\% reduction from gross).

**Clearing efficiency breakdown:**
- Bilateral netting eliminates: 51\% of gross
- Multilateral clearing adds: 34\% additional (3-member cycles: 18\%, 4-member: 9\%, 5+ member: 7\%)
- Total clearing rate: 85\%
- Residual balances (require settlement in national currency or carry-forward): 15\% of gross = EUR 127,000

**Clearing verification.** After applying Algorithm 25.1, all 50 firms have post-clearing balances within EUR 2,000 of zero — within the system's acceptable carry-forward tolerance. The maximum residual (EUR 1,890, firm 23) is within credit limit and is carried to next period.

### 25.8.3 Resilience Analysis: Firm Failure

Suppose firm 23 (maximum debit of EUR 28,600) fails at the start of the next period, unable to honor any commitments. The impact on the system:

**Direct loss:** EUR 28,600 (firm 23's net debit is written off, distributed proportionally across its creditors within the network).

**Per-creditor impact:** Firm 23 has 6 direct creditors in the network. Average loss per creditor: EUR 4,767 — absorbed within each creditor's credit limit buffer.

**System imbalance after failure:** $\sum_i b_i$ shifts by EUR 28,600 (system is no longer balanced after write-off). Resolution: redistribute the written-off amount across all remaining members as a proportional levy ($\approx$ EUR 583 per member) or absorb through a mutual credit guarantee fund (funded by a small fee on each transaction).

**Comparison with banking system contagion.** In a banking system, a EUR 28.6 million interbank obligation (scaled to economy-wide proportions) can trigger cascading defaults [C:Ch.12, 2008 crisis analysis]. In the mutual credit network, a EUR 28,600 member failure is absorbed locally with no systemic contagion — the localized loss property of Section 25.5.1.

---

## 25.9 Case Study: Swiss WIR and Sarafu Network

### 25.9.1 The Swiss WIR System (1934–Present)

The Swiss WIR Bank (Wirtschaftsring-Genossenschaft, "Economic Circle Cooperative") was founded in 1934 by a group of Swiss businessmen seeking to create a resilient B2B exchange system during the Great Depression. As of 2024, it has approximately 62,000 participating SME members across Switzerland, with annual WIR franc (CHW) turnover of approximately CHF 1.5 billion — roughly 0.25\% of Swiss GDP.

**The counter-cyclical mechanism.** WIR volume is empirically counter-cyclical: it expands when the Swiss economy contracts and contracts when the economy expands [Stodder, 2009; Stodder and Lietaer, 2016]. This is precisely what the formal model predicts:

- In contractions: CHF liquidity is scarce (banks tighten credit), increasing the relative attractiveness of WIR credit. Members shift trade toward WIR to preserve CHF for essential payments. $\Delta M^{\text{WIR}} > 0$.
- In expansions: CHF liquidity is abundant. Members prefer CHF (wider acceptance, no balance constraints). $\Delta M^{\text{WIR}} < 0$.

**Formal representation.** Let $x_{\text{WIR}}(t) \in [0,1]$ be the fraction of B2B trade settled in WIR versus CHF. The adoption dynamics:

$$\dot{x}_{\text{WIR}} = \gamma \left[\frac{\ell_{\text{WIR}}}{i_{\text{CHF}} \cdot \text{credit\_tightness}(t)} - x_{\text{WIR}}\right]$$

When credit tightness rises (recession): the WIR credit advantage (no interest, easy access) relative to CHF credit (interest rate $i_{\text{CHF}}$, tighter access) increases, pushing $x_{\text{WIR}}$ upward. This is the automatic stabilizer: the WIR system expands exactly when it is needed most.

**90-year stability.** The WIR system has operated continuously for 90 years, surviving the Great Depression, WWII, the oil shocks of the 1970s, the 2008 crisis, and the COVID-19 pandemic. Its stability reflects the Ostrom conditions of Chapter 14:

| DP | WIR Implementation | Score |
|:---|:---|:---|
| DP1 (Boundaries) | Formal membership with vetting | 2/2 |
| DP2 (Congruence) | WIR pricing reflects CHF prices | 2/2 |
| DP3 (Collective choice) | Cooperative governance, member votes | 2/2 |
| DP4 (Monitoring) | Bank-maintained transaction records | 2/2 |
| DP5 (Graduated sanctions) | Balance limits, fee escalation | 2/2 |
| DP6 (Conflict resolution) | Arbitration panel | 1/2 |
| DP7 (External recognition) | FINMA-regulated bank | 2/2 |
| DP8 (Nested enterprises) | Regional chapters within national bank | 2/2 |
| **Total** | | **15/16** |

The WIR's Ostrom score of 15/16 is identical to the Maine lobster fishery of Chapter 14 — and its 90-year longevity is the mutual credit analogue of that fishery's six-decade governance success.

### 25.9.2 The Sarafu Network (Kenya, 2010–Present)

The Sarafu Network is a commitment pooling system operating in low-income communities across Kenya, coordinated by Grassroots Economics Foundation. As of 2023, it involves approximately 55,000 registered users across more than 100 community currencies, with cumulative transaction volume exceeding USD 50 million equivalent.

**The commitment structure.** Each community issues its own token (e.g., "Kangemi" for the Kangemi neighborhood of Nairobi) representing commitments from community businesses and households. A haircut shop commits to redeem 10 Kangemi for one haircut; a vegetable seller commits to redeem 5 Kangemi for 1 kg of tomatoes. The token circulates within the community as a medium of exchange — backed by the collective productive capacity of its members.

**The Sarafu hub.** Sarafu (the network-level token) connects local community currencies: at coordinated swap events, community tokens can be exchanged for Sarafu at a managed rate, enabling inter-community trade. This creates the nested polycentric structure of the Cosmo-Local model [C:Ch.13]: local currencies for local trade, Sarafu for inter-community coordination.

**Poverty reduction outcomes.** Ruddick et al. (2021) report that Sarafu users in food-insecure households increased their food expenditure by approximately 15\% above the control group during the COVID-19 lockdown period (2020), when national currency liquidity was severely constrained. The mutual credit mechanism provided liquidity precisely when the formal monetary system contracted — the same counter-cyclical property as the WIR, now in a low-income informal economy context.

**Formal assessment.** In the commitment pooling framework (Definition 25.5):
- $\Theta_{\text{cur}}$ (Curation): Community governance committees vet new members and commitments.
- $\Theta_{\text{val}}$ (Valuation): Commitments are denominated in Sarafu at community-determined rates.
- $\Theta_{\text{lim}}$ (Limitation): Accounts are limited to 200 Sarafu negative balance; monitoring through mobile app.
- $\Theta_{\text{exc}}$ (Exchange): Sarafu tokens transferred via mobile phone using the Celo blockchain for settlement.

The system implements all four commitment pooling functions at low cost, leveraging mobile phone penetration (96\% in Kenya) and blockchain settlement (no bank account required) to operate in contexts where traditional mutual credit infrastructure would be infeasible.

---

## Chapter Summary

This chapter has developed the formal theory of mutual credit and clearing, demonstrating that decentralized liquidity provision through bilateral commitments and multilateral clearing can achieve the same real allocations as a monetary economy with perfect credit markets — without requiring any central money creation.

Mutual credit mechanics (Definition 25.1) create liquidity through double-entry commitments: each transaction credits one party and debits another, preserving the system balance constraint $\sum b_i = 0$ at all times (Proposition 25.1). Total liquidity provided equals the sum of credit limits, achievable with zero monetary stock.

The clearing problem (Definition 25.3) is formalized as a minimum-cost flow LP, solvable in polynomial time by the network simplex algorithm. Theorem 25.1 establishes clearing feasibility (requires strong connectivity); Corollary 25.1 connects clearing efficiency to algebraic connectivity. Theorem 25.2 (Optimal Clearing) proves equivalence with perfect credit markets; Theorem 25.3 proves multilateral dominates bilateral clearing whenever cycles of length ≥ 3 exist.

Lyapunov stability analysis (Proposition 25.2) shows that periodic optimal clearing bounds system imbalance and that member failures impose local, bounded losses — in sharp contrast to the cascading failures characteristic of scale-free financial networks [C:Ch.12].

Commitment Pooling (Definition 25.5) extends mutual credit to informal economy contexts through four governance functions (curation, valuation, limitation, exchange) that implement the Ostrom design principles at community scale.

The Swiss WIR demonstrates 90 years of mutual credit stability with a 15/16 Ostrom score and formal counter-cyclical stabilization properties — expanding exactly when CHF liquidity contracts. The Sarafu Network demonstrates that commitment pooling is achievable in low-income contexts with mobile technology, providing meaningful poverty-reduction impact during the COVID-19 crisis.

Chapter 26 develops the third alternative monetary architecture: resource-backed currencies, in which the money supply is formally linked to a physical resource, embedding ecological constraints directly into the monetary medium.

---

## Exercises

**25.1** Model a 5-node mutual credit network with initial obligations: A owes B 80, B owes C 60, C owes D 40, D owes E 30, E owes A 50, A owes C 20, D owes B 10.
(a) Compute the initial balance vector $\mathbf{b}$. Verify $\sum b_i = 0$.
(b) Formulate the clearing LP (Definition 25.3) with zero clearing costs. Write out the constraint matrix $B$ and the constraint vector $-\mathbf{b}$.
(c) Solve the LP (by inspection or linear programming). What is the optimal clearing flow? What fraction of gross obligations is cleared?
(d) Compare to bilateral netting: what fraction of gross obligations is cleared by bilateral netting alone? How much additional clearing does the multilateral approach achieve?

**25.2** The counter-cyclical WIR mechanism (Section 25.9.1):
(a) For Swiss GDP growth data 1995–2020 and WIR turnover data from the WIR Bank Annual Reports, compute the correlation between $\Delta \text{WIR}$ and $\Delta \text{GDP}$. Is it negative (counter-cyclical)?
(b) Using the formal model $\dot{x}_{\text{WIR}} = \gamma[\ell_{\text{WIR}}/(i_{\text{CHF}} \cdot \text{tightness}) - x_{\text{WIR}}]$, calibrate $\gamma$ and the tightness function using the 2008–2009 data. During the 2008 crisis, how quickly did WIR volume expand?
(c) Estimate the GDP stabilization value of the WIR system: if the WIR system had not existed during the 2008 recession, by how much would Swiss GDP have fallen? Use the multiplier relationship between WIR volume and GDP.

**25.3** Apply the stability analysis of Section 25.5 to the 50-firm worked example.
(a) Compute the Lyapunov function $V(\mathbf{b})$ before and after clearing. By how much does clearing reduce $V$?
(b) Firm 23 fails with net debit EUR 28,600. Compute the direct loss to each of its 6 creditors. Is any creditor's loss large enough to push them below their credit limit?
(c) Two firms fail simultaneously (firm 23 and firm 31, with net debit EUR 19,200). Compute the combined impact and whether the system can absorb both failures within credit limits. What mutual guarantee fund size (as a fraction of total credit limits) would be required to cover the 99th percentile of simultaneous failures?

**★ 25.4** Prove Theorem 25.2 (Optimal Clearing Theorem) in full generality.

(a) Define a monetary economy with perfect credit markets: each agent $i$ has income $y_i$, credit limit $\ell_i$, and can purchase any bundle of goods $(q_{ij})_{j \neq i}$ satisfying $\sum_j p_j q_{ij} \leq y_i + \ell_i$. Characterize the feasible set of allocations.
(b) Define a mutual credit economy: same agents, same preferences, same credit limits, but no prior money — only credit commitments. A feasible allocation requires that for every agent $i$, $\sum_j v(q_{ij}) \leq y_i + \ell_i$ (can purchase up to credit limit plus productive endowment).
(c) Show that the feasible allocation sets of the two economies are identical when perfect credit markets are assumed (no transaction costs, no rationing).
(d) Prove that complete multilateral clearing restores all agents' credit capacity to allow the next round of allocation — i.e., that repeated rounds of trade and clearing support an indefinitely sustained sequence of efficient allocations without money creation.

**★ 25.5** Prove Theorem 25.3 (multilateral dominates bilateral clearing) and derive the efficiency gain as a function of network cycle structure.

(a) Define bilateral netting formally: for each pair $(i,j)$, replace $\{w_{ij}, w_{ji}\}$ with $\{\max(0, w_{ij}-w_{ji}), 0\}$. Compute the gross obligations remaining after bilateral netting.
(b) Show that bilateral netting leaves all cycles of length $\geq 3$ uncleared. Give an example of a 3-cycle where bilateral netting clears 0\% of gross obligations.
(c) Show that multilateral clearing clears all cycles of any length. Compute the additional clearing achieved by multilateral vs. bilateral for a network with $k$ independent 3-cycles, each of value $v$.
(d) For a random directed graph with $n$ nodes and mean degree $\bar{k}$, estimate the expected number of 3-cycles and the expected additional clearing efficiency of multilateral over bilateral clearing. Express the result as a function of $n$ and $\bar{k}$.

**★★ 25.6** Design a commitment pooling mechanism for a 200-household community in a low-income urban setting (inspired by informal settlements in Nairobi or Lagos); prove incentive-compatibility; simulate outcomes under different valuation rules.

**Community specification:**
- 200 households, median income USD 80/month.
- 12 types of local goods/services: cooked food, childcare, laundry, haircut, phone charging, vegetable sales, crafts, transport, construction labor, medical herbs, tutoring, music/entertainment.
- Mobile phone penetration: 90\%. Bank account penetration: 25\%.
- Target: Design a commitment pooling system that increases food security and provides liquidity during income shocks.

(a) Specify all four commitment pooling functions ($\Theta_{\text{cur}}, \Theta_{\text{val}}, \Theta_{\text{lim}}, \Theta_{\text{exc}}$) for this community. For the valuation function, propose two alternative approaches: (i) community-voted prices; (ii) market-reference prices adjusted by a fairness committee. Compare their incentive properties.

(b) Prove that the commitment pooling mechanism is incentive-compatible: it is a weakly dominant strategy for each household to provide their committed goods/services upon request, rather than defecting (accepting credits without honoring commitments). Use the Folk Theorem framework [C:Ch.7] and the Ostrom monitoring/sanction conditions.

(c) Simulate the system for 50 households (a manageable sub-scale) over 100 trading periods using: (i) random trade patterns (each period, each household randomly purchases from 2 others); (ii) preferential trade patterns (households prefer trading with neighbors). Under each pattern, compute average balances, clearing efficiency, and fraction of households who remain within credit limits throughout.

(d) Introduce an income shock at period 50: 30 households lose 50\% of their national currency income for 10 periods (simulating a lockdown or layoff). How does the commitment pool respond? Compare food security (food purchases per household per period) under the shock in: (i) the commitment pool system; (ii) a baseline with no complementary currency.

---

*Chapter 26 develops the third alternative monetary architecture: resource-backed currencies — monetary systems in which the money supply is formally linked to a physical resource (energy, carbon, land), embedding ecological constraints directly into the monetary medium and creating automatic incentives for stewardship that purely institutional approaches cannot provide.*
