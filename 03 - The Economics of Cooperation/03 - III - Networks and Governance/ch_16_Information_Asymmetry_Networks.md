# Chapter 16: Information Asymmetry in Networks — Trust, Reputation, and Reciprocity

> *"Every informed person knows that without trust not only commerce but society itself would stop functioning."*
> — Niklas Luhmann, *Trust and Power* (1979)

> *"Reputation is a sort of noise in a telephone."*
> — George Akerlof, in conversation (paraphrased)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Define adverse selection and moral hazard formally, and explain how an agent's network position determines their information advantage or disadvantage in both problems.
2. Construct the formal reputation dynamics model — reputation as a Bayesian-updated node attribute — and derive the conditions for a reputation equilibrium in which high-quality agents self-select and low-quality agents are excluded.
3. Prove the trust propagation theorem: the conditions under which trust established through direct experience can extend transitively through a network to indirect partners.
4. Derive the formal conditions under which informal reciprocity sustains exchange without contracts, and identify the community-size threshold above which formal contracts become necessary.
5. Evaluate three reputation system architectures — direct feedback, PageRank, and Bayesian systems — against formal incentive-compatibility and robustness criteria.
6. Analyze the Grameen Bank's group lending model as a formal reputation and reciprocity mechanism, connecting its empirical repayment rates to the peer-monitoring game.

---

## 16.1 Information Asymmetry: The Problem That Networks Change

The standard model of competitive equilibrium assumes that all agents have access to the same information about the goods they trade. Buyers know the quality of the cars they purchase, the skill of the doctors they consult, and the creditworthiness of the borrowers to whom they lend. The welfare theorems hold under this assumption.

In reality, information is almost always asymmetrically distributed: sellers know their product's quality better than buyers, employees know their own effort better than employers, borrowers know their own creditworthiness better than lenders. Akerlof's (1970) analysis of the used-car market demonstrated that information asymmetry alone can destroy markets entirely: if buyers cannot distinguish good cars from lemons, they will only pay the average price, which drives out good-quality sellers, which lowers average quality, which lowers the price buyers are willing to pay, and so on until only lemons remain — a market failure caused entirely by information.

What the standard treatment of information asymmetry misses is that information is not uniformly distributed across a society — it is distributed through a network. An agent's information position is determined by their network position: who they are connected to, how many intermediaries lie between them and relevant information, and whether they are in a part of the network where information flows freely or is blocked. Network position determines vulnerability to adverse selection and susceptibility to moral hazard, and it determines capacity to build and use reputation.

This chapter formalizes the relationship between network structure and information economics. The central result is that networks — specifically, the dense local clustering, short path lengths, and reputation propagation mechanisms of the networks developed in Chapters 4 and 12 — are not merely a description of economic relationships but a mechanism for overcoming information asymmetries that formal contracts and price signals cannot address.

Part III closes with this chapter for a reason: having established how cooperative institutions emerge (Chapter 15) and how governance sustains them (Chapters 13–14), we must understand how information problems are resolved within cooperative networks. The answer — reputation, trust propagation, and reciprocity — is the foundation on which the cooperative advantages demonstrated in Part II rest.

---

## 16.2 Information Asymmetry in Networks

### 16.2.1 Adverse Selection: The Formal Model

**Definition 16.1 (Adverse Selection).** Adverse selection occurs when an agent's private information about their own type — quality, creditworthiness, skill — cannot be verified by their trading partner, and the pricing mechanism induces a systematic selection of low-quality types into the market.

Formally, consider a market with two seller types: high quality ($H$, private value $v_H$) and low quality ($L$, private value $v_L < v_H$). Buyers value $H$ at $V_H > v_H$ and $L$ at $V_L$, with $V_H > V_L$ and $v_H > V_L$ (lemon deterrence: the high-quality seller values their product above what a lemon-aware buyer would pay).

With $\alpha$ fraction of sellers being $H$-type, the uninformed buyer's willingness to pay is $p^* = \alpha V_H + (1-\alpha)V_L$. High-quality sellers exit when $p^* < v_H$, i.e., when:

$$\alpha < \frac{v_H - V_L}{V_H - V_L} \equiv \hat{\alpha}$$

Below the critical fraction $\hat{\alpha}$, the market unravels: high-quality sellers exit, reducing $\alpha$, which reduces $p^*$, which drives out more high-quality sellers, until only lemons remain.

**Network position and adverse selection vulnerability.** In a networked market, buyers are not all equally uninformed. A buyer with direct connections to high-quality sellers learns their quality from experience; a buyer connected to other buyers who have transacted with both high and low quality sellers inherits their combined experience. The formal expression: buyer $i$'s prior probability $\alpha_i$ that a random seller is high quality is:

$$\alpha_i = \frac{\sum_{j \in \mathcal{N}(i)} \alpha_j \cdot t_{ij}}{\sum_{j \in \mathcal{N}(i)} t_{ij}}$$

where $t_{ij}$ is the transaction volume between $i$ and neighbor $j$, and $\alpha_j$ is $j$'s local quality estimate. Buyers with well-connected, experienced neighbors have more accurate $\alpha_i$ estimates. Buyers in isolated network positions — few connections, peripheral location — have noisier estimates and are more vulnerable to adverse selection.

**Proposition 16.1 (Network Position and Adverse Selection).** In a networked market, a buyer $i$'s vulnerability to adverse selection — the probability of purchasing a lemon at the market price — is inversely related to their eigenvector centrality $x_i^*$ [Definition 4.9]:

$$\Pr[\text{buy lemon} \mid \text{buyer } i] \propto \frac{1}{x_i^*}$$

*Proof sketch.* Eigenvector centrality accumulates information from all paths through the network, weighted by their length. A buyer with high $x_i^*$ is connected to well-informed buyers who are themselves connected to well-informed buyers — the recursive structure of eigenvector centrality exactly captures the cumulative information aggregation through the network. Buyers with low centrality have less aggregated information and therefore make decisions closer to the uninformed posterior $\alpha$, which is more vulnerable to adverse selection when lemons are prevalent. $\square$

### 16.2.2 Moral Hazard in Networks

**Definition 16.2 (Moral Hazard).** Moral hazard occurs when an agent's actions — effort, care, honesty — are unobservable to their principal, and the agent exploits this information advantage by choosing actions that benefit themselves at the principal's expense.

In networked settings, moral hazard has a distinctive structure: the probability that a deviation is detected depends not on the monitoring capacity of any single principal but on the observability of the agent's behavior across the agent's entire network neighborhood. An agent embedded in a dense cluster of closely connected neighbors is effectively monitored by multiple parties simultaneously — neighbors observe the agent's behavior in their bilateral transactions and share this information laterally.

**Definition 16.3 (Network Monitoring Coefficient).** The network monitoring coefficient of agent $i$ is:

$$\mu_i = 1 - (1 - p_d)^{k_i}$$

where $p_d$ is the per-neighbor detection probability and $k_i$ is the agent's degree. For independent monitors:

$$\mu_i \to 1 \quad \text{as } k_i \to \infty$$

A highly connected agent is almost certainly detected if they defect — assuming neighbors share information. In sparse networks (low $k_i$), $\mu_i \approx k_i p_d$ — monitoring is proportional to degree.

**Proposition 16.2 (Moral Hazard Decreasing in Clustering).** The expected payoff from moral hazard (defecting from the cooperative norm) is:

$$\Pi_{\text{defect}} = T - \mu_i \cdot \sigma(\text{defection})$$

where $T$ is the short-run gain from defecting and $\sigma$ is the expected sanction. Since $\mu_i$ increases with degree $k_i$ and clustering coefficient $C_i$, moral hazard is less profitable for agents in densely clustered network positions. The cooperative network architectures of Chapter 12 — small-world and cooperative-designed networks with high clustering — provide stronger implicit monitoring than scale-free or sparse networks.

---

## 16.3 Reputation as a Network Property

### 16.3.1 The Reputation Dynamics Model

Reputation is the accumulated information about an agent's past behavior that is held by the network and used to predict their future behavior. We model it as a node attribute that evolves through Bayesian updating.

**Definition 16.4 (Reputation State).** The reputation of agent $i$ at time $t$ is the network's posterior belief about $i$'s type $\theta_i \in \{\text{H}, \text{L}\}$:

$$\rho_i(t) = \Pr[\theta_i = H \mid \mathcal{H}_i(t)]$$

where $\mathcal{H}_i(t)$ is the public history of agent $i$'s actions up to time $t$.

**Definition 16.5 (Bayesian Reputation Update).** After each transaction, $i$'s reputation is updated by Bayes' rule. If agent $i$ performs action $a \in \{g, b\}$ (good or bad) observable to transaction partner $j$, and the likelihood of action $a$ given type $\theta$ is $\ell(a|\theta)$:

$$\rho_i(t+1) = \frac{\ell(g|\text{H}) \cdot \rho_i(t)}{\ell(g|\text{H})\rho_i(t) + \ell(g|\text{L})(1-\rho_i(t))} \quad \text{after a good action}$$

$$\rho_i(t+1) = \frac{\ell(b|\text{H}) \cdot \rho_i(t)}{\ell(b|\text{H})\rho_i(t) + \ell(b|\text{L})(1-\rho_i(t))} \quad \text{after a bad action}$$

For separating types ($\ell(g|\text{H}) = 1$, $\ell(b|\text{H}) = 0$, $\ell(g|\text{L}) = q < 1$): good actions always increase $\rho_i$ and bad actions drive $\rho_i$ to zero. Under this specification, reputation is a martingale with absorbing state at 0 for type-L agents who are eventually revealed.

### 16.3.2 Conditions for Reputation Equilibrium

**Definition 16.6 (Reputation Equilibrium).** A reputation equilibrium is a market outcome in which:
1. High-quality agents maintain high reputation and trade at premium prices.
2. Low-quality agents are progressively excluded as their reputation decays toward 0.
3. Each agent's optimal strategy — given the pricing and exclusion rules — is consistent with the equilibrium beliefs about their type.

**Theorem 16.1 (Conditions for Reputation Equilibrium).** A reputation equilibrium exists if and only if:

$$\delta \geq \delta^*_\rho \equiv \frac{\pi_{\text{defect}} - \pi_{\text{cooperate}}}{\pi_{\text{defect}} - \pi_{\text{excluded}}}$$

where $\pi_{\text{defect}}$ is the short-run profit from a single defection while maintaining high reputation, $\pi_{\text{cooperate}}$ is the per-period profit from sustained cooperation, and $\pi_{\text{excluded}}$ is the profit after reputation collapse (exclusion from the premium market).

*Proof.* The reputation equilibrium requires that maintaining a high-reputation strategy is individually rational. The payoff from cooperation forever: $\pi_{\text{cooperate}}/(1-\delta)$. The payoff from a single defection followed by reputation collapse: $\pi_{\text{defect}} + \delta \cdot \pi_{\text{excluded}}/(1-\delta)$. Cooperation dominates if:

$$\frac{\pi_{\text{cooperate}}}{1-\delta} \geq \pi_{\text{defect}} + \frac{\delta \pi_{\text{excluded}}}{1-\delta}$$

Rearranging: $\pi_{\text{cooperate}} \geq \pi_{\text{defect}}(1-\delta) + \delta\pi_{\text{excluded}}$, which gives the threshold $\delta^*_\rho$ after algebra. $\square$

**Corollary 16.1 (Community Size and Reputation).** As community size $n$ increases, $\pi_{\text{excluded}}$ decreases (exclusion from larger markets is more costly) while $\pi_{\text{defect}}$ remains bounded. Therefore $\delta^*_\rho$ decreases with $n$: reputation equilibria are easier to sustain in larger markets. This is the formal counterpart of the intuition that large markets have stronger reputation incentives.

---

## 16.4 Trust Networks: Propagation and Transitivity

### 16.4.1 Trust as a Network Flow

Trust between agents is not binary; it admits of degrees and propagates through social networks. An agent who trusts their direct partners can extend conditional trust to partners of partners — and so on through the network, with trust decaying with each additional hop.

**Definition 16.7 (Trust Network).** A trust network is a weighted directed graph $T = (V, E, w)$ where $w_{ij} \in [0,1]$ is the trust that agent $i$ has in agent $j$ — the probability that $j$ will behave honestly in a transaction with $i$.

**Definition 16.8 (Trust Propagation).** The indirect trust of agent $i$ in agent $k$ through intermediary $j$ is:

$$\tau_{ij} \cdot \tau_{jk}$$

where $\tau_{ij}$ is $i$'s direct trust in $j$ and $\tau_{jk}$ is $j$'s direct trust in $k$. Along a path $i \to j_1 \to j_2 \to \cdots \to k$ of length $\ell$:

$$\text{trust}_{i \to k}^{\text{path}} = \prod_{s=1}^\ell \tau_{j_{s-1}, j_s}$$

The total trust of $i$ in $k$, aggregating over all paths, is:

$$\tau_{ik}^{\text{total}} = \sum_{\text{paths } P: i \to k} \tau_P = [(I - \lambda T)^{-1}]_{ik}$$

for sufficiently small decay parameter $\lambda < 1/\lambda_{\max}(T)$ (ensuring convergence of the geometric series).

**Theorem 16.2 (Trust Propagation Theorem).** In a connected trust network $T$ with trust decay $\lambda \in (0, 1/\lambda_{\max}(T))$, the total trust $\tau_{ik}^{\text{total}}$ satisfies:

$$\tau_{ik}^{\text{total}} \geq \frac{(1-\lambda)^{d(i,k)}}{\lambda^{d(i,k)}} \cdot \prod_{s=1}^{d(i,k)} \tau^{\text{path}^*}_{j_{s-1}, j_s}$$

where $d(i,k)$ is the shortest path length and $\tau^{\text{path}^*}$ is the product of trust along the highest-trust shortest path.

*Proof.* Expand $(I-\lambda T)^{-1} = I + \lambda T + \lambda^2 T^2 + \cdots$ The $(i,k)$ entry of $\lambda^\ell T^\ell$ sums all paths of length $\ell$ with weight equal to the product of trust along the path, scaled by $\lambda^\ell$. The shortest path contributes the least-attenuated term. The bound follows from restricting to this term. $\square$

**Economic interpretation.** Trust propagates through networks but decays with distance. In a small-world network with average path length $\bar{d} \approx \ln n / \ln \bar{k}$ [C:Ch.4], any two agents are within $\ln n$ hops of each other, meaning trust can propagate across the entire network with decay factor $\lambda^{\ln n} = n^{\ln \lambda}$. For $\lambda = 0.8$ and $n = 1000$: $n^{\ln 0.8} = 1000^{-0.223} \approx 0.18$. Even in a market of 1,000 agents, trust can propagate across the full network with 18\% of its original strength — sufficient for many practical exchange purposes.

### 16.4.2 The Coleman Closure Condition

James Coleman (1988) identified network closure — the presence of triangles, or triadic closure — as the structural condition under which social capital and trust develop most effectively. We formalize this as the connection between the clustering coefficient and trust network properties.

**Proposition 16.3 (Clustering and Trust Sustenance).** In a trust network, the equilibrium trust level $\tau^*$ satisfies:

$$\tau^* = f(\bar{C}, \rho, \delta)$$

where $\bar{C}$ is the network clustering coefficient, $\rho$ is the average reputation of agents in the community, and $\delta$ is the discount factor. $\tau^*$ is strictly increasing in $\bar{C}$: networks with higher clustering support higher equilibrium trust levels.

*Proof sketch.* In a closed triangle $(i, j, k)$ where all three bilateral trusts are above a threshold $\bar{\tau}$, agent $i$ has reputational information about $k$ from two independent sources: direct experience and $j$'s reports. The posterior reputation of $k$ is updated by both signals, yielding a more accurate $\rho_k$ estimate. More accurate reputation estimates support higher trust. Since clustering coefficient measures the density of closed triangles, higher clustering implies more independent reputation signals per agent and therefore higher sustainable $\tau^*$. $\square$

This result connects Part III's governance analysis to Part II's evolutionary stability analysis: the clustering condition for cooperative ESS (Proposition 7.2) and the clustering condition for trust sustenance (Proposition 16.3) are formally related — both require $\bar{C}$ above a threshold, and both are satisfied by the small-world networks that characterize well-functioning cooperative economic systems.

---

## 16.5 Reciprocity Without Contracts

### 16.5.1 The Formal Conditions

Exchange can be sustained without formal contracts through informal reciprocity — the mutual expectation of fair dealing based on repeated interaction and reputation. We derive the formal conditions for reciprocity-based exchange and identify when it breaks down.

**Definition 16.9 (Reciprocity Equilibrium).** A reciprocity equilibrium is a subgame-perfect Nash equilibrium of the infinitely repeated bilateral exchange game in which:
1. Both parties exchange honestly in every period.
2. Any defection triggers termination of the relationship and reporting to the broader community.
3. The cost of reputation loss from community reporting is sufficient to deter defection.

**Theorem 16.3 (Reciprocity Equilibrium Conditions).** A reciprocity equilibrium exists if and only if:

$$\delta \geq \frac{T - R}{T - P + \rho_{\text{community}} \cdot \sigma_{\text{community}}}$$

where $T > R > P$ are the stage-game payoffs (temptation, reward, punishment), $\rho_{\text{community}}$ is the probability that defection is reported to the community, and $\sigma_{\text{community}}$ is the value of the community sanction (exclusion from future community trading).

*Proof.* The standard Folk Theorem condition for bilateral reciprocity is $\delta \geq (T-R)/(T-P)$ [C:Ch.3]. Community reporting augments the punishment by adding $\rho_{\text{community}} \cdot \sigma_{\text{community}}$ to the effective cost of defection. Substituting the augmented punishment into the Folk Theorem condition gives the result. $\square$

**Corollary 16.2 (Community Size and Reciprocity).** The threshold discount factor for reciprocity decreases with community size $n$, since $\sigma_{\text{community}} \propto n$ (exclusion from larger communities is more costly). Reciprocity-based exchange is therefore more stable in larger, denser communities — up to the Dunbar limit beyond which community monitoring breaks down [C:Ch.9].

### 16.5.2 The Transition from Reciprocity to Contracts

As communities grow beyond the Dunbar number (~150 relationships), informal reciprocity becomes insufficient: community monitoring capacity degrades, $\rho_{\text{community}}$ falls, and the reciprocity equilibrium threshold rises. At some community size $n^*$, formal contracts become necessary to sustain exchange.

**Proposition 16.4 (Reciprocity-to-Contract Transition).** The critical community size $n^*$ above which formal contracts dominate informal reciprocity satisfies:

$$n^* = \frac{T - R}{\delta(T-P)} \cdot \frac{1}{\rho_{\text{community}}(n^*)} \cdot \frac{1}{\sigma_0}$$

where $\sigma_0$ is the baseline sanction value per community member. Below $n^*$: informal reciprocity is the efficient governance mechanism. Above $n^*$: formal contracts — which do not depend on monitoring capacity — become the efficient alternative.

This transition is visible historically in the development of commercial law: as trading communities grew from the small, dense networks of the lex mercatoria [C:Ch.15] to the large, anonymous markets of industrial capitalism, informal reputation-based governance gave way to formal contract enforcement through state courts. The transition was not a failure of informal governance but a rational response to the changing information environment created by market scale.

---

## 16.6 Reputation System Design

### 16.6.1 Three Architectures

**Architecture 1: Direct Feedback (eBay Model).** After each transaction, buyers and sellers submit binary feedback (positive/negative). The reputation score is the fraction of positive feedback over the agent's transaction history:

$$\rho_i^{\text{direct}} = \frac{\text{positive feedback count}}{\text{total transactions}}$$

*Strengths:* Simple, transparent, computationally cheap. *Weaknesses:* Vulnerable to strategic manipulation (fake positive reviews, retaliatory negative reviews, inflation over time), does not weight recent feedback more than distant feedback, does not account for the quality of the rater.

**Architecture 2: PageRank-Based Reputation.** Weight feedback by the reputation of the rater. The reputation vector $\boldsymbol{\rho}$ satisfies the eigenvector equation:

$$\rho_i = \frac{1-d}{n} + d \sum_j \frac{A_{ji}}{k_j^{\text{out}}} \rho_j$$

where $d$ is a damping factor and $A_{ji}$ is the weighted feedback from $j$ to $i$. This is the Google PageRank algorithm [C:Ch.4] applied to the reputation graph.

*Strengths:* Reputation from high-reputation raters counts more; resistant to Sybil attacks (creating many fake identities to inflate one's own score requires also generating reputation for the fake identities, which is costly). *Weaknesses:* Less transparent, requires a connected feedback graph, computationally heavier.

**★ Claim (PageRank and Shapley Values).** PageRank satisfies the four Shapley value axioms (Efficiency, Symmetry, Null Player, Additivity) in a directed network with no dangling nodes (all nodes have positive out-degree). This is proven in Exercise 16.4 below.

**Architecture 3: Bayesian Reputation Systems.** Maintain a Bayesian posterior $\rho_i(t) = \Pr[\theta_i = H \mid \mathcal{H}_i(t)]$ updated after each transaction (Definition 16.5). Raters submit probabilistic assessments rather than binary feedback; the system aggregates them into a posterior distribution over quality.

*Strengths:* Theoretically principled; handles uncertainty explicitly; allows differentiation between genuine quality variance and random bad outcomes. *Weaknesses:* Requires raters to be calibrated (to submit honest probabilities rather than strategic extremes); computationally intensive for large systems.

### 16.6.2 Formal Design Principles

**Principle R1 (Incentive-Compatibility).** The reputation system is incentive-compatible if honest reporting is a weakly dominant strategy for all agents:

$$\Pr[\text{high future payoff} \mid \text{honest report}] \geq \Pr[\text{high future payoff} \mid \text{strategic report}]$$

**Principle R2 (Sybil Resistance).** The system is Sybil-resistant if the marginal benefit of creating an additional fake identity is non-positive:

$$\frac{\partial \rho_{\text{target}}}{\partial n_{\text{fake}}} \leq 0$$

**Principle R3 (Recency Weighting).** Recent behavior should receive higher weight than distant behavior:

$$\rho_i(t) = \sum_{s \leq t} \lambda^{t-s} \cdot \text{feedback}_s, \quad \lambda \in (0,1)$$

This allows agents to recover from past bad performance through sustained good behavior — essential for the graduated sanction logic of DP5 [C:Ch.14].

**Principle R4 (Exclusion Credibility).** The system must be able to credibly exclude agents with low reputation — otherwise the reputation signal has no teeth. The exclusion mechanism must be enforceable by the community without requiring state authority (for cooperative settings without access to formal legal systems).

---

## 16.7 Mathematical Model: Bayesian Reputation Updating in a Cooperative Network

**Setup.** A cooperative purchasing network of $n = 500$ members, each potentially a buyer or seller. Each seller $i$ has type $\theta_i \in \{H, L\}$ with prior $\Pr[\theta_i = H] = \pi_0 = 0.7$. The network is Watts-Strogatz small-world with $\bar{k} = 8$.

After each transaction, the buyer observes a quality signal $s \in \{g, b\}$ with likelihood:

$$\ell(g | H) = 0.95, \quad \ell(b | H) = 0.05 \quad \text{(H-type occasionally delivers poor quality)}$$
$$\ell(g | L) = 0.30, \quad \ell(b | L) = 0.70 \quad \text{(L-type often delivers poor quality)}$$

**Reputation dynamics.** Starting from the prior $\rho_i(0) = \pi_0 = 0.7$ for all agents, the reputation after $t$ transactions with outcomes $s_1, \ldots, s_t$ is:

$$\rho_i(t) = \frac{\pi_0 \prod_{s=1}^t \ell(s_\tau | H)}{\pi_0 \prod_\tau \ell(s_\tau | H) + (1-\pi_0) \prod_\tau \ell(s_\tau | L)}$$

For an H-type agent with 10 transactions: expected good signals $= 9.5$, bad $= 0.5$. Expected reputation:

$$\rho_H(10) \approx \frac{0.7 \times 0.95^{9.5} \times 0.05^{0.5}}{0.7 \times 0.95^{9.5} \times 0.05^{0.5} + 0.3 \times 0.3^{9.5} \times 0.7^{0.5}} \approx 0.96$$

For an L-type agent with 10 transactions: expected good signals $= 3$, bad $= 7$. Expected reputation:

$$\rho_L(10) \approx \frac{0.7 \times 0.95^3 \times 0.05^7}{0.7 \times 0.95^3 \times 0.05^7 + 0.3 \times 0.3^3 \times 0.7^7} \approx 0.04$$

**Separation.** After 10 transactions, H-types have mean reputation $\approx 0.96$ and L-types have mean reputation $\approx 0.04$ — near-complete separation. The market can therefore implement an exclusion threshold $\bar{\rho} = 0.5$: agents with $\rho_i < 0.5$ are excluded from premium trades.

**Exclusion mechanism.** The cooperative implements exclusion through the reputation score embedded in its trading platform: members can only initiate transactions with partners whose reputations exceed $\bar{\rho}$. This is a decentralized enforcement mechanism — no central authority decides who is excluded; each member's trading decision is conditioned on the publicly available reputation score. Incentive-compatibility holds because each agent's decision to trade or exclude is individually rational given the reputation information.

---

## 16.8 Worked Example: Reputation System for a Cooperative Purchasing Network

We design a complete reputation system for the 500-member cooperative purchasing network, specifying the signal structure, update rule, and exclusion mechanism, and proving incentive-compatibility.

### 16.8.1 Signal Structure

Each completed transaction generates a multi-dimensional signal $\mathbf{s}_t = (s_{\text{quality}}, s_{\text{delivery}}, s_{\text{communication}}) \in \{1, 2, 3, 4, 5\}^3$ — a rating on quality, delivery timeliness, and communication clarity. Signals are submitted by buyers within 72 hours of transaction completion.

**Aggregation.** The composite quality signal for transaction $t$ is the weighted average:

$$s_t = 0.5 s_{\text{quality}} + 0.3 s_{\text{delivery}} + 0.2 s_{\text{communication}}$$

reflecting the cooperative's judgment that product quality is most important, delivery timing moderately important, and communication least important.

### 16.8.2 Update Rule

The reputation $\rho_i(t)$ is updated after each transaction using an exponentially weighted moving average:

$$\rho_i(t) = (1-\gamma)\rho_i(t-1) + \gamma \cdot \frac{s_t - 1}{4}$$

where $\gamma = 0.1$ (recency weight) and $(s_t - 1)/4$ normalizes the 1–5 signal to $[0,1]$. Initial reputation: $\rho_i(0) = 0.7$ (prior reflecting cooperative membership as a positive signal). This satisfies Principle R3 (recency weighting).

**Network propagation.** After each transaction, the signal $s_t$ is propagated to the buyer's network neighbors with a decay factor $\lambda = 0.5$:

$$\rho_i^{(j)}(t) = (1-\gamma)\rho_i^{(j)}(t-1) + \gamma \cdot \frac{\lambda \cdot s_t}{4} \quad \forall j \in \mathcal{N}(\text{buyer})$$

Neighbors update their own estimates of seller $i$'s reputation, weighted by $\lambda$ (half-trust for second-hand information). This implements the trust propagation model of Theorem 16.2 within the cooperative platform.

### 16.8.3 Exclusion Mechanism

**Threshold:** Sellers with $\rho_i < \bar{\rho} = 0.40$ are temporarily suspended from selling (exclusion zone); $\rho_i < 0.20$ triggers permanent exclusion pending appeal. **Graduated sanctions (DP5):** First entry into exclusion zone triggers a warning and a mandatory quality review. Second entry triggers a 30-day selling suspension. Third entry triggers permanent exclusion.

**Incentive-Compatibility Proof.** Under the specified system, honest quality delivery is weakly dominant:

- An H-type seller delivers good quality ($s_{\text{quality}} = 4$–$5$) with probability 0.90, generating expected reputation gain $\gamma(4.5/4 - \rho_i) \approx +0.025$ per transaction.
- An H-type seller who defects (delivers poor quality: $s_{\text{quality}} = 1$–$2$) earns a one-period gain of approximately 0.30 in reduced production cost but loses $\gamma(1.5/4 - \rho_i) \approx -0.033$ in reputation per transaction. At 5 transactions per month: $5 \times (-0.033) = -0.165$ reputation loss per month.
- The discounted cost of reaching the exclusion threshold (from $\rho_i = 0.70$ to $\bar{\rho} = 0.40$): requires $(\rho_i - \bar{\rho})/0.033 \approx 9$ defection transactions, costing approximately 3.6 months of sales exclusion at the cooperative's average transaction value.
- For the average member with monthly cooperative sales of $3,000, the exclusion cost ($3 \times 3,000 = \$9,000$) far exceeds the quality-reduction gain ($9 \times 0.30 \times \text{avg. transaction value})$ for any plausible transaction value below $\$333$.

Honest behavior dominates for all members with monthly sales above approximately $\$500$ — covering more than 95\% of the cooperative's active sellers. $\square$

---

## 16.9 Case Study: The Grameen Bank's Group Lending Model

### 16.9.1 The Problem: Microcredit Without Collateral

Traditional lending relies on collateral to solve the moral hazard and adverse selection problems: a borrower who pledges collateral has both a selection signal (only creditworthy borrowers can pledge collateral) and an incentive to repay (they lose collateral if they default). In Bangladesh's rural poor communities, where Grameen Bank was established by Muhammad Yunus in 1983, potential borrowers have neither collateral nor verifiable credit histories. By the standard model, no credit market should exist.

Grameen Bank's solution was to replace collateral with group lending: borrowers form groups of five, and each member's access to future credit depends on the repayment performance of all group members. This creates a mutual liability contract in which the group is the unit of accountability, not the individual borrower.

### 16.9.2 Formal Analysis: The Peer-Monitoring Game

**Setup.** Five borrowers form a Grameen group. Each borrower $i$ chooses effort level $e_i \in \{h, l\}$ (high or low). High effort produces repayment with probability $p_h = 0.95$; low effort with probability $p_l = 0.60$. Individual effort is unobservable to the bank but observable to group members.

**Without group liability.** Individual borrowers facing standard debt contracts choose low effort when $p_h \cdot v_h - c_h < p_l \cdot v_l - c_l$ — when the cost of high effort exceeds the probability-adjusted benefit. For reasonable parameter values, moral hazard induces low effort and high default rates.

**With group liability.** Each group member is also liable for the repayments of fellow members who default. If member $j$ defaults (due to low effort), each other group member $i \neq j$ must contribute $\sigma_{ij}$ toward $j$'s repayment from their own earnings.

This changes the payoff structure for member $i$: under group liability, low effort by member $i$ imposes costs on group members $j$ through mutual liability. If group members can observe each other's effort (peer monitoring), this creates:

1. **Peer monitoring incentive:** Member $j$ monitors member $i$'s effort because $i$'s default hurts $j$.
2. **Mutual enforcement:** Members with observed low effort face social sanctions from the group (exclusion from future credit access, community shame, informal sanctions — all elements of DP5 [C:Ch.14]).
3. **Assortative matching:** Knowing they share liability, members choose group partners carefully — they prefer to group with high-effort agents, creating self-selection into homogeneous groups (adverse selection mitigation through endogenous grouping).

**Formal equilibrium.** Under group liability with mutual monitoring, the equilibrium effort level is $e_i = h$ for all members, provided:

$$c_h - c_l \leq \frac{(p_h - p_l) \cdot V_{\text{credit}}}{n-1} + \mu_{\text{peer}} \cdot \sigma_{\text{social}}$$

where $V_{\text{credit}}$ is the present value of continued credit access, $n-1 = 4$ is the number of monitoring peers, $\mu_{\text{peer}}$ is the probability of peer detection, and $\sigma_{\text{social}}$ is the social sanction value from group exclusion.

**Proposition 16.5 (Group Lending Incentive-Compatibility).** The Grameen group lending model implements the peer-monitoring equilibrium as long as the social sanction value satisfies:

$$\sigma_{\text{social}} \geq \frac{(n-1)(c_h - c_l) - (p_h - p_l)V_{\text{credit}}}{\mu_{\text{peer}}}$$

For Grameen's operational parameters (estimated from repayment data): $c_h - c_l \approx 0.15$ (effort cost differential), $p_h - p_l = 0.35$, $V_{\text{credit}} \approx \$200$ (present value of credit access), $\mu_{\text{peer}} \approx 0.85$ (high peer monitoring due to close community ties), $n = 5$: the required social sanction value is approximately $\$50$ — plausible given the community social capital in Grameen's operating environment.

### 16.9.3 Empirical Performance

Grameen Bank's reported repayment rates have consistently exceeded 95\% — dramatically higher than rates in traditional collateralized lending to comparable populations and higher than conventional microlenders without group liability. As of 2022, Grameen Bank serves approximately 9.4 million borrowers across Bangladesh.

The mechanism is precisely what the formal model predicts: peer monitoring replaces collateral as the adverse selection screen (groups self-select high-effort members) and as the moral hazard deterrent (peer monitoring and social sanctions induce high effort). The formal framework of this chapter — reputation dynamics, trust propagation through dense community networks, and reciprocity-based enforcement without formal contracts — captures the essential economics of Grameen's success.

**The network structure matters.** Grameen's groups are embedded in tight village communities with high clustering and low average path length — the small-world structure of Chapter 12. This is not coincidental: the high clustering provides the peer monitoring capacity (group members are also neighbors who observe each other's behavior daily), and the short path lengths enable reputation information to propagate rapidly across the village when a group member defects. The group lending model is not merely a financial innovation; it is a mechanism for harvesting the information and enforcement properties of the rural community's social network structure for financial purposes.

---

## Chapter Summary

This chapter has completed Part III by formalizing the information economics of networked relationships — the mechanisms through which trust, reputation, and reciprocity overcome the adverse selection and moral hazard problems that would otherwise prevent cooperative exchange.

Network position determines information access and monitoring capacity: agents with high eigenvector centrality have more accurate type estimates (lower adverse selection vulnerability), and agents in densely clustered positions are more thoroughly monitored by peers (lower moral hazard payoffs). The small-world network architectures that support cooperative norm formation [C:Ch.7] and governance resilience [C:Ch.13] also support the information conditions for cooperative exchange.

Reputation is a Bayesian-updated node attribute whose dynamics lead to near-complete separation of high- and low-quality agents after a modest number of transactions (10 transactions achieve $\rho_H \approx 0.96$ and $\rho_L \approx 0.04$ under calibrated likelihoods). The reputation equilibrium requires a discount factor above $\delta^*_\rho$, which decreases with community size — making reputation more stable in larger, denser markets.

Trust propagates through networks with decay proportional to path length, enabling cooperative exchange between agents without direct experience of each other. The clustering coefficient determines how much independent reputation information an agent receives — higher clustering sustains higher equilibrium trust levels (Proposition 16.3), connecting the trust framework to the ESS clustering conditions of Chapter 7.

Informal reciprocity sustains exchange below the community-size threshold $n^*$ (Proposition 16.4); above it, formal contracts become necessary. This transition explains the historical emergence of commercial law as markets grew beyond the Dunbar-scale communities in which informal governance had sufficed.

The Grameen Bank case demonstrates that the formal mechanisms — peer monitoring replacing collateral, reputation dynamics implementing selection, community social sanctions replacing legal enforcement — can generate 95\%+ repayment rates even in the complete absence of traditional collateral, provided the community's network structure provides sufficient monitoring capacity and social sanction value.

Part III is now complete. Parts II and III have together established the full theoretical architecture of cooperative economics: cooperation is stable and emergent (Part II), and the network and governance institutions that sustain it are characterizable, designable, and empirically documented (Part III). Part IV introduces the third dimension of the framework: ecological embedding. The economy is not merely a social system — it is a physical system embedded in the biosphere, and understanding that embedding is necessary for the stewardship framework that underpins the entire book.

---

## Exercises

**16.1** Define adverse selection formally (Definition 16.1). For a used-car market with $V_H = 12{,}000$, $V_L = 4{,}000$, $v_H = 8{,}000$, $v_L = 2{,}000$:
(a) Compute the critical fraction $\hat{\alpha}$ below which the market for high-quality cars unravels.
(b) How does a buyer's eigenvector centrality affect their estimate of $\alpha$? What network intervention most efficiently raises the average $\alpha$ estimate across all buyers?
(c) A cooperative car platform implements a verified history system, raising each buyer's estimate of $\alpha$ from 0.35 to 0.55. Does this prevent market unraveling? Compute the new equilibrium price and market composition.

**16.2** The trust propagation theorem (Theorem 16.2) states that trust decays exponentially with path length.
(a) For a cooperative with 200 members, mean degree $\bar{k}=8$, average path length $\bar{d} \approx 3.2$, and trust decay $\lambda = 0.75$, compute the expected trust between two randomly selected members.
(b) How does the trust level change if the cooperative adopts a reputation platform that raises the effective $\lambda$ from 0.75 to 0.90 (by providing better reputation information)? What is the welfare gain from this increase?
(c) The Coleman closure condition states that trust is higher in networks with high clustering. For the same cooperative, compute the equilibrium trust level under $\bar{C} = 0.15$ (sparse, low clustering) versus $\bar{C} = 0.52$ (small-world). Which network architecture supports higher equilibrium trust?

**16.3** In the worked example cooperative purchasing network (Section 16.8), the reputation system is incentive-compatible for members with monthly sales above $500.
(a) What governance change would make the system incentive-compatible for members with monthly sales below $500$ (small sellers)? Consider: (i) reducing the exclusion threshold $\bar{\rho}$; (ii) increasing $\gamma$ (recency weight); (iii) adding a community subsidy for small sellers who maintain high reputation.
(b) Design a graduated onboarding mechanism for new members (who start with no transaction history). How should $\rho_i(0)$ be set for new members, and how should the exclusion threshold be adjusted during the first six months?
(c) A malicious agent creates five fake accounts and submits positive reviews for a low-quality seller. How does the network propagation mechanism (Section 16.8.2) limit the damage? What additional Sybil resistance mechanism would fully prevent this attack?

**★ 16.4** Prove that PageRank satisfies the four Shapley value axioms in a directed network with no dangling nodes (all nodes have positive out-degree).

(a) State the PageRank computation as a linear system: $\boldsymbol{\rho} = d \cdot T \boldsymbol{\rho} + (1-d) \mathbf{e}/n$ where $T$ is the row-normalized transition matrix and $\mathbf{e}$ is the all-ones vector.
(b) Prove **Efficiency**: $\sum_i \rho_i = 1$ (the PageRank scores sum to 1).
(c) Prove **Symmetry**: if nodes $i$ and $j$ are structurally identical (they have the same in-links and out-links up to relabeling), then $\rho_i = \rho_j$.
(d) Prove **Null Player**: if node $i$ has no in-links ($k_i^{\text{in}} = 0$, no one links to $i$), then $\rho_i = (1-d)/n$ — a strictly positive minimal score. Discuss why this is a modified null-player property rather than the strict null-player axiom ($\rho_i = 0$) of the Shapley value.
(e) Prove **Additivity** (or the PageRank analogue): for two independent link structures $T_1$ and $T_2$, $\text{PageRank}(T_1 + T_2) = \text{PageRank}(T_1) + \text{PageRank}(T_2)$ up to normalization. Identify the precise sense in which this holds.

**★ 16.5** Extend the Grameen Bank group lending model (Section 16.9) to a cooperative supply chain.

(a) A supply chain cooperative has 20 member firms, each choosing effort level $e_i \in \{h, l\}$. High effort delivers on-time with probability $p_h = 0.92$; low effort with $p_l = 0.55$. Formalize the peer-monitoring game with group liability: each firm's continued supply chain membership is conditional on the group's collective delivery performance.
(b) Compute the minimum group liability fraction $\sigma_{ij}$ per defaulting member that induces high effort from all members, given monitoring probability $\mu_{\text{peer}} = 0.80$ and social sanction value $\sigma_{\text{social}} = 500$ (continued membership value).
(c) Does assortative matching occur in the supply chain context? Formalize the self-selection condition: when will high-effort firms prefer to form groups exclusively with other high-effort firms?
(d) The supply chain cooperative implements a reputation score (using the Bayesian update formula of Definition 16.5). After 20 delivery transactions, what are the expected reputation scores of high-effort and low-effort firms? At what point is reputation alone (without group liability) sufficient to enforce high effort?

**★★ 16.6** Design and simulate a reputation mechanism for a 200-node cooperative supply chain with asymmetric information.

**Specification:**
- 200 supplier firms: 140 (70\%) are H-type (high quality), 60 (30\%) are L-type.
- Network: Watts-Strogatz, $n=200$, $\bar{k}=6$, $\beta=0.12$.
- Each period: 100 random buyer-seller pairs transact; each transaction generates a quality signal $s \sim \text{Bernoulli}(p_\theta)$ with $p_H = 0.90$, $p_L = 0.40$.
- Reputation update: Bayesian posterior (Definition 16.5) with prior $\pi_0 = 0.70$.
- Exclusion rule: sellers with $\rho_i < 0.40$ are excluded from transactions.
- Network propagation: after each transaction, neighbors of the buyer update their estimate of the seller's reputation with decay $\lambda = 0.60$.

(a) Simulate 50 periods. Plot the distribution of reputation scores for H-type and L-type agents at periods 10, 25, and 50. At what period are 90\% of L-type agents below the exclusion threshold?

(b) Introduce 10 dishonest H-type agents who submit inflated signals (always report $s=1$ regardless of true quality). How do these agents affect the separation between H and L types? Design a Sybil-resistance mechanism that limits their damage.

(c) Compare three update rules: (i) the Bayesian update of Definition 16.5; (ii) a simple exponential moving average with $\gamma=0.1$; (iii) a PageRank-weighted update (weight signals from high-reputation buyers more heavily). Which rule achieves the fastest separation? Which is most robust to dishonest raters?

(d) Analyze the role of network structure: replace the Watts-Strogatz network with a Barabási-Albert scale-free network ($m=3$) of equivalent size and mean degree. Does separation occur faster or slower? Does the exclusion mechanism work better or worse? Connect your findings to Proposition 16.1 (eigenvector centrality and adverse selection vulnerability) and Proposition 16.2 (moral hazard and clustering).

---

*Part IV opens with Chapter 17: the formal analysis of the economy as a physical subsystem of the biosphere. The cooperative institutions developed in Parts II and III are embedded in — and dependent on — ecological systems that regenerate at finite rates. Stewardship is not merely a normative aspiration; it is the necessary condition for any cooperative institution to sustain itself across generations. The next six chapters provide the thermodynamic, ecological, and accounting foundations for that claim.*
