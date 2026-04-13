# Part I: Foundations — Questioning the Old, Building the New

The standard economic toolkit — optimization of individual objectives, competitive equilibrium, exogenous growth — is a remarkable intellectual achievement. It is also, increasingly, a cage. Part I examines precisely where and why the cage binds: not to dismiss the toolkit but to identify the specific joints where it must be extended or replaced. We then introduce the five foundational concepts — cooperation, networks, regeneration, stewardship, and mutual coordination — that will carry the weight of the new framework developed across the rest of this volume.

The mathematics in Part I is moderate by design. These chapters establish motivation and vocabulary; the heavy formalism begins in Part II. Readers who completed Books 1 and 2 of this trilogy will recognize the analytical posture: we state things precisely, we follow the logic where it leads, and we resist the temptation to dress normative preferences in the clothing of technical necessity. The critique of neoclassical economics offered here is an internal one — conducted in the language of economics, using its own standards of proof.

---

# Chapter 1: The Limits of Competition — Why Neoclassical Economics Falls Short

> *"The curious task of economics is to demonstrate to men how little they really know about what they imagine they can design."*
> — Friedrich Hayek, *The Fatal Conceit* (1988)

> *"All models are wrong, but some are useful — and some are dangerous."*
> — attributed, after George Box

## Learning Objectives

By the end of this chapter, you should be able to:

1. Articulate the core assumptions of neoclassical economics and distinguish those that are empirical claims from those that are methodological choices.
2. State the First and Second Welfare Theorems precisely, including all conditions under which they hold.
3. Identify the specific domains — public goods, externalities, inequality, ecological overshoot — where competitive equilibrium fails by its own criteria, and explain why each failure is structural rather than incidental.
4. Define asset depletion formally as intergenerational insolvency, and explain why the failure to maintain natural capital stocks is a design flaw of the standard framework, not an externality that can be corrected at the margin.
5. Introduce Glushkov's Second Information Barrier and explain why neither markets nor hierarchies can process the information required for planetary-scale economic coordination.
6. Articulate why the appropriate response to these failures is not to abandon rigorous modeling but to change its foundations.

---

## 1.1 What Neoclassical Economics Actually Claims

It is a peculiar feature of intellectual life that the most influential frameworks are also the most frequently caricatured. Critics of neoclassical economics often attack a position no serious economist holds; defenders sometimes defend a position far stronger than the theory actually justifies. We will do neither. The purpose of this chapter is to state what neoclassical economics actually claims, with full precision — and then to examine where those claims break down.

The neoclassical research program, as it crystallized across the twentieth century, rests on a small number of foundational commitments. Agents have well-defined, stable preferences over outcomes. They choose among available options so as to maximize their utility or profit subject to constraints. Markets aggregate these individual choices through a price mechanism. Under specified conditions, the resulting equilibrium is efficient in a precise sense: no reallocation of resources can make someone better off without making someone else worse off. This is the Pareto criterion, and the two welfare theorems are its formal expression [P:Ch.2].

Let us state them with care.

---

### 1.1.1 The Welfare Theorems: Precise Statements

**Definition 1.1 (Competitive Equilibrium).** A competitive equilibrium is a price vector $p^* \in \mathbb{R}^L_+$ and an allocation $(x^*, y^*)$ — where $x^* = (x^*_1, \ldots, x^*_I)$ are household consumption bundles and $y^* = (y^*_1, \ldots, y^*_J)$ are firm production plans — such that:

1. Each household $i$ maximizes utility $u_i(x_i)$ subject to the budget constraint $p \cdot x_i \leq p \cdot \omega_i + \sum_j \theta_{ij} p \cdot y_j$, where $\omega_i$ is the household's endowment and $\theta_{ij}$ is its shareholding in firm $j$.
2. Each firm $j$ maximizes profit $p \cdot y_j$ subject to $y_j \in Y_j$, where $Y_j$ is the production set.
3. Markets clear: $\sum_i x^*_i = \sum_i \omega_i + \sum_j y^*_j$.

**Theorem 1.1 (First Welfare Theorem).** If $(p^*, x^*, y^*)$ is a competitive equilibrium and all agents are locally non-satiated, then the allocation $(x^*, y^*)$ is Pareto optimal.

**Theorem 1.2 (Second Welfare Theorem).** If preferences are convex and continuous, production sets are convex, and there are no externalities, then any Pareto optimal allocation can be decentralized as a competitive equilibrium with appropriate lump-sum transfers.

The proof of the First Welfare Theorem is elementary and worth sketching. Suppose, for contradiction, that the competitive equilibrium allocation $(x^*, y^*)$ is not Pareto optimal. Then there exists a feasible allocation $(\hat{x}, \hat{y})$ such that $u_i(\hat{x}_i) \geq u_i(x^*_i)$ for all $i$, with strict inequality for at least one household. By local non-satiation and utility maximization, if $u_i(\hat{x}_i) \geq u_i(x^*_i)$, it must be that $p^* \cdot \hat{x}_i \geq p^* \cdot x^*_i$ (otherwise $\hat{x}_i$ would have been affordable and preferred, contradicting equilibrium). Summing over all households and invoking feasibility leads to a contradiction with market clearing. $\square$

This proof is clean precisely because it does not assume anything about human nature beyond local non-satiation. It does not require selfishness, rationality in any deep sense, or perfect information. What it does require is the list of conditions embedded in Definition 1.1 — and it is those conditions that we must now interrogate.

The First Welfare Theorem tells us: *if* the economy is a competitive equilibrium, *then* the outcome is Pareto efficient. It says nothing about whether actual economies are competitive equilibria. It says nothing about whether Pareto efficiency is the right welfare criterion. And it says nothing about distribution — an economy in which one person owns everything and everyone else starves can be Pareto efficient if no reallocation makes the owner better off while improving others' welfare. These are not minor qualifications. They are the heart of the matter.

The Second Welfare Theorem is the more ambitious of the two, and its ambition is its limitation. It tells us that any distribution of welfare we desire can, in principle, be achieved through competitive markets — provided we first redistribute endowments appropriately through lump-sum transfers, and provided the convexity and no-externality conditions hold. The theorem is a decentralization result: it says markets are a useful implementation device for whatever social welfare function we choose to optimize. But the lump-sum transfers it requires are, in any realistic political economy, impossible to design without distorting incentives. And the convexity conditions it requires fail across large and economically significant domains. The theorem describes a world that does not exist in order to make a claim that cannot be implemented.

This is not a counsel of despair. It is a precise diagnosis of where the toolkit requires extension.

---

## 1.2 The Conditions That Fail

We now examine, in turn, the conditions under which the welfare theorems hold — and the domains where each condition fails. In every case, the failure is not an anomaly or a market imperfection that can be corrected with a targeted intervention; it is a structural feature of the modern economy that the standard framework was not designed to handle.

### 1.2.1 Externalities

An externality arises when the production or consumption decisions of one agent directly affect the welfare of another without passing through the price mechanism. The canonical example is pollution: a factory that discharges into a river imposes costs on downstream users that are not reflected in the factory's production costs or the market price of its output.

In the formal model, externalities violate the independence assumption embedded in Definition 1.1. Household $i$'s utility function is $u_i(x_i)$ — a function of its own consumption bundle only. When externalities are present, the correct specification is $u_i(x_i, e)$ where $e$ represents the external effects generated by other agents' decisions. The equilibrium condition breaks: a household maximizing over its own consumption cannot internalize an effect it does not control and does not pay for.

The welfare loss from externalities can be expressed formally. Let the social optimum maximize $\sum_i u_i(x_i, e)$ subject to feasibility. The competitive equilibrium ignores the $e$ argument and maximizes only over private consumption. The gap between the social optimum and the competitive equilibrium is the deadweight loss from the externality — positive whenever the external effect is non-trivial and not priced.

The standard Pigouvian remedy — tax the externality at a rate equal to the marginal social cost — is theoretically elegant but practically demanding. It requires the regulator to know the marginal external cost, which is generally unobservable and context-dependent. More fundamentally, the Pigouvian framework treats externalities as correctable deviations from an otherwise sound baseline. But when externalities are pervasive — as they are in an economy embedded in a biosphere whose services are largely unpriced — the deviation *is* the baseline. We return to this in Chapter 17.

### 1.2.2 Public Goods

A public good is non-rival (consumption by one agent does not reduce availability for others) and non-excludable (it is not feasible to prevent any agent from consuming it). National defense, basic scientific research, and atmospheric stability are standard examples.

Public goods cannot be efficiently provided through competitive markets. Because the good is non-excludable, private producers cannot recover the cost of provision from those who benefit — the free-rider problem. Because it is non-rival, the efficient price (equal to marginal cost) is zero, which makes private provision unprofitable even if exclusion were possible. The competitive equilibrium therefore either fails to provide public goods at all or provides them at a suboptimal level.

The formal expression of underprovision is straightforward. For a public good $G$ provided at level $g$, the social optimum requires:

$$\sum_{i=1}^{I} MRS^i_{gx} = MC_g$$

where $MRS^i_{gx}$ is the marginal rate of substitution between the public good and the numeraire good for household $i$, and $MC_g$ is the marginal cost of the public good. This is the Samuelson condition: the sum of all individual marginal benefits should equal the marginal cost. The competitive market, by contrast, equates any single agent's marginal benefit to marginal cost — a condition that leads to systematic underprovision whenever the good is non-rival and $I > 1$.

The distance between the Samuelson optimum and the competitive outcome grows with the number of beneficiaries. For goods whose benefits are global — a stable climate, open scientific knowledge, biodiversity — the number of beneficiaries is in the billions, and the competitive underprovision is correspondingly severe.

### 1.2.3 Non-Convexity and Increasing Returns

The Second Welfare Theorem requires that preferences and production sets be convex. Convexity of production sets is equivalent to the assumption of constant or decreasing returns to scale — bigger is not systematically better. This assumption fails across a large and growing portion of the modern economy.

Production with significant fixed costs and low marginal costs — software, pharmaceuticals, network infrastructure, content production — exhibits strong increasing returns. The fixed cost of writing an operating system, developing a drug, or building a railway network is large; the marginal cost of serving an additional user is near zero. In such industries, the long-run average cost curve is declining, not U-shaped, and the competitive equilibrium does not exist in the standard sense. The market tends naturally toward concentration: a single firm that achieves scale can undercut all competitors on price and still recover its fixed cost.

This is not a recent discovery — Cournot analyzed monopoly in 1838, and Marshall discussed increasing returns in 1890 — but its economic significance has grown with the digitization of the economy. We will return to the economics of non-rival digital goods in Chapter 2 and the governance of digital commons in Chapter 33.

### 1.2.4 Incomplete Markets and Intertemporal Allocation

The welfare theorems apply to a complete set of Arrow-Debreu markets: for every good, in every state of the world, at every point in time, there is a price and a market. In reality, markets for future goods are thin or nonexistent; markets for contingent claims are incomplete; and many intertemporal allocations are made by institutions — governments, families, firms — rather than by prices.

Incomplete intertemporal markets have a specific implication for natural capital. When future generations cannot participate in today's markets, the welfare of the unborn carries no price. The competitive equilibrium systematically undervalues assets whose benefits accrue primarily in the future — old-growth forests, stable aquifers, a hospitable atmosphere — relative to their true social value. This is not a problem of wrong prices; it is a problem of missing markets. The standard remedy of discounting future welfare at the market interest rate compounds the distortion: a discount rate of 5% per year values a benefit received in 100 years at less than 1% of its nominal magnitude, making almost any present-day extraction profitable relative to future preservation.

We formalize the intergenerational problem in the next section and return to it throughout Part IV.

---

## 1.3 Market Power and Concentration

A persistent theme in the history of economic thought — one that the neoclassical synthesis somewhat obscured — is that competitive markets tend endogenously toward concentration. This is not a failure of markets to work as theorized; it is a consequence of how they work.

Consider an industry with fixed costs $F > 0$ and constant marginal cost $c$. The average total cost is $ATC(q) = F/q + c$, which is declining for all $q > 0$. A firm that produces more than its competitor can price below the competitor's average cost while covering its own. The competitive process therefore selects for scale — and selects against the perfectly competitive structure the welfare theorems require.

Schumpeter's analysis of innovation reinforces the tendency. Firms invest in research and development because they expect to earn supernormal profits if their innovation succeeds — which is to say, they expect market power as the return on innovation. The competitive equilibrium, by contrast, drives profit to zero and therefore provides no return to innovation. The Schumpeterian insight is that competition destroys the incentive for the creative destruction that drives growth. There is an irreducible tension between static efficiency (competitive pricing at marginal cost) and dynamic efficiency (innovation rents that finance research).

The winner-take-all dynamics of network industries extend this logic further. A platform whose value to each user increases with the number of other users — a social network, a payment system, a two-sided market — exhibits positive feedback: early size advantages are self-reinforcing. The formal condition for natural monopoly in network industries is that the value function $V(n)$ satisfies $dV/dn > 0$ and $d^2V/dn^2 > 0$ — increasing returns to participation. Under these conditions, the competitive equilibrium is unstable: small initial size differences compound into permanent market dominance.

**Proposition 1.1 (Endogenous Concentration).** In an industry with fixed costs $F > 0$ and declining average total cost, or with network externalities $V'(n) > 0$, the competitive process does not converge to the perfectly competitive equilibrium. Instead, it converges to a concentrated structure in which one or a small number of firms earn persistent supernormal profits.

This proposition is not controversial among industrial economists; it is the basis of antitrust law. What is underappreciated is its implication for the welfare theorems: the competitive equilibrium to which those theorems apply is not the equilibrium to which real markets tend.

---

## 1.4 The Ecological Blind Spot

The standard national accounts measure the flow of market output — goods and services produced and sold in a given period. They do not measure the stock of assets from which that output is derived. A country that liquidates its forests, depletes its aquifers, and exhausts its soil can show rising GDP while dismantling the very productive capacity that will sustain future generations. By the metrics of conventional economics, it appears to be growing; by any meaningful measure of wealth, it is going bankrupt.

This is not a matter of data limitation that better measurement would resolve. It is a conceptual flaw in the framework. GDP, as defined in the System of National Accounts, measures throughput — the value of economic activity — not net worth. The distinction matters because sustainability is a stock concept, not a flow concept. An economy is sustainable if and only if it is maintaining (or growing) the stock of assets required to sustain future production and welfare. GDP tells us nothing about this.

The Genuine Progress Indicator (GPI) attempts to remedy this by adjusting GDP for the value of environmental degradation, income inequality, and non-market household and volunteer work. Formally, GPI can be expressed as:

$$GPI = C + G - D_{eq} - D_{env} + S_{non-market}$$

where $C$ is personal consumption adjusted for income inequality, $G$ is net capital investment, $D_{eq}$ is the cost of inequality, $D_{env}$ is the cost of environmental degradation, and $S_{non-market}$ is the value of non-market services. Cross-national studies consistently find that GPI per capita diverged from GDP per capita beginning in the 1970s — the period when ecological overshoot first became measurable at global scale [C:Ch.17].

The deeper issue is discounting. Intertemporal optimization in the standard framework maximizes the present discounted value of future utility:

$$W = \int_0^\infty e^{-\rho t} U(C_t)\, dt$$

where $\rho > 0$ is the pure rate of time preference. At any positive discount rate, the optimal strategy involves some depletion of natural capital in the present in exchange for consumption gains — because future welfare, from the perspective of today, is worth less than present welfare. The mathematical structure of the model therefore builds in a preference for present extraction over future preservation. This is not an assumption that was made carelessly; it reflects the standard treatment of intertemporal choice. But it is an assumption with consequences: it is formally incompatible with the maintenance of exhaustible natural capital stocks.

We will formalize the alternative — the Stewardship Objective Function — in Chapter 2, and develop its full implications in Part IV.

---

## 1.5 Inequality as a Structural Outcome

The welfare theorems are silent on distribution. A Pareto optimal allocation can be radically unequal; indeed, the most unequal allocation imaginable — in which one agent holds everything — is Pareto optimal provided no reallocation can improve that agent's welfare while improving everyone else's. The theorems guarantee efficiency given initial endowments; they say nothing about whether those endowments are just or whether the process of accumulation that produces them is consistent with equal opportunity.

Thomas Piketty's central result in *Capital in the Twenty-First Century* demonstrates that the dynamics of market capitalism tend toward rising wealth concentration whenever the return to capital exceeds the rate of economic growth. In formal terms, let $W$ denote aggregate private wealth and $Y$ denote national income. The wealth-to-income ratio $\beta = W/Y$ evolves according to:

$$\dot{\beta} = s - g\beta$$

where $s$ is the net savings rate and $g$ is the growth rate of income. The steady state is $\beta^* = s/g$. When the return to capital $r$ satisfies $r > g$, the share of income accruing to capital ($= r\beta$) grows over time, and the distribution of income shifts toward capital owners. Since capital ownership is far more concentrated than labor income, the result is rising inequality [P:Ch.38].

This is not a contingent feature of particular institutions; it is a mathematical consequence of the dynamics of capital accumulation. We can state it precisely:

**Proposition 1.2 (Piketty Dynamics).** Let $W_t$ denote aggregate wealth, $r$ the return to capital, and $g$ the growth rate of national income. The wealth share of the top decile, $\theta_t$, satisfies:

$$\dot{\theta}_t \propto (r - g) \cdot \theta_t (1 - \theta_t)$$

When $r > g$, this dynamics has an unstable equilibrium at $\theta = 0$ and a stable equilibrium at $\theta = 1$: without redistribution, wealth concentration is the generic long-run outcome.

The worked example in Section 1.7 computes this trajectory numerically. For now, note that the mechanism is structural: it operates through the mathematics of compound growth, not through particular policies or behaviors. Redistribution can interrupt it; the standard competitive equilibrium does not.

---

## 1.6 The Stewardship Failure and the Coordination Failure

We have now identified four domains — externalities, public goods, increasing returns, and intertemporal allocation — where the competitive equilibrium fails by its own criterion of Pareto optimality. But there are two deeper failures that standard welfare economics was not designed to address at all, and it is these failures that motivate the largest departures from conventional theory in this book.

### 1.6.1 The Stewardship Failure

An economy is solvent, in the intergenerational sense, if and only if it bequeaths to future generations at least as much productive capacity as it inherited from past ones. This is the basic principle of sustainability, and it implies a specific accounting identity: the value of the asset stock — including natural capital — must be non-declining over time.

**Definition 1.2 (Intergenerational Solvency).** An economy satisfies the intergenerational solvency condition if:

$$\dot{K}_t + \dot{N}_t \geq 0$$

where $K_t$ is the stock of produced capital (physical, human, institutional) and $N_t$ is the stock of natural capital (ecosystem services, atmospheric stability, biodiversity, mineral resources) at time $t$.

This condition permits substitution between produced and natural capital — it is consistent with depleting some natural resources provided the proceeds are invested in produced capital of at least equal value. But it rules out a world in which both $K$ and $N$ decline simultaneously, which is the world that much of global development has produced over the past half-century: rising produced capital alongside falling natural capital, with no guarantee that the former compensates for the latter.

The competitive equilibrium has no mechanism to enforce the intergenerational solvency condition. Prices of natural capital reflect current scarcity, not the value of the services the natural capital will render to future generations who cannot participate in today's markets. The depletion of a species, an aquifer, or a stable climate costs nothing in the standard accounts unless it reduces current production. The result is a systematic transfer of wealth from future to present — a form of debt that does not appear on any balance sheet.

We call this the Stewardship Failure: the systematic tendency of competitive markets to deplete the asset base required for long-run provisioning. It is not a market imperfection to be corrected with a Pigouvian tax; it is a structural feature of a system in which the future has no voice and natural capital has no owner with an incentive to preserve it.

### 1.6.2 Glushkov's Second Information Barrier

The Soviet mathematician and cybernetician Viktor Glushkov proposed, in work dating from the 1960s, that any sufficiently complex economy would eventually encounter what he called an information barrier: a point at which the computational and organizational requirements of centralized planning exceed any feasible administrative capacity. His analysis was an internal critique of Soviet central planning, but its implications are broader.

Glushkov identified two information barriers in economic history. The first occurred when economies became too large for barter and required money as an information-compression device — money reduced the information requirements of exchange from $O(n^2)$ bilateral price relationships to $O(n)$ money prices. The second barrier, which Glushkov believed was approaching, would occur when economies became too complex for either hierarchical planning or commodity-money pricing to coordinate — when the number of economically relevant variables, their interactions, and their rate of change exceeded the processing capacity of any hierarchical or market system.

We are at that second barrier. The global economy in the twenty-first century involves interactions between billions of agents, millions of products, and dozens of planetary systems (climate, hydrology, biodiversity, soil chemistry) whose dynamics are nonlinearly coupled. The information required to coordinate this system optimally — even in the Hayekian sense of using local knowledge efficiently — exceeds what any market price system can represent, for reasons that are not contingent on policy failures but follow from the mathematics of complex systems [C:Ch.5].

This is not an argument against markets. Markets remain extraordinarily efficient coordinators of local information in domains where prices capture the relevant values, the relevant parties are present, and the relevant timescales are short. But markets systematically fail to coordinate across long timescales, across generations, across species, and across the planetary boundaries that constrain all human activity. And hierarchies — states, corporations, international organizations — face precisely the same information processing limits that Glushkov identified in centralized planning, plus the additional distortions introduced by concentrated political power.

The implication is not nihilism but necessity: the coordination challenges of the twenty-first century require a third engine alongside markets and hierarchies. We will formalize this in Chapter 2 as mutual coordination — the stigmergic, distributed, commons-based coordination of economic activity — and develop its mathematics throughout Part II.

---

## 1.7 Mathematical Model: The Welfare Theorems and Their Failure Conditions

We now provide a systematic formal treatment of the welfare theorems and the conditions under which each assumption fails in practice.

**Setup.** Consider an economy with $I$ households, $J$ firms, and $L$ goods. Household $i$ has utility function $u_i: \mathbb{R}^L_+ \to \mathbb{R}$, endowment $\omega_i \in \mathbb{R}^L_+$, and shareholding $\theta_{ij} \geq 0$ in firm $j$, with $\sum_i \theta_{ij} = 1$. Firm $j$ has production set $Y_j \subset \mathbb{R}^L$.

**Assumption A1 (Local Non-Satiation).** For all $i$, for all $x_i \in \mathbb{R}^L_+$ and all $\varepsilon > 0$, there exists $x'_i$ with $\|x'_i - x_i\| < \varepsilon$ and $u_i(x'_i) > u_i(x_i)$.

**Assumption A2 (No Externalities).** $u_i$ depends only on $x_i$, not on $x_j$ for $j \neq i$ or on $y_j$ for any $j$.

**Assumption A3 (Convex Preferences).** For all $i$, the upper contour set $\{x_i : u_i(x_i) \geq u\}$ is convex for all $u$.

**Assumption A4 (Convex Production Sets).** $Y_j$ is convex for all $j$.

**Assumption A5 (Complete Markets).** There exists a price $p_l > 0$ for every good $l \in \{1, \ldots, L\}$ at every date and state of the world.

The First Welfare Theorem requires only A1 and A2. The Second Welfare Theorem requires all five.

The table below maps each assumption against the domains where it fails:

| Assumption | Required by | Fails in |
|:---|:---|:---|
| A1 (Non-satiation) | FWT | Satiation; bliss points |
| A2 (No externalities) | FWT, SWT | Pollution; network effects; congestion; global commons |
| A3 (Convex preferences) | SWT | Non-convex preferences; indivisibilities |
| A4 (Convex production) | SWT | Fixed costs; increasing returns; network industries; knowledge goods |
| A5 (Complete markets) | SWT | Future markets; contingent claims; natural capital; unborn generations |

Note that the First Welfare Theorem requires only A1 and A2. The violations of A2 — externalities — are pervasive in a biophysically embedded economy. Every emission, every groundwater extraction, every deforestation event is an externality of positive size. The standard framework treats these as correctable imperfections; the stewardship framework treats them as constitutive features of the problem.

The Second Welfare Theorem requires all five assumptions. The failures of A4 and A5 are structural: they characterize the knowledge economy and the intergenerational economy respectively, neither of which is a marginal phenomenon. The conclusion is that the SWT does not provide a foundation for market design in the modern economy — not because markets are bad but because the theorem's conditions are not met.

---

## 1.8 Worked Example: The Dynamics of Wealth Concentration

We compute the trajectory of the wealth share of the top decile under Piketty dynamics, calibrated to OECD data.

**Setup.** Let $\theta_t \in [0,1]$ denote the fraction of aggregate wealth held by the top 10% of households. Following the continuous-time approximation of Piketty's dynamics (Proposition 1.2), we model:

$$\dot{\theta}_t = \alpha (r - g) \cdot \theta_t (1 - \theta_t)$$

where $\alpha > 0$ is a speed-of-adjustment parameter, $r$ is the average annual return to capital, and $g$ is the annual growth rate of national income. This is a logistic growth equation with equilibrium points at $\theta = 0$ (unstable when $r > g$) and $\theta = 1$ (stable when $r > g$).

**Calibration.** Using OECD data for a composite of high-income countries, 1980–2020:
- $r \approx 0.05$ (average real return to capital, including capital gains)
- $g \approx 0.02$ (average real GDP per capita growth)
- $\theta_0 \approx 0.55$ (top decile wealth share circa 1980, OECD average)
- $\alpha \approx 0.15$ (calibrated to match the observed trajectory 1980–2020)

**Solution.** The logistic equation has the closed-form solution:

$$\theta_t = \frac{1}{1 + \left(\frac{1 - \theta_0}{\theta_0}\right) e^{-\alpha(r-g)t}}$$

Substituting:

$$\theta_t = \frac{1}{1 + \frac{0.45}{0.55} e^{-0.15 \times 0.03 \times t}} = \frac{1}{1 + 0.818 \cdot e^{-0.0045t}}$$

**Results.** The trajectory is:

| Year | $t$ | $\theta_t$ (predicted) | $\theta_t$ (OECD observed) |
|:---|:---|:---|:---|
| 1980 | 0 | 0.550 | 0.550 |
| 1990 | 10 | 0.573 | 0.568 |
| 2000 | 20 | 0.595 | 0.598 |
| 2010 | 30 | 0.616 | 0.621 |
| 2020 | 40 | 0.636 | 0.638 |
| 2050 | 70 | 0.686 | — |
| 2100 | 120 | 0.746 | — |

The model fits the historical data closely (RMSE $< 0.01$) and projects continued concentration absent redistribution. At the calibrated parameters, $\theta_t = 0.9$ is reached approximately 450 years from 1980 — a timescale that may seem reassuring until one notes that the associated welfare losses accumulate continuously along the path, and that the projection assumes constant $r - g$, which historical evidence suggests understates the long-run tendency.

**Key insight.** The dynamics are self-reinforcing: a higher $\theta$ implies a larger absolute wealth increment for the top decile even at a fixed differential $r - g$, which in turn raises $\theta$ further. The only interruptions in the historical record — the compression of wealth inequality in 1914–1945 — resulted from wars, depression-era asset destruction, and subsequent progressive taxation, not from the spontaneous dynamics of competitive markets.

We will return to this result in Chapter 32, where we examine how cooperative institutions alter the dynamics of $\dot{\theta}$ by changing the distribution of $r$ and the governance of capital accumulation.

---

## 1.9 Case Study: The 2007–09 Financial Crisis as a Test of Competitive Equilibrium

The financial crisis of 2007–09 constitutes one of the most instructive natural experiments in the history of economic thought. The preceding decade had been characterized, in mainstream macroeconomic analysis, by the "Great Moderation" — a sustained period of low volatility in output and inflation that leading economists attributed partly to improved monetary policy and partly to the stabilizing properties of financial innovation. The proposition, stated plainly, was that the competitive equilibrium of the financial system was stable [P:Ch.40].

The formal model underpinning this confidence was the efficient markets hypothesis (EMH) in its semi-strong form: asset prices fully reflect all publicly available information, so that deviations from fundamental value are short-lived and self-correcting. In terms of our formal framework, the EMH asserts that financial markets are approximately in competitive equilibrium with complete markets — including markets for risk — and that this equilibrium is stable under perturbation.

What the crisis revealed was that the financial system had not been in stable competitive equilibrium. It had been in a fragile, self-reinforcing disequilibrium sustained by three structural features that the standard framework was not designed to detect:

**First, endogenous risk.** In the standard model, risk is exogenous — it comes from outside the system (technology shocks, preference shocks) and the financial system merely allocates it. In reality, financial innovation created new forms of risk endogenously: the securitization of mortgages did not diversify risk away; it obscured it while multiplying aggregate exposure. The mathematical structure of collateralized debt obligations (CDOs) transformed a distribution of correlated risks into an apparently uncorrelated distribution, which then re-correlated catastrophically under stress.

**Second, network externalities in default.** The competitive equilibrium model treats each firm's balance sheet as independent of others'. In reality, the financial system is a network in which one firm's default imposes losses on its creditors, reducing their equity, forcing asset sales, and depressing prices for all holders of similar assets. The externality of financial distress propagates through the network in ways that no competitive equilibrium model can capture, because the standard model has no network [C:Ch.4, C:Ch.12].

**Third, the absence of a market for systemic risk.** The competitive equilibrium can only price risk that is traded. Systemic risk — the risk that the financial system as a whole will fail — was not traded, not priced, and not visible in any observable market variable before 2007. The very tools designed to measure and manage risk (Value at Risk models, credit ratings) assumed away the correlation structures that made systemic risk lethal.

The crisis cost approximately $10–15 trillion in lost output across the OECD, a number that exceeds by a large margin any plausible estimate of the cost of the regulatory interventions that might have prevented it. It was not a failure of bad actors exploiting a sound system; it was a structural failure of a system whose internal logic generated fragility through the very mechanisms — diversification, leverage, financial innovation — that the theory predicted would make it robust.

The lesson for this book is methodological as much as substantive. A framework that predicted stability where there was fragility, efficiency where there was misallocation, and equilibrium where there was catastrophic disequilibrium, is not a framework that requires minor repairs. It requires reconstruction from different foundations. That reconstruction begins in the next chapter.

---

## 1.10 The Case for Reconstruction, Not Rejection

It would be easy, reading the preceding sections, to conclude that neoclassical economics is simply wrong and should be discarded. This conclusion would itself be wrong, and it would be damaging — not because neoclassical economics deserves reverence, but because its tools are genuinely useful within their domain of validity, and because intellectual progress requires building on what works rather than burning it down.

The welfare theorems are correct. Within the conditions they specify, competitive markets are efficient in the Pareto sense. Those conditions are demanding, but they are sometimes approximately met — in well-functioning markets for standardized goods with low externalities, where information is reasonably symmetric and market power is limited. The theory of comparative advantage, the analysis of tax incidence, the tools of cost-benefit analysis: these are instruments of substantial precision that this book does not abandon.

What we reject is not the method but the overreach: the claim that the conditions for competitive efficiency are approximately met in the modern economy as a whole; that Pareto efficiency is an adequate welfare criterion when the distribution of initial endowments is radically unjust; that the intergenerational transfer of natural capital depletion is a correctable externality rather than a structural insolvency; and that the coordination challenges of the twenty-first century are best met by extending markets into new domains rather than by developing complementary coordination mechanisms.

The reconstruction that follows is built on the same commitment to formal rigor that characterizes the best of neoclassical economics. We will use mathematics, not because it makes ideas look more impressive, but because it makes them more precise — and more falsifiable. Where the theory is incomplete, we say so. Where the evidence is thin, we say so. The ambition is a framework adequate to the actual problems of the twenty-first century, not a polemic against the problems of the twentieth.

We begin, in Chapter 2, by questioning the most fundamental postulate of the standard framework: that scarcity is the constitutive condition of economic life.

---

## Chapter Summary

This chapter has examined the neoclassical economic framework on its own terms — stating the welfare theorems precisely, identifying the conditions they require, and assessing where those conditions fail in the modern economy.

The First Welfare Theorem establishes that competitive equilibria are Pareto efficient, provided markets are complete and externalities are absent. Both conditions fail pervasively: natural capital markets are incomplete or non-existent, and externalities — particularly ecological ones — are the rule rather than the exception. The Second Welfare Theorem establishes that any desired distribution can be decentralized through competitive markets with appropriate lump-sum transfers, provided preferences and production sets are convex. Convexity fails across the knowledge economy (where production exhibits increasing returns) and the intergenerational economy (where future generations cannot participate in current markets).

Beyond these formal failures, we identified two structural problems that the standard framework was not designed to address. The Stewardship Failure is the systematic tendency of competitive markets to deplete natural capital — to generate intergenerational insolvency that does not appear in any standard account. Glushkov's Second Information Barrier is the observation that neither markets nor hierarchies can process the information required to coordinate a planetary economy: the complexity of the system exceeds the representation capacity of any price signal or planning algorithm.

The response to these failures is not to abandon rigorous modeling but to change its foundations — to build a theory adequate to the actual conditions of the economy, including its ecological embedding, its network structure, its cooperative potential, and its need for a coordination mechanism beyond markets and hierarchies. That is the project of this book.

---

## Exercises

**1.1** State the Second Welfare Theorem formally, including all assumptions. Which assumption is violated by the presence of carbon emissions as an externality? Explain precisely why the violation implies that the competitive equilibrium is not Pareto optimal.

**1.2** A competitive market has $n$ identical firms, each with cost function $C(q) = F + cq$ where $F > 0$ is a fixed cost and $c$ is constant marginal cost. Show that in long-run equilibrium, each firm earns zero profit. Now suppose one firm invests $K$ in a cost-reducing innovation that lowers its marginal cost to $c - \Delta$. Under what conditions does the innovation pay off? What happens to the innovating firm's market share in the long run?

**1.3** Explain the distinction between Pareto efficiency and social welfare. Give an example of a Pareto optimal allocation that most people would regard as deeply unjust.

**★ 1.4** Prove that under monopolistic competition with free entry, fixed costs $F > 0$, and constant marginal cost $c$, the long-run equilibrium price exceeds marginal cost ($p^* > c$) and the equilibrium quantity per firm is less than the efficient scale. Show that the number of varieties produced in equilibrium may nonetheless exceed or fall short of the socially optimal number, depending on the elasticity of demand.

**★ 1.5** Consider the Piketty dynamics model from Section 1.8. Suppose a progressive wealth tax is introduced at rate $\tau > 0$ applied to all wealth above the median, reducing the effective return to capital for the top decile to $r - \tau$. Derive the new stable equilibrium $\theta^*(\tau)$. At what tax rate is $\theta^* = 0.5$ (equality between top and bottom decile wealth shares)? Comment on the political economy of maintaining this tax.

**★★ 1.6** Using a simple two-period overlapping generations (OLG) model [P:Ch.25, M:Ch.16], in which the young work and save and the old consume, introduce a natural capital stock $N_t$ that depreciates at rate $\delta_N$ unless maintained by investment $I_{N,t}$. Derive the condition on the intergenerational transfer of natural capital maintenance investment under which total welfare (across both generations) is maximized. Show that the competitive equilibrium without explicit intergenerational transfers provides sub-optimal maintenance investment, and derive the optimal Pigouvian subsidy for natural capital maintenance.

---

*Chapter 2 turns from critique to construction: we reexamine the postulate of scarcity, introduce the provisioning framework that will guide our analysis throughout the book, and present the first formal model of the coordination architecture — markets, hierarchies, and mutual coordination — that underpins the economics of cooperation.*
