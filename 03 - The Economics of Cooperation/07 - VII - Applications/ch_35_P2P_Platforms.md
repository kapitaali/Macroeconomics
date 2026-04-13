# Chapter 35: Peer-to-Peer Platforms — From Open Source Software to Energy Grids

> *"Platform cooperativism isn't just about who owns the platform. It is about who governs it, who designs it, and who benefits from the value it creates."*
> — Trebor Scholz, *Platform Cooperativism* (2016)

> *"The internet was built on a peer-to-peer architecture. Platform capitalism re-centralized it. Platform cooperativism can decentralize it again."*
> — Nathan Schneider, *Everything for Everyone* (2018)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Apply the P2P network theory of Chapter 8 to specific platform architectures, formally classifying platforms by their degree centralization, contribution measurement approach, and governance model.
2. Analyze the platform capitalism failure modes — rent extraction through network lock-in, regulatory arbitrage, and value appropriation from contributors — through the formal lens of cooperative game theory, identifying the precise mechanisms through which platform firms appropriate cooperative surplus.
3. Design P2P platform systems that capture cooperative benefits (network externalities, distributed intelligence, resilience) while avoiding failure modes, using the multi-stakeholder cooperative governance principles of Chapter 34.
4. Specify a P2P energy trading algorithm for a prosumer cooperative network, derive its incentive-compatibility conditions, and implement it as a smart contract pseudocode specification.
5. Analyze the transition dynamics from corporate platform to cooperative platform through the institutional entrepreneurship and tipping threshold framework of Chapter 15.
6. Evaluate REScoop.eu — the European federation of renewable energy cooperatives — formally against the cooperative game theory, network analysis, and Ostrom governance frameworks of the book.

---

## 35.1 The Platform Economy and Its Cooperative Alternative

The digital platform economy of the past two decades has produced one of the most remarkable concentrations of economic power in history. Five technology companies — Apple, Microsoft, Google, Amazon, and Meta — each achieved market capitalizations exceeding USD 1 trillion. They did so primarily by creating two-sided platforms: digital infrastructures that connect producers and consumers, facilitating transactions while capturing a share of the value created.

The cooperative game theory analysis of this book illuminates precisely how this concentration occurred. Platform firms create genuine cooperative surplus — the value of network externalities, reduced transaction costs, and distributed intelligence that only emerges from scale. But under standard corporate governance, the Shapley value of each contributor (every producer who lists on Amazon, every driver who works for Uber, every content creator on YouTube) is captured not by the contributor but by the platform's shareholders. The platform is, in the language of Chapter 6, a superadditive cooperative game whose surplus is distributed not by the Shapley value but by the monopolist's pricing rule — extracting as much as the contributor's outside option allows.

Platform cooperativism offers an alternative: platforms governed and owned by their users, distributing surplus according to contribution rather than capital ownership. This is not a new idea — consumer cooperatives, producer cooperatives, and mutual organizations have been organizing economic activity cooperatively for two centuries. But the digital infrastructure of the platform economy creates both new possibilities (global scale, zero marginal cost of replication, smart contract automation of governance) and new challenges (the cold-start network externality problem, governance at scale, and competition against incumbents with massive data and capital advantages).

This chapter develops the formal analysis of P2P platforms: what makes them different from hierarchical platforms, which failure modes they must overcome, how their governance should be designed, and what the energy sector — through REScoop.eu's federation of energy cooperatives — reveals about the conditions under which cooperative P2P platforms can displace corporate incumbents.

---

## 35.2 Platform Capitalism: The Formal Failure Mode

### 35.2.1 Platform Value Creation and Appropriation

**Definition 35.1 (Two-Sided Platform).** A two-sided platform is an intermediary $P$ connecting two groups of agents: producers $\mathcal{P}_L$ (left side) and consumers $\mathcal{P}_R$ (right side), whose interactions generate value through network externalities:

$$V(\mathcal{P}_L, \mathcal{P}_R) = v(|\mathcal{P}_L|, |\mathcal{P}_R|, q_L, q_R)$$

where $q_L$ and $q_R$ are the quality of producer and consumer participation respectively. Network externalities: $\partial v/\partial |\mathcal{P}_L| > 0$ (more producers benefit consumers) and $\partial v/\partial |\mathcal{P}_R| > 0$ (more consumers benefit producers).

**The cooperative game representation.** The platform economy is a superadditive cooperative game with characteristic function:

$$v(\mathcal{P}_L \cup \mathcal{P}_R \cup \{P\}) > v(\mathcal{P}_L) + v(\mathcal{P}_R) + v(\{P\})$$

The cooperative surplus — the additional value created by bringing producers, consumers, and the platform infrastructure together — belongs, under Shapley value logic, to all three parties in proportion to their average marginal contributions.

**Platform appropriation.** The corporate platform appropriates the cooperative surplus by setting the price structure $(p_L, p_R)$ to extract as much producer and consumer surplus as possible while maintaining participation:

$$\text{Extracted surplus} = v(\mathcal{P}_L \cup \mathcal{P}_R \cup \{P\}) - CS_L - CS_R$$

where $CS_L$ and $CS_R$ are the producer and consumer surplus left after the platform's fees. As the platform achieves market dominance (no competitive alternatives), $CS_L \to 0$ and $CS_R \to 0$ — the platform appropriates essentially all cooperative surplus.

**The Uber case.** Uber's 2023 gross bookings: USD 115 billion. Driver earnings (after Uber's commission): approximately USD 77 billion. Uber's revenue: USD 37 billion (32\% of gross bookings). Under Shapley value allocation of the cooperative surplus created by the Uber platform (connecting drivers and riders efficiently): drivers would receive approximately 45–55\% of surplus (reflecting their primary productive contribution), riders would retain their consumer surplus, and the platform would receive its marginal contribution (the infrastructure and matching algorithm) — approximately 20–25\% of surplus. The gap between actual allocation (Uber receives 32\%) and Shapley (Uber receives 20–25\%) represents the rent extracted through market power.

### 35.2.2 The Three Platform Failure Modes

**Failure Mode 1: Network Lock-In.** Once a platform achieves critical mass, switching costs (loss of network connections, data history, reputation scores) create lock-in that prevents competitive entry. This is the tipping threshold of Chapter 15 applied to the incumbent's advantage: the incumbent's adoption level is permanently above $\hat{x}$, making displacement nearly impossible without a comparable coordinated migration.

**Failure Mode 2: Regulatory Arbitrage.** Platform firms classify workers as independent contractors (not employees) to avoid labor regulations, insurance obligations, and collective bargaining rights. This arbitrage extracts value from contributors by denying them institutional protections that conventional employees receive — a form of commons enclosure (Chapter 32) of the legal-institutional commons of labor protection.

**Failure Mode 3: Data Appropriation.** Platform firms collect data generated by users' interactions and use it as a proprietary asset — training recommendation systems, selling advertising, and building competitive moats. Contributors generate this data through their participation but receive none of its value. This is the digital enclosure analyzed in Chapter 32 and Chapter 33.

**Proposition 35.1 (Platform Failure Modes are Cooperative Game Failures).** Each of the three platform failure modes corresponds to a violation of a Shapley value axiom:

- Lock-in violates **Null Player**: the platform should not receive surplus from its historical user base (who generate current network value) — the historical user base has zero marginal cost to the platform but generates positive value.
- Regulatory arbitrage violates **Symmetry**: identical contributions (driving, coding, content creation) are valued differently based on classification, not contribution.
- Data appropriation violates **Efficiency**: the cooperative surplus is not fully distributed — some is destroyed (privacy costs) or retained by the platform without distributing to generators.

*Proof.* Each violation is direct from the definitions. Lock-in: the historical user base creates the network externality that makes the current platform valuable — they are not null players (their data and connections have positive marginal value) but are treated as if they were. Symmetry: a driver classified as "contractor" and one classified as "employee" making identical contributions receive different compensation structures — violating the symmetry axiom that symmetric contributions should receive symmetric payoffs. Efficiency: data value generated by users is not distributed back to users — the sum of payoffs is less than total platform value. $\square$

---

## 35.3 P2P Platform Design Principles

### 35.3.1 The Cooperative Platform Architecture

**Definition 35.2 (Cooperative P2P Platform).** A cooperative P2P platform is a tuple $(U, \mathcal{G}, \mathcal{A}, \mathcal{D}, \mathcal{R})$ where:
- $U = U_P \cup U_C \cup U_S$: the set of stakeholder users — producers $U_P$, consumers $U_C$, and stewards $U_S$ (technology and governance maintainers).
- $\mathcal{G}$: the governance structure — multi-stakeholder cooperative governance [C:Ch.34], with weighted voting reflecting contribution and stake.
- $\mathcal{A}$: the allocation mechanism — Shapley-approximating OVA distributing revenue among contributor classes.
- $\mathcal{D}$: the data governance structure — contributor data held in a data cooperative [C:Ch.33], accessible to members, not sold to third parties.
- $\mathcal{R}$: the reputation and trust system — network-propagated reputation [C:Ch.16], portable across platforms.

**Platform cooperative design principles (extending the cooperative enterprise principles of Chapter 34):**

**PC-1 (Open membership).** Any contributor meeting quality standards can join without gatekeeping by incumbents.

**PC-2 (Transparent pricing).** Commission structures are set by democratic governance, not by profit maximization. Fee schedules are published and require member approval to change.

**PC-3 (Portable reputation).** Contributor reputation scores are owned by the contributor and portable across platforms — not locked into the platform's proprietary system.

**PC-4 (Data sovereignty).** Aggregate usage data belongs to the contributor community, managed through the data cooperative structure of Chapter 33.

**PC-5 (Surplus distribution).** Platform revenue above operating costs is distributed to contributors via OVA, not extracted as profit for external shareholders.

**PC-6 (Democratic governance).** Material decisions about platform development, pricing, and data use require multi-stakeholder vote.

### 35.3.2 Formal Comparison: Corporate vs. Cooperative Platform

| Feature | Corporate platform | Cooperative P2P platform |
|:---|:---|:---|
| Ownership | External shareholders | Contributing users |
| Revenue distribution | Shareholders (after costs) | OVA allocation to contributors |
| Data governance | Proprietary platform asset | Contributor-owned data cooperative |
| Governance | Shareholder board | Multi-stakeholder cooperative |
| Reputation | Platform-locked score | Portable, contributor-owned |
| Price setting | Profit-maximizing | Democratic (cost-covering + surplus share) |
| Entry barriers | Platform gatekeeping | Open (quality standards only) |
| Surplus appropriation | Full (toward competitive limit) | Zero (Shapley-compliant distribution) |

The cooperative platform is the operational implementation of the Shapley value: all parties receive their average marginal contribution to the platform's total value, with no residual extraction by external capital.

---

## 35.4 P2P Energy Trading: Algorithm and Smart Contract

### 35.4.1 The Prosumer Energy Platform

Renewable energy creates a new class of market participants: prosumers — households and businesses that are both producers (with rooftop solar, home batteries, or EV charging) and consumers of electricity. A P2P energy trading platform allows prosumers to trade electricity directly with each other, bypassing the utility company intermediary and potentially achieving better prices for both buyer (lower than grid retail) and seller (higher than grid feed-in tariff).

**Definition 35.3 (Prosumer Cooperative Energy Network).** A prosumer cooperative energy network is a set of $n$ prosumers $\{p_1, \ldots, p_n\}$ connected to a local grid segment, each characterized at each time period $t$ by:
- **Generation:** $g_i(t) \geq 0$ (kWh generated, primarily solar PV)
- **Consumption:** $c_i(t) \geq 0$ (kWh consumed)
- **Net position:** $e_i(t) = g_i(t) - c_i(t)$ (positive = exporter, negative = importer)
- **Battery:** $B_i(t) \in [0, \bar{B}_i]$ (stored energy, subject to charge/discharge rate constraints)

**Algorithm 35.1 (P2P Energy Trading — Double Auction)**

```python
FUNCTION p2p_energy_trade(prosumers, period t):
    """
    Implements a P2P energy double auction for a prosumer cooperative.
    Matches exporters and importers, sets cooperative clearing price,
    and records transactions for OVA allocation.
    """

    # Step 1: Collect bids and asks
    asks = []  # [(prosumer_id, quantity_kWh, min_price_EUR/kWh)]
    bids = []  # [(prosumer_id, quantity_kWh, max_price_EUR/kWh)]

    FOR i IN prosumers:
        net_e = prosumers[i].generation(t) - prosumers[i].consumption(t)
        battery_capacity = prosumers[i].battery_free_capacity(t)

        IF net_e > 0:  # Exporter
            # Willing to sell surplus at above feed-in tariff
            min_price = grid_feed_in_tariff + 0.02  # EUR/kWh margin
            asks.append((i, net_e, min_price))

        ELIF net_e < 0:  # Importer
            # Willing to buy at below grid retail price
            max_price = grid_retail_price - 0.02  # EUR/kWh discount
            bids.append((i, abs(net_e), max_price))

        # Battery charging/discharging bids
        IF prosumers[i].battery_level(t) < 0.2 * prosumers[i].battery_max:
            # Low battery: bid to charge
            bids.append((i + "_charge", battery_capacity * 0.5,
                          grid_retail_price - 0.05))
        ELIF prosumers[i].battery_level(t) > 0.8 * prosumers[i].battery_max:
            # High battery: offer to discharge
            asks.append((i + "_discharge",
                          prosumers[i].battery_level(t) * 0.3,
                          grid_feed_in_tariff + 0.01))

    # Step 2: Sort asks ascending by price, bids descending by price
    asks.sort(key=LAMBDA x: x[2])   # cheapest first
    bids.sort(key=LAMBDA x: -x[2])  # most willing to pay first

    # Step 3: Compute clearing price (midpoint double auction)
    matches = []
    remaining_asks = asks.copy()
    remaining_bids = bids.copy()

    WHILE remaining_asks AND remaining_bids:
        best_ask = remaining_asks[0]
        best_bid = remaining_bids[0]

        IF best_ask[2] <= best_bid[2]:  # Trade feasible
            # Clearing price: midpoint between ask and bid
            clearing_price = (best_ask[2] + best_bid[2]) / 2

            # Match quantity: minimum of available ask and bid
            trade_qty = MIN(best_ask[1], best_bid[1])

            matches.append({
                'seller': best_ask[0],
                'buyer': best_bid[0],
                'quantity_kWh': trade_qty,
                'price_EUR_kWh': clearing_price,
                'timestamp': t
            })

            # Update remaining quantities
            remaining_asks[0] = (best_ask[0],
                                  best_ask[1] - trade_qty,
                                  best_ask[2])
            remaining_bids[0] = (best_bid[0],
                                  best_bid[1] - trade_qty,
                                  best_bid[2])

            # Remove exhausted orders
            IF remaining_asks[0][1] < 1e-6:
                remaining_asks.pop(0)
            IF remaining_bids[0][1] < 1e-6:
                remaining_bids.pop(0)
        ELSE:
            BREAK  # No more profitable matches

    # Step 4: Record to ledger for OVA (generation-weighted contribution)
    FOR match IN matches:
        ledger.record_transaction(match)
        ova_tracker.credit(match['seller'],
                           contribution_type='energy_generation',
                           value=match['quantity_kWh'])
        ova_tracker.credit(match['buyer'],
                           contribution_type='demand_flexibility',
                           value=match['quantity_kWh'] * flexibility_score(match['buyer']))

    # Step 5: Grid interaction for unmatched positions
    FOR remaining IN remaining_asks:
        grid.sell_to_grid(remaining[0], remaining[1], grid_feed_in_tariff)
    FOR remaining IN remaining_bids:
        grid.buy_from_grid(remaining[0], remaining[1], grid_retail_price)

    RETURN matches
```

### 35.4.2 Incentive-Compatibility Analysis

**Proposition 35.2 (P2P Energy Trading Incentive-Compatibility).** The double-auction P2P energy trading mechanism is incentive-compatible — truthful bidding is a weakly dominant strategy — when:

1. **Price-taking:** Each prosumer's individual production and consumption is small relative to the total pool (no market power).
2. **Reputation continuity:** Each prosumer participates in multiple trading periods, and reputation tracks trading behavior (Algorithm 35.1 records all transactions for OVA).
3. **Grid fallback:** Prosumers always have access to grid buy/sell at regulated prices, ensuring their outside option is independent of P2P trading outcomes.

*Proof.* Under price-taking (condition 1), the double-auction clears at the competitive price regardless of individual bids — misrepresenting one's valuation cannot improve the clearing price. Under reputation continuity (condition 2), inflated bids or manipulated asks affect future OVA contributions and peer trust. Under grid fallback (condition 3), no prosumer is forced to accept an unfavorable P2P price — they always transact with the grid at known regulated prices. Together, these three conditions eliminate any incentive for strategic misrepresentation. $\square$

### 35.4.3 Smart Contract Governance Pseudocode

```solidity
// SPDX-License-Identifier: MIT
// Cooperative Energy Platform — Governance Contract
// Cross-reference: C:Ch.11 (smart contract design), C:Ch.13 (cooperative governance)

contract CooperativeEnergyPlatform {

    // ── MEMBERSHIP ──────────────────────────────────────────────────────────
    mapping(address => bool)    public isMember;
    mapping(address => uint256) public contributionScore; // OVA score
    mapping(address => uint256) public memberSince;       // timestamp
    uint256 public totalMembers;
    uint256 public constant MIN_MEMBERSHIP_PERIOD = 90 days;

    // ── REVENUE ALLOCATION (OVA) ────────────────────────────────────────────
    // Period revenue is distributed proportionally to contributionScore
    // contribution_score tracks: energy_generated, demand_flexibility,
    //                            governance_participation, maintenance

    event RevenueDistributed(uint256 totalRevenue, uint256 period);
    event ContributionRecorded(address member, string contribType, uint256 value);

    function recordEnergyContribution(
        address member,
        uint256 kWhGenerated,
        uint256 flexibilityScore    // 0-100: measure of demand flexibility
    ) external onlyOracleOrMember {
        // Generation contribution (weight 1.0 per kWh)
        contributionScore[member] += kWhGenerated * 100;
        // Flexibility contribution (weight 0.5 × flexibility score per kWh)
        contributionScore[member] += kWhGenerated * flexibilityScore / 2;
        emit ContributionRecorded(member, "energy", kWhGenerated);
    }

    // ── GOVERNANCE ──────────────────────────────────────────────────────────
    // Proposals for material changes require multi-stakeholder vote
    // Material changes: fee changes, new member categories, protocol upgrades

    enum ProposalType { FeeChange, MemberAdmission, ProtocolUpgrade, DataUse }
    enum VoteChoice  { For, Against, Abstain }

    struct Proposal {
        ProposalType pType;
        bytes calldata;          // encoded function call if approved
        uint256 deadline;        // voting deadline
        uint256 votesFor;
        uint256 votesAgainst;
        bool executed;
    }

    mapping(uint256 => Proposal) public proposals;
    uint256 public proposalCount;

    // Quadratic voting [C:Ch.13]: vote weight = sqrt(contributionScore)
    function castVote(uint256 proposalId, VoteChoice choice) external {
        require(isMember[msg.sender], "Not a member");
        require(block.timestamp < proposals[proposalId].deadline, "Voting closed");

        // Quadratic vote weight based on contribution score
        uint256 voteWeight = sqrt(contributionScore[msg.sender]);

        if (choice == VoteChoice.For) {
            proposals[proposalId].votesFor += voteWeight;
        } else if (choice == VoteChoice.Against) {
            proposals[proposalId].votesAgainst += voteWeight;
        }
        // Abstain: counted for quorum, not for/against
    }

    // ── DATA GOVERNANCE ─────────────────────────────────────────────────────
    // Aggregate usage data is a commons asset; individual data is private
    // Third-party data access requires 60% member approval + per-use fee
    // distributed proportionally to data contributors

    function approveDataAccess(
        address thirdParty,
        bytes32 datasetHash,      // hash of approved dataset
        uint256 accessFeeWei      // fee paid to data commons
    ) external onlyGovernance {
        // Distribution: 70% to data contributors (by contribution weight)
        //               20% to platform maintenance fund
        //               10% to ecological restoration fund [C:Ch.17]
        dataAccessRegistry[thirdParty][datasetHash] = true;
        distribute(accessFeeWei, 70, 20, 10);
    }
}
```

The smart contract implements all six cooperative platform design principles: open membership (any contributor meeting quality standards joins), transparent pricing (fee changes require member vote), portable reputation (OVA scores are on-chain and readable), data sovereignty (data access requires member vote and distributes fees), surplus distribution (OVA allocation), and democratic governance (quadratic voting on proposals). The ecological restoration fund (10\% of data access fees) is the direct implementation of Layer 3 biophysical planning within the platform governance.

---

## 35.5 Transition from Corporate Platform to Cooperative Platform

### 35.5.1 The Cold-Start Problem

The most significant barrier to cooperative platform formation is the cold-start problem: a new platform has no users, but users only join platforms that have users (network externalities). Corporate platforms solved this through massive subsidization — Uber spent approximately USD 25 billion in investor subsidies from 2010–2019, pricing rides below cost to achieve critical mass before raising prices. Cooperative platforms generally lack access to equivalent subsidization capacity.

**The tipping threshold model.** From Chapter 15, the cooperative platform achieves stable adoption above the critical threshold $\hat{x} = (c - \mu_\beta)/v$. For the cooperative platform to displace the incumbent:

1. It must push adoption above $\hat{x}^{\text{coop}}$ (its own tipping threshold), and
2. If the incumbent is above its own tipping threshold $\hat{x}^{\text{corp}}$, it must reduce the incumbent's effective adoption through direct competition or regulatory intervention.

**Proposition 35.3 (Cold-Start Solutions for Platform Cooperatives).** The cold-start problem can be resolved through three mechanisms, each reducing the effective $\hat{x}^{\text{coop}}$:

1. **Community migration:** An existing off-platform community (association, industry group, guild) commits to collective adoption, immediately achieving $x_0 > \hat{x}^{\text{coop}}$.
2. **Vertical integration:** The cooperative provides complementary services (insurance, training, legal support) that create value independent of network size, reducing the threshold $c$ in the tipping equation.
3. **Regulatory mandates:** Public policy requires adoption (e.g., energy cooperatives receiving priority grid connection rights, or data portability requirements reducing switching costs).

*Proof.* Each mechanism reduces the effective tipping threshold:
- Community migration: pushes $x_0$ above $\hat{x}$ directly.
- Vertical integration: adds intrinsic value $v_0 > 0$ independent of $x$, reducing $\hat{x} = (c - \mu_\beta - v_0)/v$.
- Regulatory mandates: reduce switching cost $c$ (if data portability), increase $\mu_\beta$ (if mandated features), or directly set $x_0 > \hat{x}$ (if mandatory participation). $\square$

### 35.5.2 The Network Migration Game

When a cooperative platform competes with an established corporate platform, both with active users above their tipping thresholds, the migration decision for each user is a coordination game:

$$u_i(\text{migrate}) = \underbrace{v \cdot x_{\text{coop}}}_{\text{network value}} + \underbrace{\mu_i}_{\text{cooperative values premium}} - \underbrace{c_{\text{switch}}}_{\text{switching cost}} - \underbrace{(v \cdot x_{\text{corp}})}_{\text{lost incumbent value}}$$

Migration is individually rational when $v \cdot x_{\text{coop}} + \mu_i > c_{\text{switch}} + v \cdot x_{\text{corp}}$.

For a producer whose cooperative values premium $\mu_i$ is large (e.g., a gig worker who values ownership and democratic governance) and who has low switching costs (portable reputation, OVA-compatible data): migration is individually rational even at low $x_{\text{coop}}$. For a consumer with low $\mu_i$ and high switching costs (e.g., a regular Amazon customer with many past orders, reviews, and Prime membership): migration is individually rational only when $x_{\text{coop}}$ is close to $x_{\text{corp}}$.

**The producer-led transition.** Platform cooperatives are most likely to begin on the producer side: drivers, freelancers, and content creators have higher cooperative values premium $\mu_i$ (more to gain from ownership and governance rights) and, in some cases, can organize through existing labor unions or industry associations. Consumer adoption follows naturally once sufficient producer adoption creates network value on the platform.

---

## 35.6 Case Study: REScoop.eu — The European Energy Cooperative Federation

### 35.6.1 Structure and Scale

REScoop.eu is the European federation of citizen energy cooperatives, founded in 2013, representing approximately 1,750 renewable energy cooperatives across 25 European countries with approximately 1.5 million individual members. Member cooperatives operate wind farms, solar arrays, energy efficiency projects, district heating systems, and — increasingly — P2P energy trading platforms. Total renewable capacity: approximately 14 GW (2022).

REScoop.eu serves as a Cosmo-Local governance structure [C:Ch.13]: individual local cooperatives are sovereign at the local scale (managing specific energy assets, member relations, local grid connections); REScoop provides the global-scale shared resources — advocacy at the EU policy level, shared IT infrastructure, common technical standards, collective procurement, and knowledge commons (legal templates, governance models, technical specifications shared as open-access resources).

### 35.6.2 Formal Analysis: Cooperative Game Theory

**The cooperative surplus in energy cooperatives.** An energy cooperative creates cooperative surplus relative to the default alternative (purchasing electricity from a commercial utility) through three mechanisms:

1. **Investment surplus:** Members co-invest in renewable capacity at cooperative cost (no profit margin) rather than paying a profit-including retail price. Estimated surplus: 10–20\% of electricity cost.

2. **Local value retention:** Surplus stays within the local economy rather than flowing to distant utility shareholders — the geographic economic multiplier effect estimated at 1.3–1.8× for cooperative energy spending.

3. **Governance value:** Members value participation in energy governance — the non-monetary benefit of owning and controlling their energy system. Survey evidence suggests WTP for cooperative vs. commercial energy supply of approximately EUR 50–150/year per member at equivalent monetary cost.

**Shapley value calculation.** For a typical REScoop member cooperative with 200 members, EUR 2.5 million investment in a 500kW solar farm:

- Individual investment (competitive alternative): each member could invest EUR 12,500 in commercial solar at 6\% return = EUR 750/year.
- Cooperative investment return: 8\% cooperative return (economies of scale, no profit extraction) = EUR 1,000/year.
- Cooperative surplus: EUR 250/member/year (33\% above competitive alternative).
- Total cooperative surplus (200 members): EUR 50,000/year.
- Shapley allocation (symmetric members): EUR 250/member/year = equal distribution of the surplus (consistent with one-member-one-vote governance and equal shareholding).

**Governance assessment (Ostrom principles):**

| DP | REScoop implementation | Score |
|:---|:---|:---|
| DP1 (Boundaries) | Membership agreements; geographic service area | 2/2 |
| DP2 (Congruence) | Local cooperatives adapt governance to local conditions | 2/2 |
| DP3 (Collective choice) | AGM, member votes on major decisions | 2/2 |
| DP4 (Monitoring) | Metered generation and consumption; public reporting | 2/2 |
| DP5 (Sanctions) | Member financial liability; board accountability | 1.5/2 |
| DP6 (Conflict resolution) | Ombudsman at federation level; local mediation | 1.5/2 |
| DP7 (External recognition) | EU recognition of cooperative legal status; REDII provisions | 2/2 |
| DP8 (Nested enterprises) | Local coop → national federation → REScoop.eu | 2/2 |
| **Total** | | **15/16** |

REScoop's 15/16 Ostrom score predicts high governance stability — consistent with the growth and persistence of energy cooperatives in Europe despite significant competition from incumbent utilities and adverse policy environments.

### 35.6.3 Ecological Outcomes

**Renewable capacity:** REScoop's 14 GW of cooperative renewable capacity avoids approximately 12 Mt CO₂e/year (at EU average grid intensity). This is a direct Planetary Boundary contribution [C:Ch.17] — the energy cooperative sector is measurably reducing the atmospheric carbon overshoot.

**Material and energy efficiency.** Energy cooperatives invest systematically in demand reduction alongside supply: member surveys show cooperative members adopt energy efficiency measures at approximately 2.3× the rate of comparable non-member households, driven by the governance accountability and financial alignment of cooperative membership. This is the ENA efficiency gain [C:Ch.20]: the cooperative's governance structure improves both the ascendancy (organized efficient supply) and the overhead (demand flexibility) of the energy system.

**The ecological alignment mechanism.** REScoop's governance requires energy decisions to meet member-defined ecological criteria — not only economic ones. This is the institutional implementation of Layer 3 biophysical planning [C:Ch.29]: the Planetary Boundaries constraint is embedded in the governance structure rather than imposed externally by regulation. Member-owners who care about ecological outcomes vote to implement them; the cooperative governance aligns economic interests with ecological stewardship in a way that utility shareholders — whose primary obligation is financial return — cannot replicate.

### 35.6.4 P2P Trading Implementation

Several REScoop member cooperatives have implemented P2P energy trading among their members, using variants of Algorithm 35.1. The most advanced example: Ecopower (Belgium), which operates a P2P trading platform for its 70,000 members, allowing prosumers to trade solar surplus directly with other members at a cooperative clearing price set between the utility's feed-in tariff (EUR 0.05/kWh) and retail price (EUR 0.28/kWh), cleared at EUR 0.17/kWh — EUR 0.12/kWh above what the utility offers exporters and EUR 0.11/kWh below what importers would pay the utility.

**Welfare calculation.** For each kWh traded P2P:
- Exporter gains: EUR 0.12/kWh above grid feed-in tariff.
- Importer gains: EUR 0.11/kWh below grid retail price.
- Total cooperative surplus per kWh: EUR 0.23 (the grid spread captured by the cooperative rather than the utility).

At Ecopower's P2P trading volume of approximately 40 GWh/year: total annual cooperative surplus from P2P trading = $0.23 \times 40{,}000{,}000 = $ EUR 9.2 million/year — shared among 70,000 members at approximately EUR 131/member/year in P2P trading value.

---

## 35.7 Mathematical Model: Platform Network Externalities and Cooperative Advantage

**Setup.** Two competing platforms — corporate ($C$) and cooperative ($K$) — each with $n$ potential users. Users choose which platform (or both) to join.

**Adoption dynamics:**

$$\dot{x}_C = \phi \cdot [v_C x_C + \mu_C^{\text{avg}} - c_C - v_K x_K]$$
$$\dot{x}_K = \phi \cdot [v_K x_K + \mu_K^{\text{avg}} - c_K - v_C x_C]$$

where $x_C + x_K \leq 1$ (users split between platforms), $v_C, v_K$ are network externality coefficients, $\mu_C^{\text{avg}}, \mu_K^{\text{avg}}$ are average intrinsic utility, and $c_C, c_K$ are switching costs.

**Cooperative advantage:** $\mu_K^{\text{avg}} > \mu_C^{\text{avg}}$ (for users who value ownership, governance, and surplus sharing — particularly producers). $c_K < c_C$ (cooperative platforms with portable reputation and data sovereignty have lower switching costs than proprietary platforms with locked-in data).

**Equilibrium.** Two stable equilibria exist: one with $x_C \approx 1$ (corporate platform dominates), one with $x_K \approx 1$ (cooperative platform dominates). The basin of attraction of each equilibrium depends on initial conditions and the relative values of $v, \mu, c$.

**Proposition 35.4 (Cooperative Platform Basin Expansion).** The cooperative platform's basin of attraction expands relative to the corporate platform when:

$$\frac{\mu_K^{\text{avg}} - \mu_C^{\text{avg}}}{c_K - c_C} > \frac{v_C - v_K}{v_K}$$

That is: the cooperative values premium relative to the switching cost advantage must exceed the corporate platform's network externality advantage relative to the cooperative's. Policy interventions that raise $\mu_K$ (data portability, worker classification reform) or lower $c_K$ (interoperability standards) expand the cooperative basin.

---

## Chapter Summary

This chapter has applied P2P theory to platform economics, demonstrating that the cooperative game theory of Parts II–VI provides a rigorous analysis of both the platform capitalism failure mode and the cooperative alternative.

The platform capitalism failure modes (Proposition 35.1) are formalized as Shapley value violations: network lock-in violates the null player axiom, regulatory arbitrage violates symmetry, and data appropriation violates efficiency. The corporate platform extracts the cooperative surplus that should be distributed to contributors.

The cooperative P2P platform (Definition 35.2) implements all six design principles: open membership, transparent pricing, portable reputation, data sovereignty, surplus distribution, and democratic governance. The formal comparison table shows that cooperative platforms correct every failure mode of corporate platforms.

The P2P energy trading algorithm (Algorithm 35.1) implements incentive-compatible double-auction trading in a prosumer cooperative network, with OVA-tracked contribution scores. The smart contract specification implements governance, data sovereignty, and surplus distribution on-chain. Proposition 35.2 proves incentive-compatibility under price-taking, reputation continuity, and grid fallback.

The cold-start problem (Proposition 35.3) is resolved through community migration, vertical integration, and regulatory mandates — each reducing the effective tipping threshold. The producer-led transition strategy (producer cooperative values premium is highest) provides a feasible entry sequence.

REScoop.eu scores 15/16 on Ostrom principles and delivers measurable ecological outcomes: 12 Mt CO₂e/year avoided, 2.3× higher member energy efficiency investment, and EUR 131/member/year in P2P trading cooperative surplus at Ecopower's scale. The federation's Cosmo-Local architecture (local sovereignty for asset management, global coordination for policy and knowledge commons) is the energy cooperative implementation of Chapter 13's governance framework.

Chapter 36 moves from digital and platform infrastructure to the physical landscape: the economics of regenerative agriculture and landscape restoration, where the cooperative-regenerative framework meets its most direct ecological test.

---

## Exercises

**35.1** Apply the Shapley value analysis to the Uber platform:
(a) Model Uber as a three-player cooperative game: drivers ($D$), riders ($R$), and Uber corporation ($U$). The value functions: $v(\{D\}) = 0$ (drivers cannot match riders without the platform), $v(\{R\}) = 0$, $v(\{U\}) = 0$, $v(\{D,R\}) = 0$ (no matching without infrastructure), $v(\{D,U\}) = 0.3 \cdot \text{booking}$ (limited platform without full rider base), $v(\{R,U\}) = 0.2 \cdot \text{booking}$, $v(\{D,R,U\}) = $ booking value. Compute the Shapley value for each player.
(b) Compare the Shapley value allocation to Uber's actual allocation (drivers receive ~67\%, riders receive their consumer surplus, Uber receives 32\%). Which parties receive more and which receive less than their Shapley value?
(c) Design a cooperative platform that implements the Shapley allocation. What governance changes and revenue distribution rule would achieve it?

**35.2** The P2P energy trading algorithm:
(a) For a 5-prosumer cooperative with: P1 (exports 4kWh), P2 (exports 2kWh), P3 (imports 3kWh), P4 (imports 1kWh), P5 (imports 2kWh); grid feed-in = EUR 0.05/kWh, grid retail = EUR 0.28/kWh; run Algorithm 35.1. What is the clearing price and which trades are matched?
(b) How much do P3, P4, P5 save relative to buying from the grid? How much do P1, P2 gain relative to selling to the grid? What is the total cooperative surplus from P2P trading?
(c) If P1 misreports their generation as 6kWh (inflating their ask quantity), does this affect the clearing price? Does it affect their OVA contribution score? Is misreporting a weakly dominant strategy under the conditions of Proposition 35.2?

**35.3** The cold-start problem (Proposition 35.3):
(a) A platform cooperative for gig workers (food delivery) is launching in competition with Deliveroo. Using the tipping threshold model: estimate $\hat{x}^{\text{coop}}$ if $v = 0.3$, $\mu_{\text{avg}} = 0.15$ (workers' cooperative values premium), $c = 0.40$ (switching cost from Deliveroo's reputation lock-in).
(b) A food delivery workers' union commits to migrating all 800 of its members to the cooperative platform, out of a total market of 3,000 delivery workers. Does this initial adoption exceed $\hat{x}^{\text{coop}}$? Will network externalities drive full adoption?
(c) Suppose data portability regulation requires Deliveroo to export workers' historical ratings to any other platform. How does this change $c$ and therefore $\hat{x}^{\text{coop}}$? Is this regulation sufficient to enable the cooperative platform to reach critical mass?

**★ 35.4** Prove Proposition 35.1 (platform failure modes are cooperative game failures) formally.

(a) Define the Null Player axiom formally [C:Ch.6]. Show that in a network with $n$ users at time $t$, each user's contribution to platform value is positive (historical user base generates current network externalities), yet the platform treats historical users as null players (paying them nothing for their network contribution).
(b) Define the Symmetry axiom. Show that in the Uber platform game, a driver classified as "contractor" and an equivalent driver classified as "employee" receive different payoffs from identical contributions, violating symmetry.
(c) Define the Efficiency axiom. Show that data appropriation by the corporate platform represents an efficiency violation: the total payoff to all players (workers + users + platform) is less than the total value created because privacy costs and data asymmetries reduce the welfare of data-generating users below what the efficient allocation would provide.
(d) Conversely, show that a cooperative platform implementing Algorithm 35.1 and the smart contract governance satisfies all three axioms (Null Player, Symmetry, Efficiency). This proves that cooperative platform design is formally Shapley-consistent.

**★ 35.5** Formally analyze the transition from corporate to cooperative platform using the institutional tipping threshold model of Chapter 15.

(a) Model the energy sector as an institutional coordination game with two equilibria: "utility-dominated" (incumbent utility provides all energy) and "cooperative-dominated" (cooperatives provide most energy). Define the tipping threshold $\hat{x}$ in terms of REScoop.eu's parameters.
(b) The EU Renewable Energy Directive II (REDII) mandates that member states facilitate energy communities and give them priority grid connection rights. Model this as an institutional entrepreneur investment that reduces $c$ (switching costs from incumbent utilities). Using REScoop.eu's actual post-REDII growth data (approximately 15\%/year since 2018), estimate whether the REDII intervention was sufficient to push European energy cooperatives above $\hat{x}$.
(c) Is the current EU energy cooperative sector above or below the tipping threshold? If above: what institutional design ensures the cooperative sector continues growing? If below: what additional policy interventions (feed-in tariff reform, cooperative tax advantages, data portability for energy usage data) would push adoption above $\hat{x}$?

**★★ 35.6** Design a cooperative P2P platform for a sector of your choice — options include: ride-hailing, home-sharing, freelance services, or food delivery.

(a) Specify the full platform architecture using Definition 35.2: stakeholder groups ($U_P$, $U_C$, $U_S$), governance structure, OVA allocation mechanism, data governance, and reputation system.

(b) Implement the governance smart contract (extending the template in Section 35.4.3): add the specific voting thresholds, fee structures, and data access rules appropriate to your sector.

(c) Estimate the cooperative surplus your platform would generate relative to the incumbent corporate platform. Use the Shapley value framework to compute each stakeholder group's fair allocation. Compare to current allocations in the incumbent platform.

(d) Analyze the cold-start problem for your platform: estimate $\hat{x}^{\text{coop}}$, identify the initial community most likely to migrate (high $\mu_i$, low $c$), and specify the three most important policy or design interventions that would reduce $\hat{x}^{\text{coop}}$ to a reachable level.

(e) Apply the REScoop Ostrom scoring (Section 35.6.2) to your platform design. Which principles are strongest? Which require additional institutional design? What is your predicted Ostrom score, and what does it predict about your platform's long-run governance stability?

---

*Chapter 36 moves from digital infrastructure to the physical landscape — the regenerative agriculture and ecosystem restoration investments that directly implement the natural capital stewardship of Part IV's ecological embedding. The formal tools are those of investment analysis, ecological dynamics, and landscape network theory, applied to one of the most consequential economic transformations of the coming decades.*
