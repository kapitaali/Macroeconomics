# Part IX: Capstone — Building a New Economic Science

The formal framework is now complete. Forty-two chapters have built the architecture of cooperative-regenerative economics from its game-theoretic foundations through network governance, ecological embedding, monetary alternatives, unified synthesis, stability analysis, post-growth economics, inequality dynamics, digital commons, six empirical applications, and transition theory. Part VII confirmed that this framework is not merely theoretical — it describes and explains institutions that already exist and operate successfully across multiple continents and sectors. Part VIII provided the transition theory and policy design needed to move from the margin to the mainstream.

Part IX does not summarize. Summary is the task of Chapter 45. Instead, this final Part looks forward: at the open questions the framework raises without resolving (Chapter 43), at the design problem that synthesizes the entire toolkit into applied practice (Chapter 44), and at the synthesis that situates the cooperative-regenerative framework within the longer arc of economic thought and extends an invitation to the reader to become a contributor rather than merely a student (Chapter 45).

A scientific framework is defined as much by what it does not yet know as by what it does. The most important questions in cooperative-regenerative economics are not the ones this book has answered — those are the beginning, not the end. The open questions of Chapter 43 are the research agenda for the next decade, the problems that will define the field, and the places where a graduate student, a field researcher, or a practitioner reading this book might make a decisive contribution.

The capstone project of Chapter 44 is the book's final formal exercise — not an exercise in answering questions but in asking the right ones about a specific real-world locality, and then applying the complete toolkit to design an economy that would serve that community's flourishing better than the one it currently has. It is, in miniature, what the entire trilogy has been building toward: not a description of economics as it is, but a design for economics as it could be.

---

# Chapter 43: Open Questions and Future Research Directions

> *"The most exciting phrase to hear in science, the one that heralds new discoveries, is not 'Eureka!' but 'That's funny…'"*
> — attributed to Isaac Asimov

> *"A theory that explains everything explains nothing. A good theory generates better questions than it started with."*
> — Karl Popper, paraphrased

## Learning Objectives

By the end of this chapter, you should be able to:

1. Identify the most important theoretical gaps in the cooperative-regenerative framework — the formal claims that are plausible but unproven, the models that are promising but incomplete, and the assumptions that have been made for tractability but require relaxation.
2. Specify empirical gaps — the data that does not yet exist, the natural experiments that have not been exploited, and the measurement problems that prevent rigorous testing of the framework's core predictions.
3. Map methodological frontiers — the tools from adjacent fields (machine learning, computational complexity, evolutionary biology, Earth system science) that could extend the cooperative-regenerative framework into currently intractable territory.
4. Articulate fifteen prioritized open research questions, each formulated with sufficient precision to guide empirical or theoretical research, and identify the methodological approach most likely to be productive for each.
5. Evaluate the Santa Fe Institute model of complexity-based interdisciplinary research — what it has contributed to economic thinking, where it has fallen short, and what an economics-specific equivalent would need to do differently.

---

## 43.1 The Nature of Theoretical Incompleteness

A formal framework is incomplete in two distinct ways. The first is gaps in the argument — places where the framework makes assertions that are reasonable but not formally proved. The second is boundaries of the framework — places where the framework's assumptions break down and a different or extended framework is needed. Both kinds of incompleteness are identified here, because both define research opportunities.

The cooperative-regenerative framework as presented in this book is deliberately mathematically rigorous about the claims it makes formally, and honest about what the formal results do and do not imply. The Cooperative Stewardship Theorem (Theorem 29.2) proves that the CRE yields higher IPI than the CE under three conditions — but it does not prove that those conditions are jointly satisfiable in all real economies, nor that the transition from CE to CRE can be achieved without welfare loss during the transition. The Transition Tipping Point Theorem (Theorem 40.2) derives $\hat{x}^{\text{CRE}}$ under specific model assumptions — but the calibration of the underlying parameters ($\theta$, $v_{\text{CRE}}$, $\phi_{\text{landscape}}$) relies on rough estimates rather than rigorous econometric identification.

These are not failures of the framework — they are invitations. Each gap is a research program.

---

## 43.2 Theoretical Gaps

### Gap T1: Cooperative Formation under Incomplete Information

The cooperative game theory of Chapter 6 assumes that all agents know the characteristic function $v(S)$ — the value of every coalition. In practice, agents do not know what value different coalitions would create: they do not know each other's productive capacities, preferences, or reliability. Under incomplete information about coalition values, the standard cooperative solution concepts (core, Shapley value, nucleolus) are not well-defined, and the cooperative formation process itself is a subject of strategic manipulation.

**The open problem.** Under what conditions do cooperative institutions form when agents have private information about their own productive contribution to potential coalitions? When is truthful revelation of private capacities incentive-compatible? Can the Shapley value be approximated through iterative negotiation processes that do not require full information?

**Promising approaches.** Bayesian cooperative game theory [Myerson, 1984]; mechanism design for cooperative formation [Maskin, 2008]; agent-based simulation of cooperative formation under incomplete information [C:Ch.10 extensions]. The connection to the data Shapley problem (Chapter 39, Algorithm 39.1) is direct — approximating Shapley values from partial information is both a data governance problem and a cooperative formation problem.

**Why it matters.** The principal practical barrier to cooperative enterprise formation is not that the cooperative surplus $\sigma_v$ is too small to be worth organizing — it is that potential cooperative members do not know in advance what the surplus would be. Formalizing this problem would provide the foundation for cooperative formation support programs that are currently designed by intuition rather than by theory.

### Gap T2: Cooperative Dynamics at Scale

The Cooperative Resilience Theorem (Chapter 30, Theorem 30.2) proves that the CRE's basin of attraction is larger than the CE's — but this is a local result, valid near the equilibrium. The global dynamics — the trajectory of the economy as it moves from the CE basin toward the CRE basin during transition — are not fully characterized.

**The open problem.** What is the full phase portrait of the coupled social-ecological-monetary dynamical system (equations 29.12–29.20)? Are there limit cycles (recurring instabilities)? Are there intermediate attractors that trap the transition before reaching the full CRE? Under what conditions is the transition trajectory monotone (welfare continuously improving) versus non-monotone (welfare temporarily falling before recovering)?

**Promising approaches.** Numerical phase portrait analysis of the full nine-equation system; continuation methods (AUTO, MATCONT) for bifurcation analysis; stochastic stability theory for economies with multiple attractors [Foster and Young, 1990; Young, 1993].

**Why it matters.** Non-monotone transition trajectories — in which welfare temporarily falls before the cooperative-regenerative institutions generate full benefit — are politically dangerous (the Legitimacy Trap of Chapter 41). Characterizing the conditions that produce non-monotone trajectories would allow transition policy to be designed specifically to avoid them.

### Gap T3: The Formal Connection between Ecological and Monetary Instability

Chapters 17–22 and Chapters 23–28 developed ecological dynamics and monetary dynamics as parallel frameworks, connected through the SFC-N accounting framework of Chapter 18. But the formal mathematical connection between ecological instability (approaching critical thresholds) and monetary instability (Minsky dynamics) has not been derived.

**The open problem.** Is there a formal coupling between natural capital depletion dynamics and financial instability? Does ecological overshoot increase the probability of Minsky events? Does monetary instability accelerate ecological depletion? If both answers are yes, is there a unified instability theorem that characterizes the joint social-ecological-monetary crisis?

**Promising approaches.** Coupled SFC-ecological dynamical systems; historical analysis of the correlation between financial crises and ecological stress indicators (2008: drought stress in food markets coincided with financial fragility); network analysis of the financial-ecological interdependency using the bipartite network structure of firms and ecological resources.

**Why it matters.** If ecological and monetary instability are coupled, then addressing only one — as current ecological economics and post-Keynesian monetary economics typically do in isolation — leaves the other unchecked. A unified instability theorem would provide the formal basis for integrated monetary-ecological policy that neither field currently offers.

### Gap T4: The Optimal Demurrage Rate under Ecological Uncertainty

Theorem 27.2 derived the optimal demurrage rate as $\delta^* = i^{\text{natural}} - \mathcal{R}^{\text{MSY}}$. But $\mathcal{R}^{\text{MSY}}$ is uncertain — natural capital regeneration rates are not known with precision, vary spatially and temporally, and are subject to regime shifts (Chapter 19) that can dramatically alter the regeneration function. Under uncertainty about $\mathcal{R}^{\text{MSY}}$, the optimal demurrage rate under expected value may differ substantially from the optimal demurrage rate under robustness (minimizing maximum regret) or under risk aversion.

**The open problem.** Under uncertainty about the natural capital regeneration rate $\mathcal{R}^{\text{MSY}}(N, t)$, what is the optimal demurrage rate policy — and how should it be updated in response to new ecological information? Is the optimal policy a fixed rate, a time-varying rate, or a state-contingent rule (demurrage rate responding to observed natural capital dynamics)?

**Promising approaches.** Robust monetary policy under ecological uncertainty (adapting the Brainard principle of monetary conservatism to the demurrage design problem); Bayesian updating of the optimal demurrage rate using the ESP ecological monitoring system of Chapter 20; Pontryagin optimal control for the demurrage-as-ecological-policy problem.

### Gap T5: Cooperative Platforms at Massive Scale

The cooperative platform theory of Chapter 35 was developed for platforms with thousands to tens of thousands of members. Many of the most consequential platforms operate at billions-of-users scale — Meta, Google, Amazon. The governance mechanisms appropriate for a 10,000-member data cooperative (Chapter 39) may not scale to a 3-billion-user social media platform.

**The open problem.** What are the formal governance mechanisms that can implement cooperative principles — democratic decision-making, Shapley value allocation, data sovereignty — at billion-user scale? At what scale does the one-member-one-vote principle become computationally infeasible, and what governance approximations preserve its essential properties?

**Promising approaches.** Liquid democracy and proxy voting at scale [Paulin, 2020]; random stratified governance samples (citizens' assemblies as governance mechanisms for large platforms); AI-assisted governance decision synthesis (aggregating millions of preferences into coherent policy without losing the democratic character).

### Gap T6: Cooperative Predistribution Dynamics

Chapter 41b formalised the Cooperative Predistribution Equilibrium — combining employment guarantee, working time reduction, and Shapley-consistent care valuation — and proved that it achieves strictly lower primary income inequality than the competitive equilibrium (Theorem 41b.4). But the chapter established this as a comparative statics result: it characterises the CPE's steady state without fully analysing the dynamics of reaching it. Several questions about the dynamics remain open.

**The open problems.** (i) The dynamic care Shapley game: Chapter 41b (Proposition 41b.1) established that care workers have positive Shapley values in the static cooperative game. The dynamic extension — in which care investment in period $t$ generates human capital that increases the value of coalitions in period $t+k$ — is significantly more complex. What is the optimal care investment trajectory in the dynamic cooperative game, and does it converge to the static Shapley allocation as a special case? (ii) The CPE transition path: Chapter 41b's Theorem 41b.4 characterises the CPE as a steady state but does not prove that the economy converges to it from a competitive equilibrium starting point. The coupled dynamics of the buffer stock ($L^{\text{EG}}$), working hours ($H$), labour Gini ($G_L$), and care OVA ($\ell_c^{\text{OVA}}$) form a four-dimensional nonlinear system whose phase portrait has not been fully characterised (Exercise ★★ 41b.6 poses this as an open problem). (iii) The gender dynamics of care redistribution: the 5Rs framework specifies that care must be redistributed toward gender equity, but the formal model of intra-household care bargaining was only sketched.

**Promising approaches.** Pontryagin optimal control for the dynamic care game [C:Appendix A.5]; stochastic stability theory [Foster and Young, 1990] for convergence to the CPE; intra-household bargaining models in the cooperative game tradition [Browning et al., 2014]; empirical calibration using the Nordic care cooperative datasets.

---

## 43.3 Empirical Gaps

### Gap E1: Longitudinal Data on Cooperative Enterprise Clusters

The empirical evidence on cooperative enterprise performance (Chapters 30, 34) draws primarily on cross-sectional or short-panel studies. The 50-year Mondragon comparison is exceptional in its length; most cooperative studies have panels of 5–10 years. This makes it difficult to distinguish the cooperative's short-run resilience advantage (well-documented) from its long-run productivity dynamics — which theory predicts should be positive (cooperative learning and investment) but empirical studies have not confirmed at long horizons.

**The research need.** Longitudinal panel databases of cooperative enterprises matching individual firm data over 20+ years, across multiple countries, with standardized accounting for cooperative-specific features (member capital accounts, patronage dividends, governance participation, internal labor market transfers). No such database currently exists at the needed scale. The most promising source: European cooperative federations (Cooperatives Europe, ICA Europe) could establish common reporting standards; national registries could be linked.

### Gap E2: Measuring Commons Value

The digital commons valuation methodology of Chapter 33 (Brynjolfsson et al. willingness-to-accept; Hoffmann et al. replacement cost) provides reasonable estimates but with substantial uncertainty. The natural capital valuation of Chapter 17 (ecosystem service pricing; benefits transfer) faces similar challenges. In both cases, the core difficulty is that non-market goods — whose value is the entire point — resist the revealed preference methods designed for market goods.

**The research need.** Better stated preference methods for digital commons value (combining WTA surveys with conjoint analysis and behavioral experiments); natural capital stock accounting standardized across countries (SEEA EA implementation by national statistical offices — currently at approximately 25\% global coverage); longitudinal tracking of commons value changes as governance improves or degrades (connecting changes in Ostrom principle compliance to changes in measured commons value).

### Gap E3: Endogeneity in Natural Capital Accounting

The natural capital shadow prices in the SFC-N framework (Chapter 18) are derived from the cooperative optimal control problem (Chapter 29, Proposition 29.3). But in practice, shadow prices are estimated from observed market behavior — which itself reflects current institutional distortions. The market price of carbon is not the shadow price of atmospheric carbon; it is the current cap-and-trade permit price, which reflects political compromises rather than the true social cost of carbon.

**The research need.** Identification strategies for natural capital shadow prices that are not contaminated by existing institutional distortions. Promising approaches: hedonic pricing of ecosystem services in housing markets (clean identification through discontinuity designs at habitat boundaries); experimental ecosystem service valuation (randomized conservation payment experiments); satellite-based monitoring that enables difference-in-differences designs exploiting exogenous variation in natural capital stocks.

### Gap E4: The Cooperative Network Externality

The Transition Tipping Point Theorem requires an estimate of $v_{\text{CRE}}$ — the network externality coefficient of cooperative-regenerative institutions. Chapter 40's calibration ($v_{\text{CRE}} \approx 0.35$) is based on the Quebec cooperative density regression (Chapter 32) and the Emilia-Romagna data (Chapter 29). But these estimates conflate the direct effect of cooperative institutions with the network externality specifically — the additional benefit each new cooperative enterprise brings to existing ones through supply chain integration, shared governance standards, and mutual credit networks.

**The research need.** Identification of the cooperative network externality separately from the direct cooperative enterprise effect. Instrumental variable approaches (exogenous variation in cooperative density — e.g., post-WWII cooperative formation driven by regional reconstruction policies — used to identify the spillover effect on surrounding firms). Input-output analysis to quantify supply chain integration across cooperative and conventional firms. Mutual credit network analysis to measure the increase in clearing efficiency as cooperative density rises.

### Gap E5: Measuring Deliberative Convergence and SIW Implementation

Chapter 41c's Theorem 41c.2 (Deliberative Convergence) predicts that citizens' assemblies operating under informational completeness and justificatory constraints converge in probability toward a Condorcet winner in the MPD welfare weight space. The French and Irish assembly evidence is consistent with this prediction, but neither assembly was designed to test it — preference data was collected only partially, the convergence dynamics were not tracked systematically, and the sessions were too short to observe the full trajectory. Similarly, the SIW accounting framework (Definition 41c.4) is conceptually specified but empirically thin: only 8 countries publish anything close to integrated SIW accounts, and none has run the framework for long enough to test Theorem 41c.3's divergence prediction against real-time data.

**The research needs.** (i) A systematic, pre-registered assembly study: a citizens' assembly on wellbeing priorities in which participant welfare weights ($\mathbf{w}_i$) are elicited at the beginning and end of each deliberative session, using a conjoint analysis instrument calibrated to the MPD dimensions. The convergence dynamics — whether preferences converge, at what rate, toward what distribution, and under what conditions convergence fails — can then be estimated directly. (ii) A comparative SIW implementation study: a panel of countries implementing different components of the SIW framework (New Zealand, Scotland, Wales, Iceland, Finland) with standardised measurement protocols, enabling cross-country identification of the SIW-GDP divergence predicted by Theorem 41c.3. The OECD Framework for Measuring Wellbeing and Progress provides a partial institutional infrastructure; what is missing is the time-series depth and the ecological stocks dimension (SEEA-EA) needed to test the theorem.

**Promising data sources.** The newDemocracy Foundation's systematic archive of mini-publics; OECD Better Life Index longitudinal data; Scotland's National Performance Framework (longest-running national wellbeing dashboard, since 2007); the UN Statistics Division's SEEA-EA pilot implementation programme.

---

## 43.4 Methodological Frontiers

### Frontier M1: Machine Learning and Economic Complexity

The cooperative-regenerative framework's models are highly nonlinear, high-dimensional, and analytically intractable at scale. Machine learning methods — particularly deep reinforcement learning, graph neural networks, and physics-informed neural networks — offer new tools for:

- **Agent-based simulation at scale:** Current ABMs [C:Ch.10] can simulate hundreds to thousands of agents. Deep reinforcement learning enables millions of heterogeneous agents with sophisticated behavioral repertoires, enabling simulation of economy-wide transitions with granularity currently impossible.

- **Cooperative game learning:** Machine learning algorithms that approximate Shapley values for cooperative games with thousands of players — enabling practical implementation of OVA in large-scale data cooperatives and cooperative enterprises without the $O(2^n)$ computational barrier.

- **Policy learning at scale:** Multi-agent reinforcement learning as a formal tool for policy design — discovering policy sequences that navigate the transition traps of Chapter 41 through computational search rather than analytical derivation.

**The frontier.** Physics-informed machine learning (PIML) — deep neural networks that are constrained to satisfy the SFC accounting identities and ecological conservation laws — would combine the computational power of ML with the formal discipline of the SFC-N framework. This would enable large-scale simulation of cooperative-regenerative transitions with both economic and ecological consistency.

### Frontier M2: Evolutionary Economics and Institutional Dynamics

The institutional dynamics of Chapter 15 (replicator dynamics, crystallization threshold, institutional emergence) are first steps in an evolutionary approach to cooperative-regenerative economics. A fully evolutionary framework would model:

- **Variation:** The generation of new institutional designs through experimentation, imitation, and recombination (the multi-armed bandit of Chapter 42 captures part of this).
- **Selection:** The mechanisms through which some institutions survive and spread while others fail (the MLP framework of Chapter 40 describes this at a qualitative level).
- **Inheritance:** The transmission of institutional designs across generations, jurisdictions, and cultural contexts (the policy diffusion model of Chapter 40 addresses this partially).

**The frontier.** Evolutionary game theory with cultural evolution [Henrich, 2016]; institutional evolutionary economics [Nelson and Winter, 1982; revisited with formal evolutionary models]; evolutionary graph theory applied to institutional diffusion on social networks. The key gap: evolutionary economics has rich verbal theory but limited formal models; the cooperative-regenerative framework provides the formal structures that evolutionary models could be built around.

### Frontier M3: Earth System Science and Planetary Boundaries

The Planetary Boundaries framework (Chapter 17) provides a nine-dimensional constraint set for the cooperative-regenerative economy. But the Planetary Boundaries science continues to evolve: the boundaries themselves are uncertain (Chapter 17, Table 17.1 acknowledges this), the interactions between boundaries are incompletely understood, and the tipping point dynamics of Earth system transitions are an active research frontier [Lenton et al., 2019].

**The frontier.** Integrating the latest Earth system model outputs (CMIP6 and beyond) with the economic SFC-N framework to create a fully coupled social-ecological-economic model. The International Earth System Model Intercomparison Project (ESMIP) and the Earth System Governance project are developing the institutional infrastructure for this integration; the cooperative-regenerative economic framework is the natural economic complement.

---

## 43.5 The Research Agenda: Eighteen Prioritized Questions

Drawing from the gaps and frontiers above — including the three new gaps opened by Chapters 41b and 41c (Gap T6, Gap E5) — we identify eighteen prioritized research questions for the cooperative-regenerative economics research program over the next decade.

**Priority A: Foundational theory (most urgent)**

**Q1.** *Cooperative formation under incomplete information:* Under what conditions do cooperative institutions form endogenously when agents have private information about their productive contributions to potential coalitions? Methodology: Bayesian cooperative game theory + mechanism design.

**Q2.** *Joint social-ecological-monetary instability:* Is there a formal coupling between natural capital depletion dynamics and Minsky financial instability that generates joint crises? Methodology: Coupled dynamical systems analysis; historical econometrics.

**Q3.** *Optimal demurrage under ecological uncertainty:* What is the robustly optimal demurrage rate policy under uncertainty about natural capital regeneration rates? Methodology: Robust control theory; Bayesian optimal policy.

**Q4.** *Transition trajectory characterization:* What is the full phase portrait of the cooperative-regenerative transition dynamical system? Are there limit cycles, intermediate attractors, or non-monotone welfare trajectories? Methodology: Continuation methods (AUTO); stochastic stability analysis.

**Q5.** *Cooperative governance at billion-user scale:* What governance mechanisms implement democratic cooperative principles at platform scales currently dominated by surveillance capitalism? Methodology: Mechanism design for large populations; distributed consensus protocols.

**Q6.** *Cooperative predistribution dynamics:* What is the optimal care investment trajectory in the dynamic Shapley game, and do the coupled CPE dynamics (buffer stock, working hours, care OVA, labour Gini) converge asymptotically to the CPE steady state from competitive equilibrium initial conditions? [C:Ch.41b, Gap T6] Methodology: Pontryagin optimal control; stochastic stability; intra-household bargaining theory.

**Priority B: Empirical foundations (most needed)**

**Q7.** *Longitudinal cooperative enterprise performance:* Does the cooperative productivity and resilience advantage persist over 20–50-year horizons, or does it converge to conventional firm performance as cooperatives age? Methodology: Long-panel econometrics; synthetic control methods.

**Q8.** *Commons value measurement:* What stated preference and behavioral methods produce reliable, replicable estimates of digital commons and natural capital value? Methodology: Conjoint analysis; experimental economics; large-scale WTA surveys.

**Q9.** *Cooperative network externality identification:* What is the causal magnitude of the cooperative network externality — the additional benefit each new cooperative brings to existing ones through supply chain and governance spillovers? Methodology: Instrumental variables using exogenous cooperative density variation; input-output network analysis.

**Q10.** *Demurrage adoption dynamics:* Do communities that adopt demurrage currencies exhibit the predicted velocity increase, reduced hoarding, and counter-cyclical stabilization? Large-scale quasi-experimental evidence from digital demurrage implementations? Methodology: DiD exploiting variation in demurrage adoption; velocity measurement from digital transaction records.

**Q11.** *Ecological-monetary feedback:* Is there an empirical correlation between natural capital depletion and subsequent financial instability, controlling for standard macroeconomic variables? Methodology: Panel econometrics across countries with natural capital stock data; Granger causality testing.

**Q12.** *Deliberative convergence and SIW implementation:* Do citizens' assemblies on wellbeing priorities exhibit the convergence dynamics predicted by Theorem 41c.2, and do the countries implementing SIW accounting exhibit the GDP-SIW divergence predicted by Theorem 41c.3? [C:Ch.41c, Gap E5] Methodology: Pre-registered assembly study with panel elicitation of welfare weights; comparative SIW panel across New Zealand, Scotland, Wales, Iceland, Finland.

**Priority C: Applied design (most actionable)**

**Q13.** *Optimal CR-SEZ design:* What specific institutional configurations of the Cooperative-Regenerative Special Economic Zone generate the largest welfare improvements? Methodology: Multi-armed bandit experimentation (Chapter 42, Algorithm 42.1) across CR-SEZ variants.

**Q14.** *Transition sequencing optimization:* Can the transition path be formally optimized using reinforcement learning to discover sequences that maximize IPI while respecting the three transition trap constraints? Methodology: Deep RL applied to the unified model (Chapter 29).

**Q15.** *Planetary Ledger implementation:* What technological, institutional, and governance architecture is needed to implement the Planetary Ledger (Chapter 20) at sufficient scale and precision for the GTA framework (Chapter 17) to be operationally binding? Methodology: Technology assessment; comparative institutional analysis of Earth system governance architectures.

**Q16.** *Cooperative monetary transition dynamics:* What is the welfare trajectory during the 20-year transition from debt-based to sovereign money + mutual credit + demurrage hybrid (Chapter 28)? Are there significant transition costs and how are they distributed? Methodology: Agent-based SFC modeling; historical analogy from Iceland and Denmark case study.

**Q17.** *Data cooperative scaling:* At what membership scale do data cooperatives generate sufficient per-member dividends to achieve self-sustaining adoption without regulatory mandate? Methodology: Shapley value simulation under different adoption scenarios; tipping threshold analysis calibrated to sector-specific data value functions.

**Q18.** *Democratic planning impact:* Does the adoption of wellbeing budgeting with democratically determined welfare weight vectors produce measurably higher IPI than GDP-centric budgeting, controlling for country fixed effects and income level? What institutional features of the citizens' assembly process most strongly predict convergence quality? [C:Ch.41c, Theorem 41c.5] Methodology: Comparative interrupted time-series design exploiting staggered national wellbeing budget adoption; structural equation modelling of assembly process features and convergence outcomes.

---

## 43.6 Adjacent Disciplines: The Interdisciplinary Map

The cooperative-regenerative framework is inherently interdisciplinary — it draws from economics, ecology, game theory, network science, political science, and thermodynamics. But it has not yet fully integrated with several adjacent fields that could significantly extend it.

**Complexity Biology.** Biological evolution has solved problems of cooperation, information processing, and resilience over billions of years of selection pressure. The formal models of Chapter 7 (replicator dynamics, ESS) and Chapter 5 (complex adaptive systems) are borrowed from evolutionary biology — but only at an introductory level. The full apparatus of evolutionary biology — multilevel selection theory, major evolutionary transitions, the evolution of cooperation in structured populations — offers much more. The evolution of eusociality in insects (where individual workers sacrifice reproductive fitness for colony success) provides a model for understanding how cooperative institutions emerge and stabilize that economics has not yet absorbed.

**Earth System Science.** The Planetary Boundaries science is developing rapidly: new boundaries (aerosol loading, novel entities), revised estimates of existing boundaries (the land-system change boundary was tightened in 2023), and better understanding of boundary interactions (the freshwater-climate nexus, the biodiversity-carbon nexus). The cooperative-regenerative framework's ecological embedding needs to track these developments, updating its constraint set and coupling mechanisms as Earth system science advances.

**Political Philosophy.** The normative grounding of the cooperative-regenerative framework — why should we care about IPI rather than GDP, why is the Stewardship Condition a moral obligation rather than merely a practical constraint — is asserted rather than argued in this book. Political philosophy (particularly in the capabilities tradition of Nussbaum and Sen, and the commons theory of Ostrom's extended work) provides the normative foundations. A serious engagement with political philosophy would ground the framework's normative commitments in the depth they deserve.

**Computer Science: Distributed Systems and Cryptography.** The digital infrastructure of the cooperative-regenerative economy — the Planetary Ledger (Chapter 20), the smart contract cooperative governance (Chapter 35), the blockchain-verified PES systems (Chapter 36), the data cooperative governance (Chapter 39) — all require distributed systems architecture that can provide security, verifiability, and trustlessness at the scale needed for the cooperative-regenerative economy to operate. The formal tools of distributed computing (Byzantine fault tolerance, consensus protocols, cryptographic commitment schemes) are required infrastructure for the digital cooperative-regenerative economy, and their integration with the economic framework is at an early stage.

---

## 43.7 Case Study: The Santa Fe Institute as a Model

### 43.7.1 What It Has Contributed

The Santa Fe Institute (SFI), founded in 1984, has been the world's leading center for complexity science applied to economics. Its contributions to the themes of this book include: agent-based modeling of economic systems (Axelrod's cooperation work; Epstein and Axtell's Sugarscape); network economics (Watts and Strogatz; Barabási-Albert); evolutionary game theory (Kauffman's NK fitness landscapes); and complexity-theoretic economics (Arthur's increasing returns, path dependence, and the economy as an evolving complex system).

The SFI model — small, interdisciplinary, bringing economists into sustained dialogue with physicists, biologists, and computer scientists — has demonstrably accelerated the development of complexity economics. Its influence is visible in every chapter of this book.

### 43.7.2 What Remains

The SFI model has limitations that a cooperative-regenerative economics research institute would need to address differently:

**The normative gap.** SFI-style complexity economics is predominantly descriptive — it models how economies work as complex adaptive systems. The cooperative-regenerative framework is explicitly normative: it argues for a specific institutional design. SFI-style research informs the descriptive modeling but does not provide the normative grounding.

**The policy disconnect.** SFI's research culture prizes intellectual novelty; policy application is secondary. The cooperative-regenerative economics research agenda requires deep engagement with practitioners — cooperative enterprises, commons governance institutions, municipal governments — whose problems cannot be addressed from a distance.

**The field engagement gap.** The most important empirical data for cooperative-regenerative economics is generated in the field — by cooperative enterprises, commons governance experiments, complementary currency networks. SFI's research model is not designed for long-term field engagement of this kind. The research institute that could advance cooperative-regenerative economics would combine SFI's interdisciplinary rigor with the field-engaged methodology of participatory action research.

**What an Economics of Cooperation Institute Would Look Like.** A dedicated cooperative-regenerative economics research institute would: (i) maintain a large-panel longitudinal database of cooperative enterprise performance across countries; (ii) run a rotating residency program bringing practitioners (cooperative managers, municipal officials, commons governance coordinators) into dialogue with researchers; (iii) maintain an open-source modeling commons — shared ABM frameworks, SFC-N model libraries, network analysis tools — freely available to researchers globally; (iv) publish an annual State of the Cooperative-Regenerative Economy report, tracking the fifteen research questions above against available evidence.

---

## Chapter Summary

This chapter has identified the theoretical gaps, empirical needs, methodological frontiers, and interdisciplinary connections that define the research agenda for cooperative-regenerative economics in the decade ahead.

Five theoretical gaps stand out as most fundamental: cooperative formation under incomplete information (Gap T1), cooperative dynamics at scale (Gap T2), the formal ecological-monetary instability connection (Gap T3), the optimal demurrage rate under ecological uncertainty (Gap T4), and cooperative governance at massive platform scale (Gap T5). Each is a well-posed research problem with identified methodological approaches; none is intractable with current tools.

Four empirical gaps require new data infrastructure: longitudinal cooperative enterprise databases (Gap E1), improved commons value measurement (Gap E2), causal identification of ecological shadow prices (Gap E3), and identification of the cooperative network externality (Gap E4). Addressing these gaps requires coordination across national statistical offices, cooperative federations, and research institutions — a collaborative infrastructure that does not yet exist at the needed scale.

Three methodological frontiers offer the highest potential for extending the framework: physics-informed machine learning for large-scale simulation (Frontier M1), evolutionary economics with formal evolutionary models (Frontier M2), and fully coupled social-ecological-economic modeling integrating Earth system science outputs (Frontier M3).

The fifteen prioritized research questions span foundational theory (Q1–Q5), empirical foundations (Q6–Q10), and applied design (Q11–Q15) — a decade's work for a research community of the scale and ambition that cooperative-regenerative economics requires.

Chapter 44 provides the synthesis exercise: a structured capstone project that applies the complete toolkit from this book and its two predecessors to the concrete design of a local cooperative-regenerative economy — demonstrating, in miniature, what the research program is ultimately for.

---

## Research Questions for Seminar Discussion

**R43.1** Among the five theoretical gaps (T1–T5), which do you consider most fundamental — the one whose resolution would most significantly advance the entire framework? Justify your answer by showing which other gaps depend on its resolution.

**R43.2** The empirical gaps (E1–E4) require institutional infrastructure that does not yet exist. Which of the following institutions is most positioned to develop this infrastructure: (i) national statistical offices; (ii) international cooperative federations; (iii) academic research consortia; or (iv) the cooperative enterprises themselves? What governance structure would be needed for the chosen institution to maintain the required data quality and independence?

**R43.3** The fifteen research questions (Q1–Q15) are prioritized within three categories. How would you re-prioritize them if: (i) the primary objective is to accelerate the transition (maximize speed); (ii) the primary objective is to minimize transition risks (maximize safety); or (iii) the primary objective is to build the political coalition (maximize feasibility)?

**R43.4** Map the cooperative-regenerative economics framework onto the four major 20th-century research programs in complexity and evolutionary economics: (i) the Santa Fe Institute's complexity economics; (ii) the evolutionary economics of Nelson and Winter; (iii) the institutional economics of Veblen, Commons, and Ayres; and (iv) the ecological economics of Daly, Costanza, and Odum. What does the cooperative-regenerative framework inherit from each? What does it add?

**R43.5** Design the Economics of Cooperation Institute described in Section 43.7.2. Specify: organizational structure (governance, funding, membership); research programs (which of the fifteen questions it would prioritize); field partnerships (which cooperative enterprises, municipalities, and commons institutions it would engage); and output formats (journals, policy briefs, practitioner toolkits, open-source models). What would make this institute genuinely different from existing economics research institutes?

**R43.6** The chapter identifies political philosophy as an adjacent discipline that provides normative grounding for the framework. Specify which normative claims in the book require philosophical grounding that the formal economic analysis alone cannot provide. For each claim, identify the most relevant political philosophy tradition and the key argument that would ground it.

---

*Chapter 44 is the capstone project — not a lecture but a structured design problem that applies the entire analytical toolkit of the trilogy to a real or hypothetical local economy. The project mirrors, at local scale, the design challenge that cooperative-regenerative economics ultimately addresses at global scale: how do we build economies that produce human flourishing without ecological destruction?*
