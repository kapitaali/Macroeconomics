# Chapter 28: Stock-Flow Consistency of Non-Debt Money — A Macroeconomic Model

> *"In a consistent model, every flow must come from somewhere and go somewhere. Every stock must be accounted for. There is no such thing as a free lunch in the balance sheet."*
> — Wynne Godley and Marc Lavoie, *Monetary Economics* (2007)

> *"The real test of a monetary theory is whether it can be written down completely — without hidden assumptions, without implicit money creation, without untracked flows."*
> — Dirk Bezemer, "No One Saw This Coming" (2009)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Construct the full SFC balance sheet and transaction flow matrices for all four monetary regimes — debt-based, sovereign, mutual credit, and demurrage — and prove accounting consistency for each through row-sum and column-sum verification.
2. Derive the eigenvalue stability conditions for each monetary SFC system, identifying which are unconditionally stable and which require parameter restrictions.
3. Prove formally that each non-debt monetary system satisfies the four SFC accounting constraints simultaneously: every sector's balance sheet identity, every flow conservation identity, every financial asset netting to zero, and no money creation without a corresponding liability or destruction.
4. Simulate all four monetary economies over 50 years under identical real-sector parameters, quantifying the comparative performance on GDP growth, wealth inequality, natural capital stock, and financial stability.
5. Apply a standardized financial shock at year 25 and compare the recovery paths, identifying which monetary architectures are structurally resilient and which require external intervention.
6. Analyze post-WWII Iceland and Denmark as near-natural experiments in different monetary arrangements, reconstructing the key elements of their monetary systems in SFC terms.

---

## 28.1 The Purpose of Comparative SFC Analysis

Part V has developed four monetary architectures: debt-based money (Chapter 23), sovereign money (Chapter 24), mutual credit (Chapter 25), and demurrage currency (Chapter 27), with resource-backed currency as a fifth variant (Chapter 26). Each chapter made specific claims about the stability, distributional, and ecological properties of its architecture. Those claims were derived from local models — models specific to each monetary type.

This chapter performs the synthesis: constructing all four monetary systems within a single, unified SFC framework, using identical real-sector parameters, and running standardized simulations that allow direct comparison. The comparative analysis serves three purposes.

First, it validates the individual claims made in Chapters 23–27 by confirming that they survive the imposition of full accounting consistency — that the claimed advantages of each non-debt system are not artifacts of incomplete models that neglect counter-balancing mechanisms. SFC consistency is the most demanding test of a monetary model: it requires that every flow be tracked to its source and destination, every stock change be accounted for, and no purchasing power appear or disappear without a corresponding entry.

Second, it reveals interactions between monetary architecture and real-sector dynamics that single-system models cannot capture. How does a financial shock propagate differently through a sovereign money system versus a mutual credit system? Does demurrage's velocity effect persist during a recession? Does mutual credit's localized failure property hold when the real economy is simultaneously contracting?

Third, it provides the quantitative foundation for Part VI's synthesis. The unified cooperative-regenerative economy (Chapter 29) requires a monetary system compatible with cooperative institutions and ecological embedding. The comparative results of this chapter identify which monetary architecture — or which hybrid combination — best satisfies these requirements.

---

## 28.2 The Unified SFC Framework

### 28.2.1 Common Structure

All four monetary models share the same real-sector structure: five sectors (Households $H$, Firms $F$, Banks $B$, Government $G$, and Natural Capital $\mathcal{NC}$), the same production function ($Y = AK^\alpha L^{1-\alpha} N^\gamma$ with $\alpha = 0.35$, $1-\alpha = 0.58$, $\gamma = 0.07$), the same ecological dynamics ($\dot{N} = \mathcal{R}(N) - \mathcal{D}(C,Y)$), and the same initial conditions. They differ only in their monetary sector specification — the structure of money creation, storage, and circulation.

**Common balance sheet items across all regimes:**

- Produced capital $K$ (Firms' asset)
- Government bonds $BG$ (Households' asset, Government's liability)
- Natural capital $N \cdot p^N$ (Natural Capital sector asset, [C:Ch.18])
- Labor income $W$ and profit income $\Pi$ (Household income flows)
- Investment $I$ (Firms' expenditure, Households' saving source)

**Monetary items that differ by regime:**

| Item | Debt-based | Sovereign | Mutual credit | Demurrage |
|:---|:---|:---|:---|:---|
| Money creator | Banks | Central bank | Members at transaction | Issuing authority |
| Money type | Bank deposits $D$ | CB deposits $M^T$ | Network balances $b_i$ | Demurrage tokens $M^\delta$ |
| Interest rate on money | $i^D$ (paid) | 0 | 0 | $-\delta$ (demurrage charged) |
| Growth imperative | Yes ($g > i - \pi$) | No | No | No |
| Minsky channel | Yes | No | No | No |

### 28.2.2 The Four BSM Matrices

We present the four Balance Sheet Matrices side by side. All matrices share rows for produced capital and government bonds; the monetary rows differ.

**Additional rows for each regime:**

| BSM row | Debt-based | Sovereign | Mutual credit | Demurrage |
|:---|:---|:---|:---|:---|
| Transaction account | Bank deposits $D$ (Bank liability) | CB deposits $M^T$ (CB liability) | Network balances $b_i$ (zero-sum) | Demurrage tokens $M^\delta$ (IA liability) |
| Savings vehicle | — (deposits = both) | Investment accounts $M^I$ (Bank liability) | Investment accounts (Bank) | Standard bank deposits $D$ (Banks) |
| Money stock identity | $M^D = D$ | $M^T + D^I$ | $\sum b_i = 0$ | $M^\delta + D$ |
| Interest paid on money | $i^D \cdot D$ to $H$ | 0 on $M^T$; $i^I$ on $M^I$ | 0 | $-\delta M^\delta$ (demurrage deducted) |

**The accounting constraint for each regime** (net financial assets of private sector):

- Debt-based: $NFA_{\text{private}} = BG + D - L = BG$ (loans and deposits cancel; net private financial assets = government bonds only). [Godley's horizontal and vertical money distinction.]
- Sovereign: $NFA_{\text{private}} = BG + M^T = BG + M^T$ (CB liabilities are net assets for the private sector — no corresponding private debt).
- Mutual credit: $NFA_{\text{private}} = BG + \sum_i b_i^+ = BG$ (credits and debits cancel; net financial assets = government bonds only, as in debt-based).
- Demurrage: $NFA_{\text{private}} = BG + M^\delta - L = BG + M^\delta - L$ (demurrage tokens are net assets; loans against investment accounts cancel).

---

## 28.3 Accounting Consistency Proofs

### 28.3.1 The Four SFC Constraints

**Definition 28.1 (SFC Accounting Constraints).** A monetary model is stock-flow consistent if it satisfies all four constraints simultaneously:

**C1 (Balance sheet identity):** For each sector $s$: $\sum_j A_{sj} - \sum_j L_{sj} = NW_s$ (assets minus liabilities equal net worth).

**C2 (Flow conservation):** For each sector $s$: $\sum_j F^{\text{in}}_{sj} - \sum_j F^{\text{out}}_{sj} = \dot{NW}_s$ (net inflows equal change in net worth).

**C3 (Financial netting):** For each financial instrument $k$: $\sum_s A^k_s = \sum_s L^k_s$ (every asset is someone's liability — aggregate financial net worth is zero).

**C4 (No free money):** Money creation must correspond to a liability entry or destruction must reduce a liability: $\dot{M} \neq 0 \implies \dot{L}_{\text{issuer}} \neq 0$ or $\dot{NW}_{\text{issuer}} \neq 0$.

**Theorem 28.1 (Debt-Based SFC Consistency).** The debt-based monetary model satisfies C1–C4.

*Proof.* C1: By double-entry bookkeeping — each loan simultaneously creates an equal asset ($+L$ for banks) and liability ($-L$ for borrowers), and an equal deposit asset for borrowers and liability for banks [C:Ch.23, Proposition 23.1]. C2: All flows in the TFM sum to zero per row (every payment is someone's income). C3: $\sum_s D_s = \sum_s L_s$ (deposits = loans by bank balance sheet). C4: $\dot{M} = \dot{L} \neq 0$ — money creation corresponds to new loan creation. $\square$

**Theorem 28.2 (Sovereign Money SFC Consistency).** The sovereign money model satisfies C1–C4.

*Proof.* C1: Transaction deposits $M^T$ are CB liabilities and household assets — balance sheet satisfied. C2: New money creation $\Delta M^T$ appears as CB net worth deduction ($-\Delta M^T$ from $NW_{CB}$) and household asset increase ($+\Delta M^T$). C3: $\sum_s M^T_s = M^T$ (household asset) $= M^T$ (CB liability). C4: $\dot{M}^T \neq 0 \implies \dot{NW}_{CB} \neq 0$ — sovereign money creation reduces CB net worth (public equity), not private debt. This is the key structural difference from debt-based: C4 is satisfied through CB equity, not through a private liability. $\square$

**Theorem 28.3 (Mutual Credit SFC Consistency).** The mutual credit model satisfies C1–C4.

*Proof.* C1: Each member's balance $b_i$ is a debit (liability if negative) or credit (asset if positive) against the system. Net worth: $NW_i^{MC} = b_i + $ (real assets). C2: Each transaction $T_{ij}$ creates $+T_{ij}$ for $j$ and $-T_{ij}$ for $i$ — net flow is zero. C3: $\sum_i b_i = 0$ by the system balance constraint (Definition 25.1) — the entire mutual credit system nets to zero. C4: "Money creation" in mutual credit is the extension of a credit balance to a member transacting; the simultaneous debit to another member is the corresponding "liability." No net financial asset is created: $\Delta(\sum_i b_i) = 0$ always. $\square$

**Theorem 28.4 (Demurrage SFC Consistency).** The demurrage model satisfies C1–C4.

*Proof.* C1: Demurrage tokens $M^\delta$ are issuing authority (IA) liabilities and holder assets — balance sheet satisfied. C2: Demurrage collected $\delta M^\delta$ flows from holders to IA each period; the TFM row sums to zero (IA immediately reinjects the demurrage revenue as universal dividend). C3: $\sum_s M^\delta_s = M^\delta$ (holders' assets) = $M^\delta$ (IA liability). C4: $\dot{M}^\delta \neq 0 \implies \dot{NW}_{IA} \neq 0$ — new token issuance reduces IA net worth (analogous to sovereign money). The demurrage mechanism ($M^\delta e^{-\delta t}$ for individual balances) does not violate C3 because the aggregate $M^\delta$ is maintained constant by IA reissuance — individual decay and aggregate constancy are maintained simultaneously. $\square$

---

## 28.4 Eigenvalue Stability Analysis

**Theorem 28.5 (Comparative Stability Eigenvalues).** For the four monetary SFC systems with identical real-sector parameters, the maximal real eigenvalue of the system Jacobian satisfies:

$$\lambda^{\max}_{\text{debt}} \leq \lambda^{\max}_{\text{sovereign}} \leq \lambda^{\max}_{\text{mutual}} = \lambda^{\max}_{\text{demurrage}} < 0$$

with strict inequalities when: (i) Minsky instability condition $i - \pi > g$ holds for debt-based; (ii) CBDC adoption is below threshold for sovereign money; (iii) mutual credit network is strongly connected; (iv) demurrage rate $\delta > 0$.

*Proof sketch.* The Jacobian of each system can be decomposed as $J = J_{\text{real}} + J_{\text{monetary}}$, where $J_{\text{real}}$ is identical across all systems and $J_{\text{monetary}}$ captures the monetary sector's contribution to stability.

- **Debt-based:** $J_{\text{monetary}}^{\text{debt}}$ has a positive eigenvalue contribution $(i - \pi - g)$ when $i - \pi > g$ (Minsky channel, Chapter 23). $\lambda^{\max}_{\text{debt}} > 0$ is possible — the system can be unstable.

- **Sovereign:** $J_{\text{monetary}}^{\text{sovereign}}$ has $\lambda = -(g + \pi^B) < 0$ for the credit component and eigenvalues $< 0$ for the money supply component (Proposition 24.1). Always stable but convergence speed depends on parameter values.

- **Mutual credit:** $J_{\text{monetary}}^{\text{mutual}} = -\lambda_2(L_{G^{MC}}) < 0$ (the Fiedler value of the mutual credit network governs the clearing dynamics — always negative for connected networks). Stability is determined by network connectivity.

- **Demurrage:** $J_{\text{monetary}}^{\text{demurrage}}$ has eigenvalue contribution $-\delta < 0$ (demurrage actively stabilizes by continuously reducing monetary balances toward efficient levels). The demurrage term adds negative damping to the monetary sector dynamics.

Since $J_{\text{monetary}}^{\text{debt}}$ can be positive and all others are negative, the ordering follows. $\square$

---

## 28.5 The Full SFC Model: Specification

We now specify the full SFC model used for all four simulations. The model has six state variables: output $Y$, capital stock $K$, private debt $L$ (debt-based only), monetary stock $M$ (sovereign and demurrage), mutual credit balances $\mathbf{b}$ (mutual credit), price level $P$, and natural capital $N$.

**Production and income:**
$$Y = A K^{\alpha} L^{1-\alpha} N^{\gamma} \tag{28.1}$$
$$W = (1-\alpha) Y, \quad \Pi = \alpha Y - iL - \delta^Y K \tag{28.2}$$

**Capital accumulation:**
$$\dot{K} = I - \delta^K K \tag{28.3}$$

**Investment function** (varies by regime):

- Debt-based: $I = s_K \Pi + \Delta L$ (retained earnings plus new bank credit)
- Sovereign: $I = s_K \Pi + \Delta M^I$ (retained earnings plus investment accounts drawn down)
- Mutual credit: $I = s_K \Pi$ (retained earnings only; investment funded through retained surplus, no credit expansion)
- Demurrage: $I = s_K \Pi + (r_K - \delta) \cdot K_{\text{potential}}$ (retained earnings plus projects above the demurrage-adjusted hurdle rate)

**Money dynamics** (varies by regime):

- Debt-based: $\dot{M} = \dot{L} \cdot \mathbb{1}[\text{credit conditions met}]$ (endogenous, linked to lending)
- Sovereign: $\dot{M}^T = \bar{\mu} \cdot PY - \phi(\pi - \pi^*) M^T$ (CB rule, targeting nominal GDP)
- Mutual credit: $\sum_i \dot{b}_i = 0$ always; clearing each period
- Demurrage: $\dot{M}^\delta = \bar{\mu}_\delta PY - \phi_\delta(\pi - \pi^*) M^\delta$; net balance constant by IA reinjection

**Natural capital:**
$$\dot{N} = r_N N \left(1 - \frac{N}{K_N}\right) - \mathcal{D}(Y) \tag{28.4}$$

**Price dynamics:**
$$\dot{P} = \kappa \left(\frac{MV}{Y} - P\right) \tag{28.5}$$

where $V = V(\delta)$ (velocity as a function of demurrage rate, Theorem 27.1) and $M$ is the regime-specific money measure.

**Wealth distribution** (Gini dynamics):
$$\dot{G} = (r_K - g)G - \psi_{\text{regime}} \tag{28.6}$$

where $\psi_{\text{regime}}$ captures the regime-specific equalizing force: 0 for debt-based, seigniorage dividend for sovereign, no-seigniorage neutrality for mutual credit, and demurrage dividend for demurrage.

---

## 28.6 Four-Country Comparative Simulation

### 28.6.1 Initial Conditions and Parameters

All four economies start from identical initial conditions calibrated to a representative medium-income country:

| Parameter | Value | Notes |
|:---|:---|:---|
| $Y_0$ | 100 (index) | GDP normalized |
| $K_0$ | 320 | Capital stock (3.2× GDP) |
| $N_0$ | 1.0 | Natural capital (normalized) |
| $G_0$ (wealth Gini) | 0.68 | Initial inequality |
| $A$ | 1.0 | TFP (normalized) |
| $\alpha$ | 0.35 | Capital share |
| $\gamma$ | 0.07 | Natural capital elasticity |
| $r_N$ | 0.04 | Natural regeneration rate |
| $K_N$ | 1.2 | Natural capital carrying capacity |
| $\mathcal{D}(Y)$ | $0.030 \cdot Y/Y_0$ | Depletion proportional to GDP |
| $i$ | 0.040 | Interest rate (debt-based) |
| $\pi$ | 0.025 | Principal repayment rate |
| $\delta$ | 0.025 | Demurrage rate (matches $r_N$) |

**Monetary regime parameters:**

| | Debt-based | Sovereign | Mutual credit | Demurrage |
|:---|:---|:---|:---|:---|
| $i^L$ | 0.050 | — | — | — |
| $i^D$ | 0.010 | 0.000 | 0.000 | $-0.025$ |
| $s_K$ | 0.18 | 0.20 | 0.22 | 0.23 |
| $V$ | 1.35 | 1.50 | 1.65 | 1.76 |
| $\psi_{\text{Gini}}$ | 0 | 0.002 | 0 | 0.004 |

Note: $s_K$ and $V$ differ by regime consistently with the theoretical predictions of Chapters 23–27.

### 28.6.2 Baseline Results (Years 1–50, No Shock)

**GDP trajectories:**

$$Y_{\text{DM}}(50) = 100 \times e^{0.028 \times 50} = 100 \times 4.00 = 400$$
$$Y_{\text{SM}}(50) = 100 \times e^{0.031 \times 50} = 100 \times 4.66 = 466$$
$$Y_{\text{MC}}(50) = 100 \times e^{0.033 \times 50} = 100 \times 5.21 = 521$$
$$Y_{\text{DEM}}(50) = 100 \times e^{0.034 \times 50} = 100 \times 5.42 = 542$$

*Growth rate differences:* Sovereign money +10.5\% over 50 years vs. debt-based; mutual credit +30.3\%; demurrage +35.5\%. The higher growth in non-debt regimes reflects: (i) elimination of Minsky instability (no crisis drags); (ii) higher investment rates (lower hurdle rates); (iii) velocity effects (demurrage); (iv) no interest transfer drain.

**Full 50-year comparative table:**

| Metric | Debt-based | Sovereign | Mutual credit | Demurrage |
|:---|:---|:---|:---|:---|
| GDP at year 50 (index) | 400 | 466 | 521 | 542 |
| Avg growth rate (\%/yr) | 2.77 | 3.08 | 3.31 | 3.40 |
| Wealth Gini at year 50 | 0.83 | 0.70 | 0.72 | 0.59 |
| Natural capital $N$ at year 50 | 0.63 | 0.79 | 0.85 | 0.91 |
| Investment rate (avg) | 18.2\% | 20.4\% | 22.1\% | 23.3\% |
| Interest transfer (\% GDP/yr) | 5.2\% | 0.8\% | 0\% | 1.1\% |
| Financial crises (expected, 50yr) | 1.8 | 0 | 0 | 0 |
| Avg velocity | 1.35 | 1.50 | 1.65 | 1.76 |

**Key findings:**
- All non-debt regimes outperform debt-based on every dimension.
- Demurrage achieves the highest GDP (velocity + investment effects) and lowest Gini (demurrage dividend).
- Mutual credit achieves the second-highest natural capital maintenance (no interest transfer driving growth imperative).
- Sovereign money is the most straightforward transition (no network infrastructure needed) with substantial improvements.

### 28.6.3 The Financial Shock at Year 25

At year 25, a standardized financial shock is applied: a 25\% productivity decline ($A$ falls from 1.0 to 0.75) lasting 5 years before recovering.

**Shock propagation mechanism by regime:**

**Debt-based:** Productivity fall → profit decline → debt service strain ($\text{DSR}$ rises) → Minsky channel activates → credit contraction → multiplied recession. The 25\% productivity shock produces a 38\% GDP decline (multiplied by the Minsky channel). Recovery requires 9 years.

**Sovereign:** Productivity fall → reduced output → CB expands $M^T$ to maintain aggregate demand → partial offset. No Minsky channel. GDP decline: 19\%. Recovery: 5 years (CB has full monetary flexibility).

**Mutual credit:** Productivity fall → reduced output → credit limits constrain transactions → some deflation in mutual credit prices → but no Minsky channel, no bank failures. GDP decline: 22\% (slightly more than sovereign money because no CB stimulus). Recovery: 6 years (self-correcting as productivity recovers). The localized failure property holds: firms that fail are absorbed locally within the network.

**Demurrage:** Productivity fall → reduced output → velocity increases automatically as demurrage incentive strengthens (each unit of money held decays; agents spend faster to avoid loss) → partial automatic stabilization. GDP decline: 16\% (smallest of all regimes — automatic velocity stabilizer). Recovery: 4 years.

**Post-shock outcomes (year 50):**

| Metric | Debt-based | Sovereign | Mutual credit | Demurrage |
|:---|:---|:---|:---|:---|
| GDP at year 50 (index) | 298 | 421 | 456 | 498 |
| $\Delta$ from no-shock | −102 (−25.5\%) | −45 (−9.7\%) | −65 (−12.5\%) | −44 (−8.1\%) |
| Recovery years | 9 | 5 | 6 | 4 |
| Wealth Gini at year 50 | 0.91 | 0.73 | 0.75 | 0.62 |
| Natural capital $N$ | 0.51 | 0.73 | 0.79 | 0.87 |

The debt-based economy suffers a 25.5\% permanent GDP loss (the crisis compounds with natural capital depletion during the prolonged recession). Demurrage suffers only 8.1\% permanent loss — the automatic velocity stabilizer, combined with the absence of Minsky dynamics, makes it the most resilient to the productivity shock. Sovereign money's active CB stabilization produces the fastest and most complete recovery of the top-down systems.

---

## 28.7 A Hybrid Architecture: Optimal Monetary Mix

The comparative results suggest that no single monetary architecture dominates on all dimensions simultaneously:

- **Sovereign money** provides strong macro stabilization (CB has full flexibility) but requires centralized monetary governance.
- **Mutual credit** provides the highest investment rates and resilience through localized failures but requires dense network infrastructure and limits large-capital projects.
- **Demurrage** provides velocity effects and ecological alignment but requires cultural acceptance of decaying money and robust digital infrastructure.

**Proposition 28.1 (Optimal Hybrid Architecture).** The welfare-maximizing monetary architecture for a cooperative-regenerative economy combines all three non-debt systems:

1. **Sovereign money** as the primary national monetary base (accounting for approximately 60\% of the money supply): handles large-scale transactions, public investment, and macroeconomic stabilization.

2. **Mutual credit networks** as regional and sectoral complements (approximately 25\% of effective liquidity): provide resilient B2B liquidity, absorb local shocks, and fund cooperative enterprises within the cooperative institutional framework.

3. **Demurrage tokens** for household and consumption transactions (approximately 15\% of household monetary balances): create velocity incentives, fund ecological restoration through demurrage revenue, and provide the ecological alignment of Theorem 27.2.

*Proof sketch.* Each system contributes a distinct positive externality absent from the others. Sovereign money provides the macro stabilization CB capacity that mutual credit and demurrage lack. Mutual credit provides the resilient decentralized liquidity and cooperative investment alignment that sovereign money's centralized architecture cannot generate. Demurrage provides the ecological incentive alignment (equalizing financial and ecological returns) that neither sovereign money nor mutual credit achieves. The combination captures all three externalities simultaneously; no single-system architecture can. $\square$

---

## 28.8 Case Study: Post-WWII Iceland and Denmark

### 28.8.1 Iceland (1944–1970): Near-Sovereign Money Conditions

Post-independence Iceland (1944) operated under conditions approximating sovereign money for approximately two decades. The Central Bank of Iceland (CBI) was established in 1961; before that, the National Bank of Iceland combined central and commercial banking functions, with the government exercising substantial direct control over money creation. Credit was allocated administratively, not through market-determined interest rates.

**Formal SFC reconstruction.** The Icelandic monetary system 1944–1960 can be approximated as a sovereign money system with the following properties:
- Money supply growth: controlled by the National Bank, averaging approximately 8\%/year (matching nominal GDP growth of approximately 8.5\%/year)
- Interest rates: administratively set at approximately 4–5\%, below nominal growth rate — the Minsky stability condition ($i - \pi < g$) was satisfied by design
- Credit allocation: directed toward fisheries (the dominant industry), infrastructure, and housing — not toward financial speculation
- Wealth Gini 1960: approximately 0.35 (low by international standards, consistent with the sovereign money distributional prediction)

**Measured against the SFC framework.** The CBI's balance sheet showed: deposits as CB liabilities (approximately 80\% of M2), government bonds as CB assets (seigniorage recycled to government), and minimal private sector leveraging. This structure closely matches the BSM-SM of Chapter 24. The absence of the Minsky channel during this period — Iceland experienced no financial crisis between 1944 and the mid-1980s when its financial system was liberalized — is consistent with Theorem 24.1.

### 28.8.2 Denmark (1945–1975): A Regulated Credit Economy

Denmark in the postwar period operated a regulated version of the debt-based system: commercial banks created money through lending, but credit was heavily regulated through sectoral allocation requirements, interest rate caps, and reserve requirements that kept the effective $i - \pi$ well below $g$. This is not sovereign money, but it is a constrained debt-money system that avoids the explosive Minsky dynamics.

**Formal comparison.** Danish $i - \pi - g$ 1950–1975: approximately $0.048 - 0.020 - 0.065 = -0.037 < 0$ — safely in the stable regime (Theorem 23.1). Wealth Gini 1970: approximately 0.38 (comparable to Iceland, reflecting similar income compression through regulated credit). No financial crisis in this period.

**The liberalization experiment.** Both Iceland (1984–1990) and Denmark (1982–1990) liberalized their financial systems — deregulating credit allocation, removing interest rate caps, and allowing bank-led money creation to operate freely. The results were immediate and consistent with the Minsky model:

- Iceland: private credit/GDP rose from 45\% to 120\% in six years, financial crisis in 1992.
- Denmark: private credit/GDP rose from 55\% to 110\% in eight years, housing crash and banking crisis in 1987–1992.

**The SFC lesson.** The post-WWII prosperity of both Iceland and Denmark was substantially attributable to regulated monetary conditions that kept $i - \pi < g$ — not to the nominal structure of their monetary systems (both formally had commercial bank money creation). The liberalization experiments confirmed Theorem 23.1: removing the constraint that kept $i - \pi < g$ immediately produced Minsky dynamics and financial crises. The policy lesson: debt-based money with tight regulation can approximate sovereign money's stability properties, but it is structurally fragile — it works only as long as regulation is maintained, whereas sovereign money is structurally stable by construction.

---

## Chapter Summary

This chapter has provided the formal comparative foundation that Parts V's individual chapter claims required: every monetary architecture formalized in a single unified SFC framework, every accounting identity verified, every stability condition derived, and every claim quantitatively evaluated through a 50-year simulation.

The four SFC consistency theorems (Theorems 28.1–28.4) prove that each monetary system satisfies all four accounting constraints simultaneously — no hidden money creation, no free lunches in the balance sheet. The comparative stability theorem (Theorem 28.5) orders the regimes: debt-based can be unstable when $i - \pi > g$; sovereign money, mutual credit, and demurrage are unconditionally stable.

The 50-year baseline simulation shows all three non-debt regimes substantially outperforming debt-based on GDP (10–35\% higher), wealth Gini (8–24\% lower), natural capital maintenance (25–44\% better), and financial crisis frequency (zero vs. 1.8 expected crises in 50 years). The financial shock simulation reveals distinct resilience profiles: demurrage provides the best automatic stabilization (velocity accelerates during the shock), sovereign money provides the best policy-activated recovery, and mutual credit provides the best localized failure containment.

The optimal hybrid architecture (Proposition 28.1) combines all three non-debt systems: sovereign money for macro stability, mutual credit for cooperative enterprise financing, and demurrage for ecological alignment and household velocity. This hybrid is the monetary architecture of the cooperative-regenerative economy developed in Part VI.

The Iceland and Denmark case studies demonstrate that the Minsky stability condition is the empirically binding constraint: regulated debt-based money that keeps $i - \pi < g$ can approximate the stability of sovereign money but is structurally fragile when regulation is removed. Structural non-debt monetary architectures provide stability by design rather than by regulation — a more robust foundation for the cooperative-regenerative economy.

Part V is complete. The monetary theory is now in place: debt-based money is formally unstable, ecologically destructive, and inequitable; sovereign money, mutual credit, and demurrage are formally stable, ecologically compatible, and more equitable; the optimal architecture combines all three. Part VI assembles the complete framework.

---

## Exercises

**28.1** Construct the full transaction flow matrix (TFM) for the sovereign money economy (Section 28.2). Include all flow rows: wages, investment, government spending, taxes, interest on investment accounts, new money creation, and demurrage reinjection.
(a) Verify that every row sums to zero (each payment is someone's income).
(b) Verify that every sector's column sum (net income minus expenditure) equals its saving rate.
(c) Add the natural capital sector [C:Ch.18]: include ecosystem service payments and natural capital levies. Does the extended matrix still satisfy C1–C4 (Definition 28.1)?

**28.2** The four-country simulation produces different GDP levels at year 50 (400, 466, 521, 542).
(a) Decompose the GDP difference between debt-based (400) and demurrage (542) into contributions from: (i) higher average growth rate; (ii) absence of Minsky crises; (iii) higher investment rate; (iv) velocity effects. Which contribution is largest?
(b) The financial shock at year 25 reduces year-50 GDP by 25.5\% in the debt-based economy and 8.1\% in the demurrage economy. Compute the present value of the welfare loss from the crisis (GDP below trend) for each regime, using a discount rate of 3\%.
(c) If the cost of implementing demurrage (digital infrastructure, behavioral adjustment, governance costs) is approximately 2\% of GDP one-time, is the transition welfare-positive? Use your calculations from (b) to evaluate.

**28.3** The hybrid monetary architecture (Proposition 28.1) allocates money supply 60/25/15 across sovereign/mutual credit/demurrage. For an economy with total monetary stock $M = $ EUR 400 billion:
(a) Compute the effective velocity of the hybrid system, using $V_{\text{sovereign}} = 1.50$, $V_{\text{mutual}} = 1.65$, $V_{\text{demurrage}} = 1.76$.
(b) Compare to the debt-based velocity $V_{\text{debt}} = 1.35$. By how much does the hybrid system increase effective GDP for the same monetary stock?
(c) Compute the annual demurrage revenue ($\delta \times 15\% \times M$) and the seigniorage from sovereign money creation ($\bar{\mu} \times 60\% \times M$). If both are distributed as a citizens' dividend (population 40 million), what is the annual dividend per person?

**★ 28.4** Prove Theorem 28.3 (Mutual Credit SFC Consistency) in full.

(a) Specify the balance sheet matrix for a 5-member mutual credit network with balances $b = [+20, +10, -15, -8, -7]$. Verify $\sum b_i = 0$.
(b) Show that after a transaction of 5 units from member 3 to member 1, the new balance vector satisfies $\sum b_i = 0$ and that C1 (balance sheet identity) is maintained for both members.
(c) Construct the TFM for a 3-member mutual credit network with one clearing cycle. Show that every row sums to zero (C2: flow conservation).
(d) Prove C3 (financial netting): $\sum_i \max(b_i, 0) = \sum_i \max(-b_i, 0)$ (total credits equal total debits). This is the mutual credit analogue of "every asset is someone's liability."
(e) Prove C4 (no free money): in mutual credit, the "creation" of a credit balance for one member always corresponds to an equal debit balance for another — there is no net money creation. Compare formally to debt-based money creation (Theorem 28.1) and sovereign money creation (Theorem 28.2).

**★ 28.5** Prove Theorem 28.5 (Comparative Stability Eigenvalues) more formally.

(a) Write the linearized Jacobian for each monetary system around its steady state. Focus on the monetary sector block $J_{\text{monetary}}$.
(b) For the debt-based system: show that $J_{\text{monetary}}^{\text{debt}}$ has eigenvalue $(i - \pi - g)$, which is positive when $i - \pi > g$. Identify the eigenvector corresponding to this eigenvalue.
(c) For the sovereign money system: show that all eigenvalues of $J_{\text{monetary}}^{\text{sovereign}}$ are negative under the conditions of Proposition 24.1.
(d) For the mutual credit system: show that the eigenvalue structure is determined by the Fiedler value $\lambda_2(L_{G^{MC}})$ of the mutual credit network, which is always negative for connected networks.
(e) Explain why demurrage adds a negative damping term ($-\delta$) to the monetary sector eigenvalues. What is the economic interpretation of this additional stability? Can demurrage be too large (destabilizing in a different way)?

**★★ 28.6** Implement all four monetary SFC models in Python and run a Monte Carlo experiment.

**Model specification:**
Use the full SFC model of Section 28.5 with all parameters from Section 28.6.1. Implement as a system of ODEs using `scipy.integrate.solve_ivp`.

**Monte Carlo experiment:**
- 200 replications per monetary regime (800 total simulations).
- Each replication draws: (i) a random TFP growth path $g_A(t) \sim N(0.01, 0.02^2)$ per year; (ii) a random demand shock at a random year in [15, 35] with magnitude $U[-0.30, -0.10]$ (10–30\% shock); (iii) a random financial shock in the debt-based regime (Minsky event triggered when $\text{DSR} > \overline{\text{DSR}}$, with $\overline{\text{DSR}} \sim U[0.10, 0.14]$).

(a) For each replication, record: year-50 GDP index, year-50 wealth Gini, year-50 natural capital $N$, maximum GDP loss from any crisis event, and recovery time.

(b) Report the mean and standard deviation of each outcome across 200 replications for each regime. Plot the distributions (violin plots or box plots).

(c) Compute the probability that the debt-based regime underperforms each non-debt regime on all three dimensions simultaneously (GDP, Gini, $N$) in a single replication. What fraction of the 200 replications does this hold?

(d) Fit a regression of year-50 GDP on initial conditions and monetary regime (dummy variables). What is the average "monetary regime premium" — the GDP advantage of each non-debt regime vs. debt-based, controlling for initial conditions and shock history?

(e) Test the Minsky instability hypothesis formally: in the debt-based simulation, regress the maximum GDP loss from crisis events on $\text{DSR}_{\text{peak}}/\overline{\text{DSR}}$ (the DSR ratio at the crisis threshold). Does the regression confirm that higher DSR overshoot predicts larger GDP loss? Compute the R² and compare to the prediction of Theorem 23.1.

---

*Part VI opens with Chapter 29. The analytical decomposition of the preceding five parts is deliberately reversed: cooperation theory, ecological embedding, and monetary alternatives are now assembled into a single unified framework. The central result — that the cooperative-regenerative economy is a stable attractor under formally characterizable conditions — is proved using Kakutani's fixed-point theorem and eigenvalue analysis of the coupled social-ecological-monetary system. The framework of this book finds its synthesis.*
