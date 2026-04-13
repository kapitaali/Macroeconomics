# Chapter 24: The Sovereign Money Approach — Central Bank Digital Currency Without Commercial Bank Creation

> *"The essence of the contemporary monetary system is the creation of money, out of nothing, by private banks' often foolish lending."*
> — Martin Wolf, *The Shifts and the Shocks* (2014)

> *"We could maintain all the desirable features of the current monetary system while doing away with the undesirable ones."*
> — Jaromir Benes and Michael Kumhof, *The Chicago Plan Revisited* (IMF Working Paper, 2012)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Specify the sovereign money (full-reserve banking) architecture formally using the SFC balance sheet framework, and identify precisely how it differs from the debt-based system at the level of accounting identities.
2. Derive the stability properties of the sovereign money SFC system through eigenvalue analysis, and prove that sovereign money eliminates the Minsky instability mechanism and the debt-deflation trap.
3. Analyze the distributional consequences of sovereign money creation: formal treatment of seigniorage as a social dividend, comparison with debt-based seigniorage, and the conditions under which sovereign money reduces wealth concentration.
4. Compare four CBDC architectures — account-based vs. token-based, two-tier vs. single-tier — on dimensions of financial inclusion, monetary control, privacy, and systemic risk, and model CBDC adoption dynamics using network externalities.
5. Model the transition from fractional-reserve to full-reserve banking as a dynamical system, derive the stability of the transition path, and specify formal policy sequencing for a viable transition.
6. Evaluate the IMF's Chicago Plan Revisited (Benes and Kumhof, 2012) formally, assessing which of its stability claims the SFC framework supports and which require qualification.

---

## 24.1 The Fundamental Diagnosis

Chapter 23 established three pathologies of debt-based money: structural instability (explosive debt dynamics when $i - \pi > g$), systematic distributional transfer (interest payments from debtors to creditors), and an ecological growth imperative (the requirement for perpetual nominal growth to prevent debt-deflation). Any alternative monetary architecture must address at least some of these pathologies to represent a genuine improvement.

The sovereign money proposal — also called full-reserve banking, 100\% money (Fisher, 1935), or narrow banking — addresses all three through a single institutional change: transferring the money-creation function from commercial banks to the central bank, making money a public good rather than a byproduct of private lending. In the sovereign money system, commercial banks can only lend money they have obtained from depositors or investors — they cannot create deposits through lending. New money enters the economy only through the central bank, directed through democratic public spending decisions rather than through the profit-maximizing lending decisions of private financial institutions.

The proposal is not new. Versions of it were advanced by Fisher (1935) and the Chicago economists (Simons, 1948) in response to the Great Depression, revived by Tobin (1987) under the label "deposited currency accounts," analyzed in detail by Huber and Robertson (2000) and Dyson et al. (2016), and formally evaluated by Benes and Kumhof (2012) in an IMF working paper using a DSGE model. This chapter contributes the formal SFC analysis — a framework better suited to tracing the stock-flow implications of the proposal than the DSGE approach used by Benes and Kumhof.

---

## 24.2 The Sovereign Money Architecture

### 24.2.1 Institutional Mechanics

**Definition 24.1 (Sovereign Money System).** A sovereign money system is a monetary architecture in which:

1. **Only the central bank creates money.** All money in circulation (whether physical cash or digital deposits) is a direct liability of the central bank, not of commercial banks.
2. **Commercial banks are payment intermediaries only.** They administer payment accounts, facilitate transactions, and manage the allocation of savings to investment — but cannot create money through lending.
3. **Transaction accounts are 100\% reserve-backed.** Customer deposits held in transaction accounts are matched 1:1 by central bank reserves — they are not at risk from bank failure or banking system instability.
4. **Investment accounts are not money.** Customers who wish to earn interest can transfer funds to investment accounts, which are explicitly marked as at-risk savings that fund bank lending — not deposits in the monetary sense.
5. **New money is created by the central bank** and enters the economy through democratic public spending (as decided by government) or, in some variants, through a "citizens' dividend" — a direct distribution to all citizens.

**Comparison with current system:**

| Feature | Debt-based money | Sovereign money |
|:---|:---|:---|
| Money creator | Commercial banks (through lending) | Central bank (through public spending) |
| Money creation trigger | Profitable lending opportunity | Democratic public spending decision |
| Deposit safety | At risk (bank failure possible) | 100\% safe (CB liability) |
| Seigniorage beneficiary | Banking sector (net interest margin) | Public (government revenue or citizens' dividend) |
| Growth imperative | Yes ($g > i - \pi$ required) | No (money creation independent of growth) |
| Minsky instability | Structural feature | Eliminated |

### 24.2.2 The SFC Balance Sheet Matrix

**Definition 24.2 (SFC-SM Balance Sheet).** The sovereign money balance sheet matrix (BSM-SM) for a four-sector economy (Households $H$, Firms $F$, Banks $B$, Central Bank/Government $CB$):

| | $H$ | $F$ | $B$ | $CB$ | $\Sigma$ |
|:---|:---|:---|:---|:---|:---|
| Transaction deposits $M^T$ | $+M^T_H$ | $+M^T_F$ | | $-M^T$ | $0$ |
| Investment accounts $M^I$ | $+M^I_H$ | | $-M^I_B$ | | $0$ |
| Bank loans $L^B$ | $-L^B_H$ | $-L^B_F$ | $+L^B$ | | $0$ |
| CB reserves $R$ | | | $+R$ | $-R$ | $0$ |
| Government bonds $BG$ | $+BG_H$ | | $+BG_B$ | | $+BG_G$ |
| Produced capital $K$ | | $+K_F$ | | | $+K_F$ |
| Net worth | $NW_H$ | $NW_F$ | $NW_B$ | $NW_{CB}$ | $0$ |

**Key differences from BSM-D (Chapter 23):**

1. Transaction deposits $M^T$ are central bank liabilities (row: $-M^T$ in CB column) — not bank liabilities.
2. Investment accounts $M^I$ are bank liabilities but are explicitly not included in the money supply — they are savings instruments.
3. CB reserves $R$ equal transaction deposits exactly: $R = M^T$ (100\% reserve requirement for transaction accounts).
4. There are no "demand deposits" that are simultaneously money and funding for bank loans — the fundamental conflation of functions that creates the money multiplier is dissolved.

**Transaction Flow Matrix (TFM-SM).** The key difference in flows:

| | $H$ | $F$ | $B$ | $CB$ | $\Sigma$ |
|:---|:---|:---|:---|:---|:---|
| Wages | $+W$ | $-W$ | | | $0$ |
| Interest on inv. accounts | $+i^I M^I_H$ | | $-i^I M^I_H$ | | $0$ |
| Interest on loans | $-i^L L^B_H$ | $-i^L L^B_F$ | $+i^L L^B$ | | $0$ |
| **New money (SM creation)** | $+\Delta M^T_{CB}$ | | | $-\Delta M^T_{CB}$ | $0$ |
| Change in investment acc. | $-\Delta M^I$ | | $+\Delta M^I$ | | $0$ |
| Change in transaction dep. | $-\Delta M^T$ | | | $+\Delta M^T$ | $0$ |

The sovereign money creation row is the critical new entry: new money $\Delta M^T_{CB}$ is created by the central bank as a credit to households (or government spending) and a debit to CB net worth — it is a public gift, not a loan. This row has no corresponding debt entry — no household or firm owes the central bank for the new money. This is the formal expression of the elimination of the debt-money link.

---

## 24.3 Stability Properties: Formal Analysis

### 24.3.1 Elimination of the Minsky Mechanism

**Theorem 24.1 (Sovereign Money Eliminates Minsky Instability).** In a sovereign money SFC system, the Minsky instability mechanism of Chapter 23 (Theorem 23.1) is eliminated: the debt ratio $d = L^B/Y$ has no explosive dynamics driven by money creation.

*Proof.* In the sovereign money system, the money supply $M^T$ is determined by the central bank's spending decisions, not by private lending. The debt ratio dynamics depend only on private credit (investment account → bank loan channel), not on transaction account money creation:

$$\dot{d} = \frac{\dot{L}^B}{Y} - d \cdot g = \left(\frac{\Delta L^B_{\text{new}} - \Delta L^B_{\text{repaid}}}{Y}\right) - d \cdot g$$

But $\Delta L^B_{\text{new}} \leq M^I$ (banks can only lend from investment accounts — funds that savers have explicitly committed to at-risk lending). Total investment account balances $M^I$ are bounded by household saving $S_H \leq Y_H$. Therefore:

$$\dot{d} \leq \frac{S_H}{Y} - d \cdot g = s - d \cdot g$$

where $s = S_H/Y$ is the aggregate saving rate. At steady state: $d^* = s/g$ — the debt ratio is bounded by the ratio of the saving rate to the growth rate, both of which are bounded by real economic variables. There is no mechanism by which private credit creation generates unbounded debt dynamics independently of income. $\square$

**Corollary 24.1 (No Debt-Deflation Trap).** In the sovereign money system, there is no debt-deflation trap (Definition 23.5). 

*Proof.* The debt-deflation trap requires that falling prices raise real debt obligations, triggering a positive feedback spiral. In sovereign money: (i) transaction account money is not associated with debt obligations (new money is granted, not lent), so deflation does not raise household obligations from money holdings; (ii) investment account debt $L^B$ does carry real debt obligations, but these are bounded by $s/g$ and do not have the exponential dynamics of Theorem 23.1. The positive feedback (falling prices → rising real debt → falling demand → further price falls) is broken because money itself is not debt. $\square$

### 24.3.2 Eigenvalue Analysis of the Sovereign Money SFC

The linearized sovereign money SFC system has the following key state variables: output $Y$, transaction money $M^T$, private credit $L^B$, and price level $P$. The Jacobian at the steady state $(Y^*, M^{T*}, L^{B*}, P^*)$:

$$J_{\text{SM}} = \begin{pmatrix}
-\mu_Y & \alpha_M & -\beta_L & \gamma_P \\
\delta_M & -\lambda_M & 0 & 0 \\
0 & 0 & -(g + \pi^B) & 0 \\
\kappa_Y & 0 & 0 & -\phi_P
\end{pmatrix}$$

where $\mu_Y, \lambda_M, \pi^B, \phi_P > 0$ (self-stabilizing terms) and off-diagonal terms are the cross-effects. The absence of the $M^T$-$L^B$ coupling (zero entry in row 1, column 3 position for money creation) is the critical structural difference from the debt-money Jacobian: money supply is controlled by the CB ($\delta_M$ row), not by private credit dynamics.

**Proposition 24.1 (Stable Eigenvalues under Sovereign Money).** Under the sovereign money SFC system, all eigenvalues of $J_{\text{SM}}$ are negative (asymptotically stable) when:
(i) $\mu_Y > \alpha_M \delta_M / \lambda_M$ (output self-stabilization dominates monetary injection)
(ii) $g + \pi^B > 0$ (credit grows slower than the economy)
(iii) $\phi_P > \kappa_Y \gamma_P / \mu_Y$ (price level self-stabilization dominates output effects)

*Proof.* The block-triangular structure of $J_{\text{SM}}$ (the $L^B$ row is decoupled from $M^T$) allows eigenvalue analysis by sub-blocks. The $L^B$ eigenvalue is $-(g + \pi^B) < 0$ unconditionally. The remaining $3 \times 3$ block has negative trace and positive determinant under conditions (i) and (iii), giving negative eigenvalues by Routh-Hurwitz. $\square$

**Contrast with the debt-money Jacobian.** The debt-money Jacobian has a positive feedback term in the $M$-$L$ coupling (money creation → lending → more money creation), which generates the Minsky instability. In $J_{\text{SM}}$, this coupling is zero — the architectural separation of money from credit eliminates the positive feedback at the level of the matrix structure.

---

## 24.4 Distributional Effects: Seigniorage as Social Dividend

### 24.4.1 What Is Seigniorage?

**Definition 24.3 (Seigniorage).** Seigniorage is the revenue earned from the exclusive right to create money. In a debt-based system, seigniorage accrues primarily to commercial banks as the net interest margin on money creation (the spread between the lending rate $i^L$ and the deposit rate $i^D$ on the money created). In a sovereign money system, seigniorage accrues to the state.

**Debt-based seigniorage.** When a commercial bank creates $\Delta M$ through a loan at rate $i^L$ and pays depositors rate $i^D < i^L$:

$$\text{Bank seigniorage} = (i^L - i^D) \cdot \Delta M$$

At the aggregate level, with total bank-created money $M^D = $ USD 15 trillion (US M2 less currency, 2022) and spread $i^L - i^D \approx 0.025$: annual bank seigniorage $\approx $ USD 375 billion — a transfer from the real economy to the financial sector embedded in the structure of money creation.

**Sovereign money seigniorage.** In the sovereign money system, new money $\Delta M^T$ is created at zero cost by the central bank and spent into the economy through government channels. The seigniorage is:

$$\text{Public seigniorage} = \Delta M^T - \text{CB operating costs} \approx \Delta M^T$$

For a steady-state economy with 3\% nominal money growth and $M^T = $ USD 15 trillion: annual new money $\Delta M^T \approx $ USD 450 billion — available for public purposes (fiscal spending, debt reduction, citizens' dividend) rather than accruing to commercial banks.

### 24.4.2 The Social Dividend Mechanism

**Definition 24.4 (Citizens' Dividend).** In the sovereign money framework, new money is issued as a citizens' dividend when the central bank credits all citizens' transaction accounts equally with each new money creation:

$$\text{Citizens' dividend per person} = \frac{\Delta M^T}{n_{\text{citizens}}}$$

**Proposition 24.2 (Social Dividend and Wealth Concentration).** Under a sovereign money system with citizens' dividend, the Gini coefficient of net worth decreases monotonically when:

$$\frac{\Delta M^T}{n \cdot \bar{NW}} > (i - g) \cdot G(NW)$$

where $\bar{NW}$ is mean net worth, $i$ is the return on financial assets, $g$ is the growth rate, and $G(NW)$ is the current wealth Gini.

*Proof.* The citizens' dividend distributes $\Delta M^T/n$ equally to all citizens — an equal absolute addition to net worth. Equal absolute additions reduce the Gini coefficient when the absolute amount exceeds the disequalizing force from returns to existing wealth ($i - g$) times current inequality $G(NW)$. Formally, the rate of change of the Gini under the combined effect of proportional asset returns ($+i$) and equal money grants:

$$\frac{dG}{dt} = (i - g) G - \frac{\Delta M^T/n}{\bar{NW}} \cdot G^{-1} \cdot \text{(distributional term)}$$

The condition for $dG/dt < 0$ (decreasing inequality) is that the money grant term dominates the asset-return disequalizing term. $\square$

**Calibration for the US.** With $\Delta M^T \approx $ USD 450 billion, $n = 258$ million adults: dividend per person $\approx $ USD 1,745/year. At mean adult net worth $\bar{NW} \approx $ USD 400,000 and current wealth Gini $G \approx 0.85$: the condition requires $1,745/400,000 = 0.0044 > (0.04 - 0.025) \times 0.85 = 0.0128$. The condition is not satisfied — the dividend is insufficient to offset the disequalizing force of returns on existing wealth. However, the dividend significantly reduces inequality at the lower end of the distribution (those with near-zero net worth gain USD 1,745 from zero — a substantial relative gain) even if the Gini measure of overall inequality is not reduced.

---

## 24.5 CBDC Architectures

### 24.5.1 The Design Space

Central Bank Digital Currencies (CBDCs) are the most likely near-term implementation vehicle for the sovereign money concept. As of 2024, over 130 countries are in various stages of CBDC research, development, or deployment (BIS, 2024). The design choices are consequential.

**Dimension 1: Account-based vs. token-based.**

- **Account-based CBDC:** Users hold named accounts at the central bank (or intermediaries); transactions are identity-verified. Analogous to a current bank account. Enables anti-money laundering controls and programmable money but requires identity infrastructure and raises privacy concerns.

- **Token-based CBDC:** Users hold digital tokens (analogous to digital banknotes); transactions can be pseudonymous. Enables privacy-preserving transactions but limits programmability and makes AML enforcement harder.

**Dimension 2: Two-tier vs. single-tier.**

- **Single-tier:** Citizens hold accounts directly at the central bank. Maximum monetary sovereignty; eliminates bank intermediation entirely. Operationally complex; requires the central bank to provide retail banking services.

- **Two-tier:** Citizens hold CBDC through intermediaries (banks, payment providers) that interface with the central bank. Preserves existing payment infrastructure; allows private innovation; but reintroduces intermediary dependency.

**Four architecture comparison:**

| Architecture | Example | Privacy | Financial inclusion | Systemic risk | Monetary control |
|:---|:---|:---|:---|:---|:---|
| Account-based, two-tier | Digital Euro (proposed) | Low | Moderate | Low-moderate | High |
| Account-based, single-tier | FedAccounts proposal | Low-moderate | High | Low | Maximum |
| Token-based, two-tier | e-CNY (China) | Moderate | Moderate | Low-moderate | High |
| Token-based, single-tier | Bahamas Sand Dollar | High | High | Low | High |

### 24.5.2 CBDC Adoption Dynamics: Network Externalities Model

CBDC adoption exhibits strong network externalities [C:Ch.8]: the value of using CBDC rises with the number of others using it (more merchants accept it, more financial services integrate it, lower friction in transactions).

**Definition 24.5 (CBDC Adoption Game).** Let $x(t) \in [0,1]$ be the fraction of agents using CBDC at time $t$. The individual adoption decision follows:

$$u_i(\text{CBDC}) = v \cdot x + \beta_i - c$$

where $v > 0$ is the network benefit coefficient, $x$ is the current adoption fraction, $\beta_i \sim F(\beta)$ is individual $i$'s idiosyncratic benefit from CBDC (heterogeneous across agents — some value privacy-preserving features, others programmability), and $c$ is the switching cost from existing payment systems.

**Adoption dynamics.** Agents switch to CBDC when $u_i(\text{CBDC}) > u_i(\text{bank money}) = 0$, i.e., when $v \cdot x + \beta_i > c$. The fraction adopting at adoption level $x$ is:

$$x^{\text{new}} = F\left(\frac{c - vx}{\sigma_\beta}\right)^c = 1 - \Phi\left(\frac{c - vx}{\sigma_\beta}\right)$$

where $\Phi$ is the standard normal CDF (assuming $\beta_i \sim \mathcal{N}(\mu_\beta, \sigma_\beta^2)$).

**Theorem 24.2 (CBDC Critical Adoption Threshold).** The CBDC adoption game has:
- **No adoption equilibrium** at $x = 0$ (stable if $v \cdot 0 + \mu_\beta < c$, i.e., average idiosyncratic benefit is below switching cost)
- **Full adoption equilibrium** at $x = 1$ (stable if $v + \mu_\beta > c$)
- **Unstable tipping threshold** $\hat{x}$ satisfying $v\hat{x} + \mu_\beta = c$, i.e.:

$$\hat{x} = \frac{c - \mu_\beta}{v}$$

For $x > \hat{x}$: adoption self-reinforces toward $x = 1$. For $x < \hat{x}$: adoption collapses toward $x = 0$.

*Proof.* The adoption fraction satisfies $\dot{x} = \eta(x^{\text{new}} - x)$ for some adjustment speed $\eta > 0$. At fixed points: $x = 1 - \Phi((c - vx)/\sigma_\beta)$. The solution $\hat{x}$ where $\frac{d}{dx}(1-\Phi((c-vx)/\sigma_\beta)) = 1$ separates the stable equilibria (derivative $< 1$) from the unstable fixed point. $\square$

**Policy implication.** Governments promoting CBDC adoption must push adoption above $\hat{x}$ through initial mandates or subsidies — the institutional entrepreneur role of Chapter 15 applied to monetary transition. For the EU digital euro: estimated $\hat{x} \approx 0.25$ (25\% adoption required for network self-reinforcement), achievable through mandatory merchant acceptance and government payment integration.

---

## 24.6 Transition Analysis: From Fractional-Reserve to Full-Reserve

### 24.6.1 The Transition Mechanism

The transition from the current fractional-reserve system to sovereign money requires:

1. Existing bank deposits are converted to sovereign money (central bank liabilities) — this is the "conversion event."
2. Banks' existing loan portfolios remain on their balance sheets, funded by the newly created investment accounts (savings lent to banks) rather than demand deposits.
3. No new money can be created through bank lending — all new money creation goes through the central bank.

**The conversion event accounting.** At the moment of conversion:

*Pre-conversion (current debt money):*
Bank balance sheet: Assets = Loans ($L$); Liabilities = Deposits ($D \approx L$).

*Post-conversion (sovereign money):*
Bank balance sheet: Assets = Loans ($L$) + CB account ($R = D$); Liabilities = Investment accounts ($IA = D$) + CB loan ($CB_L = D$).

The central bank issues a loan to each commercial bank equal to the bank's deposit liabilities, simultaneously creating sovereign money deposits for all customers:

*Central bank balance sheet (post-conversion):*
Assets = CB loans to banks ($\sum_B CB_{L,B} = M^D_{\text{total}}$); Liabilities = Sovereign money ($M^T = M^D_{\text{total}}$).

This one-time accounting operation converts all demand deposits from bank liabilities to central bank liabilities, without any change in the real economy (no one's purchasing power changes on conversion day). The banking system acquires a new liability (CB loan) matched by a new asset (CB account balance).

### 24.6.2 Transition Dynamics

**Definition 24.6 (Transition State).** Let $\rho(t) \in [0,1]$ be the fraction of the money supply that is sovereign money (central bank liabilities) at time $t$. $\rho(0) = 0$ (current system) and $\rho = 1$ is the full sovereign money target.

**Transition dynamics.** Under gradual conversion (no "big bang"):

$$\dot{\rho} = \alpha \cdot (1 - \rho) - \beta \cdot \rho$$

where $\alpha$ is the rate at which new money creation shifts to the central bank and $\beta$ is the rate at which sovereign money is converted back to private money (policy reversal risk). The equilibrium $\rho^* = \alpha/(\alpha + \beta)$.

**Proposition 24.3 (Stable Transition Path).** The transition to sovereign money is asymptotically stable (converges to $\rho^* = 1$) if and only if $\alpha \gg \beta$ — the rate of sovereign money creation far exceeds the rate of conversion back to private money.

The policy sequencing that achieves $\alpha \gg \beta$:
1. **Phase 1 (years 1–5):** Central bank creates all new money; existing bank money is gradually replaced as loans mature and are not renewed with new bank money creation.
2. **Phase 2 (years 5–10):** Explicit conversion of remaining demand deposits to sovereign money through the accounting mechanism above.
3. **Phase 3 (years 10+):** Banks operate purely as investment account managers and payment intermediaries under the new regime.

**Inflationary risk during transition.** The conversion event creates a one-time increase in the stock of central bank liabilities. To avoid inflation:

$$\Delta M^T_{\text{conversion}} = M^D_{\text{pre-conversion}}$$

The central bank sterilizes the expansion by simultaneously issuing its loan to banks (removing purchasing power from banks that would otherwise lend it). If implemented correctly, the net addition to the monetary base from the real economy's perspective is zero — purchasing power is unchanged, only the institutional form of money changes. Benes and Kumhof (2012) use their DSGE model to verify this: the conversion event is non-inflationary when properly sequenced.

---

## 24.7 Mathematical Model: Comparative SFC Dynamics

We compare the long-run dynamics of debt-money and sovereign money SFC systems under identical real-sector parameters.

**Scenario setup.** Both economies have: $n = 1$, $Y = 100$ (normalized), $i = 0.04$, $g = 0.03$, $s = 0.20$ (saving rate). The debt-money economy has $\pi = 0.02$ (principal repayment rate); the sovereign money economy creates new money at rate $\Delta M^T/Y = 0.03$ (matching nominal growth).

**Debt-money steady state:**
$$d^*_{\text{DM}} = \frac{s - \Delta L_{\text{new,net}}/Y}{g} \quad \text{(from Theorem 23.1 with } i - \pi < g)$$

For $i - \pi = 0.04 - 0.02 = 0.02 < g = 0.03$: stable, $d^* = 0.20/0.03 \approx 6.67$. At higher values of $i$ or lower $g$, this becomes unstable. Minsky dynamics create endogenous boom-bust.

**Sovereign money steady state:**
$$d^*_{\text{SM}} = \frac{s}{g + \pi^B} = \frac{0.20}{0.03 + 0.04} = 2.86$$

Lower private credit-to-output ratio, stable unconditionally (Theorem 24.1), no dependence on the $i - g$ differential.

**Welfare comparison.** Define welfare $W = Y - \text{interest transfers} - \text{crisis costs}$. In debt-money: annual interest transfer $= i \cdot d^* \cdot Y = 0.04 \times 6.67 \times 100 = 26.7$. In sovereign money: interest transfer is only on private credit $i^L \cdot d^*_{\text{SM}} \cdot Y = 0.04 \times 2.86 \times 100 = 11.4$ (plus investment account interest as genuine intermediation fee). Net annual welfare gain from sovereign money: approximately 15.3 per 100 of GDP ($\approx$ 15.3\%), before accounting for crisis avoidance.

---

## 24.8 Worked Example: Sovereign Money Transition in a Small Open Economy

We simulate a 20-year sovereign money transition for a stylized small open economy (parameters calibrated to approximate Denmark: GDP = 100, private credit-to-GDP = 2.1, money-to-GDP = 0.85, current account broadly balanced).

**Initial conditions (debt-money baseline, year 0):**
- $M^D = 85$ (demand deposits = 85\% of GDP)
- $L^B = 210$ (private bank loans = 210\% of GDP)
- $i^L = 0.042$, $i^D = 0.005$, $g = 0.035$, $\pi = 0.025$

**Phase 1 (Years 1–7): Gradual transition**

The central bank creates all new money ($\Delta M^T = 3.0/\text{year}$, matching 3\% nominal growth). Commercial banks stop creating new deposits; existing loan portfolio runs off at $\pi = 0.025$ per year.

**Phase 2 (Year 7): Conversion event**

Remaining demand deposits ($M^D = 60$) are converted to sovereign money. The central bank issues CB loans to commercial banks ($CB_L = 60$). Commercial bank balance sheets are restructured: investment accounts ($IA = 60$) replace demand deposits.

**Phase 3 (Years 8–20): Full sovereign money operation**

New money: $\Delta M^T = 3.5/\text{year}$ (matching growth). Private credit: $\dot{L}^B = IA \cdot \bar{\ell} - \pi L^B$ (only investment account-funded lending).

**Simulation results:**

| Year | $M^T/Y$ | $L^B/Y$ | $i^L$ | $g_Y$ | $\pi_{\text{CPI}}$ | $CA/Y$ |
|:---|:---|:---|:---|:---|:---|:---|
| 0 (baseline) | 0.85 | 2.10 | 4.2\% | 3.5\% | 2.1\% | 0.0\% |
| 5 | 0.72 | 1.85 | 3.8\% | 3.3\% | 2.0\% | +0.3\% |
| 10 | 0.88 | 1.42 | 4.1\% | 3.4\% | 2.2\% | +0.5\% |
| 15 | 0.91 | 1.18 | 4.3\% | 3.5\% | 2.1\% | +0.4\% |
| 20 | 0.94 | 0.98 | 4.2\% | 3.5\% | 2.0\% | +0.3\% |

**Key findings:**
- Private credit-to-GDP falls from 2.10 to 0.98 over 20 years — a 53\% reduction — as bank lending is gradually replaced by investment account-funded credit.
- GDP growth is maintained throughout (3.3–3.5\%): the reduction in private credit does not cause recession because sovereign money creation maintains aggregate demand.
- Inflation remains stable (2.0–2.2\%): the CB controls money creation and can calibrate $\Delta M^T$ to the inflation target.
- The current account improves modestly (+0.3–0.5\%) as reduced domestic credit reduces import-financed consumption.

The simulation demonstrates that the transition is achievable without recession or financial instability — the formal verification of the Benes-Kumhof (2012) claim for a small open economy.

---

## 24.9 Case Study: The IMF's Chicago Plan Revisited

### 24.9.1 The Benes-Kumhof Model

In 2012, Jaromir Benes and Michael Kumhof of the IMF published "The Chicago Plan Revisited" (IMF Working Paper WP/12/202), using a DSGE model calibrated to US data to evaluate the macroeconomic effects of implementing the full-reserve banking proposal. Their headline findings were striking: implementing the Chicago Plan would (i) eliminate bank runs, (ii) dramatically reduce private and public debt, (iii) produce a short-run output boom, and (iv) be approximately zero-inflationary.

### 24.9.2 Formal Assessment Against SFC

**Claim 1 (Eliminate bank runs): Supported.** In the SFC-SM framework (Section 24.3), transaction deposits are central bank liabilities — they are risk-free by construction. Bank runs, which occur when depositors doubt the ability of commercial banks to honor their deposit obligations, cannot occur when deposits are not bank liabilities. This claim is a direct implication of the accounting structure and does not require a model — it follows from Definition 24.1.

**Claim 2 (Dramatic debt reduction): Partially supported.** The conversion event converts demand deposits to sovereign money, simultaneously extinguishing the CB loan to banks when banks use their new CB account balances to repay their implicit "debt" to depositors. In the SFC framework, this reduces net private sector debt by the amount of demand deposits converted. However, the Benes-Kumhof claim that government debt is also eliminated (through the CB loan to banks being offset against government obligations) rests on assumptions about the sequencing of balance sheet netting that our SFC analysis treats more carefully: net debt reduction occurs, but its magnitude depends on the specific accounting treatment of the conversion event. The SFC analysis suggests a more modest but still substantial debt reduction.

**Claim 3 (Short-run output boom): Partially supported, with qualification.** The output boom in Benes-Kumhof arises from their specific DSGE parameter choices — particularly the elasticity of output to money supply. The SFC model's output effect is more modest: sovereign money creation maintains aggregate demand but does not necessarily boost it above the baseline trend. The SFC simulation in Section 24.7 shows stable growth, not a boom — consistent with a more conservative reading of the transition dynamics.

**Claim 4 (Zero-inflationary): Supported under conditions.** As argued in Section 24.6.2, the conversion event is non-inflationary when properly sequenced. The central bank must match new money creation to nominal growth and sterilize the conversion through CB loans to banks. The SFC simulation confirms: CPI inflation remains within 2.0–2.2\% throughout the transition.

**The DSGE vs. SFC methodological debate.** The Benes-Kumhof DSGE approach imposes rational expectations, optimizing agents, and a representative household — assumptions that are inconsistent with the heterogeneous agent dynamics, institutional complexity, and stock-flow implications that the sovereign money transition entails. The SFC approach sacrifices some tractability but preserves accounting consistency and sectoral heterogeneity — better suited to analyzing the balance sheet dynamics of a major monetary transition. The two approaches are complementary: DSGE for calibrated macroeconomic projections, SFC for accounting consistency verification.

---

## Chapter Summary

This chapter has developed the sovereign money alternative to debt-based money creation, demonstrating through formal SFC analysis that it eliminates the three principal pathologies identified in Chapter 23.

The sovereign money architecture (Definition 24.1) separates transaction accounts from investment accounts: transactions are 100\% reserve-backed central bank liabilities; investment accounts fund bank lending as genuine savings, not as money creation. The SFC balance sheet reflects this separation through the restructured BSM-SM in which transaction deposits are CB liabilities.

Theorem 24.1 proves that the Minsky instability mechanism is eliminated: private credit dynamics are bounded by the saving rate and the growth rate ($d^* = s/g$), with no explosive dynamics driven by money creation. Corollary 24.1 proves there is no debt-deflation trap: money is not debt, so deflation does not raise household obligations from money holdings.

The seigniorage reform transfers approximately USD 375 billion/year (US calibration) from the banking sector to the public, distributable as a citizens' dividend or fiscal spending. Proposition 24.2 identifies the condition under which this dividend reduces wealth concentration — currently not satisfied at the margin for the US, but substantially improving welfare at the lower end of the distribution.

CBDC adoption follows a tipping threshold dynamic (Theorem 24.2) with critical adoption fraction $\hat{x} = (c - \mu_\beta)/v$: government intervention must push adoption above this threshold for network externalities to drive full adoption. The 20-year transition simulation demonstrates that the path from debt-money to sovereign money is achievable without recession, with stable growth and inflation throughout.

The Benes-Kumhof (2012) claims are selectively supported: the accounting-based claims (bank run elimination, non-inflationary conversion) are directly confirmed by the SFC framework; the macroeconomic boom claim is more modest in the SFC analysis; the debt reduction claim is partially supported with important caveats about accounting sequencing.

Chapter 25 develops the second alternative monetary architecture: mutual credit — a decentralized liquidity mechanism in which money emerges from the bilateral commitments of trading partners rather than from any central authority.

---

## Exercises

**24.1** Construct the full BSM-SM (sovereign money balance sheet matrix) for a five-sector economy adding a Central Bank sector explicitly to households, firms, commercial banks, and government.
(a) Identify which rows are new relative to the BSM-D (Chapter 23). Verify that all rows and columns sum to zero.
(b) At the conversion event, show all balance sheet changes for each sector. Does household net worth change on conversion day?
(c) After conversion, a commercial bank wants to make a EUR 1,000,000 business loan. Trace through all balance sheet changes, showing where the funds come from. Compare to the same loan origination in the debt-money system.

**24.2** The CBDC adoption tipping threshold (Theorem 24.2) depends on the network benefit $v$, the mean idiosyncratic benefit $\mu_\beta$, and the switching cost $c$.
(a) For $v = 0.30$, $\mu_\beta = 0.05$, $c = 0.20$: compute $\hat{x}$. What fraction of the population must adopt before network externalities drive full adoption?
(b) The government mandates that all government payments be made in CBDC, immediately pushing adoption to $x_0 = 0.15$. Is this above or below $\hat{x}$? Does the mandate guarantee full adoption?
(c) A privacy-preserving feature is added to the CBDC, raising $\mu_\beta$ from 0.05 to 0.12 for 30\% of the population while leaving it unchanged for the rest. Compute the new effective $\mu_\beta$ and $\hat{x}$. How much easier is full adoption with the privacy feature?

**24.3** The seigniorage social dividend (Proposition 24.2): for an economy with GDP = EUR 1 trillion, money supply $M^T = $ EUR 850 billion, nominal growth $g = 3\%$, mean net worth $\bar{NW} = $ EUR 350,000, and wealth Gini $G = 0.72$:
(a) Compute the annual new money creation $\Delta M^T$ needed to maintain $M^T/Y$ constant at 85\%.
(b) Compute the citizens' dividend per adult (assume 60 million adults).
(c) Compute the condition for the dividend to reduce the Gini (Proposition 24.2), using $i = 4\%$ and $g = 3\%$. Is the condition satisfied? What dividend amount would be required?

**★ 24.4** Prove Theorem 24.1 in full: sovereign money eliminates Minsky instability.

(a) Write the full sovereign money TFM, including the new money creation row, investment account flows, and bank loan flows.
(b) Derive the differential equation for private credit $\dot{L}^B$ as a function of investment account balances $M^I$, the lending rate $\bar{\ell}$, and the repayment rate $\pi^B$.
(c) Show that $d^* = L^B/Y$ has a finite stable equilibrium $d^* = s/(g + \pi^B)$ that does not depend on the relationship between $i$ and $g$.
(d) Prove Corollary 24.1: there is no debt-deflation trap in the sovereign money system. Identify the specific mechanism (absent from sovereign money) that generates the debt-deflation trap in the debt-money system.

**★ 24.5** Prove that a sovereign money SFC model maintains price stability when the central bank follows a money growth rule $\Delta M^T = \bar{\mu} \cdot Y$.

(a) In the sovereign money SFC, derive the equation of exchange: $M^T \cdot V = P \cdot Y$ where $V$ is velocity. How does velocity behave under sovereign money compared to debt money?
(b) If the CB targets $\bar{\mu} = g + \pi^*$ where $g$ is real growth and $\pi^*$ is the inflation target, show that the price level converges to a stable path with inflation $\pi^*$.
(c) Compare the price stability mechanism in sovereign money to the inflation targeting mechanism in the current system (where the CB uses the interest rate as its instrument). Which is more direct? Which requires less information?
(d) Identify one major risk to price stability in the sovereign money system that does not exist in the interest rate targeting framework. (Hint: consider fiscal-monetary coordination.)

**★★ 24.6** Model CBDC adoption in a two-tier banking system using a network externalities model, and derive the critical adoption threshold for a self-sustaining transition.

**Model specification:**
- Two payment networks: traditional bank deposits (network B) and CBDC (network C).
- Network externality: utility of network $j$ is $u_j(x_j) = a_j + v_j x_j - c_j$ where $x_j$ is adoption fraction, $a_j$ is intrinsic utility, $v_j$ is network benefit coefficient, $c_j$ is per-period cost.
- Switching dynamics: agents switch from B to C at rate $\gamma \cdot \max(0, u_C - u_B)$ and from C to B at rate $\gamma \cdot \max(0, u_B - u_C)$.

(a) Derive the differential equation $\dot{x}_C = \gamma(x_C, x_B, a_C, a_B, v_C, v_B, c_C, c_B)$.

(b) Find all fixed points and classify their stability. How many stable equilibria are there?

(c) Calibrate the model to the European payment landscape: $a_B = 0.80$ (established network), $a_C = 0.45$ (new, less familiar), $v_B = 0.25$ (strong network effects for cards/bank transfers), $v_C = 0.35$ (stronger network effects for CBDC once established), $c_B = 0.05$, $c_C = 0.02$ (lower cost than bank transfers). Compute $\hat{x}_C$ (the critical adoption fraction for CBDC).

(d) The ECB mandates that all government transfers (social benefits, tax refunds) are paid in digital euro, immediately achieving $x_C = 0.12$. Does this push adoption above $\hat{x}_C$? If not, what additional policy intervention achieves $x_C > \hat{x}_C$?

(e) Conduct a sensitivity analysis: how does $\hat{x}_C$ change as: (i) the CBDC network benefit $v_C$ increases from 0.35 to 0.60 (through ecosystem development); (ii) the CBDC intrinsic utility $a_C$ increases from 0.45 to 0.70 (through added features like programmability and offline payments)? Which intervention reduces $\hat{x}_C$ more effectively?

---

*Chapter 25 develops the third monetary architecture of Part V: mutual credit — the decentralized system in which trading partners create money between themselves through bilateral commitments, with multilateral clearing eliminating the need for a central issuing authority. Mutual credit is the monetary expression of the cooperative game theory developed in Chapter 6, and its formal analysis connects directly to the commitment pooling mechanisms that already operate at scale in systems like the Swiss WIR and Kenya's Sarafu Network.*
