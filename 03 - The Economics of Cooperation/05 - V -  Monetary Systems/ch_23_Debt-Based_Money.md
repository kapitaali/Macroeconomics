# Part V: Core Theory IV — Monetary Systems Beyond Debt
 
Parts II and III built the cooperative institutional architecture: the game-theoretic case for cooperation, the network structures that sustain it, the governance frameworks that formalize it, and the information mechanisms that make it function at scale. Part IV embedded that architecture in its biophysical substrate: the flows of materials and energy that economic activity requires, the ecological systems that regenerate those flows, and the thermodynamic limits that constrain what any economy can achieve regardless of its institutional form.
 
Both halves of the analysis are incomplete without a monetary theory. Money is not merely a medium of exchange layered on top of a real economy; it is the institutional mechanism through which economic coordination — cooperative or otherwise — is articulated and financed. The form that money takes, how it is created, by whom, and at what cost, shapes the distributional outcomes of the economy, the stability properties of the financial system, and the growth imperatives that drive the relationship between the economy and the biosphere.
 
This Part constructs a rigorous monetary theory for the cooperative-regenerative economy. We begin with the existing system — the debt-based money architecture in which commercial banks create money through lending — and subject it to formal analysis: proving its instability properties, its distributional consequences, and its structural incompatibility with ecological steady-state. We then develop three alternative monetary architectures — sovereign money, mutual credit, and resource-backed currency — analyzing their stability, equity, and ecological properties using the SFC and game-theoretic tools developed in earlier parts. The Part closes with a comparative stability analysis and a macroeconomic model of the non-debt monetary system.
 
The reader will find that the monetary analysis of Part V connects backward to the cooperative institutions of Parts II–III (mutual credit is a direct application of cooperative game theory) and forward to the synthesis of Part VI (the unified cooperative-regenerative model requires a monetary framework compatible with both ecological and institutional constraints).
 
---
 
# Chapter 23: A Critical Examination of Debt-Based Money — Instabilities and Inequities
 
> *"The process by which banks create money is so simple that the mind is repelled."*
> — John Kenneth Galbraith, *Money: Whence It Came, Where It Went* (1975)
 
> *"Stability — it turns out — is destabilizing."*
> — Hyman Minsky, *Stabilizing an Unstable Economy* (1986)
 
## Learning Objectives
 
By the end of this chapter, you should be able to:
 
1. Model the mechanics of commercial bank money creation through double-entry bookkeeping, demonstrate that loans create deposits rather than intermediating prior savings, and represent this formally in the SFC balance sheet framework.
2. Formalize Minsky's financial instability hypothesis as a dynamical systems result — deriving the explosive debt dynamics that emerge when the interest rate exceeds the growth rate, and specifying the formal conditions for a Minsky moment.
3. Prove that debt-based money creation systematically transfers purchasing power from borrowers to lenders at the moment of loan origination, and derive the cumulative wealth concentration this implies over time.
4. Derive the growth imperative formally: show that a debt-money economy requires perpetual nominal growth to remain solvent, and identify the precise mechanism through which this imperative drives ecological overshoot.
5. Replicate the pre-2008 US household balance sheet dynamics in the SFC model and identify the formal indicators of Minsky instability that were observable before the crisis.
6. Analyze Japan's 1990–2010 debt deflation as a formal SFC model outcome and compare its resolution path to the post-WWII debt reduction approach.
 
---
 
## 23.1 The Institutional Reality of Money Creation
 
The dominant theory taught in introductory macroeconomics — the money multiplier model — holds that commercial banks act as intermediaries: they accept deposits from savers and lend those funds to borrowers, multiplying the initial base money created by the central bank through the reserve requirement ratio. This model is intuitively appealing and pedagogically convenient. It is also, as the Bank of England stated explicitly in its 2014 Quarterly Bulletin, "not an accurate description of how money creation works in reality."
 
The reality is more remarkable and more consequential. When a commercial bank makes a loan, it does not lend out pre-existing deposits; it simultaneously creates a new deposit (the borrower's bank account credit) and a new loan (the borrower's liability to the bank). The bank creates money ex nihilo — from nothing — bounded not primarily by reserve requirements but by capital requirements, credit demand, and the bank's own risk assessment. The money supply is endogenous to the lending process, not exogenous to it.
 
This institutional fact — which has been accepted in the heterodox economics literature since at least Kaldor (1970) and is now acknowledged by central banks themselves — has profound implications for the stability, distributional properties, and ecological consequences of the modern monetary system. These implications are the subject of this chapter.
 
We proceed systematically: first establishing the accounting mechanics of money creation (Section 23.2), then deriving the instability properties that follow from them (Section 23.3), then proving the distributional consequences (Section 23.4), and finally connecting debt-based money to the ecological overshoot analyzed in Part IV (Section 23.5).
 
---
 
## 23.2 How Banks Create Money: The Formal Accounting
 
### 23.2.1 Double-Entry Bookkeeping and Money Creation
 
**Definition 23.1 (Double-Entry Balance Sheet).** A double-entry balance sheet satisfies at all times:
 
$$\text{Assets} = \text{Liabilities} + \text{Net Worth}$$
 
Every transaction affects at least two entries simultaneously, preserving this identity.
 
**The loan origination accounting.** When bank $B$ makes a loan of amount $L$ to household $H$:
 
*Bank $B$'s balance sheet:*
 
| Assets | Liabilities |
|:---|:---|
| Loan to $H$: $+L$ | Deposit of $H$: $+L$ |
 
*Household $H$'s balance sheet:*
 
| Assets | Liabilities |
|:---|:---|
| Deposit at $B$: $+L$ | Loan from $B$: $+L$ |
 
Both balance sheets balance; no pre-existing funds have moved anywhere. The bank has created a new asset (the loan) and a new liability (the deposit) simultaneously. The deposit is the money supply — it is the purchasing power $H$ can now spend. The money supply has increased by $L$ through an accounting operation, not through the transfer of pre-existing savings.
 
**Proposition 23.1 (Loans Create Deposits).** In a modern commercial banking system, the primary direction of causation in money creation is loans → deposits, not deposits → loans. Money is created when banks extend credit and destroyed when loans are repaid.
 
*Proof.* From the balance sheet above: the deposit (money) and the loan (credit) are created simultaneously in the same transaction. The deposit cannot exist prior to the loan — it is created by the loan. Repayment reverses both entries: the deposit is extinguished ($-L$ from $H$'s assets) and the loan is cancelled ($-L$ from $B$'s assets). The money supply contracts by $L$. $\square$
 
**The SFC representation.** In the SFC framework [M:Ch.28], the money creation transaction appears in the Transaction Flow Matrix as:
 
| | Households | Banks | $\Sigma$ |
|:---|:---|:---|:---|
| Loans | $-\Delta L_H$ | $+\Delta L_H$ | 0 |
| Deposits | $+\Delta D_H$ | $-\Delta D_H$ | 0 |
| Change in net worth | 0 | 0 | 0 |
 
Every row sums to zero (accounting identity). The household's net worth is unchanged at the moment of origination ($+\Delta D - \Delta L = 0$); it changes only as the loan generates income and as interest accumulates.
 
### 23.2.2 The Endogenous Money Supply
 
**Definition 23.2 (Endogenous Money).** The money supply $M$ is endogenous if it is determined by the lending decisions of commercial banks and the borrowing decisions of non-bank agents, rather than by the central bank's control of base money. Formally:
 
$$\Delta M = \Delta L_{\text{net}} = \Delta L_{\text{originated}} - \Delta L_{\text{repaid}}$$
 
The money supply expands when new loans exceed repayments and contracts when repayments exceed new loans. The central bank's base money (reserves and currency) influences the cost of money creation (through the interest rate) but not directly its quantity.
 
**The endogenous money school.** The endogenous money view was developed by Kaldor (1970), Moore (1988), Minsky (1986), and Werner (1997), among others. It is empirically confirmed by: (i) the Bank of England's 2014 acknowledgment [McLeay et al., 2014]; (ii) Bundesbank's explicit statement in its Monthly Report (April 2017) that "the money supply in the real world is created by commercial banks"; and (iii) the empirical work of Werner (2014) using detailed bank transaction data.
 
**Implications for the money multiplier.** The money multiplier model — in which base money $H$ is multiplied by $1/rr$ (where $rr$ is the reserve requirement) to give the money supply $M$ — fails empirically and theoretically. In practice: (i) most jurisdictions have no reserve requirements (UK, Canada, Sweden, Australia); (ii) where they exist, they are not binding constraints on lending — banks lend to creditworthy borrowers and then seek reserves; (iii) the empirical correlation between base money and broad money is weak and time-varying, inconsistent with a stable multiplier. The money multiplier is a useful pedagogical device but an incorrect description of monetary reality.
 
---
 
## 23.3 Minsky's Financial Instability Hypothesis: Formal Treatment
 
### 23.3.1 The Three Financing Regimes
 
Hyman Minsky (1986) identified three financing regimes that describe the relationship between a firm's (or household's) income and its debt service obligations:
 
**Definition 23.3 (Minsky Financing Regimes).**
 
- **Hedge financing:** Income $Y_i > $ interest payment $iD_i$ + principal repayment $\dot{D}_i^{\text{repay}}$. The unit can service all debt obligations from current income.
- **Speculative financing:** Income $Y_i > iD_i$ but $Y_i < iD_i + \dot{D}_i^{\text{repay}}$. The unit can pay interest but must roll over principal — it depends on continuous access to credit markets.
- **Ponzi financing:** Income $Y_i < iD_i$. The unit cannot even pay interest from income — it must borrow to pay interest, increasing debt continuously.
 
In Minsky's account, a financial cycle moves through these regimes over time: an expansion phase dominated by hedge financing gives way to increasing speculative and Ponzi financing as confidence grows, leverage increases, and debt servicing becomes more fragile. The transition from speculative to Ponzi financing — when income can no longer cover interest payments — is the Minsky moment.
 
### 23.3.2 Formal Dynamical Model
 
**Definition 23.4 (Minsky Debt Dynamics).** Let $D(t)$ be the aggregate private debt level and $Y(t)$ be nominal GDP. Define the debt ratio $d = D/Y$. The dynamics of $d$ are:
 
$$\dot{d} = \frac{\dot{D}}{Y} - \frac{D\dot{Y}}{Y^2} = \frac{\dot{D}}{Y} - d \cdot g$$
 
where $g = \dot{Y}/Y$ is the nominal GDP growth rate.
 
Debt dynamics: $\dot{D} = i \cdot D - \pi \cdot D + \Delta L_{\text{new}}$, where $i$ is the interest rate, $\pi$ is the principal repayment rate, and $\Delta L_{\text{new}}$ is new lending. At steady state of new lending (flow of new credit equals repayments): $\dot{D} = (i - \pi) D$.
 
$$\dot{d} = (i - \pi - g) \cdot d$$
 
**Theorem 23.1 (Minsky Instability Theorem).** The debt ratio $d = D/Y$ in a debt-based monetary economy:
 
1. Is stable ($d \to d^* < \infty$) if and only if $i - \pi < g$: the net interest rate (interest minus principal repayment) is below the nominal growth rate.
2. Grows explosively ($d \to \infty$) if $i - \pi > g$: the net interest rate exceeds the nominal growth rate.
3. Reaches a Minsky moment — the point at which the aggregate debt service ratio crosses the threshold of collective insolvency — at time $T^*$ satisfying:
 
$$T^* = \frac{1}{i - \pi - g} \ln\left(\frac{d_{\text{crisis}}}{d_0}\right)$$
 
where $d_{\text{crisis}}$ is the debt-to-GDP ratio at which aggregate Ponzi financing becomes dominant.
 
*Proof.* The ODE $\dot{d} = (i - \pi - g) d$ has solution $d(t) = d_0 e^{(i-\pi-g)t}$.
 
For $i - \pi < g$: the exponent is negative, so $d(t) \to 0$ (debt ratio stabilizes and falls). For $i - \pi > g$: the exponent is positive, so $d(t) \to \infty$ (explosive debt dynamics). The Minsky moment occurs when $d(T^*) = d_{\text{crisis}}$; solving $d_0 e^{(i-\pi-g)T^*} = d_{\text{crisis}}$ gives the expression above. $\square$
 
**Calibration to the pre-2008 US economy.** For the US, 2002–2007:
- $i \approx 0.055$ (effective interest rate on private debt, FRED data)
- $\pi \approx 0.030$ (principal repayment rate, estimated from flow of funds data)
- $g \approx 0.055$ (nominal GDP growth rate including housing price inflation)
 
In 2002–2005: $i - \pi - g \approx 0.055 - 0.030 - 0.055 = -0.030 < 0$. The debt ratio was (barely) stable.
 
In 2006–2007: $g$ fell to approximately 0.040 (housing price growth slowed) while $i$ rose (Fed funds rate raised to 5.25\%). $i - \pi - g \approx 0.055 - 0.030 - 0.040 = +0.015 > 0$. The system crossed into explosive debt dynamics approximately 12–18 months before the crisis — consistent with the timing of the Minsky moment in autumn 2007.
 
### 23.3.3 The Full SFC Model of Debt Instability
 
The dynamical equation $\dot{d} = (i-\pi-g)d$ captures the essential Minsky mechanism but abstracts from the sectoral balance sheet structure. We now embed it in the full SFC framework.
 
**Sectors:** Households ($H$), Firms ($F$), Banks ($B$), Government ($G$).
 
**Key balance sheet equations (steady-state SFC):**
 
$$\text{Household net worth: } NW_H = K_H + D_H - L_H$$
$$\text{Bank equity: } E_B = L_{\text{total}} - D_{\text{total}} - L_{\text{reserve}}$$
 
**Flow equations (Transaction Flow Matrix):**
 
Household income identity:
$$Y_H = W + \Pi_H + iD_H - iL_H - T$$
 
where $W$ = wages, $\Pi_H$ = distributed profits, $iD_H$ = interest received on deposits, $iL_H$ = interest paid on loans, $T$ = taxes.
 
Bank profit identity:
$$\Pi_B = i^L L_{\text{total}} - i^D D_{\text{total}} - \text{admin costs}$$
 
where $i^L > i^D$ (banks earn a spread between lending and deposit rates).
 
**The SFC Minsky condition.** Define the household debt service ratio:
 
$$\text{DSR}_H = \frac{i^L L_H + \pi L_H}{Y_H}$$
 
The SFC Minsky moment occurs when the household debt service ratio reaches a threshold $\overline{\text{DSR}}$ above which debt repayment crowds out consumption sufficiently to reduce $Y_H$, which further raises $\text{DSR}_H$ — a positive feedback leading to the debt deflation spiral.
 
**Formal Minsky condition:**
 
$$\text{DSR}_H > \overline{\text{DSR}} \iff \frac{(i^L + \pi) L_H}{Y_H} > \overline{\text{DSR}}$$
 
Calibrating to US 2007: $i^L \approx 0.065$, $\pi \approx 0.040$, $L_H/Y_H \approx 0.97$ (household debt reached 97\% of disposable income). $\text{DSR}_H \approx (0.065 + 0.040) \times 0.97 \approx 0.102$. For $\overline{\text{DSR}} \approx 0.11$ (estimated from historical crisis thresholds), the system was approaching the Minsky threshold by mid-2007.
 
---
 
## 23.4 Distributional Effects: Proof of Systematic Wealth Transfer
 
### 23.4.1 The First-Order Transfer
 
**Theorem 23.2 (Debt-Money Distributional Transfer).** At the moment of loan origination, the debt-money system transfers purchasing power from the borrowing sector to the banking sector equal to the present value of future interest payments:
 
$$\text{Transfer} = \text{PV(interest payments)} = \int_0^T i \cdot D e^{-rt} dt = \frac{iD}{r}\left(1 - e^{-rT}\right)$$
 
where $r$ is the discount rate and $T$ is the loan term. This transfer is implicit — it does not appear in any single transaction — but accumulates as a flow from debtors to creditors over the loan lifetime.
 
*Proof.* At origination, the borrower receives purchasing power $D$ (the deposit) and commits to repay $D + \int_0^T i \cdot D \, dt$ over the loan term. The interest payments $iD$ per period are a flow from the household sector to the banking sector with no corresponding productive service — they represent a transfer from debtors to creditors determined by the loan contract, not by market competition over productive output. The present value of these transfers is $\int_0^T i D e^{-rt} dt$ as stated. $\square$
 
**Cumulative effect.** At the aggregate level, with total private debt $D$ and interest rate $i$:
 
Annual interest transfer: $iD$
 
For the US in 2019: $D \approx $ USD 16 trillion (private non-financial sector debt), $i \approx 0.045$. Annual transfer from debtors to creditors: approximately USD 720 billion — approximately 3.3\% of GDP, flowing persistently from borrowers (lower wealth, lower income) to creditors (higher wealth, higher income).
 
### 23.4.2 The Compounding Wealth Concentration
 
**Proposition 23.2 (Compounding Concentration).** Under the debt-money system, wealth concentration — measured by the Gini coefficient $G(NW)$ of net worth — increases monotonically over time as long as the return on financial assets $i$ exceeds the economic growth rate $g$:
 
$$\frac{dG(NW)}{dt} > 0 \iff i > g$$
 
*Proof sketch.* Those who hold financial assets (bonds, deposits, equity claims on bank equity) receive the interest transfer of Theorem 23.2 as income. Those who hold debt pay it. Since financial asset holdings are concentrated among high-wealth households (Pareto-distributed: the top 10\% hold approximately 70–90\% of financial assets in most OECD countries), the interest flow systematically increases the net worth of high-wealth households relative to low-wealth households. The rate of net worth increase is proportional to $i - g$ (Piketty's $r > g$ result [C:Ch.1]); when $i > g$, the distribution becomes more concentrated over time. $\square$
 
The formal mechanism connecting debt-money to inequality is therefore: debt-money creation → systematic interest transfer from borrowers to creditors → concentration of financial assets among creditors → further interest income → accelerating concentration. The Piketty finding that $r > g$ has been characteristic of most economies over most of modern history is the empirical signature of this mechanism.
 
---
 
## 23.5 The Growth Imperative: Formal Derivation
 
### 23.5.1 Debt Requires Growth
 
**Theorem 23.3 (Debt-Money Growth Imperative).** A debt-money economy in which money is created through bank lending requires positive nominal GDP growth to prevent aggregate insolvency. Specifically, the aggregate debt service capacity condition requires:
 
$$g > i - \pi$$
 
where $g$ is nominal GDP growth, $i$ is the average interest rate on private debt, and $\pi$ is the principal repayment rate. Failure of this condition produces explosive debt dynamics and eventual systemic crisis (Theorem 23.1).
 
*Proof.* From Theorem 23.1: debt ratio stability requires $i - \pi < g$, equivalently $g > i - \pi$. Nominal GDP growth must exceed the net interest burden. In an economy where money creation ceases (as in a steady-state economy with zero nominal growth), $g \to 0$ and the stability condition requires $i = \pi$ — the interest rate must equal the principal repayment rate, implying zero economic profit for the banking sector. This is incompatible with a profit-driven commercial banking system. $\square$
 
**Corollary 23.1 (Ecological Incompatibility).** The debt-money growth imperative is incompatible with the Stewardship Condition $\dot{N} \geq 0$ under current technologies, since:
 
1. Nominal GDP growth requires (under current material intensity) positive material throughput growth.
2. Positive material throughput growth on a finite planet with finite regeneration rates eventually violates $\dot{N} \geq 0$ for at least some essential natural capital stocks.
3. Therefore, the debt-money system structurally drives ecological overshoot.
 
The formal qualification "under current technologies" is important: if material productivity (GDP/DMC) grows faster than GDP, material throughput can fall while GDP grows — the absolute decoupling condition. Chapter 17 showed that absolute decoupling has not been achieved at the global scale. Until it is, the debt-money growth imperative structurally drives material expansion incompatible with planetary boundaries.
 
### 23.5.2 The Debt-Deflation Trap
 
Irving Fisher (1933) identified the debt-deflation mechanism: in a highly indebted economy, a fall in prices increases the real value of debt, triggering asset sales, price further falls, and bank failures — a self-reinforcing spiral.
 
**Definition 23.5 (Debt Deflation Dynamics).** Let $P$ be the price level, $D$ the nominal debt, and $d_{\text{real}} = D/P$ the real debt burden. If prices fall ($\dot{P} < 0$):
 
$$\dot{d}_{\text{real}} = \frac{\dot{D}}{P} - \frac{D\dot{P}}{P^2} = \frac{\dot{D}}{P} + |\pi_{\text{deflation}}| \cdot d_{\text{real}}$$
 
where $\pi_{\text{deflation}} = |\dot{P}/P|$ is the deflation rate. Even with no new borrowing ($\dot{D} = 0$), the real debt burden grows at the rate of deflation — a positive feedback: falling prices → rising real debt → debt distress → asset sales → further price falls.
 
**The Fisher trap condition:**
 
$$\frac{d\dot{d}_{\text{real}}}{d_{\text{real}}} = \pi_{\text{deflation}} > 0 \quad \text{(always in deflation)}$$
 
Once deflation begins, the debt burden spiral is self-reinforcing. Escaping the trap requires either inflation (reducing real debt) or debt restructuring (reducing nominal debt directly). In the absence of either, the economy is trapped in a low-growth, high-real-debt equilibrium — the debt-deflation trap.
 
---
 
## 23.6 Mathematical Model: Full SFC of the Debt Economy
 
**Setup.** A three-sector SFC model with Households ($H$), Banks ($B$), and Government ($G$), extended with Minsky instability dynamics.
 
**Balance Sheet Matrix (BSM-D):**
 
| | $H$ | $B$ | $G$ | $\Sigma$ |
|:---|:---|:---|:---|:---|
| Produced capital $K$ | $+K_H$ | | | $+K_H$ |
| Deposits $D$ | $+D$ | $-D$ | | $0$ |
| Government bonds $BG$ | $+BG_H$ | $+BG_B$ | $-BG$ | $0$ |
| Loans $L$ | $-L_H$ | $+L_H$ | | $0$ |
| Reserve money $H^M$ | $+H^M_H$ | $-H^M_B$ | $+H^M_G$ | $0$ |
| Net worth | $NW_H$ | $NW_B$ | $NW_G$ | $0$ |
 
**Transaction Flow Matrix (TFM-D)** (selected rows):
 
| | $H$ | $B$ | $G$ | CB | $\Sigma$ |
|:---|:---|:---|:---|:---|:---|
| Consumption | $-C$ | | | | |
| Investment | $+I_H$ | | | | |
| Government spending | $+G^{\text{sp}}$ | | $-G^{\text{sp}}$ | | $0$ |
| Wages | $+W$ | $-W$ | | | $0$ |
| Interest on deposits | $+i^D D$ | $-i^D D$ | | | $0$ |
| Interest on loans | $-i^L L$ | $+i^L L$ | | | $0$ |
| Change in loans | $+\Delta L$ | $-\Delta L$ | | | $0$ |
| Change in deposits | $-\Delta D$ | $+\Delta D$ | | | $0$ |
 
**Minsky dynamics extension.** The household lending growth rate depends on the current debt service ratio:
 
$$\Delta L_H = \bar{\ell} \cdot Y_H \cdot (1 - \text{DSR}_H/\overline{\text{DSR}}) \cdot \mathbb{1}[\text{DSR}_H < \overline{\text{DSR}}]$$
 
Lending accelerates when DSR is well below the threshold and collapses to zero when DSR approaches it — the formal representation of the Minsky cycle within the SFC framework.
 
**Stability analysis.** The linearized dynamical system around the steady-state debt ratio $d^*$ has eigenvalue:
 
$$\lambda = i^L - \pi - g - \bar{\ell}(1 - \text{DSR}^*/\overline{\text{DSR}}) \cdot s$$
 
where $s = Y_H/L_H$ is the income-to-debt ratio at steady state. Stability requires $\lambda < 0$:
 
$$i^L - \pi - g < \bar{\ell}(1 - \text{DSR}^*/\overline{\text{DSR}}) \cdot s$$
 
The right-hand side is the stabilizing force from endogenous lending reduction as DSR rises. If $\bar{\ell}$ is large (lenders respond aggressively to deteriorating DSR) and $\text{DSR}^*/\overline{\text{DSR}}$ is small (current DSR is well below the threshold), the system is stable. If lending is insensitive to DSR ($\bar{\ell}$ small) or DSR is already near the threshold, the system is unstable.
 
---
 
## 23.7 Worked Example: Pre-2008 US Household Balance Sheet Dynamics
 
We replicate the key dynamics of the US household sector's balance sheet from 2002 to 2009 using the SFC-Minsky model, demonstrating the formal approach to crisis detection.
 
### 23.7.1 Data and Calibration
 
**Key variables (Federal Reserve Flow of Funds, 2002–2009):**
 
| Year | $L_H$ (USD tn) | $Y_H$ (USD tn) | $d_H = L/Y$ | $i^L$ | $\text{DSR}_H$ | $g_{\text{nominal}}$ |
|:---|:---|:---|:---|:---|:---|:---|
| 2002 | 8.3 | 8.9 | 0.93 | 0.059 | 0.097 | 0.051 |
| 2003 | 9.4 | 9.2 | 1.02 | 0.055 | 0.101 | 0.048 |
| 2004 | 10.5 | 9.7 | 1.08 | 0.054 | 0.105 | 0.067 |
| 2005 | 11.9 | 10.3 | 1.16 | 0.057 | 0.109 | 0.066 |
| 2006 | 13.0 | 10.8 | 1.20 | 0.063 | 0.113 | 0.059 |
| 2007 | 13.8 | 11.1 | 1.24 | 0.065 | 0.117 | 0.044 |
| 2008 | 13.8 | 11.0 | 1.25 | 0.062 | 0.110 | 0.028 |
| 2009 | 13.4 | 10.5 | 1.28 | 0.055 | 0.095 | $-0.019$ |
 
**Threshold identification.** Calibrating $\overline{\text{DSR}} = 0.117$ (the 2007 peak): the DSR crossed the Minsky threshold in late 2007, consistent with the onset of the financial crisis in August 2007 (BNP Paribas fund suspensions) and the formal Minsky moment in September 2008 (Lehman collapse).
 
### 23.7.2 Model Validation
 
**Minsky instability check.** For 2004–2006: $i^L - \pi - g \approx 0.054 - 0.030 - 0.066 = -0.042 < 0$ (stable). For 2007: $i^L - \pi - g \approx 0.065 - 0.030 - 0.044 = +0.009 > 0$ (unstable). The model correctly identifies the transition from stable to unstable debt dynamics in 2007 — approximately 12 months before the formal crisis.
 
**Distributional transfer.** Annual interest transfer from households to banks, 2007: $i^L \times L_H = 0.065 \times 13.8 = $ USD 897 billion — approximately 8\% of household disposable income, paid to a banking sector capturing approximately 30\% of US corporate profits.
 
**The formal early warning.** The DSR-to-threshold ratio $\text{DSR}_H/\overline{\text{DSR}}$ reached 0.98 by end-2007 — 98\% of the critical threshold. In the SFC-Minsky model, this corresponds to $\lambda > 0$ (unstable) and implies that any further deterioration in $g$ or rise in $i^L$ would trigger the debt-deflation spiral. The 2008 credit tightening (raising $i^L$ for new borrowers) and the recession (reducing $g$) delivered both simultaneously.
 
---
 
## 23.8 Case Study: Japan's Zombie Economy, 1990–2010
 
### 23.8.1 The Balance Sheet Recession
 
Japan's real estate and equity bubble collapsed in 1990–1991, destroying approximately USD 10 trillion in asset values — roughly twice Japan's annual GDP. The aftermath was the "lost decade" (subsequently two lost decades): GDP growth averaged less than 1\% per year from 1991 to 2010, despite near-zero interest rates and aggressive fiscal stimulus.
 
Richard Koo (2008) diagnosed the mechanism as a "balance sheet recession": firms that had borrowed heavily to finance real estate investments during the bubble found their asset values far below their debt levels — negative equity. Their response was to prioritize debt repayment over investment, regardless of interest rates. Even at zero interest rates, firms refused to borrow because their primary objective was debt minimization, not profit maximization.
 
**Formal representation.** In the SFC framework, the balance sheet recession is the firm sector version of the Minsky Ponzi condition:
 
$$\text{Firm equity: } NW_F = K_F + RE_F - L_F < 0$$
 
With negative equity ($NW_F < 0$), firms face the constraint:
 
$$\dot{L}_F^{\text{new}} = 0, \quad \dot{L}_F^{\text{repay}} = \max\left(0, \frac{\Pi_F}{i^L + \pi}\right)$$
 
Firms use all profits to repay debt rather than invest — the economic equivalent of Ponzi financing in reverse: firms reduce leverage as fast as income permits, creating a demand shortage in the macroeconomy.
 
### 23.8.2 Why Conventional Policy Failed
 
**Monetary policy failure.** The Bank of Japan reduced the policy rate to zero (ZIRP) by 1999 and implemented quantitative easing from 2001. Neither was effective because the transmission mechanism was broken: firms in balance sheet recession do not borrow regardless of the interest rate. The Minsky model predicts this: $\bar{\ell} \to 0$ (lending demand collapses) regardless of $i^L$ when $NW_F < 0$. QE expands the monetary base but does not create new lending — the money sits as excess reserves.
 
**Fiscal policy partial effectiveness.** Government stimulus spending ($G^{\text{sp}} > 0$) did sustain GDP above the otherwise-depressed level — but raised government debt from approximately 60\% of GDP in 1990 to approximately 240\% by 2015. The private sector's debt reduction was largely offset by government debt accumulation — a sectoral balance shift (private deleveraging → public leveraging) consistent with the SFC sectoral balances identity.
 
### 23.8.3 Comparison with Post-WWII Debt Reduction
 
The contrast with post-WWII debt reduction is instructive. US government debt reached 118\% of GDP in 1946; it fell to 31\% by 1974 — a 87 percentage point reduction in 28 years, without a balance sheet recession. How?
 
**The post-WWII mechanism.** Five factors drove debt reduction:
1. **Inflation** ($\pi_{\text{price}} \approx 4\%$): reduced real value of fixed nominal debt.
2. **High growth** ($g \approx 4\%$ real): expanded the GDP denominator.
3. **Financial repression**: interest rates held below inflation by regulation, so $i < g$ (stable debt dynamics).
4. **Progressive taxation**: high top marginal rates limited wealth concentration and increased government revenue.
5. **No asset price collapse**: the debt was government debt; the private sector was not balance-sheet-constrained.
 
In the Minsky framework: the post-WWII path maintained $i - \pi_{\text{price}} < g$ through deliberate policy. Japan's 1990–2010 path violated this condition through the combination of deflation (raising real debt), low growth, and zero-but-positive nominal rates.
 
**The formal lesson.** The SFC-Minsky model implies that escaping a balance sheet recession requires either: (a) restoring $i - \pi < g$ through inflation or growth, or (b) directly reducing nominal debt through restructuring, debt jubilee, or debt monetization. Japan chose neither decisively — it sustained nominal debt levels through evergreening (zombie lending by banks to insolvent firms) while monetary policy failed to generate inflation. The result was two decades of stagnation, formally predicted by the Minsky dynamics.
 
---
 
## Chapter Summary
 
This chapter has subjected the debt-based monetary system to formal analysis, establishing three principal results.
 
First, the accounting mechanics of money creation: loans create deposits rather than intermediating prior savings (Proposition 23.1). The money supply is endogenous to the lending process, determined by the lending decisions of commercial banks and the borrowing decisions of non-bank agents. This is not merely a theoretical claim — it is acknowledged by the Bank of England, the Bundesbank, and the Federal Reserve in their official publications.
 
Second, the Minsky instability theorem: the debt ratio $d = D/Y$ is stable if and only if $i - \pi < g$ (Theorem 23.1). When the net interest rate exceeds the nominal growth rate, debt dynamics become explosive, converging to a Minsky moment at a predictable time $T^*$. The pre-2008 US economy crossed into the unstable regime in 2007, approximately 12 months before the crisis — a timing the formal model correctly retrodicts.
 
Third, the distributional proof: debt-money creation systematically transfers purchasing power from debtors to creditors through interest payments (Theorem 23.2), and this transfer compounds wealth concentration whenever $i > g$ (Proposition 23.2) — the monetary mechanism underlying Piketty's $r > g$ inequality dynamic.
 
The growth imperative follows directly (Theorem 23.3): a debt-money economy requires $g > i - \pi$ to avoid explosive debt dynamics, and this imperative is structurally incompatible with ecological steady-state under current material intensities (Corollary 23.1). The debt-money system drives ecological overshoot not through any deliberate choice, but through its structural requirement for continuous nominal expansion.
 
Japan's 1990–2010 balance sheet recession illustrates the debt-deflation trap: private sector deleveraging collapses investment and GDP, monetary policy fails, and the economy is trapped in a low-growth, high-real-debt equilibrium until debt is reduced through inflation, restructuring, or government absorption.
 
Chapters 24–28 develop alternative monetary architectures that eliminate or substantially reduce these pathologies: sovereign money (Chapter 24), mutual credit (Chapter 25), resource-backed currency (Chapter 26), demurrage systems (Chapter 27), and the full macroeconomic comparison of non-debt monetary alternatives (Chapter 28).
 
---
 
## Exercises
 
**23.1** Model a bank balance sheet before and after a EUR 100,000 mortgage loan origination. Show all entries for both the bank and the household. Verify that the row-sum and column-sum SFC identities are satisfied.
(a) What happens to broad money ($M2$) when the loan is originated? When it is repaid?
(b) If the household immediately spends the loan proceeds on a house (paid to the seller, who deposits at a different bank), trace the interbank settlement. Does the aggregate money supply change?
(c) The central bank raises the reserve requirement from 0\% to 10\%. Does this affect the bank's ability to originate the mortgage? Explain using the endogenous money framework.
 
**23.2** Apply the Minsky instability theorem (Theorem 23.1) to the following economies:
(a) A developed economy with $i = 4\%$, $\pi = 3\%$, $g = 3\%$. Is the debt ratio stable or unstable?
(b) The same economy after a recession reduces $g$ to 1\%. Compute the new stability condition. How long until the debt ratio doubles from its pre-recession level $d_0 = 0.90$?
(c) If the central bank raises $i$ to 5\% to combat inflation while $g = 2\%$, what is the Minsky moment timing $T^*$ if $d_{\text{crisis}} = 1.50$ and $d_0 = 0.95$?
 
**23.3** The debt-money distributional transfer (Theorem 23.2): for a USD 500,000 mortgage at $i = 5\%$ over 30 years, with discount rate $r = 3\%$:
(a) Compute the total present value of interest payments transferred from the borrower to the bank.
(b) Compare this to the principal. What fraction of total payments is interest?
(c) If the household earns the median US income of USD 70,000/year, what fraction of lifetime income is transferred to the financial sector through this single mortgage?
 
**★ 23.4** Derive the formal condition for a Minsky moment in the full SFC-Minsky model.
 
(a) Write the complete SFC balance sheet matrix (BSM-D) for the four-sector economy (households, firms, banks, government) with all asset categories. Verify row sums are zero.
(b) Derive the differential equation for the household debt-service ratio $\dot{\text{DSR}}_H$ as a function of $i^L$, $g_Y$, $g_L$, and the lending growth rate $\bar{\ell}$.
(c) Find the fixed points of this system. Show that the Minsky moment corresponds to the DSR crossing the threshold $\overline{\text{DSR}}$.
(d) Prove that the system is locally unstable at the Minsky threshold: the Jacobian eigenvalue is positive when $\text{DSR}_H = \overline{\text{DSR}}$. Interpret this economically — why does the system accelerate rather than stabilize at the threshold?
 
**★ 23.5** Analyze the ecological growth imperative (Corollary 23.1) formally.
 
(a) Define absolute decoupling formally: the condition under which GDP growth $g > 0$ is compatible with $\dot{E} \leq 0$ (falling material/energy throughput). Express the decoupling condition in terms of material productivity growth $g_{\text{MP}} = g_{\text{GDP}} - g_{\text{DMC}}$.
(b) Using Eurostat data, has any G7 economy achieved absolute decoupling on both carbon and material dimensions over any 10-year period? Compute $g_{\text{MP}}$ for each G7 country for 2010–2020.
(c) If absolute decoupling requires $g_{\text{MP}} > g$ (material productivity must grow faster than GDP), and the historical maximum sustained $g_{\text{MP}}$ is approximately 2\% per year, what is the maximum sustainable GDP growth rate compatible with absolute decoupling?
(d) Compare this maximum growth rate to the Minsky stability condition $g > i - \pi$. For a typical advanced economy with $i = 3.5\%$ and $\pi = 2.5\%$, is the decoupling-compatible growth rate sufficient for debt stability? What does this imply about the structural compatibility of debt-money and ecological sustainability?
 
**★★ 23.6** Replicate Godley and Lavoie's SIM model, extend it with Minsky instability, and simulate the 2007–09 crisis dynamics.
 
**Model specification:**
- Baseline: Godley-Lavoie SIM model [M:Ch.28] with households, firms, government, and banks.
- Extension 1: Add endogenous lending growth $\Delta L_H = \bar{\ell} \cdot Y_H \cdot (1 - \text{DSR}_H/\overline{\text{DSR}})$ with $\bar{\ell} = 0.15$, $\overline{\text{DSR}} = 0.117$.
- Extension 2: Add debt-deflation feedback: if $\text{DSR}_H > \overline{\text{DSR}}$, household consumption falls by $\mu \cdot (\text{DSR}_H - \overline{\text{DSR}}) \cdot Y_H$, reducing $g$ endogenously.
- Calibrate to US 2002 initial conditions from Section 23.7.
 
(a) Simulate the baseline model (no Minsky extension) from 2002 to 2012. Report GDP growth, debt ratio, and DSR. Does the model exhibit crisis dynamics?
 
(b) Add Extension 1. At what year does the simulated DSR reach $\overline{\text{DSR}}$? Compare to the actual 2007 crisis.
 
(c) Add Extension 2. After the DSR crosses the threshold, how does the debt-deflation feedback affect the simulated GDP trajectory? Compare to actual US GDP 2008–2010.
 
(d) Test two policy responses: (i) government fiscal stimulus of 5\% of GDP for 2 years; (ii) direct debt restructuring (reduce $L_H$ by 20\% in the crisis year). Which produces faster recovery in the model? Which better replicates the actual US experience?
 
---
 
*Chapter 24 develops the sovereign money alternative: a monetary architecture in which money creation is a public function of the central bank, commercial banks serve as payment intermediaries rather than money creators, and the Minsky instability mechanism is formally eliminated. We assess its stability properties, its distributional consequences, and the practical transition pathway from the existing debt-based system.*
