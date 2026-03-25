# Chapter 32 — Open Economy Macroeconomics: Exchange Rates and Global Imbalances

> *"The international monetary system is at the center of many of the most contentious issues in economic policy."*
> — Barry Eichengreen, *Globalizing Capital*, 1996

---

Chapter 21 introduced the foundations of open-economy macroeconomics — comparative advantage, the exchange rate, interest parity conditions, and the Dornbusch overshooting model. Chapter 26 examined the balance of payments and capital flows from the sectoral perspective. This chapter synthesizes these threads into an integrated treatment of open-economy macroeconomics at the level of modern DSGE analysis: the New Open Economy Macroeconomics (NOEM), the theory and evidence on global imbalances, the mechanics of currency crises, and the macroeconomic consequences of international financial integration.

---

## 32.1 The New Open Economy Macroeconomics

The **New Open Economy Macroeconomics (NOEM)**, developed by Obstfeld and Rogoff (1995, 1996), embeds open-economy international economics into an explicit microfounded dynamic general equilibrium framework with monopolistic competition, nominal rigidities, and forward-looking agents. It overcomes the main limitations of the Mundell–Fleming IS–LM framework: ad hoc behavioral equations, static structure, and the absence of welfare analysis.

### The Two-Country Framework

Consider a world with a home country ($H$) and a foreign country ($F$). Each produces a continuum of differentiated goods under monopolistic competition. Home households consume a Dixit–Stiglitz basket of home and foreign goods:

$$C_t = \left[a^{1/\rho}C_{H,t}^{(\rho-1)/\rho} + (1-a)^{1/\rho}C_{F,t}^{(\rho-1)/\rho}\right]^{\rho/(\rho-1)},$$

where $a \in (0,1)$ is the home country's expenditure share on domestic goods (home bias), $\rho > 0$ is the **Armington elasticity of substitution** between home and foreign goods, and $1-a$ is the expenditure share on foreign goods. The corresponding price index:

$$P_t = \left[a\, P_{H,t}^{1-\rho} + (1-a)(e_t P_{F,t}^*)^{1-\rho}\right]^{1/(1-\rho)},$$

where $P_{H,t}$ is the price of the home goods basket, $P_{F,t}^*$ is the foreign goods basket price in foreign currency, and $e_t$ is the nominal exchange rate. The real exchange rate $q_t = e_t P_t^*/P_t$ drives expenditure switching: when $q_t$ rises (domestic goods become cheaper relative to foreign), home demand shifts toward domestic goods and foreign demand shifts away from foreign goods.

### Pricing-to-Market and Exchange Rate Pass-Through

A crucial feature of the NOEM is the treatment of pricing: do firms set prices in their own currency (**producer currency pricing, PCP**) or in the consumer's currency (**local currency pricing, LCP**)?

**Definition (Exchange Rate Pass-Through).** **Exchange rate pass-through (ERPT)** is the percentage change in import prices (measured in domestic currency) for a 1% change in the exchange rate:

$$ERPT = \frac{\partial \ln P_M}{\partial \ln e},$$

where $P_M$ is the domestic-currency import price. Under full PCP (law of one price), $ERPT = 1$: a 10% depreciation raises domestic import prices by 10%. Under LCP (prices fixed in consumer currency), $ERPT = 0$ in the short run: a depreciation has no immediate effect on domestic import prices because prices are set in domestic currency.

Empirically, ERPT is incomplete and varies by country and time horizon. Goldberg and Knetter (1997) find ERPT to import prices of approximately 0.6 for the U.S. in the short run, rising toward 1 over 2–3 years. The incomplete ERPT has important policy implications: if a depreciation does not immediately reduce import prices, its expenditure-switching effect on trade volumes is delayed and muted, implying the Marshall–Lerner condition may not hold in the short run.

### The NOEM Three-Equation System

Log-linearizing the NOEM around a symmetric steady state and combining with Calvo pricing in both countries yields a system analogous to the closed-economy NK model but with cross-country linkages. The home dynamic IS equation:

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\pi_{t+1}] - r_t^n) + \alpha_{IS}\Delta\hat{q}_t,$$

where $\alpha_{IS}$ captures the expenditure-switching effect on aggregate demand via the real exchange rate. The home NKPC:

$$\hat{\pi}_{H,t} = \beta\mathbb{E}_t[\hat{\pi}_{H,t+1}] + \kappa(\hat{x}_t + \psi_q\hat{q}_t) + u_t,$$

where $\psi_q$ captures the impact of the real exchange rate on marginal costs (through imported input prices and the wage-price spiral). The exchange rate enters as an additional state variable, coupling the IS and NKPC equations and introducing additional dynamics not present in the closed-economy model.

---

## 32.2 Global Imbalances: The Current Account Puzzle

From the late 1990s through the 2007 Global Financial Crisis, a striking pattern of global current account imbalances developed: the United States ran persistent current account deficits averaging $-4$ to $-6\%$ of GDP, while China, Germany, Japan, and oil-exporting countries ran large surpluses. These imbalances — the "global savings glut" in Bernanke's (2005) phrase — were widely seen as unsustainable and potentially destabilizing.

### The Exorbitant Privilege

A key feature that makes U.S. imbalances more sustainable than those of other deficit countries is the **exorbitant privilege**: the United States earns higher returns on its foreign assets than it pays on its foreign liabilities. U.S. external assets (equity stakes in foreign firms, FDI, loans) earn high risk-adjusted returns; U.S. external liabilities (primarily U.S. Treasuries and safe assets held by foreign central banks) pay relatively low yields.

Gourinchas and Rey (2007) document this asymmetry: over 1952–2004, U.S. foreign assets earned approximately 6.8% per year while liabilities yielded approximately 3.5%, generating an annual excess return of approximately 3.3 percentage points. This privileged position arises because the U.S. provides the world's reserve currency and safe assets — a structural feature of the international monetary system rather than a short-term arbitrage opportunity.

The **Triffin dilemma** (Triffin, 1960): a country that supplies the world's reserve currency must run persistent current account deficits to provide the world with sufficient liquidity (dollar-denominated assets). But persistent deficits eventually undermine confidence in the reserve currency's value. This dilemma plagued the Bretton Woods system and eventually caused its collapse in 1971; it remains relevant to the current dollar-centric system.

### Explanations for Global Imbalances

Several complementary explanations for the U.S. deficits and Asian surpluses have been proposed:

**Global savings glut** (Bernanke, 2005): demographic aging in Japan and Germany increased saving; high precautionary saving in East Asian economies following the 1997 crisis; oil exporter revenues recycled into U.S. Treasuries. Excess global saving depressed world real interest rates and flowed into the U.S. current account deficit.

**Bretton Woods II** (Dooley, Folkerts-Landau, and Garber, 2003): Asian central banks (particularly China's PBOC) systematically undervalued their currencies to sustain export-led growth, accumulating foreign exchange reserves as a by-product. The arrangement was mutual: Asia got an undervalued exchange rate and export market access; the U.S. got cheap imports and low financing costs for its deficits.

**Safe asset shortage** (Caballero, Farhi, and Gourinchas, 2008): emerging markets demanded safe dollar-denominated assets (U.S. Treasuries) to self-insure against sudden stops, creating a structural global excess demand for U.S. safe assets. This demand suppressed U.S. interest rates, stimulated U.S. consumption, and generated the observed imbalances as a general equilibrium outcome.

---

## 32.3 Currency Crises: Three Generations of Models

The analysis of currency crises has evolved through three generations of models, each responding to the empirical puzzles that the previous generation could not explain.

### First-Generation: Krugman (1979)

The first-generation model derives a speculative attack as the rational consequence of an inconsistency between fiscal policy and an exchange rate peg. A government pegs its exchange rate at $\bar{e}$ while financing a primary deficit through money creation ($\dot{M}/M = \mu > 0$). The Purchasing Power Parity condition implies $\pi = \mu$: money creation causes inflation, eroding the real exchange rate. The central bank must sell foreign reserves to defend the peg, at rate $\dot{R} = -\mu M/P$.

Reserves follow: $R_t = R_0 - \int_0^t \mu(M_s/P_s)\,\mathrm{d}s$, declining linearly toward zero. Let $T^*$ be the date at which reserves would reach zero if the peg continued. Rational speculators anticipate the collapse and attack *before* $T^*$: if they wait until reserves actually reach zero, there will be a capital gain from the depreciation — but this profit opportunity implies that rational agents will buy foreign currency earlier, exhausting reserves sooner. The unique equilibrium has the attack occurring at date $T < T^*$, when the shadow floating-rate exchange rate (the exchange rate that would prevail without reserves) reaches $\bar{e}$.

The first-generation model makes a sharp prediction: crises occur only when fiscal-monetary policy is inconsistent with the peg; fundamentals fully determine the timing of the attack.

### Second-Generation: Multiple Equilibria (Obstfeld, 1994)

The 1992 European Monetary System crisis — in which the UK pound and Italian lira were forced out of the ERM — did not fit the first-generation model: neither country was obviously pursuing an unsustainable fiscal-monetary policy. The second-generation model adds **policy trade-offs**: the government actively chooses whether to defend the peg, weighing the cost of defense (higher interest rates, recession) against the cost of devaluation (credibility loss, inflationary consequences).

The circular structure generates multiple equilibria. Suppose markets expect a devaluation: to defend the peg, the government must raise interest rates, which deepens the recession and raises the political cost of defense, making devaluation more likely — validating the market's expectation. Alternatively, if markets believe the peg will hold: interest rates need not rise, the economy is not pushed into recession, and the peg is sustainable — also self-validating.

The key insight: crises can be **self-fulfilling**. A fundamentally solvent government can be driven to devaluation by a coordinated speculative attack, even though the peg could have been maintained had expectations been favorable. This explains why the 1992 ERM crisis spread from the UK and Italy to other ERM members (France, Spain, Portugal) through contagion — investors attacked other currencies not because of their fundamentals but because of the precedent set by the first devaluations.

### Third-Generation: Balance Sheets and Sudden Stops

The 1997 Asian financial crisis involved countries (Thailand, Korea, Indonesia) with seemingly sound fundamentals — low inflation, fiscal surpluses, high growth rates — that nonetheless experienced sudden stops and devastating currency collapses. The third-generation models focus on **balance sheet mismatches**:

**Currency mismatch**: firms and banks borrow in dollars but earn revenues in domestic currency. A depreciation raises the domestic-currency value of dollar-denominated debt, causing widespread insolvency even for previously solvent firms. This generates a "twin crisis" — currency crisis and banking crisis occurring simultaneously and mutually reinforcing.

**Maturity mismatch**: short-term foreign borrowing finances long-term domestic investment. When short-term loans are not rolled over (the sudden stop), firms cannot refinance, even profitable long-term projects must be liquidated, and the resulting fire sales amplify the initial shock.

The formal structure (Krugman, 1999): firm balance sheets determine collateral $\kappa\,(q_t K_t + B_t/e_t)$, where $K_t$ is capital valued at asset price $q_t$, $B_t/e_t$ is the dollar value of foreign borrowing, and $\kappa$ is the loan-to-value ratio. A depreciation (rise in $e_t$) reduces the dollar value of foreign debt relative to domestic assets, tightening the collateral constraint, forcing asset sales, depressing $q_t$, tightening the constraint further — a **Fisherian debt-deflation** in open-economy clothing.

---

## 32.4 Financial Integration and Contagion

International financial integration — the reduction of barriers to cross-border capital flows — generates benefits (risk diversification, efficient capital allocation, technology transfer) but also creates channels for **contagion**: the transmission of financial shocks across countries unrelated by fundamentals.

**Definition (Contagion).** **Financial contagion** is the transmission of a financial shock from one country to others through a channel that cannot be explained by direct trade links, common shocks, or shared fundamentals. It involves "wake-up calls" (a crisis in one country makes investors realize similar vulnerabilities exist elsewhere) and "portfolio rebalancing" (investors who lose wealth in one market sell assets in other markets to meet margin calls).

Forbes and Rigobon (2002) distinguish contagion (increased correlation during crises, beyond what existing fundamental links predict) from **interdependence** (high correlation that reflects genuine economic links). Many apparent contagion episodes turn out to be interdependence when properly measured; genuine pure contagion is harder to establish. Nevertheless, the 2008 Global Financial Crisis demonstrated compelling contagion: the failure of Lehman Brothers transmitted financial distress globally within days through interbank market freezes, collateral calls, and investor panic, regardless of countries' direct exposure to U.S. subprime assets.

The **Miranda-Agrippino and Rey (2020) global financial cycle**: a single global factor explains approximately 20% of variance in asset prices and capital flows across countries — a "Global Financial Cycle" driven primarily by U.S. monetary policy and risk appetite. When the Fed tightens, global risk appetite falls, capital flows reverse, EME currencies depreciate, and financial conditions tighten everywhere simultaneously. This evidence challenges the trilemma: even with floating exchange rates, countries cannot fully insulate themselves from U.S. monetary policy via the exchange rate, because the risk appetite channel bypasses the exchange rate.

---

*Next: Chapter 33 — Economic Development*
