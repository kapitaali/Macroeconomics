# Chapter 32 — Open Economy Macroeconomics: Exchange Rates and Global Imbalances

> *"The international monetary system is at the center of many of the most contentious issues in economic policy."*
> — Barry Eichengreen, *Globalizing Capital*, 1996

---

Chapter 21 introduced the foundations of open-economy macroeconomics — comparative advantage, the exchange rate, interest parity conditions, and the Dornbusch overshooting model. Chapter 26 examined the balance of payments from the sectoral perspective. This chapter synthesizes these threads into a more advanced treatment: the New Open Economy Macroeconomics (NOEM), the theory and evidence on global imbalances, the three-generation theory of currency crises, and the macroeconomic consequences of international financial integration. Together, these frameworks provide the analytical apparatus for understanding the most consequential events in international macroeconomics over the past three decades — the Asian financial crisis, the global savings glut, the European sovereign debt crisis, and the global financial cycle.

---

## 32.1 The New Open Economy Macroeconomics

The **New Open Economy Macroeconomics (NOEM)**, developed by Obstfeld and Rogoff (1995, 1996), embeds open-economy international economics into an explicit microfounded dynamic general equilibrium framework with monopolistic competition, nominal rigidities, and forward-looking agents. It overcomes the main limitations of the Mundell-Fleming IS-LM framework [Ch. 9]: ad hoc behavioral equations, static structure, and the absence of welfare analysis.

### The Two-Country Framework

Consider a world with a home country ($H$) and a foreign country ($F$), each producing a continuum of differentiated goods under monopolistic competition. Home households consume a Dixit-Stiglitz basket:

$$C_t = \left[a^{1/\rho}C_{H,t}^{(\rho-1)/\rho} + (1-a)^{1/\rho}C_{F,t}^{(\rho-1)/\rho}\right]^{\rho/(\rho-1)},$$

where $a \in (0,1)$ is the home expenditure share on domestic goods (**home bias**), $\rho > 0$ is the **Armington elasticity** of substitution between home and foreign goods, and the corresponding price index:

$$P_t = \left[a\, P_{H,t}^{1-\rho} + (1-a)(e_t P_{F,t}^*)^{1-\rho}\right]^{1/(1-\rho)}.$$

The real exchange rate $q_t = e_t P_t^*/P_t$ drives expenditure switching: when $q_t$ rises (domestic goods become cheaper relative to foreign), home demand shifts toward domestic goods and foreign demand shifts away from foreign goods. The strength of this expenditure-switching effect determines whether exchange rate depreciation improves the trade balance — and at what speed.

### Pricing-to-Market and Exchange Rate Pass-Through

A crucial feature of the NOEM is the treatment of pricing: do firms set prices in their own currency (**producer currency pricing, PCP**) or in the consumer's currency (**local currency pricing, LCP**)?

**Definition (Exchange Rate Pass-Through).** **Exchange rate pass-through (ERPT)** is the percentage change in import prices (in domestic currency) per 1% change in the nominal exchange rate:

$$ERPT = \frac{\partial \ln P_M}{\partial \ln e},$$

where $P_M$ is the domestic-currency import price. Under full PCP, $ERPT = 1$: a 10% depreciation immediately raises import prices by 10%. Under LCP, $ERPT = 0$ in the short run since prices are fixed in local currency; ERPT builds toward 1 only as contracts are renegotiated. Goldberg and Knetter (1997) find ERPT to U.S. import prices of approximately 0.6 in the short run, rising toward 1 over two to three years. The incomplete ERPT implies that the Marshall-Lerner condition may not hold immediately after depreciation — generating the J-curve: the trade balance initially worsens (import values rise but volumes don't adjust instantly) before improving as quantities adjust.

### The NOEM Three-Equation System

Log-linearizing around a symmetric steady state and combining Calvo pricing in both countries yields a system analogous to the closed-economy NK model but with cross-country linkages. The home dynamic IS equation is augmented by the real exchange rate:

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\pi_{t+1}] - r_t^n) + \alpha_{IS}\Delta\hat{q}_t,$$

where $\alpha_{IS}$ captures expenditure switching through the real exchange rate. The home NKPC:

$$\hat{\pi}_{H,t} = \beta\mathbb{E}_t[\hat{\pi}_{H,t+1}] + \kappa(\hat{x}_t + \psi_q\hat{q}_t) + u_t,$$

where $\psi_q$ captures the impact of the real exchange rate on marginal costs (imported inputs become cheaper when the currency appreciates). The exchange rate is an additional state variable coupling the IS and NKPC equations — its dynamics and the optimal policy response to exchange rate shocks are among the central policy questions in open-economy monetary theory.

### International Policy Coordination and the "Divine Coincidence"

In the closed-economy NK model, the central bank faces no trade-off under optimal policy: stabilizing inflation and the output gap simultaneously is possible (the "divine coincidence"). In the open-economy NOEM, cost-push shocks transmitted through the exchange rate break this coincidence: stabilizing domestic inflation may require allowing the output gap to widen, and stabilizing the output gap may require accepting exchange rate volatility that passes through to import prices. This creates genuine policy trade-offs that are absent in the closed-economy case and that motivate the exchange rate stabilization concerns of central banks, particularly in small open economies.

---

## 32.2 Global Imbalances: The Current Account Puzzle

From the late 1990s through the 2007 Global Financial Crisis, the world economy developed a striking pattern of current account imbalances: the United States ran persistent deficits averaging $-4$ to $-6\%$ of GDP, while China, Germany, Japan, and oil exporters ran large surpluses. Understanding these imbalances — their origin, their sustainability, and their connection to the financial crisis — is one of the central puzzles of international macroeconomics.

### The Exorbitant Privilege

A key feature making U.S. deficits more sustainable than those of other countries is the **exorbitant privilege**: the United States earns higher returns on its foreign assets than it pays on its foreign liabilities. U.S. external assets (FDI, equity stakes in foreign firms) earn risk-adjusted returns of approximately 6.8% per year; U.S. external liabilities (primarily U.S. Treasuries and safe assets held by foreign central banks) yield approximately 3.5% — an excess return of approximately 3.3 percentage points annually (Gourinchas and Rey, 2007). This asymmetry arises because the United States provides the world's reserve currency and the global supply of safe assets — a structural feature of the international monetary system.

**Definition (Exorbitant Privilege).** The **exorbitant privilege** refers to the advantage enjoyed by the issuer of the world's reserve currency: the ability to borrow at lower rates than others (since foreigners demand the reserve currency's safe assets regardless of yield), and to earn higher returns on outward investments than foreigners earn on their holdings of the reserve country's assets.

The **Triffin dilemma** (Triffin, 1960) identifies the fundamental tension in reserve currency provision: the country supplying the world's reserve currency must run persistent current account deficits to provide the world with sufficient dollar-denominated liquid assets. But persistent deficits eventually undermine confidence in the reserve currency's value. This tension plagued the Bretton Woods system (whose collapse in 1971 Triffin predicted in 1960) and remains a structural feature of the current dollar-centric system.

### Competing Explanations for Global Imbalances

Three complementary hypotheses address the causes of the pre-crisis imbalances.

**Global savings glut** (Bernanke, 2005): demographic aging in Japan and Germany increased saving; elevated precautionary saving in East Asian economies after the 1997 crisis; oil-exporter revenues recycled into U.S. Treasuries. Excess global saving depressed world real interest rates and flowed into the U.S. through the capital account, generating the current account deficit as a passive counterpart.

**Bretton Woods II** (Dooley, Folkerts-Landau, and Garber, 2003): Asian central banks deliberately undervalued their currencies to sustain export-led growth, accumulating foreign exchange reserves as a by-product. The arrangement was mutually beneficial: Asia received an undervalued exchange rate and access to U.S. markets; the United States received cheap imports and low financing costs for its twin deficits.

**Safe asset shortage** (Caballero, Farhi, and Gourinchas, 2008): emerging market economies demanded dollar-denominated safe assets (U.S. Treasuries, agency bonds) as self-insurance against sudden stops. This structural global excess demand for U.S. safe assets compressed U.S. interest rates, stimulated U.S. consumption, and generated the observed imbalances as a general equilibrium outcome — not as a policy failure but as the market outcome of a world with insufficient safe assets relative to demand.

---

## 32.3 Currency Crises: Three Generations of Models

The analysis of currency crises has evolved through three generations, each responding to puzzles the previous generation could not explain.

### First Generation: Krugman (1979)

The first-generation model derives a speculative attack as the rational consequence of an inconsistency between domestic monetary-fiscal policy and an exchange rate peg. A government pegs at $\bar{e}$ while financing a primary deficit through money creation at rate $\mu > 0$. Money creation causes inflation through PPP, eroding the real exchange rate. The central bank must sell reserves at rate $\dot{R} = -\mu M/P$ to defend the peg.

Reserves decline linearly: $R_t = R_0 - \int_0^t \mu(M_s/P_s)\,ds$. Let $T^*$ be the date when reserves reach zero under the pegged rate. Rational speculators attack *before* $T^*$: if they wait until reserves are exhausted, there is a discrete expected profit from the depreciation that would be immediately arbitraged away by earlier forward purchases of foreign currency. The unique equilibrium has the attack at date $T < T^*$, exactly when the shadow floating exchange rate — the exchange rate without reserve support — reaches $\bar{e}$.

**Key prediction**: crises occur only when fiscal-monetary fundamentals are inconsistent with the peg; fundamentals fully determine the timing of the attack.

### Second Generation: Self-Fulfilling Crises (Obstfeld, 1994)

The 1992 European Monetary System crisis — during which the UK pound and Italian lira were forced out of the ERM despite no obvious fiscal-monetary inconsistency — required a different framework. Second-generation models add **government optimization**: the government weighs the cost of defending the peg (higher interest rates, deeper recession) against the cost of abandoning it (credibility loss, inflationary consequences).

The circular structure: if markets expect devaluation, the government must raise interest rates to attract capital and defend the peg; higher rates deepen the recession, raising the political cost of defense; the government becomes more likely to devalue — validating the initial expectation. If markets expect no devaluation, interest rates need not rise, the economy is not pushed into recession, and defense is cheap — also self-validating. The result is **multiple equilibria**: the peg is simultaneously sustainable and attackable. A sunspot — a newspaper article, a speech by a foreign finance minister, contagion from another crisis — can coordinate expectations on the bad equilibrium and trigger a self-fulfilling crisis.

This model explains **crisis contagion**: after the UK and Italy devalued in 1992, investors attacked other ERM currencies (France, Spain, Portugal) not primarily because of their deteriorating fundamentals but because the precedent of successful attacks changed beliefs about the government's willingness to defend. The contagion was expectations-driven rather than fundamentals-driven.

### Third Generation: Balance Sheets and Sudden Stops

The 1997 Asian financial crisis confounded both previous models: Thailand, Korea, and Indonesia had low inflation, fiscal surpluses, and rapid growth before the crisis — yet experienced devastating currency collapses and severe recessions. Third-generation models focus on **balance sheet mismatches**.

**Currency mismatch**: firms and banks borrow in dollars but earn revenues in domestic currency. A depreciation raises the domestic-currency value of dollar debt, causing widespread insolvency even for previously profitable firms. The currency crisis triggers a banking crisis, and vice versa — the **twin crisis** structure.

**Maturity mismatch**: short-term foreign borrowing finances long-term domestic investment. When short-term loans are not rolled over — a sudden stop — firms cannot refinance profitable projects, forced liquidations begin, and fire-sale asset prices amplify the initial shock.

The formal structure (Krugman, 1999; Caballero and Krishnamurthy, 2001): a firm's collateral is $\kappa(q_t K_t + B_t^*/e_t)$, where $q_t K_t$ is the domestic asset value and $B_t^*/e_t$ is the foreign-currency debt measured in domestic units. A depreciation ($e_t$ rises) reduces the collateral ratio, forcing asset sales, depressing $q_t$, tightening the constraint further — a Fisherian debt-deflation in open-economy form. The key asymmetry: appreciation is not the mirror image of depreciation because balance sheet damage from depreciation is asymmetric (firms cannot "gain" from appreciation if they have no dollar assets).

---

## 32.4 The Global Financial Cycle and Its Implications for Policy

### The Global Financial Cycle

Miranda-Agrippino and Rey (2020) document a **global financial cycle**: a single global factor explains approximately 20% of the variance in risky asset prices and capital flows across countries. This global factor is driven primarily by U.S. monetary policy and risk appetite (measured by the VIX). When the Fed tightens, the VIX rises, global risk appetite falls, capital flows to emerging markets reverse, EME currencies depreciate, and credit conditions tighten everywhere simultaneously — regardless of local economic conditions.

This finding has a profound implication: the standard trilemma [Ch. 9] may be a **dilemma** in practice. Under floating exchange rates with free capital mobility, a country retains monetary autonomy in theory (it can set its own interest rate). But if global risk appetite — driven primarily by U.S. monetary policy — simultaneously affects domestic asset prices, credit conditions, and capital flows, then monetary autonomy is substantially constrained even under flexible rates. Exchange rate flexibility provides only partial insulation from global financial shocks; macroprudential policy and capital flow management tools are required to achieve genuine monetary independence.

### Financial Integration and Contagion

International financial integration generates benefits — risk sharing, more efficient capital allocation, technology transfer through FDI — but also creates channels for contagion: the transmission of financial distress across countries through channels unrelated to fundamentals.

**Definition (Financial Contagion).** **Financial contagion** is the transmission of a financial shock from one country to others through channels beyond direct trade links, common shocks, or shared fundamentals. Mechanisms include: "wake-up calls" (a crisis in one country alerts investors to similar vulnerabilities elsewhere); portfolio rebalancing (investors who suffer losses in one market sell assets in others to meet margin calls); and common lender effects (banks that are creditors to both countries reduce exposure to both simultaneously).

Forbes and Rigobon (2002) distinguish true contagion (increased correlation during crises, beyond fundamental links) from interdependence (persistently high correlations that reflect genuine economic integration). Many apparent contagion episodes turn out on closer examination to reflect interdependence. Nevertheless, the 2008 Global Financial Crisis demonstrated genuine contagion: the Lehman bankruptcy transmitted financial distress globally within days through interbank market freezes, collateral calls, and investor panic — affecting countries with minimal direct exposure to U.S. subprime mortgages.

---

## Chapter Summary

- The **NOEM** (Obstfeld-Rogoff 1995) embeds open-economy macro into a microfounded GE framework with monopolistic competition and Calvo pricing. The open-economy IS and NKPC are augmented by the real exchange rate $\hat{q}_t$; ERPT (0.6 for the U.S.) determines the speed of expenditure switching and the Marshall-Lerner J-curve.

- **Global imbalances** (U.S. deficits of $-4$ to $-6\%$ of GDP, 2000–07) are explained by the global savings glut, Bretton Woods II currency undervaluation, and the safe asset shortage — three complementary frameworks each capturing different aspects of the phenomenon. The **exorbitant privilege** (U.S. earns 3.3 pp more on foreign assets than it pays on liabilities) partially explains the sustainability of U.S. deficits; the **Triffin dilemma** identifies the structural tension in reserve currency provision.

- **Currency crises** have evolved through three theoretical generations: first-generation (Krugman) — fiscal-monetary inconsistency determines timing; second-generation (Obstfeld) — self-fulfilling multiple equilibria driven by government optimization; third-generation — balance sheet mismatches (currency and maturity) generating twin crises and Fisherian debt-deflation spirals.

- The **global financial cycle** (Miranda-Agrippino and Rey) — a single U.S.-driven factor explaining 20% of global asset price variance — implies that floating exchange rates provide only partial insulation from external shocks, converting the standard trilemma into an effective dilemma requiring macroprudential complements to achieve genuine monetary autonomy.

---

*Next: Chapter 33 — Economic Development: Growth, Poverty, and Structural Change*
