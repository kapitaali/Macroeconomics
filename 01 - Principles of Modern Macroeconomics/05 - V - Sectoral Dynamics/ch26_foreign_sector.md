# Chapter 26 — The Foreign Sector: Balance of Payments and Capital Flows

> *"The balance of payments always balances."*
> — Paul Krugman and Maurice Obstfeld, *International Economics*, 2018

---

Every open economy maintains continuous economic relationships with the rest of the world: it exports and imports goods and services, borrows and lends financial capital, sends and receives migrant workers, and holds (and is held by) foreign financial assets and liabilities. The accounting framework that captures all these relationships is the **balance of payments**. Introduced in Chapter 4 as an accounting identity, the balance of payments becomes in this chapter a dynamic analytical tool: by understanding what drives current account balances, what determines capital flows, and how exchange rate regimes shape the interaction between external and internal balance, we can analyze some of the most important macroeconomic questions of the modern era — from the "global imbalances" of the 2000s to the sudden stops that have triggered emerging market crises, to the macroeconomic consequences of the eurozone's institutional design.

---

## 26.1 The Current Account: Intertemporal Approach

The saving-investment identity from Chapter 4 connects the current account to domestic saving and investment:

$$CA_t = S_t - I_t + (T_t - G_t) = (S_H - I_H) + (S_F - I_F) + S_G,$$

where the last expression uses the sectoral decomposition from Chapter 4. This identity has a powerful implication: the current account is not primarily a trade phenomenon but a saving-investment phenomenon. A country runs a current account deficit not because it produces poorly or trades unfairly, but because it invests more than it saves — it is borrowing from the rest of the world to finance domestic investment (or to consume more than it produces).

The **intertemporal approach** to the current account (Sachs, 1981; Obstfeld and Rogoff, 1994) treats the current account as the optimal response of a forward-looking representative household and firm to temporary versus permanent shocks to income and investment opportunities.

Consider a small open economy with world interest rate $r^*$. The representative household maximizes lifetime utility subject to the intertemporal budget constraint:

$$\sum_{t=0}^\infty\frac{CA_t}{(1+r^*)^t} = \sum_{t=0}^\infty\frac{CA_t}{(1+r^*)^t} = 0,$$

or equivalently, the present value of the current account balances must equal zero (the economy cannot run a Ponzi scheme on the rest of the world). The **current account as intertemporal trade**: just as agents trade goods across markets at a given date, they trade goods across time — running surpluses (lending) when endowments are high and deficits (borrowing) when endowments are low.

**Transitory versus permanent shocks.** Consider a negative income shock. If it is temporary ($y_t$ falls for one period but $y_{t+1}$ returns to normal), the household borrows internationally to smooth consumption, running a current account deficit: $CA_t = S_t - I_t = Y_t - C_t^* - I_t \approx Y_t - \bar{C} - I_t < 0$, where $\bar{C}$ is the smoothed permanent-income consumption level. If the shock is permanent, the household reduces consumption to match the lower permanent income: $\Delta C_t \approx \Delta Y_t$ and $CA_t \approx 0$. This sharp prediction — temporary shocks generate current account deficits; permanent shocks do not — is broadly consistent with evidence from natural disasters, commodity price shocks, and terms-of-trade changes.

---

## 26.2 Capital Flows: Determinants and Volatility

**Definition (Capital Flow).** A **capital flow** is a cross-border transaction in financial assets: foreign direct investment (FDI, establishing or acquiring a controlling stake in a foreign enterprise), portfolio investment (buying foreign stocks or bonds below controlling threshold), bank lending, and other investment. The **financial account** of the balance of payments records net capital flows: a positive financial account balance means the country is a net recipient of capital inflows (foreigners are accumulating claims on the domestic economy).

Capital flows are driven by three broad forces.

**Return differentials.** Capital moves toward higher risk-adjusted returns. The uncovered interest parity condition (Chapter 21) is the benchmark: $i_t = i_t^* + \mathbb{E}_t[\hat{e}_{t+1}]$. Deviations from UIP — when domestic returns exceed foreign returns even after expected depreciation — generate capital inflows until the differential closes through either rate adjustment or exchange rate movement.

**Risk diversification.** International capital flows allow investors to diversify idiosyncratic country risk. In a frictionless world (no transaction costs, no political risk, no information asymmetry), optimal portfolios would hold each country's assets proportional to its share of world market capitalization. In practice, strong **home bias** is observed: U.S. investors hold far more U.S. equities than world market weights imply. The home bias reflects transaction costs, information asymmetries, and the desire to hedge domestic labor income risk.

**Push and pull factors.** Capital flows to emerging markets are driven by both "pull" factors (domestic growth prospects, commodity booms, institutional reforms) and "push" factors (low interest rates in advanced economies that make carry trades to EMEs attractive). The dominance of push factors — documented by Miranda-Agrippino and Rey (2020) — means that EME financial conditions are largely determined by the global financial cycle, which is itself driven by U.S. monetary policy. When the Fed tightens, global risk appetite falls, capital flows out of EMEs, currencies depreciate, and financial conditions tighten everywhere simultaneously.

---

## 26.3 Sudden Stops and External Crises

**Definition (Sudden Stop).** A **sudden stop** is an abrupt reversal of capital inflows to an economy — a sharp deterioration in the financial account that forces a rapid current account adjustment. The term was introduced by Calvo (1998) to describe the Mexican crisis of 1994–95, in which a sudden withdrawal of foreign capital forced a sharp reversal of the current account deficit from $-7\%$ to $+1\%$ of GDP in less than a year.

The macroeconomic consequences of a sudden stop are severe because the forced current account adjustment requires either: (i) a large reduction in domestic absorption ($C + I + G$ must fall to match the fall in external financing); (ii) a large exchange rate depreciation (to improve competitiveness and shift demand toward exports); or (iii) some combination. In the presence of **original sin** (Chapter 35 and discussed below), the depreciation worsens rather than ameliorates the crisis by raising the domestic-currency value of foreign-currency debt.

The Calvo (1998) model of a sudden stop: Suppose the economy faces a binding **collateral constraint**:

$$D_t \leq \kappa\,\theta_t\, K_t,$$

where $D_t$ is foreign-currency debt, $\kappa$ is the collateral constraint coefficient, $\theta_t$ is an exogenous "credit conditions" variable, and $K_t$ is the capital stock. When $\theta_t$ falls (a sudden tightening of global credit conditions), the constraint binds: the economy must reduce debt to satisfy $D_{t+1} \leq \kappa\theta_{t+1}K_{t+1}$, requiring a current account surplus. Investment must contract sharply to satisfy the binding constraint.

The resulting dynamics are amplified by the **Fisherian debt-deflation** mechanism in a small open economy with foreign-currency debt: the depreciation triggered by the sudden stop raises the real debt burden, tightening the collateral constraint further, causing more capital outflows and more depreciation — a self-reinforcing spiral.

---

## 26.4 The Optimal Currency Area and the Eurozone

**Definition (Optimal Currency Area).** An **optimal currency area (OCA)** is a group of economies for which a single currency and common monetary policy yields higher welfare than maintaining separate currencies with flexible exchange rates between them. The concept was developed by Mundell (1961), who identified the conditions under which labor mobility substitutes for exchange rate flexibility as an adjustment mechanism after asymmetric shocks.

The OCA criteria: (i) high factor mobility (especially labor) within the area; (ii) high economic integration (trade and financial); (iii) synchronized business cycles (so that a common monetary policy is appropriate for all members simultaneously); (iv) fiscal transfer mechanisms to compensate regions hit by asymmetric shocks. When these conditions hold, the costs of giving up exchange rate flexibility are small (synchronized cycles mean the exchange rate would rarely need to change) and the benefits of a common currency (reduced transaction costs, price transparency, elimination of exchange rate risk) are large.

Applying OCA criteria to the eurozone reveals significant tensions. The eurozone was formed primarily for political reasons — to deepen European integration and make economic interdependence irreversible — rather than because it satisfied strict OCA criteria. Labor mobility within the eurozone is low (language and cultural barriers), business cycles are imperfectly synchronized across member states (Germany and periphery countries diverged significantly in the 2000s), and the eurozone lacks a fiscal union that could redistribute resources from booming to depressed regions.

The asymmetric shock problem: when a negative demand shock hits a peripheral eurozone country (say, Spain or Greece), the standard adjustment mechanisms are: (i) exchange rate depreciation → impossible within the eurozone; (ii) fiscal expansion → constrained by SGP rules and market access; (iii) wage and price cuts → extremely slow due to nominal rigidities; (iv) labor migration → limited by cultural and linguistic barriers; (v) fiscal transfers from the EU → limited by the small EU budget and political resistance. The 2010–13 eurozone debt crisis exposed these structural weaknesses sharply, with peripheral countries facing a brutal internal devaluation (wage and price cuts) to restore competitiveness — a process that took years and generated enormous social costs.

---

## 26.5 The International Monetary System

**Definition (International Monetary System).** The **international monetary system** is the set of institutional arrangements governing exchange rates, capital flows, and the provision of international reserve assets. It determines: what serves as the international unit of account and medium of exchange; how external imbalances are financed; and what adjustment mechanisms operate when those imbalances become unsustainable.

The history of the international monetary system is a sequence of regimes, each collapsing under the strain of its inherent tensions. The **gold standard** (roughly 1870–1914): exchange rates were fixed by each country's commitment to convert its currency to gold at a fixed price. External adjustment worked through the price-specie-flow mechanism: trade deficits caused gold outflows, reducing the money supply, deflating prices, improving competitiveness, and restoring balance. The system was symmetric and automatic but required procyclical adjustment (deflation during downturns) that imposed severe social costs. The **Bretton Woods system** (1944–71): fixed but adjustable exchange rates pegged to the dollar, which was in turn pegged to gold at $35/oz. Capital controls allowed countries to maintain exchange rate pegs and independent monetary policies simultaneously — resolving the trilemma by sacrificing capital mobility. The system collapsed when the U.S. gold commitment became inconsistent with expanding dollar balances abroad (Triffin's dilemma). The **current system** (post-1973): major currencies float against each other, with periodic interventions. Capital is largely free to move; the trilemma is resolved by allowing exchange rates to float. The dollar remains the dominant reserve currency despite the end of its gold link, an institutional arrangement known as the "Bretton Woods II" system (Dooley, Folkerts-Landau, and Garber, 2003).

---

*Next: Chapter 27 — Business Cycles: Booms, Recessions, and Stabilization*
