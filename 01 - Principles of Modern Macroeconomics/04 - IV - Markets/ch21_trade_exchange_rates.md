# Chapter 21 — International Trade and Exchange Rates: Connecting Economies

> *"When goods don't cross borders, soldiers will."*
> — Attributed to Frédéric Bastiat

---

The previous chapters analyzed goods, money, labor, and financial markets as if each national economy were a closed system. But economies are deeply interconnected: goods and services cross borders continuously, financial capital flows in search of the highest risk-adjusted returns, exchange rates translate the prices of one economy into the currency of another, and monetary policy decisions in the United States ripple through bond markets in Europe and equity markets in Asia within minutes. Understanding these international linkages is not optional for macroeconomics — it is central to understanding how monetary policy works, why current-account imbalances develop, and what determines the exchange rate.

This chapter develops the basic theory of open-economy macroeconomics in three parts. First, the foundations: comparative advantage and the gains from trade, the distinction between nominal and real exchange rates, and the purchasing power parity theory. Second, the key arbitrage conditions connecting domestic and foreign financial markets: covered and uncovered interest parity. Third, the most important dynamic model of exchange rate determination: the Dornbusch overshooting model, which shows why exchange rates are far more volatile than prices or interest rates, and why they overshoot their long-run equilibrium values following monetary shocks.

---

## 21.1 Comparative Advantage and the Gains from Trade

The theory of comparative advantage, developed by David Ricardo in 1817, is the most durable proposition in international economics: countries gain from specializing in and trading goods in which they have a comparative — not necessarily absolute — advantage in production.

**Definition (Comparative Advantage).** Country $H$ has a **comparative advantage** in good $x$ relative to good $y$ if its opportunity cost of producing $x$ in terms of $y$ is lower than in country $F$:

$$\frac{a_{Hx}}{a_{Hy}} < \frac{a_{Fx}}{a_{Fy}},$$

where $a_{ij}$ is the units of labor required to produce one unit of good $i$ in country $j$. Country $H$ should specialize in $x$ (and import $y$ from $F$), even if $H$ is less efficient in absolute terms at producing everything. The gains from trade arise from the expansion of the production possibility frontier when each country concentrates on its comparative advantage.

The Ricardian model has been extended to incorporate capital (Heckscher–Ohlin) and product differentiation (Krugman's new trade theory), but the core insight survives: international specialization generates mutual gains regardless of absolute productivity differences. The Stolper–Samuelson theorem, derived from the Heckscher–Ohlin model, establishes the distributional implications: trade liberalization in a capital-abundant country benefits capital owners and harms low-skilled workers (by reducing the relative demand for the scarce factor). This theorem underpins the policy concern that trade liberalization, while raising aggregate income, can generate large distributional costs concentrated among specific groups — an issue developed in Chapter 38.

---

## 21.2 The Nominal and Real Exchange Rate

**Definition (Nominal Exchange Rate).** The **nominal exchange rate** $e_t$ is the price of one unit of foreign currency in terms of domestic currency — the number of domestic currency units required to purchase one unit of foreign currency. Under this convention, a rise in $e_t$ is a **depreciation** of the domestic currency (it now takes more domestic currency to buy the same amount of foreign currency); a fall is an **appreciation**. (Note: some texts use the inverse convention; we follow the economics convention throughout.)

**Definition (Real Exchange Rate).** The **real exchange rate** $q_t$ is the price of the foreign consumption basket in terms of the domestic consumption basket:

$$q_t = \frac{e_t\, P_t^*}{P_t},$$

where $P_t^*$ is the foreign price level and $P_t$ is the domestic price level, both expressed in their respective currencies. When $q_t$ rises, foreign goods become more expensive relative to domestic goods — domestic goods have become more **competitive**. A depreciation of the nominal exchange rate ($e_t$ rises) with constant prices increases $q_t$, improving international competitiveness. An increase in domestic prices ($P_t$ rises) with constant $e_t$ and $P_t^*$ reduces $q_t$, eroding competitiveness.

The real exchange rate is the relative price that drives trade flows: exports increase when $q_t$ rises (domestic goods are cheaper for foreigners), and imports decrease (foreign goods are more expensive for domestic residents). In the IS framework, net exports $NX(q_t, Y_t, Y_t^*)$ are increasing in $q_t$ and decreasing in domestic income $Y_t$ (higher income raises imports), with the Marshall–Lerner condition ensuring that the net effect of depreciation is positive for the trade balance.

---

## 21.3 Purchasing Power Parity

**Definition (Purchasing Power Parity — Absolute).** **Absolute purchasing power parity (PPP)** states that identical goods sell for the same price across countries when expressed in the same currency:

$$P_t = e_t\, P_t^*, \quad \text{or equivalently} \quad e_t = \frac{P_t}{P_t^*}.$$

Absolute PPP implies $q_t = 1$: the real exchange rate equals one. It is the law of one price applied to the aggregate price level.

**Definition (Purchasing Power Parity — Relative).** **Relative PPP** states that the rate of change of the exchange rate equals the difference in inflation rates:

$$\hat{e}_t \equiv \frac{\Delta e_t}{e_{t-1}} = \pi_t - \pi_t^*.$$

Relative PPP does not require price levels to be equalized across countries — only that their *rates of change* move together.

### Empirical Evidence and the Balassa–Samuelson Effect

Absolute PPP is clearly violated in the data: price levels differ enormously across countries, with rich countries having systematically higher price levels than poor countries. The most prominent explanation is the **Balassa–Samuelson effect**: productivity differences between countries are concentrated in the traded goods sector (manufacturing, technology) while non-traded services have more uniform productivity. In a country with high productivity in traded goods:

- Wages in traded goods are bid up by high productivity.
- Wages in non-traded services must rise to retain workers (labor mobility within a country).
- Higher non-traded sector wages raise the prices of non-traded services.
- The overall price level ($P$) is higher, generating a higher real exchange rate $q = eP^*/P < 1$.

Formally, let $\mathcal{A}_T$ and $\mathcal{A}_{NT}$ be traded and non-traded productivity. The Balassa–Samuelson model predicts:

$$\hat{e}_t = (\pi_t - \pi_t^*) + \left(\hat{\mathcal{A}}_{T,t}^* - \hat{\mathcal{A}}_{T,t}\right) + \left(\hat{\mathcal{A}}_{NT,t} - \hat{\mathcal{A}}_{NT,t}^*\right).$$

Countries with faster traded-sector productivity growth (like Korea and Japan in their catch-up phases) tend to experience real exchange rate appreciation over time. This explains why rapid-growth emerging markets often appear to have undervalued currencies under simple absolute PPP comparisons: their low price levels partly reflect genuinely lower non-traded goods prices, not currency undervaluation.

Relative PPP holds much better at long horizons and high inflation rates. The half-life of real exchange rate deviations from PPP is estimated at approximately 3–5 years (Rogoff, 1996), meaning departures from PPP are substantial but slowly reverting.

---

## 21.4 Interest Parity Conditions

Financial capital can flow freely across borders in a world without capital controls. This mobility creates powerful arbitrage conditions connecting domestic and foreign interest rates and the exchange rate.

### Covered Interest Parity

An investor can eliminate exchange rate risk by covering the foreign investment with a forward contract. Invest domestically at rate $i_t$ or invest abroad at rate $i_t^*$ while selling the proceeds forward at rate $f_t$ (the forward exchange rate):

- Domestic investment: \$1 grows to $(1+i_t)$ at $t+1$.
- Covered foreign investment: \$1 buys $1/e_t$ foreign currency, which grows to $(1+i_t^*)/e_t$ abroad, which is sold forward for $(1+i_t^*)\cdot f_t/e_t$ domestic currency at $t+1$.

No-arbitrage requires these two returns to be equal:

$$(1+i_t) = (1+i_t^*)\frac{f_t}{e_t},$$

or in log approximation:

$$i_t - i_t^* = f_t - e_t, \quad \text{where } f_t = \ln F_t,\; e_t = \ln E_t.$$

**Definition (Covered Interest Parity).** **Covered interest parity (CIP)** states that the domestic-foreign interest rate differential equals the forward premium on the domestic currency:

$$i_t - i_t^* = f_t - e_t.$$

CIP is a pure no-arbitrage condition that should hold exactly in frictionless markets with equal counterparty risk. In practice, it held remarkably well until 2008 but broke down systematically during the Global Financial Crisis and again after 2015, when CIP deviations became persistent (Du, Tepper, and Verdelhan, 2018). The deviations reflect constraints on bank balance sheets that prevent fully capitalizing on the arbitrage — financial intermediary balance sheets are the binding constraint on capital markets.

### Uncovered Interest Parity

If investors are risk-neutral and do not hedge their foreign investments, they require equal expected returns across currencies. An investor who buys a foreign bond rather than a domestic bond expects to receive:

$$(1+i_t^*)\frac{\mathbb{E}_t[e_{t+1}]}{e_t} = (1+i_t).$$

In log approximation, the **uncovered interest parity (UIP)** condition is:

$$i_t - i_t^* = \mathbb{E}_t[\hat{e}_{t+1}] \equiv \mathbb{E}_t[e_{t+1} - e_t].$$

**Definition (Uncovered Interest Parity).** **Uncovered interest parity (UIP)** states that the interest rate differential between two countries equals the expected change in the exchange rate:

$$i_t = i_t^* + \mathbb{E}_t[\hat{e}_{t+1}].$$

A country with a higher interest rate than its trading partner must expect its currency to depreciate — otherwise risk-neutral investors would all want to invest in the high-rate country and the rate differential could not be sustained.

UIP is weaker than CIP because it involves expectations rather than observable forward rates, and it assumes risk neutrality. The **UIP puzzle** (Fama, 1984): empirically, high-interest-rate currencies tend to *appreciate* rather than depreciate, exactly the opposite of UIP's prediction. The "carry trade" — borrowing in low-rate currencies and investing in high-rate currencies — generates positive expected profits, which is inconsistent with risk-neutral UIP. The resolution likely involves risk premia: the carry trade delivers low or negative returns in bad states (financial crises, risk-off episodes), so its positive average return is compensation for bearing this risk.

---

## 21.5 The Dornbusch Overshooting Model

The most important theoretical model of exchange rate dynamics is the **Dornbusch (1976) overshooting model**, which explains why exchange rates are far more volatile than fundamentals (money supplies, price levels, interest rates) and why they respond non-monotonically to policy shocks — initially moving too far and then gradually converging to the new equilibrium.

The model has two sectors with different speeds of adjustment. Asset markets (money and bonds) clear instantaneously; goods markets (prices) adjust slowly. This difference in speeds generates the overshooting.

### Setup

The model has five equations:

**Goods-market equilibrium (aggregate demand):**

$$y_t = \bar{y} + \delta(e_t - p_t) - \sigma r_t, \quad \delta, \sigma > 0,$$

where $e_t - p_t$ is the log real exchange rate (competitiveness) and $r_t$ the real interest rate.

**Money-market equilibrium (LM in log form):**

$$m_t - p_t = \phi y_t - \lambda i_t, \quad \phi, \lambda > 0.$$

**Uncovered interest parity (asset market clearing):**

$$i_t = i^* + \mathbb{E}_t[\dot{e}_t].$$

**Price adjustment (sticky prices in the goods market):**

$$\dot{p}_t = \pi(y_t - \bar{y}), \quad \pi > 0.$$

**Long-run equilibrium:** In the long run, $\dot{p} = 0$, which requires $y = \bar{y}$, and from the LM equation the long-run nominal exchange rate $\bar{e}$ satisfies $\bar{e} = \bar{p} + (m - \phi\bar{y} + \lambda i^*)/1$.

### Saddle-Path Dynamics and Overshooting

The model's key insight comes from combining UIP with the slow price adjustment equation. In the saddle-path equilibrium, after a monetary expansion (increase in $m$), the exchange rate must jump immediately to preserve UIP: the long-run equilibrium $\bar{e}$ moves proportionally with $m$ (one-for-one by the quantity theory), but the goods market cannot adjust instantly.

At the moment of the monetary shock:
- $p_t$ is predetermined (cannot jump — sticky prices).
- $e_t$ is free to jump immediately (asset markets clear instantaneously).
- The long-run equilibrium shifts to $\bar{e}' > \bar{e}$.

For UIP to hold at the instant of the shock with $i^*$ fixed, the domestic interest rate must fall (monetary expansion). For the interest rate to fall below $i^*$, investors must expect the exchange rate to *appreciate* (domestic currency to strengthen) from the new post-shock level. This means $e_0 > \bar{e}'$ — the exchange rate must *overshoot* the new long-run equilibrium on impact.

The exchange rate jumps to $e_0 = \bar{e}' + (\bar{e}' - \bar{e})\cdot(1/(\pi\lambda) + 1) > \bar{e}'$ on impact, and then gradually appreciates back toward $\bar{e}'$ as the price level rises. The adjustment satisfies the saddle-path condition:

$$e_t - \bar{e}' = (e_0 - \bar{e}')e^{-\mu t}, \quad \mu > 0,$$

where $\mu$ is the stable eigenvalue of the linearized system. The overshooting result:

$$e_0 - \bar{e}' = -\frac{1}{\lambda\mu}(m_t - m_{t-1}),$$

which shows that the extent of overshooting is larger when money demand is less interest-elastic ($\lambda$ small) and when goods prices adjust slowly ($\mu$ small — which happens when $\pi$ is small).

The Dornbusch model provides an elegant explanation for the **exchange rate disconnect puzzle** — the empirical observation that exchange rates are far more volatile than the fundamentals (inflation differentials, interest rate differentials, trade balances) that are supposed to drive them. If goods markets adjust slowly but asset markets adjust instantly, exchange rates will always overshoot their long-run values, generating large short-run swings around slowly-moving fundamentals.

---

## 21.6 The Mundell–Fleming Trilemma and Policy Implications

The **Mundell–Fleming trilemma** (or "impossible trinity") states that a country cannot simultaneously maintain all three of:

1. **A fixed exchange rate**: committing to peg the currency to another at a fixed rate.
2. **Perfect capital mobility**: allowing free cross-border flows of financial capital.
3. **Independent monetary policy**: setting the domestic interest rate freely.

**Definition (Mundell–Fleming Trilemma).** The trilemma arises because: if the exchange rate is fixed and capital flows freely, the domestic interest rate must equal the world rate $i^* $ (by UIP). If the central bank tries to set $i \neq i^*$, capital flows immediately to exploit the differential, forcing the central bank to buy or sell foreign exchange reserves until either the exchange rate peg breaks or the interest rate returns to $i^*$. Countries can achieve any two of the three objectives but not all three simultaneously.

Historical exchange rate regimes map onto choices along the trilemma triangle. The Bretton Woods system (1944–71) maintained fixed exchange rates and capital controls (sacrificing capital mobility) to preserve monetary independence. After the breakdown of Bretton Woods, advanced economies moved to floating rates with capital mobility, regaining monetary independence. The eurozone chose fixed exchange rates (a single currency) with capital mobility, sacrificing independent monetary policy for member states.

The trilemma has important implications for the effectiveness of macroeconomic policies. Under a fixed exchange rate with perfect capital mobility, monetary policy is powerless (any rate cut triggers capital outflows that reverse it), while fiscal policy is fully potent (no crowding out, because the interest rate cannot rise). Under a floating rate with capital mobility, monetary policy works primarily through the exchange rate channel (a rate cut depreciates the currency, stimulating net exports), while fiscal policy is partly crowded out by exchange rate appreciation.

---

*Next: Chapter 22 — The Government Sector: Fiscal Policy and Public Debt*
