# Chapter 17 — The Goods Market: Equilibrium and Adjustments

> *"Supply and demand — the centrepiece of all economic analysis — is really nothing other than a theory of how prices coordinate the decisions of self-interested individuals."*
> — Paul Seabright, *The Company of Strangers*, 2004

---

The models of Part II introduced goods-market equilibrium in the IS–LM framework, and the microfoundations of Part III derived the behavioral equations underlying aggregate consumption and investment. This chapter revisits goods-market equilibrium from a more integrated perspective, asking: what determines whether output is demand-determined or supply-determined, how do adjustment dynamics operate in real time, and what is the precise structure of the forward-looking New Keynesian IS curve that underlies the three-equation policy framework? The chapter also addresses a subtlety that the simpler treatments of Chapters 7–9 glossed over: the distinction between goods-market equilibrium as an instantaneous condition and the dynamic process through which the economy converges to it.

---

## 17.1 Goods-Market Equilibrium: The General Framework

In a closed economy, the goods market clears when the value of output produced equals the value of output demanded. This is the national income identity $Y \equiv C + I + G$ reinterpreted as an equilibrium condition: when firms produce exactly the quantity that all sectors wish to purchase at current prices, there is no pressure for output to change.

The goods-market equilibrium condition in its most general form is:

$$Y_t = \mathcal{E}(Y_t,\, r_t,\, G_t,\, T_t,\, \mathbb{E}_t[\mathbf{x}_{t+1}]),$$

where $\mathcal{E}$ is aggregate planned expenditure, $r_t$ is the real interest rate, $G_t$ government spending, $T_t$ taxes, and $\mathbb{E}_t[\mathbf{x}_{t+1}]$ is the vector of expectations of all future relevant variables. The partial derivative $\partial\mathcal{E}/\partial Y_t = C'(1-t) + m_Y^I \in (0,1)$ — the marginal propensity to spend out of income — is less than one (the stability condition). Partial derivatives with respect to $r_t$ are negative (higher rates reduce investment and durable consumption); with respect to $G_t$ positive; with respect to $T_t$ negative; with respect to expected future income positive (forward-looking households smooth consumption).

**Definition (Goods-Market Equilibrium).** The goods market is in **equilibrium** at date $t$ when planned expenditure equals actual output:

$$Y_t = \mathcal{E}(Y_t, r_t, G_t, T_t, \mathbb{E}_t[\mathbf{x}_{t+1}]).$$

This is a fixed-point condition: equilibrium output $Y_t^*$ is the level at which output and expenditure are mutually consistent. It exists and is unique under the stability condition $\partial\mathcal{E}/\partial Y < 1$, since then the function $Y \mapsto \mathcal{E}(Y, \cdot)$ is a contraction and has a unique fixed point.

The equilibrium is **not** a full-economy equilibrium in the sense of Chapter 1 — it holds conditional on a given interest rate $r_t$ and given expectations. Full general equilibrium requires all markets to clear simultaneously, including the money market (which pins down $r_t$) and the labor market (which ensures the real wage is consistent with employment). The IS–LM model addresses the first issue; Part V addresses both together.

---

## 17.2 Inventory Adjustment and the Convergence to Equilibrium

When expenditure exceeds output, firms find their inventories depleting unexpectedly; when output exceeds expenditure, inventories accumulate. These **unplanned inventory changes** provide the disequilibrium signal that drives the adjustment process.

**Definition (Unplanned Inventory Investment).** At any date, actual gross investment $I_t$ equals planned investment $I_t^p$ plus unplanned inventory accumulation $\Delta N_t^u$:

$$I_t = I_t^p + \Delta N_t^u, \quad \Delta N_t^u = Y_t - \mathcal{E}(Y_t, \cdot).$$

When $Y_t > \mathcal{E}$, inventories accumulate unintentionally ($\Delta N_t^u > 0$); when $Y_t < \mathcal{E}$, inventories are depleted ($\Delta N_t^u < 0$). Firms respond by adjusting production in the following period: rising inventories signal overproduction and prompt output cuts; falling inventories signal excess demand and prompt output expansion.

The formal adjustment process, in continuous time, is:

$$\dot{Y}_t = \phi\bigl[\mathcal{E}(Y_t, r_t, \cdot) - Y_t\bigr], \quad \phi > 0,$$

where $\phi$ is the speed of adjustment. This is a first-order linear ODE. Linearizing around equilibrium $Y^*$:

$$\dot{Y}_t \approx \phi\bigl[\partial\mathcal{E}/\partial Y - 1\bigr](Y_t - Y^*) = -\phi(1 - \mathcal{E}_Y)(Y_t - Y^*),$$

where $\mathcal{E}_Y \equiv \partial\mathcal{E}/\partial Y \in (0,1)$. The coefficient $-\phi(1-\mathcal{E}_Y) < 0$ is negative, so the system is globally stable: any deviation from $Y^*$ is corrected at rate $\phi(1-\mathcal{E}_Y)$. The equilibrium $Y^*$ is approached asymptotically from any initial condition:

$$Y_t - Y^* = (Y_0 - Y^*)e^{-\phi(1-\mathcal{E}_Y)t}.$$

The speed of convergence $\phi(1-\mathcal{E}_Y)$ decreases with $\mathcal{E}_Y$: the closer the marginal propensity to spend is to one, the slower the adjustment. At $\mathcal{E}_Y = 1$ (zero multiplier stability margin), adjustment is arbitrarily slow and the model breaks down.

This dynamic specification captures an important empirical feature of goods-market adjustment: it is not instantaneous. Following a demand shock, output adjusts over several quarters as firms revise production plans, hire workers, and rebuild or run down inventory buffers. The slow adjustment is part of why business cycles have duration measured in years rather than months.

---

## 17.3 The New Keynesian IS Curve: Deriving the Forward-Looking Relationship

The IS curves of Chapters 8 and 9 were static: they related current output to the current interest rate without explicitly modeling expectations. The New Keynesian IS curve is derived from household intertemporal optimization (the Euler equation of Chapter 11) and is inherently dynamic and forward-looking. Its derivation is one of the most important in modern macroeconomics.

Start from the household's consumption Euler equation, which we derived in Chapter 11 for the CRRA case:

$$\frac{\dot{c}_t}{c_t} = \frac{r_t - \rho}{\sigma},$$

or in discrete time, using the log approximation:

$$\ln c_{t+1} - \ln c_t \approx \frac{1}{\sigma}(r_t - \rho) + \text{uncertainty term}.$$

In a closed economy with no government and no investment, goods-market clearing requires $c_t = Y_t$ at every date. Therefore $\ln Y_{t+1} - \ln Y_t = (r_t - \rho)/\sigma$, or in terms of log-deviations from the balanced growth path (denoted $\hat{y}_t$):

$$\hat{y}_t = \mathbb{E}_t[\hat{y}_{t+1}] - \frac{1}{\sigma}(r_t - \rho).$$

With investment and government spending, the derivation requires a more careful aggregation of the household and firm Euler equations, but the structure is the same. Let $\hat{x}_t \equiv \hat{y}_t - \hat{\bar{y}}_t$ denote the output gap — the deviation of output from its flexible-price (potential) level — and $r_t^n$ the natural (Wicksellian) rate of interest, defined as the real interest rate consistent with $\hat{x}_t = 0$. The **New Keynesian IS curve** (also called the Dynamic IS equation or DIS) is:

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\pi_{t+1}] - r_t^n),$$

where $i_t$ is the nominal interest rate and $\mathbb{E}_t[\pi_{t+1}]$ is expected inflation, so that $i_t - \mathbb{E}_t[\pi_{t+1}]$ is the ex-ante real interest rate.

**Definition (Natural Rate of Interest — Revisited).** The **natural rate of interest** $r_t^n$ is the real interest rate that would prevail if prices and wages were fully flexible — the equilibrium real rate in the RCK model. It is determined by productivity growth, households' discount rate, and other supply-side factors. It is not directly observable but can be estimated from the gap between the observed real interest rate and the output gap. When $i_t - \mathbb{E}_t[\pi_{t+1}] = r_t^n$, the output gap is constant and evolves only with expectation changes.

The DIS has several critical properties that distinguish it from the static IS curve.

**Forward-looking dynamics.** The output gap today depends on the *expected* future output gap, not just the current interest rate. This is the hallmark of the New Keynesian approach: current outcomes are functions of the entire expected future path of policy. A promise of low interest rates in the future — forward guidance — affects current output even before rates actually change.

**Iterating the DIS forward.** By substituting the DIS equation one period ahead into itself and iterating:

$$\hat{x}_t = -\sigma\sum_{k=0}^{\infty}\mathbb{E}_t\bigl[i_{t+k} - \pi_{t+k+1} - r_{t+k}^n\bigr].$$

Current output is the infinite discounted sum of all future expected deviations of the real interest rate from the natural rate. This is the IS-curve analogue of the long-rate pricing formula for bonds: what matters is not just today's short rate but the entire expected path of real rates. Monetary policy therefore operates through the *entire yield curve*, not just the overnight rate, which is why central banks devote so much attention to expectations management.

**The role of $\sigma$.** The parameter $\sigma$ (the inverse of the EIS) governs the sensitivity of the output gap to interest rate deviations. With $\sigma = 1$ (log utility), a 100-basis-point deviation of the real rate from neutral for one year generates a 1-percentage-point output gap. With $\sigma = 2$, the same gap requires only a 50-basis-point interest rate deviation. Empirical estimates of $\sigma$ vary widely across specifications and data sets, and this uncertainty is a major source of uncertainty about the transmission power of monetary policy.

---

## 17.4 Demand Complementarities and the Possibility of Multiple Equilibria

The standard goods-market equilibrium has a unique fixed point under the stability condition $\mathcal{E}_Y < 1$. But several mechanisms can generate multiple equilibria — situations in which the economy may settle at different output levels depending on expectations and coordination, with no fundamental difference in technology or preferences.

**Definition (Strategic Complementarity in Demand).** Demand decisions exhibit **strategic complementarity** when the optimal choice of each agent is increasing in the aggregate choices of other agents — when higher aggregate demand makes it optimal for each individual agent to demand more. Strategic complementarities in demand create a self-reinforcing structure: optimism generates high demand, which validates the optimism.

To formalize, suppose aggregate expenditure takes the form:

$$\mathcal{E}(Y_t;\, Y_t^e) = \mathcal{E}_1(Y_t) + \mathcal{E}_2(Y_t^e),$$

where $Y_t^e$ is the expected level of aggregate demand and both $\mathcal{E}_1$ and $\mathcal{E}_2$ are increasing. The fixed-point condition for a rational expectations equilibrium requires $Y_t = Y_t^e$ (expectations are fulfilled), giving equilibrium conditions:

$$Y^* = \mathcal{E}_1(Y^*) + \mathcal{E}_2(Y^*).$$

When $\partial(\mathcal{E}_1 + \mathcal{E}_2)/\partial Y^* < 1$, the fixed point is unique. But if the combined derivative exceeds one over some range, multiple equilibria emerge — a "good" equilibrium $Y_H$ and a "bad" equilibrium $Y_L < Y_H$ — separated by an unstable interior fixed point $Y_M$.

Mechanisms that generate strategic complementarities in demand include: **increasing returns in production** (as aggregate output rises, input costs fall and productivity rises, making additional production more profitable for each firm); **thick-market externalities** (more trade means better matches, higher productivity, and more purchasing power for all participants); and **aggregate demand externalities in New Keynesian models** (when nominal rigidities prevent firms from cutting prices during downturns, the demand externality works through the real wage and employment channels, as analyzed by Blanchard and Kiyotaki, 1987).

The relevance of multiple equilibria for macroeconomic policy is significant. If the economy can coordinate on different equilibria depending on expectations, then "animal spirits" — coordinating devices like newspaper headlines, presidential speeches, or market sentiment indices — can shift the economy between equilibria without any change in fundamentals. This provides a theory of self-fulfilling recessions distinct from models that require fundamental shocks.

---

## 17.5 The Goods Market in the Open Economy

In an open economy, the goods-market equilibrium condition becomes:

$$Y_t = C_t + I_t + G_t + NX_t,$$

where net exports $NX_t = X_t(Y_t^*, q_t) - M_t(Y_t, q_t)$ depend on foreign income $Y_t^*$, the real exchange rate $q_t = e_t P_t^*/P_t$, and domestic income $Y_t$. The real exchange rate enters because depreciation ($q$ rises) makes domestic goods cheaper relative to foreign goods, increasing exports and reducing imports.

Differentiating the open-economy equilibrium condition with respect to the real exchange rate:

$$\frac{\partial Y^*}{\partial q} = \frac{X_q - M_q}{1 - C'(1-t) + M_Y} > 0 \quad \text{if } X_q > M_q,$$

where $X_q > 0$ (depreciation raises exports) and $M_q < 0$ (depreciation reduces imports, since foreign goods become more expensive). The condition $X_q - M_q > 0$ — that depreciation improves the trade balance — is the **Marshall–Lerner condition**: in levels, it requires that the sum of the price elasticities of export and import demand exceeds one.

**Definition (Marshall–Lerner Condition).** The **Marshall–Lerner condition** states that a real depreciation improves the trade balance if and only if:

$$|\varepsilon_X| + |\varepsilon_M| > 1,$$

where $\varepsilon_X > 0$ is the price elasticity of export demand and $\varepsilon_M > 0$ is the price elasticity of import demand (both expressed as positive numbers). When this condition holds, a 1% real depreciation increases the volume of exports and decreases the volume of imports by enough in total to raise the real value of net exports.

Empirically, the Marshall–Lerner condition is satisfied in the medium and long run for most countries, but not necessarily in the short run. In the short run, after a depreciation, the volumes of imports and exports adjust slowly while the prices change immediately — because trade contracts are typically denominated in the exporter's currency and take time to renegotiate. This generates the **J-curve**: the trade balance initially worsens after a depreciation (because import prices rise immediately but volumes fall only gradually) before eventually improving as volumes adjust.

---

## 17.6 Goods-Market Equilibrium and the Output Gap in DSGE Models

In a quantitative DSGE model used for policy analysis, the goods-market equilibrium condition is not a single equation but a system of resource constraints and market-clearing conditions that must hold simultaneously at every date. The representative form of the goods-market block in a medium-scale NK DSGE model (following Smets and Wouters, 2007) is:

$$Y_t = c_y C_t + i_y I_t + g_y G_t + \epsilon_t^m,$$

where $c_y$, $i_y$, $g_y$ are steady-state shares of each expenditure component in output (summing to one in a closed economy), and $\epsilon_t^m$ is a measurement error or non-modeled expenditure component. Each component $C_t$, $I_t$ is determined by its respective optimization problem (Chapters 11 and 12), connected through factor markets and the interest rate to the monetary policy block.

The estimated output gap from such a model differs from the HP-filtered output gap because the DSGE model-based gap accounts for supply-side shocks (technology, preference, markup shocks) that move potential output, while the HP filter removes only the statistical trend. During the COVID pandemic, for instance, lockdowns shifted the supply side dramatically, causing a large divergence between HP-filtered and model-based output gaps — with the model-based gap being much smaller (because potential output fell alongside actual output) and the HP-filtered gap large (because the statistical trend was unaffected by the supply shock).

---

*Next: Chapter 18 — The Money Market*
