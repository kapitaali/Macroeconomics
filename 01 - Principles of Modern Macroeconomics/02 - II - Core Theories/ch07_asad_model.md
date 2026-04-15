# Chapter 7 — The Aggregate Demand–Aggregate Supply Model: The Workhorse of Macroeconomics

> *"In the long run, prices and wages adjust to clear markets. In the short run, they do not."*
> — N. Gregory Mankiw

---

The IS–LM model (developed in Chapter 9) characterizes short-run equilibrium in the goods and money markets for a given price level. But the price level is itself an endogenous variable: it responds to the balance between aggregate demand and the economy's productive capacity. The **aggregate demand–aggregate supply (AS–AD) framework** closes this gap by combining a theory of how demand determines output (for a given price level) with a theory of how supply — the behavior of price and wage setters — determines the price level (for a given level of demand). Together, the two curves determine both real output and the price level simultaneously, and the intersection of them defines the short-run macroeconomic equilibrium from which the economy adjusts toward long-run potential.

The AS–AD framework is a workhorse rather than a frontier model. Its simplicity is its strength: it organizes thinking about short-run stabilization, supply shocks, and the distinction between demand-pull and cost-push inflation in a way that is directly applicable to policy analysis. Its limitations — the absence of intertemporal optimization, rational expectations, and a financial sector — motivate the more sophisticated models of subsequent parts.

---

## 7.1 Aggregate Demand

**Definition (Aggregate Demand Curve).** The **aggregate demand (AD) curve** is the locus of price level $P$ and real output $Y$ combinations at which the goods market (IS equilibrium) and the money market (LM equilibrium) are simultaneously satisfied, holding all other exogenous variables fixed.

### Derivation and Slope

To derive the slope, recall the IS condition ($Y = C(Y-T) + I(r) + G$) and the LM condition ($M/P = L(Y, r + \pi^e)$). When $P$ rises, the real money supply $M/P$ falls, shifting LM upward and raising the interest rate for any given $Y$, which reduces investment and hence output. Totally differentiating the LM equation:

$$\mathrm{d}r = \frac{1}{L_r}\!\left(-\frac{M}{P^2}\mathrm{d}P - L_Y\,\mathrm{d}Y\right).$$

Substituting into the differentiated IS curve and solving:

$$\left.\frac{\mathrm{d}Y}{\mathrm{d}P}\right|_{AD} = \frac{I'(r)\cdot M/P^2}{L_Y I'(r) - L_r} < 0.$$

The AD curve is **downward sloping** through three distinct channels.

The **real-balance effect** (Pigou effect): a higher price level reduces real wealth $M/P$ for holders of nominal money balances, reducing their consumption. This effect is typically small because money holdings are a small fraction of total household wealth.

The **interest rate effect** (Keynes effect): the dominant channel. A higher price level reduces the real money supply, pushing up interest rates, which reduce interest-sensitive investment spending and consumer durables purchases.

The **exchange rate effect** (Mundell-Fleming effect), operative in open economies: a higher domestic price level makes domestic goods more expensive relative to foreign goods, reducing net exports. This channel is developed fully in Chapter 21.

### Shifts in the AD Curve

Any change in an exogenous variable other than $P$ that affects IS or LM shifts the entire AD curve. The most important shifters are fiscal policy (an increase in $G$ or a cut in $T$ shifts AD right by the IS-LM fiscal multiplier times the change), monetary policy (an increase in $M$ shifts LM right and AD right), consumer or investor confidence (autonomous shifts in $C$ or $I$ shift AD), and, in the open economy, foreign income or the exchange rate. Understanding precisely which shocks shift AD versus SRAS — versus both simultaneously — is essential for diagnosing macroeconomic episodes.

---

## 7.2 Aggregate Supply

### Long-Run Aggregate Supply

In the long run, all nominal prices and wages are fully flexible. Firms produce at the level determined by the real productive capacity of the economy regardless of the price level.

**Definition (Long-Run Aggregate Supply).** The **long-run aggregate supply (LRAS) curve** is a vertical line at $Y = \bar{Y}$, where $\bar{Y}$ is potential output — the level of output consistent with the economy's capital, labor, and technology endowment, the natural rate of unemployment, and the efficient allocation of resources. It is vertical because a higher price level in the long run leads to proportionally higher nominal wages, leaving the real wage and hence employment and output unchanged. The LRAS shifts only when potential output itself changes: capital accumulation, labor force growth, improvements in education and skills, or technological change.

The distinction between SRAS and LRAS corresponds to the short-run/long-run distinction in macroeconomics. The long run is not a fixed calendar length — it is the horizon over which all prices and wages have adjusted to their market-clearing levels. In economies with flexible labor markets and limited nominal contracting, the long run may be measured in months; in economies with multi-year wage contracts and substantial nominal rigidity, it may stretch to several years.

### Short-Run Aggregate Supply: The Sticky-Wage Model

In the short run, nominal wages are typically fixed by contracts or by the social norms that make rapid nominal wage cuts psychologically and institutionally difficult — even when prices have fallen and the real wage has risen above its market-clearing level. Suppose the nominal wage is fixed at $\bar{W}$. A profit-maximizing firm hires labor until the marginal product of labor equals the real wage:

$$F_L(K, N) = \frac{\bar{W}}{P}.$$

When $P$ rises, the real wage $\bar{W}/P$ falls, making it profitable to hire more workers and produce more output. Hence output rises with $P$ in the short run. Solving for $N$ as a function of $P$ and substituting into the production function gives:

$$Y^s = F\!\left(K,\, N\!\left(\frac{\bar{W}}{P}\right)\right), \quad \frac{\partial Y^s}{\partial P} = -F_L\cdot N'\!\left(\frac{\bar{W}}{P}\right)\cdot\frac{\bar{W}}{P^2} > 0.$$

The SRAS curve is upward sloping. Its position — the price level at which output equals potential — is pinned down by the nominal wage $\bar{W}$. When nominal wages rise (a shift in the SRAS curve input), the SRAS shifts left; when wages fall, it shifts right. Supply shocks that raise input costs — an oil price increase, a drought, a supply chain disruption — also shift SRAS left by raising costs at every level of production.

### Short-Run Aggregate Supply: The Lucas Imperfect-Information Model

A different foundation for the upward-sloping SRAS, based on rational behavior with imperfect information rather than nominal rigidities, comes from Lucas (1973). Producers observe their own output prices precisely but observe the general price level only with noise. When a producer sees the price of her product rise, she cannot immediately determine whether relative demand for her product has increased (a real signal warranting more production) or whether all prices have risen uniformly (a nominal signal warranting no change). The Bayesian signal extraction problem yields:

$$Y_t = \bar{Y}_t + \alpha(P_t - P_t^e), \quad \alpha > 0,$$

where $P_t^e \equiv \mathbb{E}_{t-1}[P_t]$ is the expected price level. Output exceeds potential only when the actual price level surprises agents — only when there is genuine unanticipated inflation. Under rational expectations, this means that only unanticipated monetary shocks can raise output above potential, while anticipated policy changes are fully reflected in expectations and have no real effect. This is the foundation of the policy ineffectiveness proposition [Ch. 16].

The key empirical prediction of the Lucas model: the slope $\alpha$ should be lower in high-inflation countries, where agents are less able to distinguish relative from aggregate price changes. Lucas (1973) confirmed this cross-country prediction: the output-inflation trade-off is indeed shallower in high-inflation economies — a finding consistent with rational signal extraction and inconsistent with models in which the Phillips curve slope is a structural constant.

---

## 7.3 Short-Run and Long-Run Equilibrium

### Solving the System

Short-run equilibrium is the intersection of AD and SRAS. Using the Lucas supply function and a log-linear AD relationship:

$$\text{SRAS: } Y_t = \bar{Y}_t + \alpha(P_t - P_t^e)$$
$$\text{AD: } Y_t = Y_0 - \gamma P_t + \text{(exogenous shifts)}$$

where $\gamma > 0$ captures the combined AD slope from interest rate and exchange rate effects. These two equations simultaneously determine $Y_t$ and $P_t$ given the state of expectations and exogenous policy variables.

Long-run equilibrium adds the fulfilled-expectations condition $P_t = P_t^e$, which forces $Y_t = \bar{Y}_t$. The price level is then determined by the intersection of AD with the vertical LRAS.

### The Adjustment Mechanism

The adjustment mechanism from short-run to long-run equilibrium operates through the gradual revision of nominal wages and prices. If the economy is below potential ($Y_t < \bar{Y}_t$), excess supply of labor puts downward pressure on nominal wages; as wages fall, the SRAS curve shifts right over time, increasing supply and gradually restoring output to $\bar{Y}_t$ at a lower price level. If the economy is above potential, excess demand in labor markets drives wages up, shifting SRAS left and reducing output back to potential.

The speed of adjustment depends on how quickly nominal wages and prices respond to demand conditions — equivalently, on the degree of nominal rigidity. An important asymmetry: downward nominal wage rigidity (the empirical tendency for nominal wages to resist cuts) makes adjustment from above potential faster (wages rise readily when demand is strong) than adjustment from below potential (wages resist falling when demand is weak), generating a systematic difference in the duration and depth of inflationary versus recessionary episodes.

---

## 7.4 The Macroeconomic Effects of Shocks

One of the central uses of the AS-AD framework is organizing the economy's response to exogenous disturbances. Two shock classes with different signatures in the data deserve careful treatment.

### Demand Shocks

A **positive demand shock** — an increase in government spending, a monetary expansion, or a surge in consumer confidence — shifts the AD curve rightward. In the short run, before wages and prices have fully adjusted, output rises above potential and the price level rises: both output and prices move together. This co-movement of output and prices is the signature of a demand shock and distinguishes it from a supply shock empirically.

In the long run, wages adjust upward, shifting the SRAS left until output returns to $\bar{Y}$ at a permanently higher price level. The real effects of the demand shock are entirely temporary; the nominal effects are permanent. This is the sense in which money is neutral in the long run: a sustained monetary expansion ultimately raises the price level proportionally without affecting output.

### Supply Shocks

An **adverse supply shock** — an oil price spike, a natural disaster, a pandemic-induced lockdown, a crop failure — shifts the SRAS leftward: at any given price level, firms produce less because their costs have risen. Short-run equilibrium exhibits **stagflation**: lower output and a higher price level simultaneously. This distinctive combination — inflation and recession together — was precisely what the 1970s oil shocks produced. It was inexplicable under purely demand-oriented models that dominated policy at the time, which could not generate a simultaneous worsening of inflation and unemployment. The AS-AD framework, by distinguishing demand and supply disturbances and their different inflation-output signatures, was the theoretical advance that made stagflation comprehensible.

The appropriate policy response to a supply shock is also different from the response to a demand shock. Accommodating a supply shock with monetary expansion prevents the output decline but accepts the full inflationary impact; tightening policy to prevent inflation accepts the full output decline; intermediate responses involve a policy trade-off on the AS-AD diagram. This trade-off is formalized as the central bank's stabilization problem in the New Keynesian model [Ch. 23].

---

## 7.5 Empirical Application: The 1970s Stagflation and the 2021–22 Inflation Surge

### The 1970s

The oil price shock of 1973–74 (OPEC embargo) quadrupled the price of oil in less than a year. The AS-AD framework's prediction: adverse supply shock shifts SRAS left → stagflation (higher inflation, lower output). The data confirmed this precisely. U.S. CPI inflation rose from 3.4% in 1972 to 12.3% in 1974; real GDP fell by 0.5% in 1974 and 0.2% in 1975; unemployment rose from 4.9% in 1973 to 8.5% in 1975. The Federal Reserve faced the AS-AD trade-off directly: tightening to contain inflation would deepen the recession; accommodating to preserve output would allow inflation to become entrenched. The policy failure — accommodation — allowed inflation expectations to become dislodged, ultimately requiring the Volcker disinflation [Ch. 30] to re-anchor them at substantial output cost.

### The 2021–22 Inflation Surge

The COVID-19 inflation episode of 2021–22 involved simultaneous adverse supply shocks (supply chain disruptions, semiconductor shortages, the Russian invasion of Ukraine) and favorable demand shocks (extraordinary fiscal transfers, pent-up demand). The AS-AD framework organizes the episode: multiple leftward SRAS shifts combined with a large rightward AD shift produced inflation substantially above target. The simultaneous occurrence of supply and demand shocks complicated the real-time policy diagnosis: the Fed's early characterization of inflation as "transitory" (primarily supply-driven and self-resolving) was inconsistent with the magnitude of the AD shift from fiscal stimulus.

---

## 7.6 The Dynamic Extension: New Keynesian AS–AD

The static AS–AD framework has several important limitations for policy analysis: it treats the price level rather than the inflation rate as the key nominal variable, does not model expectations formation explicitly, and provides no propagation dynamics. A minimal dynamic extension replaces price levels with inflation rates and specifies a central bank reaction function. The three-equation **New Keynesian baseline**:

$$\hat{\pi}_t = \beta\,\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\,\hat{x}_t + u_t \quad \text{(NKPC)}$$

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\pi_{t+1}] - r_t^n) \quad \text{(NK-IS)}$$

$$i_t = \bar{r} + \pi^* + \phi_\pi(\pi_t - \pi^*) + \phi_Y\hat{x}_t \quad \text{(Taylor Rule)}$$

where $\hat{x}_t = Y_t - \bar{Y}_t$ is the output gap, $r_t^n$ is the natural rate of interest, and $u_t$ is a cost-push shock. This system — derived from first principles in Chapters 9 and 10 — is the modern analytical successor to the AS-AD model. It makes the expectations formation explicit through the forward-looking NKPC and NK-IS, introduces the policy reaction function through the Taylor rule, and generates fully dynamic impulse responses to shocks of different types.

---

## Chapter Summary

- The **AD curve** is downward sloping through three channels — real balance (Pigou), interest rate (Keynes), and exchange rate (Mundell-Fleming) effects; the interest rate channel dominates. AD shifts with fiscal policy (multiplier $\times \Delta G$ or $\times \Delta T$), monetary policy ($\Delta M$), and autonomous expenditure changes.

- The **LRAS curve** is vertical at potential output $\bar{Y}$, shifting only when underlying productive capacity changes (capital, labor, technology). The **SRAS** is upward sloping: the sticky-wage model generates it through real wage compression at fixed nominal wages; the Lucas model generates it through rational signal extraction under imperfect information.

- Short-run equilibrium intersects AD and SRAS; long-run equilibrium additionally requires $P_t = P_t^e$, forcing $Y_t = \bar{Y}_t$. Adjustment from short to long run occurs through nominal wage revision, with downward nominal wage rigidity creating asymmetry between inflationary (fast) and recessionary (slow) adjustment.

- **Demand shocks** generate co-movement of output and prices (both move together); **supply shocks** generate stagflation (output and prices move in opposite directions). The 1970s oil shocks and the 2021–22 inflation episode both illustrated supply-shock dynamics in the AS-AD framework.

- The **dynamic NK extension** replaces price levels with inflation rates, adds forward-looking expectations, and closes with a Taylor rule — making the model suitable for quantitative policy analysis.

---

*Next: Chapter 8 — The Keynesian Cross and Multiplier Effect*
