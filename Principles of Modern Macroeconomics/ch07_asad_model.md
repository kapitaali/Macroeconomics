# Chapter 7 — The Aggregate Demand–Aggregate Supply Model: The Workhorse of Macroeconomics

> *"In the long run, prices and wages adjust to clear markets. In the short run, they do not."*
> — N. Gregory Mankiw

---

The IS–LM model (developed in Chapter 9) characterizes short-run equilibrium in the goods and money markets for a given price level. But the price level is itself an endogenous variable: it responds to the balance between aggregate demand and the economy's productive capacity. The **aggregate demand–aggregate supply (AS–AD) framework** closes this gap by combining a theory of how demand determines output (for a given price level) with a theory of how supply — the behavior of price and wage setters — determines the price level (for a given level of demand). Together, these two curves determine both real output and the price level simultaneously.

The AS–AD framework is a workhorse rather than a frontier model. Its simplicity is its strength: it organizes thinking about short-run stabilization, supply shocks, and the distinction between demand-pull and cost-push inflation in a way that is directly applicable to policy analysis. Its limitations — the absence of explicit intertemporal optimization, rational expectations, and a financial sector — motivate the more sophisticated models of Parts III and beyond.

---

## 7.1 Aggregate Demand

**Definition (Aggregate Demand Curve).** The **aggregate demand (AD) curve** is the locus of price level $P$ and real output $Y$ combinations at which the goods market (IS equilibrium) and the money market (LM equilibrium) are simultaneously in equilibrium. It is derived by asking: if the price level changes, how does the equilibrium output level change, holding all other exogenous variables (government spending, money supply, taxes, expectations) fixed?

To derive the slope, recall the IS condition ($Y = C(Y-T) + I(r) + G$) and the LM condition ($M/P = L(Y, r + \pi^e)$). When $P$ rises, the real money supply $M/P$ falls. This shifts LM to the left, raising the interest rate for any given $Y$, which reduces investment and hence output. Formally, differentiating the LM equation:

$$\mathrm{d}r = \frac{1}{L_r}\!\left(-\frac{M}{P^2}\mathrm{d}P - L_Y\,\mathrm{d}Y\right).$$

Substituting into the differentiated IS curve and solving:

$$\left.\frac{\mathrm{d}Y}{\mathrm{d}P}\right|_{AD} = \frac{I'(r)\cdot M/P^2}{L_Y I'(r) - L_r} < 0.$$

The AD curve is **downward sloping**: a higher price level reduces real output through three distinct channels.

The first is the **real-balance effect** (Pigou effect): a higher price level reduces real wealth $M/P$ for holders of nominal money balances, reducing their consumption. This effect is typically small because money holdings are a small fraction of household wealth.

The second is the **interest rate effect** (Keynes effect): the dominant channel. A higher price level reduces the real money supply, pushing up interest rates, which in turn reduce interest-sensitive investment spending and consumer durable purchases.

The third is the **exchange rate effect** (Mundell–Fleming effect), operative in open economies only: a higher domestic price level makes domestic goods relatively more expensive than foreign goods, reducing net exports. This channel is developed in Chapter 21.

### Shifts in the AD Curve

Any change in an exogenous variable other than $P$ that affects IS or LM will shift the entire AD curve. The most important shifters are:

- **Fiscal policy**: an increase in $G$ or a cut in $T$ raises autonomous expenditure, shifting IS right and AD right (at any given $P$, equilibrium $Y$ is higher).
- **Monetary policy**: an increase in $M$ shifts LM right, reducing $r$ at any given $Y$, which stimulates investment and shifts AD right.
- **Consumer or investor confidence**: an autonomous increase in $C$ or $I$ shifts IS and AD right.
- **Foreign income or exchange rate**: for open economies, higher foreign income raises exports, shifting AD right.

---

## 7.2 Aggregate Supply

The aggregate supply side of the model describes how much output producers wish to supply at any given price level. This depends critically on the time horizon, because the behavior of price and wage setters differs fundamentally in the short run (when nominal rigidities bind) and the long run (when all prices are flexible).

### Long-Run Aggregate Supply

In the long run, all nominal prices and wages are fully flexible. Firms will produce at the level determined by the real productive capacity of the economy — the quantity of capital, labor, and technology available — irrespective of the price level. This is precisely the definition of potential output $\bar{Y}$ from Chapter 1.

**Definition (Long-Run Aggregate Supply).** The **long-run aggregate supply (LRAS) curve** is a vertical line at $Y = \bar{Y}$, where $\bar{Y}$ is potential output. It is vertical because in the long run, a higher price level simply leads to proportionally higher nominal wages, leaving the real wage and hence employment and output unchanged. The LRAS shifts only when the underlying determinants of potential output change: when the capital stock grows, the labor force expands or becomes more educated, or technology improves.

### Short-Run Aggregate Supply: The Sticky-Wage Model

In the short run, nominal wages are typically fixed by contracts or by social norms that make rapid nominal wage cuts psychologically and institutionally difficult. Suppose the nominal wage is fixed at $\bar{W}$ for one period. A profit-maximizing firm hires labor until the marginal product of labor equals the real wage:

$$F_L(K, N) = \frac{\bar{W}}{P}.$$

When $P$ rises, the real wage $\bar{W}/P$ falls, making it profitable for firms to hire more workers and produce more. Hence output rises with $P$ in the short run. Solving for $N$ and substituting into the production function:

$$Y^s = F\!\left(K,\, N\!\left(\frac{\bar{W}}{P}\right)\right), \quad \frac{\partial Y^s}{\partial P} > 0.$$

The short-run aggregate supply (SRAS) curve is **upward sloping**: higher prices, by eroding real wages, induce more production.

### Short-Run Aggregate Supply: The Lucas Imperfect-Information Approach

A different foundation for the upward-sloping SRAS, based entirely on rational behavior rather than nominal rigidities, comes from Lucas (1973). In the Lucas model, firms observe their own output price perfectly but observe the general price level only with a lag. When a firm sees its price rise, it cannot immediately determine whether this reflects a relative price increase (genuine demand for its product has increased, warranting more production) or a general price level increase (all prices have risen uniformly, warranting no change in production). The firm rationally attributes the price rise partly to each explanation, and since it believes part of the increase is a real relative price improvement, it increases production:

$$Y_t = \bar{Y}_t + \alpha(P_t - P_t^e), \quad \alpha > 0,$$

where $P_t^e \equiv \mathbb{E}_{t-1}[P_t]$ is the expected price level formed before period-$t$ prices are observed. This is the **Lucas supply function**: output exceeds potential only when the actual price level exceeds the expected price level — only when agents are surprised. Under rational expectations, only unanticipated changes in the price level (driven by unanticipated monetary shocks) can raise output above potential.

---

## 7.3 Short-Run and Long-Run Equilibrium

Short-run equilibrium is the intersection of AD and SRAS. Using the log-linear Lucas supply function:

$$\text{SRAS: } Y_t = \bar{Y}_t + \alpha(P_t - P_t^e)$$
$$\text{AD: } Y_t = \Phi(P_t, M_t, G_t, T_t)$$

where $\Phi$ is decreasing in $P_t$ and increasing in $M_t$ and $G_t$. These two equations simultaneously determine $Y_t$ and $P_t$.

Long-run equilibrium adds the requirement that expectations are fulfilled: $P_t = P_t^e$. This forces $Y_t = \bar{Y}_t$ — the economy is at potential. The price level is then determined by the intersection of AD with the vertical LRAS.

The **adjustment mechanism** from short-run to long-run equilibrium is the gradual revision of price expectations. If the economy is below potential ($Y_t < \bar{Y}_t$), there is excess supply of labor and goods. Wage growth slows; firms reduce prices relative to what was expected; $P_t < P_t^e$, which by the SRAS equation pushes $Y_t$ back toward $\bar{Y}_t$. The adjustment continues until $Y_t = \bar{Y}_t$. The speed of this adjustment depends on how quickly expectations and nominal contracts respond to current conditions — in other words, on the degree of nominal rigidity in the economy.

---

## 7.4 The Macroeconomic Effects of Shocks

One of the central uses of the AS–AD framework is the analysis of how the economy responds to exogenous shocks. Two classes of shocks deserve attention.

A **positive demand shock** — an increase in government spending, a monetary expansion, or a surge in consumer confidence — shifts the AD curve rightward. In the short run, before prices fully adjust, output rises above potential and the price level rises. In the long run, wages and prices adjust upward, shifting SRAS left until output returns to $\bar{Y}$ at a permanently higher price level. The real effects of the shock are entirely temporary; the nominal effects (higher price level) are permanent. This is the sense in which money is "neutral in the long run but not in the short run."

An **adverse supply shock** — an oil price spike, a natural disaster, a pandemic lockdown — shifts the SRAS curve leftward: at any given price level, firms produce less because their costs have risen. The short-run equilibrium exhibits **stagflation**: lower output and a higher price level simultaneously. This combination was precisely what the 1970s oil shocks produced, and it was inexplicable under the purely demand-oriented models that dominated policy at the time. The AS–AD framework, by distinguishing demand and supply disturbances and their different implications for the inflation-output relationship, was the theoretical advance that made sense of the stagflation episode.

---

## 7.5 The Dynamic Extension: New Keynesian AS–AD

The static AS–AD framework has several important limitations for policy analysis. It treats the price level rather than the inflation rate as the key nominal variable, which makes it awkward to analyze monetary policy that targets inflation. It does not model expectations formation explicitly. And it is entirely static, with no dynamics beyond the qualitative story about medium-run adjustment.

A minimal dynamic extension replaces $P_t$ with $\pi_t$ and specifies how the central bank sets the interest rate. The three-equation **New Keynesian baseline** system is:

$$\hat{\pi}_t = \beta\,\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\,\hat{x}_t \quad \text{(New Keynesian Phillips Curve)}$$

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\pi_{t+1}] - r_t^n) \quad \text{(New Keynesian IS Curve)}$$

$$i_t = \bar{r} + \pi^* + \phi_\pi(\pi_t - \pi^*) + \phi_Y\hat{x}_t \quad \text{(Taylor Rule)}$$

Here $\hat{x}_t = Y_t - \bar{Y}_t$ is the output gap, $r_t^n$ is the natural (Wicksellian) rate of interest, and $\pi^*$ is the inflation target. The Taylor rule specifies how the central bank adjusts the nominal interest rate in response to inflation and output deviations from target. This three-equation system is derived from first principles in Chapters 9 and 10 and analyzed as the framework for monetary policy in Chapter 23.

---

*Next: Chapter 8 — The Keynesian Cross and Multiplier Effect*
