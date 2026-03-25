# Chapter 10 — The Phillips Curve: Inflation and Unemployment

> *"There is always a temporary trade-off between inflation and unemployment; there is no permanent trade-off."*
> — Milton Friedman, Nobel Memorial Lecture, 1977

---

The aggregate supply side of the economy — the relationship between inflation and economic activity — is one of the most important and most contested areas of macroeconomics. The AS–AD model of Chapter 7 posited an upward-sloping short-run supply curve and left its slope as a parameter. The Phillips curve literature attempts to give that slope empirical content and theoretical foundations. It also has a long history of being misunderstood, misused by policymakers, and then progressively refined as the errors became apparent. Few episodes in the history of economic thought so vividly illustrate the difference between a robust structural relationship and a historically contingent correlation.

---

## 10.1 The Original Empirical Relationship

In 1958, the New Zealand economist A.W. Phillips published a paper documenting an empirical regularity in nearly a century of British data. Using annual observations on money wage inflation and unemployment from 1861 to 1957, he found a stable, negatively sloped, nonlinear relationship: when unemployment was low, wages rose rapidly; when unemployment was high, wages grew slowly or even fell.

Phillips estimated the relationship:

$$\hat{w}_t = f(u_t), \quad \hat{w}_t \equiv \frac{\mathrm{d}\ln W_t}{\mathrm{d}t}, \quad f' < 0,$$

where $\hat{w}_t$ is the rate of change of nominal wages and $u_t$ is the unemployment rate. Samuelson and Solow (1960) then translated this into a relationship between price inflation $\pi_t$ and unemployment for the United States, using the observation that prices are set as a markup over wages. The resulting **Phillips curve** appeared to offer policymakers a permanent menu of choices: accept a little more inflation and enjoy a little less unemployment, or vice versa.

The policy implications were immediately seized upon. If the inflation–unemployment trade-off were stable and exploitable, governments could choose their preferred point on the curve — perhaps 4% unemployment at the cost of 3% inflation — and maintain it indefinitely through demand management. The 1960s, during which this view dominated policy, appeared initially to confirm it: expansionary fiscal and monetary policy in the United States reduced unemployment while inflation rose only modestly. The curve seemed to be holding.

---

## 10.2 The Friedman–Phelps Critique and the Natural Rate

The appearance of stability was illusory, and two economists identified the theoretical reason why before the empirical breakdown became obvious. Milton Friedman (1968) and Edmund Phelps (1968) each published, independently, what is now called the **expectations-augmented Phillips curve**.

Their argument begins with a simple observation: the original Phillips curve relates nominal wage growth to unemployment. But workers and firms negotiate wages in anticipation of the real purchasing power those wages will deliver. A worker who demands a 5% wage increase when she expects 5% inflation has no reason to work more; she is merely trying to maintain her real wage. The original Phillips curve, by relating nominal wage growth directly to unemployment without accounting for inflation expectations, implicitly assumed that workers have static expectations — they always expect zero inflation. This assumption cannot hold indefinitely in an economy where inflation is persistently positive.

The correct relationship, Friedman and Phelps argued, is:

$$\pi_t = \pi_t^e - \alpha(u_t - u^*) + \epsilon_t, \quad \alpha > 0,$$

where $\pi_t^e$ is the **expected inflation rate** — the rate of inflation that workers and firms expect when they negotiate wages and prices — and $u^*$ is the **natural rate of unemployment** (defined in Chapter 1 as the unemployment rate consistent with stable, non-accelerating inflation). The term $\epsilon_t$ captures supply-side disturbances such as oil price shocks that move inflation independently of demand conditions. This is the **expectations-augmented Phillips curve (EAPC)**.

The EAPC has a critical implication: the inflation–unemployment trade-off is not between the level of inflation and the level of unemployment; it is between *unexpected* inflation and unemployment deviations from the natural rate. If policymakers try to permanently hold unemployment below $u^*$ by maintaining inflation higher than expected, they will find that inflation expectations rise to match actual inflation — closing the gap between $\pi_t^e$ and $\pi_t$ — and the trade-off evaporates.

To see this formally, suppose expectations are **adaptive**: agents expect this period's inflation to equal last period's, $\pi_t^e = \pi_{t-1}$. Then the EAPC becomes:

$$\pi_t - \pi_{t-1} = -\alpha(u_t - u^*) + \epsilon_t.$$

This is called the **accelerationist Phillips curve**, because it predicts that holding unemployment below $u^*$ causes inflation not merely to be elevated but to continuously *accelerate*. Any unemployment rate below $u^*$ is consistent with any constant inflation rate, as long as the central bank is willing to raise inflation each period. But the inflation rate required to maintain $u_t < u^*$ grows without bound over time. The only unemployment rate consistent with a stable, non-accelerating inflation rate is $u^* =$ the NAIRU (Non-Accelerating Inflation Rate of Unemployment).

The stagflation of 1973–75 — in which both unemployment and inflation rose simultaneously following the OPEC oil shock — provided dramatic empirical confirmation of the Friedman–Phelps framework and a decisive refutation of the original stable trade-off interpretation. A negative supply shock (captured by $\epsilon_t > 0$ in the EAPC) shifts the curve so that any given unemployment rate is associated with higher inflation. The 1960s expansion had pushed inflation expectations upward, making the apparent trade-off increasingly unfavorable; the oil shock then pushed both unemployment and inflation higher simultaneously, a combination the original Phillips curve predicted was impossible.

---

## 10.3 The New Keynesian Phillips Curve

The expectations-augmented Phillips curve of Friedman and Phelps is an improvement over the original, but it has a theoretical limitation: it derives inflation expectations from an assumed rule of thumb (adaptive expectations) rather than from optimal behavior. By the 1980s, the rational expectations revolution (Chapter 16) had made clear that any empirical relationship whose validity depends on agents making systematic forecast errors is suspect.

The **New Keynesian Phillips Curve (NKPC)** is the product of the New Keynesian research program's effort to derive the inflation–output relationship from explicit microfoundations — from the optimization problems of individual price-setting firms — while maintaining rational expectations throughout. The key micro ingredient is Calvo's (1983) model of staggered price setting, introduced in Chapter 2.

Recall the Calvo setup: in each period, a randomly chosen fraction $1-\theta$ of firms can reset their price optimally, while the remaining fraction $\theta$ must keep their price unchanged. A firm that resets at date $t$ is in effect setting a price that may persist for several periods; it therefore chooses $p_t^*$ to maximize the expected present discounted value of future profits over the spell during which the price will be in effect:

$$\max_{p_t^*} \sum_{k=0}^\infty (\beta\theta)^k\,\mathbb{E}_t\!\left[\lambda_{t+k}\!\left(\frac{p_t^*}{P_{t+k}} - mc_{t+k}\right) y_{t+k|t}\right],$$

where $\lambda_{t+k}$ is the household's marginal utility of nominal income at $t+k$ (reflecting the fact that firms are owned by households), $mc_{t+k}$ is the real marginal cost at $t+k$, and $y_{t+k|t}$ is the demand the firm faces at $t+k$ given it last reset its price at $t$. The discount factor $\beta\theta$ combines the household's discount factor $\beta$ with the probability $\theta$ that the price set today will still be in effect next period.

The first-order condition for $p_t^*$, log-linearized around a zero-inflation steady state and expressed as a deviation from the steady state (hat notation):

$$\hat{p}_t^* = (1-\beta\theta)\sum_{k=0}^\infty (\beta\theta)^k\,\mathbb{E}_t\bigl[\widehat{mc}_{t+k} + \hat{P}_{t+k}\bigr].$$

This says the optimal reset price is the discounted average of expected future nominal marginal costs — the firm is setting a price today that will be appropriate across all future periods in which it remains in effect.

The aggregate price level evolves as $\hat{P}_t = \theta\hat{P}_{t-1} + (1-\theta)\hat{p}_t^*$: a weighted average of last period's price level (carried over by the $\theta$ non-adjusters) and the new optimal price. Combining this with the expression for $\hat{p}_t^*$ and using $\hat{\pi}_t = \hat{P}_t - \hat{P}_{t-1}$, after some algebra the **New Keynesian Phillips Curve** emerges:

$$\hat{\pi}_t = \beta\,\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\,\widehat{mc}_t, \qquad \kappa = \frac{(1-\theta)(1-\beta\theta)}{\theta}.$$

The slope coefficient $\kappa$ is a decreasing function of price stickiness $\theta$: when prices are very sticky (large $\theta$), the NKPC is flat — even large movements in marginal costs generate little inflation, because most firms are locked into their existing prices. When prices are flexible (small $\theta$), $\kappa$ is large and marginal cost movements translate immediately into inflation.

Using the fact that under specific functional forms (log-linear preferences, Cobb–Douglas production), real marginal cost is proportional to the **output gap** $\hat{x}_t \equiv Y_t - \bar{Y}_t$:

$$\hat{\pi}_t = \beta\,\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\,\hat{x}_t.$$

This is the compact form of the NKPC. By iterating forward — substituting $\mathbb{E}_t[\hat{\pi}_{t+1}] = \beta\mathbb{E}_t[\hat{\pi}_{t+2}] + \kappa\mathbb{E}_t[\hat{x}_{t+1}]$ and so on — we obtain:

$$\hat{\pi}_t = \kappa\sum_{k=0}^\infty \beta^k\,\mathbb{E}_t[\hat{x}_{t+k}].$$

This is a striking result: inflation today is the discounted sum of all *future* expected output gaps. It is a **present-value formula for inflation**, analogous to the formula for stock prices as the present value of future dividends. The implication is that inflation is entirely forward-looking — it depends not on past output gaps (as in the accelerationist model) but on what the economy is expected to do in the future. A credible commitment to return the economy to potential output tomorrow should reduce inflation today, even without any current period of slack.

This forward-looking property is both the NKPC's greatest theoretical strength and its principal empirical challenge. It implies that disinflations can be achieved costlessly if they are fully credible — a prediction that seems at odds with the large output losses observed in historical disinflation episodes. It also implies that the "sacrifice ratio" — the output loss per unit of inflation reduction — should be zero under credible policy, while empirical estimates cluster at 1.4–2.8. The resolution of this tension motivates the hybrid NKPC discussed next.

---

## 10.4 The Hybrid Phillips Curve

The pure NKPC's purely forward-looking structure implies that past inflation is irrelevant to current inflation — a prediction that conflicts with the pervasive inertia in inflation observed in the data. Inflation in the United States, for instance, is highly autocorrelated: high inflation last year is a strong predictor of high inflation this year, even after controlling for current output gaps and expected future inflation.

Galí and Gertler (1999) proposed the **hybrid New Keynesian Phillips Curve** to reconcile the theoretical forward-looking structure with observed inertia. They introduce a fraction $\omega$ of firms that, instead of optimizing forward-looking price-setting, use a **backward-looking rule of thumb**: they set their price equal to the previous period's average reset price plus last period's inflation. Firms use backward-looking rules not because they are irrational but because optimization is costly and rules of thumb that track recent inflation adequately approximate the optimum in a low-inflation environment.

The hybrid NKPC is:

$$\pi_t = \frac{\omega}{\phi}\,\pi_{t-1} + \frac{(1-\omega)\beta}{\phi}\,\mathbb{E}_t[\pi_{t+1}] + \frac{\kappa(1-\omega)}{\phi}\,\hat{x}_t,$$

where $\phi = \omega + (1-\omega)(1 + \beta\omega)/(1-\omega\beta)$. When $\omega = 0$ (no backward-looking firms), the hybrid reduces to the pure NKPC. When $\omega \to 1$ (all firms are backward-looking), it approaches the expectations-augmented PC with adaptive expectations.

Galí and Gertler estimate $\omega \approx 0.3$ for the U.S., implying that forward-looking behavior dominates: roughly 70% of price-setting is forward-looking, with 30% backward-looking. This generates moderate inflation inertia without completely eliminating the forward-looking component.

---

## 10.5 The Sacrifice Ratio

Disinflation — the deliberate reduction of inflation from a higher to a lower level — typically requires a period of below-potential output and above-natural unemployment. The output cost of disinflation is summarized by the **sacrifice ratio**.

**Definition (Sacrifice Ratio).** The sacrifice ratio is the cumulative output loss, expressed as a percentage of annual GDP, required to reduce inflation by one percentage point:

$$SR = \frac{\sum_{t=0}^T (Y_t^{base} - Y_t)/\bar{Y}}{\pi_0 - \pi^*},$$

where $Y_t^{base}$ is the counterfactual output path in the absence of disinflation, $Y_t$ is actual output during the disinflation, and $\pi_0 - \pi^*$ is the total inflation reduction achieved.

The theoretical sacrifice ratio depends on the model of expectations. Under the accelerationist EAPC with adaptive expectations, a permanent one-unit reduction in inflation requires a permanent one-unit increase in unemployment above the natural rate, giving $SR = 1/\alpha$ from the EAPC equation. Under the pure NKPC with perfectly credible policy, $SR \approx 0$: a credible announcement that the central bank will tighten policy until inflation falls immediately shifts expectations down, reducing inflation without requiring any period of slack. The reality lies between these extremes.

Ball (1994) documents sacrifice ratios across 28 historical disinflation episodes in OECD countries, finding wide variation from 0.5 to 4, with the U.S. Volcker disinflation (1979–83) yielding $SR \approx 2$. The variation correlates with the speed of disinflation (faster disinflations have lower sacrifice ratios, consistent with the forward-looking component) and with institutional credibility (central banks with stronger inflation-fighting reputations face lower costs). The intermediate values — between zero and $1/\alpha$ — are consistent with the hybrid NKPC's mixture of forward- and backward-looking price-setters, combined with the observation that monetary policy announcements are only partially credible.

---

*Next: Chapter 11 — Consumption Theory*
