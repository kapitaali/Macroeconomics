# Chapter 16 — Rational Expectations and New Classical Economics: The Role of Anticipations

> *"The theory of economic policy as it has been applied in practice assumes that the government can 'fool' the public. But if people are rational, they won't be fooled."*
> — Thomas Sargent, *Rational Expectations and Inflation*, 1986

---

Chapter 15 surveyed the psychological and empirical arguments against full rationality. This chapter takes the complementary perspective: it develops the case for rational expectations as an analytical discipline, derives its most consequential macroeconomic implications, and examines those implications' empirical support and theoretical limits. The chapter proceeds through three arguments that together constitute the New Classical research program: the Lucas imperfect information model and its critique of the Phillips curve; the policy ineffectiveness proposition; and the time inconsistency problem that produces an inflationary bias in discretionary monetary policy. Each of these arguments has permanently changed how macroeconomists think about policy, regardless of how one evaluates the precise empirical magnitude of the effects they identify.

---

## 16.1 The Lucas Imperfect Information Model

Lucas (1972, 1973) derived a short-run output-inflation relationship — superficially resembling the Phillips curve — without relying on any nominal rigidity, purely from rational agents who are temporarily imperfectly informed about the aggregate price level. The key ingredient is **imperfect information**: producers observe their own output prices precisely but observe the aggregate price level only with a lag or with error.

The model's mechanism: when a firm sees the price of its product rise, it faces a signal extraction problem. The price rise might reflect an increase in relative demand for its product — a genuine signal that it should expand production — or it might reflect an economy-wide nominal inflation that has raised all prices proportionally, which contains no real signal about relative profitability. A rational firm will increase output only if it believes the price rise is partly a real relative price signal.

Under joint normality of aggregate and relative price shocks, the optimal **Bayesian signal extraction** gives the firm's estimate of the general price level:

$$\mathbb{E}_t[P_t \mid p_{it}] = \frac{\sigma_z^2}{\sigma_\eta^2 + \sigma_z^2}\,p_{it} + \frac{\sigma_\eta^2}{\sigma_\eta^2 + \sigma_z^2}\,\bar{P}_t^e,$$

where $\sigma_\eta^2$ is the variance of aggregate nominal shocks and $\sigma_z^2$ is the variance of relative price shocks. The weight $\sigma_z^2/(\sigma_z^2 + \sigma_\eta^2)$ placed on the observed price signal measures the signal-to-noise ratio: when most price variability is relative ($\sigma_z^2$ large relative to $\sigma_\eta^2$), any observed price rise is likely to be a real signal and output responds strongly. When most variability is aggregate nominal ($\sigma_\eta^2$ large), price changes are mostly inflationary noise and output responds little.

Aggregating across all producers yields the **Lucas aggregate supply curve**:

$$y_t = \bar{y}_t + b(P_t - \mathbb{E}_{t-1}[P_t]),$$

where $b = \gamma\sigma_z^2/(\sigma_z^2 + \sigma_\eta^2)$ and $\mathbb{E}_{t-1}[P_t]$ is the prior expectation of the aggregate price level. Output deviates from potential only when the actual price level surprises agents — when $P_t \neq \mathbb{E}_{t-1}[P_t]$.

The slope $b$ is not a structural constant: it **decreases with aggregate price volatility** $\sigma_\eta^2$. In high-inflation environments, where nominal aggregate disturbances dominate the price signal, firms rationally attribute any price change to nominal noise rather than real demand — the signal-to-noise ratio for real information falls — and the real output response to a nominal shock shrinks. This prediction has been confirmed empirically across countries and time periods: the apparent output-inflation trade-off is steeper in low-inflation, stable-policy regimes and shallower in high-inflation environments.

This finding is the **Lucas critique of the Phillips curve**: the observed slope of the inflation-unemployment trade-off is not a structural parameter of the economy. It depends on the variance of monetary policy, which is itself a policy choice. Any attempt to exploit a historically estimated Phillips curve by accelerating inflation will, under rational expectations, simultaneously reduce the slope of the curve — the output gain per unit of surprise inflation falls as the inflation environment becomes noisier — until the predictable inflation produces no output gain at all.

---

## 16.2 The Policy Ineffectiveness Proposition

Sargent and Wallace (1975) derived a sweeping conclusion from combining the Lucas supply curve with rational expectations: **systematic monetary policy cannot affect real output**.

The argument is elegantly simple. Suppose the money supply follows the rule $m_t = \phi\mathbf{x}_{t-1} + \epsilon_t$, where $\phi\mathbf{x}_{t-1}$ is the systematic, rule-based component depending on any variables observable at $t-1$, and $\epsilon_t$ is an unforecastable random shock. Under rational expectations, agents know the rule $\phi$: before period $t$ begins, they already anticipate $\phi\mathbf{x}_{t-1}$. This anticipated component is therefore already embedded in $\mathbb{E}_{t-1}[P_t]$. From the Lucas supply curve, only the **unanticipated** component of the price level $P_t - \mathbb{E}_{t-1}[P_t]$ affects output. But the unanticipated component is driven only by the random innovation $\epsilon_t$ — the systematic component has been fully anticipated and thus fully reflected in expectations. The policy multiplier on the systematic component $\phi\mathbf{x}_{t-1}$ is therefore zero.

This proposition is logically tight and its empirical support is genuine at the level of its main prediction: evidence following Barro (1977) consistently shows that only **unanticipated** money growth is correlated with real output, while anticipated money growth is not. The proposition's limitations are equally important, however. It holds only under the joint assumptions of flexible prices and rational expectations. In New Keynesian models with nominal rigidities — staggered price and wage setting — anticipated monetary policy does have real effects precisely because not all prices can adjust immediately to any anticipated change; the firms and workers locked into existing contracts face a real change even when the aggregate price level adjusts toward its new level. The policy ineffectiveness proposition is therefore a statement about what would be true if prices were fully flexible, not a description of a world with pervasive nominal stickiness.

---

## 16.3 Time Inconsistency and the Inflationary Bias

The most policy-consequential contribution of the rational expectations approach is the analysis of **time inconsistency** in monetary policy. This problem arises even when the central bank has entirely benevolent intentions and correctly understands the economy.

**Definition (Time Inconsistency).** A policy plan is **time inconsistent** if the plan that appears optimal at date $t$, looking forward, is no longer optimal to implement at a later date $s > t$ once expectations formed at $t$ have been incorporated into contracts and private decisions. A policymaker who re-optimizes at each date — acting under **discretion** — will systematically deviate from previously announced plans, and rational agents who anticipate this will adjust their behavior accordingly, changing the outcomes of the announcement itself.

Kydland and Prescott (1977) identified time inconsistency as a fundamental problem for macro policy. The Barro–Gordon (1983) model gives the mechanism precisely. The central bank minimizes:

$$L = \frac{1}{2}\pi^2 + \frac{\lambda}{2}(y - y^*)^2,$$

where $y^*$ is the central bank's output target (set above the natural rate $\bar{y}$ due to pre-existing labor market distortions — taxes, unions, monopoly power — that keep the laissez-faire equilibrium inefficiently below the social optimum) and $\lambda > 0$ weights the output objective against inflation. The Lucas supply curve gives $y = \bar{y} + b(\pi - \pi^e)$: output rises when inflation exceeds expectations.

Under **discretion**, the central bank takes expected inflation $\pi^e$ as given — it has already been embedded in wage contracts before the central bank acts — and minimizes $L$ over $\pi$. The first-order condition is $\pi + \lambda b(\pi - \pi^e - y^* + \bar{y}) = 0$... simplifying, in rational expectations equilibrium where $\pi^e = \pi$:

$$\pi^D_{eq} = b\lambda(y^* - \bar{y}) > 0.$$

This is the **inflationary bias**: a positive equilibrium inflation rate emerges purely from the desire to exploit the short-run Phillips curve. The central bank's output target above the natural rate creates a standing temptation to inflate. Because rational agents anticipate this, they set $\pi^e = \pi^D_{eq}$, and in equilibrium the central bank validates their expectations. Output remains at $\bar{y}$ — the inflation produces no output benefit — but inflation is permanently positive. The discretionary equilibrium has higher inflation and the same output as the zero-inflation commitment equilibrium, representing a pure welfare loss.

The **commitment equilibrium**, in which the central bank credibly commits to $\pi = 0$ and agents accordingly set $\pi^e = 0$, achieves lower welfare loss. The tension is that this commitment is not credible under discretion: once agents have set $\pi^e = 0$ in wage contracts, the central bank faces the temptation to surprise them with positive inflation and push output temporarily above $\bar{y}$. Rational agents, anticipating this temptation, will not believe $\pi = 0$ unless the commitment is somehow institutionally enforced.

### Institutional Responses to the Time-Inconsistency Problem

The theoretical identification of the inflationary bias had direct practical consequences. It provided the intellectual foundation for a set of institutional reforms to monetary policy governance that began in the 1980s and accelerated in the 1990s:

**Central bank independence** (Alesina and Summers, 1993): insulating the central bank from day-to-day political pressure reduces the effective weight $\lambda$ placed on the output objective, since politicians facing elections have shorter horizons and stronger incentives to exploit short-run trade-offs. Empirically, cross-country evidence shows a negative correlation between central bank independence indices and average inflation, consistent with the model's prediction.

**Inflation targeting** (Svensson, 1997): making the inflation objective explicit, public, and institutionally binding changes the payoff structure for the policymaker, raising the reputational and legal costs of deviating from the target. New Zealand adopted the first formal inflation targeting framework in 1990; by 2020, over 40 countries had adopted variants of it.

**Rogoff's conservative central banker** (Rogoff, 1985): appointing a central banker with a lower weight on the output objective ($\lambda^{CB} < \lambda^{society}$) reduces the inflationary bias, though at the cost of suboptimally slow responses to output shocks. The optimal appointment trades off these two distortions.

**Walsh's optimal contracts** (Walsh, 1995): designing central banker compensation to depend on inflation performance can, in principle, eliminate the inflationary bias entirely without introducing output-stabilization suboptimality, by aligning the central banker's incentives directly with the social welfare function.

Each solution addresses the commitment problem differently, and each has different practical advantages and costs. Their common logic — that institutional design can substitute for credibility that discretionary policymakers cannot sustain through announcements alone — represents one of the practical payoffs of the rational expectations research program.

---

## 16.4 The Lucas Critique: A Deeper Look

The Lucas critique, introduced conceptually in Chapter 1 and applied to the Phillips curve above, deserves a more formal treatment as the methodological capstone of the New Classical program.

The formal statement: the parameters of any estimated reduced-form econometric relationship between policy instruments and macroeconomic outcomes are valid only conditional on the policy rule in effect during the estimation period. When the rule changes, agents update their optimal decision rules — their consumption behavior, their wage-setting behavior, their investment rules — and the estimated reduced-form coefficients change with them. Policy evaluation based on reduced-form estimates is therefore invalid for any policy change large enough to affect agents' expectations about the future policy regime.

Consider the estimated consumption function $C = a + b(Y - T)$. The marginal propensity to consume $b$ is not a structural parameter of household preferences; it is the equilibrium response of household consumption to income changes given the household's expectation about the future path of taxes and income. If the government shifts from funding deficits with future tax increases (Ricardian households expect to pay eventually, reducing $b$) to funding them with future monetization (households expect real income to rise, potentially raising $b$), the same estimated $b$ cannot be used for both. The Lucas critique is precisely this: behavioral parameters estimated under one regime cannot be used to predict behavior under a different regime.

The methodological response is the **structural approach**: building models from **deep parameters** — preference parameters such as $\sigma$ (the elasticity of intertemporal substitution) and $\beta$ (the discount factor), and technology parameters such as $\alpha$ (capital's share) — that are invariant to policy regime changes because they reflect preferences and production possibilities rather than equilibrium behavioral rules. These structural parameters can be estimated or calibrated from microeconomic evidence, and the resulting structural model can then be used for policy experiments. This is the rationale for the DSGE approach and its dominance in modern quantitative macroeconomics [Chs. 27–31 in the companion *Methods* volume].

The critique has practical limits worth acknowledging. Not all behavioral parameters shift with regime changes; some may be sufficiently stable across plausible policy environments that reduced-form estimation provides a reasonable first approximation. The Lucas critique is a warning about the conditions under which reduced-form evaluation is reliable, not an assertion that it is never reliable. Empirical work in the structural VAR tradition [Ch. 27] attempts to identify regime-invariant structural shocks rather than estimating regime-dependent reduced-form relationships, precisely to navigate this tension.

---

## Chapter Summary

- The **Lucas supply curve** $y_t = \bar{y}_t + b(P_t - \mathbb{E}_{t-1}[P_t])$ derives a temporary output-inflation relationship from rational signal extraction under imperfect information, without nominal rigidities. The slope $b$ decreases with nominal aggregate price variability — the signal-to-noise critique of the Phillips curve.

- The **policy ineffectiveness proposition** (Sargent and Wallace, 1975): under flexible prices and rational expectations, only unanticipated monetary policy affects real output. Anticipated systematic policy is fully embedded in price expectations and generates no surprise. The proposition holds for flexible-price models; New Keynesian nominal rigidities restore real effects of anticipated policy.

- **Time inconsistency** generates an **inflationary bias** $\pi^D_{eq} = b\lambda(y^* - \bar{y}) > 0$ in discretionary monetary policy: the central bank's output target above the natural rate creates a standing temptation to inflate, which rational agents anticipate, producing positive equilibrium inflation with no output gain.

- Institutional responses — central bank independence, inflation targeting, conservative central bankers, performance contracts — address the commitment problem by changing either the policymaker's incentives or the feasibility of deviating from announced policy.

- The **Lucas critique** states that reduced-form econometric parameters estimated under one policy regime cannot be used for policy evaluation under a different regime, because agents will revise their optimal decision rules when the regime changes. The structural DSGE approach responds by building from preference and technology parameters that are regime-invariant.

---

*Next: Chapter 17 — The Goods Market*
