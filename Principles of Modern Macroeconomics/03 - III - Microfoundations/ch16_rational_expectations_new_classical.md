# Chapter 16 — Rational Expectations and New Classical Economics: The Role of Anticipations

---

Chapter 15 described the spectrum of assumptions about expectations formation and highlighted the empirical and psychological arguments against full rationality. This chapter takes the opposite tack: it develops the case for rational expectations as a modeling discipline, derives its most important macroeconomic implications, and examines the major critiques. Whatever one's view of the psychological realism of rational expectations, its implications for macroeconomic policy are too important to bypass.

---

## 16.1 The Lucas Imperfect Information Model

Lucas (1972, 1973) used rational expectations to derive a short-run inflation–output trade-off from a microfounded model with no nominal rigidities. The key ingredient is **imperfect information**: producers observe their own output price precisely but observe the general price level only with a lag. When a firm sees its price rise, it cannot immediately tell whether this reflects an increase in relative demand for its product (a genuine signal to produce more) or a general inflation shock (no real signal at all).

Under joint normality of aggregate and relative price shocks, **Bayesian signal extraction** gives the optimal forecast of the general price level:

$$\mathbb{E}_t[P_t \mid p_{it}] = \frac{\sigma_z^2}{\sigma_\eta^2 + \sigma_z^2}\,p_{it} + \frac{\sigma_\eta^2}{\sigma_\eta^2 + \sigma_z^2}\,\bar{P}_t^e,$$

where $\sigma_\eta^2$ is the variance of aggregate nominal shocks and $\sigma_z^2$ is the variance of relative price shocks. The ratio $\sigma_z^2/(\sigma_\eta^2 + \sigma_z^2)$ is the weight the firm places on the observed price signal versus the prior expectation. Aggregating across producers yields the **Lucas aggregate supply curve**:

$$y_t = \bar{y}_t + \frac{\gamma\sigma_z^2}{\sigma_z^2 + \sigma_\eta^2}(P_t - \mathbb{E}_{t-1}[P_t]).$$

Output deviates from potential only when the actual price level surprises agents — when $P_t \neq \mathbb{E}_{t-1}[P_t]$. The slope $b \equiv \gamma\sigma_z^2/(\sigma_z^2 + \sigma_\eta^2)$ is the key parameter: it decreases with aggregate price volatility $\sigma_\eta^2$. In high-inflation environments, firms rationally attribute any price signal to nominal disturbances rather than to real demand — the signal-to-noise ratio for real signals falls — so the real output response to a nominal shock is smaller. This is the **Lucas critique of the Phillips curve**: the slope of the trade-off is not structural but depends on the variance of monetary policy, which is itself a policy choice.

---

## 16.2 The Policy Ineffectiveness Proposition

Sargent and Wallace (1975) derived a sweeping conclusion from the Lucas model combined with rational expectations: **systematic monetary policy cannot affect real output**.

Suppose the money supply follows the rule $m_t = \phi\mathbf{x}_{t-1} + \epsilon_t$, where $\mathbf{x}_{t-1}$ is any vector of variables observable at $t-1$ and $\epsilon_t$ is an unforecastable random shock. Under rational expectations, agents know the rule $\phi$: they fully anticipate the systematic component $\phi\mathbf{x}_{t-1}$ before period $t$ begins. From the Lucas supply curve, only $P_t - \mathbb{E}_{t-1}[P_t]$ — the unanticipated component of the price level — affects output. Since the systematic component of money is fully anticipated, it is already embedded in $\mathbb{E}_{t-1}[P_t]$, and the price level surprise is driven only by $\epsilon_t$. The policy multiplier of the systematic component $\phi\mathbf{x}_{t-1}$ on real output is therefore zero.

This proposition is strong and controversial. It holds in models with flexible prices and rational expectations; it does not hold in New Keynesian models with nominal rigidities. But even in New Keynesian models, the rational expectations assumption substantially reduces the real effects of anticipated policy changes compared to adaptive expectations models — and this is the empirically relevant observation.

---

## 16.3 Time Inconsistency and the Inflationary Bias

Perhaps the most policy-relevant contribution of the rational expectations, game-theoretic approach to macroeconomics is the analysis of **time inconsistency** in monetary policy.

**Definition (Time Inconsistency).** A policy plan is **time inconsistent** if the plan that appears optimal at date $t$ (looking forward) is no longer optimal to implement at a later date $s > t$ once expectations formed at $t$ have been incorporated into contracts and decisions. A policymaker operating under discretion — able to re-optimize at each date — will deviate from the originally announced plan, and rational agents will anticipate this, affecting the outcomes of the original announcement.

Kydland and Prescott (1977) identified time inconsistency as a fundamental problem for monetary policy. The Barro–Gordon (1983) model makes the mechanism precise. The central bank minimizes the loss function:

$$L = \frac{1}{2}\pi^2 + \frac{\lambda}{2}(y - y^*)^2,$$

where $y^* > \bar{y}$ is the central bank's output target (set above the natural rate due to labor market distortions — taxes, monopoly power — that make the laissez-faire equilibrium inefficient). The Lucas supply curve gives $y = \bar{y} + b(\pi - \pi^e)$.

Under **discretion**, the central bank takes expected inflation $\pi^e$ as given (formed before the central bank acts) and minimizes $L$ over $\pi$. The first-order condition: $\pi = b\lambda(y^* - \bar{y}) + \pi^e$. In rational expectations equilibrium, $\pi^e = \pi$, which gives:

$$\pi^{D}_{eq} = b\lambda(y^* - \bar{y}) > 0.$$

This positive equilibrium inflation — the **inflationary bias** — emerges purely from the attempt to exploit the short-run Phillips curve. The central bank's desire to raise output above the natural rate leads it to choose positive inflation in equilibrium, but because agents anticipate this, output remains at $\bar{y}$ and the only result is permanently higher inflation with no output gain.

Under **commitment** to a rule $\pi = 0$, the central bank achieves $\pi^* = 0$, $y = \bar{y}$, and lower welfare loss. The commitment equilibrium Pareto dominates the discretionary equilibrium, yet under discretion the central bank always has an incentive to deviate: once $\pi^e = 0$ is embedded in wage contracts, the central bank can temporarily raise output to $y^*$ by choosing $\pi > 0$. This is the time-inconsistency trap.

The practical solutions to the inflationary bias problem include: **central bank independence** (insulating the central bank from political pressure to pursue the output objective); **inflation targeting** (making the inflation objective explicit and legally binding); **Rogoff's conservative central banker** (Rogoff, 1985: appointing a central banker with $\lambda < \lambda^{society}$, so she cares relatively more about inflation); and **Walsh's optimal incentive contracts** (Walsh, 1995: designing contracts that make central bank compensation contingent on inflation performance). Each of these institutional arrangements attempts to solve the commitment problem by changing the incentives of the policymaker or by constraining her choices.

---

## 16.4 The Lucas Critique: A Deeper Look

The Lucas critique, introduced in Chapter 1 and formalized in Chapter 6, deserves a deeper treatment here in light of the rational expectations framework.

The key insight is that any reduced-form empirical relationship between policy instruments and economic outcomes is valid only conditional on the policy rule in effect when the data were generated. When the rule changes, agents update their optimal behavior — their consumption rules, wage-setting rules, investment rules — and the reduced-form relationship shifts accordingly.

Consider the consumption function $C = a + b(Y - T)$. The MPC $b$ is not a structural parameter: it reflects the optimal response of households to income changes *given the policy environment they expect*. If households expect temporary tax cuts to be followed by tax increases (a Ricardian equivalence expectation), $b$ will be small. If households are credit-constrained and expect no such reversal, $b$ will be large. When the policy environment changes — say, the government shifts from temporary to permanent tax cuts — $b$ changes too, and any prediction based on the historically estimated $b$ will be wrong.

This is not merely a theoretical curiosity. Attempts in the 1970s to use estimated Phillips curve trade-offs to design inflationary policies that would permanently reduce unemployment failed precisely because the estimated trade-offs reflected the low-inflation regime of the 1960s and collapsed once the regime changed. The methodological lesson is that structural models — built from preference and technology parameters that do not depend on the policy regime — are required for policy analysis. The DSGE approach is the field's response to this requirement.

---

*Next: Chapter 17 — The Goods Market*
