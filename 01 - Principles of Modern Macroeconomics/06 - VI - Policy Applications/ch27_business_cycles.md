# Chapter 27 — Business Cycles: Booms, Recessions, and Stabilization

> *"The curious task of economics is to demonstrate to men how little they really know about what they imagine they can design."*
> — Friedrich Hayek, *The Fatal Conceit*, 1988

---

Business cycles — the recurrent expansions and contractions in aggregate economic activity — are among the most studied phenomena in macroeconomics. Their study combines three distinct projects: the empirical project of characterizing their statistical regularities, the theoretical project of explaining their sources and propagation, and the normative project of evaluating whether stabilization policy can make them shorter and less severe. Each project has generated significant controversy. The debates about business cycles — whether they are efficient responses to real shocks or costly departures from potential output requiring policy correction — are among the most intellectually consequential in economics, and they are far from settled.

---

## 27.1 Empirical Regularities: The Stylized Facts

Before any theoretical framework can be evaluated, the phenomenon to be explained must be characterized precisely. The standard approach, pioneered by Kydland and Prescott (1982), is to extract cyclical components using the Hodrick-Prescott filter [P:Ch.6], compute second moments — standard deviations, correlations with output, and autocorrelations — and evaluate models by comparing their implied moments to the data.

Using U.S. quarterly data from 1960Q1 to 2019Q4, the key stylized facts are:

| Variable | $\sigma_x/\sigma_y$ | $\mathrm{corr}(x,y)$ | First autocorr. |
|---|---|---|---|
| Output $Y$ | 1.00 | 1.00 | 0.87 |
| Consumption $C$ | 0.60 | 0.84 | 0.85 |
| Investment $I$ | 3.67 | 0.92 | 0.82 |
| Hours worked $N$ | 0.98 | 0.87 | 0.83 |
| Real wage $w$ | 0.45 | 0.34 | 0.74 |
| Productivity $Y/N$ | 0.58 | 0.55 | 0.68 |
| Price level $P$ | 0.58 | $-0.12$ | 0.96 |
| Money $M1$ | 0.70 | 0.32 | 0.87 |

Several patterns stand out. Investment is nearly four times more volatile than output — making it the dominant source of output fluctuations in an expenditure decomposition and the hardest aggregate for any model to match. Consumption is substantially smoother than output, consistent with the consumption smoothing theories of Chapter 11. The real wage is only weakly procyclical (correlation of 0.34), not strongly procyclical as competitive frictionless labor markets predict — a diagnostic that points toward demand-driven fluctuations and nominal rigidities. The price level is only weakly countercyclical, suggesting that supply shocks (strongly countercyclical prices) and demand shocks (procyclical prices) coexist and partially cancel. Money leads output, but the lead-lag relationship does not establish causation — it is consistent with money accommodating demand or with monetary easing preceding recovery.

### International Evidence

The stylized facts above are broadly robust across developed economies, though with important cross-country variation. Business cycles in small open economies (Australia, Canada, Sweden) are more correlated with global cycles than in large economies, reflecting trade and financial linkages. The labor market stylized facts differ: European economies exhibit less employment volatility and more real wage volatility than the U.S., consistent with stronger employment protection legislation reducing the extensive margin of adjustment and shifting fluctuations onto the wage margin. These cross-country differences are important for evaluating whether a single model can explain business cycles universally or whether institutional variation requires country-specific models.

---

## 27.2 The Real Business Cycle Model

The RBC research program, launched by Kydland and Prescott (1982) and Long and Plosser (1983), attempted to explain business cycles as the optimal response of households and firms to exogenous fluctuations in total factor productivity. The central claim: business cycle fluctuations are not inefficient deviations from potential requiring policy correction — they are the competitive equilibrium responses to real shocks.

### Model Structure

The representative household maximizes expected lifetime utility:

$$U = \mathbb{E}_0\sum_{t=0}^\infty\beta^t\!\left[\frac{c_t^{1-\sigma}-1}{1-\sigma} + \chi\frac{(1-n_t)^{1-\eta}-1}{1-\eta}\right],$$

where $c_t$ is consumption, $1-n_t$ is leisure, $\sigma$ is the coefficient of relative risk aversion, and $\eta$ governs the curvature of the disutility of labor. The Frisch elasticity of labor supply is $\varepsilon_F = (1-n)/(\eta n)$ [Ch. 13].

Production uses capital and labor under stochastic TFP $A_t$:

$$Y_t = A_t K_t^\alpha N_t^{1-\alpha}.$$

TFP follows a stationary AR(1) in logs:

$$\ln A_t = \rho_A\ln A_{t-1} + \epsilon_t^A, \quad \epsilon_t^A \sim \mathcal{N}(0, \sigma_A^2).$$

Capital accumulation: $K_{t+1} = (1-\delta)K_t + I_t$. All markets clear every period at fully flexible prices; the competitive equilibrium is Pareto efficient.

### Solution and Calibration

The model is solved by log-linearizing around the non-stochastic steady state, yielding a set of linear expectational difference equations solved by the method of undetermined coefficients. The policy functions take the form $\hat{x}_t = \psi_x^k \hat{k}_t + \psi_x^A \hat{A}_t$ for each aggregate $x$.

Standard calibration: $\alpha = 1/3$, $\beta = 0.99$ (quarterly), $\delta = 0.025$, $\sigma = 1$ (log utility), $\varepsilon_F \approx 1.3$, $\rho_A = 0.95$, $\sigma_A = 0.72\%$ (chosen to match output volatility).

The calibrated model generates the following second moments, compared to data:

| Variable | Data $\sigma_x/\sigma_y$ | RBC model $\sigma_x/\sigma_y$ |
|---|---|---|
| Consumption | 0.60 | 0.49 |
| Investment | 3.67 | 3.16 |
| Hours | 0.98 | 0.93 |
| Real wage | 0.45 | 0.55 |
| Productivity | 0.58 | 0.73 |

The match is remarkable for a single-shock model: it correctly generates more volatile investment than consumption, and roughly correct volatility ratios for hours and real wages. This was striking enough to launch a major research program.

### Criticisms and Controversies

Despite its successes, the RBC model faces deep empirical and theoretical challenges that have prevented it from achieving consensus status.

**The TFP shock puzzle.** The required TFP volatility ($\sigma_A = 0.72\%$ quarterly) is large. Solow residuals — the empirical counterpart to $\hat{A}_t$ — are contaminated by variable factor utilization (workers and machines work harder in booms, making measured productivity procyclical even with constant true TFP), measurement error, and non-constant returns to scale. Correcting for utilization (Basu, Fernald, and Kimball, 2006) reduces purified TFP volatility substantially, leaving less for the model to explain.

**The labor supply elasticity puzzle.** The model requires a Frisch elasticity of approximately 1.5–2 to match employment volatility. Microeconometric estimates for prime-age males yield $\varepsilon_F \approx 0.1$–$0.3$ [Ch. 13]. This order-of-magnitude discrepancy has no clean resolution: the indivisible labor extension (Hansen, 1985) raises the aggregate elasticity by shifting adjustment to the extensive margin, but even then calibrated extensive-margin elasticities are below what the model requires.

**Monetary non-neutrality.** The RBC model is silent on the real effects of monetary policy — by construction, money is neutral. This is not merely a modelling gap but an empirical failure: the Romer-Romer (1989) narrative evidence, the Bernanke-Blinder (1992) VAR evidence, and the natural experiment evidence from Nakamura and Steinsson (2018) all document significant real effects of monetary contractions. A model that cannot account for the Fed's largest observed effects on the economy has a limited claim to be a general theory of business cycles.

**The role of recessions.** The RBC model implies that recessions are Pareto-efficient — the social planner would choose the same allocation as the competitive equilibrium. This implies no role for stabilization policy: recessions should be embraced as the efficient response to negative technology shocks. This normative implication is contested by the robust evidence that: (i) recessions have large and persistent effects on individual workers through unemployment scarring [Ch. 31]; (ii) demand management in the 1930s could have substantially shortened the Great Depression [Ch. 40]; and (iii) monetary policy has documented real effects.

---

## 27.3 The New Keynesian Model of Business Cycles

The New Keynesian (NK) model retains the RBC framework's commitment to microfoundations and rational expectations while adding two modifications: **monopolistic competition** in goods and labor markets, and **nominal price rigidities** via Calvo pricing. These modifications restore short-run monetary non-neutrality without abandoning the DSGE methodology.

### The Three-Equation System

The log-linearized NK model [Chs. 7, 10, 23]:

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\pi_{t+1}] - r_t^n) \quad \text{(NK-IS)}$$
$$\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t + u_t \quad \text{(NKPC)}$$
$$i_t = r^n + \phi_\pi\pi_t + \phi_x\hat{x}_t + \epsilon_t^m \quad \text{(Taylor rule)},$$

where $r_t^n$ is the natural rate of interest driven by technology and preference shocks, $u_t$ is a cost-push shock capturing markup fluctuations, and $\epsilon_t^m$ is a monetary policy shock.

### Impulse Responses

The model generates qualitatively different impulse responses for different shocks, matching the sign patterns in the data.

A positive **demand shock** ($\epsilon_t^m < 0$, a monetary easing) raises both output and inflation — the central bank can temporarily move the economy up the NKPC. A positive **supply shock** (a rise in $r_t^n$ from improved TFP) raises output and reduces inflation — the divine coincidence [Ch. 23]. A positive **cost-push shock** ($u_t > 0$, a markup increase or supply disruption) raises inflation but reduces output — stagflation — creating the policy trade-off that was so prominent in the 1970s.

### Estimated NK Models

Christiano, Eichenbaum, and Evans (2005) augment the baseline NK model with capital accumulation, investment adjustment costs, variable capital utilization, habit formation in consumption, wage stickiness, and indexation — features motivated by the inability of the baseline model to generate sufficient internal propagation. Smets and Wouters (2007) estimate this medium-scale model on U.S. data using Bayesian methods [M:Ch.30] with seven structural shocks and find:

- Technology shocks explain approximately 20–30% of output variance at business cycle frequencies.
- Investment-specific technology shocks explain another 20–30%.
- Monetary policy and markup shocks together explain approximately 30–40%.

This decomposition straddles the RBC view (real shocks dominate) and the Old Keynesian view (demand management is central). Neither camp fully wins; both types of shocks matter at empirically relevant magnitudes.

---

## 27.4 Identification of Business Cycle Shocks

A central methodological challenge in business cycle research is identifying structural shocks from the reduced-form comovement of macroeconomic variables. Several strategies have been employed.

### Narrative Identification

Romer and Romer (1989, 2004) construct monetary policy shock series from the Fed's own Greenbook forecasts and FOMC meeting minutes, identifying episodes in which the Fed deliberately tightened to reduce inflation — as opposed to responding to an already-deteriorating economy. This narrative approach cleanly separates policy-induced contractions from endogenous policy responses. Their estimated response of output to a one-percentage-point monetary tightening: output falls by approximately 3% over 18 months, with unemployment rising by approximately 1.5–2 percentage points at the peak.

Blanchard and Quah (1989) identify aggregate demand and supply shocks in a bivariate VAR (output, unemployment) using the restriction that supply shocks have permanent effects on output while demand shocks have only temporary effects. Their finding: supply shocks account for roughly 70% of output variance in the long run, while demand shocks account for the majority of short-run and unemployment fluctuations.

### Sign Restrictions

Uhlig (2005) identifies monetary policy shocks by imposing that a contractionary shock should raise interest rates and not simultaneously raise output or prices for several quarters — restrictions consistent with any reasonable theory. His identified shock produces output declines consistent with the narrative approach, though his point estimate is somewhat smaller. The sign restriction approach avoids the strong a priori structure of Cholesky identification while still leveraging theoretical discipline.

### Cross-Sectional and Quasi-Experimental Evidence

Nakamura and Steinsson (2018) exploit the differential response of regions to monetary shocks — exploiting variation in the exchange rate exposure of industries across U.S. states — to estimate monetary non-neutrality in a setting where monetary policy is genuinely exogenous to local conditions. Their estimates of the real effect of monetary policy are substantially larger than those from aggregate time-series identification, suggesting that attenuation bias from endogeneity is significant in aggregate VARs.

---

## 27.5 The Great Moderation and What Ended It

From the mid-1980s to 2007, U.S. and global business cycle volatility fell sharply — the period known as the **Great Moderation**. U.S. output volatility fell by roughly half; unemployment fluctuations shrank; recessions became less frequent and less deep. Three explanations attracted support: improved monetary policy (the Fed's adherence to the Taylor principle after 1979 eliminated the self-reinforcing inflation dynamics that amplified the 1970s recessions); reduced shock volatility (commodity prices became less volatile, and inventory management improved through just-in-time supply chains); and structural change (the shift from manufacturing to services reduced the procyclical investment-intensive sector's share of GDP).

The 2007–09 financial crisis ended the Great Moderation abruptly and revealed several ways the preceding calm had been illusory. The DSGE models standard at the time lacked financial sectors capable of generating systemic crises; the housing market dynamics that created the crisis were treated as outside the model; and the ELB, treated as a theoretical curiosity by most practitioners, bound severely and for years. Chapter 40 examines the Great Recession as a detailed case study; the methodological lessons — the need for financial frictions, heterogeneous agents, and non-Gaussian tail risks — have reshaped the research frontier described in Chapter 39.

---

## Chapter Summary

- Business cycle **stylized facts** — investment four times more volatile than output, consumption smoother than output, weakly procyclical real wages, weakly countercyclical prices — provide the empirical target that all models must match.

- The **RBC model** explains cycles as efficient responses to TFP shocks. It matches investment and consumption volatility but requires implausibly large Frisch elasticities and TFP shocks, is silent on monetary non-neutrality, and implies that recessions are socially optimal.

- The **New Keynesian model** adds monopolistic competition and Calvo pricing, restoring monetary non-neutrality. Estimated medium-scale NK models attribute roughly equal shares of output variance to real shocks (technology, investment-specific) and nominal/financial shocks (monetary policy, markups).

- **Structural identification** via narrative methods (Romer-Romer, Blanchard-Quah), sign restrictions (Uhlig), and quasi-experimental variation (Nakamura-Steinsson) converges on significant real effects of monetary policy and a dominant role for demand shocks in unemployment fluctuations.

- The **Great Moderation** reflected better policy and reduced shocks; its abrupt end in 2008 revealed the limitations of models without financial frictions, fat tails, and endogenous risk.

---

*Next: Chapter 28 — Fiscal Policy in Practice: Government Spending and Taxation*
