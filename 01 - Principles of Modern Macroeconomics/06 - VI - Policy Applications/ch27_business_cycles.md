# Chapter 27 — Business Cycles: Booms, Recessions, and Stabilization

> *"The curious task of economics is to demonstrate to men how little they really know about what they imagine they can design."*
> — Friedrich Hayek, *The Fatal Conceit*, 1988

---

Business cycles — the recurrent expansions and contractions in economic activity that affect employment, investment, output, and welfare — are among the most studied phenomena in macroeconomics. Their study combines the empirical project of characterizing their statistical regularities with the theoretical project of explaining their sources and propagation mechanisms, and the normative project of evaluating whether stabilization policy can make them shorter and less severe. Each of these projects has generated significant controversy, and the debates between schools of thought about the nature of business cycles — whether they are efficient responses to real shocks or costly departures from potential output requiring policy correction — are among the most intellectually important in economics.

---

## 27.1 Empirical Regularities: The Stylized Facts

Before any theoretical framework can be evaluated, the phenomenon to be explained must be characterized precisely. The standard approach is to compute HP-filtered log-deviations from trend for each aggregate variable and compute their second moments — standard deviations, correlations with output, and autocorrelations.

Using U.S. quarterly data (1960Q1–2019Q4), the key stylized facts are:

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

Several patterns stand out. Investment is nearly four times more volatile than output, making it the dominant source of output fluctuations in an expenditure decomposition. Consumption is substantially smoother than output — consistent with consumption smoothing theories. The real wage is weakly procyclical, not strongly procyclical as competitive labor market models predict. The price level is only weakly countercyclical, suggesting that supply shocks (which would make prices strongly countercyclical) and demand shocks (which make prices procyclical) coexist. Money leads output, but the lead-lag relationship does not establish that money causes output.

---

## 27.2 The Real Business Cycle Model

The RBC research program, launched by Kydland and Prescott (1982) and Long and Plosser (1983), attempted to explain business cycles as the optimal response of households and firms to exogenous fluctuations in total factor productivity. The central proposition: business cycle fluctuations are not deviations from potential that policy should offset — they are the efficient equilibrium responses to real shocks.

### Model Structure

The representative household maximizes expected lifetime utility:

$$U = \mathbb{E}_0\sum_{t=0}^\infty\beta^t\!\left[\frac{c_t^{1-\sigma}-1}{1-\sigma} + \chi\frac{(1-n_t)^{1-\eta}-1}{1-\eta}\right],$$

where $c_t$ is consumption, $1-n_t$ is leisure (hours worked is $n_t$), $\sigma$ is risk aversion, $\eta$ governs labor supply elasticity, and $\chi$ scales the utility of leisure. The Frisch elasticity of labor supply is $\varepsilon_F = (1-n)/(\eta n)$.

The production technology with stochastic TFP $A_t$:

$$Y_t = A_t K_t^\alpha N_t^{1-\alpha}.$$

TFP follows a stationary AR(1) in logs:

$$\ln A_t = \rho_A\ln A_{t-1} + \epsilon_t^A, \quad \epsilon_t^A \sim \mathcal{N}(0, \sigma_A^2).$$

Capital accumulation: $K_{t+1} = (1-\delta)K_t + I_t$. All markets clear every period; prices are fully flexible; wages and the rental rate on capital clear labor and capital markets.

### Solution by Log-Linearization

The model is solved by log-linearizing around the deterministic steady state. Define $\hat{x}_t = \ln x_t - \ln x^*$ for each variable. The log-linearized system reduces to a set of linear expectational difference equations:

$$\mathbf{E}_t[\hat{\mathbf{x}}_{t+1}] = \Gamma\hat{\mathbf{x}}_t + \Omega\hat{\epsilon}_{t+1},$$

solved by the method of undetermined coefficients or by finding the stable eigendecomposition of $\Gamma$. The unique bounded solution gives the **policy functions**:

$$\hat{c}_t = \psi_c^k\hat{k}_t + \psi_c^A\hat{A}_t, \quad \hat{n}_t = \psi_n^k\hat{k}_t + \psi_n^A\hat{A}_t, \quad \hat{i}_t = \psi_i^k\hat{k}_t + \psi_i^A\hat{A}_t.$$

### Calibration and Results

Standard calibration: $\alpha = 1/3$, $\beta = 0.99$ (quarterly), $\delta = 0.025$, $\sigma = 1$ (log utility), $\eta = 2$ (Frisch elasticity $\approx 1.3$), $\rho_A = 0.95$, $\sigma_A = 0.0072$.

The calibrated model generates the following second moments:

| Variable | Data $\sigma_x/\sigma_y$ | RBC model $\sigma_x/\sigma_y$ |
|---|---|---|
| Consumption | 0.60 | 0.49 |
| Investment | 3.67 | 3.16 |
| Hours | 0.98 | 0.93 |
| Real wage | 0.45 | 0.55 |
| Productivity | 0.58 | 0.73 |

The RBC model matches the data reasonably well for the volatility of consumption and investment relative to output, and for the high correlation of hours with output. This was striking: a model with a single shock (technology) and no nominal rigidities could replicate key business cycle facts.

### Criticisms and Controversies

The RBC model faces several serious empirical and theoretical challenges.

**The TFP shock puzzle.** The required TFP volatility ($\sigma_A = 0.72\%$ quarterly) is large. Solow residuals — the empirical counterpart to $\hat{A}_t$ — are contaminated by variable factor utilization, measurement error, increasing returns, and other non-TFP factors. Correcting for utilization (Basu, Fernald, and Kimball, 2006) substantially reduces measured TFP volatility, leaving a smaller residual for the model to explain.

**The labor supply elasticity puzzle.** The model requires a large Frisch elasticity ($\varepsilon_F \approx 1.5$–$2$) to match employment volatility. Microeconometric estimates of the Frisch elasticity for prime-age males cluster at 0.1–0.3. The two-order-of-magnitude difference between the macro requirement and the micro estimate is the central unresolved tension in the RBC literature.

**The absence of monetary non-neutrality.** The RBC model has no role for monetary policy — money is neutral by construction. This was an intended feature (the program was partly motivated by a desire to show that business cycles could be explained without monetary frictions), but it renders the model silent on some of the most prominent empirical regularities in macroeconomics: the impact of Fed tightenings on unemployment, the relationship between money growth and inflation, and the effects of the zero lower bound.

---

## 27.3 The New Keynesian Model of Business Cycles

The New Keynesian (NK) model retains the RBC framework's commitment to microfounded, rational-expectations dynamics while adding two crucial modifications: **monopolistic competition** in goods markets and **nominal rigidities** via Calvo pricing. These modifications restore short-run monetary non-neutrality without abandoning the DSGE methodology.

The log-linearized NK model (the three-equation system from Chapters 7, 10, and 23) takes the form:

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\pi_{t+1}] - r_t^n) \quad \text{(DIS)}$$
$$\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t + u_t \quad \text{(NKPC)}$$
$$i_t = r^n + \phi_\pi\pi_t + \phi_x\hat{x}_t + \epsilon_t^m \quad \text{(Taylor rule)},$$

where $r_t^n$ is the natural rate (driven by technology and preference shocks) and $u_t$ is a cost-push shock. The model generates business cycles from three types of shocks:

**Demand shocks** ($\epsilon_t^m$, monetary policy surprises; preference shocks entering the DIS): shift AD without shifting the NKPC. Output and inflation move in the same direction; the output gap and inflation are both stable under optimal policy.

**Supply shocks** ($r_t^n$ movements from technology; cost-push shocks $u_t$): shift AS. Output and inflation move in opposite directions; there is a stabilization trade-off.

**Markup shocks**: fluctuations in the degree of monopoly power ($\kappa$ or $\mu$). These generate stagflationary episodes — high inflation and low output simultaneously.

The NK model fits U.S. business cycle data substantially better than the basic RBC model once financial shocks, investment adjustment costs, and variable capital utilization are added (Smets and Wouters, 2007). In an estimated medium-scale NK model, approximately 60% of output variance is explained by technology and preference shocks (the "real" component, consistent with the RBC view), and approximately 40% by monetary policy shocks, markup shocks, and financial shocks (the "nominal/financial" component, emphasizing the role of market frictions).

---

## 27.4 Identification of Business Cycle Shocks

The empirical study of business cycles requires identifying structural shocks — exogenous, economically interpretable disturbances — from the reduced-form comovement of macro variables. The methods of Appendix B apply here; we discuss the major identification approaches and their findings.

**Narrative identification.** Romer and Romer (1989, 2004) identify monetary policy shocks from the Fed's stated intentions (Greenbook forecasts and meeting minutes), constructing dates on which the Fed deliberately tightened to reduce inflation. These shock dates are used as instruments in VARs. The estimated response of output to a monetary tightening: output falls by approximately 3% over 1–2 years, with the peak effect at about 18 months. Unemployment rises by approximately 1.5–2 pp at the peak.

**Structural VAR with sign restrictions.** Uhlig (2005) identifies monetary policy shocks by imposing the sign restrictions that a contractionary shock (interest rate rises) should not increase output or prices for several quarters. His identified shock generates a decline in output and prices consistent with the narrative approach, but his central estimate of the output effect is somewhat smaller (approximately 1–1.5%).

**Long-run restrictions.** Blanchard and Quah (1989) identify aggregate demand shocks (no long-run effect on output) and supply shocks (permanent output effects) using long-run restrictions in a bivariate VAR with output and unemployment. Supply shocks account for roughly 70% of output variance in the long run; demand shocks are responsible for most short-run volatility and the bulk of unemployment fluctuations.

These identification exercises converge on a broadly consistent picture: both technology and demand shocks are quantitatively important; monetary policy tightening has significant short-run real effects consistent with NK nominal rigidities; and financial shocks play a large role in severe recessions (the great moderation of the 1980s–2000s was partly due to reduced shock volatility, partly to better monetary policy).

---

## 27.5 The Great Moderation and the Financial Crisis

From the mid-1980s to 2007, U.S. and global business cycle volatility fell sharply — a period called the **Great Moderation**. Output volatility halved; unemployment fluctuations became smaller; the frequency and severity of recessions declined. This coincided with improved monetary policy (better adherence to the Taylor principle), reduced shock volatility (less volatile commodity prices), and possible structural changes in the economy (better inventory management via JIT, financial deepening).

The 2007–09 financial crisis shattered the complacency that the Great Moderation had engendered. The crisis revealed: (i) DSGE models without financial sectors underestimated systemic risk; (ii) housing market dynamics and mortgage securitization created fragilities invisible to standard macro models; (iii) the zero lower bound could bind severely and persistently; (iv) hysteresis could cause temporary recessions to permanently raise the natural rate of unemployment. Part VIII examines these post-crisis developments in detail; the case study of the 2008 recession is in Chapter 40.

---

*Next: Chapter 28 — Fiscal Policy in Practice*
