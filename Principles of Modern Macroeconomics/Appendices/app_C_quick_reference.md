# Appendix C — Study Quick Reference Guide

This appendix provides compact summaries of the book's core models, equations, and results for rapid review. It is organized by topic and is intended to be used alongside — not instead of — the main text.

---

## C.1 The Core Equations at a Glance

### National Accounts Identities
$$Y \equiv C + I + G + NX \quad \text{(expenditure approach)}$$
$$(S-I) + (T-G) = NX \quad \text{(saving-investment identity)}$$
$$CA + FA + \Delta R = 0 \quad \text{(balance of payments identity)}$$

### The Solow Model (Cobb–Douglas, $f(\tilde{k}) = \tilde{k}^\alpha$)
$$\dot{\tilde{k}}_t = s\tilde{k}_t^\alpha - (n+g+\delta)\tilde{k}_t$$
$$\tilde{k}^* = \left(\frac{s}{n+g+\delta}\right)^{1/(1-\alpha)}, \quad \tilde{y}^* = \left(\frac{s}{n+g+\delta}\right)^{\alpha/(1-\alpha)}$$
$$\text{Convergence rate: } \lambda = (1-\alpha)(n+g+\delta)$$
$$\text{Golden Rule: } f'(\tilde{k}^{GR}) = n+g+\delta$$

### Ramsey–Cass–Koopmans Model
$$\text{Euler equation: } \frac{\dot{c}_t}{c_t} = \frac{r_t - \rho}{\sigma}$$
$$\text{Steady state: } f'(\tilde{k}^*) = \delta + \rho + \sigma g$$
$$\text{Transversality: } \lim_{t\to\infty} e^{-(\rho-n)t}\mu_t\tilde{k}_t = 0$$

### IS–LM Model (Linear)
$$\text{IS: } Y = \bar{A} - b_r r \quad (\text{with multiplier absorbed in }\bar{A})$$
$$\text{LM: } M/P = kY - hi$$
$$\text{Fiscal multiplier: } \frac{\Delta Y}{\Delta G} = \frac{h}{h+b_r k}$$
$$\text{Monetary multiplier: } \frac{\Delta Y}{\Delta(M/P)} = \frac{b_r}{h+b_r k}$$
$$\text{Keynesian cross multiplier: } \kappa_G = \frac{1}{1-b}, \quad \kappa_T = \frac{-b}{1-b}$$

### AS–AD and Phillips Curves
$$\text{Lucas supply: } Y_t = \bar{Y}_t + \alpha(P_t - P_t^e)$$
$$\text{EAPC: } \pi_t = \pi_t^e - \alpha(u_t - u^*) + \epsilon_t$$
$$\text{Accelerationist: } \pi_t - \pi_{t-1} = -\alpha(u_t - u^*) + \epsilon_t$$
$$\text{NKPC: } \hat{\pi}_t = \beta\,\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\,\hat{x}_t, \quad \kappa = \frac{(1-\theta)(1-\beta\theta)}{\theta}$$
$$\text{NKPC (forward solution): } \hat{\pi}_t = \kappa\sum_{k=0}^\infty \beta^k\,\mathbb{E}_t[\hat{x}_{t+k}]$$

### New Keynesian Three-Equation Model
$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\pi_{t+1}] - r_t^n) \quad \text{(NK IS)}$$
$$\hat{\pi}_t = \beta\,\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\,\hat{x}_t \quad \text{(NKPC)}$$
$$i_t = r^n + \pi^* + \phi_\pi(\pi_t-\pi^*) + \phi_y\hat{x}_t \quad \text{(Taylor rule)}$$

### Investment
$$\text{Neoclassical: } F_K(K,L) = r+\delta \equiv c^K$$
$$\text{q model: } I_t/K_t = (q_t-1)/\psi$$
$$\text{Real option trigger: } \Pi^* = \frac{\beta}{\beta-1}c^K, \quad \beta = \frac{1}{2}-\frac{\mu}{\sigma^2}+\sqrt{\left(\frac{\mu}{\sigma^2}-\frac{1}{2}\right)^2+\frac{2r}{\sigma^2}}$$

### Consumption
$$\text{Euler equation: } u'(c_t) = \beta(1+r)\,\mathbb{E}_t[u'(c_{t+1})]$$
$$\text{Log Euler (CRRA): } \mathbb{E}_t[\Delta\ln c_{t+1}] = (r-\rho)/\sigma + \tfrac{\sigma}{2}\mathrm{Var}_t[\Delta\ln c_{t+1}]$$
$$\text{Life-cycle: } c = (A_0 + T_W y)/T, \quad \text{MPC} = T_W/T$$

### Labor Market
$$\text{Natural rate: } u^* = \delta/(\delta+f(\theta^*))$$
$$\text{Efficiency wage NSC: } w \geq b + e_H(r+q_f)/q_f\cdot[1 + r/(u/(1-u))]$$
$$\text{Okun's Law: } Y_t - \bar{Y}_t = -\psi(u_t - u^*)$$

### Money
$$\text{Fisher equation: } MV = PY \implies \hat{m}+\hat{v} = \pi+\hat{y}$$
$$\text{Baumol-Tobin: } L^{BT} = \sqrt{bY/(2i)}, \quad \varepsilon_Y = 1/2, \quad \varepsilon_i = -1/2$$
$$\text{Money multiplier: } M1 = m\cdot H, \quad m = (1+c_r)/(c_r+rr+er)$$
$$\text{Seigniorage: } S = \pi\cdot(M/P)$$

### Government Budget and Debt
$$\text{Flow: } \dot{b}_t = (r_t-g_t)b_t - s_t$$
$$\text{IGBC: } b_0 = \int_0^\infty e^{-\int_0^t(r_s-g_s)\mathrm{d}s}s_t\,\mathrm{d}t$$
$$\text{Sustainability: } s_t > (r_t-g_t)b_t$$

### Barro–Gordon Inflationary Bias
$$\text{Loss: } L = \tfrac{1}{2}\pi^2 + \tfrac{\lambda}{2}(y-y^*)^2$$
$$\text{Discretionary equilibrium: } \pi^D = b\lambda(y^*-\bar{y})$$

### Asset Pricing
$$\text{SDF pricing: } p_t = \mathbb{E}_t[M_{t+1}(p_{t+1}+d_{t+1})]$$
$$\text{CCAPM SDF: } M_{t+1} = \beta(c_{t+1}/c_t)^{-\sigma}$$
$$\text{Risk premium: } \mathbb{E}_t[R^j]-R^f = -\mathrm{Cov}_t(M_{t+1},R^j)/\mathbb{E}_t[M_{t+1}]$$
$$\text{Gordon growth: } P^{eq} = D/(i^e-g)$$

### Exchange Rates
$$\text{UIP: } i_t = i_t^* + \mathbb{E}_t[\hat{e}_{t+1}]$$
$$\text{CIP: } i_t - i_t^* = f_t - e_t$$
$$\text{Relative PPP: } \hat{e}_t = \pi_t - \pi_t^*$$

---

## C.2 Key Parameters and Typical Calibrated Values

| Parameter | Symbol | Typical Value | Source |
|-----------|--------|---------------|--------|
| Capital share | $\alpha$ | 0.33 | National accounts |
| Depreciation rate (quarterly) | $\delta$ | 0.025 | Investment data |
| Population growth (annual) | $n$ | 0.01–0.02 | Demographics |
| Technology growth (annual) | $g$ | 0.015–0.02 | TFP estimates |
| Discount factor (quarterly) | $\beta$ | 0.99 | Real interest rate target |
| Risk aversion / EIS$^{-1}$ | $\sigma$ | 1–2 | Micro estimates |
| Price stickiness (Calvo) | $\theta$ | 0.75 | Frequency of price changes |
| NKPC slope | $\kappa$ | 0.1–0.2 | Estimated from inflation data |
| Taylor rule inflation | $\phi_\pi$ | 1.5 | Taylor (1993) |
| Taylor rule output | $\phi_y$ | 0.5 | Taylor (1993) |
| Interest rate smoothing | $\rho_i$ | 0.85 | Estimated |
| MPC (Keynesian cross) | $b$ | 0.75 | Consumption surveys |
| Technology shock persistence | $\rho_A$ | 0.95 | Solow residual |
| Technology shock std. dev. | $\sigma_A$ | 0.007 | Solow residual |
| Natural unemployment rate | $u^*$ | 0.04–0.05 | Structural estimates |
| Okun coefficient | $\psi$ | 2.0 | OLS regression |
| Sacrifice ratio (US) | $SR$ | 1.4–2.8 | Ball (1994) |

---

## C.3 Model Comparison Table

| Feature | Keynesian Cross | IS–LM | AS–AD (static) | NK 3-Eq. | RCK |
|---------|:-:|:-:|:-:|:-:|:-:|
| Price level | Fixed | Fixed | Endogenous | Inflation rate | Endogenous |
| Interest rate | Exogenous | Endogenous | IS-LM equilibrium | Taylor rule | $r = f'(k)-\delta$ |
| Investment | Exogenous | Endogenous | IS curve | NK IS | Euler equation |
| Expectations | Static | Static | Exogenous $P^e$ | Rational | Rational |
| Microfoundations | No | No | No | Partial | Yes |
| Long-run | Not modeled | Not modeled | LRAS vertical | Not modeled | Balanced growth |
| Policy analysis | Fiscal only | Fiscal + monetary | Both | Both + welfare | Welfare-optimal |

---

## C.4 Key Empirical Facts to Remember

**Business cycles (US, quarterly, HP-filtered):**
- $\sigma_y \approx 1.5\%$, $\sigma_c \approx 0.9\%$, $\sigma_i \approx 5.5\%$
- $\text{corr}(y,c) \approx 0.84$, $\text{corr}(y,i) \approx 0.92$, $\text{corr}(y,u) \approx -0.85$

**Inflation:**
- Boskin Commission bias: ~1.1 pp/year overstatement in CPI
- CPI income elasticity: ~0.5 (Baumol–Tobin prediction)
- Long-run: inflation ≈ money growth (quantity theory)

**Growth:**
- US long-run growth: ~2% per year (per capita, real)
- Convergence speed $\lambda \approx 0.02$–$0.04$ per year
- Capital share $\alpha \approx 1/3$ in most countries

**Monetary policy:**
- Pass-through from rate change to output: ~0.5% GDP per 1 pp rate increase (peak at 1–2 years)
- QE: ~90 bps reduction in 10-year yields per $1.75T in Fed purchases (Gagnon et al., 2011)

**Fiscal policy:**
- Spending multiplier: 0.6–1.5 (depends on regime/cycle)
- Tax multiplier: −2 to −3 (larger in absolute value than spending multiplier per dollar of revenue impact)
- ELB multiplier: potentially >1.5 (Christiano, Eichenbaum, Rebelo, 2011)

**Labor markets:**
- US Phillips curve slope $\alpha$: ~0.3–0.5
- Frisch elasticity: ~0.1–0.3 (micro estimates) vs. ~2 (macro requirement)
- Sacrifice ratio: ~1.4–2.8 for US disinflation episodes

---

## C.5 Reading a Macroeconomics Paper: A Checklist

When reading a primary research paper in macroeconomics, the following questions organize the critical evaluation:

**Model/theory papers:**
1. What is the question? What economic phenomenon is the model designed to explain?
2. What are the key assumptions? Which are standard (shared with the literature) and which are novel?
3. What is the equilibrium concept? Are prices flexible or sticky? Are expectations rational?
4. What is the main result? Is it existence, uniqueness, comparative statics, or quantitative?
5. What assumptions drive the result? Would it survive relaxing them?
6. What does the model predict that can be tested empirically?

**Empirical papers:**
1. What is the question? What causal effect is being estimated?
2. What is the identification strategy? What is the source of exogenous variation?
3. What are the identifying assumptions? Are they plausible? How are they tested?
4. What is the sample? Does it have external validity beyond the time period and countries studied?
5. How large and how precisely estimated are the effects? Are they economically significant?
6. What alternative explanations are considered and how are they ruled out?

**Structural estimation papers:**
1. What model is being estimated? What are its key equations and calibrated/estimated parameters?
2. Which parameters are estimated jointly, and which are calibrated from outside the model?
3. How is the model solved? Linear approximation? Global methods?
4. How is the likelihood evaluated? Kalman filter? Particle filter?
5. What prior distributions are used? Are they informative? Defensible?
6. How well does the estimated model fit the data? How does it compare to reduced-form benchmarks?

---

## C.6 Suggested Reading Paths by Topic

**Monetary economics:** Chapters 14, 18, 23, 29; Appendix B (Sections B.3, B.5)

**Growth theory:** Chapters 5, 25, 33; Appendix D (Section D.4)

**Business cycles:** Chapters 7, 10, 27; Appendix B (Sections B.2, B.3)

**Fiscal policy:** Chapters 8, 22, 28; Appendix B (Section B.4)

**Open economy:** Chapters 21, 26, 32, 35; Appendix J

**Financial crises:** Chapters 20, 24, 34, 40; Appendix I

**Labor markets:** Chapters 13, 19, 31; Appendix D (Section D.2)

**Inequality and distribution:** Chapters 25, 38; Appendix I

---

*For further guidance on specific topics, see Appendix H (Bibliography and Further Reading).*
