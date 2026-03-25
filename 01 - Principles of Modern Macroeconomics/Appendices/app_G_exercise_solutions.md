# Appendix G — Exercise Solutions and Hints

This appendix provides worked solutions to selected exercises from each chapter, along with hints for the remaining exercises. Starred ($\star$) exercises are more advanced; double-starred ($\star\star$) require techniques beyond the main text and are intended for doctoral-level study.

---

## Part I: Foundations

### Chapter 1

**Exercise 1.1.** Define the distinction between positive and normative macroeconomics. Give one example of each from the policy debate about central bank independence.

*Solution.* A positive statement is a testable claim about facts or causal relationships. A normative statement expresses a value judgment. Positive example: "Countries with independent central banks have experienced lower average inflation over 1990–2020." This is testable against data. Normative example: "Central banks should be independent because price stability is more important than democratic accountability over monetary policy." This involves a value judgment about institutional design that cannot be resolved by evidence alone, even if evidence can inform it.

---

**Exercise 1.2.** Explain the Lucas critique using the Phillips curve as an example. Why does a historically estimated Phillips curve provide unreliable guidance for policy?

*Solution.* The historical Phillips curve $\pi_t = \pi_{t-1} - \alpha(u_t - u^*)$ is estimated under a policy regime in which the central bank sets money growth at some rule $m_t = \phi\mathbf{x}_{t-1} + \epsilon_t$. The estimated coefficient $\alpha$ reflects both the structural price-setting behavior of firms and workers' rational responses to the prevailing monetary regime. If the central bank changes its rule — say, adopting a stricter inflation target — agents' expectations formation changes, which changes how wages are set, which changes the estimated $\alpha$. A prediction of disinflation cost based on the historically estimated $\alpha$ will therefore be wrong. This is precisely what occurred in the 1970s: the estimated Phillips curves of the 1960s implied low-cost disinflation, but the Volcker disinflation was much more expensive because the regime shift required building credibility from scratch.

---

### Chapter 3

**Exercise 3.1.** The CPI uses a Laspeyres formula. Suppose the price of apples doubles and the price of oranges halves, while consumers shift their purchases entirely from apples to oranges. Show that the CPI overstates the true increase in the cost of living in this case.

*Solution.* Let the base-period basket be $q_0 = (q_A^0, q_O^0)$ with $q_A^0 > 0$ and $q_O^0 > 0$. Suppose $p_A^0 = p_O^0 = 1$. New prices: $p_A^1 = 2$, $p_O^1 = 0.5$. Base-period expenditure: $E_0 = q_A^0 + q_O^0$.

Laspeyres CPI: $\text{CPI}_1 = (2q_A^0 + 0.5q_O^0)/(q_A^0 + q_O^0) > 1$. The index rises above one.

True cost of living: at the new prices, the consumer buys only oranges. If she maintained the same utility level by consuming oranges, her expenditure could be $0.5 \cdot q_{O}^{true}$ where $q_O^{true}$ achieves the base-period utility. In many utility specifications, she is better off at new prices with a different basket. The Laspeyres index, by keeping the base-period basket fixed, assigns positive weight to the more expensive good (apples) even though the consumer has substituted away from it, thereby overstating the cost of living.

---

### Chapter 5

**Exercise 5.1.** In the Solow model with Cobb–Douglas production ($\alpha = 1/3$), $n = 0.01$, $g = 0.02$, $\delta = 0.05$, $s = 0.25$: (a) compute $\tilde{k}^*$ and $\tilde{y}^*$; (b) compute the half-life of convergence; (c) by what percentage does $\tilde{y}^*$ increase if $s$ rises to 0.30?

*Solution.*

**(a)** $\mu = n+g+\delta = 0.08$. From the Cobb–Douglas formula:
$$\tilde{k}^* = (s/\mu)^{1/(1-\alpha)} = (0.25/0.08)^{3/2} = (3.125)^{1.5} = 5.53$$
$$\tilde{y}^* = (s/\mu)^{\alpha/(1-\alpha)} = (3.125)^{0.5} = 1.77.$$

**(b)** Convergence rate: $\lambda = (1-\alpha)\mu = (2/3)(0.08) = 0.0533$. Half-life: $\ln 2/\lambda = 0.693/0.0533 \approx 13.0$ years.

**(c)** New $s' = 0.30$: $\tilde{y}^{*\prime} = (0.30/0.08)^{0.5} = (3.75)^{0.5} = 1.936$. Percentage change: $(1.936-1.77)/1.77 \approx 9.4\%$. A 20% increase in the saving rate raises steady-state income per effective worker by about 9.4%.

---

**Exercise 5.2** ($\star$). Show that in the AK model, the saving rate $s$ permanently affects the growth rate. Contrast with the Solow model, where $s$ affects only the level of income.

*Hint.* In the AK model, the capital accumulation equation is $\dot{K} = sAK - \delta K$, so $g_K = sA - \delta$. Since $Y = AK$, $g_Y = g_K = sA - \delta$: the growth rate of output equals $sA-\delta$, which is increasing in $s$. In the Solow model with $f(\tilde{k}) = \tilde{k}^\alpha$, $\alpha < 1$, the steady state is $s\tilde{k}^{*\alpha} = \mu\tilde{k}^*$, giving $\tilde{y}^* \propto s^{\alpha/(1-\alpha)}$. A higher $s$ raises the level of $\tilde{y}^*$ but the long-run growth rate remains $g$ regardless of $s$. The difference arises from diminishing returns: the Solow model has $\alpha < 1$, generating a stable steady state; the AK model has effective $\alpha = 1$, eliminating the steady state and sustaining growth.

---

## Part II: Core Theories

### Chapter 8

**Exercise 8.1.** In the Keynesian cross with $b = 0.8$ and $t = 0.25$ (proportional tax rate): (a) compute the spending multiplier; (b) compute the effect on equilibrium income of a simultaneous $\Delta G = 100$ and $\Delta T_{\text{lump-sum}} = 100$.

*Solution.*

**(a)** $\kappa_G^{prop} = 1/[1-b(1-t)] = 1/[1-0.8\times 0.75] = 1/0.4 = 2.5$.

**(b)** This is a balanced budget expansion. The spending multiplier (2.5) applies to $\Delta G = 100$, raising income by 250. The lump-sum tax multiplier is $\kappa_T = -b\kappa_G^{lump} = -0.8 \times (1/[1-b]) = -0.8\times 5 = -4$. So $\Delta Y_T = -4\times 100 = -400$. Net: $\Delta Y = 250 - 400 = -150$. Note: with proportional taxes already in place, adding a lump-sum tax causes a larger income reduction than the spending increase creates, so the net effect is contractionary despite a balanced budget. The balanced budget multiplier is exactly one only in the simplest lump-sum-tax, no-proportional-tax version.

---

### Chapter 9

**Exercise 9.1.** In the IS–LM model with $b_r = 2$, $k = 0.5$, $h = 4$, $b = 0.75$: (a) compute the IS–LM fiscal multiplier; (b) compute the IS–LM monetary multiplier; (c) compare to the Keynesian cross multiplier and explain the difference.

*Solution.*

**(a)** Fiscal: $\Delta Y/\Delta G = h/(h+b_r k) = 4/(4+2\times 0.5) = 4/5 = 0.8$.

**(b)** Monetary: $\Delta Y/\Delta(M/P) = b_r/(h+b_r k) = 2/5 = 0.4$.

**(c)** The Keynesian cross multiplier is $\kappa_G = 1/(1-0.75) = 4$. The IS–LM multiplier (0.8) is five times smaller. The difference is crowding out: the fiscal expansion raises income, raising money demand, raising interest rates, reducing investment. The amount of crowding out is $(1-0.8/4) = 80\%$ of the gross fiscal impulse is offset by reduced investment.

---

### Chapter 10

**Exercise 10.1.** Using the NKPC $\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t$ with $\beta = 0.99$, $\kappa = 0.15$: (a) If $\hat{x}_t = 0$ for all $t \geq 1$ and $\hat{x}_0 = 1$ (a one-period output gap boom), compute $\hat{\pi}_0$. (b) Compare to the prediction of the expectations-augmented PC with adaptive expectations and the same $\hat{x}_0 = 1$.

*Solution.*

**(a)** From the NKPC forward solution: $\hat{\pi}_0 = \kappa\sum_{k=0}^\infty \beta^k \mathbb{E}_0[\hat{x}_k] = \kappa\hat{x}_0 = 0.15\times 1 = 0.15$ (since $\hat{x}_k = 0$ for $k \geq 1$). A one-period boom generates exactly $\kappa$ units of inflation contemporaneously.

**(b)** The EAPC: $\pi_0 = \pi_{-1} - \alpha(u_0-u^*)$. Since $\hat{x}_0 > 0$ implies $u_0 < u^*$: $\pi_0 = \pi_{-1} + \alpha|\hat{x}_0|/\psi$ (using Okun's law). If $\alpha = 0.3$, $\psi = 2$: $\pi_0 - \pi_{-1} = 0.3/2 = 0.15$. The level of inflation in period 0 is higher than in period $-1$ by 0.15 pp. But the key difference from the NKPC: under adaptive expectations, this 0.15 pp increase in inflation persists into period 1 as $\pi_{-1}^{new} = \pi_{-1}+0.15$, even if $\hat{x}_1 = 0$. Under the NKPC, $\hat{\pi}_1 = \beta\mathbb{E}_1[\hat{\pi}_2] + \kappa\hat{x}_1 = 0$: inflation returns to zero immediately. The NKPC implies no inflation inertia; the EAPC with adaptive expectations implies full inertia.

---

## Part III: Microfoundations

### Chapter 11

**Exercise 11.1.** A household has CRRA utility with $\sigma = 2$, discount factor $\beta = 0.98$, and faces a constant real interest rate $r = 0.03$. (a) What is the optimal consumption growth rate? (b) If $r$ increases to $0.05$, by how much does consumption growth change? (c) Interpret the result in terms of the elasticity of intertemporal substitution.

*Solution.*

**(a)** Log Euler equation: $\mathbb{E}[\Delta\ln c] = (r-\rho)/\sigma$. With $\beta = 0.98$, $\rho = -\ln(0.98) \approx 0.0202$. $(0.03-0.0202)/2 \approx 0.0049$. Consumption grows at approximately 0.49% per period.

**(b)** New growth: $(0.05-0.0202)/2 = 0.0149$. Change: $0.0149-0.0049 = 0.01$, or 1 percentage point per year.

**(c)** The EIS is $1/\sigma = 0.5$. A 1 pp increase in the real interest rate raises consumption growth by 0.5 pp — the EIS gives the percentage change in consumption growth per percentage-point change in the real rate. Higher $\sigma$ (lower EIS) means households are more reluctant to substitute consumption intertemporally, so interest rate changes have smaller effects on the time path of consumption.

---

### Chapter 12

**Exercise 12.1** ($\star$). In the $q$ model with adjustment cost parameter $\psi = 10$ and steady-state $q^* = 1$: (a) What is the steady-state investment rate $I^*/K^*$? (b) If $q_0 = 1.5$ (above steady state), what is the initial investment rate? (c) What does this imply about the speed of capital adjustment?

*Solution.*

**(a)** Steady-state $q^* = 1$: $I^*/K^* = (q^*-1)/\psi = 0/10 = 0$. In steady state, gross investment just covers depreciation: $I^* = \delta K^*$.

**(b)** $I_0/K_0 = (q_0-1)/\psi = 0.5/10 = 0.05$. The firm invests at 5% of the capital stock per period in addition to replacement investment. 

**(c)** The sluggish response — only 5% of the capital stock is installed per period despite $q$ being 50% above one — reflects the high adjustment cost. Convex adjustment costs smooth investment over time. The higher $\psi$, the slower firms adjust their capital stock to shocks. Empirically, this matches the observation that investment responds smoothly and gradually to changes in the interest rate or profitability rather than in discrete jumps.

---

## Parts IV–IX: Hints for Selected Exercises

**Chapter 22, Exercise 22.2.** To verify Ricardian equivalence in a two-period model: write the household's budget constraint in both periods, add them together to get the intertemporal constraint, substitute the government budget constraint, and show that the household's lifetime resources are unchanged by debt versus tax finance.

**Chapter 27, Exercise 27.1.** To show that the RBC model matches the variance of output, begin by log-linearizing the model around the steady state, solving for the policy functions (consumption and investment as functions of capital and productivity), and simulating with the calibrated shock process. The solution method involves either the method of undetermined coefficients or numerical eigendecomposition.

**Chapter 34, Exercise 34.2** ($\star\star$). To derive the Minsky cycle dynamics from a balance sheet model with endogenous leverage, write the model of Brunnermeier and Sannikov (2014) in which the fraction of the economy financed by expert intermediaries evolves stochastically, and show that multiple equilibria — a high-wealth stable regime and a low-wealth crisis regime — can emerge.

**Chapter 38, Exercise 38.1.** To compute the Gini coefficient for a log-normal income distribution with mean $\mu$ and variance $\sigma^2$: the Gini is $G = 2\Phi(\sigma/\sqrt{2}) - 1$, where $\Phi$ is the standard normal CDF. Verify this formula using the definition of $G$ and the properties of the log-normal distribution.

---

*Full worked solutions to all exercises are available in the accompanying solutions manual. Instructors may request access from the publisher.*
