# Chapter 11 — Consumption Theory: How Households Decide to Spend and Save

> *"It is not rational to decide on consumption by computing marginal utilities on the spot. The decision is made once and for all."*
> — Franco Modigliani

---

Private consumption is typically the largest component of GDP — around 70% in the United States — and its behavior over the business cycle is central to macroeconomic dynamics. When households cut spending, they set off the multiplier chains of Chapter 8 in reverse; when they run down saving to maintain consumption during recessions, they provide an automatic stabilizing force. Understanding what drives household consumption and saving decisions is therefore one of the most important questions in macroeconomics, with implications for fiscal policy, growth theory, and the monetary transmission mechanism.

The Keynesian consumption function $C = a + b(Y - T)$ provided the first systematic answer: consumption is a stable increasing function of current disposable income. This was initially empirically successful — the cross-sectional relationship between income and consumption looks linear and stable in household survey data. But two empirical puzzles quickly emerged. First, the time-series relationship between consumption and income is less steep than the cross-sectional relationship — the time-series MPC is lower than the cross-sectional MPC. Second, the Keynesian prediction that the saving rate should rise with income as countries develop — because higher income should allow more saving — was contradicted by Simon Kuznets's finding that the U.S. saving rate was approximately constant over long historical periods despite large income growth.

These puzzles motivated the life-cycle and permanent income theories that replaced the Keynesian consumption function as the foundation of modern consumption theory.

---

## 11.1 The Life-Cycle Hypothesis

Modigliani and Brumberg (1954) proposed the **life-cycle hypothesis (LCH)**: households plan their consumption and saving over their entire lifetimes, smoothing consumption across periods of high income (working years) and low income (childhood and retirement). The unit of analysis is not a single period's income but the household's entire lifetime resources.

Consider a household that lives $T$ periods, works and earns income $y$ for the first $T_W$ periods, then retires for $T_R = T - T_W$ periods with zero labor income. Let $A_0$ denote initial wealth. Under the simplifying assumptions that the real interest rate is zero and that the household wants to maintain constant consumption $c$ throughout life:

$$c \cdot T = A_0 + y \cdot T_W \implies c = \frac{A_0 + T_W\, y}{T}.$$

This formula has several important implications. The **MPC out of current labor income** is $T_W/T$ — the fraction of remaining life still in the working phase. For a middle-aged worker with 20 years left to work and 30 years until death, $T_W/T = 20/50 = 0.4$: only 40 cents of each additional dollar of income is consumed, the rest is saved for retirement. This is far below Keynes's MPC of $b \approx 0.8$. The **MPC out of wealth** is $1/T$: a household that inherits a lump sum of wealth spreads the windfall over its remaining lifetime.

The resolution of the two empirical puzzles follows immediately. The cross-sectional MPC is high because cross-sectional income variation is dominated by age effects and transitory factors — households in their peak earning years have high income and save heavily because $T_W/T$ is high, but this is not the marginal propensity to consume when income changes permanently. The time-series saving rate is approximately constant because, in an economy growing at a steady rate, each generation is wealthier than the previous one by the same proportion, and the aggregate saving of young high-earners exactly offsets the dissaving of retirees in steady state.

The life-cycle hypothesis also implies a **hump-shaped life-cycle saving profile**: households save during their working years (building up retirement wealth) and dissave during retirement (drawing down accumulated wealth). This prediction is broadly consistent with household survey data, though wealth in old age is higher than the pure life-cycle model predicts — suggesting that bequest motives, precautionary saving, or uncertainty about the date of death also influence saving (Kotlikoff and Summers, 1981).

---

## 11.2 The Permanent Income Hypothesis

Friedman (1957) proposed the **permanent income hypothesis (PIH)**: households base their consumption on their **permanent income** $y^P$ — a measure of their long-run average income — rather than on their current income.

**Definition (Permanent Income).** **Permanent income** $y_t^P$ is the level of income that, if received as a constant perpetuity, would generate the same present discounted value as the household's entire expected future income stream. Formally, if the real interest rate is $r$ and the household has infinite horizon:

$$y_t^P = r \cdot \left[a_t + \sum_{s=0}^\infty \frac{\mathbb{E}_t[y_{t+s}]}{(1+r)^s}\right],$$

where $a_t$ is current wealth. Permanent income is the "interest" on the household's total wealth — financial wealth plus the present value of human capital.

Friedman decomposed actual income into a permanent component and a **transitory component**:

$$y_t = y_t^P + y_t^T, \quad \mathbb{E}[y_t^T] = 0,$$

where transitory income $y_t^T$ represents unexpected, one-off fluctuations around the permanent level. The PIH hypothesis then states:

$$c_t = \alpha\, y_t^P$$

for some $\alpha$ close to one. Transitory income shocks are saved; only permanent income changes are consumed.

The PIH resolves both empirical puzzles. The cross-sectional income–consumption relationship is steeper because most of the cross-sectional variation in income reflects permanent differences (education, skills, lifetime earnings potential) that affect $y^P$ one-for-one. The time-series MPC is lower because short-run income fluctuations are partly transitory: a temporary bonus or a one-off recession should generate little change in permanent income and hence little change in consumption.

A key policy implication: the consumption response to fiscal transfers depends on whether they are perceived as permanent or transitory. A permanent tax cut raises $y^P$ significantly; a one-time rebate is largely saved because it barely moves $y^P$. This prediction is broadly consistent with empirical evidence from U.S. tax rebates, though the MPC is not zero — it is typically around 0.3–0.5 — suggesting constraints that prevent households from perfectly smoothing consumption.

---

## 11.3 The Euler Equation and Optimal Intertemporal Consumption

The life-cycle and permanent income frameworks provide economic intuition but are not derived from explicit utility maximization. The Euler equation approach, due to Hall (1978), provides the rigorous microfoundation.

A household maximizes expected lifetime utility:

$$U_t = \mathbb{E}_t\sum_{s=0}^\infty \beta^s\, u(c_{t+s}),$$

where $\beta = 1/(1+\rho) \in (0,1)$ is the discount factor (with $\rho > 0$ the subjective rate of time preference, measuring how impatient the household is), and $u(\cdot)$ is a strictly concave, increasing instantaneous utility function. The household's flow budget constraint is:

$$a_{t+1} = (1 + r_t)(a_t + y_t - c_t),$$

stating that wealth next period equals wealth this period plus saving ($a_t + y_t - c_t$) scaled by the gross return.

The Bellman equation for this dynamic program:

$$V(a_t, y_t) = \max_{c_t}\bigl\{u(c_t) + \beta\,\mathbb{E}_t\,V(a_{t+1}, y_{t+1})\bigr\}.$$

The first-order condition from the Bellman equation, combined with the envelope condition, yields the **consumption Euler equation**:

$$u'(c_t) = \beta(1+r_{t+1})\,\mathbb{E}_t[u'(c_{t+1})].$$

**Definition (Euler Equation for Consumption).** The **Euler equation** is the first-order condition for optimal intertemporal consumption allocation. It states that the household has no incentive to shift consumption between periods: the marginal utility lost by consuming one dollar less today must equal the discounted expected marginal utility gained by consuming the resulting savings next period. If $\beta(1+r) > 1$ (the return on saving exceeds the rate of impatience), the household optimally chooses rising consumption; if $\beta(1+r) < 1$, consumption declines over time.

With CRRA utility $u(c) = c^{1-\sigma}/(1-\sigma)$, the Euler equation becomes $c_t^{-\sigma} = \beta(1+r_{t+1})\mathbb{E}_t[c_{t+1}^{-\sigma}]$. Taking logarithms and applying a second-order Taylor expansion (to account for the covariance between future consumption and the interest rate):

$$\mathbb{E}_t[\Delta\ln c_{t+1}] = \frac{r_{t+1} - \rho}{\sigma} + \frac{\sigma}{2}\,\mathrm{Var}_t[\Delta\ln c_{t+1}].$$

The first term says consumption grows when the real return on saving exceeds the rate of impatience. The second term is the **precautionary saving effect**: when future consumption is uncertain (high variance), risk-averse households reduce current consumption and save more as a buffer against bad outcomes. The coefficient of relative risk aversion $\sigma$ governs both the sensitivity of consumption growth to the interest rate and the strength of the precautionary motive.

---

## 11.4 Hall's Random Walk Hypothesis and Excess Sensitivity

Hall (1978) derived a striking implication of the Euler equation under quadratic utility: if $u(c) = -(b-c)^2/2$ (which implies $u''' = 0$, eliminating the precautionary motive), the Euler equation simplifies to $c_t = \mathbb{E}_t[c_{t+1}]$. The level of consumption is a **martingale**: its best predictor for next period is its current level.

**Definition (Martingale).** A stochastic process $\{c_t\}$ is a **martingale** if $\mathbb{E}_t[c_{t+1}] = c_t$ for all $t$ — the expected change is zero. If consumption is a martingale, then $c_{t+1} - c_t = \epsilon_{t+1}$, where $\epsilon_{t+1}$ is unpredictable from any information available at $t$. No variable in the information set $\mathcal{F}_t$ should help predict the change in consumption.

This is the **random walk hypothesis for consumption**: changes in consumption should be white noise, unpredictable from past information including past income. The policy implication is that consumption is unresponsive to predictable income changes — only genuinely surprising news should affect consumption.

Hall tested this hypothesis and found that it was rejected: lagged income growth has significant predictive power for consumption growth. This finding of **excess sensitivity** — consumption responds more strongly to predictable income changes than the PIH–Euler equation predicts — has been one of the most robust findings in consumption economics for decades. The two most compelling explanations are liquidity constraints (households who cannot borrow to smooth consumption must follow current income even when future income is predictably high) and hand-to-mouth behavior (a large fraction of households are near their borrowing limit at all times).

---

## 11.5 Liquidity Constraints and Buffer-Stock Saving

Many households cannot borrow freely at the risk-free rate. They face either quantity constraints (explicit limits on borrowing) or price constraints (interest rates on consumer debt far above the risk-free rate). In either case, their consumption cannot be fully smoothed across periods.

Suppose a household faces a binding borrowing constraint $a_{t+1} \geq \underline{b}$. The Kuhn–Tucker complementary slackness condition adds a non-negative multiplier $\mu_t$ on this constraint to the Euler equation:

$$u'(c_t) = \beta(1+r)\,\mathbb{E}_t[u'(c_{t+1})] + \mu_t, \quad \mu_t \geq 0.$$

When the constraint binds ($a_{t+1} = \underline{b}$, $\mu_t > 0$), the household is consuming more than it "should" relative to the unconstrained optimum — it is consuming exactly its current income rather than saving for future uncertainty. In this regime, consumption tracks current income one-for-one, reproducing the excess sensitivity finding as a rational response to borrowing constraints rather than a violation of optimization.

The **buffer-stock model** (Deaton, 1991; Carroll, 1997) combines borrowing constraints with prudence ($u''' > 0$) to generate a theory of precautionary wealth accumulation. Households target a ratio of assets to permanent income $A^* = a_t/y_t^P$, holding enough wealth to buffer against income shocks without accumulating beyond what precautionary motives justify. When $a_t > A^* y_t^P$, the household runs down assets toward target; when $a_t < A^* y_t^P$, it saves aggressively. The buffer-stock model generates realistic patterns of wealth-income ratios at the household level and explains why the aggregate consumption function is approximately proportional to income in the long run (consistent with Kaldor facts) but exhibits excess sensitivity in the short run.

---

*Next: Chapter 12 — Investment Theory*
