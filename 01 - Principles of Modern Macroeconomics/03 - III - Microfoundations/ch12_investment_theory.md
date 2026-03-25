# Chapter 12 — Investment Theory: Firms' Decisions on Capital and Inventory

> *"Business men play a mixed game of skill and chance, the average results of which to the players are not known by those who take a hand."*
> — Keynes, *The General Theory*, 1936

---

Investment is the most volatile component of aggregate expenditure — roughly three times more variable than output over the business cycle — and understanding what drives it is essential for understanding both recessions and recoveries. Investment is also the mechanism through which the capital stock grows, linking short-run fluctuations to the long-run growth process studied in Chapter 5. And it is the primary channel through which monetary policy affects real activity: when the central bank raises interest rates, it is investment spending — on equipment, structures, and housing — that falls first and most sharply.

Yet investment is also the most conceptually challenging component of aggregate demand. Unlike consumption, which responds to current income and wealth, investment decisions are fundamentally forward-looking: a firm invests in a new machine today because it expects that machine to generate profits over its entire productive life, which may be ten or twenty years. Any theory of investment must therefore confront the problem of expectations about the future under genuine uncertainty.

---

## 12.1 The Neoclassical Theory of Investment

The neoclassical approach, developed by Jorgenson (1963), asks: what is the optimal capital stock for a profit-maximizing firm, and how quickly does the firm move toward it? The firm's objective is to maximize the present discounted value of profits:

$$V_0 = \int_0^\infty e^{-rt}\bigl[F(K_t, L_t) - w L_t - p^I I_t\bigr]\,\mathrm{d}t,$$

where $F(K, L)$ is the production function, $w$ is the real wage, $p^I$ is the real price of investment goods, and $I_t = \dot{K}_t + \delta K_t$ is gross investment (the amount needed to produce net capital accumulation plus replacement of depreciated capital). Solving the optimal control problem (with $K_t$ as the state variable and $I_t$ as the control) yields the first-order condition:

$$F_K(K_t, L_t) = r + \delta \equiv c^K.$$

**Definition (User Cost of Capital).** The **user cost of capital** $c^K = r + \delta$ is the implicit rental price of owning one unit of capital for one period: the opportunity cost of funds tied up in capital ($r$) plus the depreciation that occurs during the period ($\delta$). A profit-maximizing firm hires capital until its marginal product equals this user cost, exactly as it hires labor until the marginal product of labor equals the real wage.

In a tax environment with corporate tax rate $\tau_c$, investment tax credit $k^{ITC}$, and present value of depreciation allowances $D^{PV}$, Hall and Jorgenson (1967) derive the after-tax user cost:

$$c^K = \frac{(1-k^{ITC})(1-\tau_c D^{PV})}{1-\tau_c}(r + \delta).$$

The critical prediction of the neoclassical model is that investment adjusts instantaneously to equate the marginal product of capital to $c^K$. Whenever the interest rate changes, tax policy changes, or technology improves, the firm immediately adjusts its capital stock to the new optimal level. This implies implausibly large and rapid investment responses — the data show that investment adjusts slowly and smoothly, not in discrete jumps. The resolution of this contradiction is the adjustment cost model.

---

## 12.2 Tobin's q and the Market Value of Capital

James Tobin (1969) proposed an elegant theory: firms invest when the market value of installed capital exceeds its replacement cost, and they disinvest or allow depreciation when the market value falls short of replacement cost.

**Definition (Tobin's q).** **Tobin's q** is the ratio of the market value of a firm's installed capital to its replacement cost:

$$q_t \equiv \frac{\text{Market value of installed capital}}{\text{Replacement cost of capital}} = \frac{V_t}{p_t^I K_t},$$

where $V_t$ is the stock market value of the firm and $p_t^I K_t$ is the cost of replacing all physical capital at current investment goods prices.

When $q_t > 1$, installing a new unit of capital is worth more than it costs: the market values the installed capital above its production cost. Firms should invest until the rising capital stock drives down the marginal product (and hence $q$) back to one. When $q_t < 1$, the market values the firm's capital below what it would cost to replace — investment should cease and the capital stock should contract through depreciation.

The intuition for $q > 1$ is that capital, once installed, has some degree of irreversibility: it cannot be costlessly converted back to cash. Installed capital therefore commands a premium over uninstalled capital if the firm is profitable, because it is already in place and generating returns. This premium is captured by $q$.

---

## 12.3 Adjustment Costs and the $q$ Model

To generate smooth investment dynamics consistent with the data, the standard approach adds **convex adjustment costs** to the firm's problem. These costs capture the disruption and installation expenses of installing new capital rapidly: it is cheaper to install a small amount of new capacity per period than a large amount all at once.

Assume total adjustment costs are $C(I, K) = (\psi/2)(I/K)^2 K$, which is strictly convex in the investment rate $I/K$. This specification implies that doubling the rate of investment more than doubles the cost — there are increasing marginal costs of rapid installation. The firm maximizes:

$$V_0 = \int_0^\infty e^{-rt}\!\left[F(K_t, L_t) - wL_t - I_t - \frac{\psi}{2}\!\left(\frac{I_t}{K_t}\right)^2 K_t\right]\mathrm{d}t,$$

subject to $\dot{K}_t = I_t - \delta K_t$. Forming the current-value Hamiltonian with costate variable $q_t$ (the shadow value of an additional unit of installed capital), the optimality conditions are:

$$1 + \psi\frac{I_t}{K_t} = q_t \quad \text{(optimality in } I_t\text{)},$$

$$\dot{q}_t = rq_t - F_K(K_t, L_t) + \frac{\psi}{2}\!\left(\frac{I_t}{K_t}\right)^2 \quad \text{(costate equation)}.$$

The first condition gives the **investment function**:

$$\frac{I_t}{K_t} = \frac{q_t - 1}{\psi}.$$

This is the central result of the $q$ model: the investment rate is proportional to the excess of $q$ over one. Investment is positive only when $q > 1$ and increases with the gap. The parameter $\psi$ governs the speed of adjustment: large $\psi$ (high adjustment costs) means the firm moves slowly toward its optimal capital stock even when $q$ is far from one.

**Hayashi's theorem** (Hayashi, 1982) establishes that when both the production function and the adjustment cost function exhibit constant returns to scale, the marginal $q$ appearing in the investment function equals the *observable* average $q = V_t/(p_t^I K_t)$ — the ratio of the stock market value to the replacement cost of capital. This result is crucial for empirical work: it allows researchers to test the investment model using the observable stock market valuation rather than the unobservable shadow value.

In practice, empirical $q$ models have disappointing explanatory power — average $q$ typically explains only 5–10% of investment variation. The most important reasons are: measurement error in $q$ (accounting-based capital stock measures poorly approximate replacement cost), financial constraints that drive a wedge between the firm's borrowing rate and the risk-free rate, and the presence of irreversibilities and uncertainty, to which we now turn.

---

## 12.4 Irreversibility, Uncertainty, and Real Options

The neoclassical and $q$ models both treat investment as continuously reversible: the firm can install or remove capital freely. In reality, much capital investment is largely irreversible — a factory, a building, or specialized equipment cannot be sold for anywhere near its purchase price once installed. This irreversibility has a profound effect on optimal investment behavior, because it creates an **option value of waiting** that a reversible model ignores.

**Definition (Irreversibility).** An investment is **irreversible** if the installed capital cannot be costlessly converted back to other uses or sold for its purchase price. Irreversibility means that investing today destroys the option to invest tomorrow under potentially more favorable conditions. The value of preserving this option — not investing today so as to retain the flexibility to invest tomorrow — is the **option value of waiting**.

The analogy with financial options is exact. A firm holding an investment opportunity is in the same position as a holder of a call option on an asset: it has the right, but not the obligation, to invest (exercise the option) at any time by paying the investment cost (the exercise price). The value of this option depends on the current profitability of the project, the volatility of future profitability, and the cost of investment. Option pricing theory — specifically, the framework of Dixit and Pindyck (1994) — provides the analytical tools.

Assume future project cash flows $\Pi_t$ follow a geometric Brownian motion:

$$\mathrm{d}\Pi_t = \mu\Pi_t\,\mathrm{d}t + \sigma\Pi_t\,\mathrm{d}W_t,$$

where $\mu$ is the drift, $\sigma$ is the volatility, and $W_t$ is a standard Brownian motion. The firm invests when $\Pi_t$ exceeds a **trigger value** $\Pi^*$. Using the Bellman equation for the option value and imposing value-matching and smooth-pasting conditions, the optimal trigger satisfies:

$$\Pi^* = \frac{\beta}{\beta - 1}\, c^K, \qquad \beta = \frac{1}{2} - \frac{\mu}{\sigma^2} + \sqrt{\!\left(\frac{\mu}{\sigma^2} - \frac{1}{2}\right)^2 + \frac{2r}{\sigma^2}} > 1.$$

The factor $\beta/(\beta - 1) > 1$ means the investment trigger $\Pi^*$ exceeds the neoclassical threshold $c^K$: the firm requires a premium above the user cost before investing, because investing destroys the option to wait. Crucially, **higher volatility $\sigma$ raises $\Pi^*$**: greater uncertainty increases the value of the option to wait, making firms more reluctant to invest. This real-options mechanism explains why investment collapses sharply in recessions when uncertainty rises — measured, for instance, by the VIX or by cross-sectional dispersion of firm-level returns — even when the expected return on investment remains positive (Bloom, 2009).

---

*Next: Chapter 13 — Labor Supply and Demand*
