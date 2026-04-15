# Chapter 12 — Investment Theory: Firms' Decisions on Capital and Inventory

> *"Business men play a mixed game of skill and chance, the average results of which to the players are not known by those who take a hand."*
> — Keynes, *The General Theory*, 1936

---

Investment is the most volatile component of aggregate expenditure — roughly three to four times more variable than output over the business cycle — and understanding its determinants is essential for understanding both recessions and recoveries. Investment is also the mechanism through which the capital stock grows and productivity advances, linking short-run cyclical fluctuations to the long-run growth process studied in Chapter 5. And it is the primary channel through which monetary policy affects real activity: when the central bank raises interest rates, it is investment spending — on equipment, structures, software, and residential housing — that falls first and most sharply.

Yet investment is the most conceptually difficult component of aggregate demand. Unlike consumption, which responds predictably to current income and wealth, investment decisions are fundamentally forward-looking under genuine, irreducible uncertainty. A firm that builds a factory today is committing resources to a project whose returns will be realized over ten or twenty years, under economic and technological conditions that are unknown at the time of investment. Any adequate theory of investment must confront this intertemporal problem honestly.

This chapter develops three approaches. The neoclassical theory derives the optimal capital stock from profit maximization, yielding the user cost of capital as the key determinant. Tobin's $q$ theory connects investment to the stock market, providing a theoretically elegant and empirically testable reduced form. The adjustment cost model generates the smooth investment dynamics observed in the data. Finally, the real options approach shows that irreversibility and uncertainty — ignored by the preceding frameworks — fundamentally alter the investment decision.

---

## 12.1 The Neoclassical Theory of Investment

### The Optimal Capital Stock

The neoclassical approach, formalized by Jorgenson (1963), begins by asking what capital stock maximizes the present value of a competitive firm's profits. The firm maximizes:

$$V_0 = \int_0^\infty e^{-rt}\bigl[F(K_t, L_t) - w L_t - p^I I_t\bigr]\,\mathrm{d}t,$$

where $F(K,L)$ is the production function satisfying $F_K, F_L > 0$ and $F_{KK}, F_{LL} < 0$, $w$ is the real wage, $p^I$ is the real price of investment goods, and $I_t = \dot{K}_t + \delta K_t$ is gross investment. Solving the optimal control problem with $K_t$ as the state variable and $I_t$ as the control, the first-order condition for the optimal capital stock is:

$$F_K(K_t, L_t) = r + \delta \equiv c^K.$$

**Definition (User Cost of Capital).** The **user cost of capital** $c^K = r + \delta$ is the implicit rental price of owning one unit of capital for one period. It has two components: the opportunity cost of the funds invested ($r$, which represents what the same funds could earn in a financial investment) and the depreciation that occurs during the period ($\delta$). A profit-maximizing firm expands its capital stock until the marginal product of capital equals this user cost — exactly as it expands employment until the marginal product of labor equals the real wage.

### The Tax-Adjusted User Cost

In a realistic tax environment with corporate tax rate $\tau_c$, investment tax credit $k^{ITC}$, and present value of depreciation allowances $D^{PV}$, Hall and Jorgenson (1967) derive the after-tax user cost:

$$c^K = \frac{(1 - k^{ITC})(1 - \tau_c D^{PV})}{1 - \tau_c}(r + \delta).$$

Tax policy enters directly: an investment tax credit ($k^{ITC} > 0$) reduces the effective cost of acquiring capital; accelerated depreciation (which raises $D^{PV}$) does the same through the tax shield on depreciation. These instruments have been used deliberately as levers of investment policy — the Kennedy investment tax credit of 1962 and the Reagan accelerated depreciation schedules of 1981 are the most prominent U.S. examples — and the Jorgenson formula provides the quantitative link from tax parameters to the cost of capital.

### The Instantaneous Adjustment Problem

The critical prediction of the neoclassical model is that the capital stock adjusts instantaneously to its optimal level: whenever $r$, $\tau_c$, or technology changes, firms jump to the new optimal $K^*$ in one period. This prediction is empirically false — investment adjusts slowly and smoothly, not in discrete jumps. The resolution is the adjustment cost model.

---

## 12.2 Tobin's $q$ and the Market Value of Capital

James Tobin (1969) proposed a more direct link between financial markets and investment decisions: firms invest when the market value of installed capital exceeds its replacement cost, and disinvest (or allow depreciation to shrink the capital stock) when market value falls short of replacement cost.

**Definition (Tobin's $q$).** **Tobin's $q$** is the ratio of the market value of a firm's installed capital to its replacement cost:

$$q_t \equiv \frac{\text{Market value of installed capital}}{\text{Replacement cost of capital}} = \frac{V_t}{p_t^I K_t},$$

where $V_t$ is the firm's stock market capitalization and $p_t^I K_t$ is the cost of replacing all its physical capital at current investment goods prices.

### Interpretation

When $q_t > 1$, the market values each unit of installed capital at more than it costs to install a new one: expansion is profitable and the firm should invest. When $q_t < 1$, the market values the firm's capital below replacement cost: the firm should contract by allowing depreciation to exceed gross investment. The investment decision rule $I > 0$ iff $q > 1$ follows directly from the stock market's forward-looking valuation of the firm's future profits from capital.

The power of Tobin's $q$ as an empirical approach is that it summarizes all the information about future profitability — technology, demand, interest rates, taxes, competition — that a rational market would incorporate into the stock price. Rather than measuring each determinant separately, the researcher can use $q$ as a single sufficient statistic for investment incentives.

### Limitations of the $q$ Approach

In practice, empirical $q$ models have disappointing explanatory power: average $q$ typically explains only 5–10% of variation in investment in firm-level data. The most important reasons are measurement error in $q$ (the book value of capital, used to approximate the replacement cost $p^I K$, is a poor proxy for the true replacement value especially for intangible assets); financial constraints that drive a wedge between the social discount rate and the firm's borrowing cost; and the presence of irreversibility and adjustment costs, which mean that the standard $q$ model must be augmented.

---

## 12.3 Adjustment Costs and the $q$ Model

### Setup

To generate the smooth investment dynamics observed in the data, the standard approach introduces **convex adjustment costs**: installing new capital rapidly is more expensive than installing it slowly, due to organizational disruption, installation complexity, and the learning required to use new equipment productively. These costs are assumed to be $C(I, K) = (\psi/2)(I/K)^2 K$, strictly convex in the investment rate $I/K$. The firm maximizes:

$$V_0 = \int_0^\infty e^{-rt}\!\left[F(K_t, L_t) - wL_t - I_t - \frac{\psi}{2}\!\left(\frac{I_t}{K_t}\right)^2 K_t\right]\mathrm{d}t,$$

subject to $\dot{K}_t = I_t - \delta K_t$.

### Optimality Conditions

Forming the current-value Hamiltonian with costate variable $q_t$ (the shadow value of an additional unit of installed capital, i.e., the marginal $q$):

$$1 + \psi\frac{I_t}{K_t} = q_t \quad \text{(optimality in } I_t\text{)},$$

$$\dot{q}_t = rq_t - F_K(K_t, L_t) + \frac{\psi}{2}\!\left(\frac{I_t}{K_t}\right)^2 \quad \text{(costate equation)}.$$

### The Investment Function

Rearranging the first optimality condition gives the central result:

$$\frac{I_t}{K_t} = \frac{q_t - 1}{\psi}.$$

The investment rate is proportional to the excess of $q$ over one. When $q = 1$, the firm maintains its capital stock with zero net investment. When $q > 1$, the firm expands; when $q < 1$, it contracts. The parameter $\psi$ governs adjustment speed: high $\psi$ (steep adjustment costs) means the firm moves slowly toward its optimal capital stock even when the gap between $q$ and one is large.

This smooth investment function, derived from optimality rather than assumed, is what generates the gradual, procyclical investment dynamics observed in the data. During recessions, $q$ falls below one as stock prices fall and profit expectations worsen; the firm reduces investment but does so gradually, proportional to $\psi^{-1}(q-1)$. During booms, $q$ exceeds one and investment accelerates. The model therefore predicts that the stock market should lead investment — a prediction broadly confirmed by the data, though the quantitative relationship is weaker than the simple model suggests.

**Hayashi's theorem** (1982) establishes the crucial empirical bridge: when both the production function and the adjustment cost function exhibit constant returns to scale, the unobservable marginal $q$ equals the observable average $q = V_t/(p_t^I K_t)$ — the ratio of stock market value to replacement cost. This provides the empirical justification for using stock market valuations to test the investment model.

---

## 12.4 Irreversibility, Uncertainty, and Real Options

### The Option to Invest

The neoclassical and $q$ models both treat investment as continuously reversible: installed capital can be removed and sold for its purchase price at no cost. In reality, much capital investment is largely irreversible. A purpose-built factory, specialized equipment, or a mine has high value in its intended use but low resale value — the "lemons" problem of thinly traded second-hand capital markets means that installed capital sells for substantially below its purchase price. This irreversibility has a profound effect on optimal investment behavior that the preceding frameworks miss entirely.

**Definition (Option Value of Waiting).** When investment is irreversible, a firm that invests today permanently surrenders the **option to invest tomorrow** under potentially more favorable conditions. The **option value of waiting** is the value of retaining this flexibility. A firm should invest only when the expected return on the project exceeds not only the direct cost of investment but also the value of the option that investment destroys. This creates an **investment trigger** above the neoclassical threshold.

### The Dixit-Pindyck Model

The analogy with financial options is exact: a firm holding an investment opportunity is like a holder of a perpetual American call option on the project's value. Dixit and Pindyck (1994) formalize this using continuous-time stochastic calculus. Assume project cash flows $\Pi_t$ follow a geometric Brownian motion:

$$\mathrm{d}\Pi_t = \mu\Pi_t\,\mathrm{d}t + \sigma\Pi_t\,\mathrm{d}W_t,$$

where $\mu$ is the expected growth rate of cash flows, $\sigma > 0$ is their volatility, and $W_t$ is a standard Brownian motion. The firm invests when $\Pi_t$ first exceeds a trigger value $\Pi^*$. Applying the Bellman equation and the value-matching and smooth-pasting conditions (which ensure the option is exercised optimally):

$$\Pi^* = \frac{\beta}{\beta - 1}\, c^K, \qquad \beta = \frac{1}{2} - \frac{\mu}{\sigma^2} + \sqrt{\!\left(\frac{\mu}{\sigma^2} - \frac{1}{2}\right)^2 + \frac{2r}{\sigma^2}} > 1.$$

### The Key Predictions

The ratio $\beta/(\beta - 1) > 1$ means the investment trigger $\Pi^*$ exceeds the neoclassical user cost $c^K$. The firm requires a **premium above the user cost** before investing — the premium compensates for the destroyed option. Several comparative statics follow:

**Uncertainty raises the investment threshold.** $\partial\Pi^*/\partial\sigma > 0$: greater volatility in cash flows increases the option value of waiting (the asymmetry of the option payoff — you benefit from upside but can wait out the downside — makes waiting more valuable when volatility is high). This provides the theoretical underpinning for the empirical finding that uncertainty reduces investment even when expected returns are positive.

**Higher interest rates raise the threshold through two channels.** Directly, $c^K = r + \delta$ rises. Indirectly, higher $r$ increases $\beta$, which reduces $\beta/(\beta-1)$ — but the direct user cost effect dominates, making higher rates unambiguously contractionary for investment.

**Irreversibility amplifies investment sensitivity to demand shocks.** In a fully reversible world, firms install capital whenever $\Pi > c^K$; in an irreversible world, they wait until $\Pi > (\beta/(\beta-1))c^K$. When a negative shock lowers $\Pi$ below $c^K$, reversible firms immediately disinvest; irreversible firms hold their capital (since they cannot recover installation costs) and stop new investment. This asymmetric response — quick contraction, slow expansion — generates the observed cyclical asymmetry in investment: investment recovers more slowly from recessions than it falls into them.

---

## 12.5 The Financial Accelerator

The models above treat investment as if firms can freely borrow at the market interest rate $r$. In reality, information asymmetries between borrowers and lenders mean that external finance — borrowing — is more expensive than internal finance (retained earnings), and the wedge between the two (the **external finance premium**) varies with the borrower's balance sheet health.

**Definition (External Finance Premium).** The **external finance premium** $\rho_t$ is the spread between the cost of external borrowing and the risk-free rate, reflecting agency costs of external finance: $1 + \rho_t = f(N_t/Q_t K_t)$, where $N_t$ is the firm's net worth and $Q_t K_t$ is the value of its capital stock. When net worth is high relative to assets, lenders face lower adverse selection and moral hazard risks, and the premium is low; when net worth is low (as in recessions), the premium rises.

This creates a **financial accelerator** (Bernanke, Gertler, and Gilchrist, 1999): in a recession, falling asset prices reduce firm net worth; lower net worth raises the external finance premium; higher borrowing costs further reduce investment; lower investment reduces income and asset prices further — a self-reinforcing cycle that amplifies the initial shock. Conversely, in booms, rising asset prices raise net worth, compress the premium, stimulate investment, and amplify the expansion.

The financial accelerator is why financial crises generate such severe and persistent investment collapses: the collapse in financial market asset prices simultaneously reduces net worth (raising the premium) and collapses $q$ (reducing the incentive to invest). Both forces point in the same direction. The Great Recession of 2008–09, during which U.S. business investment fell by approximately 25%, is the canonical illustration [Ch. 40].

---

## Chapter Summary

- The **neoclassical model** yields the optimal capital stock by equating $F_K = c^K = r + \delta$: the marginal product of capital must equal the user cost. Tax policy (investment tax credits, depreciation allowances) shifts the after-tax user cost and hence the desired capital stock. The model predicts instantaneous adjustment, which fails empirically.

- **Tobin's $q$** (market value/replacement cost) provides an observable sufficient statistic for investment incentives under constant returns to scale (Hayashi's theorem). Empirical $q$ models explain only 5–10% of investment variance, reflecting measurement error, financial constraints, and irreversibility.

- The **adjustment cost model** yields the investment function $I/K = (q-1)/\psi$, generating smooth investment dynamics consistent with data. Investment is proportional to the gap between $q$ and one; the speed of adjustment decreases with the adjustment cost parameter $\psi$.

- The **real options framework** (Dixit-Pindyck) shows that irreversibility creates an option value of waiting, raising the investment trigger to $\Pi^* = (\beta/(\beta-1))c^K > c^K$. Higher uncertainty raises $\Pi^*$, reducing investment even at unchanged expected returns — the theoretically grounded explanation for the sharp investment collapse during high-uncertainty episodes (Bloom, 2009).

- The **financial accelerator** (BGG): information asymmetries create an external finance premium $\rho_t = f(N_t/Q_tK_t)$ that varies inversely with firm net worth. In recessions, falling net worth raises the premium and depresses investment, amplifying initial shocks through a self-reinforcing cycle.

---

*Next: Chapter 13 — Labor Supply and Demand*
