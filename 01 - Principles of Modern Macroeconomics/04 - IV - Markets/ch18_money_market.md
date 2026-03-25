# Chapter 18 — The Money Market: Interest Rates and Liquidity

> *"The essence of the monetary system is the confidence that underlies it."*
> — Ben S. Bernanke

---

The money market is the market in which the supply of and demand for short-term liquidity is equilibrated. It is the primary channel through which central banks transmit policy impulses to the broader economy: by controlling the quantity of reserves supplied to the banking system (or, equivalently, by setting a target for the overnight interest rate), the central bank affects the cost of liquidity throughout the financial system, which in turn propagates to longer-term interest rates, asset prices, exchange rates, and ultimately to consumption and investment. Understanding the money market in precise terms — what is traded, how demand and supply interact, how the equilibrium interest rate is determined, and what happens at the effective lower bound — is prerequisite to understanding monetary policy transmission.

This chapter develops the money market in three steps. First, we characterize money demand, drawing on both the inventory-theoretic model of Chapter 14 and the portfolio-theoretic approach. Second, we analyze the supply of money through the banking system and the monetary base. Third, we examine money-market equilibrium and its breakdown in the liquidity trap — one of the most important limiting cases in all of monetary economics.

---

## 18.1 The Demand for Money: A Multi-Motive Analysis

The demand for money arises from three distinct motives, each of which generates a separable component of total money demand.

### The Transactions Motive

The **transactions motive** arises because money eliminates the double coincidence of wants problem. Households and firms must hold money between income receipts and expenditure payments to facilitate transactions. The transactions demand for money increases with the volume of transactions — proxied by nominal income $PY$ — and decreases with the nominal interest rate, which measures the cost of holding non-interest-bearing money rather than bonds.

The Baumol–Tobin (1952, 1956) model, developed formally in Chapter 14, gives the transactions demand:

$$L^T = \sqrt{\frac{bPY}{2i}},$$

where $b$ is the fixed cost per bank trip and $i$ is the nominal interest rate. In real terms, dividing by $P$:

$$\ell^T = \frac{L^T}{P} = \sqrt{\frac{bY}{2i}} \propto Y^{1/2}\, i^{-1/2}.$$

The income elasticity of transactions money demand is $1/2$, less than one: a 1% increase in real income raises money demand by 0.5%. The interest elasticity is $-1/2$.

### The Precautionary Motive

Beyond planned transactions, households hold money to guard against unforeseen expenditures — unexpected medical bills, urgent repairs, sudden travel. The **precautionary demand** for money increases with income uncertainty and decreases with the nominal interest rate (since higher rates make it more costly to hold precautionary balances in cash rather than in interest-bearing liquid assets).

Formally, the precautionary money demand $\ell^P$ can be derived from a stochastic extension of the inventory model in which income arrives randomly rather than at fixed intervals. Under uncertainty about the timing of expenditures, the optimal cash holding exceeds the Baumol–Tobin deterministic level, and the precautionary component increases with the variance of expenditure timing.

### The Speculative (Portfolio) Motive

Tobin (1958) argued that money is held as part of an optimal financial portfolio, not merely for transactions. Even if a household has no immediate transactions needs, it may prefer to hold some money rather than risky bonds if bonds carry capital loss risk. The **speculative demand** reflects the portfolio-theoretic tradeoff between the return on bonds and the safety of money.

In Tobin's mean-variance framework, a household allocates wealth $W$ between money (zero nominal return, zero risk) and a risky bond with expected nominal return $\mu_b$ and standard deviation $\sigma_b$. With bond share $\alpha$, portfolio mean return $\mu_p = \alpha\mu_b$ and standard deviation $\sigma_p = \alpha\sigma_b$. Maximizing mean-variance utility $U = \mu_p - (\lambda/2)\sigma_p^2$:

$$\alpha^* = \frac{\mu_b}{\lambda\sigma_b^2}, \quad \ell^S = (1-\alpha^*)W = W\!\left(1 - \frac{\mu_b}{\lambda\sigma_b^2}\right).$$

Speculative money demand increases with risk aversion $\lambda$ and bond return variance $\sigma_b^2$, and decreases with the expected bond return $\mu_b$. This is consistent with the flight-to-safety observed during financial crises: when bond market volatility spikes ($\sigma_b^2$ rises), households and institutions shift out of risky assets into cash and government money-market instruments, irrespective of the level of interest rates.

### The Aggregate Money Demand Function

Combining all three motives, the aggregate real money demand function $\ell(Y, i)$ satisfies:

$$\ell(Y, i) = \ell^T(Y, i) + \ell^P(Y, \sigma^2) + \ell^S(W, i, \sigma_b^2),$$

with $\partial\ell/\partial Y > 0$ (income elasticity positive), $\partial\ell/\partial i < 0$ (interest semi-elasticity negative), and $\partial\ell/\partial\sigma^2 > 0$ (uncertainty raises money demand). For empirical estimation, it is common to use the log-linear specification:

$$\ln\ell_t = \eta_Y\ln Y_t + \eta_i\ln i_t + \text{const} + \epsilon_t,$$

where $\eta_Y > 0$ is the income elasticity and $\eta_i < 0$ the interest elasticity. Early empirical estimates (Meltzer, 1963) found $\eta_Y \approx 1$ and $\eta_i \approx -0.5$. Post-1980 financial innovation destabilized these relationships: as money market mutual funds and sweep accounts proliferated, households could economize on non-interest-bearing M1 without any change in income or the interest rate, shifting the demand function unpredictably. This instability is one reason that modern central banks target interest rates directly (via the Taylor rule) rather than targeting the money supply.

---

## 18.2 Money-Market Equilibrium and the LM Curve

In the standard framework, money-market equilibrium requires real money demand to equal the real money supply, controlled by the central bank through the monetary base $H$:

$$\frac{M_t}{P_t} = \ell(Y_t, i_t),$$

where $M_t = m\cdot H_t$ is broad money (M1 or M2) determined by the money multiplier $m = (1+c_r)/(c_r + rr + er)$ as in Chapter 14. For a given real money supply $M_t/P_t$, this condition defines the LM curve: the set of $(Y_t, i_t)$ pairs at which the money market clears.

**Definition (LM Curve).** The **LM curve** is the locus of income $Y$ and nominal interest rate $i$ combinations consistent with money-market equilibrium at a given real money supply $M/P$:

$$\text{LM:} \quad \frac{M}{P} = \ell(Y, i), \quad \text{i.e.,} \quad i = \ell^{-1}\!\left(Y, \frac{M}{P}\right),$$

where $\ell^{-1}$ denotes the inverse of $\ell$ with respect to $i$. The LM curve is upward sloping: at higher $Y$, money demand rises; to restore equilibrium at a fixed $M/P$, the interest rate must rise to reduce money demand to its supply.

Using the linear demand specification $\ell(Y,i) = kY - hi$ with $k,h > 0$:

$$i^* = \frac{kY - M/P}{h}.$$

The slope $\partial i/\partial Y|_{LM} = k/h > 0$ is steeper when money demand is less interest-sensitive (small $h$) and more income-sensitive (large $k$). An open market purchase that increases $M/P$ shifts the LM curve downward (reduces $i^*$ at any $Y$) by $\Delta i = -\Delta(M/P)/h$.

### The Interest Rate Channel of Monetary Transmission

The LM curve captures the first step in the monetary transmission mechanism: a change in the money supply changes the interest rate, which then propagates to the real economy through the IS relationship. The full IS–LM transmission chain:

$$\Delta M \to \Delta(M/P) \to \Delta i \to \Delta I \to \Delta Y.$$

Quantitatively, the monetary multiplier from Chapter 9:

$$\frac{\Delta Y}{\Delta(M/P)} = \frac{b_r}{h + b_r k},$$

where $b_r$ is investment interest sensitivity. The transmission is stronger when investment responds sharply to interest rates (large $b_r$) and weaker when money demand is highly interest-elastic (large $h$). In the liquidity trap, $h\to\infty$ and the monetary multiplier falls to zero — the subject of Section 18.3.

---

## 18.3 The Liquidity Trap: Mechanism, History, and Policy Implications

The liquidity trap is perhaps the most practically important limiting case in monetary economics. It describes a situation in which conventional monetary policy — reducing the short-term interest rate — loses its ability to stimulate the economy, trapping it at below-potential output even though the central bank is willing to print money without limit.

**Definition (Liquidity Trap).** The economy is in a **liquidity trap** when: (i) the nominal interest rate is at or near its effective lower bound ($i_t \approx i^{ELB} \approx 0$); (ii) money demand is perfectly interest-elastic at this rate — the public is willing to hold any additional money supply at the prevailing rate without the interest rate falling further; and (iii) further conventional monetary expansion (increases in $M$) fails to reduce $i$ or stimulate output. The LM curve is effectively horizontal at $i = i^{ELB}$.

The mechanism behind the liquidity trap is the **portfolio preference for liquidity over bonds at low interest rates**. When the interest rate is very close to zero, holding money rather than bonds has virtually no opportunity cost. If agents also expect bond prices to fall (interest rates to rise) in the future, the risk-adjusted return on bonds may be *below* zero, making money strictly preferable. At this point, any additional money injected into the economy is willingly absorbed into idle money holdings — the "idle balances" of Keynes's terminology — without any reduction in interest rates.

Formally, define $i^{ELB}$ as the effective lower bound (approximately zero, but potentially slightly negative if storage costs make currency holding costly). When $i \leq i^{ELB}$, the market-clearing condition becomes:

$$\frac{M_t}{P_t} = kY_t - h\cdot i^{ELB} + H_t,$$

where $H_t \geq 0$ is the quantity of excess or "hoarded" money balances held at the lower bound. The central bank can increase $M_t/P_t$ indefinitely, but this only increases $H_t$ — the idle cash pool — without lowering $i$ below $i^{ELB}$ or increasing $Y$. The LM curve has become infinitely elastic at $i^{ELB}$.

### Historical Episodes

The liquidity trap was not merely a theoretical curiosity for Keynes: it was motivated by the Great Depression experience (1929–33), in which nominal interest rates reached near-zero levels while the economy remained severely depressed. Friedman and Schwartz (1963) argued that the U.S. Fed could have broken the trap by sufficiently aggressive money supply increases, but Keynes's view was that — at the lower bound — fiscal policy was required.

The next major episode was Japan in the 1990s and 2000s. After the collapse of the asset price bubble in 1989–90, the Bank of Japan progressively reduced rates toward zero; by 1999, the uncollateralized overnight call rate was essentially at zero and the Japanese economy remained sluggish despite large fiscal expansions. The Bank of Japan introduced **quantitative easing** (QE) in 2001 — purchases of government bonds to expand bank reserves — in an attempt to escape the trap; the results were modest, and the experience confirmed that QE does not automatically generate inflation or strong output growth when expectations are depressed.

The 2008–09 Global Financial Crisis brought the liquidity trap to advanced economies more broadly: the U.S. federal funds rate reached 0–0.25% in December 2008, the ECB's deposit rate reached zero in 2012 and went negative in 2014, and the Bank of England reached 0.1% in 2020. All major central banks deployed unconventional tools in response, including QE, forward guidance, and negative policy rates.

### The Paradox of Thrift

An important macroeconomic implication of the liquidity trap is the **paradox of thrift**: an attempt by households to save more, which would normally raise the capital stock and future income in the Solow model, actually reduces current output and may leave saving unchanged or lower in a demand-constrained economy at the lower bound.

To see this, suppose households wish to reduce consumption and increase saving when $i = i^{ELB}$. The reduction in $C$ shifts the IS curve left — at any given $i$, equilibrium $Y$ is lower. But since $i$ cannot fall below $i^{ELB}$, the fall in $C$ is not offset by rising investment; instead, output falls. The fall in output reduces income, which reduces saving as well (since $S = Y - C$). Society's attempt to save more results in no increase in aggregate saving — a logical impossibility in the standard flexible-price model but a coherent equilibrium outcome in the demand-constrained Keynesian model.

Formally: in goods-market equilibrium $Y = C(Y-T) + I(r) + G$. With $r$ fixed at $r^{ELB}$, differentiate with respect to autonomous consumption $a$:

$$\frac{\mathrm{d}Y^*}{\mathrm{d}a} = \frac{1}{1-C'} = \kappa_G > 0.$$

A decline in $a$ (increased thrift) reduces $Y^*$ by a multiple of the decline. Aggregate saving $S = Y^* - C(Y^*-T) - G$: differentiating, $\mathrm{d}S/\mathrm{d}a = (1 - C')\mathrm{d}Y^*/\mathrm{d}a - 1 = (1-C')\kappa_G - 1 = 0$. Aggregate saving is unchanged despite the increased desire to save — the multiplier mechanism ensures that reduced income exactly offsets the increased propensity to save.

---

## 18.4 The Interbank Market and Modern Central Bank Operations

In the modern institutional setting, the central bank does not directly set the money supply; it targets the overnight interest rate at which banks lend reserves to one another. In the United States, this is the **federal funds rate** — the rate at which banks borrow and lend excess reserves overnight in the federal funds market.

**Definition (Federal Funds Market).** The **federal funds market** is the overnight interbank market in which commercial banks lend and borrow reserve balances held at the Federal Reserve. A bank with excess reserves (more than required to meet reserve requirements) lends to a bank with a reserve shortfall at the federal funds rate $i_{FF}$. The Fed's Open Market Desk controls reserve supply to keep $i_{FF}$ near its target.

The federal funds market equilibrium:

$$R^D(i_{FF}) = R^S,$$

where $R^D$ is aggregate reserve demand (decreasing in $i_{FF}$: higher rates make holding excess reserves more costly) and $R^S$ is the Fed-controlled reserve supply. In the **corridor system** used before 2008 (scarce reserves), $R^S$ was tightly controlled and small changes in reserve supply generated precise changes in $i_{FF}$.

After 2008, the Fed shifted to an **ample-reserves framework** (also called a floor system): it pays **interest on reserve balances (IORB)** at a rate set by policy, and maintains a large surplus of reserves in the system. In this framework:

$$i_{FF} \geq i^{IORB},$$

because no bank will lend reserves in the federal funds market at a rate below what the Fed pays to hold those same reserves. The IORB rate therefore serves as an effective floor for $i_{FF}$. The Fed also conducts **overnight reverse repurchase agreements (ON RRP)** at a rate $i^{ONRRP} \leq i^{IORB}$, which provides a floor accessible to non-bank financial institutions. Together, IORB and ON RRP create a corridor system for federal funds rates even with abundant reserves.

### The Transmission from Policy Rate to Market Rates

The federal funds rate affects broader financial conditions through several channels. The **expectations hypothesis of the term structure** links the long-term nominal interest rate to expected future short rates:

$$i_t^n = \frac{1}{n}\sum_{k=0}^{n-1}\mathbb{E}_t[i_{t+k}^1] + \phi_t^n,$$

where $\phi_t^n$ is the term premium (compensation for duration risk). A commitment to keep the federal funds rate low for an extended period — forward guidance — reduces long-term rates by lowering the expected future short rates, even without any current change in policy. This is the mechanism behind the "lower for longer" guidance deployed by the Fed after 2008 and again in 2020.

The **credit channel** operates in addition to the interest rate channel. When the central bank raises rates, bank profitability may be affected (especially if banks have floating-rate assets and fixed-rate liabilities), reducing their capacity to lend. The **bank lending channel** — the Bernanke–Blinder (1988) contribution — shows that monetary tightening reduces bank reserves, constrains bank lending (not just its price), and is particularly powerful for bank-dependent borrowers (small firms, households without access to capital markets).

---

## 18.5 Negative Interest Rates and the Reversal Rate

Several central banks have experimented with negative policy rates: the ECB deposit facility rate reached $-0.5\%$ in 2019, the Bank of Japan's policy rate balance was set to $-0.1\%$ in 2016, and the Danish, Swiss, and Swedish central banks pushed rates significantly below zero. This experience has raised the question of whether there is a lower bound below zero below which further cuts become counterproductive.

**Definition (Reversal Rate).** The **reversal rate** (Brunnermeier and Koby, 2019) is the interest rate below which further cuts reduce rather than increase bank lending, because the negative effect on bank net interest margins (and hence bank profitability and equity) outweighs the stimulative effect on borrower demand. Below the reversal rate, monetary easing tightens rather than loosens credit conditions, reversing the conventional monetary transmission mechanism.

The intuition: banks fund themselves partly with deposits whose rates cannot go negative (because depositors would switch to cash), while bank assets have rates that fall with the policy rate. As rates go increasingly negative, the net interest margin (lending rate minus deposit rate) narrows, compressing bank profits and equity. When bank equity falls, banks face tighter regulatory constraints on lending (minimum capital ratios) and reduce credit supply — exactly the opposite of the intended monetary stimulus.

The reversal rate is not the same as the ELB: the ELB is a technological lower bound set by cash storage costs; the reversal rate is an economic lower bound set by bank profitability and the institutional structure of deposit banking. Empirical estimates suggest the reversal rate is approximately $-0.5\%$ to $-1\%$ for major banking systems (Ulate, 2021).

---

*Next: Chapter 19 — The Labor Market*
