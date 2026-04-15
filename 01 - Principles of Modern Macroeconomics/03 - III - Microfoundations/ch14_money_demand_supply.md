# Chapter 14 — Money Demand and Supply: The Role of Central Banks

> *"Inflation is always and everywhere a monetary phenomenon."*
> — Milton Friedman, *A Program for Monetary Stability*, 1960

---

Money occupies a singular position in the architecture of macroeconomics. Unlike other assets — stocks, bonds, real estate — money is held not primarily for the return it generates (typically zero in nominal terms, negative in real terms during inflation) but for the services it provides: as a medium of exchange it eliminates the double-coincidence-of-wants problem of barter; as a unit of account it coordinates price-setting across a decentralized economy; as a store of value it bridges transactions over time. The demand for money is therefore the demand for these transaction services, and its relationship to income and interest rates is the crucial link between the monetary and real sides of the economy. This chapter builds that link: from the definition and measurement of money, through the inventory-theoretic theory of money demand, to the fiscal arithmetic of seigniorage, and finally to how modern central banks actually implement monetary policy — an account that differs substantially from the textbook money multiplier story.

---

## 14.1 The Functions of Money and the Definition of Monetary Aggregates

**Definition (Money).** **Money** is any asset that is generally accepted as a medium of exchange in transactions. It performs three functions: (i) **medium of exchange** — it eliminates the requirement for a double coincidence of wants; (ii) **unit of account** — prices, contracts, and debts are denominated in it; (iii) **store of value** — it transfers purchasing power across time, though less efficiently than interest-bearing alternatives.

The boundary between money and non-money is not sharp. Currency is unambiguously money; a long-term government bond is unambiguously not money. Demand deposits, money market funds, and short-term commercial paper lie along a spectrum of liquidity. Statistical agencies construct multiple monetary aggregates to capture different points along this spectrum:

$$M0 \equiv \text{Currency in circulation} + \text{Bank reserves at the central bank}$$
$$M1 \equiv \text{Currency} + \text{Demand deposits} + \text{Other checkable deposits}$$
$$M2 \equiv M1 + \text{Savings deposits} + \text{Small time deposits} + \text{Retail money market funds}$$

The choice of aggregate matters for empirical analysis. M1 is the most liquid measure; M2 better captures the full range of assets households use to finance transactions over somewhat longer horizons. For most modern central banks, however, neither aggregate is directly targeted: the collapse of stable money demand relationships in the 1980s and 1990s — as financial innovation altered the velocity of M1 and M2 — led central banks to abandon money targeting in favor of interest rate rules [Ch. 23].

### The Money Multiplier

The **money multiplier** relates the broad money supply to the monetary base $H$ (currency plus bank reserves). If $c_r = C/D$ is the currency-to-deposit ratio, $rr$ is the required reserve ratio, and $er$ is the excess reserve ratio:

$$M1 = \frac{1 + c_r}{c_r + rr + er}\, H \equiv m \cdot H.$$

The multiplier $m > 1$ when $rr + er < 1$: each dollar of base money supports multiple dollars of deposits through fractional-reserve banking, as each new deposit enables new loans that create new deposits. Following 2008, when the Federal Reserve began paying interest on excess reserves (IOER), commercial banks accumulated enormous quantities of excess reserves ($er$ rose dramatically toward 10–15% of deposits), and the money multiplier fell sharply. This episode illustrates that $m$ is not a stable structural parameter but an endogenous outcome that depends on banks' incentives to hold reserves — a point the textbook deposit expansion story obscures.

---

## 14.2 Baumol–Tobin: The Inventory-Theoretic Demand for Money

The first rigorous theory of money demand treats it as an inventory management problem. The question is precisely posed: a household receives income $Y$ at the start of the period and spends it evenly. It can hold wealth as money (yielding zero nominal interest) or as bonds (yielding nominal rate $i$). Transferring between bonds and money — a trip to the bank, a call to a broker — costs a fixed amount $b$ per transaction (representing time, effort, and any fees). How many transfers should the household make?

If the household makes $n$ transfers per period, it withdraws $Y/n$ each time and holds an average money balance of $Y/(2n)$ between withdrawals. Total cost per period is the opportunity cost of holding money plus the transaction cost:

$$TC(n) = i \cdot \frac{Y}{2n} + b \cdot n.$$

Minimizing over $n$: the first-order condition $\partial TC/\partial n = -iY/(2n^2) + b = 0$ gives $n^* = \sqrt{iY/(2b)}$. Substituting back:

$$L^{BT} = \frac{Y}{2n^*} = \sqrt{\frac{bY}{2i}}.$$

**Definition (Baumol–Tobin Money Demand).** The **Baumol–Tobin money demand function** $L^{BT} = \sqrt{bY/(2i)}$ predicts: (i) an income elasticity of $+1/2$ — doubling income less than doubles desired money balances, consistent with economies of scale in cash management; (ii) an interest elasticity of $-1/2$ — doubling the interest rate reduces money demand by about 30%. These predictions are broadly consistent with empirical estimates of long-run money demand and stand in contrast to the earlier quantity-theoretic assumption of proportionality between money and income.

The income elasticity below unity has an important implication: as economies grow richer and transactions become larger, the demand for money as a share of income declines — a prediction consistent with the secular decline in currency-to-GDP ratios observed in high-income countries over the twentieth century. The interest elasticity, meanwhile, implies that low interest rate environments (such as Japan since the 1990s or the global economy after 2008) will be associated with elevated money demand — a precondition for the liquidity trap discussed in Chapter 18.

### Portfolio Demand for Money: Tobin's Extension

Baumol's model focuses on transactions demand. James Tobin (1958) showed that a **portfolio motive** for holding money also exists: even a risk-averse investor who holds bonds for yield will want to hold some money as insurance against the variance of bond returns. In Tobin's mean-variance framework, money is the zero-variance asset, and its share in the optimal portfolio depends on risk aversion and the risk-return tradeoff offered by risky assets. This portfolio motive for money demand, rather than the transactions motive, is the more relevant consideration for the large institutional holders — pension funds, corporations — whose treasury management decisions explain a substantial part of aggregate M2 demand.

---

## 14.3 Seigniorage and the Inflation Tax

A government unable or unwilling to finance its deficit through taxes or bond sales may resort to money creation — instructing the central bank to purchase government debt with newly printed reserves. The revenue earned from this monopoly power over currency issuance is called **seigniorage**.

**Definition (Seigniorage).** **Seigniorage** is the real revenue the government earns from money creation:

$$S = \frac{\dot{M}}{P} = \frac{\dot{M}}{M} \cdot \frac{M}{P} = \hat{m} \cdot \frac{M}{P},$$

where $\hat{m} = \dot{M}/M$ is the money growth rate and $M/P$ is the real money supply. In steady state where money growth equals inflation ($\hat{m} = \pi$):

$$S = \pi \cdot \frac{M}{P}.$$

The interpretation is the **inflation tax**: holders of money balances lose real purchasing power at rate $\pi$ per period, effectively transferring resources to whoever printed the money — typically the government. The real money stock $M/P$ is the tax base, and $\pi$ is the tax rate.

Like any tax, the seigniorage revenue function has a **Laffer curve** shape. At low inflation, raising $\pi$ increases revenue because the tax rate rises faster than the base shrinks. At high inflation, the base — real money demand — collapses (households economize on money holdings through costly avoidance strategies such as barter, foreign currency, or rapid circulation), and further increases in $\pi$ reduce revenue. The **revenue-maximizing inflation rate** is:

$$\pi^{max} = \frac{1}{|\varepsilon_i|},$$

where $\varepsilon_i < 0$ is the interest semi-elasticity of money demand (the percentage change in real money demand per percentage point change in the nominal interest rate). For typical estimates of $|\varepsilon_i| \approx 0.05$–$0.1$, the revenue-maximizing inflation rate lies between 10% and 20% annually — well above the targets of modern inflation-targeting central banks.

### Hyperinflation and Fiscal Collapse

Cagan (1956) applied the seigniorage model to seven postwar hyperinflationary episodes — Germany 1922–23, Hungary 1945–46, Greece 1941–44, Poland 1923–24, Russia 1921–24, Austria 1921–22, and Hungary 1923–24 — each involving monthly inflation rates exceeding 50%. He fitted the money demand function $\ln(M/P)_t = \alpha - \beta\pi^e_t$ with adaptive expectations and showed that governments in hyperinflation consistently operated on the right-hand (downward-sloping) side of the Laffer curve: as inflation accelerated, real money demand collapsed faster than nominal money growth increased, causing seigniorage revenue to shrink even as inflation soared. This is consistent with a fiscal system that had lost access to conventional tax revenue and bond markets — one printing money as a desperate last resort while that last resort steadily destroyed itself.

The fiscal theory of price level determination [Ch. 22] provides the deeper account: hyperinflation is ultimately a fiscal phenomenon, reflecting the market's judgment that the intertemporal budget constraint of the government is unsustainable. The correct stabilization of hyperinflation always involves a credible fiscal reform alongside the monetary reform — Sargent (1982) documents this pattern in each of Cagan's four European cases.

---

## 14.4 The Quantity Theory and Its Modern Reformulation

The quantity theory of money, rooted in the works of Hume and Fisher and formalized by Friedman, asserts a proportional long-run relationship between the money supply and the price level.

**Definition (Quantity Equation).** The **quantity equation** is the accounting identity:

$$MV = PY,$$

where $M$ is the money supply, $V$ is the **velocity of circulation** (the average number of times a unit of money is used in transactions per period), $P$ is the price level, and $Y$ is real output. In growth rate form:

$$\hat{m} + \hat{v} = \pi + \hat{y}.$$

As an identity, the quantity equation is always true — it defines velocity. The quantity **theory** adds the behavioral assertion that $V$ is approximately stable (or at least predictable), so that movements in $M$ translate approximately one-for-one into movements in $P$ (or $PY$ in the short run). Friedman's monetarist position — that "inflation is always and everywhere a monetary phenomenon" — rests on this stability claim.

The empirical record is mixed. Over long horizons and across countries, the correlation between money growth and inflation is strong and close to one-for-one — Friedman's proposition holds in the long run. Over shorter horizons, however, velocity is volatile, particularly during periods of financial innovation (the 1980s) or financial stress (2008–2020), when broad money demand shifted substantially. This instability is why modern central banks have largely abandoned money supply targets in favor of inflation-forecast targeting with an interest rate instrument.

---

## 14.5 The Money Supply Process and Central Bank Operations

**Definition (Open Market Operation).** An **open market operation (OMO)** is the purchase or sale of government securities by the central bank in exchange for bank reserves. An open market purchase injects reserves into the banking system, expanding the monetary base $H$; an open market sale withdraws reserves, contracting $H$.

The textbook account has the central bank controlling the money supply by controlling $H$, with the multiplier $m$ mechanically expanding this into $M1$. Modern central banking works differently. Central banks set a **target interest rate** — the federal funds rate in the United States, the repo rate in the UK — and supply whatever reserves are needed to maintain the overnight interbank rate at the target. The money supply is demand-determined: banks create deposits by extending loans, and the central bank accommodates the resulting demand for reserves at the target rate.

**Definition (Interest Rate Corridor).** Modern central banks typically operate a **corridor system**: the policy rate is bounded below by the rate paid on excess reserves (IOER or deposit facility rate) and above by the rate charged for borrowing reserves from the central bank (the discount rate or marginal lending facility rate). By adjusting the corridor, the central bank steers the overnight rate without necessarily specifying a money supply target.

**Definition (Effective Lower Bound).** The **effective lower bound (ELB)** on nominal interest rates is the rate below which further cuts cease to ease financial conditions as intended, because holding physical currency — which yields zero nominal return — becomes attractive relative to deposits earning negative rates. Physical currency storage costs set a floor somewhat below zero; empirical experience across Europe and Japan places the ELB at approximately $-0.5\%$ to $-1.0\%$. At the ELB, conventional monetary policy — cutting the policy rate — is exhausted. This motivates the unconventional tools discussed in Chapter 29: quantitative easing, forward guidance, yield curve control, and, at the limit, proposals for abolishing physical currency to eliminate the ELB altogether.

### The Transmission of Monetary Policy

The effect of monetary policy on the real economy operates through several channels, each with different lags and magnitudes.

The **interest rate channel** is the most direct: a cut in the policy rate reduces borrowing costs, stimulating investment [Ch. 12] and durable goods consumption [Ch. 11]. The **asset price channel** works through the discount rate effect on equity and real estate valuations; rising asset prices increase household wealth, stimulating consumption (Chapter 11's wealth effect). The **exchange rate channel** — a lower domestic interest rate depreciates the currency, improving net exports [Ch. 21] — is particularly important for small open economies. The **credit channel** (Bernanke and Gertler, 1989) works through the financial accelerator: lower policy rates ease collateral constraints on firms and households, amplifying the direct interest rate effect through improved borrowing capacity [Ch. 24]. Finally, the **expectations channel** — the direct effect of central bank communications on beliefs about future rates and inflation — has gained importance as central banks have adopted more transparent, forward-looking communication frameworks [Ch. 29].

---

## Chapter Summary

- **Money** performs three functions — medium of exchange, unit of account, store of value — and is measured across a spectrum of liquidity (M0, M1, M2). The money multiplier $m = (1+c_r)/(c_r+rr+er)$ is endogenous, not structural.

- The **Baumol–Tobin model** derives money demand $L^{BT} = \sqrt{bY/(2i)}$ from inventory optimization, predicting income elasticity $+1/2$ and interest elasticity $-1/2$, consistent with empirical evidence and with economies of scale in cash management.

- **Seigniorage** $S = \pi \cdot M/P$ is the inflation tax; its Laffer curve shape implies a revenue-maximizing inflation rate $\pi^{max} = 1/|\varepsilon_i|$. Hyperinflations represent the collapse of seigniorage revenue as real money demand implodes, and are fundamentally fiscal phenomena requiring fiscal reform to resolve.

- The **quantity equation** $MV = PY$ is an identity; the quantity theory adds the behavioral assumption of stable velocity. The long-run proportionality of money growth and inflation is well established; short-run velocity instability explains why central banks abandoned money targeting.

- Modern central banks implement policy by targeting the overnight interest rate through a corridor system, supplying reserves on demand. The **effective lower bound** at approximately $-0.5\%$ to $-1\%$ marks the limit of conventional policy, motivating unconventional tools.

- Monetary policy reaches the real economy through multiple channels: the interest rate, asset prices, the exchange rate, the credit accelerator, and expectations. Their relative importance varies across countries, time periods, and the financial system's structure.

---

*Next: Chapter 15 — Expectations and Behavioral Macroeconomics*
