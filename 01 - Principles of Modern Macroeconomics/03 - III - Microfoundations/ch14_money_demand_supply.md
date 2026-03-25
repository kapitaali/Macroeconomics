# Chapter 14 — Money Demand and Supply: The Role of Central Banks

> *"Inflation is always and everywhere a monetary phenomenon."*
> — Milton Friedman, *A Program for Monetary Stability*, 1960

---

Money occupies a unique position in macroeconomics. Unlike other assets — stocks, bonds, real estate — money is held not primarily for the return it generates (typically zero in nominal terms, negative in real terms during inflation) but for the services it provides as a medium of exchange. The demand for money is therefore the demand for transaction services, and its relationship to income and interest rates is the crucial link between the monetary and real sides of the economy. Understanding how money is created, how much of it households and firms choose to hold, and how the central bank controls its supply is prerequisite to analyzing monetary policy and inflation.

---

## 14.1 The Functions of Money and the Definition of Monetary Aggregates

**Definition (Money).** **Money** is any asset that is generally accepted as a medium of exchange in transactions. It performs three functions: (i) **medium of exchange** — it eliminates the double coincidence of wants problem of barter; (ii) **unit of account** — prices are quoted in terms of it; and (iii) **store of value** — it transfers purchasing power across time, though less efficiently than interest-bearing assets.

The boundary between money and non-money is not sharp. Currency is unambiguously money; a Treasury bill is unambiguously not money. Demand deposits, money market funds, and short-term commercial paper lie on a spectrum. Statistical agencies construct multiple monetary aggregates to capture different points along this spectrum:

$$M0 \equiv \text{Currency} + \text{Bank reserves at central bank (monetary base)}$$
$$M1 \equiv \text{Currency} + \text{Demand deposits} + \text{Checkable deposits}$$
$$M2 \equiv M1 + \text{Savings deposits} + \text{Small time deposits} + \text{Retail money market funds}$$

The **money multiplier** relates the broad money supply to the monetary base. If $c_r = C/D$ is the currency-to-deposit ratio, $rr$ is the required reserve ratio, and $er$ is the excess reserve ratio:

$$M1 = \frac{1 + c_r}{c_r + rr + er}\, H \equiv m\cdot H,$$

where $H$ is the monetary base. The multiplier $m > 1$ when $rr + er < 1$: each dollar of base money supports multiple dollars of deposits through fractional-reserve banking. Following 2008, when the Federal Reserve began paying interest on excess reserves (IOER), banks accumulated enormous quantities of excess reserves ($er$ rose dramatically), and the money multiplier fell sharply — illustrating that $m$ is not a stable structural parameter but an endogenous outcome depending on banks' incentives to hold reserves.

---

## 14.2 Baumol–Tobin: The Inventory-Theoretic Demand for Money

The first rigorous theory of money demand treats it as an inventory management problem. Why do households hold money when they could hold interest-bearing bonds? Because converting bonds to money involves a transaction cost, just as converting warehouse inventory to retail floor space involves a physical cost.

A household receives income $Y$ at the start of the period and spends it evenly throughout. It can hold it as money (yielding zero interest) or as bonds (yielding nominal rate $i$). Transferring between money and bonds — a trip to the bank or a financial intermediary — costs $b$ per trip (time and bother). If the household makes $n$ trips per period, withdrawing $Y/n$ each time, average money holdings are $Y/(2n)$.

Total cost: opportunity cost of holding money plus transaction cost:

$$TC(n) = i\cdot\frac{Y}{2n} + b\cdot n.$$

Minimizing over $n$: $n^* = \sqrt{iY/(2b)}$. Substituting back, optimal average money demand is:

$$L^{BT} = \frac{Y}{2n^*} = \sqrt{\frac{bY}{2i}}.$$

**Definition (Baumol–Tobin Money Demand).** The **Baumol–Tobin money demand function** $L^{BT} = \sqrt{bY/(2i)}$ implies: (i) an income elasticity of $+1/2$ — doubling income less than doubles money demand; (ii) an interest elasticity of $-1/2$ — doubling the interest rate reduces money demand by about 30%. These predictions — income elasticity less than one and negative interest semi-elasticity — are broadly consistent with empirical estimates of long-run money demand.

---

## 14.3 Seigniorage and the Inflation Tax

A government that finances its deficit by printing money — by having the central bank create reserves to purchase government debt — earns **seigniorage**: the revenue that accrues from the monopoly right to issue currency.

**Definition (Seigniorage).** **Seigniorage** is the real revenue earned by the government from money creation:

$$S = \frac{\dot{M}}{P} = \frac{\dot{M}}{M}\cdot\frac{M}{P} = \hat{m}\cdot\frac{M}{P},$$

where $\hat{m} = \dot{M}/M$ is the money growth rate and $M/P$ is the real money supply. In steady state, when $\hat{m} = \pi$ (money growth equals inflation), seigniorage equals $\pi \cdot (M/P)$.

The economic interpretation is the **inflation tax**: households who hold money lose real purchasing power at rate $\pi$ per period, effectively transferring resources to the government. The real money balance $M/P$ is the **tax base** and $\pi$ is the **tax rate**. Like any tax, there is a Laffer curve: raising $\pi$ raises revenue at low inflation rates (higher tax rate on a given base) but reduces the base at high rates (households economize on money holdings). The revenue-maximizing inflation rate is $\pi^* = -1/\varepsilon_i$, where $\varepsilon_i < 0$ is the interest semi-elasticity of money demand.

Cagan (1956) applied this model to seven hyperinflationary episodes, fitting the money demand function $\ln(M/P) = \alpha - \beta\pi^e$ and showing that governments in hyperinflation consistently operated on the right side of the Laffer curve — collecting ever-shrinking seigniorage revenue as money demand collapsed — consistent with a fiscal system that had entirely lost access to conventional tax revenue and was printing money as a last resort.

---

## 14.4 The Money Supply Process and Central Bank Operations

**Definition (Open Market Operation).** An **open market operation (OMO)** is the purchase or sale of government securities by the central bank in exchange for bank reserves. An open market purchase injects reserves into the banking system, expanding the monetary base $H$; an open market sale withdraws reserves, contracting $H$.

Under conventional monetary policy, the central bank sets the overnight interest rate — the federal funds rate in the United States — by adjusting reserve supply to meet reserve demand at the target rate. The money market equilibrium:

$$R^D(i_{FF}) = R^S,$$

where $R^D$ is aggregate reserve demand (decreasing in the funds rate, as higher rates make reserve holdings more expensive) and $R^S$ is the supply of reserves controlled by the central bank. Before 2008, reserve supply was scarce and the funds rate was above the IOER rate (which was zero); after 2008, reserve supply was abundant and the IOER rate became the effective floor.

**Definition (Effective Lower Bound).** The **effective lower bound (ELB)** on nominal interest rates is the interest rate below which further cuts cease to stimulate the economy, because holding currency (which yields zero nominal return) becomes attractive relative to lending at negative rates. Physical currency storage costs set the ELB below zero; empirical experience in Europe and Japan suggests the ELB is approximately $-0.5\%$ to $-1\%$. At the ELB, conventional monetary policy — cutting the policy rate — is exhausted, motivating the unconventional tools discussed in Chapter 29.

---

*Next: Chapter 15 — Expectations and Behavioral Macroeconomics*
