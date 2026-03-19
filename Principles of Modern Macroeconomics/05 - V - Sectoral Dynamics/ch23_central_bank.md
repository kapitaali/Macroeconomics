# Chapter 23 — The Central Bank and Monetary Policy: Tools and Targets

> *"Price stability is not an end in itself, but a precondition for achieving sustainable growth and employment."*
> — Mario Draghi, President of the ECB, 2011–2019

---

The central bank is the institution at the heart of modern monetary policy. It controls the overnight interest rate, issues currency, acts as lender of last resort to the banking system, and — increasingly — deploys unconventional tools to influence long-term interest rates and credit conditions when the policy rate is at its effective lower bound. Understanding how the central bank operates, what it is trying to achieve, and how its decisions interact with the economy's forward-looking dynamics is indispensable for applied macroeconomics.

This chapter develops the theory of central banking in three parts. First, the objectives of monetary policy and the institutional design of modern inflation-targeting central banks. Second, optimal monetary policy in the New Keynesian framework — the theoretical benchmark against which actual policy is evaluated. Third, the tools and challenges of monetary policy at the effective lower bound, where the standard interest rate instrument is exhausted and central banks must rely on unconventional measures.

---

## 23.1 Monetary Policy Frameworks and the Taylor Rule

Modern central banks operate under **inflation targeting** frameworks, in which the central bank commits to keeping CPI inflation near a quantitative target (typically 2%) over a medium-term horizon of 18–24 months. The framework has four core features: a quantitative inflation target; operational (instrument) independence from government; transparency through regular inflation reports and press conferences; and accountability mechanisms (parliamentary hearings, published minutes, open letters when the target is missed).

Inflation targeting was introduced by the Reserve Bank of New Zealand in 1990 and has since been adopted by over 30 central banks, including the Bank of England (1992), the Bank of Canada (1991), and the Bank of Sweden (1993). The European Central Bank and the U.S. Federal Reserve use variants that emphasize both inflation and employment mandates.

### The Taylor Rule

In a seminal empirical contribution, Taylor (1993) showed that U.S. monetary policy during 1987–1992 was well described by a simple linear reaction function:

$$i_t = r^n + \pi^* + \phi_\pi(\pi_t - \pi^*) + \phi_y\hat{x}_t,$$

with $r^n = 2\%$, $\pi^* = 2\%$, $\phi_\pi = 1.5$, $\phi_y = 0.5$. The Taylor rule has two important features. First, it incorporates the **Taylor principle**: $\phi_\pi > 1$, so that when inflation rises by 1 percentage point, the nominal rate rises by more than 1 pp, causing the real rate to rise. This is the key stabilizing property: if the real rate did not rise with inflation, higher inflation would reduce real rates, stimulate demand, and generate a self-reinforcing inflationary spiral. Second, it responds to the output gap $\hat{x}_t$: when output is above potential, the central bank tightens preemptively to prevent inflationary overheating.

A more general specification incorporates **interest rate smoothing** — the empirically robust observation that central banks move rates gradually rather than in large discrete steps:

$$i_t = \rho_i\, i_{t-1} + (1-\rho_i)\bigl[r^n + \pi^* + \phi_\pi(\pi_t - \pi^*) + \phi_y\hat{x}_t\bigr] + \epsilon_t,$$

where $\rho_i \approx 0.85$ for the U.S. Fed, indicating that 85% of the desired rate change is carried over from the previous period. Smoothing is rationalized by: (i) uncertainty about the output gap (measured in real time, it is highly unreliable); (ii) the desire to avoid sending a misleadingly aggressive signal to markets; (iii) the possibility that smooth, predictable moves maintain credibility better than volatile ones.

### The Taylor Principle and Determinacy

The Taylor principle $\phi_\pi > 1$ is not merely a desirable property — it is necessary for the New Keynesian model to have a unique bounded rational expectations equilibrium. To see this, consider a simplified NK model with the IS curve and NKPC:

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\pi_{t+1}] - r^n)$$
$$\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t.$$

Substituting the Taylor rule $i_t = r^n + \phi_\pi\pi_t + \phi_y\hat{x}_t$ (ignoring smoothing and the inflation target for simplicity) and solving for the system's eigenvalues: the equilibrium is determinate (unique bounded solution) if and only if $\phi_\pi > 1$ and $\phi_y(1-\beta)/\kappa + \phi_\pi > 1$.

When $\phi_\pi < 1$ (violation of the Taylor principle), the equilibrium is **indeterminate**: multiple bounded equilibria exist, indexed by arbitrary self-fulfilling "sunspot" expectations. Any sunspot shock — random noise uncorrelated with fundamentals — can shift the economy between equilibria, generating volatility unrelated to technology or preferences. Clarida, Galí, and Gertler (2000) estimated that U.S. monetary policy before 1979 violated the Taylor principle ($\phi_\pi \approx 0.83$), consistent with the indeterminacy interpretation of the high-volatility 1970s.

---

## 23.2 Optimal Monetary Policy in the New Keynesian Model

### The Central Bank's Problem

In the canonical New Keynesian framework, the central bank minimizes a quadratic loss function:

$$\mathcal{L} = \frac{1}{2}\sum_{t=0}^\infty \beta^t\bigl[\pi_t^2 + \lambda_x\hat{x}_t^2\bigr],$$

where $\lambda_x > 0$ is the relative weight on output gap stabilization. This loss function can be derived as a second-order approximation to the representative household's welfare loss from inflation and output gap variability (Woodford, 2003). The two terms reflect: (i) the welfare cost of price dispersion generated by non-simultaneous price adjustment (Calvo model); and (ii) the welfare cost of output gap variation reflecting misallocation of resources between different goods.

The central bank minimizes this loss subject to the NK IS curve, the NKPC, and a resource constraint. Under **discretion** — optimizing period by period, taking private sector expectations as given — the first-order condition for $\pi_t$:

$$\pi_t + \frac{\lambda_x}{\kappa}\hat{x}_t = 0 \implies \hat{x}_t = -\frac{\kappa}{\lambda_x}\pi_t.$$

This **lean-against-the-wind** condition says the central bank reduces the output gap whenever inflation is positive. It is the period-by-period optimality condition under discretion, connecting the two stabilization objectives.

### The Divine Coincidence

In the standard NKPC without cost-push shocks ($\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t$), the discretionary equilibrium satisfies $\hat{\pi}_t = 0$ and $\hat{x}_t = 0$ at every date — the **divine coincidence** (Blanchard and Galí, 2007). Stabilizing inflation also stabilizes the output gap, and the two objectives are not in conflict. This holds because the only shocks in the standard model are demand shocks (which move $\hat{x}$ and $\pi$ in the same direction) and technology shocks (which move potential output but not the gap if the central bank adjusts appropriately).

The divine coincidence breaks down when cost-push shocks $u_t$ are added to the NKPC:

$$\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t + u_t.$$

Now a positive $u_t$ (supply shock) raises inflation and requires a negative $\hat{x}_t$ to stabilize inflation — a genuine trade-off between the two objectives. The optimal discretionary policy under this trade-off:

$$\hat{x}_t = -\frac{\kappa}{\lambda_x + \kappa^2/\lambda_x}\,u_t, \quad \pi_t = \frac{\lambda_x/\kappa}{\lambda_x/\kappa + \kappa}\,u_t.$$

Both inflation and the output gap move with the supply shock — some stabilization of each objective is sacrificed.

### Commitment Versus Discretion in the NK Model

Under **commitment**, the central bank can precommit to future policy rules at $t=0$, binding its future selves. The optimal commitment policy exploits the forward-looking nature of the NKPC: because current inflation depends on expected future output gaps (from the NKPC forward solution), the central bank can credibly reduce current inflation at lower cost by committing to a future recession that agents anticipate today.

The optimal commitment policy (in the presence of cost-push shocks) exhibits **history dependence**: the central bank keeps the output gap negative for longer after a positive cost-push shock than the period-by-period optimal policy would prescribe. This "promise to suffer in the future" reduces expected future inflation, which — via the forward-looking NKPC — reduces current inflation more cheaply than a purely contemporaneous approach.

The commitment policy is sometimes characterized as **price-level targeting**: instead of stabilizing inflation each period, the central bank stabilizes the *price level* around a target path, accepting above-target inflation to make up for past below-target inflation. This "make-up" strategy introduces beneficial history-dependence.

**Definition (Average Inflation Targeting).** **Average inflation targeting (AIT)** is a policy framework in which the central bank commits to keeping average inflation over a medium-term horizon (e.g., 5 years) near the target, rather than targeting inflation in each individual period. AIT is a practical approximation to optimal commitment policy under the NKPC: periods of below-target inflation are followed by above-target periods, maintaining the long-run average at target while reducing the output cost of disinflation. The Federal Reserve adopted AIT in August 2020.

---

## 23.3 Monetary Policy at the Effective Lower Bound

### The Constraints of the ELB

When the Taylor rule prescribes a negative nominal interest rate — as occurred in the United States after 2008, in the eurozone after 2012, and globally after 2020 — the ELB constrains monetary policy. The IS–LM model shows that at the ELB, conventional interest rate policy loses its traction: money-market equilibrium holds with idle cash accumulation, and output remains below potential even with the policy rate at zero. The NK model formalizes this: at the ELB, the DIS equation cannot be satisfied at $\hat{x} = 0$ when $r^n < i^{ELB}$:

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i^{ELB} - \mathbb{E}_t[\pi_{t+1}] - r_t^n) < 0.$$

The output gap is negative and the NKPC implies deflationary pressure, which raises the real interest rate further, deepening the recession — a **deflationary trap** that can be self-sustaining.

### Forward Guidance

**Definition (Forward Guidance).** **Forward guidance** is a monetary policy tool by which the central bank communicates its intended future path of interest rates to the public, with the aim of influencing longer-term interest rates and inflation expectations even when the current policy rate cannot be reduced further.

Forward guidance works through the term structure: from the DIS equation iterated forward, current output depends on the entire expected future path of real interest rates. If the central bank credibly commits to keeping rates at the ELB for longer than the market currently expects, it reduces expected future short rates, lowers long-term real rates, and stimulates current aggregate demand. The formal mechanism:

$$\hat{x}_t = -\sigma\sum_{k=0}^\infty\mathbb{E}_t[i_{t+k} - \pi_{t+k+1} - r_{t+k}^n].$$

By shifting $\mathbb{E}_t[i_{t+k}]$ downward for future $k$, forward guidance reduces $\hat{x}_t$ even when $i_t$ is constrained.

In practice, central banks have deployed two types of forward guidance. **Date-based guidance** ("rates will remain at zero until [date]") is simple to communicate but lacks the conditionality that makes it credible when the economy recovers faster than expected. **Outcome-based guidance** ("rates will remain at zero until unemployment falls below [threshold] and inflation reaches [level]") is conditioned on economic outcomes and is more robust to surprises in the recovery path. The Fed's 2012 guidance ("exceptionally low through mid-2015") was date-based; its 2020 AIT guidance ("until maximum employment has been reached and inflation has risen to 2% and is on track to moderately exceed 2% for some time") was outcome-based.

### Quantitative Easing

**Definition (Quantitative Easing).** **Quantitative easing (QE)** is a monetary policy tool in which the central bank purchases large quantities of long-duration assets — typically government bonds and mortgage-backed securities — with the aim of reducing term premia and long-term interest rates when the short-term policy rate cannot be reduced further.

QE operates through two main channels. The **portfolio balance channel** (Vayanos and Vila, 2009): central bank purchases of long-duration assets reduce the supply available to private investors, who are assumed to have preferred habitats for particular maturities. The reduced supply forces down term premia, reducing long-term yields. The **signaling channel**: QE signals the central bank's commitment to keeping rates low for an extended period, reducing expected future short rates and hence long rates via the expectations hypothesis.

The empirical evidence on QE effectiveness is broadly positive but uncertain in magnitude. Gagnon et al. (2011) estimate that the Fed's first QE program ($1.75 trillion, 2009–10) reduced 10-year Treasury yields by approximately 91 basis points. D'Amico and King (2013) estimate that QE2 ($600 billion, 2010–11) reduced 10-year yields by approximately 26 basis points per $100 billion. The effects were larger for long-duration and higher-risk assets (mortgage-backed securities) and propagated into corporate bond markets, though with significant attenuation.

### Negative Interest Rates

Several central banks pushed policy rates below zero: the ECB deposit facility rate fell to $-0.5\%$; the SNB and Danish National Bank to $-0.75\%$; the Bank of Japan to $-0.1\%$. Negative rates aim to penalize excess reserve holdings, encouraging banks to lend rather than hoard cash. The empirical evidence (Jobst and Lin, 2016; Bottero et al., 2019) finds that negative rates did transmit to lending rates in the short run, but the effects were attenuated compared to positive-rate cuts, consistent with the reversal rate framework of Chapter 18.

The combined toolkit of forward guidance, QE, and negative rates — alongside the use of AIT to generate make-up inflation — represents the modern central bank's response to the ELB constraint. Whether these tools are sufficient to prevent a prolonged deflationary trap — and whether they introduce new financial stability risks through their effects on bank profitability, asset price inflation, and risk-taking — are among the most important open questions in monetary economics.

---

*Next: Chapter 24 — The Business Sector*
