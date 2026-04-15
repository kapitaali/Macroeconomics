# Chapter 29 — Monetary Policy in Practice: Interest Rates and Quantitative Easing

> *"The first priority of a central bank is to keep the value of its currency stable... To achieve this, it must be willing to act preemptively."*
> — Janet Yellen, *A Minsky Meltdown: Lessons for Central Bankers*, 2009

---

Chapter 23 developed the theory of optimal monetary policy in the New Keynesian framework: a benevolent, fully informed central bank solving a welfare optimization problem with perfect commitment. Real-world monetary policy involves none of these luxuries. Information arrives with lags and is substantially revised; models are contested; the natural rate of interest is unobservable and shifts over time; and credibility must be earned through repeated interaction rather than simply assumed. This chapter examines how monetary policy is actually conducted — the design of interest rate rules, the record of inflation targeting, and the unconventional tools deployed after 2008 — and what the empirical evidence says about the quantitative effects of each.

---

## 29.1 The Taylor Rule in Practice and Its Variants

### Specification and Variants

The Taylor (1993) rule provides a benchmark for evaluating whether policy is tight or loose relative to economic conditions:

$$i_t = r^n + \pi^* + \phi_\pi(\pi_t - \pi^*) + \phi_y \hat{x}_t,$$

with Taylor's original calibration $\phi_\pi = 1.5$, $\phi_y = 0.5$. The **Taylor principle** requires $\phi_\pi > 1$: the nominal rate must rise more than one-for-one with inflation so that the real rate rises when inflation increases, stabilizing inflation expectations [Ch. 23].

Several variants add interest rate smoothing (reflecting central banks' preference for gradual rate changes):

$$i_t = \rho_i i_{t-1} + (1-\rho_i)[r^n + \pi^* + \phi_\pi(\pi_t - \pi^*) + \phi_y\hat{x}_t] + \epsilon_t^m,$$

with $\rho_i \approx 0.85$ estimated for the Fed over 1987–2007. Smoothing reduces market volatility from policy surprises and provides additional accommodation at the ELB (because expectations of future rate cuts follow a gradual path rather than an immediate drop).

### The Real-Time Data Problem

A persistent challenge in evaluating Taylor rule adherence is that policymakers observe real-time data subject to revision, not the revised data used in academic evaluations. Orphanides (2001) showed that the apparent deviation of 1970s Fed policy from the Taylor principle largely disappears when real-time output gap estimates (available from the Fed's own Greenbook forecasts) are used instead of ex post revised data: the Fed was roughly following the Taylor rule given what it knew at the time, but was systematically overestimating potential output — a data problem, not a preference problem.

### The Neutral Rate Estimation Problem

The Taylor rule requires $r^n$, the neutral real rate, but $r^n$ is not directly observable. The Laubach-Williams (2003) Kalman filter model estimates that U.S. $r^n$ declined from approximately 3.5% in the early 1990s to approximately $-0.5\%$ by 2016 — a decline of more than 4 percentage points over two decades. If $r^n$ is near zero or negative, the Taylor rule prescribes a near-zero nominal rate even at full employment and 2% inflation, implying that the ELB binds frequently under normal conditions. This structural shift in $r^n$ — driven by demographic trends, slowing productivity growth, and global saving glut dynamics — is one of the most important and least certain empirical findings in modern macroeconomics.

---

## 29.2 Inflation Targeting: Design and Empirical Performance

### The Framework

Inflation targeting was introduced by the Reserve Bank of New Zealand in 1990 and adopted by the Bank of England, Sweden, Canada, and Australia through the 1990s. By 2020, more than 40 countries were using formal inflation targeting frameworks. The standard design has four components: (i) a publicly announced numerical inflation target (typically 2% CPI); (ii) **instrument independence** — the central bank freely chooses its policy rate; (iii) transparency through regular Inflation Reports, published forecasts, and post-meeting press conferences; and (iv) accountability through parliamentary testimony and legally mandated open letters when the target is missed by more than a specified margin.

### The Evidence on Credibility and Sacrifice Ratios

One of the principal claims for inflation targeting is enhanced credibility: if households and wage-setters believe the central bank will maintain low inflation, the sacrifice ratio for any necessary disinflation is reduced because expectations of inflation fall quickly alongside actual disinflation, rather than lagging as under adaptive expectations.

The empirical evidence is mixed. Ball and Sheridan (2005) found that inflation fell similarly in OECD targeting and non-targeting countries during the 1990s, raising the question of whether the adoption of a formal target itself caused the disinflation or merely coincided with a global disinflationary trend. Subsequent event-study analyses exploiting the timing of individual adoptions have been more favorable, finding that adopting inflation targeting reduced long-run inflation expectations (as measured by inflation compensation embedded in indexed bond yields) by approximately 0.5–1.0 percentage points over three years.

### Average Inflation Targeting

The Federal Reserve's August 2020 framework revision introduced **Average Inflation Targeting (AIT)**: the Fed would henceforth seek inflation that averages 2% over time, meaning that periods of below-target inflation would be followed by deliberate above-target periods. This "make-up strategy" was designed to address a specific problem: if the ELB binds frequently and the Fed always aims to return to 2% as quickly as possible, average realized inflation will be systematically below 2%, keeping the neutral nominal rate (equal to $r^n + \pi^*$) low and making ELB episodes even more frequent. By committing to above-2% inflation following below-target periods, AIT raises the expected average policy rate and reduces the probability of persistent ELB binding.

AIT's practical implementation during the 2021–22 inflation surge — when the Fed delayed rate increases while inflation exceeded 5% and then 8% — raised sharp questions about whether the framework provided sufficient guidance for an inflationary shock of unexpected magnitude and whether the make-up strategy had been interpreted so loosely as to be indistinguishable from discretion. The subsequent rapid tightening and the decline of inflation back toward target by 2023–24 without a severe recession suggested that credibility was not permanently damaged, but the episode prompted significant reassessment of AIT's operational specifics.

---

## 29.3 Quantitative Easing: Mechanisms and Evidence

When the policy rate reaches the ELB, the standard interest rate channel is exhausted. **Quantitative easing (QE)** extends the toolkit by having the central bank purchase large quantities of longer-term assets, reducing their yields and putting downward pressure on a broader range of borrowing costs.

### Theoretical Channels

**Portfolio balance channel (preferred habitat theory):** Vayanos and Vila (2009) model investors as having preferred maturities (pension funds prefer long-duration bonds to match long-duration liabilities; money market funds prefer short-duration assets). When the central bank purchases long-duration bonds, preferred-habitat investors who held those bonds must reinvest in shorter-duration or alternative assets, bidding down yields across the maturity and credit spectrum. The effect is proportional to the quantity purchased relative to the outstanding supply of long-duration bonds.

**Signaling channel:** QE purchases may signal to the market that the central bank intends to keep short rates low for an extended period, reducing long rates through the expectations component of the yield (rather than the term premium). Bauer and Rudebusch (2014) decompose QE's yield effects into signaling (approximately 50%) and portfolio balance (approximately 50%) components.

**Market functioning channel:** in the acute phase of a financial crisis (March-April 2008, September-November 2008, March 2020), credit markets seize up — bid-ask spreads widen to multiples of normal levels, market-making capacity evaporates, and liquidity premia dominate yields. Central bank purchases directly restore market liquidity, reducing the liquidity premium component of yields independently of expected future rates.

### Empirical Evidence on Yield Effects

**Federal Reserve QE1** (November 2008–March 2009, $\$1.75$ trillion in MBS, agency debt, and Treasuries): Gagnon et al. (2011) estimate that the program reduced 10-year Treasury yields by approximately 90–100 basis points using event-study methods and time-series regressions. The 30-year mortgage rate fell by a comparable magnitude.

**Federal Reserve QE2** (November 2010–June 2011, $\$600$ billion in Treasuries): estimated yield reduction of approximately 20–30 bps per $\$100$ billion purchased (D'Amico and King, 2013). The smaller effect per dollar reflects both the shorter expected duration of QE2 and the faster market adjustment by the time of announcement.

**ECB PSPP (Public Sector Purchase Programme)**, launched March 2015 at €60 billion/month eventually reaching €80 billion/month: Altavilla, Carboni, and Motto (2015) estimate 10-year sovereign yield reductions of 0.5–0.7 percentage points across the eurozone in the first year, with larger effects in peripheral countries (Italy, Spain) where the portfolio balance effect was reinforced by a reduction in redenomination risk premia.

### Transmission to the Real Economy

The crucial and less certain step is the pass-through from lower bond yields to real economic activity. Lower long-term yields reduce borrowing costs for firms (corporate bond yields) and households (mortgage rates), stimulate investment and durable goods purchases, depreciate the exchange rate (improving net exports), and raise asset prices through valuation effects (the lower discount rate raises the present value of future cash flows, as in the Gordon growth model [Ch. 20]).

Weale and Wieladek (2016) estimate that the Bank of England's QE programs raised UK GDP by approximately 1.5% and inflation by approximately 1.25 percentage points relative to the no-QE counterfactual, using a suite of VAR and DSGE models. For the Federal Reserve's programs, Chung et al. (2012) at the Fed estimate slightly smaller effects: approximately 1% on GDP and 0.5 percentage points on core PCE inflation.

These estimates imply that QE is a meaningful but modest instrument relative to conventional rate cuts. The 90 bps yield reduction from QE1, for example, is equivalent to roughly 2–3 conventional rate cuts of 25 bps each (given standard estimates of pass-through from policy rate to long rates). Central banks can thus compensate partly for ELB constraints through large-scale asset purchases, but the substitution is imperfect.

---

## 29.4 Forward Guidance: Credibility and the Paradoxes

### The Theoretical Case

Forward guidance is theoretically the most potent unconventional tool. The NK-IS curve depends on the entire expected future path of the real interest rate:

$$\hat{x}_t = -\sigma\sum_{k=0}^\infty\mathbb{E}_t[i_{t+k} - \pi_{t+k+1} - r_{t+k}^n].$$

A credible commitment to keep rates at zero for five years instead of two reduces this sum substantially, stimulating current output without any immediate change in the balance sheet. Because the IS curve responds to the entire future rate path, forward guidance operates with zero implementation lag — unlike spending programs, which must be designed, appropriated, and contracted.

### The Forward-Guidance Puzzle

The puzzle: in the standard NK model, forward guidance is implausibly powerful. A credible commitment to zero rates for five extra quarters can raise current output by several percentage points — far more than VAR evidence on the actual effects of forward guidance suggests. This extreme sensitivity to distant future expectations arises because the NK model's forward solution compounds expectations over an infinite horizon: each future period's expected real rate enters the current IS curve with equal weight.

The behavioral NK model (Chapter 15) resolves this by cognitive discounting: boundedly rational agents discount future policy promises by $M^k$ over $k$ periods ahead, attenuating the effects of forward guidance on distant horizons. With $M = 0.85$ (Gabaix's estimate), a five-year forward guidance commitment has roughly one-third the effect it would have under full rationality — much closer to the empirically observed magnitudes.

### Forms and Credibility

**Date-based guidance** ("rates will remain at zero through mid-2015") is easy to communicate but creates a time-inconsistency problem: if the economy recovers faster than expected, the central bank faces pressure to raise rates before the announced date, undermining the commitment.

**Outcome-based (threshold) guidance** ("rates will remain at zero until unemployment falls below 6.5% or inflation exceeds 2.5%") conditions the rate path on economic outcomes, avoiding the calendar problem. It is harder to communicate and introduces uncertainty about the rate path — households must forecast inflation and unemployment in addition to the central bank's reaction function — but it is more theoretically coherent and credibility-preserving.

**Odyssean versus Delphic guidance** (Campbell et al., 2012): Delphic forward guidance is a forecast of the future rate path based on current information — informative but not binding. Odyssean guidance commits the central bank to a future action, tying itself to the mast like Odysseus to resist the siren call of short-run deviations. Only Odyssean guidance has the expectations-channel effects that make forward guidance valuable at the ELB; Delphic guidance merely improves information, which has smaller and less predictable effects.

---

## 29.5 Monetary Policy and Financial Stability

### The Leaning Debate

A long-running debate concerns whether central banks should use interest rates to "lean against" financial imbalances — raising rates above the Taylor rule when credit and asset price growth is elevated, to reduce the probability of a future crisis — or whether they should "clean up" after crises with aggressive monetary easing while using macroprudential tools to limit imbalances.

**The case for leaning** (Borio and White, 2004; BIS research tradition): financial imbalances build slowly and are not captured by current inflation and output gap data; by the time a crisis materializes, it is too late; a central bank that ignores financial conditions may face much larger eventual instability costs that exceed the current output cost of preemptive rate increases.

**The case against leaning** (Svensson, 2017): interest rates are a blunt instrument that affects all borrowing, not just the risky kind. The benefits of leaning (lower crisis probability, with long and uncertain lags) are small relative to the costs (higher unemployment and lower output in the near term, with certainty). Moreover, crises often originate in shadow banking or foreign currency positions that interest rate increases affect only indirectly. Macroprudential tools — countercyclical capital buffers, dynamic provisioning, LTV and DTI limits — are better targeted at the financial stability objective without the broad macroeconomic collateral damage.

### The Emerging Consensus

The post-crisis consensus, embodied in frameworks at the IMF, the BIS, and most inflation-targeting central banks, distinguishes clearly between the price stability mandate (pursued with the policy rate) and the financial stability mandate (pursued primarily with macroprudential tools, with the policy rate as a backstop when macroprudential tools are unavailable or exhausted). The 2008 crisis demonstrated that ignoring financial stability is catastrophically costly; but it also demonstrated that interest rates alone cannot contain financial imbalances without unacceptable economic costs.

---

## Chapter Summary

- The **Taylor rule** $i_t = r^n + \pi^* + 1.5(\pi_t-\pi^*) + 0.5\hat{x}_t$ is the operational benchmark for evaluating monetary policy; the Taylor principle ($\phi_\pi > 1$) is necessary for determinacy. Real-time data problems (Orphanides) and the declining neutral rate (Laubach-Williams) are the main practical challenges.

- **Inflation targeting** — with operational independence, transparency, and accountability — has been adopted by more than 40 countries since 1990. Evidence suggests it reduced inflation levels and improved anchoring of expectations; **AIT** extends this by committing to make up for past below-target periods, aiming to raise average inflation and reduce ELB frequency.

- **Quantitative easing** works through portfolio balance (removing duration from the market), signaling (committing to low rates), and market functioning (restoring liquidity). Empirical estimates place QE1's yield reduction at approximately 90 bps; pass-through to GDP is meaningful but modest (approximately 1% per $\$1.75$ trillion).

- **Forward guidance** is theoretically the most powerful ELB tool, working through the IS curve's dependence on the entire expected rate path. The forward-guidance puzzle (implausibly large NK model effects) is resolved by cognitive discounting in the behavioral NK model. Outcome-based (Odyssean) guidance is more credibility-preserving than date-based guidance.

- The **leaning versus cleaning** debate: the emerging consensus assigns financial stability primarily to macroprudential tools, with monetary policy as a backstop — while acknowledging that ignoring financial imbalances (as pre-2008 central banks largely did) is costly.

---

*Next: Chapter 30 — Inflation and Deflation: Causes, Consequences, and Cures*
