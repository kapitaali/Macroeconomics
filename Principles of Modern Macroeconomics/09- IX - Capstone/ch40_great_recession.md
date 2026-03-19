# Chapter 40 — Case Study: The Great Recession of 2008 — Causes and Policy Responses

> *"We saved the world financial system from a second Great Depression. We didn't save it for the bankers. We saved it for the workers, the pensioners, and the savers."*
> — Gordon Brown, Prime Minister of the United Kingdom, 2009

---

The Global Financial Crisis of 2007–09 and the Great Recession that followed constitute the most severe macroeconomic event in the advanced world since the 1930s. Understanding it — tracing the buildup of vulnerabilities, the cascade of failures, the policy responses, and the aftermath — is the most important applied macroeconomics exercise in living memory. This chapter integrates the frameworks developed throughout the book to provide a comprehensive account of what happened and why, grounded in data and disciplined by theory.

---

## 40.1 Chronology and Macroeconomic Magnitudes

The crisis unfolded in three phases.

**Phase 1: Housing boom and subprime buildup (2000–2006).** U.S. house prices rose approximately 90% in real terms from 2000 to their peak in Q2 2006, driven by: historically low interest rates (the Fed funds rate fell to 1% in 2003–04); relaxed lending standards in the mortgage market; rapid financial innovation in mortgage securitization; and extrapolative expectations about continued house price appreciation. Household debt-to-GDP rose from 66% (2000) to 99% (2007). The non-agency securitization market — mortgage-backed securities (MBS), collateralized debt obligations (CDOs) backed by subprime loans — grew from near zero in 2000 to approximately $2.5 trillion in outstanding balances by 2007.

**Phase 2: Financial system crisis (2007–2008).** House prices began declining in mid-2006; delinquency rates on subprime mortgages rose sharply from early 2007. In June–July 2007, two Bear Stearns hedge funds heavily invested in subprime CDOs failed. In August 2007, the interbank market froze: LIBOR–OIS spread spiked from 10 basis points to 40 basis points — banks refused to lend to each other because they could not assess each other's subprime exposures. The investment bank Bear Stearns failed in March 2008 and was acquired by JPMorgan Chase with Fed assistance. The GSEs Fannie Mae and Freddie Mac were placed in government conservatorship in September 2008. On September 15, 2008, Lehman Brothers — the fourth-largest U.S. investment bank, with $639 billion in assets — filed for bankruptcy. The LIBOR–OIS spread reached 364 basis points. Money market mutual funds "broke the buck" (net asset values fell below $1.00). Global financial markets froze.

**Phase 3: Recession and recovery (2008–2016).** U.S. real GDP fell 4.3% from Q4 2007 to Q2 2009. Unemployment rose from 4.4% (May 2007) to 10.0% (October 2009). The Case–Shiller national home price index fell 33% from its 2006 peak to its 2012 trough. The S&P 500 fell 57% from October 2007 to March 2009. Output did not return to its pre-crisis trend in most advanced economies by 2016; some economists argued it never would (hysteresis).

---

## 40.2 Origins: Vulnerabilities and Their Sources

The crisis had deep roots in structural vulnerabilities that built up over years, not months. Identifying them requires integrating financial economics, macroeconomics, and regulatory economics.

### The Originate-to-Distribute Model

The traditional mortgage model: a bank originates a mortgage, holds it on its balance sheet, and therefore has a strong incentive to screen borrowers carefully. The **originate-to-distribute (OTD) model**: a bank or non-bank lender originates a mortgage and immediately sells it into the securitization market, retaining no long-run exposure. The incentive to screen disappears — if the loan will be sold, why spend resources evaluating the borrower's creditworthiness?

Keys, Mukherjee, Seru, and Vig (2010) provide quasi-experimental evidence of this screening failure: loans with FICO scores just above 620 (a common securitization threshold) were significantly more likely to be securitized and significantly more likely to default than loans just below the threshold — consistent with the OTD model reducing screening effort for securitized loans.

### Leverage, Complexity, and Information Opacity

Investment bank leverage ($\Lambda = A/E$) rose from approximately 15 in 2000 to approximately 35 in 2007. At $\Lambda = 35$, a 3% decline in asset values wipes out all equity — the bank becomes insolvent. The risk management models used by banks (Value-at-Risk, VaR) treated subprime MBS as low-risk based on historically low mortgage default rates and their apparent diversification (mortgages across different regions). These models catastrophically underestimated correlation risk: in a nationwide house price decline, regional diversification provides no protection.

CDO structures added another layer of opacity. A CDO pools 100 different MBS tranches and issues new securities from the pool. The "super-senior" CDO tranche — with the first claim on the pool's cash flows — was rated AAA by Moody's and S&P, implying negligible default probability. This rating was based on assumed independence of the underlying mortgages: if default rates were independent, the probability that enough mortgages defaulted simultaneously to impair the super-senior tranche was vanishingly small. But independence was a false assumption: house prices fell everywhere simultaneously, correlating all the underlying mortgages. When the models' correlation assumptions failed, the "safe" AAA tranches suffered catastrophic losses.

### Regulatory Gaps and Moral Hazard

The shadow banking system — money market mutual funds, repo markets, special investment vehicles (SIVs), structured investment vehicles — performed bank-like maturity transformation (borrowing short, lending long) but was not subject to bank capital requirements, deposit insurance, or supervisory oversight. When confidence collapsed in September 2008, the shadow banking system experienced a classic run: money market funds were subject to redemption demands; repo markets froze as counterparties demanded haircuts that made overnight funding impossible; SIVs could not roll over their commercial paper.

The regulatory architecture had not kept pace with financial innovation. The Gramm-Leach-Bliley Act (1999) repealed the Glass-Steagall separation of commercial and investment banking; the Commodity Futures Modernization Act (2000) excluded OTC derivatives from regulation; Basel II capital rules allowed banks to hold less capital against "low-risk" AAA-rated securities — creating a regulatory incentive to accumulate the very assets that proved most toxic.

---

## 40.3 Propagation: The Financial Accelerator in Action

The mechanics of crisis propagation illustrate the financial accelerator (Chapter 24) operating at full force and at multiple scales simultaneously.

### The Bank Run on Shadow Banking

The run on shadow banking began when money market funds holding Lehman Brothers commercial paper "broke the buck" on September 16, 2008. Institutional investors withdrew $300 billion from prime money market funds within a week. These funds had financed themselves through overnight repo, lending to investment banks against collateral. As money market funds withdrew from the repo market, investment banks lost their primary short-term funding source, forcing asset sales at fire-sale prices that further depressed valuations.

The fire-sale mechanism: Bank A, forced to sell MBS at 70 cents on the dollar, drives down the market price of MBS held by Bank B to 70 cents (even if B's loans have not defaulted). Bank B marks its MBS to market (accounting rules requiring mark-to-market for trading books), reducing its reported equity. B's capital ratio falls below regulatory minimums, forcing it to either raise capital (impossible in a panic) or sell assets (intensifying the fire sale). The fire sale thus becomes self-reinforcing — the **cascade dynamics** of a systemic financial crisis.

### Credit Tightening and the Real Economy

The propagation from the financial sector to the real economy operated through multiple channels:

**Credit channel**: the external finance premium spiked (the Gilchrist–Zakrajšek EBP rose from near zero to approximately 200 basis points), dramatically raising borrowing costs for non-financial firms. Small and medium enterprises, highly dependent on bank credit, cut investment sharply. The Federal Reserve's Senior Loan Officer Survey showed the sharpest tightening of lending standards in its history.

**Wealth effect**: household net worth fell by approximately $13 trillion (2007–2009) as house prices and equity prices collapsed simultaneously. The MPC out of housing wealth (approximately 6 cents per dollar) and equity wealth (approximately 3 cents per dollar) implied a consumption reduction of approximately $700–900 billion — a severe drag on aggregate demand.

**Confidence collapse**: investment collapsed (–24% peak-to-trough) as uncertainty spiked (VIX reached 80 in November 2008), real options values of waiting rose, and the accelerator operated through the investment-net worth-EFP channel.

---

## 40.4 Policy Responses: Monetary, Fiscal, and Financial

### Monetary Policy

The Fed reduced the federal funds rate from 5.25% (September 2007) to 0–0.25% (December 2008), the fastest easing cycle in the post-war period. The ELB was reached 15 months before the recession ended — leaving no remaining conventional space. The Fed then deployed three rounds of quantitative easing:

- **QE1** (November 2008–March 2010): $1.75 trillion in MBS, agency bonds, and Treasuries — focused on restoring mortgage market functioning.
- **QE2** (November 2010–June 2011): $600 billion in Treasuries — aimed at further reducing long-term interest rates.
- **QE3** (September 2012–October 2014): open-ended purchases ($85 billion/month initially) conditional on labor market improvement — the first outcome-based QE program.

The Fed also introduced nine emergency credit facilities (TALF, CPFF, PDCF, AMLF, MMIFF, TSLF, ABCP MMMF, MMMF Liquidity Facility, Primary Dealer Credit Facility) to restore functioning to specific markets; extended dollar swap lines to 14 foreign central banks; and participated in the coordinated G7 central bank rate cuts of October 2008.

### Fiscal Policy

The **American Recovery and Reinvestment Act (ARRA, February 2009)** provided $787 billion over 10 years: approximately $290 billion in tax cuts and $500 billion in spending (infrastructure, education, healthcare, state fiscal relief, extended UI). The Romer–Bernstein (2009) CEA projection used spending multipliers of 1.57 and tax multipliers of 0.99, projecting that the ARRA would save or create 3.5 million jobs by end-2010.

Subsequent analysis: Wilson (2012), using cross-state variation in ARRA grants for identification, estimates a spending multiplier of approximately 1.5 over 2 years. Chodorow-Reich et al. (2012) estimate that transfers to states (ARRA fiscal relief) generated approximately 140,000 state and local government jobs, consistent with multipliers in the range 1.5–2.0.

The Troubled Asset Relief Program (TARP, October 2008): $700 billion authorized; approximately $426 billion disbursed. TARP purchased equity stakes in banks, GM, and Chrysler; supported AIG (the insurance company that had written credit default swap protection on $440 billion of CDOs). By 2013, TARP had recovered $442 billion — the program ultimately made a small profit on the banking investments (though lost substantially on auto and AIG rescue).

### Financial System Stabilization

The SCAP stress tests (May 2009) were a turning point. The Fed publicly disclosed the capital shortfalls of 19 largest banks and gave them 6 months to raise capital privately or accept government equity. The transparency — unprecedented at the time — restored market confidence: banks were able to raise $140 billion of private capital within months. The stress test methodology became the template for subsequent Basel III-required stress testing globally.

---

## 40.5 The Slow Recovery and Lessons

The U.S. recovery was the slowest of any post-WWII recession: output did not return to its 2007 peak trend until approximately 2016 by some measures, and never did by others. Hysteresis appears to have been substantial: the labor force participation rate fell from approximately 66% (2007) to approximately 62.5% (2015) and did not recover.

Several theoretical lessons follow from the 2008–09 experience:

**DSGE models without financial sectors failed to identify or predict the crisis.** The New Keynesian models in use at central banks in 2007 had complete asset markets, no default risk, and no banking sector. These models could not generate the fire-sale dynamics, credit crunches, or liquidity spirals that characterized the crisis. The post-crisis research program (financial accelerator models, HANK, intermediary asset pricing) is a direct response to this failure.

**Microprudential regulation was insufficient for macroprudential stability.** Each bank could be individually sound while the system as a whole was fragile. The regulatory framework needed macroprudential tools — countercyclical capital buffers, systemic risk surcharges, leverage limits — to address system-wide externalities.

**The ELB is not a theoretical curiosity.** The federal funds rate was at zero for 7 years (December 2008–December 2015), far longer than any pre-crisis model considered. Unconventional monetary policy (QE, forward guidance) became permanent features of the toolkit rather than emergency departures.

**Fiscal policy matters at the ELB.** The theoretical prediction that ELB fiscal multipliers substantially exceed one appears to have been borne out: ARRA spending generated employment roughly consistent with a multiplier of 1.5–2.0.

**Recovery from balance-sheet recessions is slow.** Reinhart and Rogoff (2009) documented that recoveries from financial crises take, on average, 7–10 years to restore peak employment and output — far longer than recoveries from ordinary recessions. The 2008–09 recession fit this pattern.

---

*Next: Chapter 41 — Case Study: The COVID-19 Pandemic*
