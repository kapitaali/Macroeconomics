# Appendix F — Important Economic Indicators and Their Interpretation

This appendix describes the most important macroeconomic indicators published by national statistical agencies and international organizations — what they measure, how they are constructed, their release timing, and how they are used in analysis. Understanding these indicators is essential for applying macroeconomic theory to real-world situations.

---

## F.1 Output and Activity Indicators

### Gross Domestic Product (GDP)
**Published by:** Bureau of Economic Analysis (U.S.), ONS (UK), Eurostat (EU), National Statistics offices globally.  
**Frequency:** Quarterly (advance, second, third estimates); annual revisions.  
**Timing:** Advance estimate ~30 days after quarter end; second ~60 days; third ~90 days.

GDP measures total real output. The advance estimate uses incomplete source data and is subsequently revised. The revision pattern matters for real-time policy analysis: the advance estimate of U.S. GDP has a mean absolute revision of ~1.2 pp at annual rates relative to the final. Analysts focus on: (i) the headline growth rate (quarter-on-quarter annualized); (ii) the contributions from each expenditure component ($C$, $I$, $G$, $NX$); and (iii) the implicit price deflator. A quarter of negative GDP growth is often called a "technical recession" — though the NBER Business Cycle Dating Committee uses a broader set of indicators to date recessions.

**Interpretation pitfalls:** Seasonally adjusted quarterly GDP is annualized in the U.S. (multiply by 4), so a single-quarter growth rate of 1% is reported as 4%. Non-annualized and annualized figures are sometimes confused in media reporting. Revisions to GDP can change the sign of growth: what appeared to be a positive quarter may be revised to negative months later.

### Industrial Production Index (IPI)
**Published by:** Federal Reserve (U.S.), Eurostat (EU), national central banks.  
**Frequency:** Monthly.  
**Timing:** ~15 days after month end.

The IPI measures real output in manufacturing, mining, and utilities sectors. It is a leading indicator for overall GDP because it is released earlier and covers the most cyclically sensitive sectors. A decline in IPI over several months typically precedes a GDP contraction. **Capacity utilization** — the ratio of actual output to potential output in the industrial sector — is published alongside IPI and indicates demand pressure on productive resources.

### Purchasing Managers' Indices (PMI)
**Published by:** S&P Global (formerly IHS Markit), Institute for Supply Management (ISM).  
**Frequency:** Monthly.  
**Timing:** First business day of the following month.

PMIs are surveys of purchasing managers asking whether business conditions are better, the same, or worse than the previous month. A reading above 50 indicates expansion; below 50, contraction. PMIs are **flash estimates** released before hard data and are widely followed as early cycle indicators. The composite PMI (manufacturing + services weighted) is particularly informative for overall economic conditions.

**Interpretation pitfalls:** PMIs measure the direction of change (diffusion indices), not the level or rate of change. A PMI of 52 following a PMI of 55 indicates slower expansion, not contraction.

---

## F.2 Labor Market Indicators

### Non-Farm Payrolls and Unemployment Rate
**Published by:** Bureau of Labor Statistics (U.S.), Eurostat, ILO.  
**Frequency:** Monthly.  
**Timing:** First Friday of the following month (U.S.).

The Employment Situation report contains two surveys: the **Establishment Survey** (payrolls: how many jobs were on business payrolls) and the **Household Survey** (unemployment: how many people report being unemployed and searching). These can diverge, particularly at business cycle turning points. Markets respond most strongly to payroll surprises relative to expectations.

Key sub-components: (i) average hourly earnings growth — a leading indicator for wage inflation; (ii) average weekly hours worked — a more flexible margin than employment; (iii) labor force participation rate — critical for assessing slack when unemployment rate is low; (iv) U-6 unemployment — the broadest measure of underutilization.

### Job Openings and Labor Turnover Survey (JOLTS)
**Published by:** Bureau of Labor Statistics (U.S.).  
**Frequency:** Monthly.  
**Timing:** ~35 days after month end (lagged one month behind payrolls).

JOLTS provides vacancy data essential for plotting the Beveridge curve ($v$ vs. $u$ space) and measuring labor market tightness $\theta = V/U$. The hires and separations components decompose employment changes into gross flows, revealing structural dynamics invisible in net payroll changes.

---

## F.3 Inflation Indicators

### Consumer Price Index (CPI)
**Published by:** Bureau of Labor Statistics (U.S.), Eurostat (EU).  
**Frequency:** Monthly.  
**Timing:** ~15 days after month end.

The CPI is the primary measure of consumer price inflation. Analysts distinguish: (i) **headline CPI** — all items; (ii) **core CPI** — ex food and energy; (iii) **trimmed mean CPI** (Cleveland Fed) — dropping price-change outliers from both tails. The Fed's preferred inflation measure is the **PCE deflator** (Personal Consumption Expenditures) rather than CPI, because: PCE uses chain weighting (less substitution bias), covers a broader basket including healthcare paid by employers, and has lower weight on shelter costs.

### Producer Price Index (PPI)
**Published by:** Bureau of Labor Statistics (U.S.).  
**Frequency:** Monthly.  
**Timing:** ~2 weeks after month end.

The PPI measures prices received by producers for their output — a leading indicator for CPI since input price changes eventually pass through to consumer prices. The PPI final demand index tracks prices for goods and services sold to final users.

### Inflation Expectations Indicators
Several measures track inflation expectations, which are central to Phillips curve dynamics:
- **University of Michigan Survey**: asks consumers about expected inflation 1 year and 5–10 years ahead. A breakout in 5–10 year expectations signals de-anchoring of long-run inflation expectations.
- **Survey of Professional Forecasters (SPF)**: asks professional forecasters about expected CPI and PCE inflation over various horizons. Considered the most credible survey-based measure for monetary policy.
- **Breakeven inflation rates**: the yield spread between nominal Treasury bonds and Treasury Inflation Protected Securities (TIPS) of the same maturity, $\pi^{BE} = i^{nominal} - i^{TIPS}$. Affected by liquidity premia and risk premia, so not a pure measure of expected inflation.
- **Cleveland Fed inflation expectations model**: combines surveys, TIPS, and the term structure to extract "true" inflation expectations purged of risk and liquidity premia.

---

## F.4 Financial Market Indicators

### Interest Rates
**Federal funds rate:** The overnight interbank rate that is the Fed's policy instrument. The Fed sets a target range; deviations from target are corrected through open market operations.  
**Treasury yield curve:** The relationship between yields and maturities for U.S. government bonds. The **yield spread** (10-year minus 2-year or 10-year minus 3-month) is a standard recession predictor: inversions (short rates above long rates) have preceded all U.S. recessions since 1960 with a lead time of 6–18 months.  
**LIBOR / SOFR spread over OIS:** The London Interbank Offered Rate (now replaced by SOFR — Secured Overnight Financing Rate) minus the overnight indexed swap rate measures credit risk in the interbank market. During the 2008 crisis, the LIBOR–OIS spread spiked to 364 bps, signaling frozen money markets.

### Asset Prices
**S&P 500 / equity indices:** Stock prices aggregate information about expected future corporate profits and discount rates. A standard leading indicator for GDP. The **CAPE ratio** (Shiller's cyclically adjusted price-earnings ratio using 10-year average real earnings) predicts long-horizon returns.  
**VIX:** The CBOE Volatility Index, measuring implied volatility of S&P 500 options over 30 days. A "fear index" — high VIX signals market stress and elevated uncertainty. Used by Bloom (2009) as a measure of macroeconomic uncertainty.  
**Credit spreads:** The yield premium of corporate bonds over equivalent-maturity Treasuries. High-yield (junk bond) spreads signal credit market stress and tightening financial conditions. The **Moody's Baa–Aaa spread** is a longer historical series measuring investment-grade credit risk.

---

## F.5 Monetary and Credit Indicators

### Monetary Aggregates (M1, M2)
The Fed publishes M1 and M2 weekly. Following 2020, the Fed discontinued M3 and redefined M1 to include savings deposits. The money multiplier $m = M2/H$ provides information about bank lending behavior — a declining multiplier (as after 2008) signals that reserve injections are not translating into money creation.

### Bank Credit and Lending Standards
The **Senior Loan Officer Opinion Survey (SLOOS)**, published quarterly by the Fed, asks whether banks are tightening or loosening lending standards. A net-tightening reading is a leading indicator of reduced business investment and consumer credit; it often precedes recessions.

### Federal Reserve Balance Sheet (H.4.1)
Weekly publication of Federal Reserve assets and liabilities. The growth in the Fed balance sheet from under $1T (2008) to over $8T (2022) reflects QE programs and emergency facilities. The composition of assets (Treasuries vs. MBS vs. special facilities) matters for understanding transmission channels.

---

## F.6 International and External Sector Indicators

### Current Account Balance
**Published by:** Bureau of Economic Analysis (U.S.), IMF Balance of Payments Statistics.  
**Frequency:** Quarterly.  
**Timing:** ~70 days after quarter end.

The current account balance equals net exports plus net factor income plus net transfers. A deficit means the domestic economy is spending more than it produces, financed by net capital inflows. Large sustained deficits can signal vulnerability to sudden stops. The **trade-weighted exchange rate** (published daily by the Federal Reserve) provides a broader measure than bilateral exchange rates, weighting partner currencies by trade shares.

### IMF World Economic Outlook (WEO)
**Published by:** IMF.  
**Frequency:** Twice yearly (April and October), with updates.

The WEO provides the IMF's projections for GDP growth, inflation, current accounts, and fiscal balances for ~190 countries, along with thematic analysis of global economic risks. It is widely used as the reference set of macroeconomic forecasts for cross-country analysis.

---

## F.7 Interpreting Economic Data Releases: A Framework

Economic data releases should be interpreted relative to expectations, not in absolute terms. Markets move on **surprises** — the deviation of actual data from the consensus forecast (Bloomberg or Reuters survey). A strong payroll number that exactly meets expectations causes little market reaction; a weak number that misses by 100k jobs causes a significant reaction.

A five-step framework for interpreting any data release:

1. **What did the data say?** Read the headline number and the key components.
2. **What did the market expect?** Compare to consensus forecast; compute the surprise.
3. **What are the relevant revisions?** Revised data for prior months may matter as much as current data.
4. **What does this imply for the policy path?** Does the surprise push the central bank toward tightening or easing? By how much and over what horizon?
5. **What are the caveats?** Seasonal adjustment issues? Sample changes? Likely future revisions? One month's data is always noisy — assess the trend, not the single observation.

---

*For international data sources, see Appendix J. For online data tools, see Appendix K.*
