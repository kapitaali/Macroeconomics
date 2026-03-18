# Appendix K — Online Resources and Tools

This appendix lists the primary online resources for macroeconomic research, organized by function: data retrieval, model solving, research access, and professional development.

---

## K.1 Data Retrieval and APIs

### FRED API (Federal Reserve Bank of St. Louis)
**URL:** https://fred.stlouisfed.org/docs/api/  
**Access:** Free, registration required for API key.  
**Languages:** R (`fredr` package), Python (`fredapi` package), direct JSON/XML.

The FRED API provides programmatic access to over 800,000 economic time series. A basic call in Python:
```python
from fredapi import Fred
fred = Fred(api_key='your_key')
gdp = fred.get_series('GDPC1')  # Real GDP, seasonally adjusted
```
Essential series codes: `GDPC1` (real GDP), `CPIAUCSL` (CPI), `UNRATE` (unemployment), `FEDFUNDS` (fed funds rate), `GS10` (10-year Treasury yield), `SP500` (S&P 500).

### World Bank Open Data API
**URL:** https://datahelpdesk.worldbank.org/knowledgebase/articles/889392  
**Access:** Free, no registration.  
**Python:** `wbdata` or `pandas_datareader`.

Provides access to WDI and other World Bank datasets. Indicator codes follow the pattern `NY.GDP.MKTP.KD` (GDP at constant USD) and are searchable at https://data.worldbank.org/indicator.

### IMF Data API
**URL:** https://datahelp.imf.org/knowledgebase/articles/667681  
**Access:** Free, no registration.

Provides WEO, BOP, IFS (International Financial Statistics), and other IMF databases via JSON API. The `imfpy` Python package simplifies access.

### OECD Data API
**URL:** https://data.oecd.org/api/  
**Access:** Free, no registration.

Access to OECD databases including Economic Outlook, Main Economic Indicators, and structural indicators.

### BIS Data API
**URL:** https://www.bis.org/statistics/full_data_sets.htm  
**Access:** Free download (bulk) or query via API.

Provides credit, property prices, effective exchange rates, and international banking statistics.

---

## K.2 Econometric and Modeling Software

### R
**URL:** https://www.r-project.org/  
**Primary use:** Econometrics, time-series analysis, data visualization.

The dominant language for empirical economics. Key packages:
- `dynlm`: Dynamic linear models with time-series support.
- `vars`: Vector autoregression estimation, Granger causality, impulse responses.
- `tseries`: Unit root tests (ADF, KPSS, PP).
- `urca`: Comprehensive unit root and cointegration testing (Johansen test).
- `sandwich`: HAC (Newey–West) standard errors.
- `AER`: Applied Econometrics with R — IV, 2SLS, hypothesis tests.
- `plm`: Panel data models (fixed/random effects, Hausman test).
- `bvarsv`: Bayesian VAR with stochastic volatility.
- `ggplot2`, `dplyr`, `tidyr`: Data manipulation and visualization (the tidyverse).

### Python
**URL:** https://www.python.org/  
**Primary use:** Data science, DSGE model solving, machine learning applications.

Key packages for macroeconomics:
- `pandas`: Time-series data manipulation.
- `statsmodels`: OLS, IV, time-series models (VAR, ARIMA, state-space).
- `linearmodels`: Panel data models, IV estimation.
- `pymc3` / `pymc`: Bayesian estimation via MCMC.
- `quantecon`: Tools from Ljungqvist and Sargent; dynamic programming, Markov chains.

### MATLAB / Octave
Traditional platform for DSGE model solving. The **Dynare** toolbox (below) runs on MATLAB and Octave and is the standard for DSGE estimation.

### Dynare
**URL:** https://www.dynare.org/  
**Access:** Free, open source (works with MATLAB or Octave).

Dynare is the standard tool for solving and estimating DSGE models. Given a `.mod` file specifying model equations, it: (i) linearizes around the steady state; (ii) solves the linear rational expectations model; (iii) computes impulse response functions; (iv) estimates parameters by MLE or Bayesian MCMC; (v) conducts model comparison. Used by central banks worldwide.

### Stata
**URL:** https://www.stata.com/  
**Primary use:** Microeconometrics, panel data, program evaluation.

Less common in pure macro than R or Python but standard in empirical work drawing on microeconomic data (household surveys, firm-level data). Key commands: `xtfe` (panel FE), `ivreg2` (IV), `xtivreg2` (panel IV), `xtabond2` (System GMM), `rdrobust` (RD estimation).

---

## K.3 DSGE Model Databases and Replication Archives

### MMB (Macroeconomic Model Comparison Facility)
**URL:** https://www.macromodelbase.com/  
**Description:** A collection of 100+ published DSGE models in a common Dynare framework, enabling systematic comparison of model properties and impulse responses across models.

### CEPR DP Archive
**URL:** https://cepr.org/publications/dp  
**Description:** Working papers from European macroeconomists. Subscription required for full access (many universities have access).

### NBER Working Papers
**URL:** https://www.nber.org/papers  
**Description:** The primary archive for U.S.-based macroeconomic research. Papers are freely available to the public after an embargo period.

### AEA Journals + Data/Code Archive
**URL:** https://www.aeaweb.org/journals  
**Description:** The American Economic Review and associated journals require authors to post data and code for replication. The replication archives contain the datasets and code for most published empirical papers in these journals.

### Social Science Research Network (SSRN)
**URL:** https://ssrn.com/  
**Description:** Preprint server for economics and other social sciences. Most working papers are freely available here before (and after) publication.

---

## K.4 Central Bank Research and Publications

### Federal Reserve Board Research
**URL:** https://www.federalreserve.gov/econres.htm  
Key resources:
- **FEDS Working Papers**: Staff research.
- **Finance and Economics Discussion Series**: Peer-reviewed working papers.
- **Monetary Policy Reports**: Semi-annual to Congress, with current economic assessment.
- **Beige Book**: Qualitative regional economic conditions, released 8× per year.

### Federal Reserve Banks Research
Each regional Fed bank publishes research and commentary:
- **NY Fed Liberty Street Economics**: Accessible commentary from NY Fed economists.
- **Chicago Fed**: Focus on banking and financial stability.
- **San Francisco Fed**: Focus on West Coast economy, international macroeconomics.
- **St. Louis Fed (FRED)**: Data, accessible research summaries.
- **Minneapolis Fed**: Strong research tradition in monetary economics.
- **Philadelphia Fed (RTDSM)**: Real-time macroeconomic datasets.

### European Central Bank (ECB)
**URL:** https://www.ecb.europa.eu/pub/research/  
Key resources:
- **Working Paper Series**: ECB and Eurosystem staff research.
- **Economic Bulletin**: Bi-weekly assessment of euro area conditions.
- **Economic Perspectives**: Accessible commentary for a broader audience.

### Bank for International Settlements (BIS)
**URL:** https://www.bis.org/  
Key resources:
- **BIS Working Papers**: Research on monetary policy, financial stability.
- **BIS Quarterly Review**: Financial market commentary and statistical tables.
- **Annual Economic Report**: Comprehensive annual assessment of global economy.

### IMF Research and Publications
**URL:** https://www.imf.org/en/Publications/  
Key resources:
- **IMF Working Papers**: Large archive of IMF staff research.
- **World Economic Outlook**: Semi-annual global economic assessment.
- **Global Financial Stability Report**: Semi-annual financial stability assessment.
- **Fiscal Monitor**: Semi-annual public finance assessment.
- **IMFBlog**: Accessible commentary on current economic issues.

---

## K.5 Learning Resources and Courses

### QuantEcon
**URL:** https://quantecon.org/  
Open-source lectures on quantitative economics (Python and Julia). Written by Thomas Sargent and John Stachurski. Covers dynamic programming, search models, asset pricing, Markov chains, and linear rational expectations models at graduate level.

### NBER Summer Institute Lectures
**URL:** https://www.nber.org/conferences/summer-institute  
Many sessions post lecture slides and videos. Particularly valuable: Macroeconomics and Productivity, International Finance and Macroeconomics, Monetary Economics, and Economic Fluctuations and Growth.

### AEA Continuing Education
**URL:** https://www.aeaweb.org/conference/cont-ed  
Annual intensive courses at the ASSA meetings in January. Recent topics have included machine learning in economics, empirical bayes methods, and identification in structural models.

### The Economics of Climate Change (Yale Climate Economics)
**URL:** https://sites.yale.edu/nordhaus/  
Resources on integrated assessment modeling, DICE/RICE models, and climate-economy linkages.

---

## K.6 Data Visualization and Communication

### Our World in Data
**URL:** https://ourworldindata.org/  
Accessible, carefully constructed visualizations of long-run global economic, social, and environmental data. Useful for communicating macroeconomic trends to non-specialist audiences.

### Bloomberg / Reuters / Financial Times
Professional financial news sources with real-time data. Bloomberg Terminal provides the most comprehensive real-time financial data but is expensive; the Bloomberg Economics portal (bloomberg.com/economics) provides accessible commentary.

### The Economist and Project Syndicate
**The Economist's** economics coverage and graphics are consistently high quality. **Project Syndicate** publishes short commentaries by leading economists on current policy questions.

---

*For international statistical databases, see Appendix J. For research methodology, see Appendix B.*
