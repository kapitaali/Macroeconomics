# Appendix J — International Comparisons and Data Sources

This appendix describes the principal international datasets used in macroeconomic research, the methodological issues that arise in cross-country comparisons, and guidance on accessing and interpreting these sources.

---

## J.1 Methodological Issues in Cross-Country Comparisons

Comparing economic quantities across countries requires converting national currencies into a common unit. Two approaches exist.

### Market Exchange Rate Conversion

The simplest approach converts national GDP at market exchange rates: if Finnish GDP is €280 billion and the EUR/USD rate is 1.10, Finnish GDP in USD is $308 billion. Market exchange rates reflect the prices of internationally traded goods but not non-traded goods, which are systematically cheaper in lower-income countries (the Balassa–Samuelson effect). This means that market-exchange-rate-converted GDP understates the real living standards of developing countries relative to rich countries, because a dollar of income buys more real goods and services in a low-income country than in a high-income one.

### Purchasing Power Parity Conversion

PPP conversion equates the purchasing power of currencies by computing the exchange rate at which a reference basket of goods costs the same in both countries. PPP-adjusted GDP better reflects true differences in real living standards and is the appropriate metric for welfare comparisons. The Penn World Tables (PWT) and the World Bank's World Development Indicators (WDI) both provide PPP-adjusted GDP series.

PPP exchange rates are estimated from the International Comparison Program (ICP), which collects price data for comparable goods and services across countries approximately every six years. Revisions to the ICP benchmark can substantially change PPP estimates and hence the ranking of countries by real income — the 2011 ICP revision significantly raised the estimated GDP of China and other Asian economies relative to earlier estimates.

---

## J.2 Principal International Datasets

### Penn World Tables (PWT)
**Provider:** University of Groningen, Center for Growth and Development.  
**URL:** https://www.rug.nl/ggdc/productivity/pwt/  
**Coverage:** 183 countries, 1950–present (updated annually).  
**Key variables:** Real GDP (output-side, expenditure-side, income-side), capital stock, total factor productivity, employment and hours worked, PPP exchange rates, price levels.  
**Reference:** Feenstra, Inklaar, and Timmer (2015), AER.

The PWT is the standard dataset for international growth research. Its distinguishing feature is multiple GDP concepts: the **output-side** measure (production at domestic prices) is appropriate for productivity analysis; the **expenditure-side** measure (absorption at PPP prices) is appropriate for welfare comparisons across countries.

### World Development Indicators (WDI)
**Provider:** World Bank.  
**URL:** https://datatopics.worldbank.org/world-development-indicators/  
**Coverage:** ~217 countries, 1960–present (updated annually).  
**Key variables:** GDP (current and constant prices, PPP), population, life expectancy, education indicators, poverty rates, trade, financial sector indicators.

The WDI is the most comprehensive single source for development-related macroeconomic indicators. It draws on national statistical offices, the IMF, WHO, ILO, and UNESCO. The emphasis is on development indicators rather than business cycle analysis.

### Maddison Project Database
**Provider:** University of Groningen, following the work of Angus Maddison.  
**URL:** https://www.rug.nl/ggdc/historicaldevelopment/maddison/  
**Coverage:** Up to 169 countries; historical series from 1820 (some series from year 1 CE).  
**Key variables:** GDP per capita (2011 PPP international dollars), population.  
**Reference:** Bolt, Inklaar, de Jong, and van Zanden (2018).

The Maddison Project is the essential source for very long-run historical comparisons of living standards. It documents the Great Divergence — the emergence of large income gaps between currently rich and currently poor countries — and its timing relative to industrialization.

### IMF World Economic Outlook (WEO) Database
**Provider:** International Monetary Fund.  
**URL:** https://www.imf.org/en/Publications/WEO  
**Coverage:** ~190 countries; historical data plus 5-year projections, updated April and October each year.  
**Key variables:** GDP growth, inflation (CPI), unemployment, current account, fiscal balance, government debt, trade.

The WEO database is the standard source for international macroeconomic aggregates in policy settings. Its projections are the IMF's official economic forecasts. Historical data are harmonized by IMF staff and are more comparable across countries than raw national statistics.

### OECD Economic Outlook Database
**Provider:** Organisation for Economic Cooperation and Development.  
**URL:** https://www.oecd.org/economic-outlook/  
**Coverage:** 38 OECD member countries plus major non-members (Brazil, China, India, Russia, South Africa).  
**Key variables:** GDP components, employment, wages, prices, fiscal accounts, current accounts, interest rates.

The OECD Economic Outlook provides the most detailed and methodologically consistent data for advanced economies. Its **OECD Main Economic Indicators (MEI)** database gives monthly data for a wide range of economic indicators. The **OECD Economic Policy Reforms: Going for Growth** publication provides comparable structural indicators (labor market regulation, product market regulation, tax rates) across OECD countries.

### BIS Statistical Warehouse
**Provider:** Bank for International Settlements.  
**URL:** https://stats.bis.org/  
**Coverage:** 44+ reporting countries; back to 1970s for many series.  
**Key variables:** Credit to private sector, property prices, exchange rates, interest rates, international banking, derivatives.

The BIS warehouse is the primary source for financial stability indicators. Its **total credit to private non-financial sector** series is widely used to measure credit cycles and financial vulnerabilities. The BIS **effective exchange rate indices** (narrow and broad) provide trade-weighted exchange rate measures for 60+ economies.

### World Inequality Database (WID)
**Provider:** World Inequality Lab, Paris School of Economics.  
**URL:** https://wid.world/  
**Coverage:** 100+ countries; historical series for some countries back to 1800.  
**Key variables:** Top income shares (top 1%, 10%), wealth distribution, fiscal income, national income.  
**Reference:** Piketty, Saez, Zucman, and collaborators.

WID is the definitive source for long-run distributional data, constructed from tax records, national accounts, and surveys. It provides both income inequality (shares of national income accruing to different percentiles) and wealth inequality (distribution of net worth).

### FRED (Federal Reserve Economic Data)
**Provider:** Federal Reserve Bank of St. Louis.  
**URL:** https://fred.stlouisfed.org/  
**Coverage:** United States primary; substantial international data.  
**Key variables:** Every major U.S. macroeconomic time series plus international data from BIS, World Bank, OECD.

FRED is the most user-friendly source for U.S. time-series data and is the de facto standard for economists working on U.S. macroeconomics. It includes over 800,000 series and allows direct API access for automated data retrieval. The **ALFRED** (Archival FRED) database provides vintaged data essential for real-time analysis.

---

## J.3 Bilateral and Trade Data

### UN Comtrade Database
**Provider:** United Nations Statistics Division.  
**URL:** https://comtradeplus.un.org/  
**Coverage:** 200+ reporter countries; 1962–present.  
**Key variables:** Bilateral trade flows by commodity (HS classification) and partner country.

Comtrade is the standard source for bilateral trade data. Mirror discrepancies (exports reported by one country differ from imports reported by the partner) are common and require careful handling.

### OECD Trade in Value Added (TiVA)
**Provider:** OECD.  
**URL:** https://www.oecd.org/sti/ind/trade-in-value-added.htm  
**Coverage:** 66 economies, 1995–2018.

TiVA decomposes trade flows to identify the domestic value-added content of exports, correcting the double-counting in gross trade flows that arises from global value chains.

### IMF Balance of Payments Statistics (BOP)
**Provider:** IMF.  
**URL:** https://data.imf.org/?sk=7A51304B-6426-40C0-83DD-CA473CA1FD52  
**Coverage:** ~190 countries; 1948–present for some.  
**Key variables:** Current account components, financial account components, reserve assets, international investment position.

---

## J.4 Labor Market Data

### ILO ILOSTAT
**Provider:** International Labour Organization.  
**URL:** https://ilostat.ilo.org/  
**Coverage:** ~200 countries; 1991–present for most.  
**Key variables:** Employment, unemployment (ILO definition), labor force participation, wages, hours worked, informal employment.

ILOSTAT provides internationally harmonized labor market statistics based on the ILO's standard definitions. Essential for cross-country comparisons of unemployment and labor market structure.

### OECD Labour Force Statistics
**Provider:** OECD.  
**URL:** https://stats.oecd.org/  
**Coverage:** OECD members; back to 1970s.  
**Key variables:** Unemployment rate, employment rate, participation rate, hours worked, unit labor costs, labor productivity, wage growth.

The OECD labor statistics are more methodologically consistent across countries than ILO data (because all OECD members follow common standards) and provide longer time series for advanced economies.

---

## J.5 Data Access and Reproducibility

Several best practices apply to the use of international data in macroeconomic research.

**Document the vintage.** Always record which version of a dataset was used (the PWT 10.01, the WDI April 2024 download). International data are revised and updated; results may not replicate with a different vintage.

**Understand the methodology.** Read the documentation accompanying each dataset before using it. The PWT, for instance, uses multiple GDP concepts that are not interchangeable; using the wrong one invalidates cross-country comparisons.

**Check for missing data patterns.** Missing data in international datasets is often non-random: conflict-affected countries, very small countries, and countries with weak statistical capacity are overrepresented among missing observations. Inferences based on available data may not generalize to the full population of countries.

**Use multiple sources.** When estimates from different datasets differ materially (e.g., IMF and World Bank GDP growth rates for the same country), investigate the source of discrepancy before proceeding. Discrepancies often reflect differences in methodology, coverage, or currency conversion.

**Cite primary sources.** When using PWT data, cite Feenstra, Inklaar, and Timmer (2015); when using the Maddison Project, cite Bolt et al. (2018). Do not cite FRED or World Bank as primary sources for data that originates elsewhere.

---

*For online access tools and APIs, see Appendix K.*
