# Chapter 3 — Key Concepts: GDP, Inflation, and Unemployment

> *"Not everything that can be counted counts, and not everything that counts can be counted."*
> — Attributed to William Bruce Cameron

---

Before we can analyze how the economy works, we need to agree on how to measure it. This is less straightforward than it might appear. Economic activity is an extraordinarily complex phenomenon — billions of transactions, in millions of distinct goods and services, across all sectors of a modern economy — and reducing it to a handful of summary statistics inevitably involves choices, approximations, and conceptual compromises. Understanding what those choices are, and what they imply for the interpretation of economic data, is an indispensable foundation for everything that follows.

This chapter develops three measurement concepts in depth: gross domestic product (GDP), the inflation rate, and the unemployment rate. These are not merely accounting conventions; they are the primary empirical referents for all the theoretical models developed in subsequent chapters. When a model predicts that "output rises," we need to know what output means. When a chapter discusses "the inflation rate," we need to know which price index is being measured and why the choice matters. When Chapter 31 analyzes unemployment policy, we need to know what the unemployment rate is actually counting.

---

## 3.1 Gross Domestic Product: Definition and Measurement

GDP is the most widely cited single statistic in macroeconomics. But it is also frequently misunderstood, both in terms of what it measures and in terms of its limitations as a measure of welfare.

**Definition (Gross Domestic Product).** Real **gross domestic product (GDP)** $Y_t$ is the total market value of all **final** goods and services produced within a country's geographic borders during period $t$, expressed at constant (base-year) prices. Three qualifications are embedded in this definition: "total market value" means goods are valued at transaction prices, not at social benefit; "final" means the good reaches its ultimate user without further transformation in production; and "within a country's borders" means the location of production, not the nationality of producers.

### The Three Approaches and Their Equivalence

The same aggregate can be measured in three logically distinct but numerically equivalent ways, corresponding to different points of observation in the circular flow of income and expenditure.

The **expenditure approach** decomposes GDP by the type of spending:

$$Y \equiv C + I + G + NX,$$

where $C$ is private consumption (household purchases of goods and services), $I$ is gross private investment (business fixed investment, residential investment, and inventory changes), $G$ is government consumption and gross investment, and $NX = X - M$ is net exports (exports minus imports). This is an accounting identity — it holds by definition, not as a behavioral prediction.

The **income approach** decomposes GDP by the type of factor income generated in production:

$$Y \equiv W + \Pi + R + T^{ind} + D,$$

where $W$ is compensation of employees (wages, salaries, benefits), $\Pi$ is gross operating surplus (corporate profits plus proprietors' income), $R$ is gross mixed income, $T^{ind}$ is taxes less subsidies on production, and $D$ is the consumption of fixed capital (depreciation). This identity holds because every dollar of output generates exactly one dollar of income for someone — either as labor income, capital income, or taxes.

The **production approach** sums **value added** across all industries:

$$Y \equiv \sum_{j=1}^J (Q_j - M_j),$$

where $Q_j$ is the gross output of industry $j$ and $M_j$ is its intermediate inputs. The reason for summing value added rather than gross output is to avoid double-counting: the steel that goes into a car is not counted separately from the car itself, because the value of the steel is already embedded in the value of the car. Measuring only value added at each stage ensures that each unit of final output is counted exactly once.

### Real Versus Nominal GDP

Nominal GDP at date $t$ is $P_t Y_t$ — the product of the price level and real output. It rises whenever either real output rises or prices rise. To isolate changes in productive activity from changes in prices, we construct real GDP by holding prices constant at a base-year level.

**Definition (GDP Deflator).** The **GDP deflator** is the ratio of nominal to real GDP, expressed as an index:

$$\text{GDP Deflator}_t = \frac{\text{Nominal GDP}_t}{\text{Real GDP}_t} \times 100 = \frac{\sum_i p_{it} q_{it}}{\sum_i p_{i0} q_{it}} \times 100,$$

where $p_{it}$ and $q_{it}$ are the price and quantity of good $i$ in period $t$, and $p_{i0}$ is its base-year price. The deflator is not a consumer price index — it covers all goods and services produced in the economy, not just those consumed by households — but it serves the same purpose of stripping out price changes.

The GDP deflator above is a **Laspeyres quantity index** (fixed base-year prices applied to current quantities). An alternative is the **Paasche quantity index** (fixed current-period prices applied to base-year quantities). Neither alone is ideal because both are subject to **substitution bias**: as relative prices change, agents substitute toward cheaper goods, but fixed-weight indices do not capture this substitution. The **Fisher ideal index** addresses this by taking the geometric mean of the Laspeyres and Paasche indices:

$$P_F = \sqrt{P_L \cdot P_P}.$$

The Bureau of Economic Analysis constructs U.S. real GDP using a **chain-weighted Fisher index**, which applies the Fisher formula sequentially to adjacent periods and chains the results together. This approach minimizes substitution bias and avoids the distortions introduced by using a single distant base year.

### What GDP Does Not Measure

GDP is not a measure of welfare, and treating it as one leads to systematic errors of interpretation. Several important dimensions of welfare are excluded. Household production — cooking, childcare, cleaning, and maintenance performed within the home — is not counted because it does not pass through market transactions, even though it creates genuine value. The value of leisure is not counted. The depletion of natural capital — the extraction of non-renewable resources or the degradation of ecosystems — is counted as positive production without any offsetting subtraction for the loss of the underlying asset. The distribution of income is invisible in GDP: a country in which all income accrues to 1% of the population and a country with equal distribution can have identical GDPs. These limitations do not make GDP useless — it remains the best available single-number summary of productive activity — but they define the scope of what it measures.

---

## 3.2 The Consumer Price Index and Inflation Measurement

Inflation, defined in Chapter 1 as the rate of change of the price level, is in practice always measured relative to a specific price index. The choice of index is not neutral: different indices answer different questions and are subject to different measurement biases.

The most widely reported price index in most countries is the **Consumer Price Index (CPI)**. It measures the cost, in the current period, of purchasing a fixed basket of goods and services that represents the typical consumption of households in a base period.

**Definition (Consumer Price Index).** The CPI at date $t$ with base-period basket $\mathbf{q}_0$ is a Laspeyres price index:

$$\text{CPI}_t = \frac{\sum_i p_{it} q_{i0}}{\sum_i p_{i0} q_{i0}} \times 100,$$

where $q_{i0}$ is the quantity of good $i$ in the base-period basket and $p_{it}$ is its current-period price. The inflation rate is:

$$\pi_t = \frac{\text{CPI}_t - \text{CPI}_{t-1}}{\text{CPI}_{t-1}}.$$

The Laspeyres formula uses base-period quantities as weights. This introduces a systematic upward bias in measured inflation through three channels. First, **substitution bias**: as relative prices change, consumers shift toward goods that have become relatively cheaper, but the fixed basket assumes no substitution. If apples become more expensive relative to oranges, households buy more oranges, but the CPI continues to measure the cost of the old apple-heavy basket. Second, **outlet bias**: as prices rise, consumers shift to discount retailers, but the CPI does not fully account for this. Third, **quality and new goods bias**: price increases partly reflect genuine quality improvements, and new goods (smartphones, streaming services) may not enter the basket for years after their introduction. The Boskin Commission (1996) estimated total U.S. CPI bias at approximately 1.1 percentage points per year — a substantial figure implying that measured inflation persistently overstates the true increase in the cost of living.

Because food and energy prices are highly volatile — subject to large swings due to weather, geopolitical events, and commodity market dynamics — central banks and analysts often focus on **core inflation**, which excludes these two categories:

$$\pi_t^{\text{core}} = \frac{\text{CPI}_t^{\text{core}} - \text{CPI}_{t-1}^{\text{core}}}{\text{CPI}_{t-1}^{\text{core}}}.$$

Core inflation is a better predictor of future inflation than headline inflation because it filters out transitory supply disturbances, isolating the underlying trend driven by demand conditions. However, it is not a perfect signal: if energy prices rise permanently (as during the 1973 oil shock), excluding them from the inflation measure will cause the central bank to underestimate the inflationary pressure being transmitted to wages and other prices through the cost channel.

A statistically more principled approach is the **trimmed mean** inflation measure, which removes a fraction $\alpha$ from both tails of the cross-sectional distribution of individual item price changes before averaging:

$$\pi_t^{\text{trim}}(\alpha) = \frac{1}{1-2\alpha}\int_\alpha^{1-\alpha} F_t^{-1}(p)\,\mathrm{d}p,$$

where $F_t^{-1}$ is the quantile function of the cross-sectional distribution. The trimmed mean drops outliers — both unusually large price increases and unusually large price decreases — rather than dropping specific categories like food and energy. It therefore responds to persistent supply shocks that affect many items simultaneously while filtering out idiosyncratic item-level volatility.

---

## 3.3 Unemployment: Concepts, Measurement, and the Natural Rate

Of the three core macroeconomic variables, unemployment is arguably the one whose measurement is most fraught with conceptual ambiguity. The official unemployment rate is a well-defined statistic, but it understates — sometimes by a large margin — the true extent of labor market slack.

**Definition (Unemployment Rate).** Following the International Labour Organization (ILO), a person is **unemployed** if she is: (i) without work — not in paid employment or self-employment; (ii) currently available for work; and (iii) seeking work — has taken specific steps to find employment in the recent reference period. The **labor force** $L_t$ is the sum of the employed $N_t$ and the unemployed $U_t$. The **unemployment rate** is:

$$u_t = \frac{U_t}{L_t}.$$

The critical word in the definition is "seeking." Persons who have become discouraged and stopped looking for work — because they believe no jobs are available for them — are classified as **out of the labor force**, not unemployed. During severe recessions, a substantial number of workers exit the labor force in this way, causing the measured unemployment rate to understate the degree of labor market distress.

To capture this, the U.S. Bureau of Labor Statistics publishes several alternative unemployment measures. The most comprehensive is U-6:

$$u_t^{U6} = \frac{U_t + M_t + P_t}{L_t + M_t},$$

where $M_t$ is marginally attached workers (those who want work and have searched recently but not in the past four weeks, including discouraged workers) and $P_t$ is persons working part-time for economic reasons (those who want full-time work but can only find part-time work). The U-6 rate is typically 3–5 percentage points higher than the official U-3 rate and provides a more complete picture of labor market slack.

### The Natural Rate of Unemployment

A central concept in the theory of unemployment is the distinction between the unemployment that would exist even in a well-functioning economy — because matching workers to jobs takes time and some skill mismatches are inevitable — and the unemployment that reflects macroeconomic deficiency of demand.

**Definition (Natural Rate of Unemployment).** The **natural rate of unemployment** $u^*$ is the unemployment rate that would prevail if all prices and wages were fully flexible and the economy were at its potential output $\bar{Y}_t$. It is the equilibrium rate consistent with stable, non-accelerating inflation. It is sometimes called the **NAIRU** (Non-Accelerating Inflation Rate of Unemployment) to emphasize this property.

The natural rate is not zero. It consists of two components. **Frictional unemployment** is the unemployment that arises because the process of matching workers to jobs takes time: new graduates searching for their first job, workers laid off from a declining industry who need time to find positions in a growing one, and workers voluntarily quitting to seek better opportunities. Frictional unemployment is broadly compatible with labor market efficiency — some search time is socially useful because it leads to better matches. **Structural unemployment** arises from persistent mismatches between the skills workers have and the skills employers need, often resulting from technological change or geographic shifts in industry.

The natural rate is not fixed: it changes as the structure of the economy, the generosity of unemployment insurance, the efficiency of job-matching platforms, and demographic composition all evolve. In the United States, $u^*$ was estimated at approximately 6% in the 1980s, fell to around 5% in the late 1990s, and many economists now place it closer to 4%–4.5% — reflecting improved labor market matching technology (including online job boards), changing industry composition, and demographic shifts toward older (lower-turnover) workers.

The search-and-matching model (Mortensen and Pissarides, 1994) provides the microeconomic foundation for the natural rate. With matching function $m_t = m(U_t, V_t)$ exhibiting constant returns to scale, the job-finding rate for workers is $f(\theta_t) = m_t/U_t$ and the vacancy-filling rate for firms is $q(\theta_t) = m_t/V_t$, where $\theta_t = V_t/U_t$ is **labor market tightness**. In steady state:

$$u^* = \frac{\delta}{\delta + f(\theta^*)},$$

where $\delta$ is the job-destruction rate. The natural rate rises with $\delta$ (more frequent job losses) and falls with $f(\theta^*)$ (faster job-finding). Chapter 13 develops this framework in full.

---

## 3.4 The Relationships Among Output, Inflation, and Unemployment

The three core variables are not independent: they are linked by two empirical regularities that recur throughout macroeconomics.

**Okun's Law.** Arthur Okun (1962) documented that the unemployment rate and the output gap — the deviation of actual output from potential — move closely together. The relationship, known as Okun's Law, states:

$$Y_t - \bar{Y}_t = -\psi(u_t - u^*), \quad \psi > 0.$$

Verbally: when unemployment is above its natural rate by one percentage point, output is below potential by $\psi$ percent. Empirical estimates for the United States cluster around $\psi \approx 2$: a one-percentage-point rise in unemployment is associated with approximately a 2% output shortfall. This relationship holds reliably in U.S. data across many decades but is not a structural equation — the parameter $\psi$ varies across countries and over time. It is a reduced form that reflects the underlying production technology, labor market institutions, and hours-versus-employment adjustment mechanisms of each economy.

**The Phillips Curve.** The connection between inflation and unemployment (or equivalently, between inflation and the output gap via Okun's Law) is captured by the Phillips curve, developed in depth in Chapter 10. The expectations-augmented version:

$$\pi_t = \pi_t^e - \alpha(u_t - u^*) + \epsilon_t, \quad \alpha > 0,$$

links the output side of the economy (unemployment, or equivalently the output gap) to the nominal side (inflation). Together, Okun's Law and the Phillips curve form the empirical backbone of the AS–AD framework.

---

*Next: Chapter 4 — The Circular Flow and National Income Accounting*
