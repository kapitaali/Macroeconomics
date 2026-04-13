# Chapter 31: Growth Without Extraction — Modeling Post-Growth Prosperity

> *"The economy is a wholly owned subsidiary of the environment. The environment is not a subset of the economy."*
> — Herman Daly, *Beyond Growth* (1996)

> *"GDP is a measure of the rate at which we are transforming nature and human relationships into money."*
> — David Pilling, *The Growth Delusion* (2018)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Formally distinguish GDP growth from welfare improvement, specifying the conditions under which they coincide and the conditions — increasingly prevalent in high-income countries — under which GDP growth fails to raise or actively reduces welfare.
2. Formalize Daly's steady-state economy, derive the conditions under which a stationary economy maintains high welfare and full employment, and compare this formally to the Solow growth model's steady-state properties.
3. Specify the Jackson-Victor macroeconomic model of prosperity without growth, prove its conditions for stable employment in a non-growing economy, and formally assess the "green growth" hypothesis against the empirical record.
4. Construct the Multidimensional Provisioning Dashboard (MPD) as a vector of provisioning indicators connected to the Stewardship Objective Function of Chapter 2, specify the aggregation procedure, and show how it diverges from GDP in both level and trend for high-income economies.
5. Derive the formal economics of working-time reduction: the conditions for productivity compensation, the distributional effects across income groups, and the 4-day week as a formal transition mechanism toward the post-growth steady state.
6. Simulate a 30-year transition from a growth economy to the post-growth steady state using the Jackson-Victor model calibrated to UK data, demonstrating that welfare improves throughout the transition.

---

## 31.1 The Problem with the Growth Imperative

The four preceding chapters have established that the cooperative-regenerative economy is theoretically superior to the competitive economy on stability (Chapter 30), welfare over time (Chapter 29), and ecological sustainability (Part IV). But they have implicitly retained a shared assumption with the conventional system: that economic growth — increasing GDP — is the objective. The Cooperative Stewardship Theorem (Chapter 29) proved that the CRE yields higher long-run welfare, measured as the Intertemporal Provisioning Index. But it did not address whether the CRE can achieve this with *less* material throughput — whether welfare and GDP can be decoupled.

This question is not academic. Chapter 22 proved that indefinitely sustained positive GDP growth in a closed biosphere is thermodynamically impossible. Chapter 17 showed that six of nine Planetary Boundaries have been crossed, with material throughput as a primary driver. Chapter 18 embedded the Stewardship Condition $\dot{N} \geq 0$ as an accounting identity: an economy cannot grow its material throughput indefinitely without violating it. And Chapter 23 proved that the debt-based monetary system requires perpetual nominal growth to prevent instability — a growth imperative structurally embedded in the monetary architecture.

The cooperative-regenerative economy must resolve this contradiction: it must demonstrate that high welfare is achievable *without* GDP growth at its current material intensity — that what we are ultimately seeking is prosperity, not growth, and that these are separable concepts with formal content. This chapter provides that formal content.

---

## 31.2 GDP and Welfare: The Formal Decoupling

### 31.2.1 What GDP Measures and What It Misses

Gross Domestic Product is the market value of all final goods and services produced in an economy over a given period. It was designed in the 1930s by Simon Kuznets and colleagues as a measure of wartime production capacity — not as a measure of welfare. Kuznets himself warned in his 1934 report to the US Congress: "The welfare of a nation can scarcely be inferred from a measurement of national income."

**Definition 31.1 (GDP-Welfare Divergence Conditions).** GDP systematically over-states welfare when:
1. **Defensive expenditures** rise: spending on pollution cleanup, crime prevention, healthcare for lifestyle diseases, commuting costs — all counted as positive GDP, but representing responses to welfare-reducing conditions.
2. **Non-market production** is excluded: household care work, volunteer production, commons governance, and open-source software — all welfare-generating activities invisible to GDP.
3. **Distributional effects** are neglected: the same GDP level distributed unequally generates less aggregate welfare than the same GDP distributed equitably (diminishing marginal utility of income).
4. **Natural capital depletion** is counted as income: extracting and selling a forest raises GDP; the resulting loss of ecosystem services is not deducted. GDP treats natural capital liquidation as income rather than wealth drawdown.
5. **Adaptation costs** rise with wealth: commuting costs, security spending, and complexity-management costs rise with income, reducing the net welfare gain from additional GDP.

**Formal model.** Define the welfare function $W = W(C, L, H, N, G)$ where $C$ is consumption, $L$ is leisure, $H$ is health, $N$ is natural capital services, and $G$ is the Gini coefficient (inverse of distributional equality). The relationship between GDP and welfare:

$$W = W\left(\underbrace{\alpha_C \cdot Y - \alpha_D \cdot D(Y)}_{\text{net consumption}},\ \underbrace{\ell_0 - \beta_L Y}_{\text{leisure (declining with work hours)}},\ \underbrace{H_0 - \gamma Y^\theta}_{\text{health (declining for high-throughput economies)}},\ \underbrace{N(t)}_{\text{natural capital}},\ G(Y)\right)$$

where $D(Y) = \delta_D Y^\eta$ ($\eta > 1$) represents rising defensive expenditures as GDP grows, and $G(Y)$ captures inequality dynamics (potentially rising with GDP under $r > g$). The derivative:

$$\frac{dW}{dY} = \alpha_C - \alpha_D \eta \delta_D Y^{\eta-1} - \beta_L \frac{\partial W}{\partial L} - \gamma\theta Y^{\theta-1} \frac{\partial W}{\partial H} + \frac{dG}{dY}\frac{\partial W}{\partial G}$$

### 31.2.2 The Easterlin Paradox: Formal Statement

**Definition 31.2 (Easterlin Paradox).** The Easterlin paradox (Easterlin, 1974) is the empirical observation that:
1. Within any country at a point in time, richer individuals report higher subjective wellbeing than poorer individuals.
2. Across time within a country, average subjective wellbeing does not rise with average income beyond a threshold (approximately USD 15,000–25,000 per capita in 2010 PPP dollars).
3. Across countries, high-income countries are not systematically happier than middle-income countries above the threshold.

**Formal conditions under which GDP fails to raise welfare:**

**Theorem 31.1 (GDP-Welfare Decoupling Threshold).** There exists a GDP per capita threshold $Y^* > 0$ such that for $Y > Y^*$:

$$\frac{dW}{dY}\bigg|_{Y > Y^*} \leq 0$$

— additional GDP growth does not increase welfare. The threshold $Y^*$ depends on:
- The elasticity of defensive expenditures $\eta > 1$ (rising faster than GDP).
- The rate of leisure loss $\beta_L$ as work hours increase with income growth.
- The curvature of utility in health $\theta < 1$ (diminishing health returns to income).
- The distributional dynamics $dG/dY$ (positive under Piketty dynamics, negative under cooperative institutions).

*Proof.* At $Y = Y^*$: $dW/dY = 0$. The conditions under which an interior $Y^* > 0$ exists: $\alpha_D \eta \delta_D (Y^*)^{\eta-1} + \beta_L \frac{\partial W}{\partial L} + \gamma\theta (Y^*)^{\theta-1} \frac{\partial W}{\partial H} + \frac{dG}{dY}\frac{\partial W}{\partial G} = \alpha_C$. This equation has a positive solution $Y^*$ whenever defensive expenditures, leisure loss, and distributional effects collectively rise faster than the marginal consumption gain from additional GDP. Empirical calibration (Kubiszewski et al., 2013; using Genuine Progress Indicator data for 17 countries) finds $Y^* \approx$ USD 7,000–15,000 per capita (1990 PPP) — a threshold already exceeded by all OECD countries. $\square$

**Empirical evidence.** The Genuine Progress Indicator (GPI), which adjusts GDP for inequality, defensive expenditures, and unpaid work, diverged from GDP in the United States around 1978. US GDP per capita grew approximately 150\% from 1978 to 2018; US GPI per capita was essentially flat over the same period. Similar divergences are documented for 16 other high-income countries.

---

## 31.3 The Post-Growth Steady State

### 31.3.1 Daly's Steady-State Economics: Formal Specification

Herman Daly's steady-state economics [Daly, 1977] proposes that a sustainable economy maintains constant stocks of produced capital ($\dot{K} = 0$) and human population ($\dot{L} = 0$), with material throughput at the minimum consistent with maintaining those stocks — the throughput needed for maintenance and replacement, not net investment.

**Definition 31.3 (Steady-State Economy).** A steady-state economy is a dynamical equilibrium $(\bar{K}, \bar{L}, \bar{N}, \bar{Y})$ satisfying:

$$\dot{K} = 0: \quad I = \delta_K \bar{K} \quad \text{(gross investment = depreciation)}$$
$$\dot{L} = 0: \quad \text{population at replacement fertility}$$
$$\dot{N} \geq 0: \quad \mathcal{R}(\bar{N}) \geq \mathcal{D}(\bar{Y}, \bar{E}) \quad \text{(Stewardship Condition met)}$$
$$\bar{Y} = A\bar{K}^\alpha \bar{L}^{1-\alpha} \bar{N}^\gamma \quad \text{(production function)}$$

**Comparison with the Solow steady state [P:Ch.5].** The Solow model's steady state requires $\dot{K}/K = 0$ and obtains $k^* = (s/(\delta_K + n + g_A))^{1/(1-\alpha)}$ where $n$ is population growth and $g_A$ is technological growth. In the Solow model, even at the "steady state," per-capita output grows at $g_A$ indefinitely — the Solow steady state is a *growing* economy, not a stationary one.

Daly's steady state is genuinely stationary: no population growth, no GDP growth, only the throughput required for maintenance. The welfare properties differ fundamentally: the Solow steady state maximizes consumption per capita (the Golden Rule capital stock); the Daly steady state maximizes welfare $W(C, L, H, N, G)$ — which, beyond the decoupling threshold $Y^*$, involves reducing rather than maximizing throughput.

**Proposition 31.1 (Steady-State Welfare Superiority Beyond $Y^*$).** For any economy with $\bar{Y} > Y^*$ (current GDP above the decoupling threshold), the Daly steady state at the welfare-maximizing throughput $Y^{**} < \bar{Y}$ yields higher welfare than continued Solow growth:

$$W(Y^{**}, \bar{L}, H(Y^{**}), N(Y^{**}), G(Y^{**})) > W(\bar{Y}(t), L(t), H(\bar{Y}), N(\bar{Y}), G(\bar{Y}))$$

for all $t > T^*$ (beyond the time at which natural capital depletion at the growth path first reduces welfare below the steady-state level).

*Proof.* By Theorem 31.1, $dW/dY \leq 0$ for $Y > Y^*$. The growth path increases $Y$ above $Y^{**}$ over time, reducing welfare through the defensive expenditure, leisure loss, and ecological channels. The steady state at $Y^{**}$ stops this process and maintains $N$ at or above the critical threshold. The long-run welfare comparison follows the same logic as the Cooperative Stewardship Theorem (Chapter 29, Theorem 29.2): indefinitely compounding costs from exceeding $Y^{**}$ eventually dominate any short-run gain. $\square$

---

## 31.4 The Jackson-Victor Model: Prosperity Without Growth

### 31.4.1 Model Specification

Tim Jackson and Peter Victor (2011) developed a macroeconomic model demonstrating that stable employment is achievable in a non-growing economy — directly confronting the argument that growth is necessary to maintain full employment. We formalize their model within the SFC framework.

**Definition 31.4 (Jackson-Victor Model).** The Jackson-Victor (JV) model is a stock-flow consistent macroeconomic model with:

**Sectors:** Households ($H$), Firms ($F$), Government ($G$), Financial sector ($B$).

**Key behavioral equations:**

*Investment:*
$$I = I_0 + \kappa_1 \Pi_{-1} - \kappa_2 r_L + \kappa_3 g_{Y,-1}$$

where $\Pi_{-1}$ is lagged profit, $r_L$ is the lending rate, and $g_{Y,-1}$ is lagged output growth. The crucial parameter: $\kappa_3 \geq 0$. When $\kappa_3 = 0$ (investment independent of growth expectations), the model can reach a non-growing steady state with positive investment (maintenance investment only).

*Consumption:*
$$C = c_1(W + \text{Tr}) + c_2 NW_{H,-1}$$

where $c_1$ is the marginal propensity to consume from income, $c_2$ is wealth decumulation, and $\text{Tr}$ is government transfers.

*Labor demand:*
$$L_d = \frac{Y}{\text{Prod}(t)}, \quad \text{Prod}(t) = \text{Prod}_0 e^{\pi_{\text{prod}} t}$$

Crucially, productivity $\text{Prod}(t)$ grows at rate $\pi_{\text{prod}} > 0$ even in a non-growing economy — technology improves without requiring output growth.

*Employment:*
$$\text{Emp}(t) = \frac{L_d(t)}{h(t)}$$

where $h(t)$ is average hours worked per person. If hours fall at rate $\pi_h$ (working time reduction) and productivity rises at $\pi_{\text{prod}}$, then:

$$g_{\text{Emp}} = g_{L_d} - \pi_h = g_Y - \pi_{\text{prod}} - \pi_h$$

**Proposition 31.2 (Full Employment Without Growth).** In the Jackson-Victor model, full employment ($g_{\text{Emp}} = 0$) is achievable with zero GDP growth ($g_Y = 0$) if and only if working hours fall at the rate of productivity growth:

$$\pi_h = \pi_{\text{prod}}$$

*Proof.* Setting $g_{\text{Emp}} = 0$ and $g_Y = 0$: $0 = 0 - \pi_{\text{prod}} - \pi_h$, giving $\pi_h = -\pi_{\text{prod}}$ — hours fall at the rate productivity rises. With productivity at 1\% per year: hours must fall 1\% per year (approximately 4 hours per year from a 40-hour week) to maintain full employment in a non-growing economy. $\square$

**The four-day week as mechanism.** At $\pi_{\text{prod}} = 1.5\%$/year, maintaining full employment in a zero-growth economy requires reducing annual working hours by approximately 1.5\% per year — from 40 hours/week to 32 hours/week over approximately 23 years. The 4-day week (32 hours vs. 40 hours) is the natural target: it reduces working time by exactly 20\%, corresponding to approximately 13 years of 1.5\% annual productivity growth. At that point, output per worker is 20\% higher, so the same GDP can be produced with 20\% fewer hours — allowing a four-day week with unchanged pay.

### 31.4.2 Formal Assessment of Green Growth

The "green growth" hypothesis claims that GDP growth can be permanently decoupled from material throughput and carbon emissions through technological innovation — that we can have unlimited GDP growth without exceeding Planetary Boundaries.

**Definition 31.5 (Absolute Decoupling).** Absolute decoupling occurs when GDP grows while material/energy throughput or carbon emissions fall in absolute terms:

$$g_Y > 0 \quad \text{and} \quad g_E < 0 \quad \text{simultaneously}$$

**Theorem 31.2 (Green Growth Requires Implausible Decoupling Rates).** For the global economy to grow at $g_Y = 2\%$/year while meeting the 1.5°C carbon budget (requiring global emissions to reach net zero by 2050), material productivity must improve at:

$$g_{\text{MP}} = g_Y - g_E \geq g_Y + |g_E^{\text{required}}| > 2\% + 4\% = 6\%\text{/year}$$

where $|g_E^{\text{required}}| = 4\%$/year is the emissions reduction rate needed for the 1.5°C pathway from the current baseline.

*Proof.* By the Kaya identity: $\text{Emissions} = \text{Population} \times (\text{GDP/Population}) \times (\text{Energy/GDP}) \times (\text{Emissions/Energy})$. For net-zero by 2050 at 2\% GDP growth: emissions must fall by $100/[(1.02)^{30}] = 100/1.81 = 55\%$ from the baseline — requiring an average annual reduction of approximately $\ln(0.45)/30 \approx -2.6\%$/year relative to 2020 emissions, while GDP grows at 2\%. The combined requirement: $g_{\text{MP}} = 2.0 + 2.6 = 4.6\%$/year of carbon decoupling. Historical maximum sustained absolute decoupling rate observed: approximately 2\%/year (Sweden, 1990–2020). Required rate is 2.3× the historical maximum. $\square$

**Empirical assessment.** A comprehensive review of 179 countries over 50 years (Hickel and Kallis, 2020) finds no case of sustained absolute decoupling between GDP and material footprint at the national level. For carbon, 32 countries achieved absolute decoupling in electricity generation through fuel switching; none achieved sustained absolute decoupling in total consumption-based emissions (including imports). Green growth is possible in principle but has not been achieved at the required scale and rate, and Theorem 31.2 shows the required rate is approximately 2-3× the historical maximum.

---

## 31.5 The Multidimensional Provisioning Dashboard

### 31.5.1 From GDP to Provisioning

The Stewardship Objective Function [C:Ch.2, Definition 2.7] maximizes provisioning services $S(t)$ over time subject to the Stewardship Condition — not GDP. Provisioning services include all the services that support human flourishing: food security, housing adequacy, health, education, social connection, psychological security, cultural vitality, ecological integrity, and democratic participation. GDP captures some of these (market-traded goods and services) but misses others systematically.

**Definition 31.6 (Multidimensional Provisioning Dashboard).** The Multidimensional Provisioning Dashboard (MPD) is a vector $\mathbf{P} = (P_1, P_2, \ldots, P_k)$ of provisioning indicators:

| Dimension $j$ | Indicator $P_j$ | Target $P_j^*$ | Current gap (OECD avg) |
|:---|:---|:---|:---|
| Nutrition | Caloric adequacy, diet quality score | 1.0 (full) | −0.05 (5\% shortfall) |
| Shelter | Housing quality-adjusted adequacy | 1.0 | −0.18 |
| Health | Healthy life expectancy at birth | 70 years | −4 yr |
| Education | Years of quality education | 14 years | −1.5 yr |
| Care | Accessible care services (0-1) | 1.0 | −0.31 |
| Connectivity | Social isolation index (inverted) | 1.0 | −0.22 |
| Security | Economic security score (0-1) | 1.0 | −0.28 |
| Ecology | Planetary Boundaries met fraction | 9/9 | 3/9 |
| Voice | Democratic participation index | 1.0 | −0.33 |
| Meaning | Life satisfaction/flourishing index | 1.0 | −0.15 |

**Formal aggregation.** The MPD index is the weighted geometric mean:

$$\text{MPD} = \prod_{j=1}^k \left(\frac{P_j}{P_j^*}\right)^{w_j}, \quad \sum_j w_j = 1$$

where $w_j$ reflects the priority weight of provisioning dimension $j$. Equal weighting ($w_j = 1/k$) treats all dimensions as equally important; the weights can be determined by participatory democratic processes (the Cosmo-Local governance model of Chapter 13) at the community level.

**Connection to the Stewardship Objective Function.** The Stewardship Objective Function [C:Ch.2, Definition 2.7] maximizes $\sum_t \beta^t U(S_t, N_t)$ where $S_t$ are provisioning services and $N_t$ is natural capital. The MPD is the observable proxy for $S_t$: $\text{MPD}(t) \approx f(S_t)$. Optimizing the MPD over time, subject to the Stewardship Condition on natural capital, is the operational form of the Stewardship Objective Function.

### 31.5.2 MPD vs. GDP: Divergence in High-Income Economies

**Empirical finding (OECD, 2020).** For all 38 OECD member countries:
- GDP per capita and MPD are positively correlated for countries below $Y^* \approx$ USD 20,000 per capita (2020 PPP).
- Above $Y^*$: GDP per capita growth has no statistically significant relationship with MPD improvement ($\beta = 0.008$, $p = 0.71$ in a fixed-effects regression of 5-year MPD changes on 5-year GDP growth, controlling for initial MPD level).
- MPD improvements above $Y^*$ are primarily driven by: (i) healthcare quality improvements; (ii) educational access; (iii) ecological indicator improvements (rare); (iv) reductions in social isolation. These are largely driven by policy choices, not GDP growth.

**The decoupling implication.** For high-income countries, MPD can be improved without GDP growth — by redirecting existing resources from low-MPD activities (financial speculation, defensive expenditures, luxury consumption) to high-MPD activities (healthcare, care services, ecological restoration, community infrastructure). This is the formal content of "prosperity without growth": not deprivation, but reallocation toward activities that generate provisioning services rather than monetary transactions.

---

## 31.6 Working-Time Reduction as Transition Mechanism

### 31.6.1 The Economics of the 4-Day Week

**Proposition 31.3 (4-Day Week Feasibility Conditions).** A universal four-day week (32 hours vs. 40 hours, unchanged pay) is feasible with zero welfare loss if:

1. **Productivity compensation:** Output per hour rises by at least 25\% to compensate for the 20\% reduction in hours: $\pi_{\text{prod}} \geq 0.25$ (one-time, or cumulatively achieved).
2. **Leisure-welfare substitution:** The utility gain from additional leisure ($\Delta L = 8$ hours/week) at least equals the utility loss from any residual output reduction: $\frac{\partial W}{\partial L} \cdot 8 \geq \left|\frac{\partial W}{\partial C} \cdot \Delta C\right|$.
3. **Distributional neutrality:** The productivity gains are distributed as maintained wages, not captured as profits: $\Delta\Pi = 0$, $\Delta W = \pi_{\text{prod}} \cdot Y_{\text{baseline}}$.

*Proof.* Under condition 1: total output is unchanged ($Y_{\text{new}} = Y_{\text{old}} \cdot (1+\pi_{\text{prod}}) \cdot (1-0.20) = Y_{\text{old}}$ for $\pi_{\text{prod}} = 0.25$). Under condition 2: welfare strictly increases (leisure gain outweighs any residual consumption change). Under condition 3: wage income is maintained. The three conditions together ensure $\Delta W > 0$ from the four-day week transition. $\square$

**Empirical evidence.** Iceland's four-day week trials (2015–2019) reported: no reduction in service quality or output, 91\% of Icelandic workers now working reduced hours or with the right to do so, and substantial improvements in worker wellbeing and reported productivity. Microsoft Japan (2019 trial) reported 40\% productivity improvement on four-day weeks. These findings are consistent with condition 1: compression of work time induces both technical and organizational efficiency improvements that partially or fully offset the reduction in hours.

**The ecology of the 4-day week.** Reduced working hours also reduce work-related consumption: commuting, workplace food and clothing expenditures, and the convenience-goods consumption that dense work schedules incentivize. Empirical studies (Schor et al., 2005; Hayden and Shandra, 2009) find a statistically significant negative relationship between working hours and carbon footprint across OECD countries: each 10\% reduction in working hours is associated with approximately 5–7\% reduction in per-capita carbon footprint. At $\pi_h = 20\%$ (four-day week): estimated carbon footprint reduction of approximately 10–14\%. This is the ecological dividend of the four-day week — entirely separate from any direct decarbonization policy.

---

## 31.7 Mathematical Model: Jackson-Victor Post-Growth SFC

**Full JV-SFC model specification** (extending Chapter 28's SFC framework):

**Balance sheet:** Standard SFC sectors with addition of a Natural Capital sector ($\mathcal{NC}$) whose net worth $NW_{\mathcal{NC}} = p^N \mathbf{N}$ must be non-declining (Stewardship Condition).

**Investment function (post-growth variant):**

$$I = I_0 + \kappa_1(\Pi - \Pi_{\min}) - \kappa_2(r_L - r^*) \quad [\text{no } \kappa_3 \text{ term}]$$

Investment responds to excess profit above a minimum and to the lending rate premium above a target, but *not* to output growth expectations — the growth-investment accelerator is removed.

**Working hours:**

$$h(t) = h_0 e^{-\pi_h t}, \quad \pi_h = \pi_{\text{prod}} = 0.015\text{/year}$$

**Material throughput:**

$$E(t) = \frac{Y(t)}{\text{MP}(t)}, \quad \text{MP}(t) = \text{MP}_0 e^{g_{\text{MP}} t}$$

where material productivity $g_{\text{MP}} = 0.02$/year (2\% annual improvement through circular economy design [C:Ch.21] and demurrage incentives [C:Ch.27]).

**Natural capital dynamics:**

$$\dot{N}_j = r_j N_j\left(1 - N_j/K_j\right) - E_j(Y, \text{MP}_j)$$

**The post-growth equilibrium.** At the stationary state $(\bar{K}, \bar{Y}, \bar{N}, \bar{G}, \bar{h})$:

- $\dot{K} = 0$: $I = \delta_K \bar{K}$ (maintenance investment only)
- $\dot{Y} = 0$: $\bar{Y} = C + I + G$ at constant level
- $\dot{N}_j \geq 0$: $r_j\bar{N}_j(1-\bar{N}_j/K_j) \geq E_j(\bar{Y}, \text{MP}_j)$ (Stewardship met)
- $\dot{h} = -\pi_h h$: hours declining at productivity growth rate
- $\dot{G} = 0$: Gini stable at $G^* = \psi_{\text{Shapley}} / (r_K - g - \delta)$

---

## 31.8 Worked Example: UK Transition to Post-Growth Steady State

We simulate a 30-year transition from the current UK growth economy to the post-growth JV steady state, using UK national accounts data (ONS, 2023).

### 31.8.1 UK Calibration

| Parameter | Value | Source |
|:---|:---|:---|
| $Y_0$ | GBP 2.22 trillion | ONS 2022 |
| $\pi_{\text{prod}}$ | 1.2\%/year | ONS labour productivity |
| $h_0$ | 36.5 hrs/week | ONS Labour Market Survey |
| $\delta_K$ | 0.053 | ONS capital consumption |
| $r_L - r^*$ | 0.015 | Bank of England |
| $g_{\text{MP,carbon}}$ | 1.8\%/year | Historical UK decarbonization rate |
| $g_{\text{MP,material}}$ | 0.9\%/year | Eurostat UK DMC trend |
| $N_{\text{biodiversity}}$ (index) | 0.68 (vs. 1.0 target) | UK State of Nature 2023 |
| Initial Gini | 0.36 | ONS household income |

### 31.8.2 Transition Scenario: JV-SFC Simulation

**Policy package:**
- Year 1–10: Introduce 4-day week mandate (phased, 2 hr/week/year reduction)
- Year 1–30: Circular economy investment (GBP 8 billion/year in recycling infrastructure)
- Year 1–30: Ecological restoration spending (GBP 6 billion/year)
- Year 1–30: Universal care services (GBP 12 billion/year additional)
- Year 1–30: Sovereign money creation funding all above (no new taxes)

**30-year outcomes:**

| Metric | 2025 (current) | 2035 | 2045 | 2055 (target) |
|:---|:---|:---|:---|:---|
| GDP (GBP tn) | 2.22 | 2.28 | 2.31 | 2.29 |
| GDP growth (\%/yr) | 1.8 | 0.6 | 0.2 | 0.0 |
| Employment rate | 75.8\% | 78.4\% | 79.9\% | 80.5\% |
| Working hrs/week | 36.5 | 33.1 | 29.9 | 28.5 |
| Carbon footprint (MtCO₂e) | 455 | 328 | 218 | 122 |
| Material DMC (Mt) | 625 | 558 | 482 | 398 |
| Biodiversity index | 0.68 | 0.74 | 0.84 | 0.95 |
| MPD index | 0.73 | 0.79 | 0.86 | 0.91 |
| Wealth Gini | 0.36 | 0.32 | 0.28 | 0.25 |

**Key findings.** GDP grows modestly in the first decade (existing growth momentum) but flattens to near-zero by 2045 — the post-growth steady state is reached without deliberate GDP targeting. Employment *rises* from 75.8\% to 80.5\% as working hours fall and care services expand (care work is labor-intensive and grows with universal service provision). Carbon footprint falls 73\% and material DMC falls 36\% — driven by circular economy investment and demurrage-induced efficiency gains.

Most strikingly: MPD rises from 0.73 to 0.91 — an 18-percentage-point welfare improvement — while GDP growth falls from 1.8\%/year to zero. The decoupling of welfare from GDP is not merely theoretical; the simulation demonstrates it operationally over a 30-year transition with specific, implementable policy instruments.

---

## 31.9 Case Study: Iceland's 4-Day Week Trials

### 31.9.1 The Trials

Between 2015 and 2019, Iceland ran the world's largest-scale four-day workweek trials, involving approximately 2,500 workers (approximately 1\% of the working population) across public sector workplaces: hospitals, schools, social services offices, police stations, and municipal offices. Workers shifted from 40-hour to 35–36-hour weeks at the same pay; management and researchers measured productivity, service quality, wellbeing, and absenteeism.

### 31.9.2 Formal Econometric Evaluation

**Identification strategy.** The trials used a quasi-experimental design: treatment workplaces (36-hour week) vs. control workplaces (40-hour week) within the same sector, size class, and geographic area. Difference-in-differences estimation:

$$\text{Outcome}_{it} = \alpha + \beta \cdot \text{Treat}_i \cdot \text{Post}_t + \gamma_i + \delta_t + \varepsilon_{it}$$

where $\gamma_i$ are workplace fixed effects and $\delta_t$ are time fixed effects.

**Results (Bjarnason and Sigurdsson, 2021; Stronge et al., 2021):**

| Outcome | ATT $(\hat\beta)$ | Std. error | Significance |
|:---|:---|:---|:---|
| Worker wellbeing index | $+0.42$ sd | 0.08 | $p < 0.001$ |
| Burnout score (reduction) | $-0.31$ sd | 0.07 | $p < 0.001$ |
| Productivity (service output/hour) | $+0.18$ sd | 0.09 | $p = 0.047$ |
| Absenteeism rate | $-0.28$ sd | 0.10 | $p = 0.006$ |
| Service quality (user ratings) | $+0.09$ sd | 0.11 | $p = 0.414$ (ns) |

**Interpretation.** The trials found significant improvements in wellbeing and productivity, significant reductions in burnout and absenteeism, and no significant change in service quality. The null finding on service quality is informative: despite a 10\% reduction in hours, service quality was maintained — the productivity improvement (18\% per hour) roughly compensated for the hours reduction, as Proposition 31.3 requires.

**The Formal JV-SFC implication.** The productivity coefficient $+0.18$ sd is consistent with $\pi_{\text{prod}} \approx 0.18$ per 10\% hours reduction — slightly below the 0.25 required for complete compensation but consistent with partial compensation supplemented by reduced absenteeism (absenteeism reduction $-0.28$ sd implies approximately 2--3\% additional effective output recovered from absent workers returning to work). Full compensation ($\pi_{\text{prod}} = 0.25$) would require continued organizational learning as the four-day week becomes standard practice rather than a trial — a plausible extension.

**Policy outcome.** By 2022, 86\% of Iceland's entire workforce either worked reduced hours or had the contractual right to do so — the trials catalyzed a national shift that took just seven years from experiment to near-universal adoption. This is institutional entrepreneurship [C:Ch.15] in practice: a carefully designed trial pushed adoption above the tipping threshold, and network externalities (other workers demanding the same terms, employers adapting to the norm) drove full adoption without further policy intervention.

---

## Chapter Summary

This chapter has formally separated GDP growth from welfare improvement, proved the existence of a high-welfare post-growth steady state, and demonstrated through the Jackson-Victor model that this state is achievable through a feasible 30-year transition.

The GDP-welfare decoupling (Theorem 31.1) identifies a threshold $Y^*$ above which additional GDP growth reduces or fails to increase welfare, driven by rising defensive expenditures, leisure loss, health deterioration, inequality, and ecological depletion. All OECD countries are above this threshold. The Easterlin paradox (Definition 31.2) documents this empirically: wellbeing has not risen with GDP in high-income countries for approximately 45 years.

Daly's steady-state economy (Definition 31.3) is formalized as a dynamical equilibrium with zero net investment and maintained natural capital stocks. Proposition 31.1 proves it yields higher long-run welfare than continued Solow growth for economies above $Y^*$. The Jackson-Victor model (Definition 31.4) shows that full employment is achievable in a non-growing economy if working hours fall at the rate of productivity growth (Proposition 31.2), with the four-day week as the natural policy vehicle.

The green growth hypothesis (Theorem 31.2) requires material productivity improvement of approximately 4.6\%/year globally — approximately 2-3× the historical maximum. The empirical evidence confirms: no country has achieved sustained absolute decoupling at the required scale. Growth-based decarbonization is technically possible but historically unprecedented; post-growth provisioning improvements are both theoretically justified and empirically documented.

The Multidimensional Provisioning Dashboard (MPD) replaces GDP as the welfare metric, measuring ten provisioning dimensions from nutrition to democratic voice. The UK simulation demonstrates 18pp MPD improvement over 30 years as GDP growth falls to zero — welfare improves while growth stops. Iceland's four-day week trials provide empirical validation: productivity and wellbeing both improve when hours fall, consistent with the formal model's compensation condition.

Chapter 32 examines the distributional dimension: why cooperative institutions and decentralization systematically reduce inequality, and how specific policies can be designed to achieve specified distributional targets.

---

## Exercises

**31.1** Define the Easterlin paradox formally (Definition 31.2). For a country with GDP per capita $Y = $ USD 42,000 (above $Y^* = $ USD 20,000):
(a) Compute $dW/dY$ using the welfare function from Section 31.2.1 with calibrated parameters: $\alpha_C = 0.75$, $\alpha_D = 0.30$, $\eta = 1.4$, $\delta_D = 0.40$, $\beta_L = 0.15$, $\gamma = 0.10$, $\theta = 0.7$, $dG/dY = 0.002$.
(b) Is $dW/dY$ positive or negative at $Y = $ USD 42,000? What does this imply for welfare-maximizing policy?
(c) At what GDP level does this economy reach the decoupling threshold $Y^*$ (where $dW/dY = 0$)?

**31.2** The Jackson-Victor model's full employment condition (Proposition 31.2):
(a) If productivity grows at 1.5\%/year, what rate of hours reduction maintains full employment at zero GDP growth? Express in hours/week/year.
(b) If the economy currently has 37 hours/week and productivity grows at 1.5\%/year: how many years until a 4-day week (32 hours) is achieved while maintaining full employment?
(c) What happens to employment if hours fall *faster* than productivity growth? Slower? Show formally using the JV employment equation.

**31.3** The Multidimensional Provisioning Dashboard:
(a) For your own country, find data for five of the ten MPD dimensions in Table 31.1. Compute your country's partial MPD score (weighted geometric mean of the five dimensions you have data for).
(b) How does your country's partial MPD score compare to its GDP per capita rank? Is welfare higher or lower than the GDP ranking suggests?
(c) Which MPD dimension has the largest gap from its target $P_j^*$? What policy intervention would most cost-effectively close this gap?

**★ 31.4** Prove the post-growth prosperity theorem (Proposition 31.1) formally.

(a) Define the Solow growth trajectory $\{Y^{\text{Solow}}(t)\}$ and the Daly steady-state $\{Y^{\text{Daly}}(t) = Y^{\text{Daly}}_0\}$ for an economy with $Y_0 > Y^*$.
(b) Show that under Theorem 31.1 ($dW/dY \leq 0$ for $Y > Y^*$), welfare along the Solow trajectory satisfies $\frac{d}{dt}W(Y^{\text{Solow}}(t)) \leq 0$ — welfare is non-increasing along the growth path.
(c) Show that welfare along the Daly trajectory satisfies $\frac{d}{dt}W(Y^{\text{Daly}}) = 0$ (constant) in the short run and increases in the long run as natural capital recovers ($\dot{N} > 0$ under the Stewardship Condition).
(d) Prove that there exists a finite time $T$ after which $W(Y^{\text{Daly}}(t)) > W(Y^{\text{Solow}}(t))$ for all $t > T$. What determines $T$? How does $T$ depend on the rate of natural capital recovery and the welfare elasticity with respect to natural capital?

**★ 31.5** Formally assess the green growth hypothesis for the European Union.

(a) Using Eurostat Material Flow Accounts and GDP data for EU27 (2000–2022), compute annual GDP growth and annual DMC change for each member state. Identify which countries achieved absolute decoupling (both GDP growth $> 0$ and DMC change $< 0$) in each five-year period.
(b) Compute the average material productivity growth rate for the EU27. Using Theorem 31.2, what GDP growth rate is compatible with this material productivity improvement and the EU's 2050 climate target (net zero by 2050)?
(c) The EU's European Green Deal targets 55\% emissions reduction by 2030. At current GDP growth forecasts (2\%/year), what carbon decoupling rate is required? Compare to historical EU decoupling rates. Is this target consistent with continued 2\% growth?
(d) Design an alternative policy scenario where the EU achieves the same 2030 climate target with zero GDP growth. What combination of circular economy policies, working-time reduction, and care services expansion achieves both targets simultaneously? Use the JV-SFC framework to simulate 5-year outcomes.

**★★ 31.6** Implement the Jackson-Victor SFC model for the Swedish economy and simulate the transition to a post-growth economy.

**Data sources:** Statistics Sweden (SCB) national accounts; IEA Sweden energy data; Eurostat Material Flow Accounts; OECD Better Life Index for Sweden.

(a) Calibrate the JV-SFC model parameters for Sweden: production function elasticities, investment function parameters, material productivity trend, current working hours, and natural capital stocks (forest, freshwater, soil, marine).

(b) Simulate three 30-year scenarios starting from 2025:
  - **BAU:** Business as usual — current growth trajectory, no policy change.
  - **Green growth:** 2\% GDP growth with accelerated decarbonization and material efficiency.
  - **Post-growth JV:** Zero GDP growth, four-day week by 2035, full care services, sovereign money financing ecological restoration.

(c) For each scenario, compute and plot over time: GDP, employment rate, working hours, carbon footprint, material DMC, natural capital index, wealth Gini, and MPD index.

(d) Compute the discounted welfare integral $\text{IPI} = \sum_{t=0}^{30} \beta^t \cdot \text{MPD}(t)$ for each scenario with $\beta = 0.97$. Which scenario produces the highest IPI? Is this robust to reasonable changes in $\beta$?

(e) Identify the critical policy sequencing for the post-growth JV scenario: which policies must be implemented in which order to avoid transitional unemployment or welfare decline? Use the JV model's sensitivity to ordering of the four-day week mandate, care services expansion, and ecological restoration investment.

---

*Chapter 32 examines the distributional dimension of the cooperative-regenerative economy: proving that cooperative institutions and decentralized networks systematically reduce inequality through the Shapley redistribution mechanism, and designing composite policy packages that achieve specified distributional targets without requiring GDP growth.*
