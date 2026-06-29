# Chapter 37 — Climate Change and Macroeconomics: Sustainability and Green Growth

> *"Climate change is the greatest market failure the world has ever seen."*
> — Nicholas Stern, *The Stern Review*, 2006

---

Climate change presents macroeconomics with its most difficult and consequential challenge: a global externality whose costs are concentrated in the distant future but whose causes are embedded in the current structure of production and consumption, and whose mitigation requires international coordination on an unprecedented scale. The macroeconomic dimensions are vast: they encompass the long-run growth implications of rising temperatures, the optimal carbon pricing that internalizes climate damages, the macroeconomic transition risks from decarbonization, and the distributional consequences of both climate impacts and mitigation policies. This chapter develops the formal tools for thinking about these questions — from the basic Pigouvian externality framework to the DICE integrated assessment model to the current debates about physical and transition risk in central bank stress tests.

---

## 37.1 Climate Change as an Externality

**Definition (Carbon Externality).** The burning of fossil fuels imposes a **negative externality** on current and future generations: greenhouse gas emissions raise atmospheric CO$_2$ concentration, which warms the climate, which imposes costs (reduced agricultural yields, higher sea levels, more extreme weather, ecosystem degradation) on parties who receive no compensation from the emitter. Because these costs are not reflected in the market price of fossil fuels, the market equilibrium produces too much carbon — more than is socially optimal.

The competitive market produces emissions $E^{market}$ where marginal private cost equals marginal private benefit: $MPC(E^{market}) = MPB(E^{market})$. The socially optimal level $E^*$ satisfies:

$$MPC(E^*) + MEC(E^*) = MPB(E^*),$$

where $MEC(E^*)$ is the **marginal external cost** — the present value of all future damages from one additional ton of CO$_2$ equivalent. Since $MEC > 0$, the socially optimal emissions level $E^* < E^{market}$.

A **Pigouvian carbon tax** at rate $\tau^* = MEC(E^*)$ corrects the externality by aligning the private and social cost. This is the first-best climate policy: a single uniform carbon price achieves the socially optimal allocation without any quantitative restrictions, technology mandates, or sector-specific regulations. All other climate policies (renewable energy subsidies, fuel economy standards, building codes, cap-and-trade systems) are evaluated against this benchmark — some approximate it; most involve additional distortions.

### Why Climate Change is Harder than a Standard Externality

Three features make climate change more difficult than the standard externality framework accommodates:

**Temporal asymmetry**: costs materialize over decades and centuries, while causes are contemporaneous. Standard present-value discounting of future costs is therefore enormously consequential — and deeply contested — since the discount rate determines whether future damages are nearly worthless (high discount rate) or nearly as important as present costs (low discount rate).

**Global public good character**: the climate is a global commons. Emission reductions by one country benefit all countries regardless of whether they reduce their own emissions. This creates a free-rider incentive: each country prefers that others bear the cost of mitigation while it enjoys the benefits. The Nash equilibrium of the non-cooperative climate game involves dramatically insufficient mitigation, requiring treaty-based international coordination that is historically difficult to achieve and maintain.

**Irreversibility and tipping points**: carbon removed from the atmosphere now was emitted over centuries; it cannot be cheaply reversed. More importantly, the climate system may have **tipping points** — bifurcations beyond which the system shifts to a qualitatively different and much warmer state (Arctic permafrost carbon release, Amazon dieback, West Antarctic Ice Sheet collapse). Near tipping points, the damage function becomes convex in a way that standard integrated assessment models do not capture.

---

## 37.2 Integrated Assessment Models: The DICE Framework

The formal quantitative framework for climate economics is the **integrated assessment model (IAM)**, linking an economic model to a physical climate model. The most influential is Nordhaus's DICE (Dynamic Integrated model of Climate and the Economy).

### The Structure of DICE

Output net of climate damages and abatement costs:

$$Y_t^{net} = A_t K_t^\alpha L_t^{1-\alpha}\cdot\Omega(T_t)\cdot(1 - \Lambda(\mu_t)),$$

where $\Omega(T_t) = 1/(1 + \pi_2 T_t^2)$ is the **damage function** reducing output as global mean temperature anomaly $T_t$ rises. With $\pi_2 \approx 0.00267$ (Nordhaus), a 3°C warming costs approximately 9% of GDP. The **abatement cost function** $\Lambda(\mu_t) = b_1\mu_t^{b_2}$ with $b_2 \approx 2.6$ captures rising marginal costs of emission reduction; $\mu_t \in [0,1]$ is the abatement rate.

The climate module:

$$E_t = (1-\mu_t)A_t\sigma_t K_t^\alpha L_t^{1-\alpha}, \quad M_{t+1} = M_t + \beta_E E_t - \beta_M(M_t - M^{pre}),$$
$$T_{t+1} = T_t + \zeta_T[\eta\ln(M_{t+1}/M^{pre}) - T_t],$$

where $\sigma_t$ is the carbon intensity of output, $M_t$ is atmospheric CO$_2$ concentration, $T_t$ is the temperature anomaly, and $\eta \approx 3$°C is the equilibrium climate sensitivity (warming from CO$_2$ doubling).

The social planner maximizes:

$$\max_{\{\mu_t, K_{t+1}\}}\; U = \int_0^\infty e^{-\rho t}\frac{c_t^{1-\sigma}-1}{1-\sigma}L_t\,\mathrm{d}t,$$

subject to the economic and climate modules. The solution paths for $\mu_t^*$ (optimal abatement) and $C_t^*$ are found by optimal control [M:Ch.11].

### The Damage Function Debate

The quadratic damage function $\Omega(T) = 1/(1+\pi_2 T^2)$ implies manageable, smooth damages. Several critiques challenge this:

**Burke, Hsiang, and Miguel (2015)** use panel data across countries to estimate the relationship between temperature and GDP growth non-parametrically, finding strong non-linearities and much larger damages than the Nordhaus calibration — approximately 23% loss of global income with 4°C of warming. Their estimates imply that the SCC is several times higher than the DICE model suggests.

**Tail risk**: the damage function is estimated from the historical variation in temperatures, which does not include scenarios of 5–6°C warming. Extrapolating is highly uncertain. Howard and Sterner (2017) show that higher damage function estimates, within the plausible range, raise the SCC by a factor of 3–4.

---

## 37.3 The Social Cost of Carbon

**Definition (Social Cost of Carbon).** The **social cost of carbon (SCC)** is the present discounted value of the marginal damage from an additional ton of CO$_2$ emitted today:

$$SCC_t = -\mathbb{E}_t\int_t^\infty e^{-\int_t^s r(\tau)\,\mathrm{d}\tau}\frac{\partial Y_s^{net}}{\partial E_t}\,\mathrm{d}s.$$

The SCC is the correct Pigouvian carbon tax and the central policy output of any IAM.

### The Stern-Nordhaus Disagreement

The SCC is exquisitely sensitive to the discount rate, generating one of the most important normative debates in economics.

**Nordhaus**: use market interest rates. With $\rho \approx 1.5\%$ pure time preference, $\sigma = 1.45$, and $g \approx 2\%$, the Ramsey rule gives $r = \rho + \sigma g \approx 4.4\%$. This reflects market revealed preferences for trading present against future consumption. The Nordhaus SCC is approximately $\$37$/tCO$_2$ (2010 dollars), rising to $\$100$ by 2050.

**Stern**: the pure rate of time preference should be near zero ($\rho \approx 0.1\%$) on the ethical ground that future people's welfare should not be discounted merely because they are in the future — future people are morally equivalent to present people. With $\sigma = 1$ and $g \approx 1.3\%$, the Stern discount rate is $r \approx 1.4\%$. The Stern SCC is approximately $\$300$/tCO$_2$ — roughly eight times the Nordhaus estimate — implying dramatically more aggressive near-term mitigation.

The core disagreement is normative, not empirical. $\rho$ is a value judgment about intergenerational equity that neither empirical economics nor cost-benefit theory alone can resolve. However, the disagreement has clear policy implications: the Nordhaus approach accepts 3–4°C of warming with modest mitigation costs; the Stern approach implies aggressive decarbonization to limit warming to 1.5–2°C at substantially higher current costs.

### Fat Tails and Catastrophic Risk

Weitzman (2009) argues that the standard expected utility framework understates optimal carbon prices because the distribution of climate damages has heavy tails: there is a small but non-negligible probability of catastrophic warming (5–6°C) that would destroy a large fraction of world output. When damages are catastrophic and uncertainty is large, risk aversion generates very high optimal carbon prices — potentially far above both Nordhaus and Stern estimates — as insurance against the worst-case tail scenarios.

---

## 37.4 Macroeconomic Risks from Climate Change

Beyond DICE's smooth damages, climate change generates macroeconomic risks that are more abrupt and harder to model.

### Physical Risks

Increases in the frequency and severity of extreme weather events (hurricanes, floods, wildfires, heat waves) destroy physical capital, disrupt supply chains, reduce agricultural yields, and generate large migration flows. Hsiang and Jina (2014) estimate that an average tropical cyclone reduces GDP growth for 20 years after impact, with the long-run GDP loss exceeding 50 times the immediate physical damage. These long-duration effects reflect the disruption of growth-promoting agglomeration economies, the destruction of human capital, and the reallocation of investment toward reconstruction rather than productive expansion.

**Macro-financial exposure**: large infrastructure investments (ports, coastal roads, power plants, real estate) have economic lives of 30–80 years. Assets built assuming historical climate conditions are exposed to physical damage risk under future climate conditions. The financial institutions holding mortgages on coastal properties or loans to agricultural firms face credit risk from physical climate damages — risk that is not fully priced in current markets because climate damages are still in the future.

### Transition Risks

A rapid shift to low-carbon energy systems will strand fossil fuel assets — coal mines, oil fields, gas pipelines, and power plants that become economically unviable before the end of their planned productive lives under a stringent carbon price. The IPCC estimates stranded assets at $1–4 trillion globally under a 1.5°C scenario. These losses — concentrated in the energy sector and in banks with significant energy lending — could trigger financial distress if the transition is sudden rather than gradual.

Central banks and financial regulators have begun conducting **climate stress tests**: scenarios in which a sudden carbon price increase (the "disorderly transition" scenario) strands fossil fuel assets and generates credit losses. The NGFS (Network for Greening the Financial System) develops standardized scenarios; the ECB's 2022 climate stress test found that under an orderly transition, bank losses are manageable (approximately €70 billion over 30 years), but under a disorderly transition, losses could be several times larger with significant systemic implications.

---

## 37.5 Carbon Pricing Instruments: Tax versus Cap-and-Trade

### The Policy Comparison

**Carbon tax**: the government sets the carbon price; the market determines the quantity of emission reductions. Provides price certainty, enabling long-term investment decisions by firms and households. Revenue can be recycled as dividends (the "carbon dividend" approach), reducing regressivity.

**Cap-and-trade (Emissions Trading System, ETS)**: the government caps total emissions and issues tradeable permits. Provides quantity certainty; the market price of permits equals the shadow price of the cap. Price uncertainty discourages long-term low-carbon investment if the cap is not binding.

**Welfare comparison (Weitzman, 1974)**: under certainty, the two are equivalent. Under uncertainty: a carbon tax is preferred when the marginal abatement cost curve is steep (price certainty more important than quantity certainty); cap-and-trade is preferred when the marginal damage curve is steep (quantity certainty more important). For climate change, where damages are highly non-linear near tipping points, this result shifts toward quantity instruments for stringent targets.

### The Political Economy of Carbon Pricing

Carbon taxes face stronger political opposition than cap-and-trade because they make the carbon price salient (fuel price spikes coincide with the tax). British Columbia's revenue-neutral carbon tax (2008, now C$65/tCO$_2) and Sweden's carbon tax (€130/tCO$_2$) are the highest in the world and have demonstrably reduced emissions while maintaining economic growth — but both required careful political packaging (revenue recycling, border adjustments, sectoral exemptions) to build coalitions. The EU ETS has achieved significant emission reductions but at volatile prices (ranging from €5 to €100/tCO$_2$ since 2005), creating investment uncertainty.

A **price corridor** — a floor and ceiling price in an ETS — combines price and quantity certainty, providing bounded carbon price predictability within a quantity commitment. This hybrid design is increasingly seen as the optimal institutional arrangement for a carbon pricing system that must persist through multiple election cycles.

---

## Chapter Summary

- The **carbon externality** means competitive markets overproduce emissions relative to the social optimum by $MEC > 0$ per ton. The first-best corrective is a Pigouvian tax $\tau^* = MEC(E^*)$ equal to the marginal external cost. Climate's global public good character, temporal asymmetry, and potential tipping points make it harder than standard externalities.

- The **DICE model** links a Ramsey growth economy to a physical climate module through a damage function $\Omega(T)$ and abatement cost function $\Lambda(\mu)$. The social planner optimizes abatement $\mu_t^*$ by Pontryagin's principle; the optimal carbon price equals the shadow price of the carbon state.

- The **SCC controversy** (Nordhaus $\approx$ $\$37$/tCO$_2$ vs. Stern $\approx$ $\$300$/tCO$_2$) turns almost entirely on the choice of pure time preference $\rho$ — a normative parameter with no empirically correct value. Higher damage function estimates (Burke et al.) and fat-tail catastrophic risk (Weitzman) both push the SCC above the Nordhaus baseline.

- **Physical risks** (extreme weather destroying capital and growth) and **transition risks** (stranded fossil fuel assets under rapid decarbonization) are macroeconomic risks now entering central bank stress-testing frameworks. The NGFS disorderly transition scenarios show substantially larger systemic losses than orderly scenarios.

- **Carbon tax vs. cap-and-trade**: under certainty, equivalent; under uncertainty, the choice between price and quantity instruments depends on whether the marginal damage or marginal abatement cost curve is steeper near the optimum. **Price corridors** in ETS designs provide bounded price certainty within a quantity commitment.

---

*Next: Chapter 38 — Inequality and Macroeconomics*
