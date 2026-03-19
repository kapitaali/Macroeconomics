# Chapter 37 — Climate Change and Macroeconomics: Sustainability and Green Growth

> *"Climate change is the greatest market failure the world has ever seen."*
> — Nicholas Stern, *The Stern Review*, 2006

---

Climate change presents macroeconomics with its most difficult and consequential challenge: a global externality whose costs are concentrated in the distant future but whose causes are embedded in the current structure of production and consumption, and whose mitigation requires international coordination on an unprecedented scale. The macroeconomic dimensions of climate change are vast: they encompass the long-run growth implications of rising temperatures and extreme weather events, the optimal carbon pricing that internalizes climate damages, the macroeconomic transition risks from decarbonization, and the distributional consequences of both climate impacts and mitigation policies. This chapter develops the formal tools for thinking about these questions.

---

## 37.1 Climate Change as an Externality

**Definition (Externality).** An **externality** is a cost or benefit that falls on parties not involved in a transaction. A **negative externality** is a cost imposed on third parties; a **positive externality** is a benefit conferred. The burning of fossil fuels imposes a **negative externality** on current and future generations through greenhouse gas emissions: the emitter does not pay the full social cost of her emissions, leading to overproduction of carbon-intensive goods and underinvestment in abatement.

The competitive market equilibrium without a carbon price produces a quantity of carbon emissions $E^{market}$ such that the marginal private cost of emission equals the marginal private benefit:

$$MPC(E^{market}) = MPB(E^{market}).$$

The socially optimal quantity $E^*$ satisfies:

$$MPC(E^*) + MEC(E^*) = MPB(E^*),$$

where $MEC(E^*)$ is the **marginal external cost** (the marginal climate damage imposed on all current and future people per additional ton of CO$_2$ equivalent). Since $MEC > 0$, the socially optimal emission level $E^* < E^{market}$: the market overproduces carbon emissions relative to the social optimum.

A **Pigouvian carbon tax** at rate $\tau^* = MEC(E^*)$ — equal to the marginal external cost at the optimal emission level — corrects this externality by aligning the private and social cost of emissions. In principle, this single tax is sufficient to achieve the socially optimal allocation without any quantitative restrictions, subsidies, or technology mandates. The simplicity and efficiency of the carbon tax — it is the "first-best" policy — is the theoretical benchmark against which all real-world climate policies are evaluated.

---

## 37.2 Integrated Assessment Models: The DICE Framework

The formal quantitative framework for analyzing climate economics is the **integrated assessment model (IAM)**, which links an economic model to a physical climate model to assess the welfare costs and benefits of alternative emission reduction paths. The most influential IAM is Nordhaus's DICE (Dynamic Integrated model of Climate and the Economy), developed from 1992 onward.

### The Structure of DICE

The economic module is a modified Ramsey growth model. Output net of climate damage and abatement costs:

$$Y_t^{net} = A_t K_t^\alpha L_t^{1-\alpha}\cdot\Omega(T_t)\cdot(1 - \Lambda(\mu_t)),$$

where:
- $A_t K_t^\alpha L_t^{1-\alpha}$ is gross output from the standard Cobb–Douglas production function.
- $\Omega(T_t) = \frac{1}{1 + \pi_1 T_t + \pi_2 T_t^2}$ is the **damage function**, which reduces output as global mean temperature anomaly $T_t$ rises. With $\pi_1 \approx 0$, $\pi_2 \approx 0.00267$ (Nordhaus's estimate), the damage from a 3°C warming is approximately $\Omega(3) = 1/(1 + 0.0267\times 9) \approx 0.91$ — a 9% output loss.
- $\Lambda(\mu_t) = b_1\mu_t^{b_2}$ is the **abatement cost function**, with $b_1$ the cost coefficient and $b_2 \approx 2.6$ the exponent. $\mu_t \in [0,1]$ is the emission reduction rate (fraction of uncontrolled emissions that are abated). Abatement is costly, with costs rising superlinearly in the abatement rate.

The climate module links carbon emissions to atmospheric CO$_2$ concentration to global mean temperature:

$$E_t = (1-\mu_t)A_t\sigma_t K_t^\alpha L_t^{1-\alpha}, \quad M_{t+1} = M_t + \beta_E E_t - \beta_M(M_t - M^{pre}),$$
$$T_{t+1} = T_t + \zeta_T[\eta\ln(M_{t+1}/M^{pre}) - T_t],$$

where $\sigma_t$ is the carbon intensity of output (decreasing over time as the energy mix decarbonizes exogenously), $M_t$ is atmospheric CO$_2$ concentration, and $T_t$ is global mean temperature anomaly. The climate sensitivity parameter $\eta$ determines how much warming results from a doubling of CO$_2$ concentration; the central estimate from climate science is $\eta \approx 3$°C (range 2.5–4°C).

The social planner maximizes:

$$\max_{\{\mu_t, K_{t+1}\}}\; U = \int_0^\infty e^{-\rho t}\frac{c_t^{1-\eta}-1}{1-\eta}L_t\,\mathrm{d}t,$$

subject to the economic and climate modules, resource constraints, and transversality conditions.

---

## 37.3 The Social Cost of Carbon

**Definition (Social Cost of Carbon).** The **social cost of carbon (SCC)** is the present discounted value of the marginal damage caused by an additional ton of CO$_2$ emissions today, measured in current dollars. It is the correct Pigouvian carbon tax that internalizes the full external cost of emissions:

$$SCC_t = -\mathbb{E}_t\int_t^\infty e^{-\int_t^s r(\tau)\,\mathrm{d}\tau}\frac{\partial Y_s^{net}}{\partial E_t}\,\mathrm{d}s,$$

where the integral accumulates all future marginal damages from the additional ton of emissions, discounted at the social discount rate $r(\tau)$.

The SCC is exquisitely sensitive to the discount rate. Two competing positions have generated one of the most important normative debates in economics:

**The Nordhaus approach**: use market interest rates to discount future damages — $\rho \approx 1.5\%$ pure rate of time preference, $\sigma = 1.45$ risk aversion, real growth $g \approx 2\%$, giving a social discount rate $r = \rho + \sigma g \approx 4.4\%$. This reflects the view that markets implicitly reveal social preferences for trading off present against future consumption. At this discount rate, the Nordhaus optimal carbon price (SCC) is approximately $\$37$/tCO$_2$ in 2010 dollars, rising to approximately $\$100$/tCO$_2$ by 2050.

**The Stern approach**: use a near-zero pure rate of time preference ($\rho \approx 0.1\%$) on ethical grounds that future people's welfare should not be discounted merely because they are in the future. With $\sigma = 1$ and $g \approx 1.3\%$, the Stern discount rate is $r \approx 1.4\%$ — far lower than Nordhaus's. At this rate, the SCC is approximately $\$300$/tCO$_2$ — roughly eight times the Nordhaus estimate. This generates dramatically stronger current mitigation policy.

The core of the disagreement is normative, not empirical. The pure rate of time preference $\rho$ is a value judgment about the moral equivalence of present and future generations. There is no "correct" answer from economics alone; the choice requires ethical reasoning. However, economists can clarify the implications: the Nordhaus approach implies accepting approximately 3–4°C of warming with modest mitigation costs; the Stern approach implies aggressive decarbonization to limit warming to 1.5–2°C at substantially higher current costs.

### Uncertainty and Fat Tails

A further reason to favor aggressive early action is the **fat-tailed uncertainty** about climate damages. The damage function $\Omega(T)$ is estimated from economic and physical science studies with large uncertainty. Weitzman (2009) argues that the standard expected utility framework understates the cost of extreme climate outcomes because the distribution of possible climate damages has very heavy tails: there is a small but non-negligible probability of catastrophic warming (5–6°C) that would destroy a large fraction of world output. When damages are catastrophic and tails are heavy, risk aversion alone (even without discounting) generates very high optimal carbon prices.

---

## 37.4 Macroeconomic Risks from Climate Change

Beyond the smooth damages captured by DICE's quadratic damage function, climate change generates macroeconomic risks that are more abrupt and harder to model.

**Physical risks**: increases in the frequency and severity of extreme weather events (hurricanes, floods, wildfires, heat waves) destroy physical capital, disrupt supply chains, reduce agricultural yields, and generate large migration flows. The macroeconomic modeling of these risks is in its early stages; Hsiang and Jina (2014) estimate that an average tropical cyclone reduces GDP growth for 20 years after the event, with the long-run GDP loss exceeding 50 times the immediate physical damage.

**Transition risks**: a rapid shift to low-carbon energy systems will strand fossil fuel assets — coal mines, oil and gas fields, refineries, pipelines, and power plants that become uneconomic before the end of their planned lives under a stringent carbon price. The IPCC (2018) estimates stranded fossil fuel assets at $1–4 trillion globally under a 1.5°C warming scenario. These stranded assets could generate severe losses for the financial institutions that hold them, potentially triggering a financial crisis if the transition is sudden rather than gradual.

Central banks and financial regulators (the Network for Greening the Financial System, NGFS) have begun conducting **climate stress tests** for banks: scenarios in which a sudden carbon price increase strands fossil fuel assets, generating large credit losses and potential systemic risk. The results suggest significant but manageable financial system exposure, with the caveat that the models used are linear and may understate tail risks from abrupt transitions.

---

## 37.5 Carbon Pricing Instruments: Tax versus Cap-and-Trade

Two primary instruments implement the carbon price:

**Carbon tax**: a per-unit tax on the carbon content of fuels, equal to the SCC. The government sets the price; the market determines the quantity of emissions reduction. A carbon tax provides certainty about the price of carbon, enabling long-term investment planning by firms and households.

**Cap-and-trade (emissions trading system, ETS)**: the government sets a cap on total emissions and distributes (by auction or free allocation) permits equal to the cap. Firms must surrender one permit per ton of CO$_2$ emitted; they can buy and sell permits in a market. The market price of permits equals the SCC under efficient markets. A cap-and-trade provides certainty about the quantity of emissions but uncertainty about the carbon price.

From a welfare perspective, under certainty and with no market failures, the two instruments are equivalent (Weitzman, 1974). Under uncertainty, the relative performance depends on whether the damage function is convex or concave:

- If marginal climate damages are steep (convex damage function), it is more important to limit the quantity of emissions; cap-and-trade is preferred.
- If marginal abatement costs are steep (steep MAC curve), it is more important to limit the price of carbon; a carbon tax is preferred.

In practice, cap-and-trade systems (EU ETS, California's system, RGGI) allow the price to fluctuate with economic conditions and political changes, creating uncertainty that may deter long-term low-carbon investment. Carbon taxes (British Columbia's, Sweden's) provide a stable price signal but face stronger political opposition when fuel price spikes make the tax salient. A **carbon price corridor** — a floor and ceiling price in an ETS — combines elements of both approaches, providing bounded price certainty while maintaining a quantity commitment.

---

*Next: Chapter 38 — Inequality and Macroeconomics*
