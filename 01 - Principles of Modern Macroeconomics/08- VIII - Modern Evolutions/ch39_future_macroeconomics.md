# Chapter 39 — The Future of Macroeconomics: New Theories and Methods

> *"The ideas of economists and political philosophers, both when they are right and when they are wrong, are more powerful than is commonly understood."*
> — John Maynard Keynes, *The General Theory*, 1936

---

Macroeconomics is not a finished discipline. Each generation of economists has believed it had largely solved the core problems of the field, and each generation has been proved wrong in important ways. The classical economists believed Say's Law ruled out demand-driven recessions; the Great Depression refuted them. The New Keynesian synthesizers of the 2000s believed that sophisticated inflation-targeting central banks had largely solved the business cycle; the Global Financial Crisis refuted that. The appropriate posture toward the current frontier is not complacency but engaged uncertainty: significant progress has been made, and significant open questions remain. This chapter surveys the most active research frontiers — heterogeneous-agent models, machine learning, fat-tailed risk, and the macro-finance interface — providing the conceptual vocabulary to follow a rapidly evolving literature.

---

## 39.1 Heterogeneous-Agent Macroeconomics: HANK and Beyond

The most transformative development in macroeconomics since the Global Financial Crisis has been the shift from representative-agent to heterogeneous-agent models — recognizing that the distribution of income and wealth across households matters for aggregate dynamics in ways that the representative-agent shortcut suppresses.

### Why Heterogeneity Matters

In the representative-agent New Keynesian (RANK) model, aggregate consumption satisfies:

$$\hat{C}_t = \mathbb{E}_t[\hat{C}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\pi_{t+1}] - r^n).$$

Aggregate consumption is driven primarily by the real interest rate through intertemporal substitution — the Euler equation mechanism.

In the Heterogeneous-Agent New Keynesian (HANK) model, the aggregate consumption response to monetary policy has two components:

$$\Delta C_t^{HANK} = \underbrace{\text{Direct effect (intertemporal substitution)}}_{\text{small: $\approx 10\%$ of total}} + \underbrace{\text{Indirect effect (labor income channel)}}_{\text{large: $\approx 90\%$ of total}}.$$

Kaplan, Moll, and Violante (2018) find that the indirect effect — working through the employment and wage changes that monetary policy generates — accounts for approximately 90% of the aggregate consumption response. The reason: many households are "hand-to-mouth" (they consume their entire income each period) either because they hold no liquid assets (wealthy hand-to-mouth households with illiquid wealth) or because they have neither liquid nor illiquid wealth. For these households, the intertemporal substitution channel is inactive; what matters is whether they are employed and how much they earn.

This reversal has profound implications for policy:

- **Monetary transmission** works primarily through the labor market channel, not the Euler equation channel. Policies that affect employment and wages (including fiscal policies) are the dominant monetary transmission mechanism.
- **Fiscal transfer policy** is highly effective at the ELB: transfers to hand-to-mouth households generate MPC near 1, so helicopter money has large aggregate effects even without any employment response.
- **Distributional consequences of monetary policy** are large and asymmetric, as Chapter 38 develops.

### Computational Methods for HANK

The computational challenge is formidable. The state space includes not only aggregate variables (output, inflation, capital) but the entire cross-sectional distribution of household wealth and income $\mathcal{F}_t(a, y)$ — an infinite-dimensional object.

**Reiter's (2009) linearization method**: linearize the household decision rules and the distribution jointly around the stationary equilibrium, reducing the problem to a high-dimensional but linear system solvable with standard rational expectations methods. This approach can handle distributions with several hundred grid points efficiently.

**Winberry's (2018) perturbation method**: parameterize the distribution using polynomial coefficients; perturb around the steady state in the distribution moments. Dramatically reduces the state space while maintaining reasonable accuracy for the business cycle moments of interest.

**Achdou, Han, Lasry, Lions, and Moll (2022) continuous-time methods**: reformulate the household problem as a pair of partial differential equations — the Hamilton-Jacobi-Bellman (HJB) equation for the value function and the Kolmogorov-Fokker-Planck (KFP) equation for the distribution's evolution. These PDEs can be solved efficiently on fine grids using finite-difference methods, and many analytical insights unavailable in discrete-time models become available in continuous time.

The continuous-time HANK approach has spawned a growing literature deriving closed-form or near-closed-form results: expressions for the aggregate MPC, the redistribution effects of monetary policy, and the welfare costs of business cycles — insights that complement the numerical results from the discrete-time literature.

---

## 39.2 Machine Learning and Structural Macroeconomics

Machine learning methods have entered macroeconomics in two distinct roles: solving complex economic models and identifying causal relationships in data. Both represent methodological advances over the previous state of the art.

### Deep Learning for Model Solution

DSGE models with many state variables, occasionally binding constraints (ELB, collateral constraints), and non-linear dynamics are difficult to solve with perturbation or value function iteration. Neural networks — which can approximate arbitrary continuous functions to arbitrary accuracy on compact sets — offer a scalable alternative.

The approach (Duarte, 2018; Maliar and Maliar, 2021): parameterize the value function or policy function as a neural network and minimize the Bellman equation error across simulated state-space points:

$$\mathcal{L} = \mathbb{E}\!\left[\left\|V_\theta(s) - \max_c\{u(c) + \beta\mathbb{E}[V_\theta(s')]\}\right\|^2\right],$$

where $\theta$ are the network parameters updated by stochastic gradient descent. The neural network scales to state spaces of 20–100 dimensions that are intractable for conventional projection methods — enabling DSGE models with rich heterogeneity, multiple occasionally binding constraints, and complex financial sectors.

**Reinforcement learning for policy games**: in models where the central bank and fiscal authority optimize simultaneously, multi-agent reinforcement learning algorithms learn equilibrium policy functions through iterative simulation, providing a new approach to time-inconsistency and fiscal-monetary interaction models.

### Double Machine Learning for Causal Inference

The challenge in empirical macroeconomics is identifying causal effects with many potential confounders and limited degrees of freedom. **Double Machine Learning** (Chernozhukov et al., 2018) allows researchers to control for many nuisance variables without introducing regularization bias in the coefficient of interest:

1. Regress the outcome $Y$ on controls $\mathbf{X}$ using a flexible ML method (LASSO, random forest, neural network) to obtain residuals $\tilde{Y} = Y - \hat{m}_Y(\mathbf{X})$.
2. Regress the treatment $D$ on controls $\mathbf{X}$ to obtain residuals $\tilde{D} = D - \hat{m}_D(\mathbf{X})$.
3. Regress $\tilde{Y}$ on $\tilde{D}$ by OLS to obtain the causal estimate.

This "partialing out" is valid under mild conditions on the ML estimators and controls for high-dimensional confounders that traditional IV and panel methods cannot accommodate — critical for macroeconomic applications where hundreds of potential controls are available but degrees of freedom are scarce.

**Local projections with machine learning controls** (Plagborg-Møller and Wolf, 2021): local projections for IRF estimation can be augmented with ML-selected controls to address omitted variable bias while maintaining the non-parametric flexibility of local projections over restricted SVAR lag structures.

---

## 39.3 Non-Gaussian Risks: Fat Tails, Rare Disasters, and Ambiguity

Standard DSGE models assume Gaussian shocks. This assumption is computationally convenient but empirically problematic: macroeconomic outcomes have fat tails — extreme events occur far more frequently than Gaussian probability implies.

### Rare Disasters

Barro (2006) documents that consumption disasters of 15% or more have occurred in approximately 3.5% of country-years historically. The **rare disasters model** adds disaster states to the consumption growth process:

$$\Delta\ln c_{t+1} = \mu_c + \epsilon_{t+1}^c + v_{t+1}^c\cdot\mathbf{1}\{J_{t+1}=1\},$$

where $J_{t+1} \sim \text{Bernoulli}(p)$ with $p \approx 0.035$ and $v_{t+1}^c$ is disaster severity (mean approximately $-30\%$). In disaster states, the stochastic discount factor $M_{t+1} = \beta(c_{t+1}/c_t)^{-\sigma}$ becomes very large (marginal utility is very high when consumption collapses), generating large risk premia on assets that pay poorly in disasters. With $p = 0.035$ and $\bar{v} \approx -0.3$, the model matches the historical equity premium of approximately 6% with $\sigma \approx 4$ — far more plausible than the $\sigma \approx 50$ required by the Mehra-Prescott standard model.

The rare disasters framework has also been applied to currency crashes (Farhi and Gabaix, 2016) and sovereign debt crises (Roch and Uhlig, 2018), generating risk premia and return patterns consistent with observed data.

### Ambiguity Aversion and Robust Control

Standard models assume risk — uncertainty with known probabilities. Many decisions involve **Knightian uncertainty**: situations where even the probability distribution is unknown. **Ambiguity aversion** (Ellsberg, 1961) — preferring known risks to unknown ones — is empirically robust and macroeconomically significant.

Hansen and Sargent (2008) develop a **robust control** framework: rather than optimizing under a single assumed model, the decision-maker minimizes welfare loss under the worst-case model within an uncertainty set:

$$\pi^{robust} = \arg\min_\pi\max_{P \in \mathcal{U}} L(\pi, P).$$

The optimal robust policy is more cautious and more inertial than the point-estimate optimal policy — consistent with observed central bank gradualism. Interest rate smoothing, in this framework, is rational under model uncertainty: large rate changes would be very costly if the model is wrong; gradual changes provide the opportunity to observe the economy's response and adjust.

---

## 39.4 The Macro-Finance Interface

The integration of financial economics with macroeconomics — modeling macroeconomic dynamics with financial sector balance sheets as state variables — is the most active research area in macroeconomics at the frontier.

### Intermediary Asset Pricing

He and Krishnamurthy (2013) and Brunnermeier and Sannikov (2014) model financial intermediaries as the marginal investor in asset markets. The intermediary sector's net worth $W_t$ (as a fraction of total assets) is the key state variable. The intermediary SDF prices assets based on intermediary health rather than household consumption growth, generating predictions:

- Asset risk premia are highest when intermediary net worth is low (following losses) — consistent with the flight-to-quality observed after Lehman.
- Liquidity premia move with intermediary balance sheet conditions — consistent with the breakdown of covered interest parity post-2008.
- Recovery from financial crises is slow because intermediaries must rebuild equity gradually before credit conditions normalize — consistent with the slow post-2008 recovery.

The intermediary asset pricing framework provides a macroeconomic foundation for financial market anomalies that are inexplicable in representative-agent RANK models, and it connects the financial accelerator of Chapter 24 to asset pricing in a unified equilibrium framework.

### The Natural Rate of Interest and Secular Stagnation

The decline of the natural rate $r^n$ from approximately 3–4% in the 1980s to near zero or below in the 2010s is one of the most important empirical facts in macroeconomics (Laubach and Williams, 2003; Rachel and Smith, 2017). Proposed causes include demographic aging (higher aggregate saving, lower investment demand), rising inequality (wealthy households save more than poor ones — the Mian-Straub-Sufi mechanism), a global savings glut, and a fall in the relative price of investment goods.

The natural rate determines the ELB problem's severity: when $r^n$ is near zero, the Taylor rule prescribes nearly zero nominal rates at full employment, leaving the ELB almost always binding in recessions. Whether $r^n$ will remain low or revert toward historical levels is a central unknown for monetary policy frameworks. Arguments for reversal: climate transition investment needs, defense spending increases, and demographic reversals in saving behavior. Arguments for persistence: the structural forces that drove $r^n$ down (globalization, technology, inequality) remain in place. The empirical distinction requires years of additional data and represents one of the most consequential open questions in macroeconomics.

---

## Chapter Summary

- In **HANK models** (Kaplan-Moll-Violante), the indirect income channel accounts for approximately 90% of the aggregate consumption response to monetary policy — inverting the RANK result. Hand-to-mouth households with MPCs near 1 amplify fiscal transfers; distributional consequences of monetary policy are large. Computational advances (Reiter linearization, Winberry perturbation, Achdou-Moll continuous-time PDE methods) make large-scale HANK tractable.

- **Machine learning** enters macroeconomics as a model-solution tool (neural networks for high-dimensional DSGE) and as a causal inference tool (Double Machine Learning for partialing out high-dimensional confounders). Reinforcement learning enables equilibrium computation in policy games.

- **Rare disasters** (Barro, 2006) — with $p \approx 3.5\%$ probability and $\sim 30\%$ severity — generate large equity premia with plausible risk aversion ($\sigma \approx 4$) without the $\sigma \approx 50$ of the Mehra-Prescott puzzle. **Robust control** (Hansen-Sargent) formalizes ambiguity aversion, generating gradual, inertial optimal policies consistent with observed central bank behavior.

- **Intermediary asset pricing** (He-Krishnamurthy; Brunnermeier-Sannikov) makes the financial sector's net worth a state variable for asset prices, explaining post-crisis risk premia, CIP deviations, and slow recoveries. The **declining natural rate** $r^n$ (from $\sim 3.5\%$ to near zero since the 1980s) makes ELB binding more frequent and is one of the most consequential open empirical questions for monetary policy design.

---

*Next: Chapter 40 — Case Study: The Great Recession of 2008*
