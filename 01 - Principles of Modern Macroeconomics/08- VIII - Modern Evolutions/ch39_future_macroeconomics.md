# Chapter 39 — The Future of Macroeconomics: New Theories and Methods

---

Macroeconomics is not a finished discipline. Each generation of economists has believed it had largely solved the core problems of the field, and each generation has been proved wrong in important ways. The classical economists believed Say's Law ruled out recessions; the Great Depression refuted them. The neoclassical synthesizers believed fiscal and monetary policy could permanently stabilize unemployment below its natural rate; the stagflation of the 1970s refuted them. The New Keynesian consensus of the 2000s believed that sophisticated inflation-targeting central banks had largely solved the business cycle problem; the Global Financial Crisis refuted that. The appropriate posture toward the current frontier is not "we have largely solved it" but "we have made significant progress and face significant open questions."

This chapter surveys the most active research frontiers in macroeconomics — areas where current models are known to be inadequate, where significant methodological innovation is underway, and where the answers are genuinely uncertain. The goal is not to provide a definitive summary of an unsettled literature, but to give the reader the conceptual vocabulary to follow it as it evolves.

---

## 39.1 Heterogeneous-Agent Macroeconomics: HANK and Beyond

The most transformative development in macroeconomics since the Global Financial Crisis has been the shift from representative-agent to heterogeneous-agent models. Chapter 25 introduced the HANK (Heterogeneous-Agent New Keynesian) framework; this section examines its broader significance and the computational and theoretical challenges of extending it.

### Why Heterogeneity Matters for Aggregates

In the representative-agent NK (RANK) model, the aggregate consumption Euler equation is:

$$\hat{C}_t = \mathbb{E}_t[\hat{C}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\pi_{t+1}] - r^n).$$

This equation says that aggregate consumption is driven primarily by the real interest rate. The intertemporal substitution mechanism is the dominant channel.

In the HANK model, the aggregate consumption response to a monetary policy shock is:

$$\Delta C_t^{HANK} = \underbrace{\text{Direct effect (intertemporal substitution)}}_{\text{small, as in RANK}} + \underbrace{\text{Indirect effect (income channel: employment, wages)}}_{\text{large, via hand-to-mouth households}}.$$

Kaplan, Moll, and Violante (2018) find quantitatively that the indirect (income) channel accounts for approximately 90% of the aggregate consumption response to a monetary policy shock in their HANK model — the direct intertemporal substitution channel, dominant in RANK, is almost negligible. This reversal has profound implications: it suggests that monetary policy works not primarily through the Euler equation mechanism that the RANK literature emphasized, but through the labor market and employment channel.

The HANK implication extends to fiscal policy: because liquidity-constrained households have MPCs near one, helicopter money transfers (direct payments to households, regardless of the labor market effect) generate large multipliers even without any employment effect. This makes transfer-based fiscal policy substantially more powerful than the RANK model suggests.

### Computational Methods for Solving HANK Models

The computational challenge of HANK is formidable: the state space includes not only aggregate state variables (output, inflation, the capital stock) but also the entire distribution of household wealth and income $\mathcal{F}_t(a, y)$, which is an infinite-dimensional object. Several approaches have made large-scale HANK models tractable:

**Reiter's (2009) method**: linearize the household's decision rules and the distribution around the stationary equilibrium, reducing the problem to a high-dimensional but linear system that can be solved with standard rational expectations techniques.

**Winberry's (2018) perturbation method**: use polynomial approximations to the distribution, dramatically reducing the state space while maintaining reasonable accuracy.

**Achdou, Han, Lasry, Lions, and Moll (2022) continuous-time methods**: reformulate the household problem as a system of partial differential equations (the Hamilton–Jacobi–Bellman equation and the Kolmogorov–Fokker–Planck equation for the distribution evolution) that can be solved efficiently with finite-difference methods on fine grids.

The continuous-time approach has proven particularly powerful: many analytical results that are unavailable in discrete-time HANK models can be derived in continuous time, providing theoretical insights alongside quantitative results.

---

## 39.2 Machine Learning and Structural Macroeconomics

Machine learning (ML) methods have begun entering macroeconomics in two distinct roles: as tools for solving complex economic models and as tools for identifying causal relationships in data.

### Deep Learning for DSGE Model Solution

DSGE models with many state variables and non-linear dynamics (occasionally binding constraints, occasionally binding ELB, financial frictions) are difficult to solve with conventional perturbation or value function iteration methods. Neural networks — which can approximate arbitrary functions — have shown promise as universal approximators for value functions and policy functions.

Duarte (2018) proposes training neural networks to represent the value function $V(k, A, z)$ in a model with aggregate capital $k$, productivity $A$, and a financial friction state $z$, using the Bellman equation as the loss function:

$$\mathcal{L}_{Bellman} = \left\|V(k, A, z) - \max_c\{u(c) + \beta\mathbb{E}[V(k', A', z')]\}\right\|^2.$$

The neural network is trained by stochastic gradient descent to minimize the Bellman error. Compared to conventional projection methods, the neural approach scales much better to high-dimensional state spaces, though at the cost of reduced accuracy guarantees and interpretability.

**Reinforcement learning** for economic policy: policy games in which the central bank and the fiscal authority optimize simultaneously (a dynamic game rather than a single-agent optimization) can be solved using multi-agent reinforcement learning algorithms that learn equilibrium policy functions through iterative simulation.

### Machine Learning for Causal Inference

The challenge in empirical macroeconomics is identifying causal effects in the presence of many potential confounders and limited degrees of freedom. Two ML approaches have proven particularly useful.

**Double Machine Learning** (Chernozhukov, Chetverikov, Demirer, Duflo, Hansen, Newey, and Robins, 2018): allows researchers to include many nuisance controls (variables that affect the outcome and treatment but are not the object of interest) without introducing regularization bias in the coefficient of interest. The procedure:

1. Regress the outcome $Y$ on controls $X$ using a flexible ML method (lasso, random forest, neural network) to obtain residuals $\tilde{Y}$.
2. Regress the treatment $D$ on controls $X$ using the same ML method to obtain residuals $\tilde{D}$.
3. Regress $\tilde{Y}$ on $\tilde{D}$ by OLS to obtain the causal estimate.

This "partialing out" procedure is valid under mild conditions on the nuisance estimators and controls for high-dimensional confounders that traditional IV and panel methods cannot accommodate.

**Synthetic control** (Abadie, Diamond, and Hainmueller, 2010): constructs a weighted average of comparison units that matches the treated unit's pre-treatment trajectory, providing a credible counterfactual for policy evaluation. Extended via machine learning to handle many potential controls and high-dimensional matching.

---

## 39.3 Non-Gaussian Risks: Fat Tails, Rare Disasters, and Ambiguity

Standard DSGE models assume Gaussian shocks with thin tails. This assumption is convenient (linear approximations are accurate under small Gaussian shocks) but empirically problematic: the distribution of macroeconomic outcomes has fat tails — extreme events (depressions, hyperinflations, wars) occur far more often than Gaussian probability implies.

### Rare Disasters

Barro (2006) documents that consumption disasters of 15% or more have occurred in approximately 3.5% of country-years historically — far more frequently than the tails of a normal distribution would predict. The **rare disasters model** adds disaster states to the standard consumption growth process:

$$\Delta\ln c_{t+1} = \mu_c + \epsilon_{t+1}^c + v_{t+1}^c\cdot\mathbf{1}\{J_{t+1}=1\},$$

where $J_{t+1} \sim \text{Bernoulli}(p)$ with $p \approx 0.035$ and $v_{t+1}^c \sim F_J$ (disaster severity, with mean approximately $-30\%$). The SDF in a rare disasters model:

$$M_{t+1} = \beta\left(\frac{c_{t+1}}{c_t}\right)^{-\sigma}\cdot\mathbf{1}\{J_{t+1}=0\} + \beta\left(\frac{c_{t+1}}{c_t}\right)^{-\sigma}e^{-\sigma v_{t+1}^c}\cdot\mathbf{1}\{J_{t+1}=1\}.$$

In disaster states, the SDF is very large (consumption collapsed, marginal utility is very high), so assets that pay poorly in disasters carry large risk premia. With $p = 0.035$ and $v \sim -0.3$, the model can match the 6% equity premium with $\sigma \approx 4$ — far more plausible than the $\sigma \approx 50$ required by the Mehra–Prescott standard model.

### Ambiguity Aversion

Standard models assume risk — uncertainty with known probabilities — but many economic decisions involve **ambiguity** (Knightian uncertainty): situations where even the probability distribution is unknown. **Ambiguity aversion** (Ellsberg, 1961) is the empirically robust preference for known risks over unknown ones.

Hansen and Sargent (2008) develop a **robust control** framework that gives decision-makers with ambiguity aversion policies that perform well across a range of models rather than optimizing against a single assumed model. The central bank problem becomes:

$$\min_\pi\max_{P \in \mathcal{U}} L(\pi, P),$$

where $\mathcal{U}$ is an uncertainty set of probability models. The resulting "worst-case" optimal policy is more cautious and more inertial than the point-estimate optimal policy — consistent with the observed gradualism of central bank rate changes. This provides a theoretical foundation for interest rate smoothing as a rational response to model uncertainty rather than as an ad hoc assumption.

---

## 39.4 The Macro-Finance Interface

The integration of financial economics with macroeconomics — pricing assets with SDFs derived from macroeconomic models, and modeling macroeconomic dynamics with financial sector balance sheets as state variables — is the most active intersection of research at the time of writing.

### Intermediary Asset Pricing

He and Krishnamurthy (2013) and Brunnermeier and Sannikov (2014) model financial intermediaries as the marginal investor in asset markets. The intermediary sector's equity $W_t$ (net worth as a fraction of total assets) is the key state variable. The intermediary SDF:

$$M_{t+1}^{int} = \beta\frac{W_{t+1}}{W_t}\cdot\text{(risk adjustment)},$$

prices assets based on the health of the intermediary sector rather than on consumption growth. This generates several novel predictions:

- Asset risk premia are highest when intermediary net worth is low (following losses) — consistent with the flight-to-quality observed after the Lehman failure.
- Liquidity premia move with intermediary balance sheet conditions — consistent with the breakdown of covered interest parity in 2008 and after 2015.
- Recovery from financial crises is slow because the intermediary sector must rebuild equity gradually — consistent with the slow post-2008 recovery.

### The Natural Rate of Interest and Secular Stagnation

The decline of the natural rate of interest $r^n$ from approximately 3–4% in the 1980s to near zero or below in the 2010s–2020s is one of the most important empirical facts in macroeconomics (Rachel and Smith, 2017; Laubach and Williams, 2003). Proposed causes include: demographic aging (higher saving, lower investment demand), rising inequality (wealthy households save more), a global savings glut, and a fall in the relative price of investment goods (computers, machinery are cheaper, so a given level of investment spending buys more). The natural rate determines the ELB problem: when $r^n$ is very low, the Taylor rule prescribes negative interest rates even in equilibrium, trapping the economy at the ELB regularly.

Whether $r^n$ will remain low or revert toward historical levels is the central unknown for monetary policy over the next decade. If it rises — due to climate investment needs, defense spending, and aging-reversal of saving — central banks will regain conventional policy room. If it remains depressed, AIT, QE, and forward guidance will be permanent rather than emergency features of the monetary policy toolkit.

---

*Next: Chapter 40 — Case Study: The Great Recession of 2008*
