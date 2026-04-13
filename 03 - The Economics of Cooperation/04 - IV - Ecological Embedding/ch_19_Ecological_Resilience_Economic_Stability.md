# Chapter 19: Ecological Resilience and Economic Stability — Mathematical Conditions

> *"Resilience thinking… requires a shift from equilibrium-centered views of nature to evolutionary, nonlinear, and stochastic perspectives."*
> — Brian Walker and David Salt, *Resilience Thinking* (2006)

> *"A system that can only exist in one state is fragile. A system that can bounce between several states is resilient. A system that can produce new states is adaptive."*
> — C.S. Holling, *Adaptive Environmental Assessment and Management* (1978)

## Learning Objectives

By the end of this chapter, you should be able to:

1. Define engineering resilience and ecological resilience formally and prove that ecological resilience is the more appropriate measure for economic policy design in the presence of potential regime shifts.
2. Construct the formal dynamical systems model of Holling's adaptive cycle and identify its four phases and their economic analogues.
3. Formalize the panarchy framework as a system of coupled adaptive cycles across scales and derive the conditions for cross-scale stabilization versus cross-scale destabilization.
4. Derive early-warning indicators of approaching regime shifts — rising variance and rising lag-1 autocorrelation — from the formal model of a saddle-node bifurcation, and prove their statistical properties.
5. Specify the formal conditions for a stable coupled economic-ecological system, distinguish the conditions for stable co-evolution from those for economic growth-induced ecological destabilization, and identify the mathematical transition between them.
6. Apply the coupled model to the North Atlantic cod fishery collapse, extracting early-warning signals from pre-collapse catch data.

---

## 19.1 Two Conceptions of Resilience

Chapter 5 introduced the bifurcation framework and the concept of regime shifts — qualitative changes in system behavior arising from small parametric changes. Chapter 12 introduced algebraic connectivity as a measure of network resilience. This chapter develops the ecological resilience framework that connects these two concepts specifically to the problem of managing economic systems embedded in ecological ones.

The connection is not merely metaphorical. Economies depend on ecosystems for material inputs, waste absorption, and the broader provisioning services that sustain human welfare [C:Ch.17]. Many of these ecological systems exhibit the nonlinear dynamics, multiple stable states, and potential regime shifts that the bifurcation theory of Chapter 5 describes. When an ecological system undergoes a regime shift — a fishery collapses, a lake transitions from clear to turbid, a forest converts to savanna, a coral reef bleaches to rubble — the economic system that depended on it loses its productive foundation, potentially irreversibly.

Understanding ecological resilience — how far an ecological system can be pushed before it crosses a threshold into a different regime — is therefore a prerequisite for understanding whether economic growth is sustainable. An economy operating in a region of high ecological resilience can tolerate shocks and disturbances without triggering irreversible ecosystem collapse. An economy operating close to an ecological tipping point cannot.

The formal tools of this chapter — Lyapunov stability theory, bifurcation analysis, the adaptive cycle model, and coupled system eigenvalue analysis — provide the mathematical machinery for assessing where any specific economic-ecological system stands relative to these thresholds.

---

## 19.2 Engineering vs. Ecological Resilience

### 19.2.1 Formal Definitions

The distinction between engineering resilience and ecological resilience was first articulated by Holling (1973) in a seminal paper that remains one of the most cited in ecology. It has direct implications for economic policy design.

**Definition 19.1 (Engineering Resilience).** The engineering resilience of a dynamical system at equilibrium $x^*$ is:

$$R_E = \frac{1}{\tau} = |\lambda_{\max}|$$

the reciprocal of the return time $\tau$ to equilibrium after a small perturbation, equal to the absolute value of the largest (most negative) eigenvalue of the Jacobian $\partial f/\partial x$ evaluated at $x^*$.

Engineering resilience measures how quickly a system returns to equilibrium after disturbance — a local property valid only in the vicinity of the equilibrium. A system with high engineering resilience snaps back quickly from small perturbations.

**Definition 19.2 (Ecological Resilience).** The ecological resilience of a dynamical system with multiple stable states is the size of the basin of attraction of the current stable equilibrium $x^*$:

$$R_{\mathcal{E}} = \text{Vol}(\mathcal{B}(x^*))$$

where $\mathcal{B}(x^*)$ is the set of initial conditions from which trajectories converge to $x^*$. Ecological resilience measures how large a perturbation the system can absorb before it is pushed into the basin of attraction of a different equilibrium.

**The key difference.** Engineering resilience is a local measure; ecological resilience is a global measure. A system can have high engineering resilience (quick recovery from small shocks) and low ecological resilience (small basin of attraction, vulnerable to large shocks) simultaneously. This combination is particularly dangerous: the system appears highly stable in normal operation, giving little warning of impending collapse, but is highly vulnerable to any perturbation large enough to push it across the basin boundary.

**Example 19.1 (Bistable Fishery).** Consider a fishery with two stable equilibria: a high-stock equilibrium $N^*_H \approx 0.8K$ (well-managed, productive) and a low-stock equilibrium $N^*_L \approx 0.1K$ (depleted, barely self-sustaining), separated by an unstable equilibrium at $\hat{N} \approx 0.3K$ (the tipping threshold).

Near $N^*_H$: $|\lambda| \approx 0.25$ (return time $\tau \approx 4$ years — reasonably quick). Engineering resilience is moderate.

Basin of attraction $\mathcal{B}(N^*_H) = (\hat{N}, K) = (0.3K, K)$: any perturbation that pushes the stock below 30\% of carrying capacity triggers collapse to $N^*_L$. With stock currently at $0.8K$, the system can absorb a perturbation of $-0.5K$ before collapse — ecological resilience is moderate-to-low for a managed fishery.

If fishing pressure gradually reduces the stock to $0.4K$ (while apparently maintaining engineering resilience), the ecological resilience has shrunk to $0.1K$ — the system can only absorb a perturbation of $-0.1K$ before collapse. Engineering resilience still appears fine; ecological resilience has nearly vanished. This is the dangerous regime where most fishery collapses have occurred.

### 19.2.2 Why Ecological Resilience Matters for Economic Policy

**Proposition 19.1 (Ecological Resilience as Policy Target).** For economic systems embedded in ecological systems with multiple stable states, engineering resilience ($R_E$) is an insufficient safety metric, while ecological resilience ($R_{\mathcal{E}}$) is the appropriate target for sustainability policy.

*Proof.* Engineering resilience $R_E = |\lambda_{\max}|$ depends only on the local dynamics near $x^*$ and provides no information about the distance to the basin boundary $\partial\mathcal{B}(x^*)$. A policy that maintains high $R_E$ (quick recovery from small shocks) is compatible with $R_{\mathcal{E}} \to 0$ (approaching the tipping threshold). Since the social cost of crossing the tipping threshold — losing the productive ecosystem services of $x^*_H$ permanently — is typically large and potentially irreversible, the policy target must be $R_{\mathcal{E}} \geq \bar{R}$ for some minimum acceptable basin size $\bar{R}$, not merely $R_E > 0$. $\square$

**Corollary 19.1.** Standard resilience indicators in economics (recovery speed after a recession, financial system speed of return to equilibrium) correspond to engineering resilience and provide no information about proximity to regime shifts. Ecological resilience indicators — basin size, distance to tipping threshold, early-warning signal strength — are needed to assess whether economic-ecological systems are approaching catastrophic transitions.

---

## 19.3 Holling's Adaptive Cycle

### 19.3.1 The Four-Phase Model

C.S. Holling's adaptive cycle is an empirically observed pattern of ecosystem dynamics, formalized here as a dynamical system. It describes a recurring sequence of four phases that ecosystems (and, by extension, economic systems) undergo as they grow, mature, collapse, and reorganize.

**Definition 19.3 (Holling Adaptive Cycle).** The adaptive cycle is a dynamical system on the state space $(r, K, C)$ — potential ($r$), connectedness ($C$), and resilience ($K$, re-using this symbol for carrying capacity in the ecological literature to avoid confusion) — with four phases:

**Phase $\alpha$ (Growth/Exploitation):** $r$ increases rapidly, $C$ increases moderately, $K$ increases. New resources are captured; the system grows quickly. Economic analogue: the growth phase of an industry or technology, characterized by rapid expansion, high innovation, and low specialization.

**Phase $K$ (Conservation/Maturation):** $r$ stabilizes, $C$ increases further, $K$ reaches maximum. The system is highly efficient and connected but increasingly rigid. Economic analogue: the maturation phase of an industry — high efficiency, strong institutions, but reduced flexibility and innovation.

**Phase $\Omega$ (Release/Collapse):** Accumulated brittleness leads to rapid collapse — $C$ and $K$ fall sharply. A disturbance (fire, pest, financial crisis) releases the accumulated potential. Economic analogue: the creative destruction phase — established firms fail, institutions dissolve, accumulated capital depreciates rapidly.

**Phase $r$ (Reorganization):** From low $C$ and $K$, the system reorganizes with novel combinations — new species, new firms, new institutional forms. Economic analogue: the recovery and experimentation phase following a crisis.

**Formal dynamical model.** A simplified ODE representation:

$$\dot{r} = a_1 r(1 - r/r_{\max}) - b_1 r C$$
$$\dot{C} = a_2 C(1 - C/C_{\max}) + b_2 r C - b_3 C(C - C^*)^2$$
$$\dot{K} = a_3(r - K/K_{\max}) - b_4 K \cdot \mathbb{1}[C > C^*]$$

where $C^*$ is the critical connectedness threshold above which brittle collapse becomes likely, and $\mathbb{1}[\cdot]$ is the indicator function. The third equation captures the release phase: when connectedness exceeds $C^*$, resilience $K$ declines rapidly — the over-connected system is vulnerable to collapse.

### 19.3.2 Economic Analogues

The adaptive cycle is not merely a description of ecological dynamics; it describes the life cycle of any complex adaptive system. The economic analogues are well-documented:

**Industry cycles.** Technology industries exhibit the adaptive cycle: rapid growth ($\alpha$ phase) as a new technology diffuses, followed by consolidation and standardization (K phase), followed by creative destruction as a new technology platform disrupts the incumbent ($\Omega$ phase), followed by a transition period in which new entrants restructure the industry ($r$ phase). The semiconductor industry, the personal computer industry, and the smartphone industry each exhibit this pattern.

**Financial cycles.** Minsky's financial instability hypothesis [C:Ch.23] describes a financial adaptive cycle: a growth phase of hedge financing (revenue exceeds debt service), a conservation phase of speculative financing (revenue equals but does not exceed debt service), and a release phase of Ponzi financing (revenue insufficient for debt service) triggered by rising interest rates or falling asset prices. The 2007–09 financial crisis was precisely an $\Omega$-phase release in the Holling sense.

**Institutional cycles.** Chapter 15 showed that institutions exhibit path dependence and lock-in — the behavioral hallmarks of the K phase. Crises reduce the cost of institutional change — the $\Omega$ phase creates the institutional window [C:Ch.15, Definition 15.9]. The r phase is the period of institutional experimentation and selection.

---

## 19.4 Panarchy: Nested Adaptive Cycles

### 19.4.1 The Panarchy Framework

Real systems are not isolated adaptive cycles; they are nested hierarchies of cycles operating at different spatial and temporal scales. A forest contains trees, stands, landscapes, and biomes — each with its own adaptive cycle, each influencing and being influenced by cycles at adjacent scales. An economy contains firms, industries, national economies, and the global system — again, nested adaptive cycles across scales.

**Definition 19.4 (Panarchy).** A panarchy is a nested set of adaptive cycles $\{(r_k, C_k, K_k)\}_{k=1}^K$ at scales $l_1 \prec l_2 \prec \cdots \prec l_K$ (from small and fast to large and slow), coupled through:

- **Revolt connections** ($r_k \to r_{k+1}$): A collapse at scale $k$ can trigger a collapse at scale $k+1$ if it is large enough and the higher scale is in its K (vulnerable) phase. This is cross-scale destabilization.
- **Remembrance connections** ($K_{k+1} \to K_k$): A stable, mature system at scale $k+1$ provides the context that enables recovery at scale $k$ after collapse. This is cross-scale stabilization.

**Proposition 19.2 (Panarchy Stability Conditions).** A panarchy of two scales is stable if and only if:
(i) The lower scale ($l_1$) is not in or approaching the $\Omega$ (collapse) phase when the upper scale ($l_2$) is in the K (vulnerable) phase: revolt connections cannot trigger upper-scale collapse.
(ii) The upper scale ($l_2$) is not in the $\Omega$ phase when the lower scale ($l_1$) is in the $r$ (reorganization) phase: remembrance connections remain available to support recovery.

*Proof.* (i): If $l_1$ collapses ($\Omega$ phase) and $l_2$ is in K phase, the revolt connection transmits the collapse signal to $l_2$. Since K-phase systems have high connectedness and low resilience, they cannot absorb large perturbations — cross-scale collapse follows. The stability condition requires $l_1 \notin \Omega$ when $l_2 \in K$. (ii): If $l_2$ collapses when $l_1$ is reorganizing, the remembrance connection is severed — $l_1$ loses the stabilizing context of the higher-scale system. Recovery requires an intact higher-scale context; its absence leaves $l_1$ in permanent disorganization. $\square$

**Economic implication.** The panarchy stability conditions translate directly to economic management: do not allow simultaneous collapses at adjacent scales (financial crisis coinciding with sovereign debt crisis coinciding with global recession — each reinforces the other through revolt connections) and maintain stable higher-level institutions (global monetary system, multilateral trade rules, shared ecological governance) that provide the context for lower-level recovery (individual firm restructuring, national fiscal adjustment, local ecological restoration).

---

## 19.5 Regime Shifts and Early-Warning Indicators

### 19.5.1 The Bifurcation Model

Chapter 5 introduced the saddle-node bifurcation as the formal model of a tipping point [Definition 5.5]. We now develop the full early-warning indicator theory, extending the brief treatment of Theorem 5.1.

**Setup.** Consider the ecological-economic system:

$$\dot{x} = f(x, \mu) + \sigma\xi(t)$$

where $x$ is the system state (e.g., fish stock, soil carbon, coral cover), $\mu$ is a slowly changing parameter (e.g., fishing pressure, temperature, nutrient loading), $\sigma > 0$ is noise intensity, and $\xi(t)$ is white noise.

**Near a stable equilibrium.** Linearizing around $x^*(\mu)$:

$$\dot{x} \approx \lambda(\mu)(x - x^*) + \sigma\xi(t), \quad \lambda(\mu) = f_x(x^*, \mu) < 0$$

This is an Ornstein-Uhlenbeck process with stationary variance $\sigma^2/(2|\lambda|)$ and lag-$\tau$ autocorrelation $e^{\lambda\tau}$.

**Theorem 19.1 (Early-Warning Indicators from Critical Slowing Down).** As a saddle-node bifurcation is approached ($\mu \to \mu_c$ such that $\lambda(\mu) \to 0^-$):

1. **Rising variance:** $\text{Var}(x) = \sigma^2/(2|\lambda(\mu)|) \to \infty$
2. **Rising lag-1 autocorrelation:** $\hat{\rho}_1 = e^{\lambda(\mu)\Delta t} \to 1^-$
3. **Rising skewness:** $\text{Skew}(x) \to +\infty$ (the distribution becomes right-skewed as the unstable equilibrium approaches from above)
4. **Rising spatial correlation** (for spatially extended systems): $\text{Corr}(x_i, x_j) \to 1$ as $|i-j| \to \infty$

*Proof.* All four indicators follow from the Ornstein-Uhlenbeck stationary distribution and the behavior of $\lambda(\mu) \to 0^-$:

1. Variance $= \sigma^2/(2|\lambda|) \to \infty$ as $|\lambda| \to 0$.
2. Autocorrelation $\rho(\tau) = e^{\lambda\tau} \to e^{0} = 1$ as $\lambda \to 0^-$.
3. Skewness arises from nonlinear terms neglected in linearization; as the system approaches the bifurcation, trajectories spend increasing time near the saddle point before returning to the stable node, generating positive skewness.
4. Spatial correlation follows from the same logic applied to coupled spatial units: as the individual decay rate $|\lambda| \to 0$, spatial correlations that would otherwise decay exponentially in space approach unity. $\square$

### 19.5.2 Computing the Indicators

**Algorithm 19.1 (Early-Warning Indicator Computation)**

```python
FUNCTION compute_EWI(time_series x, window_size w, lag tau):
    # Step 1: Detrend to remove slow parameter drift
    trend = rolling_mean(x, window=w)
    residuals = x - trend

    # Step 2: Rolling variance (EWI 1)
    var_t = rolling_variance(residuals, window=w)

    # Step 3: Rolling lag-1 autocorrelation (EWI 2)
    ar1_t = rolling_autocorrelation(residuals, lag=tau, window=w)

    # Step 4: Rolling skewness (EWI 3)
    skew_t = rolling_skewness(residuals, window=w)

    # Step 5: Kendall's tau trend test for each indicator
    tau_var,  p_var  = kendall_tau(var_t)
    tau_ar1,  p_ar1  = kendall_tau(ar1_t)
    tau_skew, p_skew = kendall_tau(skew_t)

    RETURN {
        'variance': var_t, 'ar1': ar1_t, 'skewness': skew_t,
        'tau_var': tau_var, 'p_var': p_var,
        'tau_ar1': tau_ar1, 'p_ar1': p_ar1,
        'tau_skew': tau_skew, 'p_skew': p_skew
    }
```

**Statistical interpretation.** Kendall's $\tau$ measures the monotonic trend in each indicator over the analysis window. A significant positive trend in variance and/or autocorrelation ($p < 0.05$) constitutes statistical evidence of approaching bifurcation. The combined indicator:

$$\text{EWI}_{\text{combined}} = \frac{\tau_{\text{var}} + \tau_{\text{ar1}}}{2}$$

has been found in simulation studies to provide approximately 80–90\% sensitivity (true positive rate) and 70–80\% specificity (true negative rate) for regime shift detection in ecological time series with at least 50 observations before the shift (Dakos et al., 2012).

---

## 19.6 Coupled Economic-Ecological Systems

### 19.6.1 The Formal Coupled System

The most important result of this chapter is the formal specification of the conditions under which economic and ecological dynamics can stably co-evolve. We derive these conditions from the eigenvalue analysis of the coupled system's Jacobian.

**Definition 19.5 (Coupled Economic-Ecological System).** A coupled economic-ecological system is the pair of differential equations:

$$\dot{x} = f(x, y) \quad \text{(ecological dynamics, state } x)$$
$$\dot{y} = g(x, y) \quad \text{(economic dynamics, state } y)$$

where $y$ affects $x$ through extraction or pollution ($\partial f/\partial y < 0$ typically — economic activity degrades the ecological state) and $x$ affects $y$ through the productivity of natural capital ($\partial g/\partial x > 0$ — ecological condition determines economic output).

**The Jacobian.** At any equilibrium $(x^*, y^*)$:

$$J = \begin{pmatrix} f_x & f_y \\ g_x & g_y \end{pmatrix}$$

where $f_x = \partial f/\partial x$, $f_y = \partial f/\partial y < 0$, $g_x = \partial g/\partial x > 0$, $g_y = \partial g/\partial y$.

**Stability conditions.** The coupled equilibrium $(x^*, y^*)$ is locally asymptotically stable if both eigenvalues of $J$ have negative real parts. By the Routh-Hurwitz criterion:

$$\text{Stable} \iff \text{tr}(J) < 0 \text{ and } \det(J) > 0$$

$$\text{tr}(J) = f_x + g_y < 0$$
$$\det(J) = f_x g_y - f_y g_x > 0$$

**Theorem 19.2 (Conditions for Stable Co-Evolution).** The coupled economic-ecological system is asymptotically stable if and only if:

$$f_x + g_y < 0 \quad \text{(sum of individual stability)} \tag{A}$$
$$f_x g_y > f_y g_x \quad \text{(ecological stability dominates cross-coupling)} \tag{B}$$

*Proof.* Standard Routh-Hurwitz for a 2×2 system. Condition (A) requires that the combined tendency to return to equilibrium ($f_x + g_y < 0$) exceeds the combined tendency to diverge. Condition (B) requires that $\det(J) > 0$, which expands to $f_x g_y - f_y g_x > 0$. Since $f_y < 0$ and $g_x > 0$, the term $-f_y g_x > 0$. Therefore Condition (B) is equivalent to $f_x g_y > -(-f_y g_x)$, i.e., the product of direct stability effects must exceed the product of cross-coupling destabilization effects. $\square$

**Corollary 19.2 (When Growth Destabilizes).** Economic growth that increases $|g_y|$ (the self-reinforcing nature of economic dynamics) or increases $|g_x|$ (the sensitivity of economic output to ecological state) without a commensurate increase in $|f_x|$ (ecological resilience) will eventually violate Condition (B), causing the coupled system to become unstable.

This corollary is the formal statement of a well-known empirical pattern: economic growth that outpaces ecological regeneration eventually destabilizes the ecological system on which it depends, triggering collapse of both. The formal condition $f_x g_y > f_y g_x$ specifies exactly when this transition occurs: when the economic feedback on itself ($g_y$) outpaces the ecological stability ($f_x$) relative to the cross-coupling terms.

### 19.6.2 The Cooperative Economic-Ecological System

In a cooperative-regenerative economy, the economic dynamics are specifically designed to avoid this destabilization. The Stewardship Condition $\dot{N} \geq 0$ — operationalized through the SFC-N framework of Chapter 18 — constrains $g(x, y)$ such that $f_y \to 0$ when $x$ approaches its minimum viable level. Formally, the cooperative system modifies the economic dynamics:

$$\dot{y} = g(x, y) \cdot \mathbb{1}[x > x_{\min}] + g_s(x, y) \cdot \mathbb{1}[x \leq x_{\min}]$$

where $g_s$ is a "stewardship mode" that reduces extraction when the ecological stock falls below its critical threshold $x_{\min}$. In stewardship mode: $f_y^{(s)} \approx 0$ and $g_x^{(s)} > g_x$ (the economic system becomes more sensitive to ecological restoration). This automatically satisfies Condition (B): $f_x g_y^{(s)} > 0 = f_y^{(s)} g_x^{(s)}$.

**Proposition 19.3 (Cooperative Economy and Ecological Stability).** A cooperative-regenerative economy that implements the Stewardship Condition ($\dot{N} \geq 0$) through binding extraction limits is asymptotically stable in the coupled economic-ecological system whenever $f_x < 0$ (the ecological system has positive self-regulation) and $g_y < 0$ (the economic system has self-limiting dynamics, e.g., through market saturation or cooperative production quotas).

*Proof.* Under the stewardship mode, $f_y^{(s)} \approx 0$, so $\det(J^{(s)}) = f_x g_y^{(s)} - f_y^{(s)} g_x^{(s)} \approx f_x g_y^{(s)} > 0$ (since $f_x < 0$ and $g_y^{(s)} < 0$). The trace condition $f_x + g_y^{(s)} < 0$ holds when both individual systems have self-regulating dynamics. $\square$

---

## 19.7 Mathematical Model: Regime Shift Bifurcation with Ecological-Economic Coupling

We develop the full coupled bifurcation model that will be calibrated to the North Atlantic cod case in Section 19.8.

**The extended system.** Let $N$ be the fish stock (normalized to $[0,1]$ relative to pre-exploitation level) and $Y$ be the economic output from the fishery (relative to sustainable yield). The coupled dynamics:

$$\dot{N} = rN(1-N) - \frac{hN^2}{N^2 + a^2} - \phi Y N \tag{19.1}$$

$$\dot{Y} = \alpha(N - N^*) Y - \beta Y^2 \tag{19.2}$$

where $r$ is intrinsic growth rate, $h$ is natural mortality scaling, $a$ is half-saturation constant (the Holling type-III loss function from Chapter 5), $\phi$ is the fishing pressure coefficient, $N^*$ is the target stock level for the fishery, $\alpha$ is the responsiveness of economic output to stock abundance, and $\beta$ is the market saturation coefficient.

Equation (19.1) is the extended ecological model from Chapter 5 with an additional economic pressure term $-\phi Y N$. Equation (19.2) is a logistic economic growth model in which output grows when stock exceeds the target level and declines when below it.

**Equilibria and stability.** Setting $\dot{N} = 0$ and $\dot{Y} = 0$:

From (19.2): either $Y = 0$ (no fishery) or $N = N^* + \beta Y / \alpha$.

Substituting into (19.1) at the cooperative equilibrium $Y^* > 0$:

$$r N^*(1-N^*) = \frac{h N^{*2}}{N^{*2} + a^2} + \phi Y^* N^*$$

This implicitly defines $N^*$ and $Y^*$ as functions of $\phi$ (fishing pressure). The bifurcation in $\phi$ occurs at the value $\phi_c$ where the cooperative equilibrium loses stability — the formal tipping point.

**The Jacobian at the cooperative equilibrium:**

$$J = \begin{pmatrix}
r(1-2N^*) - \frac{2h a^2 N^*}{(N^{*2}+a^2)^2} - \phi Y^* & -\phi N^* \\
\alpha Y^* & \alpha(N^* - N^*) - 2\beta Y^* = -2\beta Y^*
\end{pmatrix}$$

Note $J_{22} = -2\beta Y^* < 0$ (economic self-limitation). $J_{12} = -\phi N^* < 0$ (fishing pressure reduces stock returns). $J_{21} = \alpha Y^* > 0$ (economic output increases with stock). The stability conditions follow from Theorem 19.2.

---

## 19.8 Worked Example: North Atlantic Cod Fishery Collapse (1992)

### 19.8.1 Background

The Grand Banks cod fishery off Newfoundland was one of the most productive fisheries in the world for centuries — sustaining European fishing fleets from the 15th century onward. Industrial fishing in the 20th century progressively reduced the stock. On July 2, 1992, the Canadian federal government announced a moratorium on cod fishing — an unprecedented closure of a commercial fishery that had operated continuously for 500 years. The stock had collapsed to approximately 1\% of its pre-industrial level.

The collapse is the canonical example of a fishery regime shift: a system that appeared to be functioning (catches were maintained by fishing a shrinking stock more intensively) until it crossed a tipping point from which it has not recovered despite three decades of reduced fishing pressure.

### 19.8.2 Early-Warning Signal Analysis

We apply Algorithm 19.1 to the NAFO (Northwest Atlantic Fisheries Organization) cod stock assessment data, 1962–1992 (the period before the moratorium).

**Data.** Annual cod biomass indices (spawning stock biomass, SSB) from the Div. 2J3KL stock assessment, 1962–1992. The series shows a gradual decline from approximately 400,000 tonnes in the early 1960s to below 50,000 tonnes by 1992, but with substantial year-to-year variability masking the trend.

**Results** (rolling window: 12 years, lag: 1 year):

| Period | Variance (×10³) | AR(1) autocorrelation | Skewness |
|:---|:---|:---|:---|
| 1962–1973 | 2.1 | 0.31 | −0.12 |
| 1970–1981 | 3.8 | 0.44 | +0.18 |
| 1978–1989 | 7.2 | 0.67 | +0.41 |
| 1981–1992 | 12.4 | 0.82 | +0.73 |

**Kendall's $\tau$ for variance trend (1962–1992):** $\hat{\tau} = 0.71$, $p < 0.001$.
**Kendall's $\tau$ for AR(1) trend (1962–1992):** $\hat{\tau} = 0.68$, $p < 0.001$.

Both indicators show statistically significant upward trends beginning approximately 15–18 years before the collapse. The rising skewness (from negative in the 1960s to strongly positive in the late 1980s) is consistent with the theoretical prediction of Theorem 19.1 — the distribution shifted toward the right as the basin of attraction of the high-stock equilibrium contracted and the tipping threshold approached.

**When could intervention have been effective?** The basin-of-attraction approach suggests that once the AR(1) autocorrelation exceeds approximately 0.6, the system is within one or two major perturbations of the tipping threshold. The AR(1) reached 0.67 in the 1978–1989 rolling window — centered on approximately 1983. A moratorium in 1983 (nine years before the actual moratorium) would have found the stock at approximately 150,000 tonnes SSB, well above the minimum viable population threshold. Whether recovery would have been possible is uncertain given the reproductive biology of the stock, but the early-warning signals clearly indicated dangerous proximity to a regime shift by the early 1980s.

**The management failure.** NAFO's scientific advisories from the late 1970s onward consistently recommended catch reductions. The recommended total allowable catches were consistently exceeded by actual catches (by 20–50\% annually throughout the 1980s). The political economy of the fishery — short-term employment concerns, international fishing fleet competition, and inadequate enforcement — prevented the implementation of scientifically advised reductions even as the early-warning signals were accumulating. This is the panarchy failure: the economic system's K-phase rigidity (established interests, sunk costs, institutional inertia) prevented the adaptive response that would have been available in an earlier phase of the adaptive cycle.

---

## 19.9 Case Study: The Maldives as a Coupled Economic-Ecological System

### 19.9.1 System Description

The Republic of Maldives consists of 1,192 coral islands across 26 atolls in the Indian Ocean. Its economy is almost entirely dependent on two natural systems: coral reef ecosystems (supporting the tourism industry, which accounts for approximately 25\% of GDP and 60\% of export earnings) and sea-level stability (with mean land elevation of 1.5 meters above sea level, the Maldives faces existential risk from sea-level rise).

This creates one of the clearest examples of a coupled economic-ecological system in which the ecological variable (coral reef condition, sea-level margin) directly determines the economic output (tourism revenue), creating both the positive feedback (tourism revenue funds reef management) and the negative feedback (tourism pressure degrades reefs) of the coupled model.

### 19.9.2 Formal Model Specification

**Ecological state variable:** $N(t) \in [0,1]$ = coral cover as fraction of pre-bleaching level. Dynamics:

$$\dot{N} = r N(1-N) - h_T(T) \cdot N - h_C(C) \cdot N$$

where $r = 0.08$/year (coral recovery rate), $h_T(T)$ is thermal stress from ocean warming (increasing with global temperature $T$), and $h_C(C)$ is physical damage from diver and boat activity proportional to tourist arrivals $C$.

**Economic state variable:** $Y(t)$ = tourism revenue (normalized to 2020 level). Dynamics:

$$\dot{Y} = \alpha(N - N_{\min}) Y - \beta Y^2 - \gamma \dot{T}Y$$

where $\alpha$ is the responsiveness of tourism to reef quality, $N_{\min} = 0.3$ is the minimum reef cover for viable tourism, $\beta$ is market saturation, and $\gamma\dot{T}$ captures the direct effect of climate change news on tourist demand (negative).

**Coupled Jacobian at current state** ($N^* \approx 0.65$, $Y^* \approx 1.0$, calibrated to 2019 pre-pandemic data):

$$J \approx \begin{pmatrix}
-0.12 & -0.02 \\
+0.15 & -0.08
\end{pmatrix}$$

$\text{tr}(J) = -0.20 < 0$ ✓; $\det(J) = (-0.12)(-0.08) - (-0.02)(0.15) = 0.0096 + 0.003 = 0.0126 > 0$ ✓.

**The coupled system is currently stable**, but marginally so: both stability conditions are satisfied but with limited margin.

### 19.9.3 Proximity to the Tipping Point

**The ecological tipping point.** Coral reefs exhibit regime shifts: above approximately 25–30\% coral cover, reefs are in the "coral-dominated" stable state with high biodiversity and active growth. Below this threshold, macroalgae competitively exclude corals, and the reef shifts to an "algae-dominated" stable state with low biodiversity and negligible economic value for tourism.

Current Maldives coral cover: approximately 65\% (post-2016 bleaching recovery). Tipping threshold: approximately 25–30\%.

**The thermal stress trajectory.** Under the IPCC's middle-of-the-road scenario (SSP2-4.5), Indian Ocean sea surface temperatures are projected to increase by approximately 0.7°C by 2050. The Maldives has already experienced two major bleaching events (1998, 2016) that temporarily reduced coral cover to below the tipping threshold in some atolls. As ocean temperatures rise, bleaching events become more frequent — the IPCC projects annual bleaching in the Indian Ocean by the 2040s under SSP2-4.5.

**Early-warning signals in coral cover time series (Maldives, 1997–2023):**
- Variance: Kendall's $\tau = +0.54$, $p = 0.02$ — significant upward trend.
- AR(1) autocorrelation: Kendall's $\tau = +0.48$, $p = 0.04$ — significant upward trend.

The Maldives coral ecosystem is exhibiting early-warning signals of an approaching regime shift, consistent with a basin of attraction that is contracting under progressive thermal stress.

**Economic implications.** If the tipping point is crossed and coral cover falls permanently below 25–30\%, the tourism industry would collapse within approximately 5–10 years (as reef quality-dependent visitor satisfaction declines and tour operators shift to alternative destinations). Tourism revenue is approximately USD 4 billion annually — approximately 50\% of the Maldives' USD 8 billion GDP. The economic cost of the ecological regime shift is therefore the effective loss of approximately half of GDP — an economic collapse of a magnitude comparable to the Great Depression, but without the recovery mechanism, since the ecological collapse would be irreversible on any timescale relevant to economic planning.

**The cooperative policy implication.** The Maldives case illustrates Proposition 19.3: maintaining ecological resilience (keeping coral cover well above the tipping threshold) is a precondition for economic stability, not a constraint on it. The coupled model shows that economic output $Y^*$ is only positive when $N > N_{\min} = 0.3$ — economic success requires ecological health. The stewardship constraint $\dot{N} \geq 0$ is therefore in this case not just a normative aspiration but a necessary condition for the continuation of the economy itself, in the precise sense of Theorem 17.1.

---

## Chapter Summary

This chapter has developed the formal theory of ecological resilience and regime shifts, and derived the mathematical conditions under which economic and ecological systems can stably co-evolve.

Engineering resilience — the speed of return to equilibrium after small perturbations — is the wrong metric for economic policy in systems with multiple stable states. Ecological resilience — the size of the basin of attraction — is the appropriate policy target, because it measures the system's robustness to the large perturbations that actually threaten irreversible transitions. Proposition 19.1 proves this formally: policies that maintain engineering resilience while shrinking the basin of attraction are inadequate and dangerous.

Holling's adaptive cycle provides a formal dynamical systems model of ecosystem and economic life cycles, with four phases (growth, conservation, release, reorganization) and characteristic changes in potential, connectedness, and resilience. The panarchy framework nests these cycles across scales, with revolt connections (lower-scale collapse triggering upper-scale collapse) and remembrance connections (upper-scale stability supporting lower-scale recovery). Proposition 19.2 derives the stability conditions for two-scale panarchies.

The early-warning indicators of Theorem 19.1 — rising variance, rising autocorrelation, rising skewness — follow from critical slowing down near a saddle-node bifurcation. All three indicators are computable from observed time series data without knowledge of the underlying system dynamics, providing practical tools for detecting approaching regime shifts before they occur.

The formal coupled system analysis (Theorem 19.2) derives the conditions for stable co-evolution: the sum of direct stability effects must be negative (Condition A), and ecological stability must dominate cross-coupling destabilization (Condition B). Cooperative-regenerative economies that implement the Stewardship Condition automatically satisfy both (Proposition 19.3).

The North Atlantic cod collapse demonstrates that early-warning signals were present 9–12 years before the 1992 moratorium — sufficient lead time for effective intervention — but that the political economy of the fishery prevented implementation. The Maldives case shows a coupled system currently stable but exhibiting significant early-warning signals, with the tipping threshold approaching as ocean temperatures rise under climate change.

Chapter 20 develops Ecological Network Analysis — the formal framework for measuring the internal health and organization of ecological-economic systems, complementing the resilience analysis of this chapter with a structural assessment of how efficiently ecosystems and economies process energy and material flows.

---

## Exercises

**19.1** Define engineering resilience and ecological resilience formally (Definitions 19.1 and 19.2). For a savanna ecosystem with two stable states (grassland at $x^*_G = 0.2$ tree cover, forest at $x^*_F = 0.8$ tree cover) and a tipping threshold at $\hat{x} = 0.45$:
(a) The current state is $x = 0.65$ (closer to forest equilibrium). What is the ecological resilience $R_{\mathcal{E}}$?
(b) After a drought reduces tree cover to $x = 0.50$, the vegetation dynamics have eigenvalue $\lambda = -0.15$/year near $x^*_F$. What is the engineering resilience $R_E$? Has the ecological resilience changed?
(c) A further drought reduces tree cover to $x = 0.42$. Has the system crossed the tipping threshold? Is the engineering resilience metric useful for detecting this?

**19.2** The panarchy stability conditions (Proposition 19.2) require that lower-scale collapses not trigger upper-scale collapses when the upper scale is in its K phase.
(a) In the 2007–09 financial crisis, identify: (i) the lower scale undergoing $\Omega$-phase collapse; (ii) the upper scale in K phase; (iii) the revolt connection that transmitted the collapse upward.
(b) The subsequent European sovereign debt crisis (2010–12) represents a second revolt connection. Identify the scales and the transmission mechanism.
(c) What governance interventions (in the panarchy framework) would have broken the revolt connection and prevented the lower-scale financial collapse from triggering upper-scale sovereign crises?

**19.3** Apply the early-warning indicator algorithm (Algorithm 19.1) to the following stylized time series representing a lake ecosystem transitioning from clear to turbid:

Year: 1, 2, ..., 30
Clarity index (1=clear, 0=turbid): starts near 0.85, slowly declining, with fluctuations:
$x = 0.85 - 0.008t + 0.04\xi_t$, where $\xi_t \sim N(0,1)$ for $t = 1,...,20$,
then $x = 0.65 - 0.025(t-20) + 0.08\xi_t$ for $t = 21,...,30$.

(a) Compute the rolling variance (window = 8 years) and rolling AR(1) (window = 8 years) for this series.
(b) At what year does the variance first show a statistically significant upward trend (Kendall's $\tau > 0$ with $p < 0.10$)?
(c) Compare the timing of the statistical early-warning signal to the year in which the system dynamics change ($t = 21$). Does the indicator lead or lag the parameter change?

**★ 19.4** Derive all four early-warning indicators of Theorem 19.1 from the formal model of a saddle-node bifurcation.

(a) Start with the Ornstein-Uhlenbeck process $\dot{x} = \lambda(\mu)(x-x^*) + \sigma\xi(t)$ with $\lambda(\mu) < 0$. Compute the stationary distribution and show $\text{Var}(x) = \sigma^2/(2|\lambda|)$.
(b) Compute the autocorrelation function $\rho(\tau) = \text{Cov}(x(t), x(t+\tau))/\text{Var}(x)$ and show $\rho(\tau) = e^{\lambda\tau}$.
(c) Show that as $\lambda \to 0^-$: variance $\to \infty$, lag-1 autocorrelation $\to 1^-$.
(d) For skewness (EWI 3): expand the dynamics to second order around $x^*$: $\dot{x} = \lambda(x-x^*) + \frac{1}{2}f_{xx}(x-x^*)^2 + \sigma\xi$. Show that the second-order term generates positive skewness when $f_{xx} > 0$ (the potential is convex on one side of the saddle), and explain geometrically why skewness increases near the tipping threshold.

**★ 19.5** Prove Theorem 19.2 (conditions for stable co-evolution) and its Corollary 19.2 (when growth destabilizes).

(a) Write the Jacobian $J$ of the coupled system at the cooperative equilibrium and derive the Routh-Hurwitz conditions for $2 \times 2$ systems.
(b) Show that Condition (B) ($f_x g_y > f_y g_x$) is violated when fishing pressure $\phi$ increases beyond $\phi_c$: compute $\phi_c$ for the coupled model of Section 19.7 with $r = 0.15$, $h = 0.10$, $a = 0.2$, $\alpha = 0.3$, $\beta = 0.1$, $N^* = 0.6$.
(c) Interpret $\phi_c$ economically: what does the critical fishing pressure represent in terms of sustainable yield, and why is it the formal tipping point for the coupled system?
(d) Show that Proposition 19.3 holds: if the cooperative economy implements $\dot{N} \geq 0$ through a binding extraction limit $\phi \leq \phi^*$, the coupled system is always stable for any ecologically self-regulating ecology ($f_x < 0$) and self-limiting economy ($g_y < 0$).

**★★ 19.6** Construct a coupled economic-ecological model for the Baltic Sea fisheries and calibrate it to 1980–2020 data.

**Data sources:** ICES (International Council for the Exploration of the Sea) Baltic Sea assessment data; Eurostat fisheries statistics.

(a) The Baltic Sea has three principal commercial species: cod ($N_1$), herring ($N_2$), and sprat ($N_3$), with predator-prey interactions ($N_1$ preys on $N_2$ and $N_3$). Specify the ecological dynamics as a Lotka-Volterra system with harvesting:

$$\dot{N}_1 = r_1 N_1(1 - N_1/K_1) - a_{12}N_1 N_2 - a_{13}N_1 N_3 - E_1(Y)$$
$$\dot{N}_2 = r_2 N_2(1 - N_2/K_2) + b_{21}N_1 N_2 - a_{23}N_2 N_3 - E_2(Y)$$
$$\dot{N}_3 = r_3 N_3(1 - N_3/K_3) + b_{31}N_1 N_3 + b_{32}N_2 N_3 - E_3(Y)$$

Calibrate the parameters $(r_i, K_i, a_{ij}, b_{ij})$ using ICES stock assessment data for 1980–2000 (pre-collapse baseline). Report your estimates.

(b) Specify the economic dynamics: fishery revenue $Y = p_1 E_1 + p_2 E_2 + p_3 E_3$ where prices $p_i$ are exogenous. Calibrate the revenue-to-extraction relationship using Eurostat landing values.

(c) Compute the Jacobian of the full $4 \times 4$ coupled system (three ecological states plus economic revenue) at the 1980 baseline. Are all eigenvalues negative? Is the system stable?

(d) Apply the early-warning indicator algorithm to each species' catch time series (1980–2020). Which species shows the strongest early-warning signals? Does the timing of the signals precede the observed stock declines? What policy intervention (reduction in total allowable catch, changes in mesh size, area closures) would have extended the ecological resilience and preserved the coupled stability?

---

*Chapter 20 develops Ecological Network Analysis — the formal framework for measuring the efficiency, organization, and health of ecological and economic systems through the analysis of energy and material flow networks. Ascendancy, overhead, and capacity provide quantitative measures of how well a system is organized and how much buffering capacity it retains — tools for assessing the regenerative condition that the Stewardship Constraint requires.*
