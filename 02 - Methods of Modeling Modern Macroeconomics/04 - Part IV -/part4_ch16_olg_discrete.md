# Chapter 16: The Discrete-Time Overlapping Generations Model

*Social Security and Capital Accumulation*

> *"The overlapping-generations model is the essential tool for thinking about intergenerational issues: pension reform, public debt, and the distribution of welfare across generations."*
> — Peter Diamond, Nobel Lecture, 2010

**Cross-reference:** *Principles* Ch. 25.1 (Diamond OLG: dynamic efficiency, PAYG social security); Ch. 22.2 (Ricardian equivalence failure) **[P:Ch.25.1, P:Ch.22.2]**

---

## 16.1 The Diamond (1965) OLG Model: Structure and Motivation

Chapter 12 analyzed the **continuous-time** OLG model of Blanchard and Yaari, which achieves tractability by assuming a constant mortality rate. The **discrete-time Diamond (1965) model** is the more common benchmark in the macroeconomics literature. It features two-period-lived households — a young generation and an old generation coexist at each date — and a competitive production sector with capital and labor.

The OLG structure introduces a feature absent from the Ramsey model: **generational heterogeneity**. Different generations have different intertemporal budget constraints because they were born at different times and face different wage and interest rate environments. This heterogeneity creates scope for markets to fail in a specific way: the competitive equilibrium can be dynamically inefficient, accumulating more capital than is socially optimal. The government can improve welfare by taxing away excess saving — a policy dimension entirely absent from the representative-agent framework.

The Diamond model is also the foundation for the analysis of PAYG social security [P:Ch.25.1] and the fiscal theory of the price level in an OLG context [P:Ch.22.2]. This chapter derives the model completely, proves the dynamic efficiency conditions, and analyzes social security reform.

---

## 16.2 Households: Two-Period Optimization

Each household lives for two periods: young (period $t$) and old (period $t+1$). There is no uncertainty.

**Endowments:** The young earn wage $w_t$ in period $t$; the old earn nothing (they retire).

**Budget constraints:**

$$c^y_t + s_t = w_t \quad \text{(young: consume and save)}$$
$$c^o_{t+1} = (1+r_{t+1})s_t \quad \text{(old: consume all savings plus returns)}$$

where $c^y_t$ is young consumption, $c^o_{t+1}$ is old consumption, $s_t \geq 0$ is saving (borrowing is ruled out by the terminal age), and $r_{t+1}$ is the interest rate from period $t$ to $t+1$.

**Lifetime utility maximization:**

$$\max_{s_t \geq 0} u(w_t - s_t) + \beta u\bigl((1+r_{t+1})s_t\bigr).$$

**First-order condition (assuming interior solution):**

$$-u'(w_t - s_t) + \beta(1+r_{t+1})u'\bigl((1+r_{t+1})s_t\bigr) = 0,$$

$$\boxed{u'(c^y_t) = \beta(1+r_{t+1})u'(c^o_{t+1}).}$$

This is the discrete-time Euler equation — the same condition as in the infinite-horizon Ramsey model, now for a two-period problem.

**Definition 16.1 (Saving Function).** The **saving function** $s(w_t, r_{t+1})$ is the solution to the household's optimization problem. Under standard assumptions (strictly concave utility, intertemporal normality), $s$ is:
- Increasing in $w_t$ (higher income → more saving, if the income effect dominates).
- Ambiguously signed in $r_{t+1}$: higher returns increase the return to saving (substitution effect: save more) but also make any given saving target achievable with less saving (income effect). The sign depends on the elasticity of intertemporal substitution.

### 16.2.1 Closed-Form Solution: Log Utility

For $u(c) = \ln c$ and $\beta \in (0,1)$:

$$\frac{1}{c^y_t} = \beta(1+r_{t+1})\frac{1}{c^o_{t+1}} = \beta(1+r_{t+1})\frac{1}{(1+r_{t+1})s_t} = \frac{\beta}{s_t}.$$

So $c^y_t = s_t/\beta$. Combined with the young budget constraint $c^y_t + s_t = w_t$:

$$\frac{s_t}{\beta} + s_t = w_t \implies s_t\left(\frac{1+\beta}{\beta}\right) = w_t \implies \boxed{s_t = \frac{\beta}{1+\beta}w_t.}$$

Two remarkable properties of the log-utility saving function: (1) it is independent of $r_{t+1}$ (income and substitution effects cancel exactly); (2) it is a constant fraction $\beta/(1+\beta)$ of wage income.

---

## 16.3 Firms: Competitive Production

A representative firm produces output using capital and labor with a CRS Cobb–Douglas technology:

$$Y_t = K_t^\alpha (A_t L_t)^{1-\alpha}, \quad A_t = (1+g)^t A_0.$$

In intensive form ($\tilde{k}_t = K_t/(A_tL_t)$): $\tilde{y}_t = \tilde{k}_t^\alpha$.

Competitive factor markets imply:

$$r_t = f'(\tilde{k}_t) - \delta = \alpha\tilde{k}_t^{\alpha-1} - \delta, \qquad \tilde{w}_t \equiv \frac{w_t}{A_t} = f(\tilde{k}_t) - \tilde{k}_tf'(\tilde{k}_t) = (1-\alpha)\tilde{k}_t^\alpha.$$

---

## 16.4 General Equilibrium: The Capital Accumulation Map

**Market clearing:** The capital stock at $t+1$ equals the aggregate savings of the young generation at $t$:

$$K_{t+1} = L_t s_t,$$

since the old dissave entirely ($s_t$ was their only saving, now returned as $(1+r_{t+1})s_t$). In intensive form ($\tilde{k}_{t+1} = K_{t+1}/(A_{t+1}L_{t+1})$):

$$\tilde{k}_{t+1} = \frac{L_t s_t}{A_{t+1}L_{t+1}} = \frac{s_t}{(1+g)(1+n)A_t} = \frac{s(\tilde{w}_t, r_{t+1})}{(1+g)(1+n)}.$$

**Definition 16.2 (Equilibrium Map).** With log utility: $s(\tilde{w}_t, r_{t+1}) = \frac{\beta}{1+\beta}\tilde{w}_t \cdot A_t$. Therefore:

$$\tilde{k}_{t+1} = \frac{\beta}{(1+\beta)(1+g)(1+n)}\tilde{w}_t = \frac{\beta(1-\alpha)}{(1+\beta)(1+g)(1+n)}\tilde{k}_t^\alpha \equiv \phi(\tilde{k}_t).$$

**Theorem 16.1 (Equilibrium Dynamics — Log Utility).** Under log utility and Cobb–Douglas production, the dynamics of the normalized capital stock are governed by the map:

$$\tilde{k}_{t+1} = \phi(\tilde{k}_t) = \frac{\beta(1-\alpha)}{(1+\beta)\mu_0}\tilde{k}_t^\alpha, \quad \mu_0 \equiv (1+g)(1+n).$$

This is a one-dimensional increasing map on $\mathbb{R}_+$.

### 16.4.1 Steady State

The steady state satisfies $\tilde{k}^* = \phi(\tilde{k}^*)$:

$$\tilde{k}^* = \frac{\beta(1-\alpha)}{(1+\beta)\mu_0}\tilde{k}^{*\alpha} \implies \tilde{k}^{*1-\alpha} = \frac{\beta(1-\alpha)}{(1+\beta)\mu_0}.$$

$$\boxed{\tilde{k}^* = \left[\frac{\beta(1-\alpha)}{(1+\beta)(1+g)(1+n)}\right]^{1/(1-\alpha)}.}$$

**Stability:** $\phi'(\tilde{k}^*) = \alpha\phi(\tilde{k}^*)/\tilde{k}^* = \alpha\tilde{k}^{*\alpha-1}\cdot\phi(1) = \alpha \in (0,1)$. Since $\phi'(\tilde{k}^*) = \alpha < 1$, the steady state is **locally asymptotically stable** — the economy converges monotonically to $\tilde{k}^*$ from any initial condition. (The map $\phi$ is globally stable on $\mathbb{R}_{++}$ because it is increasing and concave with a unique crossing of the 45-degree line.)

### 16.4.2 Multiplicity with Non-Log Utility

For general CRRA utility with $\sigma \neq 1$, the saving function $s(w, r)$ depends on $r$, which itself depends on $\tilde{k}_{t+1}$ (a forward-looking dependence). The equilibrium map is then an implicit equation:

$$\tilde{k}_{t+1} = \frac{s(\tilde{w}(\tilde{k}_t), r(\tilde{k}_{t+1}))}{(1+g)(1+n)}.$$

This implicit equation can have multiple solutions for $\tilde{k}_{t+1}$ given $\tilde{k}_t$ when $\sigma < 1$ (high EIS, strong substitution effect) — multiple equilibria. The Diamond model can exhibit **multiplicity of steady states** and **indeterminate dynamics** for certain parameter configurations, unlike the Ramsey model which has a unique globally stable steady state.

---

## 16.5 Dynamic Efficiency and Inefficiency

**Definition 16.3 (Dynamic Efficiency).** An allocation $\{c^y_t, c^o_t, \tilde{k}_t\}$ is **dynamically efficient** if there is no feasible reallocation that makes every generation better off. The competitive equilibrium is dynamically inefficient iff the economy is **over-accumulating capital** — investing more than the Golden Rule.

**Theorem 16.2 (Dynamic Inefficiency Criterion).** The steady state $\tilde{k}^*$ is dynamically inefficient (over-accumulated capital) if and only if:

$$r^* < n + g, \quad \text{i.e.,} \quad f'(\tilde{k}^*) - \delta < n + g.$$

Equivalently: the net marginal product of capital is less than the population-plus-technology growth rate.

*Proof.* The Golden Rule capital stock $\tilde{k}^{GR}$ maximizes $\tilde{c}^* = f(\tilde{k}) - (n+g+\delta)\tilde{k}$ (steady-state consumption per effective worker), satisfying $f'(\tilde{k}^{GR}) = n+g+\delta$.

If $\tilde{k}^* > \tilde{k}^{GR}$: the economy is to the right of the Golden Rule. Reducing capital (reducing $s$) would raise steady-state consumption and be Pareto-improving across all generations. Hence $\tilde{k}^* > \tilde{k}^{GR}$ implies dynamic inefficiency.

$\tilde{k}^* > \tilde{k}^{GR}$ iff $f'(\tilde{k}^*) < f'(\tilde{k}^{GR}) = n+g+\delta$, i.e., $r^* = f'(\tilde{k}^*)-\delta < n+g$. $\square$

**Is the Diamond OLG model dynamically inefficient in practice?** Abel, Mankiw, Summers, and Zeckhauser (1989) test dynamic efficiency by comparing aggregate profits ($\approx r^*K^*$) to aggregate investment ($\approx (n+g+\delta)K^*$). For all OECD economies, aggregate profits substantially exceed investment — implying $r^* > n+g$ and dynamic efficiency. The Diamond model's over-accumulation is a theoretical possibility that empirical data suggest is not realized.

**Why can the Diamond OLG over-accumulate when the Ramsey model cannot?** In the Ramsey model, the competitive equilibrium coincides with the social planner's solution (First Welfare Theorem), which never accumulates beyond the Golden Rule because impatient households optimize. In the Diamond OLG model, the First Welfare Theorem does not apply because markets are incomplete: the young and the unborn cannot trade with each other before they are born. This market incompleteness allows the equilibrium to reach the over-accumulation region.

---

## 16.6 PAYG Social Security and Welfare Analysis

A **Pay-As-You-Go (PAYG) social security system** taxes the current young at rate $\tau$ and transfers the proceeds to the current old:

$$d_t = \tau w_t L_t / L_{t-1} = \tau w_t (1+n).$$

The household's budget constraints with social security:

$$c^y_t + s_t = (1-\tau)w_t \quad \text{(young: after-tax wage)}$$
$$c^o_{t+1} = (1+r_{t+1})s_t + d_{t+1} = (1+r_{t+1})s_t + \tau w_{t+1}(1+n) \quad \text{(old: savings plus transfer)}$$

With log utility, the optimal saving becomes:

$$s_t = \frac{\beta}{1+\beta}(1-\tau)w_t - \frac{\tau w_{t+1}(1+n)}{(1+\beta)(1+r_{t+1})}.$$

The second term is the present value of the social security transfer discounted at $r_{t+1}$. The PAYG system reduces private saving dollar-for-dollar if $r^* = n$ (where the implicit return on PAYG contributions equals the market return), less than dollar-for-dollar if $r^* > n$ (which the system crowds out more valuable private saving), and more than dollar-for-dollar if $r^* < n$ (dynamically inefficient case).

**Capital accumulation with PAYG:** The equilibrium map becomes:

$$\tilde{k}_{t+1} = \frac{s_t}{(1+g)(1+n)} < \frac{\beta(1-\alpha)\tilde{k}_t^\alpha}{(1+\beta)\mu_0} \quad \text{(for } \tau > 0\text{)}.$$

Social security reduces the steady-state capital stock. **If the economy is dynamically inefficient** ($r^* < n$), this reduction is welfare-improving: the economy was over-accumulating, and PAYG corrects the inefficiency by transferring resources from young savers to old consumers. **If the economy is dynamically efficient** ($r^* > n$), PAYG reduces welfare-enhancing private capital formation.

**Generation-by-generation welfare analysis:**

| Generation | Effect of introducing PAYG |
|---|---|
| First old generation | Gains: receives transfer $d$ without having contributed |
| Middle generations (steady state) | Ambiguous: lose $r^*$ on reduced saving, gain $n$ implicit return |
| Future generations (new SS) | Gain iff $r^* < n$ (dynamic inefficiency corrected); lose iff $r^* > n$ |

---

## 16.7 The Recursive Competitive Equilibrium

**Definition 16.4 (Recursive Competitive Equilibrium).** A **recursive competitive equilibrium** (RCE) for the Diamond OLG model consists of:

1. A **value function** $V(\tilde{k}, \text{age})$ for each household type.
2. **Policy functions**: $c^y(\tilde{k})$, $s(\tilde{k})$, $c^o(\tilde{k})$.
3. **Price functions**: $r(\tilde{k})$, $w(\tilde{k})$.
4. A **transition function**: $\tilde{k}' = \Phi(\tilde{k})$.

such that: (i) households optimize taking prices as given; (ii) firms optimize (factor markets clear); (iii) goods market clears: $Y_t = c^y_t L_t + c^o_t L_{t-1} + K_{t+1} + \delta K_t$; (iv) the transition function $\Phi$ is consistent with savings behavior: $\tilde{k}' = s(\tilde{k})/\mu_0$.

The RCE reduces the economy's dynamics to the single map $\tilde{k}_{t+1} = \Phi(\tilde{k}_t)$, which we can analyze with the tools of Chapters 14 and this chapter.

---

## 16.8 Worked Example: Full Algebraic Solution and Phase Diagram

*Cross-reference: Principles Ch. 25.1 (Diamond model)* **[P:Ch.25.1]**

**Log utility, Cobb–Douglas, Calibration:**

$\alpha = 0.36$, $\beta = 0.5$ (one period = 30 years, so $\beta = 0.5$ implies annual discount of about $0.978^{30}$), $n = 0.3$ (30-year population growth), $g = 0.5$ (30-year technology growth), $\delta = 1$ (capital fully depreciates in 30 years).

$\mu_0 = (1+n)(1+g) = 1.3\times1.5 = 1.95$.

**Steady state:**

$\tilde{k}^* = \left[\frac{0.5\times0.64}{1.5\times1.95}\right]^{1/0.64} = \left[\frac{0.32}{2.925}\right]^{1.5625} = [0.1094]^{1.5625} = 0.0359$

$r^* = 0.36\times 0.0359^{-0.64} - 1 = 0.36\times 28.6 - 1 = 9.3 - 1 = 8.3$ per 30 years → annual rate $\approx 7.7\%$

$n+g = 0.3 + 0.5 = 0.8$ per 30 years → annual rate $\approx 2\%$

Since $r^* = 8.3 \gg n+g = 0.8$, the economy is dynamically efficient. PAYG social security reduces welfare in this calibration.

```apl
⍝ APL — Diamond OLG: map simulation and phase diagram
⎕IO←0 ⋄ ⎕ML←1

alpha←0.36  ⋄  beta←0.5  ⋄  n←0.3  ⋄  g←0.5  ⋄  delta←1
mu0 ← (1+n)×1+g

⍝ Factor prices
wage   ← {(1-alpha)×⍵*alpha}        ⍝ w̃(k̃) = (1-α)k̃^α
rate   ← {alpha×⍵*alpha-1) - delta} ⍝ r(k̃) = αk̃^(α-1) - δ

⍝ Saving function (log utility: independent of r)
saving ← {(beta÷1+beta)×wage ⍵}     ⍝ s̃(k̃) = β/(1+β) × w̃(k̃)

⍝ Equilibrium map
phi ← {(saving ⍵)÷mu0}              ⍝ k̃' = s̃(k̃)/μ₀

⍝ Steady state: k̃^(1-α) = β(1-α)/((1+β)μ₀)
kstar_formula ← ((beta×1-alpha)÷(1+beta)×mu0)*÷1-alpha
kstar_formula    ⍝ steady state

⍝ Verify: phi(kstar) = kstar
phi kstar_formula    ⍝ should equal kstar_formula

⍝ Simulate from k0 = 0.01 (below steady state)
k0 ← 0.01
path ← phi \ 30 ⍴ k0    ⍝ scan: 30 iterations of the map
path[0] path[14] path[29]    ⍝ converges monotonically to k*

⍝ Phase diagram: plot phi(k) vs k on [0, 0.1]
k_grid ← (⍳100)×0.001
k_next ← phi¨ k_grid
⍝ Fixed point: where k_next = k_grid
⍝ Dynamic efficiency: r* vs n+g
rstar  ← rate kstar_formula
rstar (n+g)    ⍝ rstar >> n+g → dynamically efficient

⍝ Effect of PAYG social security (tau=0.1)
tau ← 0.1
⍝ With PAYG, saving changes: s_payg = β/(1+β) × (1-τ)w̃
saving_payg ← {(beta÷1+beta)×(1-tau)×wage ⍵}
phi_payg    ← {(saving_payg ⍵)÷mu0}
kstar_payg  ← (phi_payg ⍣ (1e¯8∘>|⊢-phi_payg)) ⊢ kstar_formula × 0.9
kstar_payg    ⍝ lower than kstar_formula (PAYG crowds out capital)
```

```python
import numpy as np; import matplotlib.pyplot as plt
from scipy.optimize import fixed_point

alpha, beta, n, g, delta = 0.36, 0.5, 0.3, 0.5, 1.0
mu0 = (1+n)*(1+g)

wage  = lambda k: (1-alpha)*k**alpha
rate  = lambda k: alpha*k**(alpha-1) - delta
phi   = lambda k: (beta/(1+beta))*wage(k)/mu0  # equilibrium map (log utility)
kstar = ((beta*(1-alpha))/((1+beta)*mu0))**(1/(1-alpha))

print(f"k̃* = {kstar:.5f}")
print(f"r* = {rate(kstar):.3f}  vs  n+g = {n+g:.3f}  → {'efficient' if rate(kstar)>n+g else 'inefficient'}")

# Phase diagram
k_grid = np.linspace(0.001, 0.08, 300)
k_next = phi(k_grid)

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(k_grid, k_next, 'b-', linewidth=2, label='k̃\' = φ(k̃)')
ax.plot(k_grid, k_grid, 'k--', label='45° line')
ax.axvline(kstar, color='red', linestyle=':', label=f'k̃*={kstar:.4f}')

# PAYG experiment
for tau in [0.0, 0.1, 0.2]:
    phi_tau = lambda k, t=tau: (beta/(1+beta))*(1-t)*wage(k)/mu0
    kstar_tau = fixed_point(phi_tau, kstar*0.9)
    ax.plot([kstar_tau], [phi_tau(kstar_tau)], 'o', markersize=8, label=f'τ={tau}: k̃*={kstar_tau:.4f}')

ax.set_xlabel('k̃ₜ'); ax.set_ylabel('k̃ₜ₊₁')
ax.set_title('Diamond OLG: Phase Diagram and PAYG Effects')
ax.legend(); plt.tight_layout(); plt.show()

# Simulate convergence
path = [0.01]; [path.append(phi(path[-1])) for _ in range(30)]
print(f"\nPath: {[round(x,5) for x in path[:5]]} ... {round(path[-1],5)} → k̃*={kstar:.5f}")
```

```julia
alpha, beta, n, g, delta = 0.36, 0.5, 0.3, 0.5, 1.0
mu0 = (1+n)*(1+g)

wage(k) = (1-alpha)*k^alpha
phi(k)  = (beta/(1+beta))*wage(k)/mu0
kstar   = ((beta*(1-alpha))/((1+beta)*mu0))^(1/(1-alpha))
rstar   = alpha*kstar^(alpha-1) - delta

println("k̃* = $(round(kstar,digits=5))")
println("r* = $(round(rstar,digits=3)) vs n+g = $(n+g) → $(rstar>n+g ? "efficient" : "INEFFICIENT")")
println("\nStability: φ'(k*) = α = $alpha < 1 → stable ✓")

# Simulate
k = 0.01; path = [k]
for _ in 1:40; k = phi(k); push!(path,k); end
println("Converges to: $(round(path[end],digits=5)) (k*=$(round(kstar,digits=5)))")
```

```r
alpha<-0.36; beta<-0.5; n<-0.3; g<-0.5; delta<-1.0; mu0<-(1+n)*(1+g)
wage <- function(k) (1-alpha)*k^alpha
phi  <- function(k) (beta/(1+beta))*wage(k)/mu0
kstar <- ((beta*(1-alpha))/((1+beta)*mu0))^(1/(1-alpha))
rstar <- alpha*kstar^(alpha-1) - delta
cat(sprintf("k̃*=%.5f, r*=%.3f, n+g=%.1f → %s\n",
            kstar, rstar, n+g, ifelse(rstar>n+g,"efficient","INEFFICIENT")))

# Path simulation via Reduce
path <- Reduce(function(k,.) phi(k), 1:40, kstar*0.3, accumulate=TRUE)
cat(sprintf("Path: %.5f → %.5f → %.5f → ... → %.5f (k*=%.5f)\n",
            path[1],path[2],path[3],tail(path,1),kstar))
```

---

## 16.9 Programming Exercises

### Exercise 16.1 (APL — Stability Analysis)

Write a dfn `olg_stability ← {alpha beta n g delta ← ⍵ ⋄ ...}` that returns: (a) $\tilde{k}^*$; (b) $r^*$; (c) $\phi'(\tilde{k}^*) = \alpha$ (always, confirming stability for log utility); (d) a Boolean `dynamic_efficient` flag; (e) the effect of $\tau = 0.1$ PAYG on $\tilde{k}^*$ as a percentage change. Test with the calibration of Section 16.8 and with a hypothetical economy where $r^* < n+g$ (dynamically inefficient — try very high $\beta$, low $r$).

### Exercise 16.2 (Python — CRRA Utility, Numerical Map)

For CRRA utility with $\sigma \in \{0.5, 1, 2, 5\}$: (a) compute the saving function $s(\tilde{w}, r)$ numerically (by solving the household FOC); (b) write the implicit map $\tilde{k}_{t+1} = s(\tilde{w}(\tilde{k}_t), r(\tilde{k}_{t+1}))/\mu_0$; (c) solve for $\tilde{k}^*$ using fixed-point iteration; (d) check whether multiple steady states can arise for $\sigma < 1$. Plot the map for each $\sigma$ value.

### Exercise 16.3 (Julia — PAYG Welfare Calculation)

For the log-utility calibration, compute the welfare effect of introducing a small PAYG system ($\tau = 0.01$). (a) For the current old (initial generation): compute $\Delta U = \ln(d_0) - \ln(0) = +\infty$ (they get a windfall — the first generation always gains). (b) For a middle generation (in the new steady state): $\Delta U = \ln(c^{y,new}) + \beta\ln(c^{o,new}) - [\ln(c^{y,old}) + \beta\ln(c^{o,old})]$. Is this positive or negative for the calibration (where $r^* > n$)? (c) Plot $\Delta U$ as a function of the PAYG tax rate $\tau$ for middle generations and identify the welfare-maximizing $\tau$.

### Exercise 16.4 — Fully Funded vs. PAYG ($\star$)

A **fully funded (FF)** social security system collects contributions $\tau w_t$ from the young, invests them at the market rate $r_{t+1}$, and returns $(1+r_{t+1})\tau w_t$ to the old at retirement. (a) Show that a fully funded system is equivalent to a mandated private saving account — it does not affect the equilibrium capital stock or welfare. (b) By contrast, show that PAYG provides an implicit return of $n+g$ (not $r^*$). (c) When is PAYG welfare-superior to FF? Show this is exactly when $n+g > r^*$ — the dynamic inefficiency condition.

---

## 16.10 Chapter Summary

**Key results:**

- In the **Diamond OLG model**, young households save $s(w_t, r_{t+1})$ from their wage income; old households dissave entirely. Capital accumulation is $K_{t+1} = L_t s_t$.
- With **log utility and Cobb–Douglas production**: $s_t = [\beta/(1+\beta)]w_t$ (independent of $r$); the equilibrium map is $\tilde{k}_{t+1} = [\beta(1-\alpha)/((1+\beta)\mu_0)]\tilde{k}_t^\alpha$, with unique stable steady state $\tilde{k}^* = [\beta(1-\alpha)/((1+\beta)\mu_0)]^{1/(1-\alpha)}$.
- The steady state is **dynamically efficient** iff $r^* > n+g$ (net MPK exceeds population-plus-growth rate); empirically, all OECD economies satisfy this.
- **PAYG social security** reduces $\tilde{k}^*$ by crowding out private saving. It improves welfare iff the economy is dynamically inefficient ($r^* < n+g$).
- The **Recursive Competitive Equilibrium** is the map $\Phi: \tilde{k}_t \mapsto \tilde{k}_{t+1}$ consistent with household optimization and market clearing.
- In APL: the map is `phi ← {(saving ⍵)÷mu0}`; simulation is `phi \ T ⍴ k0` (scan); the steady state is `phi ⍣ (1e¯8∘>|⊢-phi) ⊢ k0_init`.

*Next: Chapter 17 — Recursive Methods for Real Business Cycle Models*
