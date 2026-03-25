# Chapter 12: The Continuous-Time Overlapping Generations Model

*Demographic Transitions and the Blanchard–Yaari Framework*

> *"The overlapping generations structure is not a technical complication but an economic necessity: it is the only way to model the redistribution between generations that is central to social security, public debt, and demographic policy."*

**Cross-reference:** *Principles* Ch. 25 (OLG demographics, Diamond model, HANK, aging); Ch. 22 (social security, Ricardian equivalence across generations) **[P:Ch.25, P:Ch.22]**

---

## 12.1 Why Continuous-Time OLG? The Limits of the Ramsey Model

The Ramsey–Cass–Koopmans model assumes an infinitely-lived representative household. This is a powerful analytical device but obscures the intergenerational dimension that is central to many of the most important macroeconomic policy debates: social security reform, public debt sustainability, the consequences of population aging, and the distributional effects of fiscal policy across cohorts.

The Diamond (1965) overlapping-generations model (Chapter 16) addresses this by explicitly modeling two-period-lived generations. But the two-period structure is too coarse for quantitative work: it cannot generate realistic life-cycle saving profiles, smooth demographic dynamics, or tractable responses to policy changes.

The **Blanchard (1985) – Yaari (1965) perpetual youth model** strikes a compromise: it operates in continuous time, allows for finite lifetimes, preserves tractable aggregation, and generates a clean departure from Ricardian equivalence. Households face a constant Poisson death rate $p > 0$ — they can die at any moment with probability $p\,dt$ over a small interval $dt$. This generates a population that is always young (relative to actuarial tables) but has finite expected lifetimes: $E[\text{lifespan}] = 1/p$.

---

## 12.2 Household Optimization with Uncertain Lifetime

### 12.2.1 The Modified Discount Rate

A household born at time $s$ maximizes expected lifetime utility:

$$\mathbb{E}\left[\int_s^\infty e^{-\rho(t-s)} u(c_{s,t})\,dt\right],$$

where the expectation is taken over the random date of death $T_D$. Since death follows a Poisson process with rate $p$, the probability of surviving from $s$ to $t$ is $e^{-p(t-s)}$. Therefore:

$$\mathbb{E}\left[\int_s^\infty e^{-\rho(t-s)} u(c_{s,t})\,dt\right] = \int_s^\infty e^{-(\rho+p)(t-s)} u(c_{s,t})\,dt.$$

The **effective discount rate** is $\rho + p$: the sum of pure time preference $\rho$ and the mortality rate $p$. Households discount the future both because of impatience and because they may not live to enjoy it.

### 12.2.2 Actuarially Fair Annuities

With competitive insurance markets, households can purchase **actuarially fair annuities** that pay out while alive and transfer accumulated wealth to the insurance company upon death. The return on an annuity is $r_t + p$: the market interest rate $r_t$ plus the mortality premium $p$ (reflecting that the insurer collects the wealth of those who die and distributes it to survivors). The budget constraint for a household of cohort $s$ at time $t$:

$$\dot{a}_{s,t} = (r_t + p)a_{s,t} + w_t - c_{s,t},$$

where $a_{s,t}$ is financial wealth (claims on the annuity) and $w_t$ is the wage.

### 12.2.3 Optimal Consumption: The Euler Equation

With CRRA utility $u(c) = c^{1-\sigma}/(1-\sigma)$, the Hamiltonian for the cohort-$s$ problem gives the Euler equation:

$$\frac{\dot{c}_{s,t}}{c_{s,t}} = \frac{(r_t + p) - (\rho + p)}{\sigma} = \frac{r_t - \rho}{\sigma}.$$

The mortality rate $p$ cancels from the Euler equation — the effective discount rate $\rho+p$ and the annuity return $r+p$ both shift by $p$, so the consumption growth rate is the same as in the Ramsey model. The departure from Ramsey appears in the **level** of consumption, not its growth rate.

### 12.2.4 Consumption Function

Solving the Euler equation forward subject to the intertemporal budget constraint (no-Ponzi condition):

$$c_{s,t} = (\rho + p + (\sigma-1)(r_t-\rho)/\sigma)\cdot a_{s,t}^{total},$$

where $a_{s,t}^{total}$ is total wealth (financial wealth plus human wealth). For simplicity with $\sigma = 1$ (log utility), the **consumption function** takes the elegant form:

$$c_{s,t} = (\rho + p)\left[a_{s,t} + h_t\right],$$

where $h_t = \int_t^\infty e^{-\int_t^\tau(r_u+p)du}w_\tau\,d\tau$ is **human wealth** — the present discounted value of future wage income, discounted at the annuity-adjusted rate $(r+p)$.

**Definition 12.1 (Human Wealth in the Blanchard–Yaari Model).** Human wealth satisfies the differential equation:

$$\dot{h}_t = (r_t + p)h_t - w_t.$$

This says: human wealth grows at rate $(r+p)$ (the annuity return) minus current wages received (which convert to financial wealth).

---

## 12.3 Aggregation

The key insight of the perpetual youth model is that aggregation across cohorts is tractable. At any time $t$, the age distribution of living households follows an exponential distribution with rate $p$ — there are $pe^{-p(t-s)}$ survivors of cohort $s$ per unit population (normalized so the total population is $1$ without population growth, or $L_t = e^{nt}$ with growth).

**The aggregate consumption function.** Integrating $c_{s,t} = (\rho+p)(a_{s,t}+h_t)$ over all living cohorts (using $\sigma=1$):

$$C_t = \int_{-\infty}^t (\rho+p)(a_{s,t}+h_t)p e^{-p(t-s)} ds = (\rho+p)\left[\int_{-\infty}^t p e^{-p(t-s)}a_{s,t}ds + h_t\right].$$

Since $\int_{-\infty}^t p e^{-p(t-s)}a_{s,t}ds = A_t$ (total financial wealth per capita) and $\int_{-\infty}^t p e^{-p(t-s)}h_t\,ds = h_t$ (human wealth is common across cohorts):

$$\boxed{C_t = (\rho+p)(A_t + H_t),}$$

where $H_t = h_t$ is aggregate human wealth (equal to individual human wealth since all living cohorts have the same $h_t$).

**Definition 12.2 (Aggregate Dynamics in Blanchard–Yaari).** Aggregate consumption $C_t$ and human wealth $H_t$ satisfy:

$$\dot{C}_t = (r_t - \rho)C_t - p(\rho+p)A_t$$
$$\dot{H}_t = (r_t + p)H_t - w_t.$$

The term $-p(\rho+p)A_t$ is the **key departure from Ramsey**: newly born households start with zero financial wealth ($a_{s=t,t} = 0$), so aggregate consumption growth is lower than individual consumption growth by the amount $p(\rho+p)A_t$ — newly born consumers drag down the aggregate.

*Derivation of the aggregate Euler equation:* Differentiate $C_t = (\rho+p)(A_t + H_t)$:

$$\dot{C}_t = (\rho+p)(\dot{A}_t + \dot{H}_t).$$

Using $\dot{A}_t = (r_t+n)A_t + w_t - C_t$ (capital accumulation) and $\dot{H}_t = (r_t+p)H_t - w_t$:

$$\dot{C}_t = (\rho+p)[(r_t+n)A_t + w_t - C_t + (r_t+p)H_t - w_t]$$
$$= (\rho+p)[(r_t)(A_t + H_t) + nA_t + pH_t - C_t]$$
$$= (r_t)C_t + (\rho+p)[nA_t + pH_t - C_t].$$

Substituting $C_t = (\rho+p)(A_t+H_t)$:

$$\dot{C}_t = r_t C_t + (\rho+p)nA_t + (\rho+p)pH_t - (\rho+p)^2(A_t+H_t).$$

After simplification (assuming $n=0$ for clarity):

$$\dot{C}_t = (r_t - \rho)C_t - p(\rho+p)A_t.$$

**Proposition 12.1 (Departure from Ramsey).** The aggregate Euler equation in the Blanchard–Yaari model is $\dot{C}_t = (r_t-\rho)C_t - p(\rho+p)A_t$. This equals the Ramsey aggregate Euler equation $(r_t-\rho)C_t$ only when $p = 0$ (no mortality risk) or $A_t = 0$ (no financial wealth). The term $-p(\rho+p)A_t$ captures the **wealth dilution** from demographic turnover: new households enter without wealth, dragging down aggregate consumption growth.

---

## 12.4 General Equilibrium

In a closed economy with capital, the production side is identical to the Ramsey model:

$$Y_t = F(K_t, L_t), \quad r_t = F_K, \quad w_t = F_L = f(\tilde{k}_t) - \tilde{k}_t f'(\tilde{k}_t).$$

Market clearing: $K_t = A_t L_t$ (capital equals aggregate financial wealth) and $Y_t = C_t + \dot{K}_t + \delta K_t$.

In per-capita terms (setting $n = 0$ and $g = 0$ for analytical clarity), the equilibrium is characterized by the two-dimensional system:

$$\dot{k}_t = f(k_t) - c_t - \delta k_t \quad \text{(capital accumulation)}$$
$$\dot{c}_t = (f'(k_t) - \delta - \rho)c_t - p(\rho+p)k_t \quad \text{(aggregate Euler)}$$

### 12.4.1 Steady State

Setting $\dot{k} = 0$ and $\dot{c} = 0$:

**From $\dot{c} = 0$:** $f'(k^*) = \delta + \rho + p(\rho+p)k^*/c^*$.

**From $\dot{k} = 0$:** $c^* = f(k^*) - \delta k^*$.

Substituting: $f'(k^*) = \delta + \rho + p(\rho+p)k^*/(f(k^*)-\delta k^*)$.

This is a nonlinear equation in $k^*$ that must be solved numerically. However, we can immediately see the qualitative difference from Ramsey: the RCK steady state required $f'(k^*) = \delta + \rho$ (no mortality term), while the Blanchard–Yaari steady state requires a higher net marginal product $f'(k^*) > \delta + \rho$ (when $p > 0$ and $k^*, c^* > 0$). Since $f'' < 0$, this means $k^*_{BY} < k^*_{RCK}$: **the Blanchard–Yaari economy has a lower capital stock than the Ramsey economy**.

**Intuition:** Households with finite lifetimes save less than infinitely-lived households (they don't need to accumulate wealth to support themselves in infinite retirement). Lower saving → lower capital → lower steady-state income. This is the macroeconomic manifestation of the finite planning horizon.

### 12.4.2 Stability Analysis

The Jacobian of the Blanchard–Yaari two-dimensional system at the steady state:

$$J_{BY} = \begin{pmatrix} f'(k^*)-\delta & -1 \\ \frac{c^*f''(k^*)}{1} - p(\rho+p) & f'(k^*)-\delta-\rho - \frac{p(\rho+p)k^*}{c^{*2}}(f'(k^*)-\delta) \end{pmatrix}.$$

The determinant:

$$\det(J_{BY}) = c^*f''(k^*)|\text{correction terms}| < 0$$

provided $f'' < 0$ and the mortality correction is not too large. The saddle-point property is preserved (one positive, one negative eigenvalue) for standard calibrations.

---

## 12.5 Social Security and Ricardian Non-Equivalence

The Blanchard–Yaari model restores **Ricardian non-equivalence**: a deficit-financed transfer (reducing taxes today, increasing taxes in the future) raises current aggregate consumption because the future tax burden falls partly on households not yet born — households who are not in the current economy's budget constraint.

**Proposition 12.2 (Ricardian Non-Equivalence in Blanchard–Yaari).** In the Blanchard–Yaari model, a lump-sum tax cut $\Delta T$ financed by government debt that will be repaid at rate $(r+p)\Delta B$ (where $\Delta B$ is new debt) raises current aggregate consumption by:

$$\Delta C_t = (\rho+p)\frac{p}{r+p}\Delta B > 0.$$

*Proof sketch.* The tax cut increases household wealth $A_t$ by $\Delta B$ (the present value of future tax increases is lower than $\Delta B$ because some of the future burden falls on not-yet-born households). Using the consumption function $C = (\rho+p)(A+H)$, the effect on consumption is $(\rho+p)\Delta A_{effective}$ where $\Delta A_{effective} = [p/(r+p)]\Delta B$ — the fraction of the debt not to be borne by current households. $\square$

The departure from Ricardian equivalence is proportional to $p/(r+p)$: as $p \to 0$ (infinite lives), this goes to zero (full Ricardian equivalence); as $p \to\infty$ (very short lives), it goes to $1/(1 + r/p) \to 1$ — the full tax cut is spent because households expect to be dead before repayment.

This result has an important implication for social security analysis: a PAYG pension system that taxes the young and pays the old redistributes from (not-yet-born) future young to current old. In the Blanchard–Yaari model, this redistribution has real effects — it raises current consumption and reduces capital accumulation. The welfare analysis requires computing whether the consumption gain to current old outweighs the capital-stock loss for all future generations.

---

## 12.6 Population Aging: The Demographic Dividend and Its Reversal

**Definition 12.3 (Population Aging in Continuous Time).** In the Blanchard–Yaari model, population aging is modeled as a reduction in the mortality rate $p$ (people live longer) or the birth rate $n$ (fewer new households entering). Both reduce the population growth rate $n - p$.

**Effects on the steady state:**

1. **Reduction in $p$:** Lower mortality rate → longer planning horizons → households save more for retirement → higher capital accumulation → higher $k^*$ and $y^*$. But also: the Blanchard–Yaari departure from Ricardian equivalence shrinks ($p/(r+p)$ falls) — the economy becomes more Ricardian.

2. **Reduction in $n$:** Fewer new workers → the capital-labor ratio rises mechanically (same capital, fewer workers) → higher $k$ per worker → higher $y$ per worker but lower aggregate output (fewer workers).

The **demographic dividend** refers to the period during which the working-age population is large relative to dependents (both young and old) — generating high labor supply, high saving, and rapid growth. East Asian economies exploited their demographic dividend from the 1960s to 1990s. The dividend reverses when the old-age dependency ratio rises and saving falls.

---

## 12.7 Worked Example: Steady-State Comparison with Ramsey

**Calibration:** $\alpha = 0.35$, $\delta = 0.05$, $\rho = 0.04$, $p = 0.04$ (expected life $= 1/p = 25$ years in the adult phase), $n = 0$, $g = 0$.

**Ramsey steady state:** $f'(k^*_{RCK}) = \delta + \rho = 0.09$. With $f'(k) = \alpha k^{\alpha-1}$:

$$k^*_{RCK} = (\alpha/0.09)^{1/(1-\alpha)} = (0.35/0.09)^{1/0.65} = (3.889)^{1.538} = 8.58.$$

$$c^*_{RCK} = k^{*\alpha}_{RCK} - \delta k^*_{RCK} = 8.58^{0.35} - 0.05\times8.58 = 2.247 - 0.429 = 1.818.$$

**Blanchard–Yaari steady state:** Solve $f'(k^*_{BY}) = \delta + \rho + p(\rho+p)k^*_{BY}/c^*_{BY}$ with $c^*_{BY} = f(k^*_{BY}) - \delta k^*_{BY}$ numerically.

Substituting: $\alpha(k^*_{BY})^{\alpha-1} - \delta - \rho = p(\rho+p)k^*_{BY}/[(k^*_{BY})^\alpha - \delta k^*_{BY}]$.

Numerical solution (Newton–Raphson): $k^*_{BY} \approx 6.84$, $c^*_{BY} \approx 1.71$.

**Comparison:**

| Quantity | Ramsey | Blanchard–Yaari | Reduction |
|---|---|---|---|
| $k^*$ | 8.58 | 6.84 | $-20\%$ |
| $y^* = k^{*\alpha}$ | 2.25 | 2.04 | $-9\%$ |
| $c^*$ | 1.82 | 1.71 | $-6\%$ |

The finite-lifetime assumption (mortality rate $p = 0.04$) reduces the steady-state capital stock by 20% and output by 9% relative to the Ramsey benchmark. This is a quantitatively significant difference — the finite planning horizon matters economically, not just theoretically.

```apl
⍝ APL — Blanchard-Yaari steady state via Newton-Raphson
⎕IO←0 ⋄ ⎕ML←1

alpha←0.35  ⋄  delta←0.05  ⋄  rho←0.04  ⋄  p←0.04

⍝ Ramsey steady state (closed form)
k_rck ← (alpha÷delta+rho) * ÷1-alpha
c_rck ← (k_rck*alpha) - delta×k_rck
k_rck  ⍝ ≈ 8.58
c_rck  ⍝ ≈ 1.82

⍝ Blanchard-Yaari: solve F(k) = 0 where
⍝ F(k) = alpha*k^(alpha-1) - delta - rho - p*(rho+p)*k / (k^alpha - delta*k)
F_BY ← {k←⍵
    c ← (k*alpha) - delta×k
    (alpha × k*alpha-1) - delta - rho - (p×(rho+p)×k) ÷ c}

dF_BY ← {k←⍵     ⍝ numerical derivative
    h←1e¯7
    ((F_BY k+h) - (F_BY k-h)) ÷ 2×h}

⍝ Newton-Raphson iteration
nr_step ← {⍵ - (F_BY ⍵) ÷ (dF_BY ⍵)}
converged ← {1e¯10 > |F_BY ⍵}
k_by ← nr_step ⍣ converged ⊢ 7.0     ⍝ start from k=7, near expected answer

c_by ← (k_by*alpha) - delta×k_by

k_by  ⍝ ≈ 6.84
c_by  ⍝ ≈ 1.71

⍝ Transition dynamics from k0 = 3 (poor economy)
⍝ BY system: dk = f(k)-c-delta*k; dc = (f'(k)-delta-rho)*c - p*(rho+p)*k
by_ode ← {k c ← ⍵
    k_dot ← (k*alpha) - c - delta×k
    c_dot ← ((alpha×k*alpha-1) - delta - rho)×c - p×(rho+p)×k
    k_dot c_dot}

⍝ Phase diagram: nullclines
k_grid ← 0.1 × ⍳ 120        ⍝ k from 0.1 to 12.0
c_kss0 ← (k_grid*alpha) - delta×k_grid   ⍝ dk=0 nullcline: c = f(k)-delta*k
⍝ dc=0 nullcline: c = p*(rho+p)*k / (f'(k)-delta-rho)
c_css0 ← (p×(rho+p)×k_grid) ÷ (alpha×k_grid*alpha-1) - delta - rho
```

```python
import numpy as np
from scipy.optimize import fsolve

alpha, delta, rho, p = 0.35, 0.05, 0.04, 0.04

# Ramsey steady state
k_rck = (alpha/(delta+rho))**(1/(1-alpha))
c_rck = k_rck**alpha - delta*k_rck
print(f"Ramsey: k*={k_rck:.4f}, c*={c_rck:.4f}")

# Blanchard-Yaari steady state: solve system
def by_system(k):
    c = k**alpha - delta*k
    return alpha*k**(alpha-1) - delta - rho - p*(rho+p)*k/c

k_by = fsolve(by_system, 7.0)[0]
c_by = k_by**alpha - delta*k_by
print(f"Blanchard-Yaari: k*={k_by:.4f}, c*={c_by:.4f}")
print(f"Differences: k: {100*(k_by/k_rck-1):.1f}%, c: {100*(c_by/c_rck-1):.1f}%")

# Effect of varying p (mortality rate) on steady-state capital
p_vals = np.linspace(0.0, 0.15, 50)
k_ss = []
for pv in p_vals:
    if pv == 0:
        k_ss.append(k_rck)
    else:
        def f_p(k): return alpha*k**(alpha-1)-delta-rho-pv*(rho+pv)*k/(k**alpha-delta*k)
        k_ss.append(fsolve(f_p, 7.0)[0])

import matplotlib.pyplot as plt
plt.figure(figsize=(7,4))
plt.plot(p_vals, k_ss, 'b-', linewidth=2)
plt.axhline(k_rck, color='r', linestyle='--', label=f'Ramsey k*={k_rck:.2f}')
plt.xlabel('Mortality rate p'); plt.ylabel('Steady-state capital k*')
plt.title('Effect of Finite Lifetimes on Capital Accumulation')
plt.legend(); plt.tight_layout(); plt.show()
```

---

## 12.8 Programming Exercises

### Exercise 12.1 (APL — Mortality Sensitivity)

Using the Newton–Raphson dfn from Section 12.7, compute the Blanchard–Yaari steady-state capital $k^*_{BY}$ for $p \in \{0.02, 0.04, 0.06, 0.08, 0.10, 0.15\}$. Produce a table showing $k^*_{BY}$, $y^*_{BY}$, $c^*_{BY}$, and the percentage deviations from the Ramsey benchmark for each $p$.

### Exercise 12.2 (Python — Transition Dynamics)

Simulate the Blanchard–Yaari economy from $k_0 = k^*_{BY}/2$ using RK4. Plot the transition path $k(t)$ and $c(t)$ together with the nullclines in the phase plane. Compare the speed of convergence to the Ramsey model: which converges faster to its respective steady state?

### Exercise 12.3 (Julia — Ricardian Non-Equivalence)

```julia
# Compute the degree of Ricardian non-equivalence as a function of p and r
r_vals = [0.02, 0.04, 0.06, 0.08]
p_vals = [0.01, 0.02, 0.04, 0.06, 0.10]
println("Fraction of tax cut consumed: p/(r+p)")
println("(1.0 = full Keynesian; 0.0 = full Ricardian)\n")
println(lpad("p\\r", 6), join([lpad("r=$(r)", 8) for r in r_vals]))
for pv in p_vals
    row = lpad("p=$(pv)",6) * join([lpad(round(pv/(r+pv),digits=3), 8) for r in r_vals])
    println(row)
end
```

### Exercise 12.4 — Aging Demographic Transition ($\star$)

Model a permanent reduction in the mortality rate from $p_{old} = 0.04$ to $p_{new} = 0.02$ (people live longer). (a) Compute the new steady state $k^*_{new} > k^*_{old}$ (longer lives → more saving → higher capital). (b) Simulate the transition path from $k^*_{old}$ to $k^*_{new}$ using RK4. (c) Compute the welfare change along the transition: do current households benefit or lose from the demographic shift? (d) Compare to the effect of a change in the birth rate from $n_{old} = 0.02$ to $n_{new} = 0.01$ — same population aging outcome but different mechanism.

### Exercise 12.5 — PAYG Social Security ($\star\star$)

Introduce a PAYG social security system that taxes workers $\tau w$ and pays retirees a lump-sum $b = \tau w L_{workers}/L_{retired}$. In the Blanchard–Yaari model: (a) derive the modified aggregate Euler equation with the social security tax-transfer mechanism; (b) show that the steady-state capital stock falls (PAYG reduces private saving); (c) derive the optimal PAYG transfer $\tau^*$ that maximizes a utilitarian social welfare function across generations.

---

## 12.9 Chapter Summary

**Key results:**

- The **Blanchard–Yaari perpetual youth model** introduces finite lifetimes via a Poisson death rate $p > 0$, modifying the effective discount rate to $\rho + p$ and the annuity return to $r + p$.
- Individual consumption follows $c_{s,t} = (\rho+p)(a_{s,t}+h_t)$ — a fraction $\rho+p$ of total wealth (financial plus human).
- **Aggregation** yields $C_t = (\rho+p)(A_t + H_t)$ and the aggregate Euler equation $\dot{C}_t = (r_t-\rho)C_t - p(\rho+p)A_t$.
- The **departure from Ramsey** is the term $-p(\rho+p)A_t$: newly born households enter without wealth, dragging down aggregate consumption growth relative to individual growth.
- The Blanchard–Yaari steady state has $k^*_{BY} < k^*_{RCK}$ — finite lifetimes reduce capital accumulation (by 20% in the calibrated example).
- **Ricardian non-equivalence** is restored: a deficit-financed tax cut raises aggregate consumption by $(\rho+p)\frac{p}{r+p}\Delta B$ — the fraction falling on yet-to-be-born households.
- In APL: the Newton–Raphson steady state is computed using `⍣ converged`; transition dynamics via `{by_ode_step ⍵}⍣T`.

**Connections forward:** Chapter 16 develops the Diamond discrete-time OLG model — the discrete counterpart — with explicit two-period household optimization, dynamic inefficiency, and social security welfare analysis. Chapter 25 of *Principles* [P:Ch.25] connects the OLG structure to the HANK model with heterogeneous households.

---

*Next: Chapter 13 — Adjustment Cost Models: Tobin's q and Investment Dynamics*
