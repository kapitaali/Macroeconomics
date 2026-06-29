# Part III: Continuous-Time Dynamic Models in Macroeconomics

*Connects to: Principles Ch. 5 (Solow, RCK), Ch. 12 (investment), Ch. 25 (OLG)*

---

Continuous time is the natural habitat of growth theory and optimal control. The differential equation $\dot{\tilde{k}} = sf(\tilde{k}) - \mu\tilde{k}$ is a single, clean object with a well-developed analytical theory; the phase plane of the RCK model reveals the saddle-path geometry of optimal saving at a glance; and the Hamiltonian framework of Pontryagin unifies all intertemporal optimization problems into a single set of necessary conditions.

This part derives the full mathematical solutions to the growth models of *Principles* Part I, using the ODE and phase-diagram tools of Chapter 3. Chapter 3 was preparation; this part is the application. Every model here was introduced verbally and graphically in *Principles*; here we solve it.

**Chapter 10** solves the Solow growth model as a Bernoulli ODE, derives the closed-form transition path, proves existence and uniqueness of the steady state from the Inada conditions, and derives the Golden Rule algebraically. **Chapter 11** develops the full apparatus of optimal control theory — Hamiltonians, Pontryagin's maximum principle, costate equations, the transversality condition — and applies it to the Ramsey–Cass–Koopmans model to derive the Euler equation and saddle-path solution. **Chapter 12** analyzes the continuous-time Blanchard–Yaari overlapping-generations model, the tractable alternative to the Diamond (1965) discrete-time OLG. **Chapter 13** returns to the Tobin's $q$ investment model, deriving the investment dynamics from an optimal control problem and using real options theory to analyze irreversibility.

The level of mathematical detail in this part is high. Every proof is complete; every algorithm is stated in language-neutral pseudocode before any code is shown. Readers who want only the results without the derivations can follow the boxed equations and worked examples; readers who want to understand *why* each result holds should work through each proof step by step.

---

# Chapter 10: The Solow Growth Model

*Solving the Fundamental Differential Equation*

> *"The Solow model is the starting point for almost all analyses of economic growth. It is simple, tractable, and delivers surprisingly rich insights."*
> — Robert Barro and Xavier Sala-i-Martin

**Cross-reference:** *Principles* Ch. 5.2 (Solow–Swan model and the fundamental ODE); Ch. 5.5 (convergence, conditional and unconditional) **[P:Ch.5.2, P:Ch.5.5]**

---

## 10.1 Setting Up the Fundamental ODE from First Principles

*Principles* Chapter 5 presented the Solow ODE as a given. Here we derive it from the underlying production and saving assumptions, making explicit every step and every assumption.

**The model's building blocks:**

- A single good is produced by capital $K$ and effective labour $AL$ using a neoclassical production function $Y = F(K, AL)$.
- Population grows at rate $n$: $\dot{L}/L = n$.
- Technology grows at rate $g$: $\dot{A}/A = g$.
- A constant fraction $s \in (0,1)$ of output is saved: $\dot{K} = sY - \delta K$.
- Capital depreciates at rate $\delta > 0$.

**Definition 10.1 (Intensive Form).** Define capital per unit of effective labour $\tilde{k} \equiv K/(AL)$ and output per unit of effective labour $\tilde{y} \equiv Y/(AL) = f(\tilde{k})$, where $f(\tilde{k}) \equiv F(\tilde{k}, 1)$ is the intensive-form production function.

**Deriving the ODE.** Differentiate $\tilde{k} = K/(AL)$ with respect to $t$:

$$\dot{\tilde{k}} = \frac{\dot{K}}{AL} - \frac{K}{AL}\cdot\frac{\dot{A}L + A\dot{L}}{AL} = \frac{\dot{K}}{AL} - \tilde{k}\left(\frac{\dot{A}}{A} + \frac{\dot{L}}{L}\right) = \frac{sY - \delta K}{AL} - \tilde{k}(g+n).$$

Substituting $Y = ALf(\tilde{k})$ and $\delta K = \delta\tilde{k}AL$:

$$\dot{\tilde{k}} = sf(\tilde{k}) - \delta\tilde{k} - (n+g)\tilde{k} = sf(\tilde{k}) - \underbrace{(n+g+\delta)}_{\equiv\,\mu}\tilde{k}.$$

**Definition 10.2 (Fundamental Solow ODE).** The **fundamental equation of the Solow model** is:

$$\boxed{\dot{\tilde{k}} = sf(\tilde{k}) - \mu\tilde{k}, \quad \mu \equiv n + g + \delta.}$$

Every element of this equation has a clear interpretation: $sf(\tilde{k})$ is actual investment per effective worker; $\mu\tilde{k} = (n+g+\delta)\tilde{k}$ is the **break-even investment** — the investment required to keep $\tilde{k}$ constant as population grows ($n$), technology improves ($g$), and capital depreciates ($\delta$).

---

## 10.2 The Cobb–Douglas Special Case: Closed-Form Solution

With $f(\tilde{k}) = \tilde{k}^\alpha$, $\alpha \in (0,1)$, the Solow ODE becomes:

$$\dot{\tilde{k}} = s\tilde{k}^\alpha - \mu\tilde{k}.$$

### 10.2.1 The Steady State

**Definition 10.3 (Steady State).** The **steady state** $\tilde{k}^*$ satisfies $\dot{\tilde{k}} = 0$:

$$s\tilde{k}^{*\alpha} = \mu\tilde{k}^* \implies \tilde{k}^{*\alpha-1} = \frac{\mu}{s} \implies \tilde{k}^* = \left(\frac{s}{\mu}\right)^{1/(1-\alpha)}.$$

Steady-state output and consumption per effective worker:

$$\tilde{y}^* = \tilde{k}^{*\alpha} = \left(\frac{s}{\mu}\right)^{\alpha/(1-\alpha)}, \qquad \tilde{c}^* = (1-s)\tilde{y}^* = (1-s)\left(\frac{s}{\mu}\right)^{\alpha/(1-\alpha)}.$$

### 10.2.2 The Bernoulli ODE: Closed-Form Transition Path

The Cobb–Douglas Solow ODE is a **Bernoulli equation** — a nonlinear first-order ODE of the form $\dot{z} + P(t)z = Q(t)z^\alpha$. The substitution $v = \tilde{k}^{1-\alpha}$ linearizes it.

Let $v \equiv \tilde{k}^{1-\alpha}$. Differentiating:

$$\dot{v} = (1-\alpha)\tilde{k}^{-\alpha}\dot{\tilde{k}} = (1-\alpha)\tilde{k}^{-\alpha}(s\tilde{k}^\alpha - \mu\tilde{k}) = (1-\alpha)(s - \mu v).$$

This is a linear first-order ODE in $v$! It has the form $\dot{v} = -(1-\alpha)\mu v + (1-\alpha)s$, with steady state $v^* = s/\mu = \tilde{k}^{*1-\alpha}$ and convergence rate $\lambda = (1-\alpha)\mu$.

**Theorem 10.1 (Closed-Form Solow Transition Path).** The transition path of the Cobb–Douglas Solow model from initial condition $\tilde{k}_0 > 0$ is:

$$\tilde{k}(t) = \left[v^* + (v_0 - v^*)e^{-\lambda t}\right]^{1/(1-\alpha)},$$

where $v_0 = \tilde{k}_0^{1-\alpha}$, $v^* = s/\mu = \tilde{k}^{*1-\alpha}$, and $\lambda = (1-\alpha)\mu$.

*Proof.* The solution to $\dot{v} = -\lambda v + (1-\alpha)s$ with initial condition $v_0$ is:

$$v(t) = v^* + (v_0 - v^*)e^{-\lambda t}, \quad v^* = \frac{(1-\alpha)s}{\lambda} = \frac{(1-\alpha)s}{(1-\alpha)\mu} = \frac{s}{\mu}.$$

Converting back via $\tilde{k}(t) = v(t)^{1/(1-\alpha)}$ gives the result. $\square$

This is the **exact** transition path — not a linearization. It is valid for any initial condition $\tilde{k}_0$, not just those near $\tilde{k}^*$.

---

## 10.3 Linearization and the Convergence Rate

For empirical work — particularly convergence regressions — we need the **local** approximation of the transition path near the steady state.

**Theorem 10.2 (Linear Convergence Rate).** Near the steady state $\tilde{k}^*$ of the Cobb–Douglas Solow model, the log-deviation $\hat{k}(t) \equiv \ln\tilde{k}(t) - \ln\tilde{k}^*$ evolves approximately as:

$$\hat{k}(t) \approx \hat{k}_0 e^{-\lambda t}, \quad \lambda = (1-\alpha)\mu = (1-\alpha)(n+g+\delta).$$

*Proof.* From the exact solution, near $v_0 \approx v^*$:

$$\frac{v(t) - v^*}{v^*} = \frac{v_0 - v^*}{v^*}e^{-\lambda t}.$$

Since $v = \tilde{k}^{1-\alpha}$ and $\ln v \approx (1-\alpha)\ln\tilde{k}$ for $\tilde{k}$ near $\tilde{k}^*$, we have $\hat{v}(t) \approx (1-\alpha)\hat{k}(t)$, so $\hat{k}(t) \approx \hat{k}_0 e^{-\lambda t}$. $\square$

**Definition 10.4 (Half-Life of Convergence).** The **half-life** is the time for the gap to halve: $\tau_{1/2} = \ln 2 / \lambda$.

**Calibration:** With $\alpha = 1/3$, $n = 0.01$, $g = 0.02$, $\delta = 0.05$: $\mu = 0.08$, $\lambda = (2/3)(0.08) = 0.053$, $\tau_{1/2} = \ln 2/0.053 \approx 13$ years. This matches the empirical estimates of cross-country convergence rates of approximately 2–3% per year — which corresponds to $\lambda = 0.02$–$0.03$, somewhat below the model's prediction with these parameters. The discrepancy is partly explained by international capital flows and human capital omitted from the basic model.

---

## 10.4 Existence and Uniqueness: The Inada Conditions

For general production functions, the existence and uniqueness of $\tilde{k}^*$ requires the **Inada conditions**.

**Definition 10.5 (Inada Conditions).** The production function $f: \mathbb{R}_+ \to \mathbb{R}_+$ satisfies the **Inada conditions** if:

$$f(0) = 0, \quad f'(\tilde{k}) > 0, \quad f''(\tilde{k}) < 0, \quad \lim_{\tilde{k}\to 0}f'(\tilde{k}) = +\infty, \quad \lim_{\tilde{k}\to+\infty}f'(\tilde{k}) = 0.$$

**Theorem 10.3 (Existence and Uniqueness of the Steady State).** If $f$ satisfies the Inada conditions and $s \in (0,1)$, $\mu > 0$, then the equation $sf(\tilde{k}) = \mu\tilde{k}$ has a unique positive solution $\tilde{k}^* > 0$, and this steady state is globally asymptotically stable.

*Proof.* Define $g(\tilde{k}) \equiv sf(\tilde{k}) - \mu\tilde{k}$. The Solow ODE is $\dot{\tilde{k}} = g(\tilde{k})$.

**Existence:** $g(0) = sf(0) - 0 = 0$... but we need $g'(0) > 0$ and $g'(\infty) < 0$ to guarantee a crossing above zero. At $\tilde{k} = 0^+$: $sf'(0^+) = +\infty > \mu$, so $g(\tilde{k}) > 0$ for small $\tilde{k} > 0$. As $\tilde{k} \to \infty$: $f'(\tilde{k}) \to 0$, and by L'Hôpital-style argument $f(\tilde{k})/\tilde{k} \to 0$ (since $f$ is concave with $f' \to 0$), so $g(\tilde{k}) = \tilde{k}[sf(\tilde{k})/\tilde{k} - \mu] \to -\infty$. By the intermediate value theorem, $g$ crosses zero at some $\tilde{k}^* > 0$.

**Uniqueness:** $g'(\tilde{k}) = sf'(\tilde{k}) - \mu$. Since $f'' < 0$, $g'$ is strictly decreasing. At the crossing, $g'(\tilde{k}^*) = sf'(\tilde{k}^*) - \mu$. The sign of $g'(\tilde{k}^*)$ determines local stability: the crossing is from above (stable) iff $g'(\tilde{k}^*) < 0$, i.e., $sf'(\tilde{k}^*) < \mu$. By convexity of the $\mu\tilde{k}$ line and concavity of $sf(\tilde{k})$, they can cross only once, and at the crossing the slope of $sf$ is below the slope of $\mu\tilde{k}$. Hence $\tilde{k}^*$ is unique and stable. $\square$

---

## 10.5 The Golden Rule: Optimal Saving

**Definition 10.6 (Golden Rule Capital Stock).** The **Golden Rule capital stock** $\tilde{k}^{GR}$ maximizes steady-state consumption per effective worker $\tilde{c}^* = f(\tilde{k}^*) - \mu\tilde{k}^*$ over $s$ (equivalently over $\tilde{k}^*$, since $\tilde{k}^*$ is a monotone function of $s$):

$$\max_{\tilde{k}^*} \tilde{c}^* = f(\tilde{k}^*) - \mu\tilde{k}^* \implies f'(\tilde{k}^{GR}) = \mu = n + g + \delta.$$

**Theorem 10.4 (Golden Rule).** The Golden Rule capital stock satisfies $f'(\tilde{k}^{GR}) = \mu$, which corresponds to the saving rate:

$$s^{GR} = \mu\tilde{k}^{GR}/f(\tilde{k}^{GR}) = \alpha \quad \text{(for Cobb–Douglas)}.$$

*Proof (Cobb–Douglas):* With $f(\tilde{k}) = \tilde{k}^\alpha$, the Golden Rule condition $f'(\tilde{k}^{GR}) = \alpha\tilde{k}^{GR\,\alpha-1} = \mu$ gives $\tilde{k}^{GR} = (\alpha/\mu)^{1/(1-\alpha)}$. The corresponding steady state has $s^{GR}\tilde{k}^{GR\,\alpha} = \mu\tilde{k}^{GR}$, so $s^{GR} = \mu\tilde{k}^{GR\,1-\alpha} = \mu(\alpha/\mu) = \alpha$. $\square$

The Golden Rule saving rate equals the capital share $\alpha$ for Cobb–Douglas production. Countries with $s > s^{GR} = \alpha$ save too much (dynamically inefficient); countries with $s < \alpha$ save too little. For the U.S.: $\alpha \approx 0.33$ and $s \approx 0.20 < \alpha$, suggesting the U.S. is below the Golden Rule — a common finding for advanced economies [P:Ch.5.3].

**Dynamic efficiency:** The RCK model (Chapter 11) shows that optimizing households choose $s < s^{GR}$ in general, so market economies are dynamically efficient. The Golden Rule is relevant when policy interventions (mandatory saving, social security) push $s$ above $\alpha$.

---

## 10.6 Numerical Solution: The Shooting Algorithm

For production functions that do not yield a closed-form Bernoulli ODE solution, we solve numerically. The **shooting algorithm** finds the transition path from $\tilde{k}_0$ to (near) $\tilde{k}^*$ by integrating the ODE forward using RK4.

**Algorithm 10.1 (Solow Transition Path).**

Input: $f$, $s$, $\mu$, $\tilde{k}_0$, time horizon $T$, step size $h$.
1. Set $t = 0$, $\tilde{k}[0] = \tilde{k}_0$.
2. For $n = 0, 1, \ldots, T/h - 1$:
   - $F(\tilde{k}) \leftarrow sf(\tilde{k}) - \mu\tilde{k}$
   - Apply RK4 step: $\tilde{k}[n+1] \leftarrow \text{RK4}(F, \tilde{k}[n], h)$
   - $t \leftarrow t + h$
3. Return $\{\tilde{k}[n]\}$.

For the RCK model (Chapter 11), the shooting is more complex because the initial value of the **control variable** $c$ is not given — only the initial value of the **state variable** $k$ is. The shooting algorithm must find $c_0$ such that the trajectory starting at $(k_0, c_0)$ converges to the saddle-point steady state rather than diverging.

---

## 10.7 Worked Example: Mankiw–Romer–Weil Convergence Regression

*Cross-reference: Principles Ch. 5.5 (convergence evidence)* **[P:Ch.5.5]**

Mankiw, Romer, and Weil (1992) test the Solow model's conditional convergence prediction by regressing long-run average growth rates on initial income, saving rates, and population growth. The model predicts the log-linearized growth equation:

$$\ln y(T) - \ln y(0) = (1-e^{-\lambda T})\left[\frac{\alpha}{1-\alpha}\ln s - \frac{\alpha}{1-\alpha}\ln(n+g+\delta) - \ln y(0)\right],$$

where $y = Y/L$ is output per worker (not per effective worker). Setting $T = 25$ years, $g = 0.02$, $\delta = 0.03$ (so $g+\delta = 0.05$), and estimating $\lambda \approx 0.02$ from the data (implying $\alpha \approx 0.60$, suggesting physical plus human capital together):

$$\ln y(T) - \ln y(0) = 0.393\left[\frac{0.60}{0.40}\ln s - \frac{0.60}{0.40}\ln(n+0.05) - \ln y(0)\right].$$

The predicted coefficients: on $\ln s$ approximately $+0.59$; on $\ln(n+g+\delta)$ approximately $-0.59$; on $\ln y(0)$ approximately $-0.39$. MRW report empirical estimates very close to these, providing strong support for the augmented Solow model.

```apl
⍝ APL — Solow model: closed-form path and MRW regression simulation
⎕IO←0 ⋄ ⎕ML←1

⍝ Parameters
alpha ← 1÷3   ⋄   s ← 0.20   ⋄   mu ← 0.08

⍝ Steady state
kstar ← (s÷mu)*÷1-alpha
ystar ← kstar*alpha
cstar ← (1-s)×ystar
kstar ystar cstar    ⍝ ≈ 3.95 1.58 1.26

⍝ Convergence rate and half-life
lambda ← (1-alpha)×mu    ⍝ ≈ 0.053
halflife ← (⍟2)÷lambda   ⍝ ≈ 13 years

⍝ Closed-form transition path: k(t) for t=0..100
vstar ← kstar*1-alpha
k0    ← 0.5              ⍝ initial capital (below steady state)
v0    ← k0*1-alpha
T     ← 100

t ← ⍳T
v_path ← vstar + (v0-vstar)×(*-lambda×t)   ⍝ v(t) = v* + (v0-v*)e^{-λt}
k_path ← v_path*÷1-alpha                   ⍝ k(t) = v(t)^{1/(1-α)}

⍝ Gap closes by half every halflife years
gap_t ← |k_path - kstar
gap_t[0] gap_t[13] gap_t[26]   ⍝ initial, 1 halflife, 2 halflives

⍝ MRW regression: simulate 100 countries
N ← 100
⎕RL ← 42                         ⍝ set random seed (Dyalog)
s_i   ← 0.05 + 0.35 × ? N⍴0    ⍝ saving rates U[0.05, 0.40]
n_i   ← 0.00 + 0.03 × ? N⍴0    ⍝ population growth U[0, 0.03]
mu_i  ← n_i + 0.05               ⍝ n+g+delta, g+delta=0.05
k0_i  ← 0.5 + 2 × ? N⍴0        ⍝ initial capital variation

⍝ Steady state for each country
kstar_i ← (s_i÷mu_i)*÷1-alpha
lny0_i  ← ⍟(k0_i*alpha)          ⍝ log initial output per worker
lnyT_i  ← ⍟(kstar_i*alpha)       ⍝ log steady-state output (approximate long-run)
growth_i ← lnyT_i - lny0_i       ⍝ 25-year growth

⍝ MRW regressors: ln(s), ln(n+g+delta), ln(y0)
X ← ↑(⍟s_i)(⍟mu_i)(lny0_i)(N⍴1)   ⍝ 4×N design matrix
⍝ OLS: b = (X'X)^-1 X'y
Xmat ← ⍉X         ⍝ N×4
b_ols ← (⍉Xmat ⌹ Xmat) +.× Xmat +.× growth_i
b_ols   ⍝ coefficients: should be ≈ (+, -, -, +intercept)
```

```python
import numpy as np
import matplotlib.pyplot as plt

# Parameters
alpha, s, mu = 1/3, 0.20, 0.08
lambda_ = (1-alpha)*mu
kstar = (s/mu)**(1/(1-alpha))
vstar = kstar**(1-alpha)

# Closed-form transition path
k0 = 0.5
T = 100
t = np.arange(T)
v0 = k0**(1-alpha)
v_path = vstar + (v0 - vstar)*np.exp(-lambda_*t)
k_path = v_path**(1/(1-alpha))

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
ax.plot(t, k_path, 'b-', linewidth=2, label='k̃(t)')
ax.axhline(kstar, color='red', linestyle='--', label=f'k̃* = {kstar:.2f}')
ax.set_xlabel('Time (years)'); ax.set_ylabel('k̃')
ax.set_title('Solow Transition Path (Closed Form)')
ax.legend()

# MRW simulation
np.random.seed(42)
N = 100
s_i   = np.random.uniform(0.05, 0.40, N)
n_i   = np.random.uniform(0.00, 0.03, N)
mu_i  = n_i + 0.05
kstar_i = (s_i/mu_i)**(1/(1-alpha))
lns   = np.log(s_i); lnmu = np.log(mu_i)

# Simulate 25-year growth: use linearized convergence
lny0  = np.random.normal(np.log(kstar_i**alpha), 0.3)   # random initial incomes
conv  = 1 - np.exp(-lambda_*25)
lnyT  = lny0 + conv*(np.log(kstar_i**alpha) - lny0)
growth_25 = lnyT - lny0

# OLS regression: growth on ln(s), ln(n+g+d), ln(y0)
X = np.column_stack([lns, lnmu, lny0, np.ones(N)])
b = np.linalg.lstsq(X, growth_25, rcond=None)[0]
print("MRW regression coefficients:")
print(f"  ln(s):         {b[0]:.3f}  (predicted: +{alpha/(1-alpha)*conv:.3f})")
print(f"  ln(n+g+δ):     {b[1]:.3f}  (predicted: {-alpha/(1-alpha)*conv:.3f})")
print(f"  ln(y₀):        {b[2]:.3f}  (predicted: {-conv:.3f})")

ax2 = axes[1]
ax2.scatter(lny0, growth_25, alpha=0.5, s=20)
ax2.set_xlabel('ln(y₀)'); ax2.set_ylabel('25-year growth')
ax2.set_title('Conditional convergence (MRW simulation)')
plt.tight_layout(); plt.show()
```

```julia
using LinearAlgebra, Statistics, Random

Random.seed!(42)
alpha, mu = 1/3, 0.08
lambda_ = (1-alpha)*mu

# Closed-form path
kstar(s) = (s/mu)^(1/(1-alpha))
vstar(s) = kstar(s)^(1-alpha)

function solow_path(k0, s, T)
    v0 = k0^(1-alpha); vs = vstar(s)
    t = 0:T-1
    v_t = vs .+ (v0-vs).*exp.(-lambda_.*t)
    return v_t.^(1/(1-alpha))
end

k_path = solow_path(0.5, 0.20, 100)
println("k̃* = $(round(kstar(0.20), digits=3))")
println("k̃(50) = $(round(k_path[51], digits=3))  (should be ≈ k̃*)")

# MRW simulation
N = 100
s_i  = 0.05 .+ 0.35*rand(N)
n_i  = 0.00 .+ 0.03*rand(N)
mu_i = n_i .+ 0.05
ks_i = kstar.(s_i ./ mu_i .* mu_i)   # simplify: kstar per country
lny0 = log.(kstar.(s_i).*alpha) .+ 0.3*randn(N)
conv = 1 - exp(-lambda_*25)
lnyT = lny0 .+ conv.*(log.(kstar.(s_i).^alpha) .- lny0)
growth = lnyT .- lny0

X = hcat(log.(s_i), log.(mu_i), lny0, ones(N))
b = X \ growth
println("\nMRW OLS: coefficients = ", round.(b, digits=3))
```

```r
alpha <- 1/3; mu <- 0.08; lambda_ <- (1-alpha)*mu

# Closed-form path
kstar <- function(s) (s/mu)^(1/(1-alpha))
solow_path <- function(k0, s, T) {
  v0 <- k0^(1-alpha); vs <- kstar(s)^(1-alpha)
  t <- 0:(T-1)
  v_t <- vs + (v0-vs)*exp(-lambda_*t)
  v_t^(1/(1-alpha))
}

k_path <- solow_path(0.5, 0.20, 100)
cat(sprintf("k̃* = %.3f,  k̃(50) = %.3f\n", kstar(0.20), k_path[51]))

# MRW regression simulation
set.seed(42); N <- 100
s_i  <- runif(N, 0.05, 0.40); n_i <- runif(N, 0, 0.03); mu_i <- n_i + 0.05
lny0 <- log(kstar(s_i)^alpha) + rnorm(N, 0, 0.3)
conv <- 1 - exp(-lambda_*25)
lnyT <- lny0 + conv*(log(kstar(s_i)^alpha) - lny0)
growth <- lnyT - lny0

df <- data.frame(growth, lns=log(s_i), lnmu=log(mu_i), lny0=lny0)
fit <- lm(growth ~ lns + lnmu + lny0, data=df)
cat("\nMRW OLS coefficients:\n"); print(round(coef(fit), 3))
```

---

## 10.8 Programming Exercises

### Exercise 10.1 (APL — Phase Diagram)

Write a dfn `solow_phase ← {s mu alpha ← ⍵ ⋄ ...}` that takes parameters and returns: (a) the steady state $\tilde{k}^*$, convergence rate $\lambda$, and half-life $\tau_{1/2}$; (b) vectors `k_grid`, `invest`, `breakeven` of length 100 over $\tilde{k} \in [0.01, 3\tilde{k}^*]$; (c) the crossing point. Plot using `]Plot` or export to CSV. Verify the crossing is at $\tilde{k}^*$.

### Exercise 10.2 (Python — Parameter Sensitivity)

For $\alpha \in \{0.2, 0.33, 0.5\}$ and $s \in \{0.10, 0.15, \ldots, 0.40\}$, compute: (a) the steady-state capital $\tilde{k}^*$ and output $\tilde{y}^*$; (b) the Golden Rule saving rate $s^{GR} = \alpha$; (c) the fraction $(s - s^{GR})/s^{GR}$. Plot a $3\times6$ table of the dynamic efficiency gap for each $(\alpha, s)$ combination.

### Exercise 10.3 (Julia — Non–Cobb–Douglas)

Implement the RK4 solver for the Solow model with the CES production function $f(\tilde{k}) = [\gamma\tilde{k}^\rho + (1-\gamma)]^{1/\rho}$ (which nests Cobb–Douglas as $\rho \to 0$). For $\gamma = 0.4$, $\rho \in \{-1, -0.5, -0.1, 0.1, 0.5, 1\}$, $s = 0.25$, $\mu = 0.08$: (a) find $\tilde{k}^*$ numerically using Newton–Raphson; (b) compare the convergence rate to the Cobb–Douglas approximation. Does higher substitution elasticity ($\rho \to 1$) speed or slow convergence?

### Exercise 10.4 (R — The Augmented Solow Model)

Mankiw, Romer, and Weil (1992) augment the Solow model with human capital $H$: $Y = K^\alpha H^\beta (AL)^{1-\alpha-\beta}$. The fundamental equations become a 2D ODE system in $(\tilde{k}, \tilde{h})$. (a) Derive the two ODEs. (b) Show the steady state satisfies $\tilde{k}^* = (s_k^{1-\beta}s_h^\beta/\mu)^{1/(1-\alpha-\beta)}$ and similarly for $\tilde{h}^*$. (c) Simulate 100 countries with $s_k \in U[0.05, 0.30]$, $s_h \in U[0.02, 0.15]$, $n \in U[0, 0.03]$, $\alpha = 1/3$, $\beta = 1/3$, and run the MRW regression. Do the coefficients match the theoretical predictions?

### Exercise 10.5 — The AK Model Failure ($\star$)

For the AK production function $f(\tilde{k}) = A\tilde{k}$ (linear, no diminishing returns): (a) show the Solow ODE becomes $\dot{\tilde{k}} = (sA - \mu)\tilde{k}$, a linear ODE with solution $\tilde{k}(t) = \tilde{k}_0 e^{(sA-\mu)t}$; (b) there is no finite steady state — the economy grows forever at rate $g^* = sA - \mu$ if $sA > \mu$; (c) the Inada conditions fail (the lower Inada condition $f'(0) = +\infty$ is violated since $f'(\tilde{k}) = A$ is constant); (d) explain what this implies for the convergence result: show that cross-country income differences are permanent in the AK model but transitory in the Solow model.

### Exercise 10.6 — Dynamic Efficiency Test ($\star\star$)

Abel, Mankiw, Summers, and Zeckhauser (1989) test dynamic efficiency by comparing aggregate profits $\Pi = Y - wL$ to aggregate investment $I = \dot{K} + \delta K$ across countries. Show that the economy is dynamically efficient iff $\Pi > I$ (net factor payments to capital exceed net investment). (a) Express this condition in terms of the capital share $\alpha$, the investment-to-output ratio $I/Y = s$, and the capital-output ratio $K/Y$. (b) For the U.S.: $\alpha \approx 0.33$, $s \approx 0.20$, $K/Y \approx 2.5$; evaluate the AMSZ test. (c) Simulate 50 OECD countries (with calibrated parameters from the Penn World Tables documentation) and identify which are dynamically inefficient.

---

## 10.9 Chapter Summary

**Key results:**

- The **Solow ODE** $\dot{\tilde{k}} = sf(\tilde{k}) - \mu\tilde{k}$ is derived by differentiating $\tilde{k} = K/(AL)$ and substituting capital accumulation and factor growth.
- For **Cobb–Douglas** $f(\tilde{k}) = \tilde{k}^\alpha$: the substitution $v = \tilde{k}^{1-\alpha}$ transforms the nonlinear ODE into a linear one, yielding the closed-form path $\tilde{k}(t) = [v^* + (v_0-v^*)e^{-\lambda t}]^{1/(1-\alpha)}$ with $\lambda = (1-\alpha)\mu$.
- The **convergence rate** $\lambda = (1-\alpha)\mu$ and half-life $\tau_{1/2} = \ln 2/\lambda \approx 13$ years (for standard calibration) match empirical cross-country convergence estimates.
- The **Inada conditions** guarantee existence and uniqueness of the steady state; a formal existence-uniqueness proof uses the intermediate value theorem and concavity.
- The **Golden Rule** saving rate equals $s^{GR} = \alpha$ for Cobb–Douglas — countries saving more than $\alpha$ are dynamically inefficient; most OECD countries save less than $\alpha$ and are efficient.
- In APL: the closed-form path is `v_path ← vstar + (v0-vstar) × * (-lambda) × ⍳T`; `k_path ← v_path * ÷ 1-alpha` — the entire transition path in two lines.

**Connections forward:** Chapter 11 shows why optimizing households (RCK) choose $s < s^{GR}$, avoiding dynamic inefficiency. Chapter 17 implements the stochastic RBC extension of the Solow model by adding TFP shocks to the capital accumulation equation and solving via value function iteration.

---

*Next: Chapter 11 — The Ramsey–Cass–Koopmans Model: Optimal Control and the Euler Equation*
