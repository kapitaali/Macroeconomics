# Chapter 35: Integrated Assessment Models

*Climate-Economy Modeling, DICE, and the Social Cost of Carbon*

> *"The social cost of carbon is the most important price in the world. It should reflect the present discounted value of all future damages from one additional ton of CO₂."*
> — Nicholas Stern

**Cross-reference:** *Principles* Ch. 37 (climate change and macroeconomics, DICE model, social cost of carbon, Stern vs. Nordhaus); Ch. 40 (capstone: climate policy) **[P:Ch.37, P:Ch.40]**

---

## 35.1 The Integrated Assessment Modeling Framework

An **Integrated Assessment Model (IAM)** combines an economic model of production, consumption, and investment with a physical climate model of emissions, atmospheric concentrations, and temperature — solving for the socially optimal emission path and the corresponding **social cost of carbon (SCC)**.

The DICE model (Dynamic Integrated model of Climate and the Economy; Nordhaus, 1992, 2008, 2018) is the most widely used IAM. It is the quantitative underpinning of U.S. federal regulatory policy — the SCC from DICE variants is used by the EPA to value carbon regulations.

Mathematically, DICE is a **continuous-time optimal control problem**: the social planner chooses the abatement fraction $\mu(t)$ (emissions reduction relative to business-as-usual) to maximize global welfare subject to economic and climate dynamics. The optimal abatement path and the SCC emerge from the Hamiltonian conditions of Chapter 11.

---

## 35.2 The DICE Model Structure

### 35.2.1 Economic Module

The economy produces output $Y(t)$ using capital $K(t)$ and exogenous labor $L(t)$:

$$Y(t) = \Omega(T_{AT})\cdot[1-\Lambda(\mu(t))]\cdot A(t)K(t)^\gamma(A_L(t)L(t))^{1-\gamma},$$

where:
- $\Omega(T_{AT}) = 1/(1 + \pi_1 T_{AT} + \pi_2 T_{AT}^2)$ is the **damage function** — output fraction lost due to temperature increase $T_{AT}$ above pre-industrial. Nordhaus calibrates $\pi_1 \approx 0$, $\pi_2 \approx 0.00236$.
- $\Lambda(\mu) = \theta_1\mu^{\theta_2}$ is the **abatement cost function** — output fraction spent on emission reduction.
- $A(t)$ is total factor productivity (exogenous).

**Capital accumulation:** $\dot{K} = Y(t) - C(t) - \delta_K K(t)$.

**Emissions:** $E(t) = \sigma(t)(1-\mu(t))Y(t)$, where $\sigma(t)$ is emission intensity (declining exogenously with energy efficiency).

**Welfare objective:**

$$\max_{\mu(\cdot),C(\cdot)}\int_0^\infty e^{-\rho t}L(t)\frac{[C(t)/L(t)]^{1-\sigma}-1}{1-\sigma}\,dt.$$

### 35.2.2 Climate Module

**Carbon cycle:** Three carbon reservoirs (atmosphere, upper ocean, lower ocean) exchange carbon with flow rates $\phi_{ij}$:

$$\dot{M}_{AT} = E(t) - \phi_{10}M_{AT} + \phi_{01}M_{UO}, \qquad \text{etc.}$$

This is a linear system $\dot{\mathbf{M}} = \Phi\mathbf{M} + \mathbf{e}(t)$ (where $\mathbf{e}$ contains $E(t)$ in the first entry).

**Temperature dynamics:**

$$C_{AT}\dot{T}_{AT} = F(M_{AT}) - \lambda T_{AT} - \frac{C_{AT}}{C_{OC}}(T_{AT} - T_{OC}),$$
$$C_{OC}\dot{T}_{OC} = \frac{C_{AT}}{C_{OC}}(T_{AT} - T_{OC}),$$

where $F(M_{AT}) = \eta\log_2(M_{AT}/M_{AT,1750})$ is the radiative forcing (logarithm of atmospheric carbon concentration ratio).

---

## 35.3 The Social Cost of Carbon as a Costate Variable

The SCC is derived formally from the Hamiltonian of the DICE optimal control problem. The key insight is that the SCC is the **shadow price** of atmospheric carbon — exactly the costate variable of the Hamiltonian, analogous to Tobin's $q$ in Chapter 13.

**Setting up the Hamiltonian:** The state variables are $(K, M_{AT}, T_{AT}, \ldots)$. The control variables are $(\mu, C)$. The current-value Hamiltonian:

$$\mathcal{H} = L\frac{(C/L)^{1-\sigma}-1}{1-\sigma} + q_K\dot{K} + q_{AT}\dot{M}_{AT} + q_{T}\dot{T}_{AT} + \ldots$$

where $q_K$ is the shadow price of capital (Tobin's $q$), $q_{AT}$ is the shadow price of atmospheric carbon (**the SCC**), and $q_T$ is the shadow price of temperature.

**Definition 35.1 (Social Cost of Carbon).** The **social cost of carbon** is the negative of the costate variable associated with atmospheric CO₂:

$$SCC(t) \equiv -\frac{\partial\mathcal{H}}{\partial E}\bigg/\frac{\partial\mathcal{H}}{\partial C} = -\frac{q_{AT}(t)}{q_K(t)}.$$

This is the marginal cost of emitting one additional ton of CO₂, measured in units of current consumption.

**Pontryagin conditions for $\mu$:**

$$\frac{\partial\mathcal{H}}{\partial\mu} = 0 \implies \text{optimal abatement: MAC} = SCC,$$

where MAC is the marginal abatement cost ($\partial\Lambda/\partial\mu \cdot Y$). The optimal abatement fraction equates marginal abatement cost to the SCC.

**The SCC dynamics (costate equation for $q_{AT}$):**

$$\dot{q}_{AT} = \rho q_{AT} - \frac{\partial\mathcal{H}}{\partial M_{AT}} = \rho q_{AT} + q_T \cdot \frac{\partial F}{\partial M_{AT}}\cdot\frac{1}{C_{AT}}.$$

The SCC grows (in absolute value) over time at a rate related to the discount rate $\rho$ and the marginal impact of carbon on temperature.

**Theorem 35.1 (Ramsey Formula for the SCC).** In the DICE model, the steady-state growth rate of the SCC satisfies:

$$\frac{\dot{SCC}}{SCC} \approx \rho + \sigma g_c - g_{damages},$$

where $g_c$ is consumption growth, $g_{damages}$ is the growth rate of climate damages, and $\rho + \sigma g_c$ is the Ramsey discount rate [P:Ch.37.2].

*Derivation.* The costate equation gives $\dot{q}_{AT}/q_{AT} = \rho - \partial\ln\mathcal{H}/\partial M_{AT}$. The Ramsey discount rate $\rho + \sigma g_c$ adjusts for the declining marginal utility of consumption as the economy grows. $\square$

---

## 35.4 Sensitivity Analysis: Stern vs. Nordhaus

The Stern Review (2006) and Nordhaus (2008) reached dramatically different conclusions about optimal climate policy. The difference is almost entirely driven by the discount rate parameters $\rho$ (pure time preference) and $\sigma$ (CRRA curvature).

**Table 35.1: Stern vs. Nordhaus parameters.**

| Parameter | Stern | Nordhaus | Effect on SCC |
|---|---|---|---|
| $\rho$ (pure time preference) | 0.001 | 0.015 | Lower $\rho$ → higher SCC |
| $\sigma$ (CRRA) | 1.0 | 2.0 | Lower $\sigma$ → higher SCC |
| Ramsey rate ($\rho+\sigma g$) | ≈1.4% | ≈5.5% | Much lower for Stern |
| SCC (2020 estimate) | ~$300/t | ~$40/t | 7.5× difference |

The Stern–Nordhaus debate is fundamentally a disagreement about the appropriate **intergenerational discount rate** — how much we should discount welfare of future generations relative to the present [P:Ch.37.3].

---

## 35.5 Numerical Solution: Discretizing DICE

The continuous-time DICE problem is discretized to 5-year time steps, converting the optimal control problem to a **nonlinear programming (NLP)** problem.

**Discrete DICE system (each period $t = 0, 1, \ldots, T$):**

State transition: $K_{t+1} = (1-\delta_K)K_t + Y_t - C_t$, $M_{AT,t+1} = \sum_j\phi_{1j}M_{j,t} + E_t$, $T_{t+1}$ from climate equations.

Objective: $\max_{\mu_t, C_t}\sum_t R_t L_t u(C_t/L_t)$, where $R_t = e^{-\rho\Delta t}/(1+R_t^{market})$ (discount factor).

This is a **$2T$-dimensional NLP** solved by `scipy.optimize.minimize` (Python) or `Ipopt` (Julia). The Hamiltonian provides analytical gradients, accelerating convergence.

```python
import numpy as np
from scipy.optimize import minimize

# Simplified DICE-like model (discrete time, 10-period simulation)
T_steps = 20      # 20 periods of 5 years each = 100 years
dt = 5.0          # years per period

# Parameters
rho   = 0.015     # Nordhaus: pure time preference
sigma_u = 2.0     # CRRA curvature (Nordhaus)
gamma = 0.30      # capital share
delta_K = 0.10    # depreciation (per 5-year period)
pi2   = 0.00236   # damage function coefficient
theta1, theta2 = 0.0318, 2.8  # abatement cost parameters

# Exogenous paths (simplified)
A_path = 5.0 * (1.04)**np.arange(T_steps)   # TFP growth
L_path = 7e9 * (1.005)**np.arange(T_steps)  # population
sigma_E = 0.6 * (0.98)**np.arange(T_steps)  # emission intensity

def dice_welfare(controls, rho=rho, sigma_u=sigma_u):
    """Evaluate DICE welfare given control sequence."""
    n = T_steps
    mu = controls[:n]   # abatement fraction [0,1]
    C  = controls[n:]   # consumption path (trillion $)
    mu = np.clip(mu, 0, 1); C = np.maximum(C, 1.0)
    
    # State trajectories
    K = np.zeros(n); T_AT = np.zeros(n); M_AT = np.zeros(n)
    K[0] = 200.0; T_AT[0] = 1.0; M_AT[0] = 800.0  # initial conditions
    
    welfare = 0.0
    for t in range(n):
        # Damage and abatement
        damage = 1.0 / (1 + pi2 * T_AT[t]**2)
        abate  = 1.0 - theta1 * mu[t]**theta2
        Y_t = damage * abate * A_path[t] * K[t]**gamma
        
        # Welfare contribution
        c_pc = C[t] / (L_path[t] / 1e9)  # per capita consumption
        welfare += L_path[t]/1e9 * np.exp(-rho*t*dt) * (c_pc**(1-sigma_u)-1)/(1-sigma_u) * dt
        
        # Update states
        E_t = sigma_E[t] * (1-mu[t]) * Y_t
        if t < n-1:
            K[t+1] = (1-delta_K)*K[t] + Y_t - C[t]
            K[t+1] = max(K[t+1], 1.0)
            # Simplified climate
            M_AT[t+1] = 0.88*M_AT[t] + E_t*dt*3.67  # mass balance
            T_AT[t+1] = T_AT[t] + 0.01*(np.log(M_AT[t]/590)/np.log(2)*3.8 - 1.2*T_AT[t])*dt
    
    return -welfare  # minimize negative welfare

# Optimize
mu0 = np.ones(T_steps) * 0.3  # initial abatement guess
C0  = np.ones(T_steps) * 50.0  # initial consumption guess
x0 = np.concatenate([mu0, C0])

result = minimize(dice_welfare, x0, method='SLSQP',
                  bounds=[(0,1)]*T_steps + [(1, 1e5)]*T_steps,
                  options={'maxiter': 500, 'ftol': 1e-8})

mu_opt = result.x[:T_steps]; C_opt = result.x[T_steps:]
print(f"Optimal abatement path (first 5 periods): {np.round(mu_opt[:5]*100,1)}%")
print(f"Optimal consumption (first 5 periods): {np.round(C_opt[:5],1)} T$")

# Compute SCC via numerical shadow price (finite difference)
def dice_welfare_delta(delta_E=0.01, t_shock=0):
    """Perturb emission at t_shock by delta_E, measure welfare change."""
    sigma_E_perturb = sigma_E.copy()
    sigma_E_perturb[t_shock] += delta_E / (A_path[t_shock] * 200**gamma * 0.7)
    # Re-run with perturbed path... (simplified)
    return dice_welfare(result.x)  # placeholder

# Stern comparison
result_stern = minimize(lambda x: dice_welfare(x, rho=0.001, sigma_u=1.0), x0,
                        method='SLSQP',
                        bounds=[(0,1)]*T_steps + [(1, 1e5)]*T_steps,
                        options={'maxiter': 300, 'ftol': 1e-6})
mu_stern = result_stern.x[:T_steps]
print(f"\nStern abatement (period 1): {mu_stern[0]*100:.1f}%  vs Nordhaus: {mu_opt[0]*100:.1f}%")
print("(Stern advocates near-100% abatement immediately due to low discount rate)")
```

```julia
using Optim

# DICE-like welfare calculation in Julia
function dice_welfare(controls; rho=0.015, sigma_u=2.0)
    T = 20; n = T; dt = 5.0
    mu = clamp.(controls[1:n], 0, 1)
    C  = max.(controls[n+1:2n], 1.0)
    
    # Exogenous paths
    A = 5.0 .* (1.04).^(0:n-1); L = 7e9 .* (1.005).^(0:n-1)
    sig_E = 0.6 .* (0.98).^(0:n-1)
    
    K = 200.0; T_AT = 1.0; M_AT = 800.0; welfare = 0.0
    for t in 1:n
        damage = 1/(1 + 0.00236*T_AT^2); abate = 1 - 0.0318*mu[t]^2.8
        Y = damage*abate*A[t]*K^0.30
        c_pc = C[t]/(L[t]/1e9)
        welfare += L[t]/1e9 * exp(-rho*(t-1)*dt) * (c_pc^(1-sigma_u)-1)/(1-sigma_u) * dt
        K = max((1-0.10)*K + Y - C[t], 1.0)
        M_AT = 0.88*M_AT + sig_E[t]*(1-mu[t])*Y*dt*3.67
        T_AT += 0.01*(log(M_AT/590)/log(2)*3.8 - 1.2*T_AT)*dt
    end
    return -welfare
end

x0 = vcat(0.3*ones(20), 50.0*ones(20))
res = optimize(dice_welfare, x0, LBFGS(), Optim.Options(iterations=500))
mu_opt = clamp.(Optim.minimizer(res)[1:20], 0, 1)
println("Nordhaus optimal abatement (period 1): $(round(mu_opt[1]*100,digits=1))%")

# Stern comparison
dice_stern(x) = dice_welfare(x; rho=0.001, sigma_u=1.0)
res_stern = optimize(dice_stern, x0, LBFGS(), Optim.Options(iterations=500))
mu_stern = clamp.(Optim.minimizer(res_stern)[1:20], 0, 1)
println("Stern optimal abatement (period 1): $(round(mu_stern[1]*100,digits=1))%")
```

---

## 35.6 Programming Exercises

### Exercise 35.1 (APL — DICE Scan Simulation)

Implement the DICE state update as a dfn `dice_step ← {K T_AT M_AT mu C ← ⍵ ⋄ ...}` and simulate the model forward using APL's scan operator: `states ← dice_step \ T ⍴ ⊂init_state`. (a) Simulate 20 periods (100 years) with fixed $\mu = 0.3$ and $C = 50$ T$. (b) Plot temperature trajectory and cumulative welfare. (c) Use the outer product `∘.f` over a grid of $(\rho, \sigma)$ values to compute the SCC sensitivity surface.

### Exercise 35.2 (Python — SCC via Shadow Price)

Compute the SCC numerically: (a) solve DICE at the optimal $(\mu^*, C^*)$; (b) perturb $E_0$ by +1 ton of CO₂ (which changes $M_{AT,0}$ by the carbon cycle coefficient); (c) re-solve with the perturbed carbon path and compare welfare; (d) $SCC_0 = -\Delta V / \Delta E_0$ in dollars per ton. Compare to Nordhaus (2018) estimate of $\approx$\$40/t and Rennert et al. (2022) estimate of $\approx$\$185/t.

### Exercise 35.3 (Julia — Regional RICE Model)

```julia
# Multi-region DICE (RICE): 5 regions
# Each region has own TFP, population, emission intensity
n_regions = 5
A_regional = [5.0, 8.0, 3.0, 1.5, 2.0]   # TFP by region
L_regional = [1.2, 4.5, 0.7, 0.4, 0.2].*1e9  # population
sigma_regional = [0.3, 0.8, 0.6, 0.5, 0.4]   # emission intensity
damage_regional = [0.00236, 0.00354, 0.00472, 0.00590, 0.00472] # higher in poor regions

println("Multi-region (RICE) SCC differs by region:")
println("Rich regions (low damage coeff): lower SCC from domestic perspective")
println("Poor regions (high damage coeff): higher SCC from domestic perspective")
println("Social planner sets global SCC = weighted average = Pareto optimum")
println("Nash bargaining leads to under-provision of abatement (prisoner's dilemma)")
```

### Exercise 35.4 — Uncertainty and the SCC ($\star$)

**Fat-tailed risks and Weitzman's Dismal Theorem:** Weitzman (2009) argues that if climate damages have a fat-tailed distribution (Pareto or similar), the expected SCC is infinite — because the marginal utility of a small probability of catastrophic outcomes dominates the calculation. (a) Modify the damage function to include a fat-tailed catastrophic scenario: with probability $p = 0.01$, temperature rises by $10°C$ and damages are $90\%$ of GDP. (b) Compute the expected SCC with and without the catastrophic scenario. (c) Show that with CRRA utility ($\sigma > 1$) and a fat upper tail on damages, the SCC grows without bound as the catastrophe probability $p \to 0$ but the tail becomes heavier.

---

## 35.7 Chapter Summary

**Key results:**

- The **DICE model** integrates an economic growth model with a carbon-cycle and temperature model; the planner maximizes welfare over consumption and abatement paths, subject to climate dynamics.
- The **social cost of carbon** is the costate variable of the atmospheric CO₂ stock: $SCC = -q_{AT}/q_K$. The optimal abatement rule equates marginal abatement cost to the SCC (Pontryagin condition).
- The **Ramsey formula for SCC growth** (Theorem 35.1): $\dot{SCC}/SCC \approx \rho + \sigma g_c - g_{damages}$ — the SCC grows at the Ramsey discount rate minus the growth of damages.
- The **Stern–Nordhaus gap** ($\approx 7.5\times$) arises almost entirely from the discount rate: Stern uses $\rho = 0.001$ (near-zero pure time preference); Nordhaus uses $\rho = 0.015$ (calibrated to observed market rates).
- The discrete DICE NLP is solved by `scipy.optimize.minimize` (SLSQP) or Julia's `Optim.jl` (L-BFGS); the APL implementation uses `dice_step \ T ⍴ init` for the state dynamics.
- In APL: `SCC_grid ← rho_grid ∘.{SCC_function ⍺ ⍵} sigma_grid` — sensitivity surface via outer product.

*Next: Chapter 36 — Agent-Based Computational Macroeconomics*
