# Part VII: Computational Methods for DSGE Models

*Connects to: Principles Ch. 7, Ch. 9–10, Ch. 23, Ch. 27, Appendix B*

---

Dynamic Stochastic General Equilibrium models are the workhorse of modern central bank policy analysis. The Federal Reserve, the ECB, the Bank of England, the Riksbank, and the IMF all use medium-scale DSGE models — descendants of Smets and Wouters (2007) — as the primary quantitative framework for conditional forecasting, policy simulation, and communication with the public.

This part covers the complete DSGE pipeline: from the nonlinear model equations to the policy functions used for simulation and estimation. **Chapter 27** develops log-linearization — the approximation that transforms the nonlinear DSGE into a tractable linear system, with a complete worked example for both the NK and RBC models. **Chapter 28** solves the linearized system using the Blanchard–Kahn conditions and Sims' `gensys` algorithm, proving the determinacy condition formally. **Chapter 29** goes beyond first-order approximations to second-order perturbation, which is needed for welfare analysis, risk premia, and ELB dynamics. **Chapter 30** develops Bayesian estimation of DSGE models via the Metropolis–Hastings algorithm, with the Kalman filter providing the likelihood. **Chapter 31** integrates everything into working software — Dynare, APL, Python, Julia, and R — with a complete, reproducible implementation of the same NK model in every language.

Every chapter connects to the three-equation NK model of *Principles*. The reader who completes this part will be able to take any log-linearized DSGE model, solve it numerically, estimate it on data, and interpret the results — the complete toolkit of a central bank research economist.

---

# Chapter 27: Linearization and Log-Linearization

*Approximating DSGE Models Around the Steady State*

> *"The first-order approximation is the simplest thing that works. For most policy questions, it is also the right thing."*
> — Michael Woodford

**Cross-reference:** *Principles* Ch. 7 (AS–AD derivation from IS–LM); Ch. 9 (IS–LM as linear system); Ch. 23 (NK three-equation model); Appendix B (DSGE estimation) **[P:Ch.7, P:Ch.9, P:Ch.23, P:AppB]**

---

## 27.1 Why Linearize?

DSGE models are systems of nonlinear equations. The household's Euler equation is $u'(C_t) = \beta\mathbb{E}_t[(1+R_{t+1})u'(C_{t+1})]$; the firm's pricing equation involves a nonlinear aggregator over price dispersion; the capital accumulation equation has multiplicative terms. Solving these nonlinear systems exactly is computationally tractable only for very simple models (via VFI, Chapter 17) or in special cases with closed-form solutions (Chapter 10).

**The linearization strategy:** Approximate the nonlinear model by a linear one near the steady state. The approximation error is second-order: $O(\varepsilon^2)$ where $\varepsilon$ measures the deviation from steady state. For business-cycle analysis — where fluctuations are typically 1–5% of steady-state values — the linearization error is $O(0.05^2) = O(0.0025)$, small relative to the first-order effects we care about.

The linearized model is a system of linear expectational difference equations — exactly the kind solved by the methods of Chapters 18 and 28. This is why linearization is the gateway to the entire DSGE estimation pipeline.

---

## 27.2 Linear Approximation: Deviation Variables

The simplest linearization uses absolute deviations $\tilde{x}_t = x_t - x^*$ from the steady state $x^*$.

**Definition 27.1 (Linear Approximation).** For a smooth function $F(x_1, \ldots, x_n) = 0$, the first-order Taylor expansion around the steady state $(x_1^*, \ldots, x_n^*)$:

$$F(x_1, \ldots, x_n) \approx F(x_1^*, \ldots, x_n^*) + \sum_{i=1}^n F_{x_i}(x^*)\tilde{x}_i = \sum_{i=1}^n F_{x_i}(x^*)\tilde{x}_i,$$

using $F(x^*) = 0$ (steady-state condition).

The result is a linear equation in deviations $\{\tilde{x}_i\}$. For dynamic models with $\mathbb{E}_t[\cdot]$ operators, the linearization is applied similarly, treating expectations as linear operators:

$$\mathbb{E}_t[F(x_t, x_{t+1}, \ldots)] \approx \sum_i F_{x_i}(x^*)\mathbb{E}_t[\tilde{x}_{i,t}].$$

---

## 27.3 Log-Linearization: The $\hat{x}_t$ Convention

For macroeconomic variables that are always positive (output, consumption, prices, the capital stock), **log-linearization** is typically preferred over level linearization. The log-deviation:

$$\hat{x}_t \equiv \ln x_t - \ln x^* = \ln(x_t/x^*)$$

has a direct percentage interpretation: $\hat{x}_t \approx (x_t - x^*)/x^*$ is the percentage deviation from the steady state.

**Definition 27.2 (Log-Linearization Rules).** For variables $X, Y, Z$ near their steady states $X^*, Y^*, Z^*$, with $\hat{X} = \ln X - \ln X^*$, $\hat{Y} = \ln Y - \ln Y^*$:

| Expression | Log-linearized form | Rule |
|---|---|---|
| $X_t Y_t$ | $X^*Y^*(\hat{X}_t + \hat{Y}_t)$ | Log of product = sum of logs |
| $X_t/Y_t$ | $(X^*/Y^*)(\hat{X}_t - \hat{Y}_t)$ | Log of ratio = difference of logs |
| $X_t^\alpha$ | $X^{*\alpha}\alpha\hat{X}_t$ | Log of power = multiply by exponent |
| $X_t + Y_t$ | $X^*\hat{X}_t + Y^*\hat{Y}_t + X^* + Y^*$ | Requires care (see below) |
| $\mathbb{E}_t[X_{t+1}]$ | $X^*(1 + \mathbb{E}_t[\hat{X}_{t+1}])$ | Expectation passes through |

**Handling sums:** For $Z_t = X_t + Y_t$, the log-linearization is:

$$Z^*\hat{Z}_t = X^*\hat{X}_t + Y^*\hat{Y}_t,$$

i.e., $\hat{Z}_t = (X^*/Z^*)\hat{X}_t + (Y^*/Z^*)\hat{Y}_t$ — a weighted average of percentage deviations, with steady-state shares as weights.

**Proof of the product rule.** Write $X_t = X^*e^{\hat{X}_t}$ and $Y_t = Y^*e^{\hat{Y}_t}$. Then $X_tY_t = X^*Y^*e^{\hat{X}_t+\hat{Y}_t}$. Log-linearizing $e^{\hat{X}_t+\hat{Y}_t} \approx 1 + \hat{X}_t + \hat{Y}_t$ gives $X_tY_t \approx X^*Y^*(1 + \hat{X}_t + \hat{Y}_t)$. The deviation form: $X^*Y^*(\hat{X}_t + \hat{Y}_t)$. $\square$

---

## 27.4 Log-Linearizing the New Keynesian Model

*Cross-reference: Principles Ch. 23 (NK three-equation model derivation)* **[P:Ch.23]**

We log-linearize the NK model from microfoundations, connecting the algebraic steps to the equations of *Principles* Ch. 23.

### 27.4.1 Household Euler Equation

The household's CRRA Euler equation:

$$\frac{1}{C_t} = \beta\mathbb{E}_t\!\left[\frac{1+R_{t+1}}{C_{t+1}}\frac{1}{\Pi_{t+1}}\right],$$

where $R_t$ is the nominal interest rate and $\Pi_t = P_t/P_{t-1}$ is gross inflation. The steady state: $1/C^* = \beta(1+R^*)/(C^*\Pi^*)$, giving $(1+R^*)/\Pi^* = 1/\beta$.

Log-linearize each side. Left side: $-\hat{C}_t$ (using $1/C_t = (1/C^*)e^{-\hat{C}_t} \approx (1/C^*)(1-\hat{C}_t)$, so deviation is $-\hat{C}_t$).

Right side, use $1+R_{t+1} \approx (1+R^*)(1 + \hat{R}_{t+1}/(1+R^*)) \approx (1+R^*)(1 + \sigma_r\hat{R}_{t+1})$ where $\sigma_r = R^*/(1+R^*)$... more cleanly: define $\hat{r}_t = \ln(1+R_t) - \ln(1+R^*)$. Then:

$$-\hat{C}_t = \mathbb{E}_t[\hat{r}_{t+1} - \hat{\Pi}_{t+1} - \hat{C}_{t+1}].$$

Rearranging:

$$\hat{C}_t = \mathbb{E}_t[\hat{C}_{t+1}] - (\hat{r}_{t+1} - \mathbb{E}_t[\hat\Pi_{t+1}]) + \text{noise}.$$

In terms of the output gap $\hat{x}_t \equiv \hat{C}_t - \hat{C}^{flex}_t$ (deviation of consumption from its flexible-price level) and letting $\sigma = 1$ (log utility):

$$\boxed{\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - (\hat{r}_t - \mathbb{E}_t[\hat\pi_{t+1}] - r^n_t),}$$

the **Dynamic IS (DIS) equation**. Here $r^n_t$ is the natural rate of interest (the real interest rate that would prevail with flexible prices). This is precisely the DIS equation of *Principles* Ch. 23 [P:Ch.23.1].

### 27.4.2 The New Keynesian Phillips Curve

Firms set prices using Calvo (1983) contracts: each firm can adjust its price with probability $(1-\theta)$ each period. Log-linearizing the aggregate price-level equation and the firm's optimal reset-price condition gives:

$$\hat{\pi}_t = \beta\mathbb{E}_t[\hat\pi_{t+1}] + \kappa\hat{x}_t + u_t, \quad \kappa = \frac{(1-\theta)(1-\beta\theta)}{\theta}\cdot\frac{1}{\sigma + \phi},$$

where $\phi$ is the inverse Frisch elasticity of labor supply and $u_t$ is the cost-push shock [P:Ch.10.3].

**Derivation of $\kappa$.** The firm's optimal reset price satisfies (in log deviations):

$$\hat{p}_t^* = (1-\beta\theta)\sum_{j=0}^\infty(\beta\theta)^j\mathbb{E}_t[\hat{mc}_{t+j} + \hat{p}_{t+j}],$$

where $\hat{mc}_t$ is the log real marginal cost. The aggregate price equation:

$$\hat{p}_t = (1-\theta)\hat{p}_t^* + \theta\hat{p}_{t-1} \implies \hat\pi_t = (1-\theta)(\hat{p}_t^* - \hat{p}_{t-1}).$$

Under standard assumptions $\hat{mc}_t = (\sigma + \phi)\hat{x}_t$, combining gives the NKPC with $\kappa = (1-\theta)(1-\beta\theta)/\theta \cdot 1/(\sigma+\phi) / (1 + ...)$. The precise formula depends on whether labor is firm-specific or economy-wide; the standard form is as stated. $\square$

### 27.4.3 Taylor Rule

The Taylor rule requires no log-linearization — it is already linear:

$$\hat{r}_t = r^n + \phi_\pi\hat\pi_t + \phi_y\hat{x}_t + v_t.$$

The three-equation NK system — DIS, NKPC, Taylor rule — is now the linearized representation of the full nonlinear model, valid to first order in deviations from the steady state.

---

## 27.5 Log-Linearizing the RBC Model

*Cross-reference: Principles Ch. 27 (RBC model)* **[P:Ch.27]**

The RBC model has four key equilibrium conditions. We log-linearize each around the deterministic steady state.

### 27.5.1 Household Euler Equation

$U_C(t) = \beta(1+R_{t+1})U_C(t+1)$. With $U_C = C^{-\sigma}$ and $1+R_t = \alpha A_t K_t^{\alpha-1}n_t^{1-\alpha} + (1-\delta)$ (capital rental rate):

$$-\sigma\hat{C}_t = \mathbb{E}_t[-\sigma\hat{C}_{t+1} + \hat{r}_{t+1}],$$

where $\hat{r}_t = \ln(1+R_t) - \ln(1+R^*) \approx R_t/(1+R^*)$. Rearranging:

$$\hat{C}_t = \mathbb{E}_t[\hat{C}_{t+1}] - \frac{1}{\sigma}\hat{r}_{t+1}.$$

### 27.5.2 Labor Market Clearing

From $-U_N/U_C = W_t$ (intratemporal condition) with $U = C^\mu(1-N)^{(1-\mu)(1-\sigma)}/(1-\sigma)$ (KPR):

$$\frac{(1-\mu)}{\mu}\frac{C_t}{1-N_t} = W_t = (1-\alpha)A_tK_t^\alpha n_t^{-\alpha}.$$

Log-linearizing (using steady-state labor share $n^* = \mu(1-\alpha)/(...)$):

$$\hat{C}_t + \frac{n^*}{1-n^*}\hat{n}_t = \hat{w}_t = \hat{A}_t + \alpha\hat{K}_t - \alpha\hat{n}_t.$$

### 27.5.3 Capital Accumulation and Resource Constraint

$$\hat{K}_{t+1} = (1-\delta)\hat{K}_t + \delta\hat{I}_t, \quad (C^*/Y^*)\hat{C}_t + (I^*/Y^*)\hat{I}_t = \hat{Y}_t.$$

### 27.5.4 Production Function

$$\hat{Y}_t = \hat{A}_t + \alpha\hat{K}_t + (1-\alpha)\hat{n}_t.$$

The complete log-linearized RBC system is five equations in five variables $(\hat{C}_t, \hat{n}_t, \hat{K}_{t+1}, \hat{Y}_t, \hat{r}_t)$ plus the exogenous $\hat{A}_t = \rho_A\hat{A}_{t-1} + \varepsilon_t$.

---

## 27.6 Approximation Error and Second-Order Terms

**Theorem 27.1 (First-Order Approximation Error).** The first-order log-linearization of a smooth model is accurate to order $O(\varepsilon^2)$ where $\varepsilon = \max_i|\hat{x}_{i,t}|$ measures the maximum deviation from steady state. Specifically, for a twice-differentiable function $F$:

$$F(\mathbf{x}^* + \varepsilon\mathbf{v}) = F(\mathbf{x}^*) + \varepsilon\nabla F(\mathbf{x}^*)'\mathbf{v} + \frac{\varepsilon^2}{2}\mathbf{v}'H_F(\mathbf{x}^*)\mathbf{v} + O(\varepsilon^3).$$

The first-order approximation discards the $O(\varepsilon^2)$ Hessian term. For a model with $\varepsilon \approx 0.02$ (2% fluctuations), the relative approximation error is $\varepsilon/2 \approx 1\%$ — acceptable for most policy questions.

**When is the approximation poor?**

1. **Large shocks:** Crisis episodes (2008-09: -9% output gap, 2020: -10%) push far from the linearization point. Second-order approximations (Chapter 29) are more accurate.
2. **Nonlinear mechanisms:** The ELB $i_t \geq 0$ is a binding inequality constraint that the linear approximation cannot handle. Requires global methods or piecewise-linear approximation.
3. **Welfare analysis:** Welfare depends on the variance of consumption ($\approx \varepsilon^2$), which is a second-order effect — invisible in the first-order approximation.

---

## 27.7 Worked Example: Medium-Scale NK Model with Investment

*Cross-reference: Principles Ch. 23, Ch. 27 (medium-scale model)* **[P:Ch.23, P:Ch.27]**

We add investment to the NK model, following Christiano, Eichenbaum, and Evans (2005). The additional equation is the investment Euler equation (from Tobin's $q$ theory, Chapter 13):

$$\hat{q}_t = \mathbb{E}_t[\hat{q}_{t+1}] + \hat{r}_{t+1} - \hat{\pi}_{t+1} - \hat{rk}_{t+1},$$

$$\hat{I}_t = \frac{1}{\psi}\hat{q}_t,$$

where $\hat{q}_t$ is Tobin's $q$ deviation and $\hat{rk}_t$ is the log-deviation of the return on capital. The capital accumulation:

$$\hat{K}_{t+1} = (1-\delta)\hat{K}_t + \delta\hat{I}_t.$$

**Full system (8 equations, 8 variables):**

```
DIS:      x̂_t = E_t[x̂_{t+1}] - σ(î_t - E_t[π̂_{t+1}] - r^n_t)
NKPC:     π̂_t = βE_t[π̂_{t+1}] + κx̂_t
Taylor:   î_t = ρ_i î_{t-1} + (1-ρ_i)(φ_π π̂_t + φ_y x̂_t) + ε^i_t
q-eq:     q̂_t = E_t[q̂_{t+1}] + E_t[r̂k_{t+1}] - (î_t - E_t[π̂_{t+1}])
Invest:   Î_t = (1/ψ) q̂_t
K-accum:  K̂_{t+1} = (1-δ)K̂_t + δÎ_t
rk eq:    r̂k_t = (1-α-δ)(Ŷ_t - K̂_{t-1}) + ... (MPK relation)
Output:   Ŷ_t = (C*/Y*)ĈY_t + (I*/Y*)Î_t
```

This is the prototype system for the CEE (2005) and Smets–Wouters (2007) models. Chapter 28 solves it using the Blanchard–Kahn machinery.

```python
import numpy as np

# Log-linearization verification: check that linearized model
# matches nonlinear model near steady state for small shocks
alpha, delta, beta, sigma = 0.36, 0.025, 0.99, 1.0
kappa_NK, phi_pi, phi_y = 0.15, 1.5, 0.5

# Steady state (normalized: Y*=1, C*=1-delta*K*, etc.)
r_star = 1/beta - 1
K_star = (alpha/(r_star+delta))**(1/(1-alpha))
Y_star = K_star**alpha
C_star = Y_star - delta*K_star
I_star = delta*K_star

print("NK model steady state:")
print(f"  Y* = {Y_star:.4f}, K* = {K_star:.4f}, C* = {C_star:.4f}, I* = {I_star:.4f}")
print(f"  r* = {r_star*100:.2f}% (quarterly), K/Y = {K_star/Y_star:.2f}")

# Log-linearize: verify DIS, NKPC hold
# Given a unit TFP shock hat_A = 1%:
hat_A = 0.01
# Natural rate: r^n = α*rho_A*hat_A (approximately, from the Euler eq. at flex prices)
rho_A = 0.95
r_n = alpha * rho_A * hat_A  # approximately

# MSV solution from Chapter 18
G0 = np.array([[1, -kappa_NK], [sigma*phi_pi, 1+sigma*phi_y]])
G1 = np.array([[beta, 0], [-sigma, 1]])
A = np.linalg.inv(G0) @ G1
eigs = np.linalg.eigvals(A)
print(f"\nEigenvalues of A: {np.abs(eigs).round(3)} — determinacy: {np.all(np.abs(eigs)>1)}")

# Response to natural rate shock
C_A = np.linalg.inv(G0) @ np.array([0, sigma])  # loading on r_n shock
Omega = np.linalg.solve(np.eye(2) - rho_A*A, C_A.reshape(2,1))
print(f"\nMSV solution (response to 1% TFP shock via natural rate):")
print(f"  π̂_t = {Omega[0,0]:.4f} * û_t")
print(f"  x̂_t = {Omega[1,0]:.4f} * û_t")
```

```julia
# Julia — automated log-linearization via symbolic differentiation
using ForwardDiff, LinearAlgebra

# Define nonlinear NK equations F(y_t, y_{t-1}, E_y_{t+1}, z_t) = 0
# y = [pi, x], z = [r_n]
function nk_equations(yt, yt_1, Eyt1, zt; beta=0.99, kappa=0.15, sigma=1.0, phi_pi=1.5, phi_y=0.5)
    pi_t, x_t = yt
    pi_t1, x_t1 = Eyt1
    r_n = zt[1]
    r_t = phi_pi*pi_t + phi_y*x_t  # Taylor rule gives nominal rate
    
    F1 = pi_t - beta*pi_t1 - kappa*x_t             # NKPC
    F2 = x_t - x_t1 + sigma*(r_t - pi_t1 - r_n)    # DIS
    return [F1, F2]
end

# Compute Jacobians at steady state
y_ss = [0.0, 0.0]; z_ss = [0.0]
J_yt  = ForwardDiff.jacobian(y -> nk_equations(y,y_ss,y_ss,z_ss), y_ss)
J_Eyt = ForwardDiff.jacobian(Ey -> nk_equations(y_ss,y_ss,Ey,z_ss), y_ss)
J_zt  = ForwardDiff.jacobian(z -> nk_equations(y_ss,y_ss,y_ss,z), z_ss)

println("Γ₀ (J_yt) = ", round.(J_yt, digits=3))
println("Γ₁ (-J_Eyt) = ", round.(-J_Eyt, digits=3))
println("Ψ (J_zt) = ", round.(J_zt, digits=3))
# These are the Γ₀, Γ₁, Ψ matrices for Chapter 28's gensys
```

```r
# R — log-linearization rules demonstration
alpha <- 0.36; delta <- 0.025; beta <- 0.99

# Steady state
r_star <- 1/beta - 1
K_star <- (alpha/(r_star+delta))^(1/(1-alpha))
Y_star <- K_star^alpha
C_star <- Y_star - delta*K_star

# Log-linearization of production: Ŷ = Â + α*K̂ + (1-α)*n̂
# Verify: at small deviations, log-lin ≈ exact percentage change
hat_A_vec <- seq(-0.05, 0.05, 0.01)  # -5% to +5% TFP shocks

# Exact nonlinear response (holding K, n fixed at SS)
Y_exact <- Y_star * exp(hat_A_vec)  # Y = A * K^α * n^(1-α), only A changes

# Log-linearized response: Ŷ ≈ Â (with K, n fixed)
Y_linear <- Y_star * (1 + hat_A_vec)

cat("Approximation error (log-linear vs exact):\n")
for(i in c(1,3,5,7,9,11)){
  err <- abs(Y_exact[i] - Y_linear[i]) / Y_star
  cat(sprintf("  Â=%.2f: exact=%.4f, linear=%.4f, error=%.4e\n",
              hat_A_vec[i], Y_exact[i], Y_linear[i], err))
}
```

---

## 27.8 Programming Exercises

### Exercise 27.1 (APL — Log-Linearization of the Euler Equation)

Implement the log-linearization of the household Euler equation $C_t^{-\sigma} = \beta(1+R_{t+1})C_{t+1}^{-\sigma}/\Pi_{t+1}$ step by step in APL. (a) Express each term as `(1+hat_C)^{-sigma} ≈ 1-sigma×hat_C` using `1+(-sigma)×hat_C`. (b) Equate deviations on both sides to derive the DIS equation. (c) Verify numerically: generate 100 periods of model data with $\varepsilon=0.01$ deviations and check that the log-linearized DIS holds to $O(\varepsilon^2) \approx 10^{-4}$.

### Exercise 27.2 (Python — Calvo Pricing Derivation)

Derive the NKPC slope $\kappa = (1-\theta)(1-\beta\theta)/\theta$ for the standard Calvo model with linear disutility of labor. (a) Start from the firm's optimal reset price $p_t^* = (1-\beta\theta)\sum_{j=0}^\infty(\beta\theta)^j\mathbb{E}_t[mc_{t+j}+p_{t+j}]$. (b) Log-linearize and use the law of motion for the price level to obtain the NKPC. (c) Compute $\kappa$ for $\theta \in \{0.5, 0.66, 0.75, 0.85\}$ and $\beta = 0.99$.

### Exercise 27.3 — Approximation Error ($\star$)

For the RBC Euler equation $C_t^{-\sigma} = \beta\mathbb{E}_t[(R_{t+1}+1-\delta)C_{t+1}^{-\sigma}]$: (a) log-linearize to get the linearized Euler equation; (b) simulate the exact nonlinear model for 500 periods with TFP shock std dev $\sigma_A \in \{0.007, 0.05, 0.10\}$; (c) compute the Euler equation residual in both the linearized and exact versions; (d) show the linearized residual is $O(\sigma_A^2)$.

---

## 27.9 Chapter Summary

**Key results:**

- **Log-deviation** $\hat{x}_t = \ln x_t - \ln x^* \approx (x_t-x^*)/x^*$ is the natural variable for percentage deviations; log-linearization uses the rules: product → sum, ratio → difference, power $\alpha$ → multiply.
- The **DIS equation** $\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(\hat{r}_t - \mathbb{E}_t[\hat\pi_{t+1}] - r^n_t)$ and **NKPC** $\hat\pi_t = \beta\mathbb{E}_t[\hat\pi_{t+1}] + \kappa\hat{x}_t$ are the first-order approximations of the household Euler equation and Calvo pricing, valid to $O(\varepsilon^2)$.
- The **approximation error** is $O(\varepsilon^2)$: for 2% fluctuations ($\varepsilon = 0.02$), the relative error is $\approx 0.02\%$ — acceptable for most policy analysis.
- The **RBC log-linearized system** has five equations in $(\hat{C}, \hat{n}, \hat{K}', \hat{Y}, \hat{r})$; the NK system has three in $(\hat\pi, \hat{x}, \hat{r})$; medium-scale models augment these with investment, habits, price and wage stickiness.
- In APL: log-linearization of a product `X*Y` is `(X_ss*hat_X) + (Y_ss*hat_Y)`; of a ratio `X/Y` is `hat_X - hat_Y`; ForwardDiff.jl automates the Jacobian computation to give $\Gamma_0, \Gamma_1, \Psi$ directly.

*Next: Chapter 28 — Solving Linear DSGE Models: Blanchard–Kahn and Sims*
