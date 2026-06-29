# Chapter 4: Difference Equations in Discrete-Time Macro Models

*Business Cycle Dynamics*

> *"In discrete time, every model is an iterated map. The question is whether that map has a fixed point, a cycle, or chaos."*

**Cross-reference:** *Principles* Ch. 6 (time-series properties, HP filter); Ch. 10 (NKPC as a forward difference equation); Ch. 16 (rational expectations, forward solutions); Ch. 27 (RBC model log-linearized as a discrete system) **[P:Ch.6, P:Ch.10, P:Ch.16, P:Ch.27]**

---

## 4.1 Why Discrete Time? Data, Policy, and Computation

The choice between continuous and discrete time is not merely aesthetic. Several compelling reasons push modern macroeconomics toward discrete-time formulations.

**Data observability.** GDP is measured quarterly, inflation monthly, the federal funds rate daily. The natural time unit for policy analysis — the period between Federal Open Market Committee meetings — is six weeks. A continuous-time model that maps cleanly to quarterly data requires a discretization step anyway; it is simpler to start in discrete time.

**Computational tractability.** Value function iteration, the Kalman filter, the Metropolis–Hastings sampler for Bayesian DSGE estimation — all are inherently discrete. The state space must be discretized, time must be measured in steps, and probability distributions are approximated by finite grids.

**Forward-looking equations.** The New Keynesian Phillips Curve $\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t$ [P:Ch.10] and the Dynamic IS curve $\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\pi_{t+1}] - r^n_t)$ are difference equations with expectations of future values. These arise naturally from discrete-time optimization but have no direct continuous-time analogue for policy analysis.

This chapter develops the theory of difference equations — the discrete-time counterpart to Chapter 3's ODEs. Every result has a direct parallel in the continuous-time theory; the key substitution is $e^{at} \leftrightarrow \lambda^t$ and the stability condition changes from $a < 0$ to $|\lambda| < 1$.

---

## 4.2 First-Order Linear Difference Equations

### 4.2.1 The Homogeneous Case

The equation $x_{t+1} = \lambda x_t$ has the solution:
$$x_t = \lambda^t x_0.$$

Stability is determined by $|\lambda|$:
- $|\lambda| < 1$: $x_t \to 0$ (stable).
- $|\lambda| = 1$: $|x_t| = |x_0|$ constant (neutrally stable).
- $|\lambda| > 1$: $|x_t| \to \infty$ (unstable).

**Theorem 4.1 (Stability of First-Order Difference Equations).** The equilibrium $x^* = 0$ of $x_{t+1} = \lambda x_t$ is **globally asymptotically stable** if and only if $|\lambda| < 1$.

### 4.2.2 The Nonhomogeneous Case

The equation $x_{t+1} = \lambda x_t + c$ (with $\lambda \neq 1$) has the general solution:
$$x_t = \lambda^t\left(x_0 - x^*\right) + x^*, \quad x^* = \frac{c}{1-\lambda}.$$

For $|\lambda| < 1$, $x_t \to x^*$ as $t \to \infty$.

**Convergence speed:** The gap $x_t - x^*$ shrinks by factor $\lambda$ each period. After $T$ periods, the gap is $\lambda^T(x_0 - x^*)$.

### 4.2.3 Backward vs. Forward Solutions

A crucial distinction in macroeconomics is between **backward-looking** (predeterminate) variables, whose current value is determined by past behavior, and **forward-looking** (free) variables, whose current value is determined by expected future behavior.

**Backward solution** for $x_{t+1} = \lambda x_t + c$: start from $x_0$ and iterate forward. Requires $|\lambda| < 1$ for stability.

**Forward solution** for $x_{t+1} = \lambda x_t + c$: rearrange as $x_t = \lambda^{-1}(x_{t+1} - c)$ and iterate forward:
$$x_t = \lambda^{-T}x_{t+T} - c\lambda^{-1}\sum_{j=0}^{T-1}\lambda^{-j}.$$

For the forward solution to be bounded (a **transversality condition**), we need $|\lambda^{-T}x_{t+T}| \to 0$ as $T \to \infty$, which requires $|\lambda| > 1$ (so $|\lambda^{-T}| \to 0$). If $|\lambda| > 1$, the forward solution is:
$$x_t = -\frac{c}{\lambda - 1} + \sum_{j=0}^\infty \lambda^{-(j+1)}\mathbb{E}_t[c_{t+j}] = \sum_{j=0}^\infty \lambda^{-(j+1)}\mathbb{E}_t[\text{forcing variable}_{t+j}].$$

This is the **present-value formula** structure that appears throughout the New Keynesian model. The NKPC forward solution $\hat{\pi}_t = \kappa\sum_{j=0}^\infty \beta^j\mathbb{E}_t[\hat{x}_{t+j}]$ is exactly this form with $\lambda = 1/\beta > 1$ [P:Ch.10].

---

## 4.3 Second-Order Difference Equations

Second-order equations $x_{t+2} + px_{t+1} + qx_t = c$ appear in the multiplier–accelerator model of Samuelson (1939) [P:Ch.8] and in continuous-time models after discretization. The general solution is:

$$x_t = C_1\lambda_1^t + C_2\lambda_2^t + x^*,$$

where $\lambda_1, \lambda_2$ are the roots of the **characteristic equation**:
$$\lambda^2 + p\lambda + q = 0 \implies \lambda_{1,2} = \frac{-p \pm \sqrt{p^2 - 4q}}{2},$$

and $x^* = c/(1 + p + q)$ is the particular solution (assumed $1 + p + q \neq 0$). The constants $C_1$, $C_2$ are determined by two initial conditions $x_0$ and $x_1$.

**Classification by roots:**

| Discriminant $\Delta = p^2 - 4q$ | Roots | Path Type |
|---|---|---|
| $\Delta > 0$ | Two real roots $\lambda_1 \neq \lambda_2$ | Monotone (same sign) or S-shaped |
| $\Delta = 0$ | Repeated real root $\lambda_1 = \lambda_2 = -p/2$ | Monotone or boundary case |
| $\Delta < 0$ | Complex conjugate pair $\lambda = re^{\pm i\theta}$ | Oscillatory |

For complex roots $\lambda = re^{i\theta}$ (where $r = \sqrt{q}$ and $\theta = \arctan(\sqrt{4q-p^2}/(-p))$):
$$x_t = x^* + r^t(C_1\cos(\theta t) + C_2\sin(\theta t)).$$

The system oscillates with period $2\pi/\theta$ and amplitude decaying like $r^t$. Stability requires $r = \sqrt{q} < 1$.

### 4.3.1 The Multiplier–Accelerator Model

Samuelson's (1939) model of endogenous business cycles [P:Ch.8] sets:
$$Y_t = C_t + I_t + G, \quad C_t = bY_{t-1}, \quad I_t = v(C_t - C_{t-1}),$$

where $b$ is the MPC and $v$ is the accelerator. Substituting:
$$Y_t = b(1+v)Y_{t-1} - bvY_{t-2} + G.$$

This is exactly a second-order difference equation with $p = -b(1+v)$ and $q = bv$. The steady state is $Y^* = G/(1-b)$.

**Stability:** Requires both roots $|\lambda_{1,2}| < 1$. Equivalently, the necessary and sufficient conditions are $|p| < 1 + q$ and $q < 1$, i.e., $b(1+v) < 1 + bv$ (always true) and $bv < 1$. Stability requires $bv < 1$.

**Business cycles arise from complex roots:** When $p^2 < 4q$, i.e., $b^2(1+v)^2 < 4bv$, the roots are complex and the system oscillates. The amplitude $r = \sqrt{bv}$ determines whether cycles decay ($bv < 1$), persist ($bv = 1$), or explode ($bv > 1$).

---

## 4.4 Systems of Linear Difference Equations

The discrete-time state-space form:
$$\mathbf{x}_{t+1} = A\mathbf{x}_t + B\mathbf{u}_t,$$

where $\mathbf{x}_t \in \mathbb{R}^n$ is the state vector and $\mathbf{u}_t \in \mathbb{R}^m$ is an exogenous input (shock or policy variable). This is the linearized DSGE model after log-linearization — the starting point of Parts VII.

**General solution:** Iterating the state equation:
$$\mathbf{x}_t = A^t\mathbf{x}_0 + \sum_{j=0}^{t-1}A^j B\mathbf{u}_{t-1-j}.$$

Using the eigendecomposition $A = PDP^{-1}$:
$$\mathbf{x}_t = PD^tP^{-1}\mathbf{x}_0 + \sum_{j=0}^{t-1}PD^jP^{-1}B\mathbf{u}_{t-1-j}.$$

**Stability:** Requires all eigenvalues of $A$ to satisfy $|\lambda_i| < 1$.

**Impulse response:** The response of $\mathbf{x}$ to a unit shock $\mathbf{u}_0 = \mathbf{e}_j$ (with $\mathbf{u}_t = \mathbf{0}$ for $t > 0$) at horizon $h$ is $A^{h-1}B\mathbf{e}_j$ — the $j$-th column of $A^{h-1}B$. This is computed efficiently as $A^h\mathbf{v}$ using repeated matrix multiplication or via the eigendecomposition.

```apl
⍝ APL — state-space simulation and impulse response
⎕IO←0 ⋄ ⎕ML←1

⍝ A simple 2×2 system (linearized IS-LM style)
A ← 2 2 ⍴ 0.85 0.10 ¯0.05 0.90    ⍝ transition matrix
B ← 2 1 ⍴ 0.5 ¯0.3                  ⍝ shock loadings

⍝ Impulse response: response of state to unit shock at t=0
shock ← ,0.1                          ⍝ unit shock (1×1 vector)
H ← 20                                ⍝ horizons

⍝ IRF: collect A^h × B × shock for h = 0, ..., H-1
irf ← {(A⍣⍵) +.× B +.× shock} ¨ ⍳H

⍝ Stack into matrix H×2
⊃ irf

⍝ Simulate with AR(1) shocks
T ← 100
shocks ← 0.9 × {⍵+0.1×1∘⊂(?⍵⍴0)} ⍣T ⊢ 0   ⍝ AR(1) shock sequence

⍝ State simulation
sim_path ← {A +.× ⍵} \ T ⍴ ⊂(2⍴0)   ⍝ homogeneous (no shock for now)
```

---

## 4.5 Forward-Looking Difference Equations and the Minimum State Variable Solution

The distinctive challenge of modern macroeconomics is that many key equations involve **expected future values** of endogenous variables. The NKPC is $\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t$; the DIS is $\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(\hat{i}_t - \mathbb{E}_t[\hat{\pi}_{t+1}] - \hat{r}^n_t)$. These are **expectational difference equations** with no unique solution without additional restrictions.

**Definition 4.1 (Minimum State Variable Solution).** The **minimum state variable (MSV) solution** is the unique solution of a linear expectational difference equation that expresses each endogenous variable as a linear function of the minimal set of state variables — typically the exogenous forcing variables and any predetermined (backward-looking) variables.

McCallum's (1983) algorithm for finding the MSV solution:

**Algorithm 4.1 (McCallum's Undetermined Coefficients).**

Given the expectational system $\mathbf{y}_t = A\mathbb{E}_t[\mathbf{y}_{t+1}] + C\mathbf{z}_t$ where $\mathbf{z}_t$ is a vector of exogenous state variables with law of motion $\mathbf{z}_{t+1} = \Phi\mathbf{z}_t + \bm{\varepsilon}_{t+1}$:

1. **Guess** an MSV solution: $\mathbf{y}_t = \Omega\mathbf{z}_t$ (endogenous variables are linear functions of exogenous states).
2. **Compute** $\mathbb{E}_t[\mathbf{y}_{t+1}] = \Omega\mathbb{E}_t[\mathbf{z}_{t+1}] = \Omega\Phi\mathbf{z}_t$.
3. **Substitute** into the original equation: $\Omega\mathbf{z}_t = A\Omega\Phi\mathbf{z}_t + C\mathbf{z}_t$.
4. **Match** coefficients: $\Omega = A\Omega\Phi + C$.
5. **Solve** this matrix equation for $\Omega$. This is a discrete **Sylvester equation** of the form $\Omega - A\Omega\Phi = C$.

The Sylvester equation $\Omega - A\Omega\Phi = C$ can be vectorized: $\text{vec}(\Omega)(I - \Phi' \otimes A) = \text{vec}(C)$, which is a linear system in $\text{vec}(\Omega)$. It has a unique solution when $(I - \Phi' \otimes A)$ is nonsingular — equivalently, when no eigenvalue of $\Phi$ equals the reciprocal of any eigenvalue of $A$.

### 4.5.1 Application: The NK Model Under a Taylor Rule

*Cross-reference: Principles Ch. 10, Ch. 23* **[P:Ch.10, P:Ch.23]**

The three-equation NK model:
$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(\hat{i}_t - \mathbb{E}_t[\hat{\pi}_{t+1}] - \hat{r}^n_t)$$
$$\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t$$
$$\hat{i}_t = \phi_\pi\hat{\pi}_t + \phi_y\hat{x}_t + \hat{r}^n_t$$

Substituting the Taylor rule into the DIS:
$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma\phi_\pi\hat{\pi}_t - \sigma\phi_y\hat{x}_t + \sigma\hat{r}^n_t - \sigma\mathbb{E}_t[\hat{\pi}_{t+1}] + \sigma\hat{r}^n_t.$$

Wait — let us use the compact matrix form. Stack $\mathbf{y}_t = (\hat{\pi}_t, \hat{x}_t)'$ and note that both are jump variables (free to adjust at each date). The system can be written:

$$\underbrace{\begin{pmatrix} 1 & 0 \\ \sigma\phi_\pi & 1+\sigma\phi_y \end{pmatrix}}_{\equiv\Gamma_0} \mathbf{y}_t = \underbrace{\begin{pmatrix} \beta & \kappa \\ \sigma & 1 \end{pmatrix}}_{\equiv\Gamma_1} \mathbb{E}_t[\mathbf{y}_{t+1}] + \underbrace{\begin{pmatrix} 0 \\ \sigma \end{pmatrix}}_{\equiv\Psi}\hat{r}^n_t.$$

This is the canonical form $\Gamma_0\mathbf{y}_t = \Gamma_1\mathbb{E}_t[\mathbf{y}_{t+1}] + \Psi z_t$, which we solve via $\mathbf{y}_t = \Gamma_0^{-1}\Gamma_1\mathbb{E}_t[\mathbf{y}_{t+1}] + \Gamma_0^{-1}\Psi z_t$. Let $A = \Gamma_0^{-1}\Gamma_1$ and $C = \Gamma_0^{-1}\Psi$. The MSV solution is $\mathbf{y}_t = \Omega z_t$ where $\Omega$ solves the scalar Sylvester equation $\Omega - A\Omega\rho_r = C$ (since $z_t = \hat{r}^n_t$ follows AR(1) with $\phi = \rho_r$).

**Determinacy:** The system has a unique bounded solution (MSV solution is the only equilibrium) iff the number of eigenvalues of $A = \Gamma_0^{-1}\Gamma_1$ outside the unit circle equals the number of free variables — here, both $\hat{\pi}_t$ and $\hat{x}_t$ are free, so we need exactly 2 eigenvalues of $A$ outside the unit circle. This is the **Blanchard–Kahn condition** for this model (developed fully in Chapter 28).

The Taylor principle $\phi_\pi > 1$ ensures this: in the 2×2 NK model, $\phi_\pi > 1$ implies both eigenvalues of $A$ have modulus greater than 1, satisfying the Blanchard–Kahn condition. This is the formal derivation of the result stated in *Principles* Ch. 23 [P:Ch.23.1] that the Taylor principle is necessary and sufficient for determinacy.

---

## 4.6 Stability Analysis and the Unit Circle

For a system $\mathbf{x}_{t+1} = A\mathbf{x}_t$, stability is determined by the eigenvalues of $A$ relative to the unit circle in the complex plane.

**Definition 4.2 (Stable, Unstable, and Centre Subspaces).** Decompose $\mathbb{R}^n$ into three subspaces based on eigenvalues:
- **Stable subspace** $E^s$: spanned by generalised eigenvectors corresponding to $|\lambda_i| < 1$.
- **Unstable subspace** $E^u$: spanned by generalised eigenvectors corresponding to $|\lambda_i| > 1$.
- **Centre subspace** $E^c$: spanned by generalised eigenvectors corresponding to $|\lambda_i| = 1$.

For the DSGE model with $n_s$ predetermined variables and $n_f$ free (jump) variables, the Blanchard–Kahn condition requires:

> **Number of eigenvalues of $A$ outside unit circle = number of free variables $n_f$.**

When this holds, the model has a unique bounded rational expectations equilibrium. When there are too few eigenvalues outside the unit circle, there are multiple equilibria (indeterminacy — sunspot equilibria exist). When there are too many, no bounded equilibrium exists (the model "explodes").

**Checking stability in APL:**

```apl
⍝ APL — eigenvalue check for Blanchard-Kahn condition
⎕IO←0 ⋄ ⎕ML←1

⍝ Use Python interoperability for eigenvalues
⎕PY.Import 'numpy as np'
eigs ← {⍵ ⎕PY.Call 'np.linalg.eigvals'}

A ← 2 2 ⍴ 1.2 0.3 ¯0.1 0.85   ⍝ 2×2 transition matrix

eig_vals ← eigs A
moduli   ← |eig_vals            ⍝ absolute values
n_outside← +/ moduli > 1        ⍝ count eigenvalues outside unit circle
n_outside  ⍝ should equal number of free variables for determinacy
```

---

## 4.7 Worked Example: Solving the Two-Equation NK System

*Cross-reference: Principles Ch. 10, Ch. 23* **[P:Ch.Ch.10, P:Ch.23]**

Consider the simplified NK model with only a cost-push shock $u_t$:
$$\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t + u_t, \quad u_t = \rho_u u_{t-1} + \varepsilon_t$$
$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma\phi_\pi\hat{\pi}_t - \sigma\phi_y\hat{x}_t.$$

**Calibration:** $\beta = 0.99$, $\kappa = 0.15$, $\sigma = 1$, $\phi_\pi = 1.5$, $\phi_y = 0.5$, $\rho_u = 0.7$.

**Step 1: Canonical form.** Stack $\mathbf{y}_t = (\hat{\pi}_t, \hat{x}_t, u_t)'$ (the shock is a state variable since it is predetermined). Write:

$$\Gamma_0\mathbf{y}_t = \Gamma_1\mathbb{E}_t[\mathbf{y}_{t+1}] + \text{noise},$$

but since $u_t$ is predetermined, it belongs in the state vector rather than as a free variable.

**Step 2: MSV guess.** Since the only state is $u_t$, guess:
$$\hat{\pi}_t = \omega_\pi u_t, \quad \hat{x}_t = \omega_x u_t.$$

**Step 3: Solve for coefficients.** Substituting $\hat{\pi}_t = \omega_\pi u_t$, $\hat{x}_t = \omega_x u_t$, and $\mathbb{E}_t[\hat{\pi}_{t+1}] = \omega_\pi\rho_u u_t$, $\mathbb{E}_t[\hat{x}_{t+1}] = \omega_x\rho_u u_t$ into the NKPC:

$$\omega_\pi u_t = \beta\omega_\pi\rho_u u_t + \kappa\omega_x u_t + u_t \implies \omega_\pi(1 - \beta\rho_u) = \kappa\omega_x + 1.$$

And into the DIS (with Taylor rule substituted):

$$\omega_x u_t = \omega_x\rho_u u_t - \sigma\phi_\pi\omega_\pi u_t - \sigma\phi_y\omega_x u_t \implies \omega_x(1 - \rho_u + \sigma\phi_y) = -\sigma\phi_\pi\omega_\pi.$$

**Step 4: Solve the 2×2 system.** Let $a = 1 - \beta\rho_u$, $b = 1 - \rho_u + \sigma\phi_y$, $c = \sigma\phi_\pi$:

$$\omega_\pi = \frac{\kappa\omega_x + 1}{a}, \quad \omega_x = \frac{-c\omega_\pi}{b}.$$

Substituting the second into the first:

$$\omega_\pi = \frac{-c\kappa\omega_\pi/b + 1}{a} \implies \omega_\pi\left(a + \frac{c\kappa}{b}\right) = 1 \implies \omega_\pi = \frac{b}{ab + c\kappa}.$$

Then $\omega_x = -c\omega_\pi/b = -c/(ab + c\kappa)$.

**Step 5: Numerical evaluation.**

$a = 1 - 0.99 \times 0.7 = 0.307$

$b = 1 - 0.7 + 1 \times 0.5 = 0.8$

$c = 1 \times 1.5 = 1.5$

$ab + c\kappa = 0.307 \times 0.8 + 1.5 \times 0.15 = 0.2456 + 0.225 = 0.4706$

$\omega_\pi = 0.8/0.4706 = 1.700$

$\omega_x = -1.5/0.4706 = -3.188$

**Interpretation:** A positive cost-push shock of 1 unit raises inflation by 1.700 and lowers the output gap by 3.188. The large output gap response reflects the central bank's aggressive tightening (high $\phi_\pi = 1.5$) in response to the inflation increase.

```apl
⍝ APL — MSV solution for 2-variable NK system
⎕IO←0 ⋄ ⎕ML←1

beta←0.99  ⋄  kappa←0.15  ⋄  sigma←1
phi_pi←1.5  ⋄  phi_y←0.5  ⋄  rho_u←0.7

a ← 1 - beta×rho_u
b ← 1 - rho_u + sigma×phi_y
c ← sigma×phi_pi

denom ← (a×b) + c×kappa
omega_pi ← b÷denom
omega_x  ← (-c)÷denom

omega_pi  ⍝ ≈ 1.700
omega_x   ⍝ ≈ ¯3.188

⍝ Impulse response: response to unit shock at t=0
T ← 20
irf_pi ← omega_pi × rho_u * ⍳T   ⍝ omega_pi × rho_u^t
irf_x  ← omega_x  × rho_u * ⍳T

irf_pi    ⍝ inflation response path
irf_x     ⍝ output gap response path
```

---

## 4.8 The Cobweb Model and Expectational Stability

The classic **cobweb model** illustrates how adaptive versus rational expectations change the dynamics of a market with production lags.

Supply: $Q_t^s = a + b P_t^e$ (firms produce based on expected price)
Demand: $Q_t^d = c - d P_t$ (consumers respond to actual price)

Market clearing: $Q_t^s = Q_t^d$, so:
$$P_t = \frac{c - a}{d} - \frac{b}{d}P_t^e.$$

Under **adaptive expectations** $P_t^e = P_{t-1}$ (naive: expect last period's price):
$$P_t = \frac{c-a}{d} - \frac{b}{d}P_{t-1}.$$

This is a first-order difference equation with coefficient $-b/d$. Stability requires $|{-b/d}| < 1$, i.e., supply is more elastic than demand ($b < d$). If supply is more elastic ($b > d$), the cobweb explodes.

Under **rational expectations** $P_t^e = \mathbb{E}[P_t]$ (firms correctly anticipate the equilibrium price), the market always clears at the rational expectations equilibrium price $P^* = (c-a)/(b+d)$ — no cobweb dynamics at all. Rational expectations "flatten" the cobweb by eliminating the lag between expectation and realization.

This example, while simple, illustrates the fundamental point of *Principles* Chapter 16 [P:Ch.16]: rational expectations can eliminate systematic forecast errors and with them the price dynamics that drive business cycles. Whether prices and wages are actually set as rationally as this model implies is an empirical question — but it sets the theoretical benchmark.

---

## 4.9 Programming Exercises

### Exercise 4.1 (APL — Multiplier–Accelerator)

Implement the Samuelson multiplier–accelerator dynamics as a scan operation.

```apl
⎕IO←0 ⋄ ⎕ML←1

b ← 0.8    ⍝ MPC
v ← 1.2    ⍝ accelerator
G ← 100    ⍝ government spending

⍝ Recursion: Y_t = b(1+v)*Y_{t-1} - bv*Y_{t-2} + G
step ← {b×(1+v)×⊃⌽⍵) - (b×v×⊃⍵) + G}   ⍝ ⍵ is pair [Y_{t-2}, Y_{t-1}]

⍝ Collect path: start from steady state Y*=G/(1-b) with small perturbation
Yss ← G÷1-b
Y0  ← Yss × 1.05   ⍝ 5% above steady state
path ← step \ 50 ⍴ ⊂ Yss Y0   ⍝ scan: collect all pairs

⍝ Extract output sequence
Y_seq ← {⊃⌽ ⍵}¨ path     ⍝ second element of each pair
Y_seq
```

### Exercise 4.2 (Python — Stability Diagram)

```python
import numpy as np
import matplotlib.pyplot as plt

b_vals = np.linspace(0, 1, 100)
v_vals = np.linspace(0, 3, 100)
B, V = np.meshgrid(b_vals, v_vals)

# Characteristic roots
p = -B*(1+V)
q = B*V
disc = p**2 - 4*q
lambda_modulus = np.where(disc >= 0,
    np.maximum(np.abs((-p + np.sqrt(np.abs(disc)))/2),
               np.abs((-p - np.sqrt(np.abs(disc)))/2)),
    np.sqrt(q))  # for complex roots, modulus = sqrt(q)

stable = lambda_modulus < 1
oscillatory = (disc < 0) & stable

fig, ax = plt.subplots()
ax.contourf(B, V, stable.astype(int), levels=[0.5,1.5], colors=['lightblue'], alpha=0.5)
ax.contour(B, V, oscillatory.astype(int), levels=[0.5], colors=['red'])
ax.set_xlabel('MPC (b)'); ax.set_ylabel('Accelerator (v)')
ax.set_title('Stability diagram: blue=stable, red boundary=oscillatory')
plt.show()
```

### Exercise 4.3 (Julia — MSV Solution)

```julia
# MSV solution for arbitrary NK calibration
function msv_nk(; beta=0.99, kappa=0.15, sigma=1.0,
                  phi_pi=1.5, phi_y=0.5, rho_u=0.7)
    a = 1 - beta*rho_u
    b = 1 - rho_u + sigma*phi_y
    c = sigma*phi_pi
    denom = a*b + c*kappa
    omega_pi = b / denom
    omega_x  = -c / denom
    return (omega_pi=omega_pi, omega_x=omega_x)
end

sol = msv_nk()
println("ω_π = $(round(sol.omega_pi, digits=3))")
println("ω_x = $(round(sol.omega_x, digits=3))")

# Sweep over phi_pi to show determinacy
phi_pis = 0.5:0.1:3.0
for φ in phi_pis
    s = msv_nk(phi_pi=φ)
    println("φ_π = $φ: ω_π = $(round(s.omega_pi,digits=3)), ω_x = $(round(s.omega_x,digits=3))")
end
```

### Exercise 4.4 (R — Forward Solution)

```r
# Forward solution of NKPC: pi_t = kappa * sum_{j>=0} beta^j E_t[x_{t+j}]
# Given AR(1) output gap: x_t = rho_x * x_{t-1} + eps_t
kappa <- 0.15; beta <- 0.99; rho_x <- 0.8

# Closed-form: pi_t = kappa/(1-beta*rho_x) * x_t
omega_pi_forward <- kappa / (1 - beta*rho_x)
cat("π_t = ", round(omega_pi_forward, 4), "* x_t\n")

# Simulate and verify numerically
T <- 1000
set.seed(42)
x <- filter(rnorm(T, sd=0.1), rho_x, method="recursive")
pi_formula <- omega_pi_forward * x

# Compare to direct present-value calculation
pi_pv <- sapply(1:(T-50), function(t) {
  kappa * sum(sapply(0:49, function(j) beta^j * x[t+j]))
})
cat("Max discrepancy:", max(abs(pi_formula[1:(T-50)] - pi_pv)), "\n")
```

### Exercise 4.5 — Blanchard–Kahn ($\star$)

For the 2×2 NK system in Section 4.5.1, write the matrix $A = \Gamma_0^{-1}\Gamma_1$ explicitly and compute its eigenvalues for the calibration in Worked Example 4.7. Verify that both eigenvalues lie outside the unit circle when $\phi_\pi = 1.5$ and $\phi_y = 0.5$, and find the critical value of $\phi_\pi$ below which one eigenvalue moves inside the unit circle (violating determinacy). Show this threshold is consistent with the Taylor principle $\phi_\pi > 1$.

### Exercise 4.6 — Indeterminacy and Sunspots ($\star\star$)

When $\phi_\pi < 1$, the NK model has indeterminate equilibria. The general solution becomes $\mathbf{y}_t = \Omega z_t + \Pi\eta_t$ where $\eta_t$ is an arbitrary martingale difference sequence (the "sunspot"). (a) Explain why this generates volatility unrelated to fundamentals. (b) Simulate the model under $\phi_\pi = 0.8$ with a sunspot $\eta_t \sim \mathcal{N}(0,1)$ and compare the resulting volatility of $\hat{\pi}_t$ and $\hat{x}_t$ to the fundamental-only case $\phi_\pi = 1.5$. (c) Show that the Taylor principle eliminates the sunspot component by making $\Pi = \mathbf{0}$.

---

## 4.10 Chapter Summary

**Key results:**

- First-order linear difference equations $x_{t+1} = \lambda x_t + c$ are stable iff $|\lambda| < 1$, with convergence rate $|\lambda|$ and steady state $x^* = c/(1-\lambda)$.
- **Backward solutions** require $|\lambda| < 1$ (predetermined variables); **forward solutions** require $|\lambda| > 1$ (free/jump variables). The present-value formula structure of the NKPC arises from the forward solution with $|\lambda| = 1/\beta > 1$.
- Second-order equations $x_{t+2} + px_{t+1} + qx_t = c$ generate oscillatory dynamics when the discriminant $p^2 - 4q < 0$; stability requires $\sqrt{q} < 1$.
- The **MSV solution** of an expectational difference equation is found by guessing a linear function of state variables and matching coefficients — the Sylvester equation.
- The **Blanchard–Kahn condition** for a unique bounded rational expectations equilibrium: number of eigenvalues outside the unit circle equals number of free variables.
- In APL: matrix powers `A⍣h` compute impulse responses efficiently; the scan `\` collects path histories; `rho_u * ⍳T` generates geometric decay sequences directly.

**Connections forward:** Chapter 14 revisits difference equations in the multiplier–accelerator context. Chapter 18 develops the full rational expectations solution methodology (undetermined coefficients, the Sims algorithm). Chapter 28 provides the definitive treatment of the Blanchard–Kahn conditions and the QZ decomposition algorithm for DSGE models.

---

*Next: Chapter 5 — Stochastic Processes for Aggregate Shocks*
