# Chapter 3: Differential Equations in Continuous-Time Macro Models

*The Solow and Ramsey Foundations*

> *"The solution of a differential equation is a function — not a number but a whole trajectory through time."*

**Cross-reference:** *Principles* Ch. 5.2 (Solow ODE); Ch. 5.3 (RCK Hamiltonian system and phase plane); Ch. 19 (labour market dynamics); Ch. 3 (GDP convergence) **[P:Ch.5.2, P:Ch.5.3, P:Ch.19, P:Ch.3]**

---

## 3.1 Why Continuous Time? The Case for ODEs in Growth Theory

Two modeling conventions compete in macroeconomics: continuous time (differential equations) and discrete time (difference equations). The choice is largely one of mathematical convenience, but the two approaches have different strengths.

Continuous-time models are analytically cleaner for growth theory and optimal control. The differential equation $\dot{k} = sf(k) - \delta k$ is a single object with a well-developed theory. Phase plane analysis — drawing the $(k, c)$ plane, finding nullclines, identifying stable manifolds — provides geometric intuition that is harder to obtain in discrete time. The Pontryagin maximum principle, developed in Chapter 11, is most naturally stated in continuous time.

Discrete-time models are preferred for business cycle analysis, DSGE estimation, and numerical work. All macroeconomic data are observed at discrete intervals (quarterly GDP, monthly inflation), and most numerical algorithms — value function iteration, Kalman filtering, particle filtering — are discrete-time procedures.

This chapter develops continuous-time differential equations systematically, with every result illustrated by the two key models of *Principles* Part I: the Solow growth model and the RCK optimal growth model. Chapter 4 develops the discrete-time analogue.

---

## 3.2 First-Order Linear Ordinary Differential Equations

**Definition 3.1 (Ordinary Differential Equation).** An **ordinary differential equation (ODE)** of order $n$ is an equation relating a function $x(t)$ and its derivatives up to order $n$:
$$F(t, x, \dot{x}, \ddots x, \ldots, x^{(n)}) = 0.$$
It is **ordinary** (as opposed to partial) because $x$ depends on a single variable $t$.

We begin with the first-order linear case, which covers the linearized Solow dynamics and the Wicksellian adjustment process.

### 3.2.1 Homogeneous Equation

The equation $\dot{x}(t) = ax(t)$ has the general solution:

$$x(t) = x_0 e^{at},$$

where $x_0 = x(0)$ is the initial condition. The sign of $a$ determines behavior:
- $a < 0$: $x(t) \to 0$ (stable, exponential decay).
- $a = 0$: $x(t) = x_0$ (constant).
- $a > 0$: $x(t) \to \pm\infty$ (unstable, exponential growth).

**Theorem 3.1 (Stability of Linear First-Order ODE).** The equilibrium $x^* = 0$ of $\dot{x} = ax$ is **globally asymptotically stable** if and only if $a < 0$.

This is the continuous-time analogue of the discrete-time condition $|\lambda| < 1$. The connection is: if we discretize $\dot{x} = ax$ with step size $h$, we get $x_{t+1} = (1 + ah)x_t$, which is stable iff $|1 + ah| < 1$, iff $ah \in (-2, 0)$ — which for small $h$ requires $a < 0$.

### 3.2.2 Nonhomogeneous Equation and Particular Solutions

The equation $\dot{x} = ax + b$ (with $a \neq 0$) has equilibrium $x^* = -b/a$. Let $\tilde{x} = x - x^*$; then $\dot{\tilde{x}} = a\tilde{x}$, which is the homogeneous equation. Therefore:

$$x(t) = x^* + (x_0 - x^*)e^{at} = -\frac{b}{a} + \left(x_0 + \frac{b}{a}\right)e^{at}.$$

For $a < 0$, $x(t) \to x^*$ as $t \to \infty$: the system converges to the steady state at rate $|a|$. The **half-life** — time for the gap $x(t) - x^*$ to halve — is $\tau_{1/2} = \ln 2 / |a|$.

*Application:* The linearized Solow model near the steady state $\tilde{k}^*$ gives:
$$\dot{\tilde{k}} - \dot{\tilde{k}}^* \approx -\lambda(\tilde{k} - \tilde{k}^*),$$
where $\lambda = (1-\alpha)\mu > 0$ [P:Ch.5.2]. This is exactly $\dot{\tilde{x}} = -\lambda\tilde{x}$ with $a = -\lambda < 0$, giving convergence at rate $\lambda$ and half-life $\ln 2/\lambda \approx 17$ years for $\lambda \approx 0.04$ — the conditional convergence result of *Principles* Ch. 5.2.

### 3.2.3 The Integrating Factor Method

For the general first-order linear ODE $\dot{x} + p(t)x = q(t)$, the solution uses the **integrating factor** $\mu(t) = e^{\int p(t)\,dt}$:

$$x(t) = \frac{1}{\mu(t)}\left[C + \int_0^t \mu(s)q(s)\,ds\right].$$

This formula is useful when $p(t)$ or $q(t)$ are time-varying — as in models with time-varying interest rates or technology growth rates.

---

## 3.3 First-Order Nonlinear ODEs: Phase Line Analysis

Most continuous-time macro models generate nonlinear ODEs. Analytical solutions often do not exist, but the **phase line** (for scalar equations) and **phase plane** (for 2D systems) provide a complete qualitative picture of the dynamics.

### 3.3.1 Autonomous Equations and the Phase Line

An **autonomous** ODE $\dot{x} = f(x)$ does not depend explicitly on $t$. The dynamics are fully determined by the sign of $f(x)$:

**Definition 3.2 (Steady State).** A point $x^*$ is a **steady state** (or **equilibrium**) of $\dot{x} = f(x)$ if $f(x^*) = 0$. The system is stationary at $x^*$.

**Definition 3.3 (Phase Line).** The **phase line** is a one-dimensional diagram plotting $f(x)$ on the vertical axis against $x$ on the horizontal axis. The steady states are the zeros of $f$. The system moves rightward ($\dot{x} > 0$) wherever $f(x) > 0$ and leftward wherever $f(x) < 0$.

**Theorem 3.2 (Local Stability via Linearization).** Let $x^*$ be a steady state of $\dot{x} = f(x)$ and suppose $f$ is differentiable at $x^*$. Then $x^*$ is **locally asymptotically stable** if $f'(x^*) < 0$ and **locally unstable** if $f'(x^*) > 0$.

*Proof.* Let $\tilde{x} = x - x^*$. Taylor-expanding: $\dot{\tilde{x}} = f(x^* + \tilde{x}) \approx f(x^*) + f'(x^*)\tilde{x} = f'(x^*)\tilde{x}$. By Theorem 3.1, this converges iff $f'(x^*) < 0$. $\square$

### 3.3.2 The Solow Equation as a Nonlinear ODE

The fundamental equation of the Solow model [P:Ch.5.2]:
$$\dot{\tilde{k}} = sf(\tilde{k}) - \mu\tilde{k},$$
where $f(\tilde{k}) = \tilde{k}^\alpha$ (Cobb–Douglas), $s$ is the saving rate, $\mu = n + g + \delta$ is the effective depreciation rate.

**Phase line analysis:**

1. **Steady state:** $s\tilde{k}^{*\alpha} = \mu\tilde{k}^*$, giving $\tilde{k}^* = (s/\mu)^{1/(1-\alpha)}$.
2. **For $\tilde{k} < \tilde{k}^*$:** $sf(\tilde{k})/\tilde{k} > s\tilde{k}^{*\alpha}/\tilde{k}^* = \mu$ (since $f(\tilde{k})/\tilde{k}$ is decreasing in $\tilde{k}$), so $\dot{\tilde{k}} > 0$.
3. **For $\tilde{k} > \tilde{k}^*$:** $sf(\tilde{k})/\tilde{k} < \mu$, so $\dot{\tilde{k}} < 0$.

Therefore $\tilde{k}^*$ is globally stable: the economy converges to it from any initial $\tilde{k}_0 > 0$.

**Local convergence rate:**
$$f'(\tilde{k}^*) = \frac{d}{d\tilde{k}}[sf(\tilde{k}) - \mu\tilde{k}]\Big|_{\tilde{k}^*} = sf'(\tilde{k}^*) - \mu = \alpha s\tilde{k}^{*\alpha-1} - \mu = \alpha\mu - \mu = -\mu(1-\alpha) = -\lambda < 0,$$
confirming stability and giving the convergence rate $\lambda = \mu(1-\alpha)$.

**The Cobb–Douglas analytical solution.** For $f(\tilde{k}) = \tilde{k}^\alpha$, the Solow ODE is a **Bernoulli equation**, solvable in closed form. Let $z = \tilde{k}^{1-\alpha}$; then $\dot{z} = (1-\alpha)\tilde{k}^{-\alpha}\dot{\tilde{k}} = (1-\alpha)(s - \mu z)$, which is a linear ODE with solution:

$$z(t) = z^* + (z_0 - z^*)e^{-\lambda t}, \quad z^* = s/\mu.$$

Converting back: $\tilde{k}(t) = [\tilde{k}^{*1-\alpha} + (\tilde{k}_0^{1-\alpha} - \tilde{k}^{*1-\alpha})e^{-\lambda t}]^{1/(1-\alpha)}$.

This is the **exact transition path** of the Solow model — no approximation. It is used in Chapter 10 for numerical verification and to replicate the Mankiw–Romer–Weil convergence regressions.

---

## 3.4 Systems of ODEs and Phase Plane Analysis

The RCK model involves two coupled ODEs — the Euler equation for consumption and the capital accumulation equation — which together constitute a 2D autonomous system. Understanding 2D systems requires phase plane analysis.

### 3.4.1 Linear 2D Systems

For the system $\dot{\mathbf{x}} = A\mathbf{x}$ with $A \in \mathbb{R}^{2\times 2}$, the solution is:

$$\mathbf{x}(t) = e^{At}\mathbf{x}_0 = P e^{Dt} P^{-1}\mathbf{x}_0,$$

where $A = PDP^{-1}$ is the eigendecomposition (assuming $A$ is diagonalisable). The qualitative behavior near the origin is completely determined by the eigenvalues $\lambda_1, \lambda_2$.

**Definition 3.4 (Classification of 2D Equilibria).** For $\dot{\mathbf{x}} = A\mathbf{x}$ with eigenvalues $\lambda_1, \lambda_2$:

| Eigenvalues | Equilibrium Type | Behavior |
|---|---|---|
| $\lambda_1 < \lambda_2 < 0$ | Stable node | Monotone convergence |
| $0 < \lambda_1 < \lambda_2$ | Unstable node | Monotone divergence |
| $\lambda_1 < 0 < \lambda_2$ | **Saddle point** | Stable in one direction, unstable in other |
| $\text{Re}(\lambda) < 0$, $\text{Im}(\lambda) \neq 0$ | Stable spiral | Oscillatory convergence |
| $\text{Re}(\lambda) > 0$, $\text{Im}(\lambda) \neq 0$ | Unstable spiral | Oscillatory divergence |
| $\lambda_{1,2} = \pm bi$, purely imaginary | Centre | Closed orbits (neutral stability) |

The **saddle point** is the most important case in macroeconomics. The RCK model has a saddle point: there is a one-dimensional stable manifold (the **saddle path**) along which the system converges to the steady state, and convergence occurs along this path. Any trajectory starting off the saddle path either violates the transversality condition (too much capital accumulation) or runs into a corner (negative consumption). The unique optimal trajectory is thus the saddle path, selected by the appropriate initial condition on the jump variable $c_0$.

### 3.4.2 The Nullcline Method

**Definition 3.5 (Nullclines).** For a 2D autonomous system $\dot{x}_1 = f_1(x_1, x_2)$, $\dot{x}_2 = f_2(x_1, x_2)$:
- The **$\dot{x}_1 = 0$ nullcline** is the curve $\{(x_1, x_2): f_1(x_1, x_2) = 0\}$.
- The **$\dot{x}_2 = 0$ nullcline** is the curve $\{(x_1, x_2): f_2(x_1, x_2) = 0\}$.

Steady states occur at nullcline intersections. The phase plane is divided into four regions by the nullclines; in each region, the signs of $\dot{x}_1$ and $\dot{x}_2$ are constant, determining the direction of motion with arrows.

**Example — RCK Phase Plane.** The system in effective-labor units $(\tilde{k}, \tilde{c})$:

$$\dot{\tilde{k}} = f(\tilde{k}) - \tilde{c} - \mu\tilde{k}, \qquad \frac{\dot{\tilde{c}}}{\tilde{c}} = \frac{1}{\sigma}[f'(\tilde{k}) - \delta - \rho - \sigma g].$$

**$\dot{\tilde{k}} = 0$ nullcline:** $\tilde{c} = f(\tilde{k}) - \mu\tilde{k}$. This is a hump-shaped curve, peaking at the Golden Rule capital stock $\tilde{k}^{GR}$ where $f'(\tilde{k}^{GR}) = \mu$.

**$\dot{\tilde{c}} = 0$ nullcline:** $f'(\tilde{k}^*) = \delta + \rho + \sigma g$. This is a vertical line at the unique $\tilde{k}^*$.

The steady state $(\tilde{k}^*, \tilde{c}^*)$ is at the intersection. The RCK steady state satisfies $\tilde{k}^* < \tilde{k}^{GR}$ because impatience ($\rho > 0$) pushes $f'(\tilde{k}^*) = \delta + \rho + \sigma g > \delta + n + g = f'(\tilde{k}^{GR})$, implying $\tilde{k}^* < \tilde{k}^{GR}$ since $f'' < 0$.

### 3.4.3 Linearization Around the Steady State

To analyze local dynamics, linearize the system around $(\tilde{k}^*, \tilde{c}^*)$. Define deviations $x_1 = \tilde{k} - \tilde{k}^*$, $x_2 = \tilde{c} - \tilde{c}^*$. The Jacobian:

$$J = \begin{pmatrix} \partial\dot{\tilde{k}}/\partial\tilde{k} & \partial\dot{\tilde{k}}/\partial\tilde{c} \\ \partial\dot{\tilde{c}}/\partial\tilde{k} & \partial\dot{\tilde{c}}/\partial\tilde{c} \end{pmatrix}\Bigg|_{(\tilde{k}^*, \tilde{c}^*)} = \begin{pmatrix} f'(\tilde{k}^*) - \mu & -1 \\ \frac{\tilde{c}^*}{\sigma}f''(\tilde{k}^*) & 0 \end{pmatrix}.$$

**Theorem 3.3 (Saddle-Point Property of the RCK Steady State).** The Jacobian $J$ of the RCK system at the steady state has one negative and one positive real eigenvalue (i.e., the steady state is a saddle point).

*Proof.* $\det(J) = 0\cdot(f'(\tilde{k}^*) - \mu) - (-1)\cdot\frac{\tilde{c}^*}{\sigma}f''(\tilde{k}^*) = \frac{\tilde{c}^*}{\sigma}f''(\tilde{k}^*) < 0$, since $f'' < 0$ (diminishing returns) and $\tilde{c}^* > 0$. Since $\det(J) = \lambda_1\lambda_2 < 0$, the eigenvalues have opposite signs, confirming a saddle point. $\square$

The product of eigenvalues is negative: one is positive, one is negative. The negative eigenvalue $\lambda_-$ governs the speed of convergence along the saddle path. From the quadratic formula applied to $\lambda^2 - \text{tr}(J)\lambda + \det(J) = 0$:

$$\lambda_\pm = \frac{\text{tr}(J) \pm \sqrt{\text{tr}(J)^2 - 4\det(J)}}{2}.$$

Since $\det(J) < 0$, the discriminant $\text{tr}(J)^2 - 4\det(J) > \text{tr}(J)^2 > 0$, ensuring real eigenvalues.

---

## 3.5 Second-Order ODEs: The Multiplier–Accelerator

A second-order ODE $\ddot{x} + p\dot{x} + qx = r$ arises in the continuous-time multiplier–accelerator model [P:Ch.27]. It can always be converted to a 2D first-order system by letting $x_1 = x$, $x_2 = \dot{x}$:

$$\dot{x}_1 = x_2, \qquad \dot{x}_2 = -qx_1 - px_2 + r.$$

The matrix of this system is $A = \begin{pmatrix} 0 & 1 \\ -q & -p \end{pmatrix}$ with $\text{tr}(A) = -p$ and $\det(A) = q$.

The characteristic equation of the second-order ODE: $\lambda^2 + p\lambda + q = 0$. The roots:
$$\lambda = \frac{-p \pm \sqrt{p^2 - 4q}}{2}.$$

**Stability (both roots with negative real part) requires $p > 0$ and $q > 0$.**

The nature of the roots determines whether approach to equilibrium is monotone ($p^2 > 4q$, real roots) or oscillatory ($p^2 < 4q$, complex roots with non-zero imaginary part). Oscillatory convergence corresponds to the business-cycle-like fluctuations in the continuous-time multiplier–accelerator model.

---

## 3.6 Numerical Solution of ODEs

When analytical solutions are unavailable — which is the rule for nonlinear systems — numerical methods are essential. The two most important are the Euler method and the Runge–Kutta 4th-order method.

**Algorithm 3.1 (Euler Method).**

Given $\dot{x} = f(t, x)$, $x(0) = x_0$, step size $h$:
1. Set $t_0 = 0$, $x_0 = x(0)$.
2. For $n = 0, 1, \ldots, N-1$: $x_{n+1} = x_n + h \cdot f(t_n, x_n)$, $t_{n+1} = t_n + h$.

The Euler method has local truncation error $O(h^2)$ and global error $O(h)$. It is simple but requires small $h$ for accuracy.

**Algorithm 3.2 (Runge–Kutta 4th Order, RK4).**

Given $\dot{x} = f(t, x)$, $x(0) = x_0$, step size $h$:
1. Compute four slopes:
   - $k_1 = f(t_n, x_n)$
   - $k_2 = f(t_n + h/2,\; x_n + (h/2)k_1)$
   - $k_3 = f(t_n + h/2,\; x_n + (h/2)k_2)$
   - $k_4 = f(t_n + h,\; x_n + h\,k_3)$
2. Update: $x_{n+1} = x_n + (h/6)(k_1 + 2k_2 + 2k_3 + k_4)$.

RK4 has local truncation error $O(h^5)$ and global error $O(h^4)$ — dramatically more accurate than Euler for the same step size.

```apl
⍝ APL — RK4 for the Solow model
⎕IO←0 ⋄ ⎕ML←1

⍝ Parameters
alpha ← 1÷3  ⋄  s ← 0.25  ⋄  mu ← 0.08   ⍝ n+g+delta

⍝ Solow ODE: dk̃/dt = s*k̃^α - μ*k̃
solow ← {s × ⍵*alpha) - mu × ⍵}   ⍝ scalar ODE, ⍵ = k̃

⍝ One RK4 step: f is the ODE, h is step size, x is current value
rk4step ← {f h x ← ⍺⍺ ⍺ ⍵
    k1 ← f x
    k2 ← f x + (h÷2)×k1
    k3 ← f x + (h÷2)×k2
    k4 ← f x + h×k3
    x + (h÷6) × k1 + 2×k2 + 2×k3 + k4}

⍝ Solve over T periods with step h
h ← 0.1  ⋄  T ← 200  ⋄  k0 ← 0.5   ⍝ initial capital below steady state

⍝ Generate trajectory: apply rk4step T times
steps ← ⍳ T
trajectory ← {(solow rk4step h) ⍵}⍣T ⊢ k0    ⍝ iterate T times

⍝ Or collect all intermediate values using scan
all_k ← {(solow rk4step h) ⍵} \ T ⍴ k0    ⍝ collect path
```

```python
# Python — RK4 for the Solow model
import numpy as np
import matplotlib.pyplot as plt

alpha, s, mu = 1/3, 0.25, 0.08

def solow(k): return s * k**alpha - mu * k

def rk4_step(f, x, h):
    k1 = f(x)
    k2 = f(x + (h/2)*k1)
    k3 = f(x + (h/2)*k2)
    k4 = f(x + h*k3)
    return x + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

h, T, k0 = 0.1, 200, 0.5
k_ss = (s / mu) ** (1/(1-alpha))
print(f"Steady state: k* = {k_ss:.4f}")

path = [k0]
for _ in range(T):
    path.append(rk4_step(solow, path[-1], h))

plt.plot(np.arange(T+1)*h, path, label='RK4 trajectory')
plt.axhline(k_ss, linestyle='--', color='red', label='Steady state k*')
plt.xlabel('Time'); plt.ylabel('k̃'); plt.legend(); plt.show()
```

```julia
# Julia — RK4 for the Solow model
alpha, s, mu = 1/3, 0.25, 0.08
solow(k) = s * k^alpha - mu * k

function rk4_step(f, x, h)
    k1 = f(x); k2 = f(x + (h/2)*k1)
    k3 = f(x + (h/2)*k2); k4 = f(x + h*k3)
    x + (h/6)*(k1 + 2k2 + 2k3 + k4)
end

h, T, k0 = 0.1, 200, 0.5
path = accumulate((k,_) -> rk4_step(solow, k, h), 1:T; init=k0)
println("Final k: ", round(last(path), digits=4), "  Steady state: ", round((s/mu)^(1/(1-alpha)), digits=4))
```

```r
# R — RK4 for the Solow model
alpha <- 1/3; s <- 0.25; mu <- 0.08
solow <- function(k) s * k^alpha - mu * k
rk4_step <- function(f, x, h) {
  k1 <- f(x); k2 <- f(x + (h/2)*k1)
  k3 <- f(x + (h/2)*k2); k4 <- f(x + h*k3)
  x + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
}
h <- 0.1; T <- 200; k0 <- 0.5
path <- Reduce(function(k, .) rk4_step(solow, k, h), seq_len(T), accumulate=TRUE, init=k0)
plot(seq(0, T*h, h), path, type='l', xlab='Time', ylab='k̃',
     main='Solow convergence'); abline(h=(s/mu)^(1/(1-alpha)), lty=2, col='red')
```

---

## 3.7 Worked Example: Full Phase Diagram of the Solow Model

*Cross-reference: Principles Ch. 5.2* **[P:Ch.5.2]**

Calibration: $\alpha = 1/3$, $s = 0.20$, $\delta = 0.05$, $n = 0.01$, $g = 0.02$, so $\mu = 0.08$.

**Steady state:** $\tilde{k}^* = (0.20/0.08)^{3/2} = (2.5)^{1.5} = 3.953$.

**Convergence rate:** $\lambda = (1-\alpha)\mu = (2/3)(0.08) = 0.0533$. Half-life: $\ln 2/0.0533 \approx 13.0$ years.

**Phase line:** At $\tilde{k} = 1$ (below steady state): $\dot{\tilde{k}} = 0.20(1)^{1/3} - 0.08(1) = 0.12 > 0$. At $\tilde{k} = 6$ (above steady state): $\dot{\tilde{k}} = 0.20(6)^{1/3} - 0.08(6) = 0.362 - 0.480 = -0.118 < 0$. Both confirm convergence to $\tilde{k}^*$.

**Transition path (analytical):** With $z = \tilde{k}^{2/3}$ and $z^* = (s/\mu) = 2.5$:
$$\tilde{k}(t) = [z^* + (z_0 - z^*)e^{-\lambda t}]^{3/2} = [2.5 + (z_0 - 2.5)e^{-0.0533t}]^{3/2}.$$

For $\tilde{k}_0 = 1$: $z_0 = 1^{2/3} = 1$, so $\tilde{k}(t) = [2.5 - 1.5e^{-0.0533t}]^{3/2}$.

At $t = 13$ years (one half-life): gap should halve. Current gap: $z_0 - z^* = -1.5$. Predicted $z(13) = 2.5 - 1.5e^{-0.693} = 2.5 - 0.75 = 1.75$, giving $\tilde{k}(13) = (1.75)^{1.5} = 2.315$. The gap $\tilde{k}^* - \tilde{k}(13) = 3.953 - 2.315 = 1.638$, which is indeed half of the initial gap $3.953 - 1 = 2.953$. ✓

---

## 3.8 Programming Exercises

### Exercise 3.1 (APL)

Write a dfn `solow_phase ← {s delta n g alpha ← ⍵ ⋄ ...}` that takes parameters as a 5-vector and returns the steady state $\tilde{k}^*$, convergence rate $\lambda$, and half-life $\tau_{1/2}$. Test with the calibration from Section 3.7.

### Exercise 3.2 (Python — Phase Diagram)

Plot the Solow phase diagram: the curve $sf(\tilde{k})$ and the line $\mu\tilde{k}$ on the same axes for $\tilde{k} \in [0, 8]$, with the steady state marked. Add three trajectory arrows at $\tilde{k} \in \{0.5, 2, 6\}$ showing the direction of motion.

### Exercise 3.3 (Julia — RCK Saddle Path)

Implement the **reverse-shooting algorithm** for the RCK model: start from the steady state $(\tilde{k}^*, \tilde{c}^*)$, perturb slightly along the eigenvector corresponding to $\lambda_-$, and integrate backward in time to trace the saddle path. Parameters: $\alpha = 1/3$, $\rho = 0.04$, $\sigma = 2$, $\delta = 0.05$, $n = 0.01$, $g = 0.02$.

### Exercise 3.4 (R — Convergence Regression)

Using the analytical Solow transition path $\tilde{k}(t) = [z^* + (z_0 - z^*)e^{-\lambda t}]^{1/(1-\alpha)}$, generate simulated data for 100 "countries" with:
- $z_0 \sim U[0.5, 3]$ (different initial conditions)
- $s \sim U[0.1, 0.4]$ (different saving rates)
- $\mu = 0.08$ (common effective depreciation)

Compute 40-year average growth rates and regress on $\ln\tilde{k}_{i,0}$ with and without $\ln(s_i/\mu)$ as a control. Verify that conditional convergence ($\beta < 0$ conditional on the saving rate) holds but unconditional convergence does not.

### Exercise 3.5 — Stability Classification ($\star$)

For the linear system $\dot{\mathbf{x}} = J\mathbf{x}$ with Jacobian $J$ from Section 3.4.3:

(a) Express $\text{tr}(J)$ and $\det(J)$ in terms of $f'(\tilde{k}^*)$, $f''(\tilde{k}^*)$, $\tilde{c}^*$, $\sigma$, $\mu$.

(b) Show that $\det(J) < 0$ for any $f$ satisfying $f'' < 0$ and any $\tilde{c}^* > 0$, confirming the saddle-point property without reference to any specific functional form.

(c) In APL, implement a function that takes $(f', f'', \tilde{c}^*, \sigma, \mu)$ as arguments and returns both eigenvalues of $J$, verifying they have opposite signs.

### Exercise 3.6 — Endogenous Growth ($\star\star$)

The AK model [P:Ch.5.4] has $f(\tilde{k}) = A\tilde{k}$ (no diminishing returns). The growth equation becomes $\dot{K}/K = sA - \delta$. Show that: (a) there is no interior steady state of the Solow type; (b) all trajectories grow at rate $g^* = sA - \delta$ (for $sA > \delta$); (c) the phase line technique fails here because $f'(k) = A$ is constant. Explain what this implies for the convergence analysis.

---

## 3.9 Chapter Summary

This chapter developed the theory of ordinary differential equations for continuous-time macroeconomic models.

**Key results:**

- The first-order linear ODE $\dot{x} = ax + b$ has solution $x(t) = x^* + (x_0 - x^*)e^{at}$ with steady state $x^* = -b/a$; stable iff $a < 0$.
- The **linearization theorem** allows local stability analysis of nonlinear ODEs: the stability of $x^*$ for $\dot{x} = f(x)$ is determined by the sign of $f'(x^*)$.
- The **Solow ODE** is a Bernoulli equation with the exact solution $\tilde{k}(t) = [z^* + (z_0 - z^*)e^{-\lambda t}]^{1/(1-\alpha)}$, giving convergence at rate $\lambda = (1-\alpha)\mu$.
- For 2D systems $\dot{\mathbf{x}} = A\mathbf{x}$, the equilibrium type (node, saddle, spiral, center) is determined by $(\text{tr}(A), \det(A))$; the RCK steady state is a saddle point because $\det(J) < 0$.
- **RK4** is the standard numerical ODE solver, with global error $O(h^4)$ far superior to Euler's $O(h)$.

**Connections forward:** Chapter 10 develops the full Solow model analytically and numerically, including the closed-form transition path and the MRW convergence regression. Chapter 11 applies the phase plane analysis to the RCK model, using the saddle-path structure to motivate the shooting algorithm and the transversality condition.

---

*Next: Chapter 4 — Difference Equations in Discrete-Time Macro Models*
