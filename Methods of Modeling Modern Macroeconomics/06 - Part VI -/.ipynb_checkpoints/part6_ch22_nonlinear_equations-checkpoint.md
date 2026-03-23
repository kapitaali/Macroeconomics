# Part VI: Numerical Methods for Solving Macroeconomic Models

*Connects to: All preceding analytical results*

---

When analytical solutions exist, we use them. The Solow model has a closed-form transition path; the log-utility Diamond OLG has an algebraic steady state; the scalar NK MSV solution is a single formula. But these are the exceptions. In medium-scale DSGE models with ten or more state variables, financial frictions, occasionally-binding constraints, and heterogeneous agents, closed-form solutions are impossible. Numerical methods are not a fallback — they are the primary tool.

This part covers the five fundamental numerical techniques that underlie every computational macroeconomics software package, from Dynare to QuantEcon. **Chapter 22** develops root-finding: bisection, Newton–Raphson (with a proof of quadratic convergence), and quasi-Newton methods for finding DSGE steady states. **Chapter 23** covers numerical integration — trapezoidal rules, Gaussian quadrature, and Monte Carlo — for the expected-utility integrals that appear in dynamic programming and option pricing. **Chapter 24** develops numerical optimization, from gradient descent through BFGS, applied to the central bank's optimal policy problem. **Chapter 25** addresses linear system solution: LU decomposition, sparse matrix methods, and QR — the linear algebra backbone of the Kalman filter (Chapter 20), the DSGE solution (Chapter 28), and the Leontief model (Chapter 2). **Chapter 26** synthesizes everything into Monte Carlo simulation of macroeconomic models: the Tauchen discretization, Markov chain simulation, bootstrap confidence bands, and the particle filter for nonlinear state-space models.

The mathematical depth of this part is intentionally applied: every algorithm is derived from first principles (not just stated), and every derivation connects to a specific step in the macroeconomic workflow. The reader who works through these five chapters will understand not just how to call a numerical solver but what it is doing inside — and when it will fail.

---

# Chapter 22: Solving Nonlinear Equations

*Newton–Raphson for Steady-State Calculations*

> *"Newton's method is the most important algorithm in numerical analysis. Everything else is a variation on it."*

**Cross-reference:** *Principles* Ch. 5.2 (Solow steady state as a root); Ch. 5.3 (RCK modified Golden Rule); Ch. 23 (Taylor rule fixed point) **[P:Ch.5.2, P:Ch.5.3, P:Ch.23]**

---

## 22.1 The Steady State as a Root-Finding Problem

Every dynamic model has a **steady state** — a fixed point of the system's law of motion where all variables remain constant over time. Finding the steady state is the first step in any DSGE analysis: the log-linearization (Chapter 27) is performed around the steady state, and the Blanchard–Kahn solution (Chapter 28) requires knowing it precisely.

In simple models, steady states have closed-form solutions. The Solow model gives $\tilde{k}^* = (s/\mu)^{1/(1-\alpha)}$; the log-utility Diamond OLG gives a direct algebraic expression. But in medium-scale DSGE models with multiple nonlinear equilibrium conditions — price-setting optimization, labor market search-and-matching, financial frictions — no closed form exists. We must solve a system of nonlinear equations:

$$\mathbf{F}(\mathbf{x}) = \mathbf{0}, \quad \mathbf{F}: \mathbb{R}^n \to \mathbb{R}^n.$$

The methods of this chapter — bisection, Newton–Raphson, quasi-Newton — address this problem with different trade-offs between reliability and speed.

---

## 22.2 The Bisection Method

**Algorithm 22.1 (Bisection).**

*Input:* Continuous function $f: \mathbb{R} \to \mathbb{R}$; bracket $[a, b]$ with $f(a)f(b) < 0$; tolerance $\varepsilon > 0$.

1. Set $a_0 = a$, $b_0 = b$.
2. For $n = 0, 1, 2, \ldots$:
   - $m_n = (a_n + b_n)/2$ (midpoint).
   - If $|b_n - a_n| < 2\varepsilon$ or $|f(m_n)| < \varepsilon$: **return** $m_n$.
   - If $f(a_n)f(m_n) < 0$: set $a_{n+1} = a_n$, $b_{n+1} = m_n$.
   - Else: set $a_{n+1} = m_n$, $b_{n+1} = b_n$.

**Theorem 22.1 (Bisection Convergence).** The bisection method converges to a root $x^*$ with:

$$|m_n - x^*| \leq \frac{b-a}{2^{n+1}}.$$

The number of iterations to achieve tolerance $\varepsilon$ is $n \geq \lceil\log_2((b-a)/\varepsilon)\rceil$.

*Proof.* At each step the bracket $[a_n, b_n]$ halves in length: $b_n - a_n = (b-a)/2^n$. The root $x^*$ lies in the bracket, so $|m_n - x^*| \leq (b_n - a_n)/2 = (b-a)/2^{n+1}$. $\square$

**Convergence rate:** The error halves at each step — **linear convergence** with rate $1/2$. To gain one decimal digit of accuracy requires about $\log_2(10) \approx 3.3$ iterations. For 10-digit accuracy from a bracket of width 1: $n \approx 3.3 \times 10 = 33$ iterations. Bisection is slow but **guaranteed** to converge whenever a bracket exists.

**When to use bisection:** When the function is continuous but not differentiable; when a bracket is available but convergence of faster methods is uncertain; as a fallback when Newton–Raphson fails to converge.

---

## 22.3 Newton–Raphson for Scalar Equations

Newton–Raphson is the workhorse of root-finding. It exploits derivative information to converge far faster than bisection.

**Algorithm 22.2 (Newton–Raphson, Scalar).**

*Input:* Differentiable $f: \mathbb{R} \to \mathbb{R}$; initial guess $x_0$; tolerance $\varepsilon$.

1. For $n = 0, 1, 2, \ldots$:
   - Compute $x_{n+1} = x_n - f(x_n)/f'(x_n)$.
   - If $|x_{n+1} - x_n| < \varepsilon$: **return** $x_{n+1}$.

**Derivation from Taylor expansion.** Expand $f$ around $x_n$:

$$0 = f(x^*) \approx f(x_n) + f'(x_n)(x^* - x_n) \implies x^* \approx x_n - \frac{f(x_n)}{f'(x_n)} \equiv x_{n+1}.$$

The update $x_{n+1} = x_n - f(x_n)/f'(x_n)$ is the $x$-intercept of the tangent line to $f$ at $(x_n, f(x_n))$. Geometrically, Newton–Raphson replaces the nonlinear $f$ with its linear approximation and solves that.

**Theorem 22.2 (Quadratic Convergence of Newton–Raphson).** Suppose $f$ is twice continuously differentiable, $f(x^*) = 0$, $f'(x^*) \neq 0$, and the initial guess $x_0$ is sufficiently close to $x^*$. Then:

$$|x_{n+1} - x^*| \leq C|x_n - x^*|^2, \quad C = \frac{|f''(x^*)|}{2|f'(x^*)|}.$$

The error **squares** at each step — **quadratic convergence**. If $|x_0 - x^*| < 10^{-1}$, then after one step $|x_1 - x^*| \lesssim C\cdot 10^{-2}$; after two steps $|x_2 - x^*| \lesssim C^3 \cdot 10^{-4}$; after three steps $|x_3 - x^*| \lesssim C^7 \cdot 10^{-8}$. Three iterations from 1-digit accuracy gives 8-digit accuracy.

*Proof.* Let $e_n = x_n - x^*$. Taylor-expand $f(x_n)$ around $x^*$:

$$f(x_n) = f(x^*) + f'(x^*)e_n + \frac{1}{2}f''(\xi_n)e_n^2 = f'(x^*)e_n + \frac{1}{2}f''(\xi_n)e_n^2$$

for some $\xi_n$ between $x_n$ and $x^*$. The Newton step:

$$e_{n+1} = x_{n+1} - x^* = e_n - \frac{f(x_n)}{f'(x_n)} = e_n - \frac{f'(x^*)e_n + \frac{1}{2}f''(\xi_n)e_n^2}{f'(x_n)}.$$

Since $f'(x_n) \to f'(x^*)$ and $f'(x^*) \neq 0$:

$$e_{n+1} = e_n - \frac{f'(x^*)e_n}{f'(x^*)} - \frac{f''(\xi_n)e_n^2}{2f'(x_n)} = -\frac{f''(\xi_n)}{2f'(x_n)}e_n^2.$$

Taking absolute values: $|e_{n+1}| \leq C|e_n|^2$ with $C = |f''(\xi^*)|/(2|f'(x^*)|)$. $\square$

In APL, Newton–Raphson via the `⍣≡` fixed-point operator:

```apl
⍝ APL — Newton–Raphson via ⍣≡
⎕IO←0 ⋄ ⎕ML←1

⍝ Find steady state of Solow: sf(k) - μk = 0
⍝ f(k) = s*k^α - μ*k; we want the zero
alpha←1÷3  ⋄  s←0.25  ⋄  mu←0.08

f   ← {s×⍵*alpha) - mu×⍵}          ⍝ ODE right-hand side
df  ← {s×alpha×⍵*alpha-1) - mu}    ⍝ analytical derivative

⍝ One Newton step: x → x - f(x)/f'(x)
nr_step ← {⍵ - (f ⍵) ÷ df ⍵}

⍝ Iterate until convergence: ⍣≡ applies nr_step until two consecutive results agree
kstar ← nr_step ⍣ ≡ ⊢ 1            ⍝ start from k=1

kstar                                ⍝ ≈ 3.953
f kstar                              ⍝ should be ≈ 0 (residual check)
```

---

## 22.4 Newton–Raphson for Systems

For a system $\mathbf{F}(\mathbf{x}) = \mathbf{0}$ with $\mathbf{F}: \mathbb{R}^n \to \mathbb{R}^n$, Newton–Raphson generalizes by replacing the scalar derivative with the **Jacobian matrix**.

**Definition 22.1 (Jacobian).** The **Jacobian** of $\mathbf{F}$ at $\mathbf{x}$ is the $n\times n$ matrix:

$$J_F(\mathbf{x}) = \frac{\partial\mathbf{F}}{\partial\mathbf{x}} = \begin{pmatrix}\partial F_1/\partial x_1 & \cdots & \partial F_1/\partial x_n \\ \vdots & & \vdots \\ \partial F_n/\partial x_1 & \cdots & \partial F_n/\partial x_n\end{pmatrix}.$$

**Algorithm 22.3 (Newton–Raphson, System).**

*Input:* $\mathbf{F}: \mathbb{R}^n \to \mathbb{R}^n$; initial guess $\mathbf{x}_0$; tolerance $\varepsilon$.

1. For $n = 0, 1, 2, \ldots$:
   - Compute $J_F(\mathbf{x}_n)$.
   - Solve the linear system $J_F(\mathbf{x}_n)\Delta\mathbf{x} = -\mathbf{F}(\mathbf{x}_n)$ for $\Delta\mathbf{x}$.
   - Set $\mathbf{x}_{n+1} = \mathbf{x}_n + \Delta\mathbf{x}$.
   - If $\|\mathbf{x}_{n+1} - \mathbf{x}_n\| < \varepsilon$ or $\|\mathbf{F}(\mathbf{x}_{n+1})\| < \varepsilon$: **return** $\mathbf{x}_{n+1}$.

The key step is solving $J_F\Delta\mathbf{x} = -\mathbf{F}$ — a linear system in $\Delta\mathbf{x}$. This is the step that uses LU decomposition (Chapter 25). In APL, it is `delta_x ← ((-F_x) ⌹ J_x)`.

**Theorem 22.3 (Quadratic Convergence — System Newton).** If $\mathbf{F}$ is twice continuously differentiable in a neighborhood of $\mathbf{x}^*$, $J_F(\mathbf{x}^*)$ is invertible, and $\mathbf{x}_0$ is sufficiently close to $\mathbf{x}^*$, then:

$$\|\mathbf{x}_{n+1} - \mathbf{x}^*\| \leq C\|\mathbf{x}_n - \mathbf{x}^*\|^2,$$

with $C$ depending on the second derivatives of $\mathbf{F}$ and $\|J_F(\mathbf{x}^*)^{-1}\|$.

**Sufficient condition for convergence (Kantorovich theorem):** Let $h = \|J_F(\mathbf{x}_0)^{-1}\mathbf{F}(\mathbf{x}_0)\|\cdot\|J_F(\mathbf{x}_0)^{-1}\|\cdot L/2$ where $L$ is the Lipschitz constant of $J_F$. If $h \leq 1/2$, Newton–Raphson converges from $\mathbf{x}_0$.

---

## 22.5 Quasi-Newton Methods: Broyden's Method

Newton–Raphson requires computing the Jacobian $J_F$ at each iteration — expensive when $\mathbf{F}$ is costly to evaluate or when $J_F$ has no closed form. **Broyden's method** maintains an approximate Jacobian $B_n$ and updates it cheaply using only function evaluations.

**Algorithm 22.4 (Broyden's Method).**

*Input:* $\mathbf{F}$; initial guess $\mathbf{x}_0$; initial Jacobian approximation $B_0$ (often $B_0 = I$ or a finite-difference Jacobian).

For $n = 0, 1, 2, \ldots$:
1. Solve $B_n\Delta\mathbf{x}_n = -\mathbf{F}(\mathbf{x}_n)$; update $\mathbf{x}_{n+1} = \mathbf{x}_n + \Delta\mathbf{x}_n$.
2. Compute $\Delta\mathbf{F}_n = \mathbf{F}(\mathbf{x}_{n+1}) - \mathbf{F}(\mathbf{x}_n)$.
3. **Broyden update:** $B_{n+1} = B_n + \frac{(\Delta\mathbf{F}_n - B_n\Delta\mathbf{x}_n)\Delta\mathbf{x}_n'}{\Delta\mathbf{x}_n'\Delta\mathbf{x}_n}$.

The Broyden update is the **rank-1 correction** that enforces the secant condition $B_{n+1}\Delta\mathbf{x}_n = \Delta\mathbf{F}_n$, while changing $B_n$ as little as possible (in Frobenius norm). This generalizes the scalar secant method to systems.

**Convergence:** Broyden's method achieves **superlinear convergence** — faster than linear but slower than quadratic. It requires only one function evaluation per iteration (vs. $n+1$ for finite-difference Jacobian Newton), making it much faster per iteration for large systems.

---

## 22.6 Worked Example: Steady State of a Medium-Scale NK Model

*Cross-reference: Principles Ch. 23 (NK steady state)* **[P:Ch.23]**

A medium-scale NK model has 7 steady-state conditions in 7 unknowns $(C^*, K^*, Y^*, w^*, r^*, \pi^*, mc^*)$:

$$\mathbf{F}(\mathbf{x}) = \begin{pmatrix} C^* - (1-s_K)Y^* \\ K^* - s_K Y^* / \delta \\ Y^* - K^{*\alpha}n^{*1-\alpha} \\ w^* - (1-\alpha)Y^*/n^* \\ r^* - \alpha Y^*/K^* - (1-\delta) \\ \pi^* - \bar\pi \\ mc^* - w^*/[(1-\alpha)Y^*/n^*] \end{pmatrix} = \mathbf{0}$$

plus optimality conditions from household and price-setting behavior. Newton–Raphson solves this in 3–5 iterations from a reasonable starting point.

```python
import numpy as np
from scipy.optimize import fsolve

# Medium-scale NK steady state
alpha, delta, beta, theta_p = 0.36, 0.025, 0.99, 0.75
pi_bar, phi_pi, phi_y = 1.005, 1.5, 0.5  # gross inflation target

def nk_steady_state(x):
    C, K, Y, w, r, pi, mc = x
    n = 1/3  # fixed labor (for simplicity)
    
    F = np.zeros(7)
    F[0] = Y - K**alpha * n**(1-alpha)              # production
    F[1] = 1 - beta*(1+r)                            # Euler equation (r* = 1/β - 1)
    F[2] = r - (alpha*Y/K - delta)                   # MPK = r + δ
    F[3] = w - (1-alpha)*Y/n                         # MPL = w
    F[4] = mc - w/((1-alpha)*Y/n)                    # MC = w/MPL
    F[5] = pi - pi_bar                               # inflation target
    F[6] = C - (Y - delta*K)                         # goods market clearing
    return F

# Initial guess from approximate values
r_star = 1/beta - 1
K_star_approx = (alpha/(r_star+delta))**(1/(1-alpha)) * (1/3)**(1-alpha/(1-alpha))
x0 = [0.7, K_star_approx, 1.0, 2.0, r_star, pi_bar, 1.0]

x_star = fsolve(nk_steady_state, x0, full_output=False)
C_star, K_star, Y_star, w_star, r_star_sol, pi_star, mc_star = x_star
print(f"Steady state:")
print(f"  Y* = {Y_star:.4f}, K* = {K_star:.4f}, C* = {C_star:.4f}")
print(f"  r* = {r_star_sol*100:.2f}% (quarterly), mc* = {mc_star:.4f}")
print(f"  Residuals: {np.max(np.abs(nk_steady_state(x_star))):.2e}")

# Newton-Raphson from scratch (educational)
def newton_system(F, x0, tol=1e-12, max_iter=50):
    x = np.array(x0, dtype=float)
    for i in range(max_iter):
        Fx = F(x)
        if np.max(np.abs(Fx)) < tol: return x, i
        # Finite-difference Jacobian
        h = 1e-7
        J = np.zeros((len(x), len(x)))
        for j in range(len(x)):
            xp = x.copy(); xp[j] += h
            J[:,j] = (np.array(F(xp)) - Fx) / h
        dx = np.linalg.solve(J, -Fx)
        x = x + dx
    return x, max_iter

x_nr, iters = newton_system(nk_steady_state, x0)
print(f"\nNewton-Raphson converged in {iters} iterations")
print(f"Max residual: {np.max(np.abs(nk_steady_state(x_nr))):.2e}")
```

```julia
using NLsolve, LinearAlgebra

alpha, delta, beta = 0.36, 0.025, 0.99
pi_bar = 1.005; n = 1/3

function nk_ss!(F, x)
    C, K, Y, w, r, pi_, mc = x
    F[1] = Y - K^alpha * n^(1-alpha)
    F[2] = 1 - beta*(1+r)
    F[3] = r - (alpha*Y/K - delta)
    F[4] = w - (1-alpha)*Y/n
    F[5] = mc - w/((1-alpha)*Y/n)
    F[6] = pi_ - pi_bar
    F[7] = C - (Y - delta*K)
end

r_app = 1/beta - 1
K_app = (alpha/(r_app+delta))^(1/(1-alpha))*(n)^((1-alpha)/(1-alpha))
x0 = [0.7, K_app, 1.0, 2.0, r_app, pi_bar, 1.0]

sol = nlsolve(nk_ss!, x0, method=:newton, ftol=1e-12)
println("Converged: $(sol.f_converged), iterations: $(sol.iterations)")
println("Y*=$(round(sol.zero[3],digits=4)), K*=$(round(sol.zero[2],digits=4))")
```

```r
library(nleqslv)

alpha<-0.36; delta<-0.025; beta<-0.99; n<-1/3; pi_bar<-1.005

nk_ss <- function(x) {
  C<-x[1]; K<-x[2]; Y<-x[3]; w<-x[4]; r<-x[5]; pi_<-x[6]; mc<-x[7]
  c(Y - K^alpha*n^(1-alpha), 1 - beta*(1+r), r-(alpha*Y/K-delta),
    w-(1-alpha)*Y/n, mc-w/((1-alpha)*Y/n), pi_-pi_bar, C-(Y-delta*K))
}

r_app<-1/beta-1; K_app<-(alpha/(r_app+delta))^(1/(1-alpha))
x0<-c(0.7, K_app, 1.0, 2.0, r_app, pi_bar, 1.0)
sol<-nleqslv(x0, nk_ss, method="Newton")
cat(sprintf("Y*=%.4f, K*=%.4f, residual=%.2e\n",sol$x[3],sol$x[2],max(abs(sol$fvec))))
```

---

## 22.7 Programming Exercises

### Exercise 22.1 (APL — Broyden's Method)

Implement Broyden's method in APL as a dfn `broyden ← {F B0 x0 ← ⍵ ⋄ ...}`. The update step uses the rank-1 formula: `B ← B + ((dF - B+.×dx) ∘.× dx) ÷ dx+.×dx`. Test on the 2-equation Solow/RCK system and compare iteration counts to Newton–Raphson. Show that Broyden requires fewer function evaluations per solution.

### Exercise 22.2 (Python — Convergence Rate Measurement)

For Newton–Raphson applied to $f(x) = x^3 - 2$ starting from $x_0 = 2$: (a) run 10 iterations recording $e_n = |x_n - 2^{1/3}|$; (b) compute the convergence ratio $e_{n+1}/e_n^2$; (c) verify this ratio converges to $C = |f''(x^*)|/(2|f'(x^*)|) = (2\cdot 2^{-2/3})/(2\cdot 3\cdot 2^{-2/3}) = 1/3$; (d) plot $\log|e_n|$ vs. $n$ and verify the slope approximately doubles each step on the log scale.

### Exercise 22.3 (Julia — Ill-Conditioning)

```julia
# Effect of ill-conditioning on Newton convergence
using LinearAlgebra

function ill_conditioned_system(x, kappa=1000.0)
    # Near-singular system: condition number ≈ kappa
    A = [1.0 1.0; 1.0 1.0+1/kappa]
    b = [2.0, 2.0+1/kappa]
    A*x - b  # F(x) = Ax - b, solution x* = [1, 1]
end

for kappa in [1, 10, 100, 1000, 1e6]
    F!(f,x) = (f .= ill_conditioned_system(x, kappa))
    sol = nlsolve(F!, [0.5, 0.5], method=:newton)
    cond_J = cond([1.0 1.0; 1.0 1.0+1/kappa])
    println("κ=$(kappa): cond(J)=$(round(cond_J,digits=1)), " *
            "iters=$(sol.iterations), max_err=$(round(maximum(abs.(sol.zero .- 1)),digits=8))")
end
```

### Exercise 22.4 — RCK Steady State ($\star$)

The RCK model steady state satisfies: (1) $f'(k^*) = \delta + \rho + \sigma g$; (2) $c^* = f(k^*) - (n+g+\delta)k^*$. With $f(k) = k^\alpha$ these have analytical solutions. Now use $f(k) = [\gamma k^\rho + (1-\gamma)]^{1/\rho}$ (CES). (a) Write the two conditions as a $2\times2$ system $\mathbf{F}(k^*, c^*) = \mathbf{0}$. (b) Apply Newton–Raphson with finite-difference Jacobian. (c) Compare the CES steady state to the Cobb–Douglas approximation for $\rho \in \{-0.5, -0.1, 0, 0.1, 0.5\}$.

---

## 22.8 Chapter Summary

**Key results:**

- **Bisection**: guaranteed convergence at linear rate $1/2$; needs $n \approx 3.3\log_{10}(1/\varepsilon)$ iterations; requires a bracket $[a,b]$ with $f(a)f(b)<0$.
- **Newton–Raphson (scalar)**: $x_{n+1} = x_n - f(x_n)/f'(x_n)$; **quadratic convergence** proved in Theorem 22.2; error satisfies $|e_{n+1}| \leq C|e_n|^2$ near $x^*$.
- **Newton–Raphson (system)**: solve $J_F(\mathbf{x}_n)\Delta\mathbf{x} = -\mathbf{F}(\mathbf{x}_n)$ at each step via LU; quadratic convergence (Theorem 22.3); Kantorovich theorem gives sufficient conditions.
- **Broyden's method**: rank-1 Jacobian update $B_{n+1} = B_n + (\Delta\mathbf{F} - B_n\Delta\mathbf{x})\Delta\mathbf{x}'/(\Delta\mathbf{x}'\Delta\mathbf{x})$; superlinear convergence; one function evaluation per iteration.
- In APL: Newton–Raphson is `{⍵ - (f ⍵) ÷ df ⍵}⍣≡`; the system version is `{x - ((-F_x) ⌹ J_x)}⍣≡`.

*Next: Chapter 23 — Numerical Integration and Differentiation*
