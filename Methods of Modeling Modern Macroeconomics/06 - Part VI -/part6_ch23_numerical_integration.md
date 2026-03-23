# Chapter 23: Numerical Integration and Differentiation

*Computing Expected Utility and Marginal Effects*

> *"Integration and differentiation are inverses of each other analytically. Numerically, they pose completely different challenges."*

**Cross-reference:** *Principles* Ch. 11 (expected utility, Euler equation); Ch. 20 (SDF and asset pricing); Ch. 37 (DICE model, SCC integration) **[P:Ch.11, P:Ch.20, P:Ch.37]**

---

## 23.1 Why Numerical Integration?

Two fundamental operations appear constantly in quantitative macroeconomics:

**Expected utility:** The Bellman equation $V(a) = \max_c\{u(c) + \beta\mathbb{E}[V(a')]\}$ requires evaluating $\mathbb{E}[V(a')] = \int V(a')p(a'|a)da'$ where the expectation is over the continuous distribution of next-period assets. With Gaussian income shocks: $a' = (1+r)(a-c) + y'$, $y' \sim \mathcal{N}(\mu, \sigma^2)$, this is an integral over the normal density that has no closed form for general $V$.

**The SDF moment:** Asset pricing requires $\mathbb{E}[M_{t+1}R_{t+1}] = 1$ where $M_{t+1} = \beta(c_{t+1}/c_t)^{-\sigma}$ and $R_{t+1}$ has a continuous distribution. Computing this moment is a numerical integration problem.

**The SCC:** The social cost of carbon [P:Ch.37.2] is $SCC_t = -\mathbb{E}_t\int_t^\infty e^{-\int_t^s r_u\,du}\partial Y_s/\partial E_t\,ds$ — a stochastic integral over an infinite horizon, requiring numerical quadrature.

This chapter develops three classes of methods: Newton–Cotes rules (simple but limited accuracy), Gaussian quadrature (highly accurate for smooth integrands), and Monte Carlo (flexible but slow).

---

## 23.2 Newton–Cotes Rules: Trapezoidal and Simpson's

Newton–Cotes rules approximate $\int_a^b f(x)dx$ by replacing $f$ with a polynomial interpolant on a uniform grid.

**Algorithm 23.1 (Composite Trapezoidal Rule).**

Partition $[a,b]$ into $n$ equal subintervals of width $h = (b-a)/n$ with nodes $x_i = a + ih$:

$$\int_a^b f(x)\,dx \approx T_n(f) = h\left[\frac{f(x_0)}{2} + f(x_1) + f(x_2) + \cdots + f(x_{n-1}) + \frac{f(x_n)}{2}\right].$$

**Theorem 23.1 (Trapezoidal Error Bound).** If $f \in C^2[a,b]$:

$$\left|\int_a^b f\,dx - T_n(f)\right| \leq \frac{(b-a)h^2}{12}\max_{x\in[a,b]}|f''(x)| = O(h^2).$$

**Algorithm 23.2 (Composite Simpson's Rule).**

With $n$ even and $h = (b-a)/n$:

$$\int_a^b f\,dx \approx S_n(f) = \frac{h}{3}\left[f(x_0) + 4f(x_1) + 2f(x_2) + 4f(x_3) + \cdots + 4f(x_{n-1}) + f(x_n)\right].$$

Simpson's rule uses quadratic interpolation on each pair of subintervals.

**Theorem 23.2 (Simpson's Error Bound).** If $f \in C^4[a,b]$:

$$\left|\int_a^b f\,dx - S_n(f)\right| \leq \frac{(b-a)h^4}{180}\max_{x\in[a,b]}|f^{(4)}(x)| = O(h^4).$$

Simpson's converges twice as fast as trapezoidal in terms of the error exponent. For the same number of function evaluations, Simpson's is almost always preferred.

In APL, both rules exploit the `+/` reduce and the `×` weighting:

```apl
⍝ APL — Trapezoidal and Simpson's rules
⎕IO←0 ⋄ ⎕ML←1

⍝ Trapezoidal: weights are 1/2, 1, 1, ..., 1, 1/2 × h
trap ← {f a b n ← ⍵
    h ← (b-a)÷n
    x ← a + h × ⍳n+1
    w ← h × 0.5, (n-1)⍴1, 0.5    ⍝ endpoint weights 1/2, interior 1
    w +.× f¨x}

⍝ Simpson's: weights 1, 4, 2, 4, 2, ..., 4, 1 × h/3
simp ← {f a b n ← ⍵    ⍝ n must be even
    h ← (b-a)÷n
    x ← a + h × ⍳n+1
    inner_w ← (n-1) ⍴ 4 2    ⍝ alternating 4, 2
    w ← (h÷3) × 1, inner_w, 1
    w +.× f¨x}

⍝ Test: ∫₀¹ x² dx = 1/3
f_sq ← {⍵*2}
trap f_sq 0 1 100    ⍝ should be ≈ 0.3333
simp f_sq 0 1 100    ⍝ should be exactly 1/3 (quadratic integrand, Simpson is exact)
```

---

## 23.3 Gaussian Quadrature

Newton–Cotes rules use equally spaced nodes — which is suboptimal. **Gaussian quadrature** chooses nodes and weights to maximize the degree of polynomial exactness for a fixed number of function evaluations.

**Definition 23.1 (Gaussian Quadrature).** An $n$-point Gaussian quadrature rule approximates:

$$\int_a^b w(x)f(x)\,dx \approx \sum_{i=1}^n \omega_i f(x_i),$$

where $w(x)$ is a **weight function**, the nodes $\{x_i\}$ are the roots of the $n$-th orthogonal polynomial with respect to $w$, and the weights $\{\omega_i\}$ are chosen so the rule is exact for all polynomials of degree $\leq 2n-1$.

**Theorem 23.3 (Exactness of Gaussian Quadrature).** The $n$-point Gauss rule is exact for all polynomials of degree $\leq 2n-1$. No $n$-point rule can be exact for all polynomials of degree $\leq 2n$.

*Proof sketch.* The $2n$ free parameters (nodes + weights) are chosen to satisfy $2n$ moment conditions: exactness for $1, x, x^2, \ldots, x^{2n-1}$. For any rule with fixed nodes, only $n$ weight conditions are available (exactness up to degree $n-1$). Gauss quadrature uses both nodes and weights as free parameters, achieving degree $2n-1$ exactly. $\square$

### 23.3.1 Gauss–Hermite Quadrature for Normal Integrals

For integrals over the real line against the standard normal density — the most common case in macroeconomics with Gaussian shocks:

$$\int_{-\infty}^\infty f(x)\phi(x)\,dx \approx \sum_{i=1}^n \omega_i f(x_i),$$

where $\phi(x) = e^{-x^2/2}/\sqrt{2\pi}$ and the nodes $x_i$ are roots of the $n$-th Hermite polynomial $H_n(x)$.

**Key property:** 5–10 quadrature points achieve the accuracy of 1000+ Monte Carlo draws for smooth integrands — a 100–200x speedup. This is why Gauss–Hermite quadrature is used in the inner loop of dynamic programming (Chapter 15) and in continuous-state Bellman equations.

For the general normal $y \sim \mathcal{N}(\mu, \sigma^2)$: transform $x = (y-\mu)/\sigma$, then $\mathbb{E}[f(y)] = \int f(\mu+\sigma x)\phi(x)dx \approx \sum_i\omega_i f(\mu + \sigma x_i)$.

```apl
⍝ APL — Gauss-Hermite quadrature for E[f(Y)] where Y~N(μ,σ²)
⎕IO←0 ⋄ ⎕ML←1

⍝ Pre-computed GH nodes and weights for n=5 points
⍝ (standard normal basis: nodes are roots of H5(x)/√(2π))
gh_nodes   ← ¯2.02018 ¯0.95858 0 0.95858 2.02018
gh_weights ← 0.01995 0.39362 0.94531 0.39362 0.01995
⍝ Weights already normalized for standard normal

⍝ E[f(Y)] where Y~N(mu, sig²): change of variables y = mu + sig*x
E_gauss_hermite ← {f mu sig ← ⍵
    y_nodes ← mu + sig × gh_nodes    ⍝ transform to general normal
    gh_weights +.× f¨ y_nodes}        ⍝ weighted sum

⍝ Test: E[Y²] = mu² + sig² for Y~N(mu, sig)
f_sq ← {⍵*2}
E_gauss_hermite f_sq 1.5 0.5    ⍝ should be 1.5²+0.5² = 2.5
```

```python
import numpy as np
from numpy.polynomial.hermite import hermgauss

def gauss_hermite_expectation(f, mu, sigma, n=7):
    """E[f(Y)] where Y ~ N(mu, sigma^2) via Gauss-Hermite quadrature."""
    # hermgauss returns nodes/weights for ∫f(x)exp(-x²)dx
    # Transform: Y = mu + sqrt(2)*sigma*x
    nodes, weights = hermgauss(n)
    y_nodes = mu + np.sqrt(2) * sigma * nodes
    return np.sum(weights * f(y_nodes)) / np.sqrt(np.pi)

# Test: E[Y^2] = mu^2 + sigma^2
f = lambda y: y**2
mu, sigma = 1.5, 0.5
result = gauss_hermite_expectation(f, mu, sigma, n=5)
exact = mu**2 + sigma**2
print(f"GH(n=5): {result:.8f}  Exact: {exact:.8f}  Error: {abs(result-exact):.2e}")

# Compare to MC for CRRA expected utility
sigma_crra = 2.0
f_crra = lambda c: c**(1-sigma_crra) / (1-sigma_crra)
mu_c, sig_c = 1.0, 0.1

print("\nExpected CRRA utility (varying n):")
for n_pts in [3, 5, 7, 10, 15]:
    result_gh = gauss_hermite_expectation(f_crra, mu_c, sig_c, n=n_pts)
    print(f"  n={n_pts:2d}: {result_gh:.10f}")

# MC benchmark
np.random.seed(42)
c_draws = np.random.normal(mu_c, sig_c, 100000)
mc_result = np.mean(f_crra(c_draws[c_draws > 0]))
print(f"  MC(100k): {mc_result:.10f}")
```

```julia
using FastGaussQuadrature

function gh_expectation(f, mu, sigma, n=7)
    nodes, weights = gausshermite(n)
    # Transform: integrate f(mu + sqrt(2)*sigma*x) * exp(-x^2) dx / sqrt(π)
    y = mu .+ sqrt(2)*sigma.*nodes
    return sum(weights .* f.(y)) / sqrt(π)
end

# Test
println("E[Y^2] via GH(7): ", round(gh_expectation(y->y^2, 1.5, 0.5, 7), digits=8))
println("Exact:            ", 1.5^2 + 0.5^2)
```

```r
library(statmod)

gh_expectation <- function(f, mu, sigma, n=7) {
  gh <- gauss.quad.prob(n, dist="normal", mu=mu, sigma=sigma)
  sum(gh$weights * f(gh$nodes))
}

cat(sprintf("E[Y^2] via GH(7): %.8f\n", gh_expectation(function(y) y^2, 1.5, 0.5, 7)))
cat(sprintf("Exact:            %.8f\n", 1.5^2 + 0.5^2))
```

---

## 23.4 Monte Carlo Integration

**Definition 23.2 (Monte Carlo Estimator).** For $I = \int f(x)p(x)dx$ where $p$ is a density, the **Monte Carlo estimator** based on $N$ draws $\{x_i\}_{i=1}^N \sim p$ is:

$$\hat{I}_N = \frac{1}{N}\sum_{i=1}^N f(x_i).$$

**Theorem 23.4 (Monte Carlo Central Limit Theorem).** If $\mathbb{E}[f(x)^2] < \infty$:

$$\sqrt{N}(\hat{I}_N - I) \xrightarrow{d} \mathcal{N}(0, \sigma_f^2), \quad \sigma_f^2 = \text{Var}[f(x)].$$

The Monte Carlo standard error is $\sigma_f/\sqrt{N}$ — decreasing at rate $N^{-1/2}$ regardless of the dimension of the integration domain. **This dimension-independence is Monte Carlo's main advantage:** for high-dimensional integrals (many state variables), Newton–Cotes and quadrature methods suffer from the curse of dimensionality, while Monte Carlo converges at the same rate.

### 23.4.1 Importance Sampling

When the integrand $f(x)p(x)$ is concentrated in a small region relative to $p$, direct Monte Carlo is inefficient. **Importance sampling** draws from an alternative distribution $q(x)$ that better covers the support of the integrand:

$$I = \int f(x)p(x)dx = \int f(x)\frac{p(x)}{q(x)}q(x)dx \approx \frac{1}{N}\sum_{i=1}^N f(x_i)\frac{p(x_i)}{q(x_i)},$$

where $x_i \sim q$. The ratio $p(x)/q(x)$ is the **importance weight**. The optimal $q^*(x) \propto |f(x)|p(x)$ minimizes variance, but requires knowing $|f|p$ up to a constant.

### 23.4.2 Quasi-Monte Carlo

**Quasi-Monte Carlo (QMC)** replaces pseudo-random draws with **low-discrepancy sequences** (Halton, Sobol) that cover the domain more uniformly. For smooth integrands in $d$ dimensions, QMC achieves error $O((\log N)^d/N)$ — much faster than MC's $O(N^{-1/2})$ for moderate $d$.

---

## 23.5 Numerical Differentiation

**Finite differences** approximate derivatives when analytical formulas are unavailable or costly to derive.

**Definition 23.3 (Finite Difference Approximations).**

Forward difference: $f'(x) \approx \frac{f(x+h)-f(x)}{h}$ — error $O(h)$.

Backward difference: $f'(x) \approx \frac{f(x)-f(x-h)}{h}$ — error $O(h)$.

**Central difference:** $f'(x) \approx \frac{f(x+h)-f(x-h)}{2h}$ — error $O(h^2)$.

**Complex step:** $f'(x) \approx \text{Im}[f(x+ih)]/h$ for complex-valued extension of $f$ — error $O(h^2)$ but without cancellation errors (unlike real finite differences where $h$ too small causes catastrophic cancellation).

**Step-size selection:** The optimal $h$ balances truncation error ($O(h^p)$ for the method) against rounding error ($O(\varepsilon_M/h)$ where $\varepsilon_M \approx 10^{-16}$ is machine epsilon). For the central difference, the optimal step is $h^* = (\varepsilon_M)^{1/3} \approx 10^{-5\cdot\frac{1}{3}} \approx 10^{-5}$ (for double precision).

**Automatic differentiation (AD):** Modern software (PyTorch, JAX, Zygote.jl, ForwardDiff.jl) computes exact derivatives (to machine precision) of arbitrary programs via the chain rule, without finite differences. AD is increasingly preferred in DSGE estimation because it enables gradient-based optimizers for the likelihood.

---

## 23.6 Worked Example: Expected CRRA Utility Under Income Uncertainty

*Cross-reference: Principles Ch. 11.2 (stochastic Euler equation)* **[P:Ch.11.2]**

A household has CRRA utility $u(c) = c^{1-\sigma}/(1-\sigma)$, $\sigma = 2$. Income next period is $y' = \bar{y}e^{\varepsilon'}$ where $\varepsilon' \sim \mathcal{N}(0, 0.04)$ (log-normal income shocks with $\sigma_\varepsilon = 0.2$). Savings $s$ generate next-period wealth $w' = (1+r)s$, and consumption $c' = w' + y'$. We wish to compute:

$$\mathbb{E}[u(c')] = \mathbb{E}\!\left[\frac{(w' + \bar{y}e^{\varepsilon'})^{1-\sigma}}{1-\sigma}\right].$$

**Method 1 — Gauss–Hermite (7 points):** Transform $\varepsilon' = \sqrt{2}\sigma_\varepsilon z$; compute at 7 nodes.

**Method 2 — Monte Carlo (100,000 draws):** Direct simulation.

**Method 3 — Trapezoidal rule:** Integrate over a finite range $[-3\sigma_\varepsilon, 3\sigma_\varepsilon]$.

```python
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.stats import norm

sigma_crra, r, w_prime_val = 2.0, 0.04, 1.0
y_bar, sigma_eps = 1.0, 0.20

def u(c): return c**(1-sigma_crra)/(1-sigma_crra) if c > 0 else -1e20
def u_vec(c): return np.where(c > 1e-10, c**(1-sigma_crra)/(1-sigma_crra), -1e20)

def E_u_gauss_hermite(w_prime, n=7):
    nodes, weights = hermgauss(n)
    # eps' = sqrt(2)*sigma_eps*z (change of variables for hermgauss)
    eps = np.sqrt(2)*sigma_eps*nodes
    c_prime = w_prime + y_bar*np.exp(eps)
    return np.sum(weights * u_vec(c_prime)) / np.sqrt(np.pi)

def E_u_montecarlo(w_prime, N=100000):
    eps = np.random.normal(0, sigma_eps, N)
    c_prime = w_prime + y_bar*np.exp(eps)
    return np.mean(u_vec(c_prime))

def E_u_trapezoidal(w_prime, n=1000):
    eps_range = 3*sigma_eps
    eps_grid = np.linspace(-eps_range, eps_range, n)
    c_prime = w_prime + y_bar*np.exp(eps_grid)
    integrand = u_vec(c_prime) * norm.pdf(eps_grid, 0, sigma_eps)
    h = 2*eps_range/(n-1)
    return h*(0.5*integrand[0] + np.sum(integrand[1:-1]) + 0.5*integrand[-1])

np.random.seed(42)
print("Expected utility comparison (w'=1.0, σ=2, σ_eps=0.2):")
print(f"  Gauss-Hermite (n=7):     {E_u_gauss_hermite(w_prime_val, 7):.8f}")
print(f"  Gauss-Hermite (n=15):    {E_u_gauss_hermite(w_prime_val, 15):.8f}")
print(f"  Monte Carlo (N=100k):    {E_u_montecarlo(w_prime_val):.8f}")
print(f"  Trapezoidal (n=1000):    {E_u_trapezoidal(w_prime_val):.8f}")

# Efficiency comparison: time for 1% accuracy
import time
for method, func in [("GH-7", lambda: E_u_gauss_hermite(1.0, 7)),
                     ("GH-15", lambda: E_u_gauss_hermite(1.0, 15)),
                     ("MC-1k", lambda: E_u_montecarlo(1.0, 1000)),
                     ("Trap-100", lambda: E_u_trapezoidal(1.0, 100))]:
    t0 = time.perf_counter()
    [func() for _ in range(100)]
    elapsed = (time.perf_counter()-t0)/100
    print(f"  {method}: {elapsed*1e6:.1f} μs per evaluation")
```

---

## 23.7 Programming Exercises

### Exercise 23.1 (APL — Gauss–Hermite Integration)

Store the $n=7$ Gauss–Hermite nodes and weights as APL vectors. (a) Write a dfn `gh_E ← {f mu sig ← ⍵ ⋄ gh_weights +.× f¨ mu+sig×sqrt(2)×gh_nodes) ÷ sqrt(⍨○1)}` that computes $\mathbb{E}[f(Y)]$ for $Y \sim \mathcal{N}(\mu, \sigma^2)$. (b) Verify: `gh_E ({⍵*2}) 0 1` should give 1.0 (since $\mathbb{E}[Z^2]=1$). (c) Use it to evaluate the stochastic continuation value in a simplified buffer-stock Bellman equation.

### Exercise 23.2 (Python — Quadrature Convergence)

For the integral $I = \mathbb{E}[e^{-\sigma Y}]$ where $Y \sim \mathcal{N}(0, 1)$ (which equals $e^{\sigma^2/2}$ analytically), compare the absolute error of: (a) Gauss–Hermite with $n = 3, 5, 7, 10, 15$ points; (b) Trapezoidal with $n = 10, 50, 100, 1000$ points; (c) Monte Carlo with $N = 10^3, 10^4, 10^5, 10^6$ draws. Plot all three error curves on a log-log scale. Confirm GH achieves near-machine-precision accuracy with $n = 15$ while MC needs $10^6$ draws for 3-digit accuracy.

### Exercise 23.3 (Julia — Automatic Differentiation)

```julia
using ForwardDiff

# Jacobian of NK steady-state conditions via AD
alpha, delta, beta = 0.36, 0.025, 0.99; n = 1/3; pi_bar = 1.005

function nk_ss(x)
    C,K,Y,w,r,pi_,mc = x
    [Y - K^alpha*n^(1-alpha),
     1 - beta*(1+r),
     r - (alpha*Y/K - delta),
     w - (1-alpha)*Y/n,
     mc - w/((1-alpha)*Y/n),
     pi_ - pi_bar,
     C - (Y - delta*K)]
end

r_star = 1/beta - 1
K_star = (alpha/(r_star+delta))^(1/(1-alpha))
x_star = [0.65, K_star, 1.0, 1.5, r_star, pi_bar, 0.95]

# AD Jacobian (exact, no finite differences!)
J_ad = ForwardDiff.jacobian(nk_ss, x_star)
println("Condition number of Jacobian: ", round(cond(J_ad), digits=2))
println("Max eigenvalue:               ", round(maximum(abs.(eigvals(J_ad))), digits=4))

# Compare to finite-difference Jacobian
h = 1e-7
J_fd = zeros(7, 7)
F0 = nk_ss(x_star)
for j in 1:7
    xp = copy(x_star); xp[j] += h
    J_fd[:,j] = (nk_ss(xp) - F0) / h
end
println("Max error (AD vs FD): ", maximum(abs.(J_ad - J_fd)))
```

### Exercise 23.4 — Importance Sampling ($\star$)

For the tail integral $I = P(Y > k) = \int_k^\infty\phi(y)dy$ where $Y \sim \mathcal{N}(0,1)$ and $k = 3$ (a rare event): (a) direct MC estimate with $N = 10^5$ draws and note the variance; (b) importance sampling with proposal $q(y) = \text{Exp}(1)$ translated to $[k,\infty)$; (c) show the IS variance is dramatically lower; (d) compute the effective sample size $N_{eff} = (\sum w_i)^2/\sum w_i^2$ for both methods.

---

## 23.8 Chapter Summary

**Key results:**

- **Trapezoidal rule**: $O(h^2)$ error; simple but slow. **Simpson's rule**: $O(h^4)$ error; $4\times$ faster convergence for smooth $f$.
- **Gauss–Hermite quadrature**: $n$-point rule is exact for polynomials of degree $\leq 2n-1$; achieves near-machine-precision for smooth integrands with $n \approx 7$–15 points — replacing 1000+ Monte Carlo draws.
- **Monte Carlo**: convergence rate $O(N^{-1/2})$ independent of dimension; preferred for high-dimensional integrals or non-smooth integrands.
- **Importance sampling** reduces variance by drawing from a proposal close to the integrand; **QMC** with low-discrepancy sequences achieves $O((\log N)^d/N)$ for smooth $d$-dimensional integrands.
- **Central finite difference**: $O(h^2)$ error; optimal step $h \approx \varepsilon_M^{1/3} \approx 10^{-5}$. **Automatic differentiation** (ForwardDiff, JAX) computes exact derivatives of code without finite differences.
- In APL: Gauss–Hermite is `gh_weights +.× f¨ y_nodes` — weighted inner product; trapezoidal is `h × 0.5×(f a)+f b + +/ f¨ interior_nodes`.

*Next: Chapter 24 — Numerical Optimization: Finding Optimal Policy Rules*
