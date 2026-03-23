# Part I: Foundations of Mathematical Macroeconomics

*Connects to: Principles Parts I–II*

---

This part builds the mathematical scaffolding that every subsequent part uses. Readers who have taken a year of calculus and a semester of linear algebra will find it a targeted review with macroeconomic applications — not a repetition of what they already know, but a reframing of familiar tools in the language of economic models. Each chapter introduces its tools in the context of a specific macroeconomic question, so that the mathematics never floats free of motivation.

Chapter 1 develops optimization theory — the mathematical heart of all macroeconomic modeling, since macroeconomics is the study of agents who optimize subject to constraints. Chapter 2 applies linear algebra to economic systems, deriving the Leontief inverse and the IS–LM solution as matrix problems. Chapters 3 and 4 develop differential and difference equations respectively, covering the analytical toolkit for the dynamic models of Parts III–IV. Chapter 5 introduces the stochastic processes — particularly the AR(1) — that appear as the shock processes in every model from Chapter 17 onward.

A note on difficulty: Chapter 1 through Chapter 5 are sequenced so that each builds on the previous, but they are also largely self-contained. A reader who wants only the stochastic methods of Part V can read Chapter 5 directly after Chapter 2.

---

# Chapter 1: Calculus and Optimization in Macroeconomics

*From Utility Maximization to Lagrange Multipliers*

> *"Optimization is the heart of economics. Almost every economic model, at its core, is an optimization problem."*
> — Thomas Sargent

**Cross-reference:** *Principles* Ch. 1 (scope of macroeconomics); Ch. 11 (consumption theory, Euler equation); Ch. 12 (investment theory, Tobin's q); Ch. 13 (labour supply and demand) **[P:Ch.1, P:Ch.11, P:Ch.12, P:Ch.13]**

---

## 1.1 Why Calculus? The Optimization Imperative

Macroeconomics, at every level of sophistication, is a discipline about constrained choice. A household chooses a consumption path to maximize lifetime utility subject to a budget constraint. A firm chooses investment to maximize the present value of profits subject to an adjustment cost technology. A central bank chooses a path for the interest rate to minimize a loss function subject to the constraints imposed by the Phillips curve and the IS relationship. In each case, the economic content is carried by the optimization problem, and the calculus is the language in which that problem is solved.

The models of *Principles* frequently invoke optimization without deriving the conditions in full. Chapter 11 of *Principles* states the Euler equation $u'(c_t) = \beta(1+r_{t+1})\mathbb{E}_t[u'(c_{t+1})]$ and then explains it verbally; this chapter derives it, together with every other first-order condition that appears in *Principles*, from the underlying optimization problem. The derivations are not decorative. They matter because:

1. They reveal what assumptions the result depends on — and what happens when those assumptions are relaxed.
2. They provide the algebraic form needed for numerical computation (Part VI).
3. They are the basis for comparative statics — understanding how the equilibrium changes when parameters change.

We begin with unconstrained optimization, move to equality-constrained problems (the Lagrangian), then to inequality-constrained problems (Kuhn–Tucker conditions), and finally to the implicit function theorem, which turns comparative statics from a graphical exercise into an algebraic one.

---

## 1.2 Unconstrained Optimization

### 1.2.1 Single-Variable Problems

> **Definition 1.1 (Local and Global Maxima).** Let $f: \mathbb{R} \to \mathbb{R}$. A point $x^* \in \mathbb{R}$ is a **local maximum** of $f$ if there exists $\varepsilon > 0$ such that $f(x^*) \geq f(x)$ for all $x$ with $|x - x^*| < \varepsilon$. It is a **global maximum** if $f(x^*) \geq f(x)$ for all $x \in \mathbb{R}$.

The fundamental result connecting calculus to optimization:

> **Theorem 1.1 (First-Order Necessary Condition).** If $f$ is differentiable and $x^*$ is a local maximum (or minimum) of $f$, then $f'(x^*) = 0$.

*Proof sketch.* If $f'(x^*) > 0$, then by continuity of $f'$ there exists $\delta > 0$ such that $f'(x) > 0$ for all $x \in (x^* - \delta, x^* + \delta)$. Then $f$ is strictly increasing on this interval, so $f(x^* + \delta/2) > f(x^*)$, contradicting $x^*$ being a local maximum. Similarly for $f'(x^*) < 0$. $\square$

The first-order condition (FOC) $f'(x^*) = 0$ is necessary but not sufficient. Points where $f'(x) = 0$ are called **critical points** or **stationary points**; they may be maxima, minima, or inflection points.

> **Theorem 1.2 (Second-Order Sufficient Conditions).** Let $f$ be twice differentiable and let $x^*$ satisfy $f'(x^*) = 0$.
> - If $f''(x^*) < 0$, then $x^*$ is a strict local maximum.
> - If $f''(x^*) > 0$, then $x^*$ is a strict local minimum.
> - If $f''(x^*) = 0$, the test is inconclusive.

*Application to macroeconomics:* The Baumol–Tobin money demand model of *Principles* Ch. 14 minimizes total cost $TC(n) = iY/(2n) + bn$ over the number of bank trips $n$. Setting $TC'(n) = 0$:

$$-\frac{iY}{2n^2} + b = 0 \implies n^* = \sqrt{\frac{iY}{2b}}.$$

The second derivative $TC''(n) = iY/n^3 > 0$ confirms this is a minimum.

### 1.2.2 Multivariate Unconstrained Optimization

Let $f: \mathbb{R}^n \to \mathbb{R}$. The gradient is:

$$\nabla f(\mathbf{x}) = \left(\frac{\partial f}{\partial x_1}, \frac{\partial f}{\partial x_2}, \ldots, \frac{\partial f}{\partial x_n}\right)'.$$

> **Theorem 1.3 (Multivariate FOC).** If $\mathbf{x}^*$ is a local maximum or minimum of a differentiable $f: \mathbb{R}^n \to \mathbb{R}$, then $\nabla f(\mathbf{x}^*) = \mathbf{0}$.

The second-order conditions require the **Hessian matrix**:

$$H_f(\mathbf{x}) = \begin{pmatrix} \frac{\partial^2 f}{\partial x_1^2} & \frac{\partial^2 f}{\partial x_1 \partial x_2} & \cdots \\ \frac{\partial^2 f}{\partial x_2 \partial x_1} & \frac{\partial^2 f}{\partial x_2^2} & \cdots \\ \vdots & \vdots & \ddots \end{pmatrix}.$$

By Young's theorem (symmetry of mixed partials), $H_f$ is symmetric.

> **Definition 1.2 (Positive and Negative Definite Matrices).** A symmetric matrix $A$ is **positive definite** if $\mathbf{v}'A\mathbf{v} > 0$ for all nonzero $\mathbf{v} \in \mathbb{R}^n$; equivalently, all eigenvalues of $A$ are strictly positive. It is **negative definite** if $\mathbf{v}'A\mathbf{v} < 0$ for all nonzero $\mathbf{v}$; equivalently, all eigenvalues are strictly negative. It is **indefinite** if it has both positive and negative eigenvalues.

> **Theorem 1.4 (Multivariate Second-Order Conditions).** At a critical point $\mathbf{x}^*$:
> - $H_f(\mathbf{x}^*)$ negative definite $\Rightarrow$ strict local maximum.
> - $H_f(\mathbf{x}^*)$ positive definite $\Rightarrow$ strict local minimum.
> - $H_f(\mathbf{x}^*)$ indefinite $\Rightarrow$ saddle point.

A practical check for negative definiteness: the leading principal minors of $-H_f$ must all be positive (i.e., $H_f$ is negative definite iff $-H_f$ is positive definite iff all leading principal minors of $-H_f$ are positive).

### 1.2.3 The Envelope Theorem

The envelope theorem answers the question: how does the optimized value of the objective function change when a parameter changes, without re-solving the optimization problem?

> **Theorem 1.5 (Envelope Theorem, Unconstrained).** Let $f(x; \alpha)$ be differentiable in both $x$ and the parameter $\alpha$. Define the value function $V(\alpha) = \max_x f(x; \alpha)$, with maximizer $x^*(\alpha)$ (assuming it is differentiable in $\alpha$). Then:
>
> $$\frac{dV}{d\alpha} = \frac{\partial f(x^*(\alpha); \alpha)}{\partial \alpha}.$$

That is, the total derivative of the value function with respect to a parameter equals the partial derivative of the objective with respect to that parameter, evaluated at the optimum. The channel through which $\alpha$ changes the optimal $x^*$ — the indirect effect — contributes nothing to the change in $V$ because $f_x(x^*;\alpha) = 0$ at the optimum.

*Proof:*

$$\frac{dV}{d\alpha} = \frac{d}{d\alpha} f(x^*(\alpha);\alpha) = \underbrace{f_x(x^*(\alpha);\alpha)}_{=0 \text{ by FOC}} \cdot \frac{dx^*}{d\alpha} + f_\alpha(x^*(\alpha);\alpha) = f_\alpha(x^*(\alpha);\alpha). \quad \square$$

*Application:* In the firm's investment problem of *Principles* Ch. 12, the Hamiltonian is the value function and the shadow price $q_t$ (Tobin's $q$) is the envelope of the profit function with respect to the capital stock. The costate equation $\dot{q} = rq - F_K$ is the envelope theorem applied to the continuous-time optimal control problem — we derive this fully in Chapter 11.

---

## 1.3 Constrained Optimization: Equality Constraints

Most optimization problems in macroeconomics involve constraints. The household maximizes utility subject to a budget constraint; the firm maximizes profit subject to a production technology. Equality constraints are handled by the method of Lagrange multipliers.

### 1.3.1 The Lagrangian Method

> **Definition 1.3 (Lagrangian).** For the problem $\max_{\mathbf{x}} f(\mathbf{x})$ subject to $g(\mathbf{x}) = 0$, the **Lagrangian function** is:
>
> $$\mathcal{L}(\mathbf{x}, \lambda) = f(\mathbf{x}) - \lambda\, g(\mathbf{x}),$$
>
> where $\lambda \in \mathbb{R}$ is the **Lagrange multiplier** or **shadow price** associated with the constraint.

> **Theorem 1.6 (Lagrange Necessary Conditions).** If $\mathbf{x}^*$ is a local maximum of $f$ subject to $g(\mathbf{x}) = 0$, and if the constraint qualification $\nabla g(\mathbf{x}^*) \neq \mathbf{0}$ holds, then there exists $\lambda^*$ such that:
>
> $$\nabla_{\mathbf{x}} \mathcal{L}(\mathbf{x}^*, \lambda^*) = \mathbf{0} \quad \text{and} \quad g(\mathbf{x}^*) = 0.$$
>
> In components: $\partial f/\partial x_i = \lambda^* \partial g/\partial x_i$ for all $i$, plus the constraint.

*Derivation from first principles.* Consider $\max_{x_1, x_2} f(x_1, x_2)$ subject to $g(x_1, x_2) = 0$. Since the constraint implicitly defines $x_2$ as a function of $x_1$ near the optimum (by the IFT, provided $g_{x_2} \neq 0$), we can write $x_2 = h(x_1)$ where $h'(x_1) = -g_{x_1}/g_{x_2}$. Substituting into the objective and taking the total derivative with respect to $x_1$:

$$\frac{d}{dx_1}f(x_1, h(x_1)) = f_{x_1} + f_{x_2} h'(x_1) = f_{x_1} - f_{x_2}\frac{g_{x_1}}{g_{x_2}} = 0.$$

Rearranging: $f_{x_1}/g_{x_1} = f_{x_2}/g_{x_2} \equiv \lambda^*$. This gives $f_{x_i} = \lambda^* g_{x_i}$ for $i=1,2$, which are exactly the FOCs of the Lagrangian. $\square$

### 1.3.2 Economic Interpretation of the Lagrange Multiplier

The Lagrange multiplier $\lambda^*$ measures the marginal value of relaxing the constraint. Formally, if the constraint is $g(\mathbf{x}) = c$ (parameterized by $c$) and $V(c) = \max_{\mathbf{x}: g(\mathbf{x})=c} f(\mathbf{x})$, then $dV/dc = \lambda^*$. Each unit increase in $c$ allows the agent to reach a higher objective by exactly $\lambda^*$.

In the household problem with budget constraint $p_1 x_1 + p_2 x_2 = m$, the Lagrange multiplier is the marginal utility of income: $\lambda^* = dV/dm$ where $V(m)$ is the indirect utility function. In the firm's problem with a capital constraint, $\lambda^*$ is the shadow price of capital — Tobin's $q$ [P:Ch.12].

### 1.3.3 Multiple Equality Constraints

For $\max f(\mathbf{x})$ subject to $g_j(\mathbf{x}) = 0$, $j = 1, \ldots, m$ (with $m < n$), the Lagrangian is:

$$\mathcal{L}(\mathbf{x}, \bm{\lambda}) = f(\mathbf{x}) - \sum_{j=1}^m \lambda_j g_j(\mathbf{x}),$$

and the FOCs are $\nabla_{\mathbf{x}} \mathcal{L} = \mathbf{0}$ and $g_j(\mathbf{x}) = 0$ for all $j$, giving $n + m$ equations in $n + m$ unknowns.

**Worked Example 1.1: Household Utility Maximization**

> A household maximizes $u(c_1, c_2) = c_1^{1-\sigma}/(1-\sigma) + \beta c_2^{1-\sigma}/(1-\sigma)$ subject to $c_1 + c_2/(1+r) = w$, where $\beta \in (0,1)$, $\sigma > 0$, and $w$ is lifetime wealth. Find the optimal consumption plan and the Lagrange multiplier.

*Solution.* The Lagrangian:

$$\mathcal{L} = \frac{c_1^{1-\sigma}}{1-\sigma} + \beta\frac{c_2^{1-\sigma}}{1-\sigma} - \lambda\left(c_1 + \frac{c_2}{1+r} - w\right).$$

First-order conditions:

$$\frac{\partial\mathcal{L}}{\partial c_1} = c_1^{-\sigma} - \lambda = 0 \implies \lambda = c_1^{-\sigma},$$

$$\frac{\partial\mathcal{L}}{\partial c_2} = \beta c_2^{-\sigma} - \frac{\lambda}{1+r} = 0 \implies \lambda = \beta(1+r) c_2^{-\sigma}.$$

Equating: $c_1^{-\sigma} = \beta(1+r) c_2^{-\sigma}$, which gives the **Euler equation**:

$$\frac{c_2}{c_1} = [\beta(1+r)]^{1/\sigma}.$$

This is exactly equation (11.14) in *Principles* Ch. 11 [P:Ch.11.2], now derived explicitly. Substituting into the budget constraint and solving:

$$c_1^* = \frac{w}{1 + \frac{[\beta(1+r)]^{1/\sigma}}{1+r}}, \quad c_2^* = [\beta(1+r)]^{1/\sigma} c_1^*, \quad \lambda^* = (c_1^*)^{-\sigma}.$$

The shadow price $\lambda^*$ is the marginal utility of wealth: gaining one more unit of $w$ raises lifetime utility by $(c_1^*)^{-\sigma}$. $\square$

---

## 1.4 Constrained Optimization: Inequality Constraints

When constraints may or may not bind at the optimum — as in the borrowing constraint of *Principles* Ch. 11.4 — the Kuhn–Tucker conditions replace the simple Lagrange conditions.

### 1.4.1 The Kuhn–Tucker Conditions

> **Theorem 1.7 (Kuhn–Tucker Necessary Conditions).** Consider $\max_{\mathbf{x}} f(\mathbf{x})$ subject to $g_j(\mathbf{x}) \leq 0$ for $j = 1, \ldots, m$. Define the Lagrangian $\mathcal{L}(\mathbf{x}, \bm{\mu}) = f(\mathbf{x}) - \sum_j \mu_j g_j(\mathbf{x})$. If $\mathbf{x}^*$ is a local maximum and a constraint qualification holds, then there exist $\mu_j^* \geq 0$ such that:
>
> 1. **Stationarity:** $\nabla_{\mathbf{x}} \mathcal{L}(\mathbf{x}^*, \bm{\mu}^*) = \mathbf{0}$
> 2. **Primal feasibility:** $g_j(\mathbf{x}^*) \leq 0$ for all $j$
> 3. **Dual feasibility:** $\mu_j^* \geq 0$ for all $j$
> 4. **Complementary slackness:** $\mu_j^* g_j(\mathbf{x}^*) = 0$ for all $j$

The complementary slackness condition is the key: either the constraint is binding ($g_j(\mathbf{x}^*) = 0$) or the multiplier is zero ($\mu_j^* = 0$), or both. If a constraint is not binding, it plays no role in the optimum.

*Application to the liquidity-constrained consumer [P:Ch.11.4]:* The household maximizes $u(c_t) + \beta\mathbb{E}[u(c_{t+1})]$ subject to $a_{t+1} = (1+r)(a_t + y_t - c_t)$ and the borrowing constraint $a_{t+1} \geq \underline{b}$. Writing the constraint as $-a_{t+1} + \underline{b} \leq 0$, the KKT conditions add a multiplier $\mu_t \geq 0$:

$$u'(c_t) = \beta(1+r)\mathbb{E}[u'(c_{t+1})] + \mu_t, \quad \mu_t(a_{t+1} - \underline{b}) = 0.$$

When the constraint binds ($a_{t+1} = \underline{b}$): $\mu_t = u'(c_t) - \beta(1+r)\mathbb{E}[u'(c_{t+1})] > 0$ — the household is "constrained" in the sense that the Euler equation holds with a positive wedge. When the constraint is slack: $\mu_t = 0$ and the standard Euler equation holds. This is the formal derivation of *Principles* equation (11.16).

---

## 1.5 The Implicit Function Theorem and Comparative Statics

Comparative statics asks: when a parameter changes, how does the equilibrium change? The implicit function theorem (IFT) provides the algebraic tool for answering this question systematically.

> **Theorem 1.8 (Implicit Function Theorem — Scalar Version).** Let $F(x, \alpha): \mathbb{R}^2 \to \mathbb{R}$ be continuously differentiable. Suppose $F(x^*, \alpha^*) = 0$ and $F_x(x^*, \alpha^*) \neq 0$. Then there exists a neighborhood of $\alpha^*$ and a unique continuously differentiable function $x = x(\alpha)$ such that $x(\alpha^*) = x^*$ and $F(x(\alpha), \alpha) = 0$. Moreover:
>
> $$\frac{dx^*}{d\alpha} = -\frac{F_\alpha(x^*, \alpha^*)}{F_x(x^*, \alpha^*)}.$$

The multivariate version: if $\mathbf{F}(\mathbf{x}, \bm{\alpha}) = \mathbf{0}$ defines $\mathbf{x}$ implicitly as a function of $\bm{\alpha}$ near $(\mathbf{x}^*, \bm{\alpha}^*)$, and if the Jacobian $\partial\mathbf{F}/\partial\mathbf{x}$ is nonsingular at $(\mathbf{x}^*, \bm{\alpha}^*)$, then:

$$\frac{d\mathbf{x}^*}{d\bm{\alpha}} = -\left(\frac{\partial\mathbf{F}}{\partial\mathbf{x}}\right)^{-1} \frac{\partial\mathbf{F}}{\partial\bm{\alpha}}.$$

This is the fundamental formula for comparative statics in any system of equilibrium conditions.

**Worked Example 1.2: Comparative Statics for the IS Curve**

> Derive the slope of the IS curve (how $Y$ changes with $r$ in goods-market equilibrium) using the IFT.

The IS curve is the locus satisfying $F(Y, r) \equiv Y - C(Y - T) - I(r) - G = 0$. Then:

$$F_Y = 1 - C'(Y - T), \quad F_r = -I'(r).$$

By the IFT:

$$\left.\frac{dY}{dr}\right|_{\text{IS}} = -\frac{F_r}{F_Y} = -\frac{-I'(r)}{1 - C'} = \frac{I'(r)}{1 - C'} < 0,$$

since $I'(r) < 0$ and $1 - C' > 0$. The IS curve slopes downward, with steeper slope when the MPC $C'$ is large (strong multiplier) or investment is insensitive to the interest rate (small $|I'|$). This is the algebraic derivation of the IS slope stated qualitatively in *Principles* Ch. 9.1 [P:Ch.9.1]. $\square$

---

## 1.6 Concavity, Convexity, and Global Optima

The FOCs and SOCs identify local optima. For global optima, we need concavity or convexity.

> **Definition 1.4 (Concave and Convex Functions).** A function $f: C \to \mathbb{R}$ defined on a convex set $C \subseteq \mathbb{R}^n$ is **concave** if for all $\mathbf{x}, \mathbf{y} \in C$ and $\theta \in [0,1]$:
>
> $$f(\theta\mathbf{x} + (1-\theta)\mathbf{y}) \geq \theta f(\mathbf{x}) + (1-\theta) f(\mathbf{y}).$$
>
> It is **strictly concave** if the inequality is strict for $\mathbf{x} \neq \mathbf{y}$ and $\theta \in (0,1)$.

> **Proposition 1.1.** If $f$ is strictly concave, any local maximum is also a global maximum. If $f$ is twice differentiable, it is concave iff $H_f$ is negative semi-definite everywhere.

Most utility functions used in macroeconomics are strictly concave: $u(c) = \ln c$, $u(c) = c^{1-\sigma}/(1-\sigma)$ for $\sigma > 0$, and $u(c) = -(c - c^*)^2$ are all strictly concave in $c$. This guarantees that the Lagrange FOCs identify a global maximum.

For convex constraint sets (which are the rule in consumer and producer theory), the constraint qualification is automatically satisfied, and the Lagrange/KKT conditions are both necessary and sufficient under concavity of the objective.

---

## 1.7 Worked Example: Deriving the Labour Supply Curve

*Cross-reference: Principles Ch. 13 (labour supply and demand)* **[P:Ch.13.1]**

Consider a household with preferences $u(c, \ell) = \ln c + \chi \ln \ell$, where $c$ is consumption and $\ell = 1 - n$ is leisure ($n$ is hours worked), and $\chi > 0$ is the relative weight on leisure. The budget constraint is $c = w n + a$ (wage income plus non-labour income $a$), and time is normalized so $\ell + n = 1$.

**Step 1: Write the constrained problem.**

$$\max_{c, n} \ln c + \chi \ln(1-n) \quad \text{s.t.} \quad c = wn + a.$$

Substitute the constraint: $\max_n \ln(wn + a) + \chi\ln(1-n)$.

**Step 2: First-order condition.**

$$\frac{d}{dn}\left[\ln(wn+a) + \chi\ln(1-n)\right] = \frac{w}{wn+a} - \frac{\chi}{1-n} = 0.$$

**Step 3: Solve for $n^*$.**

$$w(1-n) = \chi(wn+a) \implies w - wn = \chi wn + \chi a \implies n^*(w,a) = \frac{w - \chi a}{w(1+\chi)}.$$

**Step 4: Comparative statics.**

$$\frac{\partial n^*}{\partial w} = \frac{1}{1+\chi} + \frac{\chi a}{w^2(1+\chi)} > 0 \quad \text{(substitution effect dominates for } a > 0\text{)}.$$

$$\frac{\partial n^*}{\partial a} = \frac{-\chi}{w(1+\chi)} < 0 \quad \text{(labour supply falls with non-labour income)}.$$

**Step 5: Verify second-order condition.**

$$\frac{d^2}{dn^2}\left[\ln(wn+a)+\chi\ln(1-n)\right] = -\frac{w^2}{(wn+a)^2} - \frac{\chi}{(1-n)^2} < 0.$$

The objective is strictly concave in $n$, confirming $n^*$ is a global maximum.

The upward-sloping labour supply curve ($\partial n^*/\partial w > 0$) is consistent with the competitive benchmark of *Principles* Ch. 13.1 [P:Ch.13.1]. The formula $n^* = (w - \chi a)/[w(1+\chi)]$ also shows that the Frisch elasticity of labour supply (the elasticity of $n$ with respect to $w$ holding the marginal utility of wealth constant) depends on the preference parameter $\chi$ — a result used in Chapter 13 of this volume when calibrating the RBC model.

---

## 1.8 Programming Exercises

### Exercise 1.1 (APL)

Implement Newton's method for finding the root of $f'(x) = 0$ where $f(x) = -x^4 + 3x^3 - 2$ using the APL `⍣` power operator with a convergence guard.

```apl
⍝ Dyalog APL — Newton's method via ⍣ (power operator)
⎕IO←0 ⋄ ⎕ML←1

⍝ Define f and its derivatives
f  ← {-⍵*4 + 3×⍵*3 - 2}
df ← {-4×⍵*3 + 9×⍵*2}      ⍝ f'
ddf← {-12×⍵*2 + 18×⍵}       ⍝ f''

⍝ One Newton step on df: find critical point of f
step ← {⍵ - (df ⍵) ÷ (ddf ⍵)}

⍝ Iterate to fixed point (convergence when |step(x)-x| < 1e-10)
converged ← {1e¯10 > |⍵ - step ⍵}
x_star ← step ⍣ converged ⊢ 0.5   ⍝ start from x=0.5

x_star          ⍝ display result
f x_star        ⍝ display f(x*)
```

Note the idiom `f ⍣ converged ⊢ start`: apply `f` repeatedly until the predicate `converged` is satisfied. `⍣` with a dyadic right argument iterates until convergence — this pattern appears throughout the book for fixed-point computations.

### Exercise 1.2 (Python)

```python
# Python — Newton's method for critical point of f(x) = -x**4 + 3*x**3 - 2
from scipy.optimize import minimize_scalar
import numpy as np

f   = lambda x: -x**4 + 3*x**3 - 2
df  = lambda x: -4*x**3 + 9*x**2    # first derivative
ddf = lambda x: -12*x**2 + 18*x     # second derivative

def newton_critical(x0, tol=1e-10, max_iter=100):
    x = x0
    for i in range(max_iter):
        x_new = x - df(x) / ddf(x)
        if abs(x_new - x) < tol:
            return x_new, i+1
        x = x_new
    raise RuntimeError("Did not converge")

x_star, iters = newton_critical(0.5)
print(f"x* = {x_star:.10f}, f(x*) = {f(x_star):.10f}, converged in {iters} iterations")
```

### Exercise 1.3 (Julia)

```julia
# Julia — Newton's method for critical point
f(x)   = -x^4 + 3x^3 - 2
df(x)  = -4x^3 + 9x^2
ddf(x) = -12x^2 + 18x

function newton_critical(x0; tol=1e-10, max_iter=100)
    x = x0
    for i in 1:max_iter
        x_new = x - df(x) / ddf(x)
        abs(x_new - x) < tol && return x_new, i
        x = x_new
    end
    error("Did not converge")
end

x_star, iters = newton_critical(0.5)
println("x* = $(round(x_star, digits=10)),  f(x*) = $(round(f(x_star), digits=10)),  iters = $iters")
```

### Exercise 1.4 (R)

```r
# R — Newton's method for critical point
f   <- function(x) -x^4 + 3*x^3 - 2
df  <- function(x) -4*x^3 + 9*x^2
ddf <- function(x) -12*x^2 + 18*x

newton_critical <- function(x0, tol = 1e-10, max_iter = 100) {
  x <- x0
  for (i in seq_len(max_iter)) {
    x_new <- x - df(x) / ddf(x)
    if (abs(x_new - x) < tol) return(list(x = x_new, iters = i))
    x <- x_new
  }
  stop("Did not converge")
}

res <- newton_critical(0.5)
cat(sprintf("x* = %.10f,  f(x*) = %.10f,  iters = %d\n", res$x, f(res$x), res$iters))
```

### Exercise 1.5 — Household Problem ($\star$)

Extend Worked Example 1.1 to three periods. The household maximizes $\sum_{t=1}^3 \beta^{t-1} u(c_t)$ with $u(c) = \ln c$ subject to the lifetime budget constraint $\sum_{t=1}^3 c_t/(1+r)^{t-1} = w$. (a) Derive the three Euler equations. (b) Solve analytically for $c_t^*$ as a function of $w$, $r$, and $\beta$. (c) In APL, solve the system using `⌹` (matrix divide) for a range of parameter values and produce a table of consumption shares $c_t^*/w$.

### Exercise 1.6 — Comparative Statics ($\star$)

Using the IS condition $Y = C(Y-T) + I(r) + G$ with $C(y) = a + by$, compute $\partial Y^*/\partial G$, $\partial Y^*/\partial T$, and $\partial Y^*/\partial r$ using the IFT. Verify these match the IS–LM multipliers derived in Chapter 6. Implement the computation in APL as a dfn that takes the parameters $(a, b, I', G, T)$ as a vector and returns the three partial derivatives.

### Exercise 1.7 — Kuhn–Tucker ($\star\star$)

Consider a household that maximizes $u(c_1, c_2) = c_1^{1/2} + \beta c_2^{1/2}$ subject to $c_1 + c_2/(1+r) = w$ and $c_1 \geq \underline{c}$ (a consumption floor). (a) Write the Lagrangian with two multipliers. (b) Identify the two regimes (constraint binding vs. not binding) and solve for the optimum in each. (c) Find the threshold wealth $\bar{w}$ below which the floor binds.

---

## 1.9 Chapter Summary

This chapter developed the optimization tools that underlie every macroeconomic model in *Principles* and in this volume.

**Key results:**

- The **first-order necessary condition** $\nabla f(\mathbf{x}^*) = \mathbf{0}$ must hold at any interior maximum or minimum of a differentiable function.
- The **second-order sufficient condition** requires the Hessian to be negative definite (for a maximum) or positive definite (for a minimum) at the critical point.
- The **envelope theorem** states that the marginal value of a parameter equals its direct effect on the objective, evaluated at the optimum; indirect effects through the optimal policy function contribute nothing.
- The **Lagrange method** converts equality-constrained problems into unconstrained problems over an augmented objective; the multiplier is the shadow price of the constraint.
- The **Kuhn–Tucker conditions** extend the Lagrange method to inequality constraints, with complementary slackness ensuring that a slack constraint has a zero multiplier.
- The **implicit function theorem** provides the formula for comparative statics: $d\mathbf{x}^*/d\bm{\alpha} = -(F_{\mathbf{x}})^{-1} F_{\bm{\alpha}}$.

**Connections forward:** Chapter 3 extends these tools to infinite-horizon continuous-time problems (the Hamiltonians of optimal control). Chapter 11 applies them to the RCK model. Chapter 15 develops the discrete-time analogue (Bellman equations and dynamic programming). Chapter 18 applies the Lagrange method to the household's problem under rational expectations.

---

*Next: Chapter 2 — Linear Algebra for Macroeconomic Systems*
