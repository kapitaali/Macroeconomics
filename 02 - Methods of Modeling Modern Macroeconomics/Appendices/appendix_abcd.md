# Appendix A: Mathematical Formulas and Tables

---

## A.1 Differentiation Rules

**Basic rules.** For differentiable functions $f, g$ and constant $c$:

| Rule | Formula |
|---|---|
| Constant | $(c)' = 0$ |
| Power | $(x^n)' = nx^{n-1}$ |
| Sum | $(f+g)' = f' + g'$ |
| Product | $(fg)' = f'g + fg'$ |
| Quotient | $(f/g)' = (f'g - fg')/g^2$ |
| Chain | $(f(g(x)))' = f'(g(x))\cdot g'(x)$ |
| Inverse | If $y = f^{-1}(x)$: $dy/dx = 1/f'(y)$ |

**Common functions:**

| Function | Derivative |
|---|---|
| $e^x$ | $e^x$ |
| $a^x$ | $a^x\ln a$ |
| $\ln x$ | $1/x$ |
| $\log_a x$ | $1/(x\ln a)$ |
| $\sin x$ | $\cos x$ |
| $\cos x$ | $-\sin x$ |
| $\arctan x$ | $1/(1+x^2)$ |

**Implicit differentiation.** If $F(x, y) = 0$ defines $y$ as a function of $x$: $dy/dx = -F_x/F_y$.

**Logarithmic differentiation.** For $y = f(x)^{g(x)}$: $\ln y = g(x)\ln f(x)$, then differentiate both sides.

## A.2 Integration Formulas

| Integrand | Antiderivative |
|---|---|
| $x^n$ ($n\neq-1$) | $x^{n+1}/(n+1)$ |
| $1/x$ | $\ln|x|$ |
| $e^{ax}$ | $e^{ax}/a$ |
| $\ln x$ | $x\ln x - x$ |
| $\sin ax$ | $-\cos(ax)/a$ |
| $\cos ax$ | $\sin(ax)/a$ |
| $1/(1+x^2)$ | $\arctan x$ |

**Integration by parts:** $\int u\,dv = uv - \int v\,du$.

**Gaussian integral:** $\int_{-\infty}^\infty e^{-ax^2}dx = \sqrt{\pi/a}$ for $a > 0$.

**Present value integral:** $\int_0^\infty e^{-\rho t}f(t)dt = \mathcal{L}\{f\}(\rho)$ (Laplace transform at $\rho$).

## A.3 Taylor Series Expansions

$$f(x) = \sum_{n=0}^\infty\frac{f^{(n)}(a)}{n!}(x-a)^n \quad \text{(Taylor series around }a\text{)}$$

**Standard expansions around 0:**

| Function | Expansion |
|---|---|
| $e^x$ | $1 + x + x^2/2! + x^3/3! + \cdots$ |
| $\ln(1+x)$ | $x - x^2/2 + x^3/3 - \cdots$ ($|x|<1$) |
| $(1+x)^\alpha$ | $1 + \alpha x + \frac{\alpha(\alpha-1)}{2}x^2 + \cdots$ |
| $1/(1-x)$ | $1 + x + x^2 + x^3 + \cdots$ ($|x|<1$) |
| $\sin x$ | $x - x^3/6 + x^5/120 - \cdots$ |
| $\cos x$ | $1 - x^2/2 + x^4/24 - \cdots$ |

**Key approximations (small $x$):** $\ln(1+x) \approx x$, $(1+x)^\alpha \approx 1 + \alpha x$, $e^x \approx 1 + x + x^2/2$.

## A.4 Matrix Identities

**Woodbury identity:** $(A + UCV)^{-1} = A^{-1} - A^{-1}U(C^{-1}+VA^{-1}U)^{-1}VA^{-1}$.

**Sherman–Morrison:** $(A + \mathbf{u}\mathbf{v}')^{-1} = A^{-1} - \frac{A^{-1}\mathbf{u}\mathbf{v}'A^{-1}}{1+\mathbf{v}'A^{-1}\mathbf{u}}$.

**Matrix determinant lemma:** $\det(A + \mathbf{u}\mathbf{v}') = (1+\mathbf{v}'A^{-1}\mathbf{u})\det(A)$.

**Trace and determinant:** $\text{tr}(AB) = \text{tr}(BA)$. $\det(AB) = \det(A)\det(B)$. $\det(A^{-1}) = 1/\det(A)$.

**Kronecker product:** $(A\otimes B)(C\otimes D) = (AC)\otimes(BD)$. $\text{vec}(AXB) = (B'\otimes A)\text{vec}(X)$.

**Differentiation:** $\partial(A\mathbf{x})/\partial\mathbf{x} = A$. $\partial(\mathbf{x}'A\mathbf{x})/\partial\mathbf{x} = (A+A')\mathbf{x}$. $\partial\ln\det(A)/\partial A = (A^{-1})'$.

## A.5 Special Functions in Macroeconomics

**Log-normal:** If $\ln X \sim \mathcal{N}(\mu,\sigma^2)$, then $X \sim \text{LogNormal}(\mu,\sigma^2)$ with $\mathbb{E}[X] = e^{\mu+\sigma^2/2}$, $\text{Var}[X] = (e^{\sigma^2}-1)e^{2\mu+\sigma^2}$.

**Gamma function:** $\Gamma(n) = (n-1)!$ for integer $n$; $\Gamma(1/2) = \sqrt{\pi}$; $\Gamma(n+1) = n\Gamma(n)$.

**Normal CDF:** $\Phi(z) = P(Z\leq z)$ for $Z\sim\mathcal{N}(0,1)$; $\Phi(-z) = 1-\Phi(z)$.

---

# Appendix B: Linear Algebra Review

---

## B.1 Vector Spaces

**Definition.** A **vector space** over $\mathbb{R}$ is a set $V$ with addition and scalar multiplication satisfying 8 axioms (closure, associativity, commutativity, identity, inverses, distributivity).

**Basis and dimension.** A set $\{\mathbf{v}_1,\ldots,\mathbf{v}_n\}$ is a **basis** if it is linearly independent and spans $V$. The **dimension** $\dim V$ is the cardinality of any basis.

**Standard basis** of $\mathbb{R}^n$: $\mathbf{e}_i$ has 1 in position $i$ and 0 elsewhere.

## B.2 Matrix Operations

**Rank.** $\text{rank}(A)$ = number of linearly independent rows (= columns). $\text{rank}(A) \leq \min(m,n)$ for $A\in\mathbb{R}^{m\times n}$. Full column rank: $\text{rank}(A) = n$ (columns independent). Full row rank: $\text{rank}(A) = m$.

**Null space.** $\mathcal{N}(A) = \{\mathbf{x}: A\mathbf{x}=\mathbf{0}\}$. Dimension: $n - \text{rank}(A)$ (rank-nullity theorem).

**Norms.** $\|\mathbf{x}\|_1 = \sum|x_i|$; $\|\mathbf{x}\|_2 = \sqrt{\sum x_i^2}$; $\|\mathbf{x}\|_\infty = \max|x_i|$. Matrix norm: $\|A\|_2 = \sigma_{\max}(A)$ (largest singular value); $\|A\|_F = \sqrt{\text{tr}(A'A)}$ (Frobenius).

## B.3 Determinants and Cramer's Rule

For $A\in\mathbb{R}^{2\times2}$: $\det(A) = a_{11}a_{22} - a_{12}a_{21}$.

**Cramer's rule.** For $A\mathbf{x} = \mathbf{b}$ with $\det(A)\neq0$: $x_i = \det(A_i)/\det(A)$, where $A_i$ is $A$ with column $i$ replaced by $\mathbf{b}$.

**Properties.** $\det(AB) = \det(A)\det(B)$. $\det(A') = \det(A)$. $\det(\alpha A) = \alpha^n\det(A)$ (for $n\times n$). $\det(A) = \prod_i\lambda_i$ (product of eigenvalues).

## B.4 Eigenvalues and Eigenvectors

$A\mathbf{v} = \lambda\mathbf{v}$, $\mathbf{v}\neq\mathbf{0}$. Eigenvalues are roots of $\det(A-\lambda I) = 0$ (characteristic polynomial).

**Spectral decomposition** (symmetric $A=A'$): $A = Q\Lambda Q'$ where $Q$ is orthonormal and $\Lambda = \text{diag}(\lambda_1,\ldots,\lambda_n)$.

**Stability:** Discrete: $A^t \to 0$ iff $|\lambda_i| < 1$ for all $i$. Continuous: $e^{At} \to 0$ iff $\text{Re}(\lambda_i) < 0$ for all $i$.

**Generalized eigenvalues.** $A\mathbf{v} = \lambda B\mathbf{v}$ has generalized eigenvalues $\lambda_i = T_{ii}/S_{ii}$ from the QZ decomposition $QAZ = S$, $QBZ = T$.

## B.5 The Perron–Frobenius Theorem

**Theorem (Perron–Frobenius).** Let $A$ be a square matrix with all positive entries. Then:
1. $A$ has a unique largest real eigenvalue $\lambda_1 > 0$ (the **Perron root**).
2. The corresponding eigenvector $\mathbf{v}_1$ has all positive entries (**Perron vector**).
3. $|\lambda_i| < \lambda_1$ for all other eigenvalues.

**Application in input-output analysis (Chapter 2).** For the Leontief technical coefficient matrix $A$ (with $A_{ij} \geq 0$), the Perron root $\lambda_1(A) < 1$ guarantees $(I-A)^{-1}$ exists and is non-negative — the Leontief inverse has all non-negative entries (backward linkage multipliers are non-negative).

**Application in Markov chains.** For a stochastic matrix $P$ (row sums = 1, all entries $\geq 0$), Perron–Frobenius guarantees a unique stationary distribution $\bm\pi$ with $\bm\pi'P = \bm\pi'$ and $\pi_i > 0$ for all $i$ (ergodic chains).

## B.6 QR and LU Decompositions

**LU decomposition.** $PA = LU$: $P$ permutation, $L$ unit lower triangular, $U$ upper triangular. Solves $A\mathbf{x}=\mathbf{b}$ in $O(n^3/3)$ flops (Chapter 25).

**QR decomposition.** $A = QR$: $Q$ orthonormal ($Q'Q=I$), $R$ upper triangular. Used for numerically stable OLS (Chapter 25). $\kappa(R) = \kappa(A)$ vs. $\kappa(A'A) = \kappa(A)^2$ for normal equations.

**Cholesky.** For positive definite $A$: $A = LL'$ ($L$ lower triangular). Fastest symmetric system solver. Used for sampling from multivariate normal (Chapter 26).

**SVD.** $A = U\Sigma V'$: $U, V$ orthonormal, $\Sigma = \text{diag}(\sigma_1,\ldots,\sigma_r,0,\ldots)$. Rank = number of positive singular values. Used for numerical rank determination (Chapter 41).

---

# Appendix C: Calculus Review

---

## C.1 Limits and Continuity

**Limit.** $\lim_{x\to a}f(x) = L$: for every $\varepsilon>0$ there exists $\delta>0$ such that $|x-a|<\delta \Rightarrow |f(x)-L|<\varepsilon$.

**L'Hôpital's rule.** If $\lim f = \lim g = 0$ (or $\pm\infty$): $\lim f/g = \lim f'/g'$ (when the right-hand limit exists).

**Useful limits:** $\lim_{x\to0}(1+x)^{1/x} = e$. $\lim_{n\to\infty}(1+r/n)^n = e^r$. $\lim_{x\to0}\frac{\sin x}{x} = 1$.

## C.2 Multivariate Differentiation

**Partial derivative.** $\partial f/\partial x_i$ = derivative of $f$ treating all $x_j$ ($j\neq i$) as constants.

**Gradient.** $\nabla f = (\partial f/\partial x_1, \ldots, \partial f/\partial x_n)' \in \mathbb{R}^n$.

**Hessian.** $H_f = [\partial^2f/\partial x_i\partial x_j]$ — symmetric $n\times n$ matrix of second partial derivatives. $H_f \succ 0$ (positive definite) iff $f$ is strictly convex.

**Jacobian.** For $F:\mathbb{R}^n\to\mathbb{R}^m$: $J_F = [\partial F_i/\partial x_j]$ — the $m\times n$ matrix of first partial derivatives.

**Chain rule (vector form).** For $h = f\circ g$ ($h:\mathbb{R}^m\to\mathbb{R}^p$, $g:\mathbb{R}^n\to\mathbb{R}^m$, $f:\mathbb{R}^m\to\mathbb{R}^p$): $J_h = J_f(g(\mathbf{x}))\cdot J_g(\mathbf{x})$.

## C.3 Optimization Conditions

**Unconstrained.** FOC: $\nabla f(\mathbf{x}^*) = \mathbf{0}$. SOC (minimum): $H_f(\mathbf{x}^*) \succ 0$.

**Equality constraints (Lagrange).** Maximize $f(\mathbf{x})$ s.t. $g(\mathbf{x}) = 0$: Lagrangian $\mathcal{L} = f - \lambda g$; FOCs: $\nabla f = \lambda\nabla g$.

**Inequality constraints (KKT).** Maximize $f$ s.t. $g_j(\mathbf{x}) \leq 0$: KKT conditions: $\nabla f = \sum_j\mu_j\nabla g_j$; $\mu_j \geq 0$; $\mu_jg_j = 0$ (complementary slackness).

**Envelope theorem.** For $V(\alpha) = \max_x f(x,\alpha)$ s.t. $g(x,\alpha)=0$: $dV/d\alpha = \partial\mathcal{L}/\partial\alpha|_{x=x^*(\alpha)}$.

## C.4 Integration

**Fundamental theorem of calculus.** $\frac{d}{dx}\int_a^x f(t)dt = f(x)$. $\int_a^b f'(x)dx = f(b) - f(a)$.

**Fubini's theorem.** For integrable $f$: $\int\!\int f(x,y)dxdy = \int\!\left[\int f(x,y)dx\right]dy$.

**Change of variables.** $\int_{\phi(a)}^{\phi(b)}f(x)dx = \int_a^b f(\phi(t))\phi'(t)dt$.

**Leibniz rule.** $\frac{d}{d\alpha}\int_{a(\alpha)}^{b(\alpha)}f(x,\alpha)dx = f(b,\alpha)b'(\alpha) - f(a,\alpha)a'(\alpha) + \int_a^b\frac{\partial f}{\partial\alpha}dx$.

**Dominated convergence theorem.** If $f_n \to f$ pointwise and $|f_n| \leq g$ (integrable), then $\int f_n \to \int f$.

---

# Appendix D: Probability Distributions and Statistical Tables

---

## D.1 Core Distributions

### Normal Distribution: $\mathcal{N}(\mu, \sigma^2)$

PDF: $f(x) = \frac{1}{\sigma\sqrt{2\pi}}\exp\!\left[-\frac{(x-\mu)^2}{2\sigma^2}\right]$.

Mean: $\mu$. Variance: $\sigma^2$. MGF: $M(t) = e^{\mu t + \sigma^2t^2/2}$.

**Standard normal** $\mathcal{N}(0,1)$: CDF $\Phi(z)$. Key quantiles: $\Phi^{-1}(0.025) = -1.960$, $\Phi^{-1}(0.05) = -1.645$, $\Phi^{-1}(0.10) = -1.282$.

### Log-Normal Distribution: $\text{LogN}(\mu, \sigma^2)$

$X \sim \text{LogN}(\mu,\sigma^2)$ iff $\ln X \sim \mathcal{N}(\mu,\sigma^2)$.

Mean: $e^{\mu+\sigma^2/2}$. Variance: $(e^{\sigma^2}-1)e^{2\mu+\sigma^2}$. Median: $e^\mu$.

**Role in macroeconomics:** Income and productivity shocks ($A_t = e^{z_t}$, $z_t \sim \mathcal{N}$); asset prices; firm sizes (Chapter 33).

### Gamma Distribution: $\text{Gamma}(\alpha, \beta)$

PDF: $f(x) = \frac{\beta^\alpha}{\Gamma(\alpha)}x^{\alpha-1}e^{-\beta x}$, $x > 0$.

Mean: $\alpha/\beta$. Variance: $\alpha/\beta^2$.

**Role:** Prior for positive parameters (CRRA $\sigma$, shock persistence); the chi-squared is $\chi^2(k) = \text{Gamma}(k/2, 1/2)$.

### Beta Distribution: $\text{Beta}(\alpha, \beta)$

PDF: $f(x) = \frac{x^{\alpha-1}(1-x)^{\beta-1}}{B(\alpha,\beta)}$, $x\in(0,1)$.

Mean: $\alpha/(\alpha+\beta)$. Variance: $\alpha\beta/[(\alpha+\beta)^2(\alpha+\beta+1)]$.

**Role:** Prior for parameters in $[0,1]$ — Calvo probability $\theta$, AR persistence $\rho$, capital share $\alpha$.

### Inverse-Gamma Distribution: $\text{IG}(\alpha, \beta)$

$X \sim \text{IG}(\alpha,\beta)$ iff $1/X \sim \text{Gamma}(\alpha,\beta)$.

Mean: $\beta/(\alpha-1)$ ($\alpha>1$). Variance: $\beta^2/[(\alpha-1)^2(\alpha-2)]$ ($\alpha>2$).

**Role:** Prior for variance parameters ($\sigma^2_\varepsilon$) — conjugate prior for the normal variance.

## D.2 Key Statistical Tables

### Normal Distribution Quantiles

| $p$ | $\Phi^{-1}(p)$ |
|---|---|
| 0.90 | 1.282 |
| 0.95 | 1.645 |
| 0.975 | 1.960 |
| 0.99 | 2.326 |
| 0.995 | 2.576 |

### ADF Critical Values (asymptotic, constant + trend)

| Significance | No trend | With trend |
|---|---|---|
| 1% | −3.43 | −3.96 |
| 5% | −2.86 | −3.41 |
| 10% | −2.57 | −3.12 |

### Johansen Trace Statistic Critical Values (5%)

| $r$ (null: rank $\leq r$) | $n=2$ | $n=3$ | $n=4$ |
|---|---|---|---|
| 0 | 15.5 | 29.7 | 47.2 |
| 1 | 3.8 | 15.4 | 29.7 |
| 2 | — | 3.8 | 15.4 |

## D.3 Moment-Generating Functions

| Distribution | MGF $M(t) = \mathbb{E}[e^{tX}]$ |
|---|---|
| $\mathcal{N}(\mu,\sigma^2)$ | $e^{\mu t + \sigma^2t^2/2}$ |
| $\text{Bernoulli}(p)$ | $1-p+pe^t$ |
| $\text{Poisson}(\lambda)$ | $e^{\lambda(e^t-1)}$ |
| $\text{Gamma}(\alpha,\beta)$ | $(1-t/\beta)^{-\alpha}$ ($t<\beta$) |
| $\text{Exponential}(\lambda)$ | $\lambda/(\lambda-t)$ ($t<\lambda$) |

**Key property:** $\mathbb{E}[X^k] = M^{(k)}(0)$ (the $k$-th derivative of MGF at 0).
