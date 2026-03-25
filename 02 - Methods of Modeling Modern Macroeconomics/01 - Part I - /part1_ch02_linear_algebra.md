# Chapter 2: Linear Algebra for Macroeconomic Systems

*Input-Output Matrices, Eigenvalues, and Economic Structure*

> *"The beauty of mathematics is that the same abstract structure — a system of linear equations — describes both the flow of goods between industries and the stability of a business cycle model."*

**Cross-reference:** *Principles* Ch. 4 (circular flow, input-output analysis); Ch. 6 (data and econometric methods); Ch. 9 (IS–LM as a simultaneous system); Ch. 27 (RBC model, matrix solution) **[P:Ch.4, P:Ch.6, P:Ch.9, P:Ch.27]**

---

## 2.1 Why Linear Algebra? Economic Systems as Matrix Problems

A recurring structure in macroeconomics is the **simultaneous equation system**: a set of conditions that must hold jointly, involving multiple endogenous variables that are all determined together. The IS–LM model [P:Ch.9] is a 2×2 system determining $(Y^*, i^*)$ from the intersection of two curves. The input-output model [P:Ch.4] is an $n \times n$ system determining production levels across all sectors simultaneously. The log-linearized DSGE model of Part VII is a $k \times k$ system of expectational difference equations. In each case, the right language is linear algebra.

Beyond solving systems, linear algebra provides the tools for:

- **Stability analysis:** whether a dynamic system returns to equilibrium after a shock depends on the eigenvalues of its transition matrix — a central concept in Chapters 4, 14, and 28.
- **Dimensionality reduction:** eigendecomposition reveals which directions in the state space matter most, which underlies the spectral analysis of time series and the solution algorithms for DSGE models.
- **Efficient computation:** the matrix operations `+.×` (multiply) and `⌹` (divide/solve) in APL, and their analogues in NumPy and Julia's LinearAlgebra, are the workhorses of every numerical algorithm in this book.

This chapter develops the linear algebra toolkit with these applications in mind. Every definition is illustrated with an economic example; every theorem is used somewhere in the book.

---

## 2.2 Vectors and Matrices as Economic Objects

### 2.2.1 Economic Vectors

A **vector** $\mathbf{x} \in \mathbb{R}^n$ is an ordered list of $n$ real numbers. In macroeconomics, vectors appear as:

- **State vectors** in dynamic models: $\mathbf{s}_t = (k_t, A_t)$ in the RBC model, where $k_t$ is the capital stock and $A_t$ is the technology level.
- **Output vectors** in input-output analysis: $\mathbf{x} = (x_1, x_2, \ldots, x_n)'$ where $x_i$ is gross output of industry $i$.
- **Impulse vectors** in VAR models: $\mathbf{e}_j$ is the $j$-th standard basis vector, representing a unit shock to variable $j$.

**Definition 2.1 (Inner Product).** The **inner product** of two vectors $\mathbf{u}, \mathbf{v} \in \mathbb{R}^n$ is:
$$\langle \mathbf{u}, \mathbf{v} \rangle = \mathbf{u}'\mathbf{v} = \sum_{i=1}^n u_i v_i.$$

The inner product $\mathbf{p}'\mathbf{q}$ where $\mathbf{p}$ is a price vector and $\mathbf{q}$ is a quantity vector gives the value of the basket at those prices — the fundamental operation underlying both GDP measurement and the computation of price indices [P:Ch.3.1].

### 2.2.2 Matrices as Linear Transformations

A matrix $A \in \mathbb{R}^{m \times n}$ represents a **linear map** $T: \mathbb{R}^n \to \mathbb{R}^m$ via $T(\mathbf{x}) = A\mathbf{x}$. The action of $A$ on a vector transforms it — rotating, scaling, projecting — without any nonlinear distortion. This linearity is exactly why linearized DSGE models are tractable: every period's state vector is a linear function of the previous state and the current shocks.

**Definition 2.2 (Matrix Multiplication).** For $A \in \mathbb{R}^{m \times k}$ and $B \in \mathbb{R}^{k \times n}$, the product $C = AB \in \mathbb{R}^{m \times n}$ has elements:
$$c_{ij} = \sum_{l=1}^k a_{il} b_{lj} = \mathbf{a}_i' \mathbf{b}_j,$$
where $\mathbf{a}_i'$ is the $i$-th row of $A$ and $\mathbf{b}_j$ is the $j$-th column of $B$. Matrix multiplication is associative $(AB)C = A(BC)$ and distributive $A(B+C) = AB + AC$, but **not commutative** in general: $AB \neq BA$.

In APL, matrix multiplication is the inner product `A +.× B` — read as "sum of products of corresponding elements, by column." This primitive is the single most-used operation in the entire book.

```apl
⍝ APL — matrix multiplication
⎕IO←0 ⋄ ⎕ML←1

A ← 2 2 ⍴ 1 2 3 4     ⍝ 2×2 matrix: [[1 2][3 4]]
B ← 2 2 ⍴ 5 6 7 8     ⍝ 2×2 matrix: [[5 6][7 8]]
C ← A +.× B            ⍝ matrix product: [[19 22][43 50]]
```

---

## 2.3 Matrix Inversion and the Solution of Linear Systems

### 2.3.1 The Inverse Matrix

**Definition 2.3 (Matrix Inverse).** A square matrix $A \in \mathbb{R}^{n \times n}$ is **invertible** (or **nonsingular**) if there exists a matrix $A^{-1}$ such that $A^{-1}A = AA^{-1} = I_n$, where $I_n$ is the $n \times n$ identity matrix. $A$ is invertible if and only if $\det(A) \neq 0$, equivalently, if and only if $A$ has full rank $n$.

The system $A\mathbf{x} = \mathbf{b}$ has a unique solution $\mathbf{x}^* = A^{-1}\mathbf{b}$ when $A$ is invertible.

In APL, $A^{-1}\mathbf{b}$ is written `b ⌹ A`. The monadic form `⌹A` computes $A^{-1}$ directly. The `⌹` primitive (called "domino" or "matrix divide") calls LAPACK's LU solver internally, so it is both convenient and numerically stable.

```apl
⍝ APL — solving Ax = b using ⌹
A ← 2 2 ⍴ 2 1 5 3      ⍝ coefficient matrix
b ← 8 5                  ⍝ right-hand side (as a vector)
x ← b ⌹ A               ⍝ solution: x = A⁻¹b
x                        ⍝ should give 19 ¯30 (check: 2×19 + 1×(¯30) = 8 ✓)

⍝ Explicit inverse (use with care — solve directly when possible)
Ainv ← ⌹A
Ainv +.× b               ⍝ same result via explicit inverse
```

**Note on numerical practice:** Computing $A^{-1}$ explicitly and then multiplying is less numerically stable than solving $A\mathbf{x} = \mathbf{b}$ directly. Always prefer `b ⌹ A` over `(⌹A) +.× b` in production code.

### 2.3.2 Cramer's Rule

For small systems (2×2, 3×3), Cramer's rule gives explicit closed-form solutions that are useful for comparative statics.

**Theorem 2.1 (Cramer's Rule).** For the system $A\mathbf{x} = \mathbf{b}$ with $A \in \mathbb{R}^{n \times n}$ invertible, the $i$-th component of the solution is:
$$x_i^* = \frac{\det(A_i)}{\det(A)},$$
where $A_i$ is the matrix obtained by replacing the $i$-th column of $A$ with $\mathbf{b}$.

*Proof sketch for 2×2.* With $A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$, $\mathbf{b} = \begin{pmatrix} e \\ f \end{pmatrix}$:
$$\det(A) = ad - bc, \quad \det(A_1) = \begin{vmatrix} e & b \\ f & d \end{vmatrix} = ed - bf, \quad \det(A_2) = \begin{vmatrix} a & e \\ c & f \end{vmatrix} = af - ce.$$
Direct calculation: $x_1^* = (ed - bf)/(ad - bc)$, $x_2^* = (af - ce)/(ad - bc)$, which one can verify satisfies $A\mathbf{x}^* = \mathbf{b}$. $\square$

---

## 2.4 The Leontief Input-Output Model

The input-output model of Wassily Leontief (1941) is both an important economic framework and an elegant illustration of linear algebra in macroeconomics. It is the foundation of the production-side national accounts [P:Ch.4.5].

**Definition 2.4 (Technical Coefficient Matrix).** For an economy with $n$ industries, the **technical coefficient matrix** $A = [a_{ij}]$ has element $a_{ij}$ equal to the dollar value of input from industry $i$ required to produce one dollar of gross output in industry $j$. Each column $j$ represents the production recipe of industry $j$.

The accounting identity for gross output: each industry's output equals its deliveries to other industries plus its deliveries to final demand $\mathbf{d}$:
$$\mathbf{x} = A\mathbf{x} + \mathbf{d} \implies (I - A)\mathbf{x} = \mathbf{d} \implies \mathbf{x} = (I-A)^{-1}\mathbf{d}.$$

**Definition 2.5 (Leontief Inverse).** The matrix $L = (I - A)^{-1}$ is the **Leontief inverse** or **total requirements matrix**. Its $(i,j)$ element gives the total output of industry $i$ — direct plus all indirect upstream requirements — needed to deliver one dollar of final demand for industry $j$'s product.

**Theorem 2.2 (Existence of the Leontief Inverse).** $(I - A)^{-1}$ exists and has all non-negative elements if and only if all eigenvalues of $A$ have modulus strictly less than 1. When this holds:
$$(I-A)^{-1} = I + A + A^2 + A^3 + \cdots = \sum_{k=0}^\infty A^k.$$

*Proof of the series representation.* For any matrix $B$ with all eigenvalues inside the unit circle, $B^k \to \mathbf{0}$ as $k \to \infty$. The partial sums $S_N = \sum_{k=0}^N A^k$ satisfy $(I-A)S_N = I - A^{N+1} \to I$ as $N \to \infty$. The limit is therefore $(I-A)^{-1}$. $\square$

The economic interpretation of the series is the **multiplier chain**: the first term $I$ represents direct demand; $A\mathbf{d}$ captures first-round intermediate input requirements; $A^2\mathbf{d}$ captures second-round requirements; and so on. This is precisely the input-output analogue of the Keynesian spending multiplier [P:Ch.8].

In APL, the Leontief inverse of a 3-sector economy:

```apl
⍝ APL — Leontief inverse
⎕IO←0 ⋄ ⎕ML←1

⍝ Technical coefficient matrix (3 sectors)
A ← 3 3 ⍴ 0.1 0.2 0.0
           0.3 0.1 0.2
           0.0 0.2 0.1

⍝ Identity matrix
I ← (=⍨ ⍳ ≢A)          ⍝ =⍨⍳n generates n×n identity: compare ⍳n with itself

⍝ Leontief inverse: (I-A)⁻¹ = ⌹(I-A)
L ← ⌹ I - A

⍝ Final demand vector
d ← 100 80 60

⍝ Gross output required
x ← L +.× d
x    ⍝ total output by sector

⍝ Verify: (I-A)x = d
(I - A) +.× x    ⍝ should equal d
```

The APL expression `=⍨⍳n` deserves explanation: `⍳n` generates the vector `0 1 2 ... n-1` (with `⎕IO←0`); `=⍨` applies the outer product of equality to this vector with itself, yielding an $n \times n$ identity matrix. This is a characteristic APL idiom — generating structured matrices from primitive operations.

---

## 2.5 Determinants, Rank, and Singularity

### 2.5.1 The Determinant

**Definition 2.6 (Determinant).** The **determinant** $\det(A)$ of a square matrix $A \in \mathbb{R}^{n \times n}$ is a scalar that measures the signed volume scaling factor of the linear transformation represented by $A$. For $n = 2$:
$$\det\begin{pmatrix} a & b \\ c & d \end{pmatrix} = ad - bc.$$

The determinant has several key properties:
- $\det(AB) = \det(A)\det(B)$.
- $\det(A^{-1}) = 1/\det(A)$ when $A$ is invertible.
- $\det(A') = \det(A)$.
- $A$ is invertible iff $\det(A) \neq 0$.

In macroeconomics, $\det(A)$ appears in Cramer's rule (Section 2.3.2) and in the IS–LM multiplier formulas of Chapter 6. The sign of $\det(I - A)$ in the Leontief model determines whether the economy is productive.

### 2.5.2 Rank and Linear Independence

**Definition 2.7 (Rank).** The **rank** of a matrix $A$, denoted $\text{rank}(A)$, is the dimension of the column space of $A$ — the number of linearly independent columns. Equivalently, it is the number of nonzero singular values of $A$.

For the DSGE model identification problem (Chapter 41), rank conditions are central: a model is identified only if the Jacobian matrix mapping structural parameters to model predictions has full column rank.

---

## 2.6 Eigenvalues and Eigenvectors

Eigenvalues and eigenvectors are the most important concepts in linear algebra for dynamic macroeconomics. The stability of any linear dynamic system — whether it converges to a steady state, oscillates, or explodes — is determined entirely by the eigenvalues of its transition matrix.

**Definition 2.8 (Eigenvalue and Eigenvector).** Let $A \in \mathbb{R}^{n \times n}$. A scalar $\lambda \in \mathbb{C}$ is an **eigenvalue** of $A$ with corresponding **eigenvector** $\mathbf{v} \neq \mathbf{0}$ if:
$$A\mathbf{v} = \lambda\mathbf{v}.$$

The set of all eigenvalues is the **spectrum** of $A$, denoted $\sigma(A)$. The eigenvalues are the roots of the **characteristic polynomial**:
$$p(\lambda) = \det(\lambda I - A) = 0.$$

For $n = 2$, this gives the quadratic $\lambda^2 - \text{tr}(A)\lambda + \det(A) = 0$, with roots:
$$\lambda_{1,2} = \frac{\text{tr}(A) \pm \sqrt{\text{tr}(A)^2 - 4\det(A)}}{2}.$$

Note two useful identities: $\lambda_1 + \lambda_2 = \text{tr}(A)$ and $\lambda_1 \lambda_2 = \det(A)$.

### 2.6.1 Eigenvalues and Stability

**Theorem 2.3 (Stability of Linear Discrete-Time Systems).** The system $\mathbf{x}_{t+1} = A\mathbf{x}_t$ converges to $\mathbf{0}$ from any initial condition if and only if all eigenvalues of $A$ satisfy $|\lambda_i| < 1$. The system diverges if any $|\lambda_i| > 1$.

**Theorem 2.4 (Stability of Linear Continuous-Time Systems).** The system $\dot{\mathbf{x}} = A\mathbf{x}$ converges to $\mathbf{0}$ from any initial condition if and only if all eigenvalues of $A$ satisfy $\text{Re}(\lambda_i) < 0$.

These two theorems underlie all stability analysis in macroeconomics:
- The Solow model's convergence to the steady state (Chapter 10) corresponds to the continuous-time condition $\text{Re}(\lambda) < 0$ for the linearized ODE.
- The Blanchard–Kahn condition for a unique stable solution to a DSGE model (Chapter 28) requires counting eigenvalues inside vs. outside the unit circle.
- The determinacy condition $\phi_\pi > 1$ for the Taylor rule (Chapter 28, [P:Ch.23.1]) translates into a requirement on the eigenvalues of the NK system's transition matrix.

### 2.6.2 Diagonalization and Matrix Powers

If $A$ has $n$ linearly independent eigenvectors $\mathbf{v}_1, \ldots, \mathbf{v}_n$ with eigenvalues $\lambda_1, \ldots, \lambda_n$, form the matrix $P = [\mathbf{v}_1 | \cdots | \mathbf{v}_n]$. Then:
$$A = PDP^{-1}, \quad D = \text{diag}(\lambda_1, \ldots, \lambda_n).$$

Powers are trivial: $A^k = PD^kP^{-1}$ where $D^k = \text{diag}(\lambda_1^k, \ldots, \lambda_n^k)$.

This is the basis for computing impulse response functions in VAR and DSGE models: the IRF at horizon $h$ is $A^h\mathbf{e}_j$, which reduces to $PD^hP^{-1}\mathbf{e}_j$ — the $j$-th column of $PD^hP^{-1}$.

```apl
⍝ APL — eigenvalue decomposition and matrix powers
⍝ (Dyalog APL does not have a built-in eigendecomposition primitive;
⍝  use the power method for dominant eigenvalue, or call numpy via ⎕PY)

⎕IO←0 ⋄ ⎕ML←1

⍝ For IRFs: compute A^h directly via repeated multiplication
A ← 2 2 ⍴ 0.9 0.1 0.0 0.8   ⍝ transition matrix

⍝ Impulse response of variable 1 to shock in variable 1, horizons 0..19
shock ← 1 0                    ⍝ unit shock to variable 1
irf ← {(A⍣⍵) +.× shock} ¨ ⍳ 20

⍝ Stack into matrix (20×2) for plotting
⊃ irf    ⍝ mix: convert nested vector to matrix
```

The APL idiom `(A⍣h) +.× shock` computes $A^h \mathbf{s}$: `⍣h` applies the matrix multiply `h` times. With `¨` (each) this generates the full IRF sequence in one expression.

---

## 2.7 The Spectral Decomposition and Symmetric Matrices

Many matrices in macroeconomics are symmetric: covariance matrices, Hessians, and the matrices arising from certain DSGE structures. Symmetric matrices have a particularly clean spectral structure.

**Theorem 2.5 (Spectral Decomposition of Symmetric Matrices).** If $A \in \mathbb{R}^{n \times n}$ is symmetric, then:
1. All eigenvalues are real.
2. Eigenvectors corresponding to distinct eigenvalues are orthogonal.
3. $A$ has an orthogonal eigendecomposition $A = Q\Lambda Q'$ where $Q$ is orthogonal ($Q'Q = I$) and $\Lambda = \text{diag}(\lambda_1, \ldots, \lambda_n)$.

**Definition 2.9 (Quadratic Form).** For a symmetric matrix $A$ and vector $\mathbf{x}$, the expression $\mathbf{x}'A\mathbf{x}$ is a **quadratic form**. The sign of this form for all nonzero $\mathbf{x}$ is determined entirely by the signs of the eigenvalues of $A$ — which is why positive/negative definiteness of the Hessian determines whether a critical point is a minimum/maximum (Theorem 1.4).

In APL, the quadratic form $\mathbf{x}'A\mathbf{x}$:

```apl
⍝ APL — quadratic form x'Ax
quadratic ← {⍺ +.× ⍵ +.× ⍺}    ⍝ ⍺ is x, ⍵ is A: x +.× (A +.× x)
x ← 1 2
A ← 2 2 ⍴ 3 1 1 2              ⍝ positive definite (eigenvalues both positive)
x quadratic A                    ⍝ should be 14 > 0
```

---

## 2.8 The IS–LM Model as a Linear System

*Cross-reference: Principles Ch. 9 (IS–LM model)* **[P:Ch.9]**

The IS–LM model determines equilibrium output $Y^*$ and the nominal interest rate $i^*$. Using the linear specifications from *Principles* Ch. 9.3:

**IS curve:** $Y = \bar{A} - b_r i$ (where $\bar{A}$ captures autonomous spending and $b_r > 0$ is investment interest sensitivity)

**LM curve:** $M/P = kY - hi$ (where $k > 0$ is income elasticity of money demand and $h > 0$ is interest elasticity)

Writing as a linear system $A\mathbf{y} = \mathbf{b}$:

$$\underbrace{\begin{pmatrix} 1 & b_r \\ k & -h \end{pmatrix}}_{A} \underbrace{\begin{pmatrix} Y \\ i \end{pmatrix}}_{\mathbf{y}} = \underbrace{\begin{pmatrix} \bar{A} \\ M/P \end{pmatrix}}_{\mathbf{b}}.$$

The solution by Cramer's rule:

$$\det(A) = -h - b_r k,$$

$$Y^* = \frac{1}{\det(A)}\det\begin{pmatrix} \bar{A} & b_r \\ M/P & -h \end{pmatrix} = \frac{-h\bar{A} - b_r(M/P)}{-h - b_r k} = \frac{h\bar{A} + b_r(M/P)}{h + b_r k},$$

$$i^* = \frac{1}{\det(A)}\det\begin{pmatrix} 1 & \bar{A} \\ k & M/P \end{pmatrix} = \frac{M/P - k\bar{A}}{-h - b_r k} = \frac{k\bar{A} - M/P}{h + b_r k}.$$

The fiscal multiplier $\partial Y^*/\partial G$ (noting $\bar{A}$ increases one-for-one with $G$):

$$\frac{\partial Y^*}{\partial G} = \frac{h}{h + b_r k}.$$

This is the IS–LM fiscal multiplier derived in *Principles* Ch. 9.3 [P:Ch.9.3] — now from the matrix inverse rather than from graphical reasoning. Chapter 6 of this book develops this analysis fully, adding Cramer's rule derivations for all policy multipliers and extending to the open economy (Mundell–Fleming).

```apl
⍝ APL — IS-LM solution via ⌹
⎕IO←0 ⋄ ⎕ML←1

⍝ Parameters
br ← 2    ⍝ investment-interest sensitivity
k  ← 0.5  ⍝ income elasticity of money demand
h  ← 4    ⍝ interest elasticity of money demand

⍝ Build coefficient matrix
islm_matrix ← {br k h ← ⍵
    2 2 ⍴ 1 br k (-h)}     ⍝ [[1, br], [k, -h]]

A ← islm_matrix br k h

⍝ Exogenous variables: Abar = 200, M/P = 500
Abar ← 200
MP   ← 500
b    ← Abar MP

⍝ Solve for equilibrium
Y_star i_star ← b ⌹ A
Y_star    ⍝ equilibrium output
i_star    ⍝ equilibrium interest rate

⍝ Fiscal multiplier: ∂Y*/∂G = h / (h + br×k)
dY_dG ← h ÷ h + br × k     ⍝ should match second element of gradient
```

---

## 2.9 The Jordan Normal Form and Defective Matrices

Not every matrix is diagonalisable. When a matrix has repeated eigenvalues and insufficient eigenvectors, we need the **Jordan normal form**.

**Definition 2.10 (Jordan Block).** A **Jordan block** of size $m$ with eigenvalue $\lambda$ is the $m \times m$ matrix:
$$J_m(\lambda) = \begin{pmatrix} \lambda & 1 & 0 & \cdots & 0 \\ 0 & \lambda & 1 & \cdots & 0 \\ \vdots & & \ddots & \ddots & \vdots \\ 0 & \cdots & 0 & \lambda & 1 \\ 0 & \cdots & 0 & 0 & \lambda \end{pmatrix}.$$

**Theorem 2.6 (Jordan Normal Form).** Every square matrix $A \in \mathbb{R}^{n \times n}$ (over $\mathbb{C}$) is similar to a block-diagonal matrix $J = \text{diag}(J_{m_1}(\lambda_1), \ldots, J_{m_r}(\lambda_r))$, where $\sum m_i = n$.

For the stability analysis of dynamic systems, Jordan blocks with eigenvalue $|\lambda| < 1$ still converge to zero, but more slowly than the diagonal case — the presence of the superdiagonal 1s means the $k$-th power of $J_m(\lambda)$ contains terms like $\binom{k}{j}\lambda^{k-j}$, which still go to zero as $k \to \infty$ when $|\lambda| < 1$.

In practice, most matrices arising in macroeconomic models are either diagonalisable or can be handled with the generalized Schur (QZ) decomposition (Chapter 28), so the Jordan form is primarily of theoretical importance.

---

## 2.10 Worked Example: Three-Sector Leontief Economy

*Cross-reference: Principles Ch. 4.5 (input-output analysis)* **[P:Ch.4.5]**

Consider a three-sector economy (manufacturing, services, agriculture) with the following technical coefficient matrix:

$$A = \begin{pmatrix} 0.20 & 0.15 & 0.05 \\ 0.25 & 0.10 & 0.20 \\ 0.05 & 0.08 & 0.15 \end{pmatrix},$$

and final demand vector $\mathbf{d} = (120, 80, 50)'$.

**Step 1:** Verify $A$ is productive. Compute $\rho(A) = \max_i|\lambda_i|$, the spectral radius of $A$. If $\rho(A) < 1$, the Leontief inverse exists.

The eigenvalues of $A$ can be found from its characteristic polynomial. For a rough check, note that the maximum column sum of $A$ is $\max_j \sum_i a_{ij} = 0.50 + 0.33 + 0.40 = 0.50 < 1$. Since this column-sum norm bounds the spectral radius, $\rho(A) \leq 0.50 < 1$, confirming productivity.

**Step 2:** Compute the Leontief inverse $L = (I-A)^{-1}$.

$$I - A = \begin{pmatrix} 0.80 & -0.15 & -0.05 \\ -0.25 & 0.90 & -0.20 \\ -0.05 & -0.08 & 0.85 \end{pmatrix}.$$

Computing $(I-A)^{-1}$ numerically (shown below):

$$L \approx \begin{pmatrix} 1.327 & 0.247 & 0.123 \\ 0.404 & 1.214 & 0.328 \\ 0.094 & 0.122 & 1.222 \end{pmatrix}.$$

**Step 3:** Compute gross output $\mathbf{x} = L\mathbf{d}$.

$$\mathbf{x} = L\mathbf{d} \approx \begin{pmatrix} 1.327 & 0.247 & 0.123 \\ 0.404 & 1.214 & 0.328 \\ 0.094 & 0.122 & 1.222 \end{pmatrix} \begin{pmatrix} 120 \\ 80 \\ 50 \end{pmatrix} \approx \begin{pmatrix} 185.8 \\ 130.0 \\ 84.3 \end{pmatrix}.$$

**Step 4:** Interpretation. To deliver \$120 of manufacturing to final demand, the economy must produce \$185.8 of manufacturing gross output in total — the additional \$65.8 supplies intermediate inputs to all three sectors through the full chain of upstream requirements.

**Step 5:** Multiplier. The total output multiplier for manufacturing final demand is $\sum_i L_{i1} = 1.327 + 0.404 + 0.094 = 1.825$: one dollar of final demand for manufacturing generates \$1.825 of total gross output across all sectors.

```apl
⍝ APL — Three-sector Leontief model
⎕IO←0 ⋄ ⎕ML←1

A ← 3 3 ⍴ 0.20 0.15 0.05
           0.25 0.10 0.20
           0.05 0.08 0.15

d ← 120 80 50              ⍝ final demand

I3 ← =⍨ ⍳ 3               ⍝ 3×3 identity matrix
L  ← ⌹ I3 - A             ⍝ Leontief inverse

x ← L +.× d               ⍝ gross output
x                          ⍝ ≈ 185.8 130.0 84.3

⍝ Output multipliers (column sums of L)
+⌿ L                      ⍝ total output multiplier per dollar of final demand by sector

⍝ Verify: (I-A)x = d
(I3 - A) +.× x            ⍝ should recover d
```

---

## 2.11 Programming Exercises

### Exercise 2.1 (APL)

Implement a function `leontief ← {⌹ (=⍨⍳≢⍵) - ⍵}` that takes a technical coefficient matrix and returns the Leontief inverse in one line. Test it on the 3-sector example above. Then implement the full output multiplier calculation `multipliers ← {+⌿ leontief ⍵}` and verify that the multiplier for sector $j$ equals the $j$-th column sum of $L$.

```apl
⎕IO←0 ⋄ ⎕ML←1
leontief  ← {⌹ (=⍨⍳≢⍵) - ⍵}      ⍝ (I-A)⁻¹ in one token
multipliers ← {+⌿ leontief ⍵}      ⍝ column sums = output multipliers

A ← 3 3 ⍴ 0.20 0.15 0.05 0.25 0.10 0.20 0.05 0.08 0.15
multipliers A                        ⍝ total output multipliers by sector
```

### Exercise 2.2 (Python)

```python
import numpy as np

A = np.array([[0.20, 0.15, 0.05],
              [0.25, 0.10, 0.20],
              [0.05, 0.08, 0.15]])
d = np.array([120, 80, 50])

n = A.shape[0]
L = np.linalg.inv(np.eye(n) - A)   # Leontief inverse
x = L @ d                            # gross output
print("Gross output:", x.round(2))
print("Output multipliers:", L.sum(axis=0).round(4))
print("Spectral radius:", max(abs(np.linalg.eigvals(A))).round(4))
```

### Exercise 2.3 (Julia)

```julia
using LinearAlgebra
A = [0.20 0.15 0.05;
     0.25 0.10 0.20;
     0.05 0.08 0.15]
d = [120.0, 80.0, 50.0]
L = inv(I - A)
x = L * d
println("Gross output: ", round.(x, digits=2))
println("Multipliers: ", round.(sum(L, dims=1), digits=4))
println("Spectral radius: ", maximum(abs.(eigvals(A))) |> x -> round(x, digits=4))
```

### Exercise 2.4 (R)

```r
A <- matrix(c(0.20,0.25,0.05, 0.15,0.10,0.08, 0.05,0.20,0.15), 3, 3)
d <- c(120, 80, 50)
n <- nrow(A)
L <- solve(diag(n) - A)
x <- L %*% d
cat("Gross output:", round(x, 2), "\n")
cat("Multipliers:", round(colSums(L), 4), "\n")
cat("Spectral radius:", round(max(abs(eigen(A)$values)), 4), "\n")
```

### Exercise 2.5 — IS–LM Parameter Sweep ($\star$)

Using the IS–LM matrix system from Section 2.8, write an APL dfn `islm_multipliers ← {br k h ← ⍵ ⋄ ...}` that returns the fiscal and monetary multipliers $(\partial Y^*/\partial G, \partial Y^*/\partial(M/P))$ for given parameters. Generate a 10×10 grid of fiscal multipliers over $(b_r, h) \in [0.5, 4] \times [1, 8]$ using `∘.f` outer product syntax and plot as a heat map.

### Exercise 2.6 — Eigenvalue Stability ($\star$)

For the 2×2 transition matrix $A = \begin{pmatrix} 0.9 & \phi \\ 0 & 0.8 \end{pmatrix}$, find the range of $\phi$ for which the system $\mathbf{x}_{t+1} = A\mathbf{x}_t$ is stable (all eigenvalues inside the unit circle). Note: the eigenvalues of a triangular matrix are its diagonal entries, so stability holds for any $\phi$. Now perturb to $A = \begin{pmatrix} 0.9 & 0 \\ \phi & 0.8 \end{pmatrix}$ and repeat. What does this tell you about how off-diagonal elements affect stability?

### Exercise 2.7 — Perron–Frobenius ($\star\star$)

The Perron–Frobenius theorem states that a non-negative irreducible matrix $A$ has a unique largest real eigenvalue $\lambda_{PF} > 0$ (the Perron root) with a corresponding non-negative eigenvector. In the Leontief context, if $\lambda_{PF} < 1$, the economy is productive. (a) Verify the Perron root numerically for the 3-sector matrix $A$ above. (b) Find by bisection the largest value of a scalar multiplier $\alpha$ such that $\alpha A$ remains productive (i.e., $\lambda_{PF}(\alpha A) < 1$). (c) Interpret economically.

---

## 2.12 Chapter Summary

This chapter developed the linear algebra toolkit for macroeconomic modeling.

**Key results:**

- The **Leontief inverse** $(I-A)^{-1} = \sum_{k=0}^\infty A^k$ exists when the spectral radius $\rho(A) < 1$ and gives total (direct plus indirect) output requirements per unit of final demand.
- **Cramer's rule** provides explicit closed-form solutions for small linear systems — the foundation for IS–LM multiplier derivations in Chapter 6.
- **Eigenvalues** determine stability: a discrete system $\mathbf{x}_{t+1} = A\mathbf{x}_t$ is stable iff all $|\lambda_i| < 1$; a continuous system $\dot{\mathbf{x}} = A\mathbf{x}$ is stable iff all $\text{Re}(\lambda_i) < 0$.
- **Diagonalization** $A = PDP^{-1}$ makes powers $A^k = PD^kP^{-1}$ trivial, enabling analytical and numerical computation of impulse response functions.
- In APL: `⌹` solves linear systems and computes matrix inverses; `+.×` is matrix multiplication; `=⍨⍳n` generates the identity matrix; `+⌿` computes column sums (output multipliers).

**Connections forward:** Chapter 3 uses eigenvalue analysis to classify the equilibria of differential equation systems. Chapter 4 applies it to difference equations and previews the Blanchard–Kahn condition. Chapter 6 uses Cramer's rule to derive all IS–LM multipliers. Chapter 28 uses the generalized Schur (QZ) decomposition — an extension of eigendecomposition — to solve linear DSGE models.

---

*Next: Chapter 3 — Differential Equations in Continuous-Time Macro Models*
