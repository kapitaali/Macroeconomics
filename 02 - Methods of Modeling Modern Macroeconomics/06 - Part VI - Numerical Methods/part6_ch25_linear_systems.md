# Chapter 25: Solving Systems of Equations

*From Linear Systems to Sparse Matrices*

> *"Every time you call a linear solver, LU decomposition is happening — whether you know it or not."*

**Cross-reference:** *Principles* Ch. 4 (Leontief input-output, matrix inversion); Ch. 9 (IS–LM as linear system); Ch. 27 (DSGE equilibrium conditions as large sparse system) **[P:Ch.4, P:Ch.9, P:Ch.27]**

---

## 25.1 The Ubiquity of Linear Systems in Macroeconomics

Linear systems $A\mathbf{x} = \mathbf{b}$ appear at every stage of the macroeconomic modeling workflow:

- **IS–LM:** The $2\times2$ system from Chapter 6, solved by `b ⌹ A` in APL.
- **Leontief inverse:** $(I-A)^{-1}\mathbf{d}$, a dense $n\times n$ system for $n$ sectors.
- **Kalman filter:** Prediction step $P_{\text{pred}} = FPF' + Q$ requires matrix multiplication; the innovation variance $F_t = HPH' + R$ and gain $K = PH'F_t^{-1}$ require matrix inversion.
- **DSGE log-linearization:** The linearized first-order conditions form a large sparse system; the Sylvester equation for the MSV solution is a vectorized linear system.
- **Value function iteration:** The Howard policy improvement step solves $(I-\beta P)V = u$, where $P$ is the $N\times N$ Markov transition matrix (sparse).

Understanding how linear systems are solved — not just that they are — enables better debugging, performance tuning, and recognition of when standard methods will fail.

---

## 25.2 LU Decomposition: Direct Method for Dense Systems

**Definition 25.1 (LU Decomposition).** The **LU decomposition** of an invertible matrix $A \in \mathbb{R}^{n\times n}$ is the factorization $PA = LU$, where:
- $P$ is a permutation matrix (row reordering for numerical stability).
- $L$ is unit lower triangular: $l_{ii} = 1$, $l_{ij} = 0$ for $j > i$.
- $U$ is upper triangular: $u_{ij} = 0$ for $j < i$.

**Algorithm 25.1 (LU Decomposition with Partial Pivoting).**

For $k = 1, \ldots, n-1$:
1. Find the pivot: $p = \arg\max_{i\geq k}|a_{ik}|$ (largest absolute value in column $k$ at or below the diagonal).
2. Swap rows $k$ and $p$ in $A$.
3. For $i = k+1, \ldots, n$: $l_{ik} = a_{ik}/a_{kk}$ (elimination multipliers).
4. For $i = k+1, \ldots, n$: $a_{i,k:n} \leftarrow a_{i,k:n} - l_{ik}\cdot a_{k,k:n}$ (row elimination).

**Theorem 25.1 (LU Complexity).** The LU decomposition of an $n\times n$ matrix requires $O(n^3/3)$ flops for elimination plus $O(n^2)$ for forward/back substitution to solve $A\mathbf{x} = \mathbf{b}$.

*Proof.* In the $k$-th elimination step, the inner loop performs $(n-k)^2$ operations. Total: $\sum_{k=1}^{n-1}(n-k)^2 \approx n^3/3$. $\square$

**Solving $A\mathbf{x} = \mathbf{b}$ via LU:**

1. Factor: $PA = LU$ (cost $O(n^3/3)$).
2. Forward substitution: solve $L\mathbf{y} = P\mathbf{b}$ for $\mathbf{y}$ (cost $O(n^2)$).
3. Back substitution: solve $U\mathbf{x} = \mathbf{y}$ for $\mathbf{x}$ (cost $O(n^2)$).

For multiple right-hand sides $\mathbf{b}_1, \ldots, \mathbf{b}_k$ with the same $A$: factor once ($O(n^3)$), then solve each in $O(n^2)$. This is the key efficiency insight for the Kalman filter: the innovation covariance $F_t$ changes at each time step, so each inversion is a separate $O(p^3)$ cost; but if $F_t$ is constant (stationary model), it is factored once.

**Condition number and numerical stability:**

**Definition 25.2 (Condition Number).** The **condition number** of a matrix $A$ is:

$$\kappa(A) = \|A\|\cdot\|A^{-1}\|.$$

The relative error in the computed solution $\hat{\mathbf{x}}$ satisfies:

$$\frac{\|\hat{\mathbf{x}} - \mathbf{x}^*\|}{\|\mathbf{x}^*\|} \leq \kappa(A)\cdot\varepsilon_M + O(\varepsilon_M^2),$$

where $\varepsilon_M \approx 10^{-16}$ is machine epsilon. A well-conditioned matrix ($\kappa(A) \approx 1$) introduces negligible numerical error; a poorly conditioned matrix ($\kappa(A) \sim 10^{12}$) can lose $12$ digits of accuracy in the solution.

**Partial pivoting** (choosing the largest pivot) reduces condition number amplification: it ensures that the elements of $L$ satisfy $|l_{ij}| \leq 1$, bounding the growth of rounding errors.

In APL, `⌹A` calls LAPACK's `DGETRF/DGETRS` (LU with partial pivoting) internally. The solution `b ⌹ A` is equivalent to LU-factoring $A$ and solving — but APL's `⌹` handles the complete pipeline automatically.

---

## 25.3 Iterative Methods for Large Systems

When $n$ is very large (hundreds of thousands of variables in heterogeneous-agent models), direct LU is impractical ($O(n^3)$ is too slow). **Iterative methods** compute approximate solutions by successive refinement.

### 25.3.1 Gauss–Seidel

**Algorithm 25.2 (Gauss–Seidel).**

For a system $A\mathbf{x} = \mathbf{b}$, decompose $A = D + L + U$ (diagonal, strict lower, strict upper triangular).

Initialize $\mathbf{x}^{(0)}$. For $k = 0, 1, 2, \ldots$:

$$x_i^{(k+1)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j<i}a_{ij}x_j^{(k+1)} - \sum_{j>i}a_{ij}x_j^{(k)}\right), \quad i = 1, \ldots, n.$$

Each element is updated using the most recent values of all other elements.

**Convergence:** Gauss–Seidel converges (to the true solution) iff the spectral radius $\rho(-(D+L)^{-1}U) < 1$. A sufficient condition is that $A$ is strictly diagonally dominant: $|a_{ii}| > \sum_{j\neq i}|a_{ij}|$ for all $i$.

### 25.3.2 Conjugate Gradient for SPD Systems

For **symmetric positive definite (SPD)** systems — which arise in the normal equations of least squares, the Howard policy improvement step with symmetric $P$, and certain DSGE Hessians — the **Conjugate Gradient (CG)** method achieves faster convergence.

**Algorithm 25.3 (Conjugate Gradient).**

Initialize $\mathbf{x}^{(0)} = \mathbf{0}$, $\mathbf{r}^{(0)} = \mathbf{b}$, $\mathbf{p}^{(0)} = \mathbf{r}^{(0)}$.

For $k = 0, 1, 2, \ldots$:
1. $\alpha_k = (\mathbf{r}^{(k)})'(\mathbf{r}^{(k)}) / [(\mathbf{p}^{(k)})'A\mathbf{p}^{(k)}]$
2. $\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + \alpha_k\mathbf{p}^{(k)}$
3. $\mathbf{r}^{(k+1)} = \mathbf{r}^{(k)} - \alpha_kA\mathbf{p}^{(k)}$
4. $\beta_k = (\mathbf{r}^{(k+1)})'(\mathbf{r}^{(k+1)}) / [(\mathbf{r}^{(k)})'(\mathbf{r}^{(k)})]$
5. $\mathbf{p}^{(k+1)} = \mathbf{r}^{(k+1)} + \beta_k\mathbf{p}^{(k)}$

**Theorem 25.2 (CG Convergence).** CG converges to the exact solution in at most $n$ steps. The error after $k$ steps satisfies:

$$\|\mathbf{x}^* - \mathbf{x}^{(k)}\|_A \leq 2\left(\frac{\sqrt{\kappa}-1}{\sqrt{\kappa}+1}\right)^k\|\mathbf{x}^*-\mathbf{x}^{(0)}\|_A,$$

where $\|\mathbf{v}\|_A = \sqrt{\mathbf{v}'A\mathbf{v}}$ and $\kappa = \kappa(A)$.

---

## 25.4 Sparse Matrix Methods

Many macroeconomic linear systems are **sparse**: most elements are zero, and the non-zero structure has economic meaning.

**Leontief model:** The $n\times n$ technical coefficient matrix $A$ for an $n$-sector economy has at most $n^2$ entries, but in practice each sector uses inputs from only a few others — $A$ is sparse. For $n = 500$ sectors, storing $A$ densely requires $500^2 = 250{,}000$ entries; storing sparsely requires only the non-zero entries (perhaps $5{,}000$–$10{,}000$).

**DSGE Jacobian:** The gradient and Hessian of the DSGE log-likelihood have a block structure corresponding to the model's timing restrictions — many parameters affect only specific equations and specific periods.

**Sparse storage formats:**

- **CSR (Compressed Sparse Row):** Store all non-zero values, their column indices, and row start pointers. Matrix-vector products $A\mathbf{v}$ cost $O(\text{nnz})$ where $\text{nnz}$ is the number of non-zeros.
- **CSC (Compressed Sparse Column):** Same structure transposed. Preferred for column-oriented operations.
- **COO (Coordinate format):** List of $(i, j, v)$ triples. Convenient for construction, converted to CSR/CSC for computation.

In APL, sparsity is simulated via Boolean masking:

```apl
⍝ APL — Sparse Leontief computation via Boolean masking
⎕IO←0 ⋄ ⎕ML←1

⍝ For large sparse A: mask selects non-zero entries
⍝ A_sparse = A × (|A > threshold)  — zero out small entries
sparse_mask ← {1e¯10 < |⍵}              ⍝ Boolean: 1 where non-zero
A_sparse    ← {⍵ × sparse_mask ⍵}      ⍝ retain only significant entries

⍝ Leontief inverse: still computed via ⌹ (LAPACK handles sparsity internally)
leontief_sparse ← {⌹ (=⍨⍳≢⍵) - ⍵}    ⍝ same formula; LAPACK exploits zeros

⍝ For very large systems: use compressed column storage
⍝ via ⎕PY interface to scipy.sparse
⎕PY.Import 'scipy.sparse as sp'
⎕PY.Import 'scipy.sparse.linalg as sla'

⍝ Build sparse matrix in Python from APL arrays
make_sparse ← {data rows cols shape ← ⍵
    ⎕PY.Call 'sp.csr_matrix' ((⊂data)(⊂rows, ⍨cols) shape)}
```

---

## 25.5 Overdetermined Systems: Least Squares

Many estimation problems yield **overdetermined** systems $A\mathbf{x} \approx \mathbf{b}$ with $m > n$ equations and $n$ unknowns. The OLS solution minimizes $\|\mathbf{b} - A\mathbf{x}\|_2^2$.

### 25.5.1 Normal Equations

The OLS solution satisfies the **normal equations**:

$$A'A\mathbf{x}^* = A'\mathbf{b} \implies \mathbf{x}^* = (A'A)^{-1}A'\mathbf{b}.$$

In APL: `x_ols ← (⌹ X) +.× y` — already derived in Chapter 19 for VAR estimation. This uses the pseudo-inverse $A^+ = (A'A)^{-1}A'$.

### 25.5.2 QR Decomposition

The QR decomposition is numerically superior to forming $A'A$ (which squares the condition number):

**Definition 25.3 (QR Decomposition).** For $A \in \mathbb{R}^{m\times n}$ ($m \geq n$, full column rank): $A = QR$ where $Q \in \mathbb{R}^{m\times n}$ is orthonormal ($Q'Q = I_n$) and $R \in \mathbb{R}^{n\times n}$ is upper triangular.

The OLS solution via QR: $A\mathbf{x}^* = \mathbf{b} \Rightarrow QR\mathbf{x}^* = \mathbf{b} \Rightarrow R\mathbf{x}^* = Q'\mathbf{b}$, solved by back substitution.

**Advantage of QR:** $\kappa(R) = \kappa(A)$, whereas $\kappa(A'A) = \kappa(A)^2$. For a moderately ill-conditioned system with $\kappa(A) = 10^8$, normal equations have $\kappa(A'A) = 10^{16}$ — at the limit of double precision. QR avoids this squaring.

---

## 25.6 Worked Example: 100-Sector Leontief Economy

*Cross-reference: Principles Ch. 4.5 (input-output analysis)* **[P:Ch.4.5]**

Consider a 100-sector U.S. economy with technical coefficient matrix from the 2019 BEA input-output tables. The system $(I-A)\mathbf{x} = \mathbf{d}$ (gross output from final demand) is a dense $100\times100$ linear system.

**Performance comparison:**

| Method | Cost | Time (100×100) | Notes |
|---|---|---|---|
| Dense LU (`⌹` in APL) | $O(n^3)$ | $< 1$ ms | Direct; exact to machine precision |
| Dense inverse `⌹A` then `+.×d` | $O(n^3)$ | $< 1$ ms | Less stable; avoid when possible |
| Gauss–Seidel | $O(n^2)$ per iter | ~50 ms (100 iters) | Only for diag-dominant; slow here |
| Sparse LU (scipy) | $O(\text{nnz}^{1.5})$ | $< 0.1$ ms | Much faster for sparse $A$ |

For $n = 500$ sectors: dense LU takes $O(500^3) \approx 1.25\times10^8$ flops, approximately $0.1$ seconds. Sparse LU for a matrix with $\text{nnz} = 5{,}000$ takes $O(5000^{1.5}) \approx 3.5\times10^5$ flops — nearly 400x faster.

```python
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import time

np.random.seed(42)
n = 100

# Dense random technical coefficient matrix (sparse in practice)
A_dense = np.random.uniform(0, 0.1, (n,n))  # small coefficients
# Ensure productive: column sums < 1
A_dense /= (A_dense.sum(axis=0) + 0.5)
d = np.ones(n) * 10  # final demand

# Dense solution
t0 = time.perf_counter()
x_dense = np.linalg.solve(np.eye(n) - A_dense, d)
t_dense = time.perf_counter() - t0
print(f"Dense LU:    {t_dense*1000:.2f} ms, max residual: {np.max(np.abs((np.eye(n)-A_dense)@x_dense - d)):.2e}")

# Sparse solution (keep only entries > 0.05)
A_sparse = A_dense.copy(); A_sparse[A_sparse < 0.05] = 0
A_sp = sp.csr_matrix(np.eye(n) - A_sparse)
t0 = time.perf_counter()
x_sparse = sla.spsolve(A_sp, d)
t_sparse = time.perf_counter() - t0
nnz = A_sp.nnz
print(f"Sparse LU:   {t_sparse*1000:.2f} ms, nnz={nnz}/{n**2} ({100*nnz/n**2:.1f}%), max residual: {np.max(np.abs(A_sp@x_sparse - d)):.2e}")

# Output multipliers (column sums of Leontief inverse)
L = np.linalg.inv(np.eye(n) - A_dense)
multipliers = L.sum(axis=0)
print(f"\nOutput multipliers: mean={multipliers.mean():.3f}, max={multipliers.max():.3f}, min={multipliers.min():.3f}")
```

```julia
using LinearAlgebra, SparseArrays

n = 100; Random.seed!(42)
A = rand(n,n) .* 0.1; A ./= (sum(A,dims=1) .+ 0.5)
d = ones(n) .* 10.0
I_n = I(n)

# Dense
t_dense = @elapsed x_dense = (I_n - A) \ d
println("Dense: $(round(t_dense*1000,digits=2)) ms, residual: $(maximum(abs.((I_n-A)*x_dense-d)))")

# Sparse
A_sp = A .* (A .> 0.05)
A_mat = sparse(I_n - A_sp)
t_sparse = @elapsed x_sparse = A_mat \ d
println("Sparse: $(round(t_sparse*1000,digits=2)) ms, nnz=$(nnz(A_mat))")
```

```r
n <- 100; set.seed(42)
A <- matrix(runif(n*n,0,0.1),n,n); A <- A/outer(rep(1,n),colSums(A)+0.5)
d <- rep(10,n); I_n <- diag(n)

# Dense
t_dense <- system.time(x_dense <- solve(I_n-A, d))[3]
cat(sprintf("Dense: %.2f ms, residual: %.2e\n", t_dense*1000, max(abs((I_n-A)%*%x_dense-d))))

# Sparse (using Matrix package)
library(Matrix)
A_sp <- A; A_sp[A_sp < 0.05] <- 0
M_sp <- as(I_n-A_sp, "sparseMatrix")
t_sparse <- system.time(x_sparse <- solve(M_sp, d))[3]
cat(sprintf("Sparse: %.2f ms, nnz=%d\n", t_sparse*1000, nnz(M_sp)))
```

---

## 25.7 Programming Exercises

### Exercise 25.1 (APL — Condition Number Check)

Write a dfn `check_system ← {A b ← ⍵ ⋄ ...}` that: (a) computes the condition number $\kappa(A)$ via `⌹` (approximated as $\|A\|_\infty\cdot\|A^{-1}\|_\infty$ where the $\infty$-norm is `⌈/+/|M`); (b) computes the solution `x ← b ⌹ A`; (c) computes the backward error `\|(Ax-b)/b\|_\infty`; (d) warns if `kappa(A) × machine_epsilon > 1e-8`. Test on the IS-LM system from Chapter 6 and on a near-singular matrix.

### Exercise 25.2 (Python — Howard Improvement as Sparse Linear System)

In the buffer-stock consumption model (Chapter 15), the Howard policy improvement step solves $(I-\beta P_{c^n})V^n = u_{c^n}$ where $P_{c^n}$ is the $N\times N$ sparse Markov transition matrix under policy $c^n$. (a) Build $P_{c^n}$ as a scipy sparse matrix (each row has at most 2 non-zeros for linear interpolation). (b) Solve using `scipy.sparse.linalg.spsolve`. (c) Compare the number of VFI iterations to PFI iterations (using Howard improvement). Verify PFI converges 10–20× faster.

### Exercise 25.3 (Julia — QR vs. Normal Equations)

```julia
using LinearAlgebra

# Demonstrate numerical superiority of QR over normal equations
function test_ols(m, n, cond_number)
    # Generate a matrix with specified condition number
    U,_,V = svd(randn(m,n))
    s = range(1, cond_number, length=n)
    A = U[:,1:n] * diagm(s) * V'
    x_true = randn(n)
    b = A*x_true + 0.01*randn(m)
    
    # Normal equations
    x_ne = (A'*A) \ (A'*b)
    err_ne = norm(x_ne - x_true)
    
    # QR
    Q,R = qr(A)
    x_qr = R \ (Matrix(Q)'*b)
    err_qr = norm(x_qr - x_true)
    
    return err_ne, err_qr, cond(A'*A)
end

println("Condition number comparison:")
for κ in [1e4, 1e6, 1e8, 1e10]
    err_ne, err_qr, cond_AtA = test_ols(200, 50, κ)
    println("  κ(A)=$(Int(κ)): κ(A'A)=$(round(cond_AtA,sigdigits=2)), " *
            "err(NE)=$(round(err_ne,digits=6)), err(QR)=$(round(err_qr,digits=6))")
end
```

### Exercise 25.4 — Sparse DSGE Jacobian ($\star$)

The DSGE equilibrium conditions form a system $\mathbf{F}(\mathbf{x}, \mathbf{x}_{-1}, \mathbf{x}_{+1}, \bm\varepsilon) = \mathbf{0}$ with $n_y$ endogenous variables and $n_z$ exogenous shocks. (a) For the 3-equation NK model, identify the $6\times 6$ Jacobian $\partial\mathbf{F}/\partial(\mathbf{x}, \mathbf{x}_{-1}, \mathbf{x}_{+1})$. (b) Show the Jacobian has a block structure with many zeros — quantify the sparsity ratio (zeros/total). (c) For a medium-scale model with 20 equations, estimate the memory savings from sparse storage.

---

## 25.8 Chapter Summary

**Key results:**

- **LU decomposition** $PA = LU$ solves $A\mathbf{x} = \mathbf{b}$ in $O(n^3/3) + O(n^2)$ flops; APL's `⌹` calls LAPACK internally.
- The **condition number** $\kappa(A)$ bounds the relative solution error: $\|\hat{\mathbf{x}}-\mathbf{x}^*\|/\|\mathbf{x}^*\| \leq \kappa(A)\varepsilon_M$; partial pivoting controls $\kappa$ amplification.
- **Gauss–Seidel** converges for diagonally dominant systems; **conjugate gradient** achieves error $\leq 2[(\sqrt\kappa-1)/(\sqrt\kappa+1)]^k$ for SPD systems — much faster than Gauss–Seidel for well-conditioned problems.
- **Sparse storage** (CSR/CSC) reduces memory and matrix-vector cost from $O(n^2)$ to $O(\text{nnz})$; sparse direct solvers cost $O(\text{nnz}^{1.5})$ vs. dense $O(n^3)$.
- **QR decomposition** for OLS is numerically superior to normal equations: $\kappa(R) = \kappa(A)$ vs. $\kappa(A'A) = \kappa(A)^2$.
- In APL: `b ⌹ A` for dense systems; for sparse, use `⎕PY` to call scipy.sparse; the Leontief inverse is `⌹(I-A)` regardless of sparsity (LAPACK handles zeros efficiently).

*Next: Chapter 26 — Monte Carlo Methods: Simulating Macroeconomic Models Under Uncertainty*
