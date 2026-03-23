# Part V: Stochastic Methods for Macroeconomic Modeling

*Connects to: Principles Ch. 6, Ch. 15–16, Ch. 20, Appendix B*

---

The models of Parts III and IV are deterministic or, in the case of the RBC model, stochastic in a fully numerical sense — we solved the Bellman equation on a grid without exploiting any linear structure. This part develops the analytical and statistical machinery for handling stochastic models *algebraically*: rational expectations solution methods, time-series econometrics, the Kalman filter, and structural estimation. These are the tools that connect the theoretical models of Parts III–IV to the data and to the estimation exercises of Part VII.

**Chapter 18** solves linear models with forward-looking expectations — the defining feature of New Keynesian macroeconomics. The method of undetermined coefficients (UMC) and McCallum's minimum-state-variable (MSV) procedure reduce infinite-horizon expectational difference equations to finite-dimensional algebraic problems. **Chapter 19** develops the time-series toolkit: unit roots, cointegration, VARs, structural identification, impulse responses, and forecast-error variance decompositions. **Chapter 20** derives the Kalman filter from scratch as the MMSE estimator for state-space models, applies it to the output gap and natural rate, and shows how it generates the likelihood function for DSGE estimation. **Chapter 21** completes the estimation toolkit with GMM and maximum likelihood, including full asymptotic theory, the delta method, and the application to the consumption Euler equation.

By the end of this part, the reader can take a log-linearized DSGE model (from Part VII), evaluate its likelihood using the Kalman filter, and estimate its parameters via GMM or MLE — the complete estimation pipeline used by every modern central bank model.

---

# Chapter 18: Rational Expectations

*Solving Linear Expectational Difference Equations*

> *"The assumption of rational expectations is not that agents are geniuses. It is that they are not systematically fooled by the same mistake, period after period."*
> — Robert Lucas

**Cross-reference:** *Principles* Ch. 15 (rational expectations definition and the Lucas critique); Ch. 16 (policy ineffectiveness, Barro–Gordon); Ch. 10 (NKPC forward solution as present-value formula) **[P:Ch.15, P:Ch.16, P:Ch.10]**

---

## 18.1 The Rational Expectations Problem

Chapter 4 introduced the issue of forward-looking difference equations — equations in which $\mathbb{E}_t[x_{t+1}]$ appears as an argument. Chapter 4's treatment was informal: we conjectured the MSV solution form and solved for coefficients by matching. This chapter develops the full formal theory.

A **rational expectations equilibrium** is a stochastic process $\{x_t\}$ such that: (i) households and firms form expectations $\mathbb{E}_t[x_{t+1}]$ as the mathematical conditional expectation given their information set; and (ii) the equilibrium conditions of the model are satisfied. The challenge is that $\mathbb{E}_t[x_{t+1}]$ depends on the equilibrium process for $\{x_t\}$, which in turn depends on expectations — a circular structure that requires a fixed-point argument to solve.

For **linear** models, this circular structure can be broken algebraically, yielding closed-form solutions. The tools developed here — undetermined coefficients, the MSV criterion, Blanchard–Kahn counting — are the foundation for Part VII's DSGE solution algorithms.

---

## 18.2 The Scalar Expectational Difference Equation

The simplest expectational difference equation:

$$x_t = a\mathbb{E}_t[x_{t+1}] + bz_t, \quad z_t = \rho z_{t-1} + \varepsilon_t, \quad \varepsilon_t \sim \text{WN}(0,\sigma^2).$$

Here $a, b$ are parameters, $z_t$ is an exogenous AR(1) forcing variable, and $x_t$ is the endogenous variable. We seek a solution $\{x_t\}$ consistent with rational expectations.

### 18.2.1 The Forward Solution

**Case 1: $|a| < 1$ (saddle-path stable).** Iterate the equation forward:

$$x_t = b\sum_{j=0}^\infty a^j\mathbb{E}_t[z_{t+j}] = b\sum_{j=0}^\infty a^j\rho^j z_t = \frac{b}{1-a\rho}z_t,$$

using $\mathbb{E}_t[z_{t+j}] = \rho^j z_t$ (the AR(1) forecast). This is the **unique bounded solution** when $|a| < 1$ and $|a\rho| < 1$.

**Case 2: $|a| > 1$ (indeterminate / forward-explosive).** With $|a| > 1$, the forward iteration diverges unless $x_t$ jumps to eliminate the explosive component. Rearrange as $\mathbb{E}_t[x_{t+1}] = x_t/a - bz_t/a$ and substitute the MSV guess $x_t = \omega z_t$:

$$\omega\rho z_t = \omega z_t/a - bz_t/a \implies \omega(1/a - \rho) = b/a \implies \omega = \frac{b}{1-a\rho}.$$

Remarkably, the same formula $\omega = b/(1-a\rho)$ applies regardless of whether $|a| < 1$ or $|a| > 1$! The distinction is in **uniqueness**: when $|a| < 1$ the solution is the unique bounded one; when $|a| > 1$ there are additional "sunspot" solutions.

### 18.2.2 Sunspot Solutions and Multiplicity

When $|a| > 1$, the general solution is:

$$x_t = \frac{b}{1-a\rho}z_t + \xi_t,$$

where $\xi_t$ is any process satisfying $\mathbb{E}_t[\xi_{t+1}] = \xi_t/a$. For example, $\xi_t = C\cdot a^{-t}$ for any constant $C$ is a solution (it dies out monotonically). But also any process with $\mathbb{E}_t[\xi_{t+1}] = \xi_t/a$ works — including **sunspot processes** driven by extraneous randomness.

**Definition 18.1 (Sunspot Equilibrium).** A **sunspot equilibrium** is a rational expectations equilibrium in which the endogenous variable responds to extraneous "sunspot" variables $\xi_t$ that have no fundamental economic content (they do not affect preferences, technology, or endowments). The sunspot variable satisfies $\xi_t = a^{-1}\xi_{t-1} + \eta_t$ where $\eta_t$ is an arbitrary martingale difference sequence.

When $|a| > 1$: multiple equilibria (MSV plus any sunspot). When $|a| < 1$: unique bounded equilibrium (MSV only). This is the scalar analogue of the Blanchard–Kahn counting rule.

---

## 18.3 The Method of Undetermined Coefficients (UMC)

The UMC provides a systematic procedure for finding the MSV solution of any linear rational expectations model.

**Algorithm 18.1 (Method of Undetermined Coefficients).**

1. **Identify the state variables** — the minimal set $\mathbf{s}_t$ such that the solution can be written $\mathbf{x}_t = \Omega\mathbf{s}_t$.
2. **Conjecture the MSV form** — $\mathbf{x}_t = \Omega\mathbf{s}_t$ where $\Omega$ is an unknown matrix.
3. **Compute $\mathbb{E}_t[\mathbf{x}_{t+1}]$** — using the law of motion of $\mathbf{s}_t$: $\mathbb{E}_t[\mathbf{x}_{t+1}] = \Omega\mathbb{E}_t[\mathbf{s}_{t+1}] = \Omega\Phi\mathbf{s}_t$.
4. **Substitute** into the model equations.
5. **Match coefficients** — equate coefficients of $\mathbf{s}_t$ on both sides to get a matrix equation for $\Omega$.
6. **Solve** for $\Omega$. Uniqueness of $\Omega$ is the MSV uniqueness condition.

The matrix equation obtained in step 5 is a **Sylvester equation** of the form:

$$\Omega = A\Omega\Phi + C,$$

where $A$, $\Phi$, and $C$ are known matrices from the model. Vectorizing: $\text{vec}(\Omega) = (I + \Phi'\otimes A)^{-1}\text{vec}(C)$. A unique solution $\Omega$ exists iff $(I + \Phi'\otimes A)$ is invertible — iff no eigenvalue of $\Phi$ equals the negative of any eigenvalue of $A$, which is generically satisfied.

---

## 18.4 The New Keynesian Model: Forward Solutions

*Cross-reference: Principles Ch. 10 (NKPC), Ch. 23 (NK policy)* **[P:Ch.10, P:Ch.23]**

### 18.4.1 The NKPC as a Forward Equation

The New Keynesian Phillips Curve:

$$\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t + u_t,$$

where $u_t$ is a cost-push shock. This is an expectational difference equation with $a = \beta$, which satisfies $|a| = \beta < 1$. Iterating forward:

$$\hat{\pi}_t = \kappa\sum_{j=0}^\infty\beta^j\mathbb{E}_t[\hat{x}_{t+j}] + \sum_{j=0}^\infty\beta^j\mathbb{E}_t[u_{t+j}].$$

**Definition 18.2 (Present-Value Inflation Formula).** The NKPC forward solution expresses current inflation as the present discounted value of all future expected output gaps and cost-push shocks:

$$\hat{\pi}_t = \kappa\sum_{j=0}^\infty\beta^j\mathbb{E}_t[\hat{x}_{t+j}] + \sum_{j=0}^\infty\beta^j\mathbb{E}_t[u_{t+j}].$$

This is exactly the formula in *Principles* Ch. 10.4 [P:Ch.10.4], now derived rigorously. Inflation today depends on the entire expected future path of the output gap — a forward-looking property that makes disinflation potentially costless if it is credible (agents immediately revise down all future inflation expectations).

### 18.4.2 The NK Three-Equation System Under a Taylor Rule

Combine the DIS, NKPC, and Taylor rule from *Principles* Ch. 23 [P:Ch.23]:

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(\hat{i}_t - \mathbb{E}_t[\hat{\pi}_{t+1}] - r^n_t) \quad \text{(DIS)}$$
$$\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t \quad \text{(NKPC)}$$
$$\hat{i}_t = \phi_\pi\hat{\pi}_t + \phi_y\hat{x}_t \quad \text{(Taylor rule)}$$

with exogenous shocks: $r^n_t = \rho_r r^n_{t-1} + \varepsilon^r_t$ (demand shock).

**Step 1: Substitute the Taylor rule into the DIS:**

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma\phi_\pi\hat{\pi}_t - \sigma\phi_y\hat{x}_t + \sigma\mathbb{E}_t[\hat{\pi}_{t+1}] + \sigma r^n_t.$$

Collecting $\hat{x}_t$:

$$(1+\sigma\phi_y)\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] + \sigma\mathbb{E}_t[\hat{\pi}_{t+1}] - \sigma\phi_\pi\hat{\pi}_t + \sigma r^n_t.$$

**Step 2: Stack the system.** Let $\mathbf{y}_t = (\hat{\pi}_t, \hat{x}_t)'$ and $z_t = r^n_t$. Write in the canonical form:

$$\underbrace{\begin{pmatrix}1 & -\kappa \\ \sigma\phi_\pi & 1+\sigma\phi_y\end{pmatrix}}_{\equiv\Gamma_0}\mathbf{y}_t = \underbrace{\begin{pmatrix}\beta & 0 \\ -\sigma & 1\end{pmatrix}}_{\equiv\Gamma_1}\mathbb{E}_t[\mathbf{y}_{t+1}] + \underbrace{\begin{pmatrix}0 \\ \sigma\end{pmatrix}}_{\equiv\Psi}z_t.$$

**Step 3: Reduce to standard form.** Pre-multiply by $\Gamma_0^{-1}$:

$$\mathbf{y}_t = \underbrace{\Gamma_0^{-1}\Gamma_1}_{\equiv A}\mathbb{E}_t[\mathbf{y}_{t+1}] + \underbrace{\Gamma_0^{-1}\Psi}_{\equiv C}z_t.$$

**Step 4: Apply UMC.** Guess $\mathbf{y}_t = \Omega z_t$. Then $\mathbb{E}_t[\mathbf{y}_{t+1}] = \Omega\rho_r z_t$. Substituting:

$$\Omega z_t = A\Omega\rho_r z_t + Cz_t \implies \Omega = A\Omega\rho_r + C.$$

This is the Sylvester equation for a $2\times1$ matrix $\Omega$ (since $z_t$ is scalar):

$$\Omega(I - A\rho_r) = C \implies \Omega = C(I - A\rho_r)^{-1}.$$

This requires $(I - A\rho_r)$ to be invertible — satisfied generically.

**Step 5: Verify determinacy.** The solution $\Omega$ is the unique MSV solution iff the NK system has a unique bounded equilibrium. For a system with two free (forward-looking) variables $(\hat\pi_t, \hat{x}_t)$, the Blanchard–Kahn condition requires **two eigenvalues of $A = \Gamma_0^{-1}\Gamma_1$ outside the unit circle**.

**Theorem 18.1 (Taylor Principle and Determinacy).** In the NK three-equation model, the Blanchard–Kahn condition is satisfied — the equilibrium is unique — if and only if:

$$\phi_\pi + \frac{1-\beta}{\kappa}\phi_y > 1.$$

For $\phi_y = 0$: the condition reduces to $\phi_\pi > 1$, the **Taylor principle**: the central bank must raise the nominal rate by more than one-for-one with inflation.

*Proof sketch.* The determinant and trace of $A = \Gamma_0^{-1}\Gamma_1$ can be computed explicitly. The eigenvalue analysis shows that both eigenvalues of $A$ have modulus greater than 1 iff the stated condition holds. See Chapter 28 for the full proof using the QZ decomposition. $\square$

---

## 18.5 Sunspots and Indeterminacy in the NK Model

When $\phi_\pi < 1$ (Taylor principle violated), the NK model has **indeterminate equilibria**: the MSV solution is one equilibrium, but there are infinitely many others driven by sunspot shocks. The general solution is:

$$\hat{\pi}_t = \omega_\pi z_t + \xi_t, \quad \hat{x}_t = \omega_x z_t + \varphi\xi_t,$$

where $\xi_t$ is an arbitrary sunspot process satisfying $\mathbb{E}_t[\xi_{t+1}] = \xi_t/\lambda_+$ (where $\lambda_+ < 1$ is the eigenvalue of $A$ inside the unit circle) and $\varphi$ is a coefficient determined by the eigenvector.

**Economic implication:** Sunspot equilibria generate volatility unrelated to fundamentals — the economy can fluctuate even when technology, preferences, and government policy are perfectly stable. The Taylor principle ($\phi_\pi > 1$) rules out sunspots by ensuring both eigenvalues of $A$ lie outside the unit circle, leaving the MSV solution as the unique equilibrium [P:Ch.23.1].

---

## 18.6 Worked Example: Impulse Response to a Cost-Push Shock

*Cross-reference: Principles Ch. 10.3, Ch. 23.2* **[P:Ch.10.3, P:Ch.23.2]**

**Calibration:** $\beta = 0.99$, $\kappa = 0.15$, $\sigma = 1$, $\phi_\pi = 1.5$, $\phi_y = 0.5$, $\rho_r = 0$ (i.i.d. shock, so $r^n_t = \varepsilon^r_t$).

Add a cost-push shock $u_t = \rho_u u_{t-1} + \varepsilon^u_t$ with $\rho_u = 0.5$.

**State:** $\mathbf{s}_t = (r^n_t, u_t)'$. Augment the NKPC: $\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t + u_t$.

**MSV guess:** $\mathbf{y}_t = \Omega\mathbf{s}_t$ where $\Omega$ is $2\times2$.

The Sylvester equation $\Omega = A\Omega\Phi + C$ where $\Phi = \text{diag}(\rho_r, \rho_u) = \text{diag}(0, 0.5)$ and $C$ now has a second column for $u_t$.

**Solving numerically:**

```apl
⍝ APL — NK model MSV solution
⎕IO←0 ⋄ ⎕ML←1

⍝ Parameters
beta←0.99  ⋄  kappa←0.15  ⋄  sigma←1  ⋄  phi_pi←1.5  ⋄  phi_y←0.5
rho_r←0.0  ⋄  rho_u←0.5

⍝ Γ₀ and Γ₁ matrices (stacked system: [NKPC; DIS with Taylor substituted])
G0 ← 2 2 ⍴ 1 (-kappa) (sigma×phi_pi) (1+sigma×phi_y)
G1 ← 2 2 ⍴ beta 0 (-sigma) 1
A  ← (⌹G0) +.× G1           ⍝ A = Γ₀⁻¹Γ₁

⍝ Eigenvalues of A (using ⎕PY for proper eigendecomposition)
⎕PY.Import 'numpy as np'
eigs ← A ⎕PY.Call 'np.linalg.eigvals'
|eigs    ⍝ should both be > 1 (Taylor principle satisfied)

⍝ C matrix: columns for [r^n shock, u shock]
⍝ From [NKPC: u enters as Ψ₂; DIS: r^n enters as Ψ₁]
G0_inv ← ⌹G0
Psi_r ← 0 sigma     ⍝ r^n enters DIS only (column vector)
Psi_u ← 1 0         ⍝ u enters NKPC only
C ← G0_inv +.× (2 2 ⍴ Psi_u , Psi_r)   ⍝ C is 2×2: columns are shocks

⍝ State transition
Phi ← 2 2 ⍴ rho_u 0 0 rho_r    ⍝ diagonal: [rho_u, rho_r]

⍝ Sylvester equation: Ω = A Ω Φ + C
⍝ Vectorised: vec(Ω) = (I - Φ'⊗A)⁻¹ vec(C)
⍝ For 2×2: vec gives 4-element vector; Φ'⊗A is 4×4
kron ← {⍺ ∘.× ⍵}                    ⍝ Kronecker product via outer
Phi_t ← ⍉Phi
KronPA ← Phi_t kron A                ⍝ 4×4 Kronecker product
I4 ← =⍨ ⍳ 4
vec_C ← , C                          ⍝ vectorise C (column-major)
vec_Omega ← vec_C ⌹ I4 - KronPA     ⍝ solve linear system
Omega ← 2 2 ⍴ vec_Omega             ⍝ reshape to 2×2

Omega    ⍝ [ω_π_u  ω_π_r; ω_x_u  ω_x_r]
⍝ Row 1: inflation responses to [u, r^n]
⍝ Row 2: output gap responses to [u, r^n]

⍝ Impulse responses: response to unit u shock (column 1 of Omega)
H ← 20
irf_u ← {Omega[;0]+.×Phi⍣⍵+.×(1 0)} ¨ ⍳H    ⍝ Omega × Phi^h × e_u
⍝ Stack into H×2 matrix
pi_irf ← {⊃⍵}¨ irf_u
x_irf  ← {⊃⌽⍵}¨ irf_u
pi_irf    ⍝ inflation impulse response
x_irf     ⍝ output gap impulse response
```

```python
import numpy as np

beta, kappa, sigma, phi_pi, phi_y = 0.99, 0.15, 1.0, 1.5, 0.5
rho_u, rho_r = 0.5, 0.0

# Structural matrices
G0 = np.array([[1, -kappa], [sigma*phi_pi, 1+sigma*phi_y]])
G1 = np.array([[beta, 0], [-sigma, 1]])
A  = np.linalg.inv(G0) @ G1

# Check Taylor principle
eigs = np.linalg.eigvals(A)
print(f"Eigenvalues of A: {np.abs(eigs)}  (both > 1: {np.all(np.abs(eigs)>1)})")

# Shock loadings (C matrix)
G0_inv = np.linalg.inv(G0)
Psi = np.array([[1, 0], [0, sigma]])  # [u in NKPC, r^n in DIS]
C = G0_inv @ Psi

# State transition matrix
Phi = np.diag([rho_u, rho_r])

# Sylvester equation: Omega = A @ Omega @ Phi + C
# Vectorized: (I - Phi' ⊗ A) vec(Omega) = vec(C)
from numpy import kron
KronPA = np.kron(Phi.T, A)
I4 = np.eye(4)
vec_Omega = np.linalg.solve(I4 - KronPA, C.flatten(order='F'))
Omega = vec_Omega.reshape(2, 2, order='F')
print(f"\nMSV solution matrix Ω:\n{np.round(Omega,3)}")
print(f"Interpretation: π̂_t = {Omega[0,0]:.3f}·u_t + {Omega[0,1]:.3f}·r^n_t")
print(f"                x̂_t = {Omega[1,0]:.3f}·u_t + {Omega[1,1]:.3f}·r^n_t")

# Impulse responses to cost-push shock u (column 0 of Omega × Phi^h × e_u)
H = 20
e_u = np.array([1, 0])  # unit shock to u
irf_pi = np.zeros(H); irf_x = np.zeros(H)
state = e_u.copy()
for h in range(H):
    response = Omega @ state
    irf_pi[h] = response[0]; irf_x[h] = response[1]
    state = Phi @ state

import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,4))
ax1.bar(range(H), irf_pi, color='steelblue')
ax1.axhline(0, color='k', lw=0.5); ax1.set_title('Inflation IRF (cost-push shock)')
ax1.set_xlabel('Quarters')
ax2.bar(range(H), irf_x, color='tomato')
ax2.axhline(0, color='k', lw=0.5); ax2.set_title('Output gap IRF (cost-push shock)')
ax2.set_xlabel('Quarters')
plt.tight_layout(); plt.show()
```

```julia
using LinearAlgebra

beta, kappa, sigma, phi_pi, phi_y = 0.99, 0.15, 1.0, 1.5, 0.5
rho_u, rho_r = 0.5, 0.0

G0 = [1 -kappa; sigma*phi_pi 1+sigma*phi_y]
G1 = [beta 0; -sigma 1]
A  = inv(G0) * G1

eigs = eigvals(A)
println("Eigenvalues of A: $(round.(abs.(eigs),digits=3)) (>1: $(all(abs.(eigs).>1)))")

G0_inv = inv(G0)
Psi = [1 0; 0 sigma]
C = G0_inv * Psi
Phi = diagm([rho_u, rho_r])

# Sylvester: Ω = A*Ω*Φ + C → (I - Φ'⊗A)vec(Ω) = vec(C)
KronPA = kron(Phi', A)
I4 = I(4)
vec_Omega = (I4 - KronPA) \ vec(C)
Omega = reshape(vec_Omega, 2, 2)

println("\nΩ (MSV solution):\n", round.(Omega, digits=3))

# IRF to unit cost-push shock
H = 20
irf = [Omega * (Phi^h * [1.0, 0.0]) for h in 0:H-1]
pi_irf = [r[1] for r in irf]
x_irf  = [r[2] for r in irf]
println("\nInflation IRF (first 5): $(round.(pi_irf[1:5], digits=3))")
println("Output IRF   (first 5): $(round.(x_irf[1:5],  digits=3))")
```

```r
beta<-0.99; kappa<-0.15; sigma<-1; phi_pi<-1.5; phi_y<-0.5
rho_u<-0.5; rho_r<-0.0

G0 <- matrix(c(1,sigma*phi_pi,-kappa,1+sigma*phi_y),2,2)
G1 <- matrix(c(beta,-sigma,0,1),2,2)
A  <- solve(G0) %*% G1

eigs <- eigen(A)$values
cat(sprintf("Eigenvalue moduli: %.3f, %.3f  (Taylor principle: %s)\n",
            Mod(eigs[1]), Mod(eigs[2]), all(Mod(eigs)>1)))

Psi <- matrix(c(1,0,0,sigma),2,2)
C   <- solve(G0) %*% Psi
Phi <- diag(c(rho_u, rho_r))

# Sylvester equation (using vec and Kronecker)
KronPA <- kronecker(t(Phi), A)
I4 <- diag(4)
vec_Omega <- solve(I4 - KronPA, as.vector(C))
Omega <- matrix(vec_Omega, 2, 2)
cat("\nΩ (MSV solution):\n"); print(round(Omega,3))

# IRF
H <- 20
state <- c(1, 0)
pi_irf <- x_irf <- numeric(H)
for(h in 1:H) {
  resp <- Omega %*% state; pi_irf[h] <- resp[1]; x_irf[h] <- resp[2]
  state <- Phi %*% state
}
cat("\nInflation IRF (first 5):", round(pi_irf[1:5],3))
```

---

## 18.7 Programming Exercises

### Exercise 18.1 (APL — Determinacy Region)

Compute the determinacy region in $(\phi_\pi, \phi_y)$ space for the NK model. For each pair $(\phi_\pi, \phi_y)$ on a $20\times20$ grid over $[0, 3]\times[0, 1.5]$: (a) compute $A = \Gamma_0^{-1}\Gamma_1$; (b) compute the eigenvalue moduli; (c) flag as determinate if both moduli exceed 1. Display as an APL Boolean matrix and verify the boundary matches $\phi_\pi + (1-\beta)\phi_y/\kappa = 1$.

### Exercise 18.2 (Python — Forward Guidance)

"Forward guidance" is a commitment to keep rates at zero for $T$ periods beyond when liftoff would normally occur. In the NK model, this means setting $\hat{i}_t = 0$ for $t = 0, \ldots, T-1$ and returning to the Taylor rule at $t = T$. (a) Solve the model recursively backward from $t = T$ (where the Taylor rule applies) to $t = 0$ (where the peg applies). (b) Plot the impulse response of inflation and the output gap at $t = 0$ as a function of $T$ for $T \in \{1, 2, 4, 8\}$. (c) Observe the "forward guidance puzzle": the effect grows exponentially in $T$ for standard calibrations.

### Exercise 18.3 (Julia — Sunspot Volatility)

```julia
# Compare fundamental vs. sunspot volatility
beta, kappa, sigma = 0.99, 0.15, 1.0

# Indeterminate calibration: phi_pi = 0.8 (below Taylor principle)
phi_pi_ind = 0.8; phi_y_ind = 0.0
G0_ind = [1 -kappa; sigma*phi_pi_ind 1+sigma*phi_y_ind]
G1_ind = [beta 0; -sigma 1]
A_ind  = inv(G0_ind) * G1_ind
eig_ind = eigvals(A_ind)
println("Indeterminate eigenvalues: $(round.(abs.(eig_ind),digits=3))")
# One eigenvalue inside unit circle → sunspots possible
lam_inside = eig_ind[argmin(abs.(eig_ind))]
println("Inside eigenvalue λ = $(round(lam_inside,digits=3))")
println("Sunspot process: ξ_t = $(round(1/abs(lam_inside),digits=3))·ξ_{t-1} + η_t")
```

### Exercise 18.4 — Hybrid NKPC ($\star$)

The hybrid NKPC (Galí and Gertler, 1999) has both forward- and backward-looking components:

$$\hat\pi_t = \omega_f\beta\mathbb{E}_t[\hat\pi_{t+1}] + \omega_b\hat\pi_{t-1} + \kappa\hat{x}_t, \quad \omega_f + \omega_b = 1.$$

Now $\hat\pi_{t-1}$ is a predetermined state variable. (a) Write the system in the form $\mathbf{y}_t = A\mathbb{E}_t[\mathbf{y}_{t+1}] + B\mathbf{y}_{t-1} + C\mathbf{z}_t$ with $\mathbf{y}_t = (\hat\pi_t, \hat{x}_t)'$. (b) The system now has one predetermined variable ($\hat\pi_{t-1}$) and two jump variables; the Blanchard–Kahn condition requires exactly one eigenvalue inside the unit circle. (c) Solve using the MSV approach and compare IRFs to the pure NKPC case.

---

## 18.8 Chapter Summary

**Key results:**

- The scalar expectational equation $x_t = a\mathbb{E}_t[x_{t+1}] + bz_t$ has MSV solution $x_t = [b/(1-a\rho)]z_t$, valid whenever $|a\rho| \neq 1$; unique bounded solution iff $|a| < 1$; multiple solutions (sunspots) iff $|a| > 1$.
- The **method of undetermined coefficients**: conjecture $\mathbf{y}_t = \Omega\mathbf{s}_t$, substitute, match coefficients to get the Sylvester equation $\Omega = A\Omega\Phi + C$, solve via $\text{vec}(\Omega) = (I-\Phi'\otimes A)^{-1}\text{vec}(C)$.
- The **NK model** in canonical form $\mathbf{y}_t = A\mathbb{E}_t[\mathbf{y}_{t+1}] + C\mathbf{z}_t$ is determinate iff both eigenvalues of $A = \Gamma_0^{-1}\Gamma_1$ lie outside the unit circle — the Blanchard–Kahn condition with $n_f = 2$ free variables.
- The **Taylor principle** ($\phi_\pi + (1-\beta)\phi_y/\kappa > 1$) is necessary and sufficient for determinacy in the NK model.
- In APL: the Kronecker product is `Phi_t kron A` via outer product `∘.×`; the Sylvester system is solved by `vec_C ⌹ I4 - KronPA`; IRFs are `{Omega +.× Phi⍣⍵ +.× e_shock}¨⍳H`.

*Next: Chapter 19 — Time Series Methods for Macroeconomic Data*
