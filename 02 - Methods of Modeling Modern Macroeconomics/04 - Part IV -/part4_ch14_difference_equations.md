# Part IV: Discrete-Time Dynamic Models in Macroeconomics

*Connects to: Principles Ch. 11 (consumption), Ch. 25 (OLG), Ch. 27 (RBC)*

---

Most policy models, DSGE models, and estimation exercises operate in discrete time. The reasons are practical: macroeconomic data arrive at quarterly or monthly intervals; calibration targets are sample moments computed from those data; and the numerical algorithms used for estimation — particle filters, MCMC samplers, value function iteration on finite grids — are inherently discrete.

This part develops the discrete-time parallel to Part III. Where Part III deployed differential equations, phase planes, and Pontryagin's maximum principle, this part uses difference equations, dynamic programming, and Bellman's principle of optimality. Every result has a direct continuous-time analogue; spotting the correspondence — $e^{at} \leftrightarrow \lambda^t$, stability condition $a < 0$ $\leftrightarrow$ $|\lambda| < 1$, the Hamiltonian $\leftrightarrow$ the Bellman equation — is the intellectual thread connecting the two parts.

**Chapter 14** develops the theory of difference equations systematically, solving the cobweb model and the Samuelson multiplier–accelerator in closed form, and deriving the full stability diagram in parameter space. **Chapter 15** builds dynamic programming from the ground up — the Bellman equation, the contraction mapping theorem, value function iteration (VFI), and the Howard policy improvement algorithm — applying them to the buffer-stock consumption model. **Chapter 16** solves the Diamond (1965) OLG model in discrete time, deriving dynamic efficiency conditions and analyzing PAYG social security. **Chapter 17** assembles the full stochastic RBC model as a recursive competitive equilibrium and solves it via VFI over the $(k, A)$ state space, computing the second moments used for calibration evaluation.

These four chapters are the computational foundation of Part VII. The reader who understands VFI, the Blanchard–Kahn counting rule, and the recursive competitive equilibrium definition has the tools to understand every step of the DSGE pipeline.

---

# Chapter 14: Difference Equations in Macroeconomics

*The Cobweb Model and Inventory Cycles*

> *"In discrete time, every model is an iterated map. The question is whether that map converges, oscillates, or explodes — and the eigenvalues of the linearized system give the answer."*

**Cross-reference:** *Principles* Ch. 6 (time-series properties of macro data); Ch. 17 (goods market, inventory adjustment model) **[P:Ch.6, P:Ch.17]**

---

## 14.1 Motivating Examples: Why Difference Equations Arise in Macroeconomics

Three canonical macroeconomic models generate difference equations naturally.

**The inventory adjustment model** [P:Ch.17.1]: firms target an inventory stock $N^* = \nu Y^d$ proportional to expected demand. If actual inventories $N_{t-1}$ differ from the target, firms partially adjust: $\Delta Y_t^s = \beta(N^* - N_{t-1})$, $\beta \in (0,1)$. With output determining income and demand: $Y_t = \mu \bar{A} - \beta\nu(1-b)Y_{t-1}$, a first-order linear difference equation.

**The cobweb model**: agricultural supply decisions take one period to implement. Farmers plant based on last period's price; consumers buy at the current clearing price. The resulting price dynamics: $P_{t+1} = -(S'/D')P_t + \text{const}$, where $S' > 0$ and $D' < 0$ are supply and demand slopes.

**The multiplier–accelerator** [P:Ch.8]: Samuelson's (1939) model combines the Keynesian income multiplier with an investment accelerator, generating a second-order difference equation that can produce oscillatory business cycles.

These three models span the key structural features of macroeconomic difference equations: first-order monotone dynamics, first-order oscillatory dynamics (cobweb with $|S'/D'| > 1$), and second-order complex dynamics (multiplier-accelerator with complex eigenvalues). The theory developed below characterizes all three.

---

## 14.2 First-Order Linear Difference Equations: Complete Solution

**Definition 14.1 (First-Order Linear Difference Equation).** A **first-order linear difference equation** has the form:

$$x_{t+1} = ax_t + b, \quad t = 0, 1, 2, \ldots$$

with given initial condition $x_0$.

### 14.2.1 Homogeneous Solution

For $b = 0$: $x_{t+1} = ax_t$ has solution $x_t = a^t x_0$.

**Stability:** $x_t \to 0$ iff $|a| < 1$; $x_t$ diverges iff $|a| > 1$; $x_t$ cycles between $\pm x_0$ iff $a = -1$; $x_t = x_0$ for all $t$ iff $a = 1$.

### 14.2.2 Particular Solution and General Solution

For $a \neq 1$: the steady state is $x^* = b/(1-a)$. The general solution:

$$\boxed{x_t = a^t(x_0 - x^*) + x^*, \quad x^* = \frac{b}{1-a}.}$$

The term $a^t(x_0 - x^*)$ is the **transient** (decays iff $|a| < 1$); $x^*$ is the **steady state** (particular solution).

**Convergence:**

| Condition | Behavior |
|---|---|
| $0 < a < 1$ | Monotone convergence to $x^*$ |
| $-1 < a < 0$ | Oscillatory convergence (alternating signs) to $x^*$ |
| $a > 1$ | Monotone divergence |
| $a < -1$ | Oscillatory divergence |
| $a = 1$ | No finite steady state; $x_t = x_0 + bt$ (drift) |
| $a = -1$ | Oscillation between $x_0$ and $b - x_0$ |

### 14.2.3 The Cobweb Model: First-Order Dynamics

In the cobweb model, supply is based on last period's price: $Q^s_t = \alpha + \beta P_{t-1}$, $\beta > 0$. Demand: $Q^d_t = \gamma - \delta P_t$, $\delta > 0$. Market clearing $Q^s_t = Q^d_t$:

$$P_t = \frac{\gamma - \alpha}{\delta} - \frac{\beta}{\delta}P_{t-1}.$$

This is a first-order difference equation with $a = -\beta/\delta$ and $b = (\gamma-\alpha)/\delta$. Stability requires $|a| < 1$, i.e., $\beta < \delta$: supply must be more inelastic than demand.

**Rational expectations correction:** If farmers correctly anticipate the equilibrium price $P_t^e = P^* = b/(1+\beta/\delta)$ for all $t$, there are no cobweb dynamics — markets clear immediately at the steady state. The cobweb arises only under naive expectations ($P^e_t = P_{t-1}$). This motivates the rational expectations assumption: if agents understand the model, they eliminate self-fulfilling oscillations [P:Ch.16].

---

## 14.3 Second-Order Linear Difference Equations

**Definition 14.2 (Second-Order Linear Difference Equation).** A **second-order linear difference equation** has the form:

$$x_{t+2} + px_{t+1} + qx_t = c,$$

with initial conditions $x_0$ and $x_1$. The **characteristic equation** is:

$$\lambda^2 + p\lambda + q = 0, \quad \lambda_{1,2} = \frac{-p \pm \sqrt{p^2 - 4q}}{2}.$$

**Theorem 14.1 (General Solution — Second-Order).** Let $\lambda_1, \lambda_2$ be the roots of the characteristic equation and $x^* = c/(1+p+q)$ (assuming $1+p+q \neq 0$).

**Case 1: Distinct real roots ($p^2 > 4q$).**
$$x_t = C_1\lambda_1^t + C_2\lambda_2^t + x^*.$$

**Case 2: Repeated root ($p^2 = 4q$, $\lambda_1 = \lambda_2 = \lambda$).**
$$x_t = (C_1 + C_2 t)\lambda^t + x^*.$$

**Case 3: Complex conjugate roots ($p^2 < 4q$), $\lambda_{1,2} = r e^{\pm i\theta}$.**
$$x_t = r^t(C_1\cos\theta t + C_2\sin\theta t) + x^*,$$
where $r = \sqrt{q}$ and $\theta = \arctan\left(\sqrt{4q-p^2}/(-p)\right) \in (0,\pi)$.

**Proof of Case 3.** Write $\lambda = r(\cos\theta + i\sin\theta)$ and $\bar\lambda = r(\cos\theta - i\sin\theta)$. Then $\lambda^t = r^t e^{i\theta t} = r^t(\cos\theta t + i\sin\theta t)$, and the real general solution combines $\lambda^t$ and $\bar\lambda^t$ with conjugate coefficients, giving the stated form. $\square$

**Stability for second-order equations:** The system is stable (both roots with $|\lambda| < 1$) iff:

$$|p| < 1 + q \quad \text{and} \quad q < 1.$$

These are necessary and sufficient. For complex roots ($r = \sqrt{q}$), stability reduces to $r < 1$, i.e., $q < 1$.

---

## 14.4 The Multiplier–Accelerator: Second-Order Dynamics in Full

**Definition 14.3 (Samuelson Multiplier–Accelerator).** Samuelson's (1939) model:

$$C_t = bY_{t-1}, \quad I_t = v(C_t - C_{t-1}), \quad Y_t = C_t + I_t + G.$$

Substituting: $Y_t = b(1+v)Y_{t-1} - bvY_{t-2} + G$, or:

$$Y_{t+2} - b(1+v)Y_{t+1} + bv\cdot Y_t = G.$$

This is a second-order difference equation with $p = -b(1+v)$ and $q = bv$. The steady state $Y^* = G/(1-b)$.

**Characteristic equation:** $\lambda^2 - b(1+v)\lambda + bv = 0$.

**Discriminant:** $\Delta = b^2(1+v)^2 - 4bv = b^2(1-v)^2 + 2bv(b-2)$... Let me compute directly: $\Delta = [b(1+v)]^2 - 4bv = b^2 + 2bv + b^2v^2 - 4bv = b^2(1-v)^2 + 2b(bv-2v) = b^2(1+v)^2 - 4bv$.

**Case analysis for business cycle character:**

| Region | Eigenvalue type | Dynamic path | 
|---|---|---|
| $b^2(1+v)^2 > 4bv$ | Two distinct real | Monotone approach |
| $b^2(1+v)^2 < 4bv$ | Complex conjugates | Oscillatory (cyclical) |
| $bv < 1$ | $r = \sqrt{bv} < 1$ | Convergent (dampened) cycles |
| $bv = 1$ | $r = 1$ | Neutral oscillation (limit cycle) |
| $bv > 1$ | $r > 1$ | Explosive oscillation |

**Theorem 14.2 (Stability Conditions for Multiplier–Accelerator).** The Samuelson multiplier-accelerator is stable iff $bv < 1$, and generates business cycles (oscillatory path) iff $b^2(1+v)^2 < 4bv$.

*Proof.* Stability requires $q = bv < 1$ (necessary and sufficient for the complex-root modulus to be less than 1). Oscillation requires $p^2 < 4q$, i.e., $b^2(1+v)^2 < 4bv$. $\square$

**The period of oscillation:** For complex roots with modulus $r$ and argument $\theta$, the period of one oscillation is $T = 2\pi/\theta$ periods. With $b = 0.8$ and $v = 0.9$: $r = \sqrt{0.72} = 0.849$, $\theta = \arccos(-b(1+v)/(2r)) = \arccos(-0.72/(0.849\times2\times1.9)) ≈ \arccos(-0.224) ≈ 1.796$ radians, period $\approx 2\pi/1.796 \approx 3.5$ periods — a business cycle of roughly 3–4 years for quarterly data.

---

## 14.5 Systems of Linear Difference Equations

**Definition 14.4 (State-Space Form).** A system of $n$ first-order linear difference equations:

$$\mathbf{x}_{t+1} = A\mathbf{x}_t + B\mathbf{u}_t, \quad \mathbf{x}_0 \text{ given},$$

where $\mathbf{x}_t \in \mathbb{R}^n$ is the **state vector**, $\mathbf{u}_t \in \mathbb{R}^m$ is an exogenous input, $A \in \mathbb{R}^{n\times n}$ is the **transition matrix**, and $B \in \mathbb{R}^{n\times m}$ is the **input matrix**.

**Theorem 14.3 (Solution of Linear State-Space System).** The solution of $\mathbf{x}_{t+1} = A\mathbf{x}_t$ is:

$$\mathbf{x}_t = A^t\mathbf{x}_0 = PD^tP^{-1}\mathbf{x}_0,$$

where $A = PDP^{-1}$ is the eigendecomposition. The general solution with inputs is:

$$\mathbf{x}_t = A^t\mathbf{x}_0 + \sum_{j=0}^{t-1}A^{t-1-j}B\mathbf{u}_j.$$

**Stability:** The system is stable iff all eigenvalues of $A$ satisfy $|\lambda_i| < 1$.

**Impulse response:** The response at horizon $h$ to a unit shock $\mathbf{u}_0 = \mathbf{e}_j$ is $A^{h-1}B\mathbf{e}_j = A^{h-1}B_j$ (the $j$-th column of $A^{h-1}B$). The full IRF sequence $\{B_j, AB_j, A^2B_j, \ldots\}$ decays to zero iff all $|\lambda_i(A)| < 1$.

In APL, the IRF at all horizons 0..H-1 is generated in one expression:

```apl
⍝ APL — Impulse response function from state-space system
⎕IO←0 ⋄ ⎕ML←1

A ← 2 2 ⍴ 0.9 0.1 0.0 0.8    ⍝ transition matrix
B ← 2 1 ⍴ 0.5 0.3             ⍝ shock loading vector
H ← 20                         ⍝ horizons

⍝ IRF: {A^h × B} for h=0,1,...,H-1
irf ← {(A⍣⍵) +.× B} ¨ ⍳ H   ⍝ apply A^h each, using ⍣h

⍝ Stack into H×2 matrix
⊃ irf
```

---

## 14.6 Worked Example: Full Solution of the Multiplier–Accelerator

**Setup:** $b = 0.8$, $v = 0.8$, $G = 100$.

**Steady state:** $Y^* = G/(1-b) = 100/0.2 = 500$.

**Characteristic equation:** $\lambda^2 - 0.8(1.8)\lambda + 0.8\times0.8 = \lambda^2 - 1.44\lambda + 0.64 = 0$.

**Roots:** $\lambda = (1.44 \pm \sqrt{1.44^2 - 4\times0.64})/2 = (1.44 \pm \sqrt{2.0736 - 2.56})/2$.

Discriminant: $2.0736 - 2.56 = -0.4864 < 0$ → complex roots.

$r = \sqrt{q} = \sqrt{0.64} = 0.8$, $\theta = \arccos(-p/(2r)) = \arccos(1.44/(2\times0.8)) = \arccos(0.9) \approx 0.4510$ radians.

Period: $2\pi/0.4510 \approx 13.9$ periods.

**General solution:**
$$Y_t = 500 + 0.8^t(C_1\cos(0.451t) + C_2\sin(0.451t)).$$

**Initial conditions:** $Y_0 = 480$ (4% below steady state), $Y_1 = 490$ (2% below).

From $Y_0 = 500 + C_1 = 480$: $C_1 = -20$.

From $Y_1 = 500 + 0.8(C_1\cos0.451 + C_2\sin0.451) = 490$:
$0.8(-20\times0.900 + C_2\times0.436) = -10 \Rightarrow -14.4 + 0.349C_2 = -10 \Rightarrow C_2 = 4.4/0.349 = 12.6$.

**Verified:** The economy oscillates with period ≈ 14 periods (3.5 years quarterly) and dampens at rate $r^t = 0.8^t$ — Kydland–Prescott business cycles are calibrated to match such dynamics.

```apl
⍝ APL — Multiplier-accelerator simulation and stability diagram
⎕IO←0 ⋄ ⎕ML←1

b ← 0.8  ⋄  v ← 0.8  ⋄  G ← 100
Ystar ← G ÷ 1-b    ⍝ steady state = 500

⍝ One-step recursion: Y_t = b(1+v)*Y_{t-1} - bv*Y_{t-2} + G
⍝ State = [Y_{t-1}, Y_t]; step advances to [Y_t, Y_{t+1}]
step ← {
    Yprev Ycurr ← ⍵
    Ynext ← b×(1+v)×Ycurr - b×v×Yprev + G
    Ycurr Ynext}

⍝ Simulate 80 periods from initial conditions
Y0 ← 480  ⋄  Y1 ← 490
path_pairs ← step \ 80 ⍴ ⊂ Y0 Y1   ⍝ scan collects all pairs
Y_path ← {⊃⌽ ⍵}¨ path_pairs          ⍝ extract second element (current Y)

⍝ Verify convergence to steady state
Y_path[0] Y_path[39] Y_path[79]       ⍝ should trend toward 500

⍝ Stability diagram: 20×20 grid over (b, v) ∈ [0,1]×[0,3]
b_grid ← (⍳20) ÷ 20
v_grid ← (⍳20) × 3 ÷ 20
⍝ Stability: bv < 1
stable_mat ← b_grid ∘.{⍺×⍵ < 1} v_grid     ⍝ 20×20 Boolean stability matrix
⍝ Complex roots: b²(1+v)² < 4bv
complex_mat ← b_grid ∘.{(⍺*2)×(1+⍵)*2) < 4×⍺×⍵} v_grid
⍝ Oscillatory and stable = complex_mat AND stable_mat
osc_stable ← stable_mat ∧ complex_mat
⍝ Each region characterised: 
⍝ stable_mat ∧ ~complex_mat  → monotone convergence
⍝ osc_stable                  → dampened oscillation (business cycles)
⍝ ~stable_mat ∧ complex_mat  → explosive oscillation
⍝ ~stable_mat ∧ ~complex_mat → monotone divergence
```

```python
import numpy as np; import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

b_vals = np.linspace(0.01, 0.99, 200)
v_vals = np.linspace(0.01, 2.99, 200)
B, V = np.meshgrid(b_vals, v_vals)

stable   = B*V < 1
complex_ = B**2*(1+V)**2 < 4*B*V

region = np.zeros_like(B, dtype=int)
region[stable & ~complex_]  = 1  # monotone convergence
region[stable & complex_]   = 2  # oscillatory convergence (cycles)
region[~stable & complex_]  = 3  # explosive oscillation
region[~stable & ~complex_] = 4  # monotone divergence

cmap = ListedColormap(['white','lightblue','steelblue','tomato','salmon'])
plt.figure(figsize=(8,6))
plt.contourf(B, V, region, levels=[-0.5,0.5,1.5,2.5,3.5,4.5], cmap=cmap)
plt.contour(B, V, stable.astype(int), levels=[0.5], colors='black', linewidths=2)
plt.xlabel('MPC (b)'); plt.ylabel('Accelerator (v)')
plt.title('Stability diagram: multiplier-accelerator\n(dark blue = oscillatory convergence = business cycles)')
plt.colorbar(ticks=[0,1,2,3,4], label='Region')
plt.tight_layout(); plt.show()

# Time path simulation
b, v, G = 0.8, 0.8, 100
Y = np.zeros(80); Y[0], Y[1] = 480, 490
for t in range(2, 80):
    Y[t] = b*(1+v)*Y[t-1] - b*v*Y[t-2] + G
plt.figure(); plt.plot(Y); plt.axhline(G/(1-b), c='r', ls='--', label='Y*=500')
plt.xlabel('Period'); plt.ylabel('Output'); plt.title('Multiplier-Accelerator Path (b=0.8, v=0.8)'); plt.legend(); plt.show()
```

```julia
b, v, G = 0.8, 0.8, 100.0
Ystar = G / (1-b)

# Characteristic roots
p_coef, q_coef = -b*(1+v), b*v
disc = p_coef^2 - 4*q_coef
if disc < 0
    r = sqrt(q_coef); theta = acos(-p_coef/(2r))
    println("Complex roots: r=$(round(r,digits=3)), θ=$(round(theta,digits=3)), period=$(round(2π/theta,digits=1))")
else
    lam = [(-p_coef + sqrt(disc))/2, (-p_coef - sqrt(disc))/2]
    println("Real roots: λ₁=$(round(lam[1],digits=3)), λ₂=$(round(lam[2],digits=3))")
end

Y = zeros(80); Y[1], Y[2] = 480.0, 490.0
for t in 3:80
    Y[t] = b*(1+v)*Y[t-1] - b*v*Y[t-2] + G
end
println("Y[1]=$(Y[1]), Y[40]=$(round(Y[40],digits=1)), Y[80]=$(round(Y[80],digits=1)), Y*=$(Ystar)")
```

```r
b <- 0.8; v <- 0.8; G <- 100
Ystar <- G/(1-b)
Y <- numeric(80); Y[1] <- 480; Y[2] <- 490
for(t in 3:80) Y[t] <- b*(1+v)*Y[t-1] - b*v*Y[t-2] + G
cat(sprintf("Steady state: %.0f; Y[40]=%.1f; Y[80]=%.1f\n", Ystar, Y[40], Y[80]))
# Period of oscillation
disc <- (b*(1+v))^2 - 4*b*v
if(disc < 0) {
  r <- sqrt(b*v); theta <- acos(b*(1+v)/(2*r))
  cat(sprintf("Period ≈ %.1f quarters\n", 2*pi/theta))
}
plot(Y, type='l', xlab='Period', ylab='Output', main='Multiplier-Accelerator')
abline(h=Ystar, lty=2, col='red')
```

---

## 14.7 Programming Exercises

### Exercise 14.1 (APL — Cobweb Stability)

Write a dfn `cobweb ← {alpha beta gamma delta T ← ⍵ ⋄ ...}` that simulates $T$ periods of the cobweb model starting from an arbitrary price $P_0$. (a) Verify convergence when $|\beta/\delta| < 1$ and divergence when $|\beta/\delta| > 1$. (b) Use `∘.{|⍺÷⍵} < 1` to generate a $10\times10$ Boolean stability grid over $(\beta, \delta)$ pairs and display it.

### Exercise 14.2 — Analytical Period Calculation ($\star$)

For the multiplier-accelerator with $b = 0.75$ and $v \in \{0.5, 1.0, 1.5, 2.0, 2.5\}$: (a) compute the oscillation period $2\pi/\theta$ analytically for each $v$; (b) determine the range of $v$ for which oscillations occur; (c) show that as $v \to \infty$ the period approaches $2\pi/\arccos(-b(1+v)/(2\sqrt{bv})) \to 2$ — period-2 oscillations.

### Exercise 14.3 — Companion Form Eigenvalues ($\star$)

Convert the second-order equation $Y_{t+2} - b(1+v)Y_{t+1} + bvY_t = G$ to a state-space system with state $\mathbf{x}_t = (Y_t, Y_{t-1})'$ and companion matrix $A = \begin{pmatrix} b(1+v) & -bv \\ 1 & 0 \end{pmatrix}$. (a) Verify that the eigenvalues of $A$ equal the roots of the characteristic equation. (b) Show that stability of the companion system ($|\lambda_i(A)| < 1$) is equivalent to the stability condition $bv < 1$.

### Exercise 14.4 — Inventory Cycle ($\star\star$)

Derive the inventory adjustment model from first principles [P:Ch.17.1]. Let firms target $N^* = \nu Y$, partially adjust inventories: $N_t = N_{t-1} + \beta(N^* - N_{t-1})$, and produce $Y_t^s = Y^d_{t-1} + \beta(N^* - N_{t-1})$. With $Y^d = \bar{A} + bY$: (a) derive the difference equation for $Y_t$; (b) find the stability condition; (c) calibrate with $b = 0.8$, $\beta = 0.3$, $\nu = 0.5$ and simulate the impulse response to a 10-unit increase in autonomous demand $\bar{A}$.

---

## 14.8 Chapter Summary

**Key results:**

- First-order equations $x_{t+1} = ax_t + b$ have steady state $x^* = b/(1-a)$ and general solution $x_t = a^t(x_0-x^*)+x^*$; stable iff $|a| < 1$.
- The **cobweb model** is a first-order equation with $a = -\beta/\delta$; stable iff supply is more inelastic than demand ($\beta < \delta$); rational expectations eliminate cobweb dynamics entirely.
- **Second-order equations** $x_{t+2}+px_{t+1}+qx_t = c$ have solutions characterized by their characteristic roots — real monotone, complex oscillatory, or explosive — with stability requiring $|p| < 1+q$ and $q < 1$.
- The **multiplier–accelerator** $Y_t = b(1+v)Y_{t-1} - bvY_{t-2} + G$ is stable iff $bv < 1$ and oscillatory (business cycles) iff $b^2(1+v)^2 < 4bv$; the oscillation period is $2\pi/\arccos[-b(1+v)/(2\sqrt{bv})]$.
- State-space form $\mathbf{x}_{t+1} = A\mathbf{x}_t$: stable iff all eigenvalues of $A$ satisfy $|\lambda_i| < 1$; IRFs are $A^h$ applied to the shock vector.
- In APL: simulation is a scan `step \ T ⍴ ⊂ initial_state`; IRFs are `{(A⍣⍵)+.×B}¨⍳H`; the stability diagram is `b_grid ∘.{⍺×⍵ < 1} v_grid`.

*Next: Chapter 15 — Dynamic Programming for Consumption and Saving*
