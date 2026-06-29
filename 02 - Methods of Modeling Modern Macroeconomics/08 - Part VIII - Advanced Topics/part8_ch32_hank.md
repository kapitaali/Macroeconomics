# Part VIII: Advanced Topics in Macroeconomic Modeling

*Connects to: Principles Parts VII–VIII (HANK, climate, inequality, digital economy, future methods)*

---

The representative-agent DSGE of Parts IV–VII is powerful but has limits. It cannot address distributional questions: who gains and who loses from a rate cut? It cannot model systemic financial risk: the failure of one bank spreading to others through the network. It cannot account for the bounded rationality and coordination failures that characterize financial crises and technology adoptions. And it cannot represent the planetary-scale externalities of climate change.

This part surveys the frontier methods that take macroeconomics beyond the representative agent. **Chapter 32** develops heterogeneous-agent models — the Aiyagari (1994) model and the HANK framework — using the endogenous grid method for fast household problem solution, the Young (2010) non-stochastic simulation for the stationary distribution, and the Krusell–Smith algorithm for aggregate dynamics. **Chapter 33** brings the mathematical tools of inequality analysis: Pareto distributions, the Gini coefficient in closed form, and Auclert's (2019) URE framework for decomposing policy effects across the wealth distribution. **Chapter 34** develops network models for financial contagion, formalizing the Eisenberg–Noe clearing vector and the percolation threshold for cascade propagation. **Chapter 35** builds the DICE integrated assessment model from its Hamiltonian foundations, deriving the social cost of carbon as a costate variable and analyzing Stern–Nordhaus parameter sensitivity. **Chapter 36** introduces agent-based computational macroeconomics — heterogeneous, boundedly rational agents generating emergent business cycles — and situates it relative to DSGE.

Each chapter is self-contained and connects explicitly to the relevant chapters of *Principles* (HANK in Ch. 25; inequality in Ch. 38; systemic risk in Ch. 34; climate in Ch. 37; bounded rationality in Ch. 15 and Ch. 39).

---

# Chapter 32: Heterogeneous Agent Models

*The Aiyagari Model, EGM, and the Krusell–Smith Algorithm*

> *"The HANK model changes the monetary transmission mechanism: the response of aggregate consumption to a rate cut depends on the entire wealth distribution, not just the representative household's Euler equation."*

**Cross-reference:** *Principles* Ch. 25 (HANK model, Aiyagari–Bewley economy, heterogeneous MPC); Ch. 38 (wealth inequality, distributional consequences of policy) **[P:Ch.25, P:Ch.38]**

---

## 32.1 Beyond the Representative Agent

The representative-agent (RA) models of Parts III–VII implicitly assume all households are identical: same wealth, same income, same marginal propensity to consume. This is counterfactual. The U.S. wealth distribution has a Gini coefficient of approximately 0.87 — the top 1% hold about 38% of wealth; the bottom 50% hold less than 2%. Hand-to-mouth households (those with essentially zero liquid assets) have MPCs near 1; wealthy households have MPCs near 0. A fiscal stimulus that reaches the bottom 50% will have very different aggregate effects than one that reaches the top 10%.

**Definition 32.1 (Heterogeneous Agent New Keynesian — HANK — Model).** A **HANK model** is a general equilibrium model in which: (i) households differ in their wealth $a_i$ and earnings $y_i$; (ii) markets are incomplete (households cannot fully insure against idiosyncratic income risk); (iii) prices ($r$, $w$) are determined in equilibrium and affect all households; (iv) monetary policy operates through the distribution of wealth and income, not just the Euler equation of a representative household.

The Aiyagari (1994) model is the canonical static HANK — no aggregate uncertainty, but idiosyncratic income risk with a borrowing constraint. Krusell and Smith (1998) add aggregate (business cycle) shocks.

---

## 32.2 The Aiyagari Model: Setup

**Households:** A continuum of households $i \in [0,1]$ maximize lifetime utility subject to a borrowing constraint:

$$\max_{\{c_{i,t}\}}\mathbb{E}_0\sum_{t=0}^\infty\beta^t u(c_{i,t})$$

subject to:
$$a_{i,t+1} = (1+r)a_{i,t} + y_{i,t} - c_{i,t}, \quad a_{i,t+1} \geq \underline{a},$$

where $a_{i,t}$ is wealth (savings), $y_{i,t}$ is idiosyncratic labor income following a Markov chain $\Pi^y$, and $\underline{a} \geq 0$ is the borrowing limit.

**Firms:** A representative firm produces $Y = K^\alpha L^{1-\alpha}$ with $L = 1$ (inelastic labor supply), paying $r = \alpha K^{\alpha-1} - \delta$ and $w = (1-\alpha)K^\alpha$.

**Definition 32.2 (Aiyagari Stationary Equilibrium).** A **stationary competitive equilibrium** consists of: (1) value function $V(a, y)$ and policy function $c^*(a, y)$; (2) stationary wealth distribution $\mu^*(a, y)$ (a probability measure over $(a, y)$); (3) aggregate capital $K^* = \int a\,d\mu^*(a,y)$ and prices $(r^*, w^*)$; such that: (i) $c^*$ solves the household's problem given $(r^*, w^*)$; (ii) $\mu^*$ is invariant under $c^*$ and $\Pi^y$; (iii) markets clear: $K^* = \int a\,d\mu^*$.

The key challenge: finding $r^*$ and $K^*$ simultaneously, since household savings depend on $(r^*, w^*)$ but prices depend on $K^*$.

---

## 32.3 Solving the Household Problem: The Endogenous Grid Method

Standard VFI solves the Bellman equation by searching over consumption choices for each asset grid point. The **Endogenous Grid Method** (EGM; Carroll, 2006) inverts the Euler equation, dramatically accelerating the solution.

**Definition 32.3 (EGM Idea).** Instead of finding the optimal $c$ for a given $a$ (which requires a search), EGM:

1. Fixes a grid of **next-period asset values** $\{a'_j\}$.
2. Computes the implied consumption from the Euler equation (inverted).
3. Recovers the **current assets** $a = c + a'/(1+r) - y$ that are consistent with this choice.

The resulting $(a, c)$ pairs define the policy function on an **endogenous** current-asset grid.

**Algorithm 32.1 (Endogenous Grid Method).**

Initialize: $a'$-grid $\{a'_1 < \cdots < a'_N\}$, initial guess $c_0(a', y) = r\cdot a' + w\cdot y$.

For each iteration $n = 0, 1, 2, \ldots$:

For each income state $y_j$ and each $a'_k$:
1. Compute expected marginal utility: $\mathbb{E}_{y'|y_j}[u'(c^n(a'_k, y'))] = \sum_{l}\Pi^y_{jl}u'(c^n(a'_k, y'_l))$.
2. Invert Euler equation: $u'(c^*) = \beta(1+r)\mathbb{E}[u'(c^n(a'_k, y'))]$, giving $c^*(a'_k, y_j) = (u')^{-1}[\beta(1+r)\mathbb{E}[u'(c^n(a'_k, y'))]]$.
3. Recover endogenous grid: $a^* = c^*(a'_k, y_j) + a'_k/(1+r) - w\cdot y_j$.
4. Enforce constraint: if $a^* < \underline{a}$, set $c = (1+r)\underline{a} + w\cdot y_j - \underline{a}$ (constrained).

**For CRRA utility** $u'(c) = c^{-\sigma}$: $(u')^{-1}(x) = x^{-1/\sigma}$, so step 2 becomes:

$$c^*(a'_k, y_j) = \left[\beta(1+r)\sum_l\Pi^y_{jl}c^n(a'_k, y'_l)^{-\sigma}\right]^{-1/\sigma}.$$

No search required — this is a direct formula. The EGM reduces VFI's inner loop from an $O(N^2)$ search to $O(N)$ computation per grid point, providing a 100× speedup for $N = 100$.

---

## 32.4 The Stationary Distribution

Given the policy function $c^*(a, y)$ (equivalently, the savings function $g^*(a, y) = (1+r)a + wy - c^*(a,y)$), the stationary wealth distribution $\mu^*$ satisfies:

$$\mu^*(a', y') = \sum_y \Pi^y_{yy'}\int_{\{a: g^*(a,y) = a'\}}\mu^*(a, y)\,da.$$

This is a **fixed-point equation in the space of distributions** — computationally intensive to solve directly. Two practical approaches:

**Method 1 (Histogram simulation):** Simulate $M = 50{,}000$ households forward for $T = 1{,}000$ periods, starting from an arbitrary distribution. Average the last $T_{burn}$ cross-sections.

**Method 2 (Young, 2010 non-stochastic simulation):** Replace stochastic draws with deterministic mass flows. Maintain an $N_a \times N_y$ array of masses $\mu_{jk}$ on the grid. Update each period by mapping masses through the policy function using linear interpolation:

$$\mu'_{j,k} = \sum_{j'}\omega_{jj'}^{(k)} \mu_{j'k_{\text{src}}} \cdot \Pi^y_{k_{\text{src}},k},$$

where $\omega$ are interpolation weights. Young's method is deterministic, faster, and avoids sampling noise.

In APL:

```apl
⍝ APL — Stationary distribution via Young method
⎕IO←0 ⋄ ⎕ML←1

⍝ Given: policy g_star (N_a × N_y matrix: next-period assets)
⍝        transition Pi_y (N_y × N_y income transition matrix)
⍝        asset grid a_grid (N_a vector)

young_step ← {mu g_star Pi_y a_grid ← ⍵
    N_a ← ≢a_grid  ⋄  N_y ← ≢Pi_y
    mu_new ← (N_a, N_y) ⍴ 0
    
    ⍝ For each (j,k) current mass: map to next period via interpolation
    :For k :In ⍳N_y        ⍝ current income states
        :For j :In ⍳N_a    ⍝ current asset states
            a_next ← g_star[j;k]
            ⍝ Find interpolation index and weight
            idx ← 0 ⌈ (N_a-2) ⌊ ⌊ (a_next - a_grid[0]) × (N_a-1) ÷ a_grid[N_a-1]-a_grid[0]
            w   ← (a_next - a_grid[idx]) ÷ a_grid[idx+1] - a_grid[idx]
            ⍝ Distribute mass to income transitions
            :For l :In ⍳N_y    ⍝ next income states
                mu_new[idx;l]   +← (1-w) × mu[j;k] × Pi_y[k;l]
                mu_new[idx+1;l] +← w     × mu[j;k] × Pi_y[k;l]
            :EndFor
        :EndFor
    :EndFor
    mu_new}

⍝ Iterate until stationary: mu⍣≡ converges
converged ← {1e¯8 > ⌈/⌈/|⍺-⍵}
mu0 ← ⊃(N_a×N_y) ⍴ 1        ⍝ uniform initial distribution (normalized)
mu0 ÷← +/+/ mu0
mu_star ← young_step ⍣ converged ⊢ mu0
```

---

## 32.5 The Outer Equilibrium Loop

Given policy functions and the stationary distribution, the equilibrium interest rate $r^*$ is found by the condition $K^* = \int a\,d\mu^* = $ Firms' capital demand. This is a **single equation in one unknown** $r^*$:

$$\Psi(r) \equiv \int a\,d\mu^*(a, y; r) - \left(\frac{\alpha}{r+\delta}\right)^{1/(1-\alpha)} = 0,$$

where $\mu^*(a,y;r)$ is the stationary distribution for interest rate $r$. The function $\Psi(r)$ is monotone (higher $r$ reduces capital demand but increases household savings — typically $\Psi'(r) > 0$), so bisection or Newton–Raphson finds the root efficiently.

```apl
⍝ APL — Outer Aiyagari equilibrium loop
aiyagari_excess_demand ← {r_try ←⍵
    w_try ← (1-alpha) × (alpha÷r_try+delta)*alpha÷1-alpha
    ⍝ Solve household problem at (r_try, w_try) using EGM
    c_policy ← egm r_try w_try           ⍝ EGM returns N_a×N_y policy
    g_policy ← (1+r_try)×⊂a_grid ∘.+w_try×y_grid  ⍝ savings
    g_policy -← c_policy                  ⍝ next period assets
    ⍝ Compute stationary distribution
    mu ← young_step ⍣ converged ⊢ mu_uniform
    ⍝ Aggregate capital supply
    K_supply ← +/+/ a_grid ∘.× (+⌿ mu)   ⍝ E[a] under mu
    ⍝ Firms' capital demand
    K_demand ← (alpha÷r_try+delta)*÷1-alpha
    K_supply - K_demand}    ⍝ excess demand: = 0 at equilibrium

⍝ Find equilibrium r* via bisection on [r_lo, r_hi]
r_lo ← 0.001  ⋄  r_hi ← 0.04
bisection ← {a b ← ⍵
    m ← (a+b)÷2
    (aiyagari_excess_demand m) > 0: m b
    a m}
r_star ← {⍣ (0.0001∘>|⍺-⍵) ⊢ r_lo r_hi} bisection
r_star    ⍝ equilibrium interest rate
```

---

## 32.6 The Krusell–Smith Algorithm for Aggregate Shocks

When aggregate shocks (TFP shocks) are added to the Aiyagari model, the **state** of the economy includes the entire wealth distribution $\mu$ — an infinite-dimensional object. The **Krusell–Smith (1998)** algorithm approximates this by a small number of moments (typically just the mean $K = \int a\,d\mu$).

**Algorithm 32.2 (Krusell–Smith).**

1. **Parametric law of motion:** Guess that $\ln K_{t+1} = a_0 + a_1\ln K_t + a_2\ln A_t$ (log-linear in aggregate capital and TFP).
2. **Household problem:** Given the approximate law of motion, solve each household's consumption problem using EGM. The state now includes $(a_i, y_i, K, A)$.
3. **Simulate:** Simulate $M$ households for $T$ periods, drawing aggregate shocks from the TFP process.
4. **Update coefficients:** Regress realized $\ln K_{t+1}$ on $(1, \ln K_t, \ln A_t)$ to update $(a_0, a_1, a_2)$.
5. **Iterate:** Repeat steps 2–4 until the regression $R^2 > 0.9999$ (Krusell–Smith accuracy criterion).

**Why it works:** The key insight of Krusell and Smith is that **aggregate capital** is an approximately sufficient statistic for the wealth distribution, for the purposes of forecasting next period's prices. The regression $R^2 > 0.9999$ verifies this numerically.

---

## 32.7 Reiter's Method: Linearizing Around the Stationary Distribution

**Reiter (2009)** provides an alternative to the simulation-based Krusell–Smith approach. Instead of simulating the model, Reiter linearizes the full system (including the distribution) around the stationary distribution and solves the resulting large-dimensional DSGE system using the methods of Chapters 27–28.

The state vector includes the discretized distribution $\bm\mu \in \mathbb{R}^{N_a N_y}$ plus aggregate variables. The linearized system has dimension $N_a N_y + n_{agg}$ — typically $N_a N_y \approx 200\times3 = 600$ plus a handful of aggregate variables. The gensys algorithm handles systems of this size in seconds.

**Reiter's advantage over Krusell–Smith:** Full second moments (not just means) of the distribution are tracked; the method is suitable for models where the KS approximate aggregation fails. **Reiter's cost:** The system size is large; solving the $600\times600$ eigenvalue problem in Chapter 28 takes longer than the $4\times4$ standard DSGE.

---

## 32.8 Continuous-Time HANK: HJB and KFP Equations

The **Achdou et al. (2022)** continuous-time HANK model replaces discrete Bellman equations with differential equations, enabling finite-element methods (HACT, Huggett, 2011 framework).

**Definition 32.4 (Hamilton–Jacobi–Bellman Equation).** In continuous time with income process $\{y_t\}$ (Poisson jumps with intensity $\lambda$), the value function $V(a, y)$ satisfies:

$$\rho V(a,y) = \max_c\left\{u(c) + V_a(a,y)\cdot\dot{a} + \lambda\sum_{y'}\Pi^y_{yy'}[V(a,y')-V(a,y)]\right\},$$

where $\dot{a} = ra + wy - c$ (continuous-time budget constraint).

**Definition 32.5 (Kolmogorov–Fokker–Planck Equation).** The stationary density $\mu(a, y)$ satisfies:

$$0 = -\frac{\partial}{\partial a}[\dot{a}^*(a,y)\mu(a,y)] + \lambda\sum_{y'}\Pi^y_{y'y}\mu(a,y') - \lambda\mu(a,y),$$

where $\dot{a}^*(a,y) = ra + wy - c^*(a,y)$ is the optimal asset drift.

The HJB and KFP equations are coupled: HJB gives the policy function $c^*$; KFP gives the distribution $\mu$ consistent with $c^*$. Achdou et al. (2022) solve them simultaneously using a finite-difference scheme, achieving 100× speedups over discrete-time methods.

---

## 32.9 Worked Example: Aiyagari Model Calibration and Wealth Gini

*Cross-reference: Principles Ch. 25 (Aiyagari–Bewley), Ch. 38 (wealth distribution)* **[P:Ch.25, P:Ch.38]**

```python
import numpy as np
from scipy.interpolate import interp1d

# Aiyagari model: minimal implementation
alpha, delta, beta, sigma_crra = 0.36, 0.025, 0.96, 2.0
# Income process: 2-state Markov (employed/unemployed)
y_grid  = np.array([0.1, 1.0])    # income in each state
Pi_y    = np.array([[0.1, 0.9],
                    [0.05, 0.95]]) # transition: rows=current, cols=next

# Asset grid
a_min, a_max, N_a = 0.0, 40.0, 200
a_grid = np.linspace(a_min, a_max, N_a)

def egm_step(c_old, r, w, a_grid, y_grid, Pi_y, beta, sigma):
    """One EGM iteration."""
    N_a, N_y = len(a_grid), len(y_grid)
    c_new = np.zeros((N_a, N_y))
    for j_y, y_curr in enumerate(y_grid):
        # Expected marginal utility at each a' grid point
        EMU = np.zeros(N_a)
        for l_y, y_next in enumerate(y_grid):
            # Interpolate c_old(a', y_next) on a_grid
            c_interp = interp1d(a_grid, c_old[:, l_y], fill_value='extrapolate')
            c_next = np.maximum(c_interp(a_grid), 1e-10)
            EMU += Pi_y[j_y, l_y] * c_next**(-sigma)
        # Invert Euler equation: c* = [beta*(1+r)*EMU]^{-1/sigma}
        c_endo = (beta * (1+r) * EMU)**(-1/sigma)
        # Endogenous current assets: a = c + a'/(1+r) - w*y
        a_endo = c_endo + a_grid / (1+r) - w * y_curr
        # Interpolate back to exogenous grid (clamp to borrowing limit)
        above_bc = a_endo >= a_min
        if above_bc.sum() < 2:
            c_new[:, j_y] = (1+r)*a_grid + w*y_curr - a_min  # all constrained
        else:
            a_valid = a_endo[above_bc]; c_valid = c_endo[above_bc]
            c_interp2 = interp1d(a_valid, c_valid, bounds_error=False,
                                 fill_value=(c_valid[0], c_valid[-1]))
            c_new[:, j_y] = np.maximum(c_interp2(a_grid), 1e-10)
            # Constrained region
            constrained = a_grid < a_valid[0]
            c_new[constrained, j_y] = (1+r)*a_grid[constrained] + w*y_curr - a_min
    return c_new

def find_equilibrium(alpha, delta, beta, sigma_crra, y_grid, Pi_y, a_grid,
                     r_lo=0.001, r_hi=0.038, tol=1e-4, max_iter=50):
    """Find Aiyagari equilibrium r* via bisection."""
    N_a = len(a_grid); N_y = len(y_grid)
    
    def excess_supply(r):
        w = (1-alpha) * (alpha/(r+delta))**(alpha/(1-alpha))
        # Solve household problem via EGM
        c = np.outer(a_grid, np.ones(N_y)) * 0.05 + 0.5  # initial guess
        for _ in range(500):
            c_new = egm_step(c, r, w, a_grid, y_grid, Pi_y, beta, sigma_crra)
            if np.max(np.abs(c_new - c)) < 1e-7: break
            c = c_new
        # Savings function
        g = np.zeros((N_a, N_y))
        for j_y, y_j in enumerate(y_grid):
            g[:, j_y] = (1+r)*a_grid + w*y_j - c[:, j_y]
        # Stationary distribution via simulation
        np.random.seed(42); M = 5000; T = 2000; burnin = 500
        a_sim = np.ones(M) * 5.0
        y_idx = np.zeros(M, dtype=int)
        mu_cdf = np.cumsum(Pi_y, axis=1)
        dist = []
        for t in range(T):
            a_sim_new = np.zeros(M)
            for i in range(M):
                j = y_idx[i]
                g_interp = interp1d(a_grid, g[:, j], fill_value='extrapolate')
                a_sim_new[i] = max(a_min, float(g_interp(a_sim[i])))
                u = np.random.rand()
                y_idx[i] = np.searchsorted(mu_cdf[j], u)
            a_sim = a_sim_new
            if t >= burnin: dist.append(a_sim.copy())
        a_dist = np.concatenate(dist)
        K_supply = np.mean(a_dist)
        K_demand = (alpha/(r+delta))**(1/(1-alpha))
        return K_supply - K_demand, a_dist
    
    for _ in range(max_iter):
        r_mid = (r_lo + r_hi) / 2
        excess, a_dist = excess_supply(r_mid)
        if abs(excess) < tol: break
        if excess > 0: r_hi = r_mid
        else: r_lo = r_mid
    
    return r_mid, a_dist

r_star, a_dist = find_equilibrium(alpha, delta, beta, sigma_crra, y_grid, Pi_y, a_grid)
K_star = (alpha/(r_star+delta))**(1/(1-alpha))
w_star = (1-alpha) * K_star**alpha
print(f"Equilibrium: r*={r_star*100:.2f}%, K*={K_star:.3f}, w*={w_star:.3f}")

# Wealth Gini
a_sorted = np.sort(a_dist)
n = len(a_sorted)
lorenz = np.cumsum(a_sorted) / np.sum(a_sorted)
gini = 1 - 2*np.trapz(lorenz, np.linspace(0, 1, n))
print(f"Wealth Gini: {gini:.3f}  (data: ~0.87 for U.S.)")
print(f"Top 10% share: {np.sum(a_sorted[int(0.9*n):])/np.sum(a_sorted)*100:.1f}%")
print(f"Bottom 50% share: {np.sum(a_sorted[:n//2])/np.sum(a_sorted)*100:.1f}%")
```

---

## 32.10 Programming Exercises

### Exercise 32.1 (APL — EGM Implementation)

Implement the EGM in APL as a dfn `egm ← {r w ← ⍵ ⋄ ...}`. The core step: `c_endo ← (beta×(1+r)×EE)*÷-sigma` where `EE ← Pi_y +.× (c_old*-sigma)` (matrix multiply of transition with marginal utilities). Verify the policy function converges and matches the Python implementation.

### Exercise 32.2 (Julia — Stationary Distribution via Matrix Iteration)

In Julia, implement Young's (2010) non-stochastic simulation as a matrix iteration: build the $N_aN_y \times N_aN_y$ transition matrix $\mathbf{Q}$ from the policy function and the income transition $\Pi^y$, then find the stationary distribution as the left eigenvector: `mu = eigvecs(Q')[:,1]` (normalized). Compare to the histogram simulation for the same model.

### Exercise 32.3 (Python — Krusell–Smith Regression)

```python
# Validate the KS approximate aggregation hypothesis
# Test whether K_{t+1} can be accurately predicted from (K_t, A_t) alone
# Simulate 10,000 periods of the Aiyagari model with aggregate TFP shocks
# Add TFP: Y_t = A_t * K_t^alpha, A_t = exp(z_t), z_t = rho*z_{t-1} + eps
import numpy as np; np.random.seed(42)
T = 5000; rho_A = 0.95; sig_A = 0.007
z = np.zeros(T)
for t in range(1,T): z[t] = rho_A*z[t-1] + sig_A*np.random.randn()
A = np.exp(z)

# Simulate K path (placeholder: K near steady state with AR1 dynamics + TFP loading)
K_ss = 5.0; lam = 0.9  # persistence
K = np.zeros(T); K[0] = K_ss
for t in range(1,T):
    K[t] = K_ss*(1-lam) + lam*K[t-1] + 0.3*(A[t]-1)*K_ss + 0.01*np.random.randn()

# KS regression
X = np.column_stack([np.ones(T-1), np.log(K[:-1]), np.log(A[:-1])])
y_reg = np.log(K[1:])
b = np.linalg.lstsq(X, y_reg, rcond=None)[0]
y_hat = X @ b
R2 = 1 - np.var(y_reg - y_hat)/np.var(y_reg)
print(f"KS regression: a0={b[0]:.4f}, a1={b[1]:.4f}, a2={b[2]:.4f}")
print(f"R² = {R2:.6f}  (KS threshold: R² > 0.9999)")
```

### Exercise 32.4 — HJB Finite Differences ($\star$)

For the simple one-income-state Aiyagari model (no idiosyncratic risk), the HJB equation reduces to $\rho V(a) = u(c^*(a)) + V'(a)(ra + w - c^*(a))$. (a) Discretize $V'(a)$ using upwind finite differences: $V'(a_j) \approx [V(a_{j+1})-V(a_j)]/\Delta a$ if $\dot{a}^* > 0$, else $[V(a_j)-V(a_{j-1})]/\Delta a$. (b) Write the resulting linear system $(\rho I - B)V = u$ where $B$ is the finite-difference matrix. (c) Solve with `scipy.sparse.linalg.spsolve` and compare to the discrete-time solution.

---

## 32.11 Chapter Summary

**Key results:**

- The **Aiyagari (1994) equilibrium** is a triple $(c^*, \mu^*, r^*)$ satisfying household optimization, distribution invariance, and market clearing; finding $r^*$ requires an outer bisection loop.
- The **Endogenous Grid Method** inverts the Euler equation to give $c^*(a'_k, y_j) = [\beta(1+r)\mathbb{E}[u'(c^n(a'_k, y'))]^{-1/\sigma}$ — avoiding search entirely; $O(N)$ per iteration vs. $O(N^2)$ for VFI.
- **Young's (2010) non-stochastic simulation** propagates a mass distribution through the policy function using linear interpolation — deterministic, no sampling noise, faster than histogram simulation.
- The **Krusell–Smith algorithm** approximates the law of motion for aggregate capital by a log-linear regression on $(K_t, A_t)$, validated by $R^2 > 0.9999$.
- **Reiter's method** linearizes the full HANK system (including distribution dynamics) around the stationary distribution, reducing to a large-dimensional standard DSGE solvable by gensys.
- The **HJB and KFP equations** (Achdou et al., 2022) are the continuous-time counterparts of the Bellman equation and the stationarity condition; finite-difference methods solve them simultaneously.
- In APL: EGM is `c_endo ← (beta×(1+r)×Pi_y+.×c_old*(-sigma))*÷(-sigma)`; stationary distribution via `young_step ⍣ converged ⊢ mu0`; equilibrium via `aiyagari ⍣ (tol∘>|excess) ⊢ r0`.

*Next: Chapter 33 — Income Inequality and Macroeconomic Dynamics*
