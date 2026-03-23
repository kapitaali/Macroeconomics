# Chapter 40: Policy Analysis with a New Keynesian Model

*Commitment, Discretion, ELB, and Welfare*

> *"The Taylor principle is not a recommendation — it is a necessary condition for the price level to be determined."*

**Cross-reference:** *Principles* Ch. 23 (optimal monetary policy, commitment vs. discretion, Taylor rule); Ch. 28 (fiscal policy, zero lower bound); Ch. 29 (forward guidance, QE) **[P:Ch.23, P:Ch.28, P:Ch.29]**

---

## 40.1 The Policy Analysis Pipeline

Chapter 31 showed how to solve the NK model. This chapter uses that solution for policy analysis — the ultimate purpose of DSGE modeling at central banks. The standard pipeline:

1. **Solve the model** under a benchmark policy rule (Taylor rule with standard coefficients).
2. **Compute welfare** under the benchmark.
3. **Solve the optimal policy problem** (commitment or discretion).
4. **Compute the welfare gain** from optimal relative to benchmark.
5. **Sensitivity analysis**: how does the welfare comparison depend on the model parameters?

We develop each step formally, using the NK three-equation model as the vehicle.

---

## 40.2 Welfare Measurement in the NK Model

*Cross-reference: Chapter 29 (second-order approximation required for welfare)* **[M:Ch.29]**

The central bank's loss function is the second-order approximation to household welfare. For CRRA utility and Calvo pricing, Woodford (2003) shows:

**Theorem 40.1 (NK Welfare Loss Function).** To second order, the deviation of household welfare from the first-best (flexible-price) allocation is:

$$\mathcal{W} = -\frac{1}{2}\mathbb{E}\sum_{t=0}^\infty\beta^t\left[\hat\pi_t^2 + \frac{\kappa}{\varepsilon}\hat{x}_t^2\right] + O(\varepsilon^3),$$

where $\varepsilon$ is the order of the approximation and $\kappa/\varepsilon = \lambda_x$ is the weight on the output gap (proportional to the slope of the NKPC $\kappa$ divided by the demand elasticity $\varepsilon_{NK}$).

*Proof sketch.* The second-order expansion of household utility around the zero-inflation flexible-price steady state yields cross-terms involving $\hat\pi$ and $\hat x$. The inflation term arises from price dispersion across Calvo firms (firms with different prices produce different quantities, creating output loss). The output gap term arises from the wedge between the flexible-price output gap and zero. The cross-term vanishes at the optimum (Woodford, 2003, Ch. 6). $\square$

**The welfare loss is:**

$$\mathcal{L} = \text{Var}(\hat\pi_t) + \lambda_x\text{Var}(\hat{x}_t),$$

where $\lambda_x = \kappa/(varepsilon_{NK}) \approx 0.05$ for standard calibrations. Inflation variance is heavily weighted relative to output gap variance — consistent with the "divine coincidence" [P:Ch.23.3]: stabilizing inflation and stabilizing the output gap are usually aligned.

---

## 40.3 Optimal Policy Under Commitment

Under **commitment**, the central bank pre-announces and commits to a complete sequence of future actions $\{i_t\}_{t=0}^\infty$ at time 0. This commitment is credible — agents believe the sequence will be followed.

**The commitment problem:**

$$\min_{\{i_t, \hat\pi_t, \hat{x}_t\}}\mathbb{E}_0\sum_{t=0}^\infty\beta^t\left[\hat\pi_t^2 + \lambda_x\hat{x}_t^2\right]$$

subject to:
$$\hat\pi_t = \beta\mathbb{E}_t[\hat\pi_{t+1}] + \kappa\hat{x}_t + u_t \quad \forall t \quad \text{(NKPC)}$$
$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(\hat{i}_t - \mathbb{E}_t[\hat\pi_{t+1}] - r^n_t) \quad \forall t \quad \text{(DIS)}$$

**Theorem 40.2 (Optimal Commitment Policy — Targeting Rule).** The optimal commitment policy satisfies the **targeting rule**:

$$\hat\pi_t = -\frac{\lambda_x}{\kappa}(\hat{x}_t - \hat{x}_{t-1}).$$

*Proof.* Form the Lagrangian with multipliers $\varphi_t$ for the NKPC constraint:

$$\mathcal{L}^{Lagr} = \mathbb{E}\sum_t\beta^t\left[\hat\pi_t^2 + \lambda_x\hat{x}_t^2 + \varphi_t(\hat\pi_t - \beta\hat\pi_{t+1} - \kappa\hat{x}_t - u_t)\right].$$

The first-order conditions:
- w.r.t. $\hat\pi_t$: $2\hat\pi_t + \varphi_t - \beta^{-1}\varphi_{t-1} = 0 \Rightarrow \varphi_t = \beta^{-1}\varphi_{t-1} - 2\hat\pi_t$.
- w.r.t. $\hat{x}_t$: $2\lambda_x\hat{x}_t - \kappa\varphi_t = 0 \Rightarrow \varphi_t = 2\lambda_x\hat{x}_t/\kappa$.

Substituting: $2\lambda_x\hat{x}_t/\kappa = \beta^{-1}\cdot 2\lambda_x\hat{x}_{t-1}/\kappa - 2\hat\pi_t$, giving $\hat\pi_t = -(\lambda_x/\kappa)(\hat{x}_t - \beta^{-1}\hat{x}_{t-1})$. For $\beta\approx1$: $\hat\pi_t \approx -(\lambda_x/\kappa)(\hat{x}_t - \hat{x}_{t-1})$. $\square$

**Interpretation:** Under commitment, the central bank allows some initial inflation (when hit by a cost-push shock) but then actively deflates to restore the price level — **price level targeting** (PL targeting). This contrasts with **inflation targeting** (IT), which targets $\hat\pi_t = 0$ at all times.

The targeting rule $\hat\pi_t = -(\lambda_x/\kappa)(\hat{x}_t - \hat{x}_{t-1})$ implements **price-level targeting**: because $\hat\pi_t = \hat{p}_t - \hat{p}_{t-1}$, the rule implies $\hat{p}_t = \hat{p}_{t-1} - (\lambda_x/\kappa)(\hat{x}_t - \hat{x}_{t-1})$, keeping the price level near a target path rather than just the inflation rate.

---

## 40.4 Optimal Policy Under Discretion

Under **discretion**, the central bank re-optimizes each period, taking private sector expectations as given (but knowing they are formed under rational expectations). This is a game between the central bank and the private sector.

**The equilibrium under discretion** satisfies:

$$\hat\pi_t = -\frac{\lambda_x}{\kappa}\hat{x}_t, \quad \text{(discretion targeting rule)}$$

combined with the NKPC $\hat\pi_t = \beta\mathbb{E}_t[\hat\pi_{t+1}] + \kappa\hat{x}_t + u_t$.

The discretion rule is **simpler**: $\hat\pi_t/\hat{x}_t = -\lambda_x/\kappa$ — a constant ratio of inflation to output gap, independent of history (unlike the commitment rule which depends on the lagged output gap $\hat{x}_{t-1}$).

**The commitment-discretion gap:** Under a cost-push shock $u_t = \rho_u u_{t-1} + \varepsilon_t$:

$$\mathcal{L}^{disc} > \mathcal{L}^{comm}, \quad \text{with gap: } \Delta\mathcal{L} = \frac{\lambda_x\rho_u^2}{(\lambda_x+\kappa^2)(1-\beta\rho_u)}\sigma_u^2 > 0.$$

The commitment policy achieves lower welfare loss by exploiting forward-looking inflation expectations: announcing to keep output low in the future reduces current inflation expectations, lowering the sacrifice ratio.

---

## 40.5 The ELB and Forward Guidance

**Definition 40.1 (Effective Lower Bound).** The **effective lower bound** (ELB) binds when the natural rate $r^n_t < 0$ and the central bank cannot lower the policy rate below zero. At the ELB: $\hat{i}_t = 0$ (rather than following the Taylor rule).

**The ELB problem:** At the ELB, the DIS equation gives:

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] + \sigma\mathbb{E}_t[\hat\pi_{t+1}] + \sigma r^n_t.$$

With $r^n_t < 0$ and no ability to cut rates, the economy falls into a recession ($\hat{x}_t < 0$) and deflation ($\hat\pi_t < 0$) — the **liquidity trap**.

**Forward guidance** — committing to keep rates low even after the ELB ceases to bind — provides stimulus at the ELB by raising $\mathbb{E}_t[\hat\pi_{t+1}]$ (expectation of future inflation lowers the real rate today).

**Theorem 40.3 (Forward Guidance Multiplier).** In the NK model at the ELB for $T$ periods, a commitment to keep rates at zero for $\tau$ periods beyond the natural liftoff date generates a GDP multiplier at the ELB:

$$\text{FG multiplier}(\tau) = \frac{\sigma\kappa}{\kappa^2+\lambda_x}\cdot\frac{1-(beta\rho)^{\tau+1}}{1-\beta\rho},$$

which grows with $\tau$ — potentially exponentially for $\beta\rho \to 1$ (the **forward guidance puzzle**).

*Derivation.* In the NK model with ELB for $T$ periods and rate kept at zero for $T+\tau$ periods, the output gap at the ELB satisfies a recursion backward from date $T+\tau$. Each additional period of forward guidance adds the term $(\beta\rho)^\tau\sigma\kappa/(\kappa^2+\lambda_x)$ to the impact on current output. Summing gives the stated formula. $\square$

The forward guidance puzzle: standard NK calibrations imply unrealistically large effects of forward guidance far in the future. Resolution: habit formation, finite planning horizons, or agent inattention (Gabaix, 2020) attenuate the effect of distant promises.

---

## 40.6 Worked Example: Welfare Cost of Suboptimal Policy

```python
import numpy as np
from scipy.linalg import solve_discrete_lyapunov

# NK model: welfare analysis under different policy rules
beta, kappa, sigma = 0.99, 0.15, 1.0
lambda_x = kappa / 6.0  # Woodford (2003): ε = 6 → λ_x = κ/ε ≈ 0.025
rho_u, sigma_u = 0.5, 0.01

def welfare_loss(phi_pi, phi_y, rule='taylor'):
    """Compute welfare loss under a given policy rule."""
    G0 = np.array([[1, -kappa], [sigma*phi_pi, 1+sigma*phi_y]])
    G1 = np.array([[beta, 0], [-sigma, 1]])
    
    try:
        A = np.linalg.inv(G0) @ G1
    except: return np.inf
    
    eigs = np.linalg.eigvals(A)
    if not np.all(np.abs(eigs) > 1): return np.inf  # indeterminate
    
    # MSV solution for cost-push shock
    C_u = np.linalg.inv(G0) @ np.array([1.0, 0.0])
    Omega = np.linalg.solve(np.eye(2) - rho_u*A, C_u)
    
    # Theoretical variance of [pi, x] from cost-push shock
    var_u = sigma_u**2 / (1 - rho_u**2)
    var_pi = Omega[0]**2 * var_u
    var_x  = Omega[1]**2 * var_u
    
    return var_pi + lambda_x * var_x

# Compare policy rules
rules = {
    'Taylor (1.5, 0.5)': (1.5, 0.5),
    'Aggressive (3.0, 0.5)': (3.0, 0.5),
    'Output-focused (1.5, 2.0)': (1.5, 2.0),
    'IT (φ_y=0)': (1.5, 0.0),
}

print("Welfare loss under alternative Taylor rules:")
baseline = welfare_loss(1.5, 0.5)
for name, (phi_pi, phi_y) in rules.items():
    L = welfare_loss(phi_pi, phi_y)
    print(f"  {name}: L = {L*1e4:.4f}×10⁻⁴  (ratio to baseline: {L/baseline:.3f})")

# Compute optimal Taylor rule
from scipy.optimize import minimize
result = minimize(lambda p: welfare_loss(p[0], p[1]),
                  [1.5, 0.5], method='Nelder-Mead',
                  bounds=[(1.01, 10), (0, 5)])
phi_opt = result.x
L_opt = result.fun
print(f"\nOptimal Taylor rule: φ_π={phi_opt[0]:.3f}, φ_y={phi_opt[1]:.3f}")
print(f"Welfare gain vs. baseline: {100*(baseline-L_opt)/baseline:.2f}%")

# Commitment vs. discretion
# Discretion: pi_t = -(lambda_x/kappa)*x_t
# From NKPC: pi*(1+lambda_x/kappa^2)*(1-rho_u*beta) = (rho_u*sigma_u/...)
# Simplified: welfare under discretion
var_pi_disc = (rho_u**2*sigma_u**2/(1-rho_u**2)) * (kappa/(kappa**2+lambda_x))**2
var_x_disc  = (rho_u**2*sigma_u**2/(1-rho_u**2)) * (lambda_x/kappa/(kappa**2+lambda_x))**2 * kappa**2
L_disc = var_pi_disc + lambda_x * var_x_disc

# Optimal commitment: lower by the commitment-discretion gap
delta_L = lambda_x * rho_u**2 / ((lambda_x+kappa**2)*(1-beta*rho_u)) * sigma_u**2
L_comm = L_disc - delta_L

print(f"\nCommitment vs. Discretion:")
print(f"  Discretion: L = {L_disc*1e4:.4f}×10⁻⁴")
print(f"  Commitment: L = {L_comm*1e4:.4f}×10⁻⁴")
print(f"  Gap: {(L_disc-L_comm)*1e4:.4f}×10⁻⁴ ({100*(L_disc-L_comm)/L_disc:.1f}%)")
```

```julia
using Optim, LinearAlgebra

beta, kappa, sigma = 0.99, 0.15, 1.0
lambda_x = kappa/6.0; rho_u, sigma_u = 0.5, 0.01

function welfare_nk(phi_pi, phi_y)
    G0 = [1 -kappa; sigma*phi_pi 1+sigma*phi_y]
    G1 = [beta 0; -sigma 1]; A = inv(G0)*G1
    all(abs.(eigvals(A)).>1) || return Inf
    Omega = (I(2)-rho_u*A)\(inv(G0)*[1;0])
    var_u = sigma_u^2/(1-rho_u^2)
    Omega[1]^2*var_u + lambda_x*Omega[2]^2*var_u
end

res = optimize(p->welfare_nk(p[1],p[2]), [1.5,0.5], NelderMead())
phi_opt = Optim.minimizer(res)
println("Optimal φ_π=$(round(phi_opt[1],digits=3)), φ_y=$(round(phi_opt[2],digits=3))")
println("Welfare gain: $(round(100*(welfare_nk(1.5,0.5)-welfare_nk(phi_opt[1],phi_opt[2]))/welfare_nk(1.5,0.5),digits=2))%")
```

---

## 40.7 Programming Exercises

### Exercise 40.1 (APL — Optimal Commitment Policy)

The commitment targeting rule $\hat\pi_t = -(\lambda_x/\kappa)(\hat{x}_t - \hat{x}_{t-1})$ can be implemented as an additional equation in the NK system. (a) Augment the $\Gamma_0, \Gamma_1$ matrices to include the commitment targeting rule and $\hat{x}_{t-1}$ as a predetermined state variable. (b) Verify the Blanchard–Kahn conditions: now $\hat{x}_{t-1}$ is predetermined and $\hat\pi_t, \hat{x}_t$ are jump variables — 2 jump variables need 2 unstable eigenvalues. (c) Compare the IRF to a cost-push shock under commitment vs. the Taylor rule.

### Exercise 40.2 (Python — ELB Dynamics)

Model the ELB regime as a Markov-switching model: in normal times, the Taylor rule applies; at the ELB, $\hat{i}_t = \min(\text{Taylor}, 0) = 0$. (a) Solve the model recursively backward from the date when the ELB stops binding. (b) Compute the sequence of $(\hat\pi_t, \hat{x}_t)$ during a 4-quarter ELB episode driven by a demand shock ($r^n_t < 0$). (c) Add forward guidance: keep rates at zero for 2 additional quarters after the natural rate returns to positive. Plot the difference in the IRF.

### Exercise 40.3 (Julia — Ramsey Optimal Policy)

```julia
# Ramsey-optimal fiscal policy: optimal tax-spending mix
# CB minimizes pi^2 + lambda_x * x^2 + lambda_i * i^2 (instrument cost)
function ramsey_optimal(lambda_x, lambda_i; rho_u=0.5, sigma_u=0.01)
    # Grid search over (phi_pi, phi_y, phi_i) -- include interest rate smoothing
    best_L = Inf; best_p = (1.5, 0.5, 0.0)
    for phi_pi in 1.1:0.2:4.0
        for phi_y in 0:0.2:2.0
            # Compute model-implied variances (simplified)
            G0 = [1 -0.15; sigma*phi_pi 1+sigma*phi_y]
            G1 = [0.99 0; -sigma 1]; A = inv(G0)*G1
            all(abs.(eigvals(A)).>1) || continue
            Omega = (I(2)-rho_u*A)\(inv(G0)*[1;0])
            var_u = sigma_u^2/(1-rho_u^2)
            L = Omega[1]^2*var_u + lambda_x*Omega[2]^2*var_u
            L += lambda_i*(phi_pi*Omega[1]+phi_y*Omega[2])^2*var_u  # instrument cost
            L < best_L && (best_L = L; best_p = (phi_pi, phi_y, 0.0))
        end
    end
    best_p, best_L
end

params, L = ramsey_optimal(0.025, 0.02)
println("Ramsey optimal: φ_π=$(params[1]), φ_y=$(params[2]), L=$(round(L*1e4,digits=4))×10⁻⁴")
```

### Exercise 40.4 — Price-Level Targeting vs. Inflation Targeting ($\star$)

Compare PLT (commitment) to IT (discretion or Taylor rule with $\phi_\pi > 1$): (a) compute the IRF to a cost-push shock under both regimes; (b) show that under PLT, inflation rises on impact but then deflates to below target — consistent with the targeting rule; (c) compute the welfare loss under both and verify PLT delivers lower loss; (d) implement price-level targeting in Dynare by adding $\hat{p}_t$ as a variable and replacing the Taylor rule with $\hat{i}_t = \phi_p\hat{p}_t + \phi_y\hat{x}_t$.

---

## 40.8 Chapter Summary

**Key results:**

- The **NK welfare loss function** $\mathcal{L} = \text{Var}(\hat\pi) + \lambda_x\text{Var}(\hat{x})$ is the second-order welfare approximation; $\lambda_x = \kappa/\varepsilon_{NK} \approx 0.025$ (Theorem 40.1, proved from the second-order approximation).
- The **optimal commitment targeting rule** $\hat\pi_t = -(\lambda_x/\kappa)(\hat{x}_t - \hat{x}_{t-1})$ implies price-level targeting — the central bank partially reverses past inflation deviations (Theorem 40.2, proved from the Lagrangian FOCs).
- The **commitment-discretion gap** $\Delta\mathcal{L} = \lambda_x\rho_u^2/[({\lambda_x+\kappa^2})(1-\beta\rho_u)]\sigma_u^2 > 0$: commitment achieves lower welfare loss by exploiting forward expectations.
- The **forward guidance multiplier** (Theorem 40.3) grows with commitment duration $\tau$; the forward guidance puzzle arises when $\beta\rho_u \to 1$.
- **Optimal Taylor rule** computed by minimizing the Lyapunov welfare loss over $(\phi_\pi, \phi_y)$; gains of 10–30% vs. standard (1.5, 0.5) coefficients for aggressive anti-inflation stance.
- In APL: welfare is `var_pi + lambda_x × var_x` where variances come from `Sigma_y ← Omega +.× Sigma_z +.× ⍉Omega`; optimal policy via `⍣≡` Newton on the welfare gradient.

*Next: Chapter 41 — Model Validation and Sensitivity Analysis*
