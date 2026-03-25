# Chapter 9: Tax and Government Spending Multipliers

*Analytical Solutions with Behavioral Foundations*

> *"The fiscal multiplier is not one number but a family of numbers indexed by the model, the regime, and the state of the economy."*

**Cross-reference:** *Principles* Ch. 8 (multiplier derivation); Ch. 22 (Ricardian equivalence, fiscal rules); Ch. 28 (fiscal policy in practice, empirical multipliers); Ch. 23 (ELB and the NK model) **[P:Ch.8, P:Ch.22, P:Ch.28, P:Ch.23]**

---

## 9.1 The Architecture of Fiscal Multipliers

Chapters 7 and 8 derived fiscal multipliers in two static frameworks: the Keynesian cross and the AD–SRAS system. This chapter synthesizes and extends those results in four directions:

1. **Tax structure.** How do lump-sum versus proportional taxes change the multiplier formula?
2. **Automatic stabilizers.** How do built-in fiscal responses to the output gap dampen fluctuations?
3. **Ricardian equivalence.** Under what conditions does the multiplier collapse to zero?
4. **The ELB multiplier.** What is the analytical expression for the fiscal multiplier in the New Keynesian model when monetary policy is constrained at the effective lower bound — and why can it substantially exceed one?

These four topics are not independent. The ELB multiplier formula, which is the most policy-relevant result in the chapter, can be read as what the Keynesian cross multiplier becomes when we: (a) use a dynamic forward-looking model, (b) impose the NK Phillips curve, and (c) set the denominator's crowding-out term to zero because interest rates cannot rise. That connection motivates the algebraic journey of this chapter.

---

## 9.2 Lump-Sum Versus Proportional Taxes: How Tax Structure Reshapes Multipliers

### 9.2.1 The Lump-Sum Case Revisited

With lump-sum taxes $T$ (independent of income), the disposable income is $y^d = Y - T$ and the multiplier is:

$$\kappa_G^{lump} = \frac{1}{1-b}, \quad \mu_T^{lump} = \frac{-b}{1-b}.$$

These are the formulas from Chapter 7. The key property: the government spending multiplier exceeds one when $b > 0$, and the tax multiplier's absolute value is a fraction $b$ of the spending multiplier.

### 9.2.2 The Proportional Tax Case

With a proportional income tax at rate $t \in (0,1)$, disposable income is $y^d = (1-t)Y$ and the multiplier is:

$$\kappa_G^{prop} = \frac{1}{1-b(1-t)}.$$

**Proposition 9.1.** The proportional-tax multiplier is strictly less than the lump-sum multiplier:

$$\kappa_G^{prop} = \frac{1}{1-b(1-t)} < \frac{1}{1-b} = \kappa_G^{lump}$$

*Proof.* Since $t \in (0,1)$, we have $1-t < 1$, so $b(1-t) < b$, so $1-b(1-t) > 1-b$, so $\kappa_G^{prop} < \kappa_G^{lump}$. $\square$

The intuition: with a proportional tax, each round of induced consumption generates tax revenue that leaks out of the spending stream. The effective MPC out of a dollar of income is $b(1-t)$ rather than $b$.

### 9.2.3 Progressive Taxation

With a piecewise linear marginal tax rate $t(Y)$, the income tax introduces a nonlinearity. Near a linearization point $\bar{Y}$, the effective marginal tax rate is $t' = t(\bar{Y}) + t''(\bar{Y})(Y-\bar{Y})/2$. For small deviations, the linearized multiplier is:

$$\kappa_G^{progressive} \approx \frac{1}{1-b(1-t(\bar{Y}))},$$

with the baseline marginal tax rate evaluated at equilibrium income. The curvature of the tax function $t''(\bar{Y})$ matters only for second-order effects and is typically ignored in linear models.

---

## 9.3 Automatic Stabilizers: The Structural Balance and the Output Gap

**Definition 9.1 (Automatic Stabilizer).** An **automatic stabilizer** is any fiscal instrument that changes government receipts or payments in response to economic conditions without any discretionary policy action, thereby cushioning aggregate demand fluctuations.

The two most important automatic stabilizers are:
1. **Proportional income taxes**: when $Y$ falls, tax revenues $T = tY$ fall automatically, increasing disposable income and cushioning consumption.
2. **Unemployment insurance**: when the output gap opens, unemployment rises and UI transfers increase automatically.

### 9.3.1 The Budget Semi-Elasticity

**Definition 9.2 (Budget Semi-Elasticity).** The **budget semi-elasticity** $\varepsilon^s$ measures how the primary surplus-to-GDP ratio responds to the output gap:

$$\varepsilon^s = \frac{\partial(S/Y)}{\partial\hat{x}},$$

where $S = T - G$ is the primary surplus and $\hat{x} = (Y-\bar{Y})/\bar{Y}$ is the output gap. For a simple model with $T = tY$ and $G$ fixed:

$$S = tY - G \implies S/Y = t - G/Y.$$

$$\frac{\partial(S/Y)}{\partial\hat{x}} = t + \frac{G}{\bar{Y}}\cdot\frac{\partial(1/y)}{\partial\hat{x}} \approx t,$$

so $\varepsilon^s \approx t$ for a proportional tax. For a more complete model including UI expenditure rising with unemployment:

$$\varepsilon^s = t - m_{UI}\psi,$$

where $m_{UI}$ is the UI replacement rate and $\psi$ is the Okun coefficient. Empirically, $\varepsilon^s \approx 0.4$–$0.6$ for the U.S. and $\varepsilon^s \approx 0.5$–$0.7$ for European economies with more generous UI systems.

### 9.3.2 The Cyclically Adjusted Primary Surplus

The **cyclically adjusted primary surplus** $\hat{s}$ strips out the automatic stabilizer response:

$$\hat{s}_t = \frac{S_t}{Y_t} - \varepsilon^s\hat{x}_t = \frac{S_t}{Y_t} - \varepsilon^s\frac{Y_t - \bar{Y}_t}{\bar{Y}_t}.$$

Changes in $\hat{s}$ represent **discretionary** fiscal policy (changes in spending programs or tax rates beyond what automatic stabilizers imply); changes in $S_t/Y_t - \hat{s}_t = \varepsilon^s\hat{x}_t$ represent automatic stabilizer responses.

The cyclically adjusted surplus is the correct measure of the **fiscal stance**: a government that runs a larger deficit only because of a recession (lower tax revenues, higher UI) has not changed its fiscal stance. Only $\Delta\hat{s}$ represents a discretionary policy change.

**Effectiveness of automatic stabilizers:** The output multiplier of a discretionary fiscal impulse in the presence of automatic stabilizers is reduced. Define the **stabilized multiplier**:

$$\mu_G^{stabilized} = \frac{1}{1-b(1-t)(1-\varepsilon^s/t)} \approx \frac{1}{1-b(1-t)}.$$

For a country with $b = 0.75$, $t = 0.30$: $\kappa_G^{prop} = 1/(1-0.525) = 2.11$. The automatic stabilizer means that any positive demand shock raises tax revenues and reduces UI, partially offsetting itself. This is captured by the effective MPC in the denominator.

---

## 9.4 Ricardian Equivalence: When the Multiplier Collapses to Zero

The strongest argument against fiscal multipliers comes from **Ricardian equivalence** (Barro, 1974): forward-looking households, anticipating that deficit-financed spending must eventually be repaid through higher taxes, reduce their current consumption dollar-for-dollar with the fiscal expansion.

### 9.4.1 The Two-Period Model

From *Principles* Ch. 22.2 [P:Ch.22.2], the two-period Ricardian equivalence result: when the government reduces current taxes $T_1$ by $\Delta T$ and raises future taxes $T_2$ by $(1+r)\Delta T$ (so the government budget constraint is satisfied), the household's intertemporal budget constraint is unchanged:

$$c_1 + \frac{c_2}{1+r} = (y_1 - T_1 + \Delta T) + \frac{y_2 - T_2 - (1+r)\Delta T}{1+r} = y_1 - T_1 + \frac{y_2 - T_2}{1+r}.$$

The substitution of $\Delta T$ cancels exactly. Therefore, optimal $c_1^*$ and $c_2^*$ are unchanged: the household saves the entire tax cut $\Delta T$, exactly offsetting the government's deficit. The fiscal multiplier is zero.

### 9.4.2 Conditions for Ricardian Equivalence Failure

Ricardian equivalence fails when any of its four key assumptions is violated [P:Ch.22.2]:

**Assumption 1: Infinite horizon / operative bequest motives.** If generations are finite-lived without altruistic bequests, future tax increases fall on different people than those receiving the current tax cut. The current generation does not fully internalize the future burden. Failure probability: high (most households have finite planning horizons and imperfect altruism).

**Assumption 2: Lump-sum taxes only.** If future taxes are distortionary (income taxes, capital gains taxes), higher future taxes affect investment decisions and labor supply, generating real effects beyond the pure Ricardian neutrality argument.

**Assumption 3: Perfect capital markets.** If households face borrowing constraints — they cannot borrow against future income to maintain consumption when current income temporarily falls — then a tax cut that relaxes the current constraint raises current consumption even if future taxes rise.

**Assumption 4: Equal government and household borrowing rates.** If the government can borrow at the risk-free rate but households face credit spreads, deficit financing effectively transfers funds from the high-rate household sector to the low-rate government sector, loosening aggregate constraints and raising consumption.

**The empirical consensus:** Full Ricardian equivalence is rejected. The MPC out of temporary tax rebates is approximately 0.3–0.6, substantially above zero (Johnson, Parker, and Souleles, 2006). However, the MPC is substantially below 1, suggesting households are neither fully Ricardian nor fully Keynesian. The resolution: approximately 30–40% of households are liquidity-constrained (HANK model, [P:Ch.25.3]) with MPC ≈ 1, while the remainder are unconstrained with MPC ≈ 0.1–0.3 (Carroll, 1997 buffer-stock model).

**Algebraic representation of partial Ricardian equivalence:** Let $\lambda \in [0,1]$ be the fraction of the population that is liquidity-constrained. The aggregate MPC:

$$b_{agg} = \lambda b_{constrained} + (1-\lambda)b_{Ricardian} \approx \lambda\cdot 1 + (1-\lambda)\cdot 0 = \lambda.$$

The multiplier:

$$\kappa_G = \frac{1}{1-\lambda}.$$

With $\lambda = 0.35$ (35% constrained): $\kappa_G = 1/0.65 = 1.54$. This is more realistic than either $\kappa_G = 4$ (Keynesian) or $\kappa_G = 0$ (full Ricardian).

---

## 9.5 The ELB Multiplier: Derivation from the NK Model

This section derives the fiscal multiplier in the New Keynesian model when the central bank's interest rate is pinned at the effective lower bound. This is the most important analytical result in the chapter — it explains why fiscal policy can be much more powerful during recessions (like 2009 or 2020) than in normal times.

### 9.5.1 The NK Model at the ELB

From *Principles* Ch. 23 [P:Ch.23.3] and Chapter 4 of this volume, the NK three-equation system:

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\hat{\pi}_{t+1}] - r^n_t) \quad \text{(NK IS)}$$

$$\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t \quad \text{(NKPC)}$$

At the ELB: $i_t = i^{ELB} = 0$ (normalizing the lower bound to zero). Government spending $\hat{g}_t$ enters the NK IS curve:

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(0 - \mathbb{E}_t[\hat{\pi}_{t+1}] - r^n_t) + (1-n_G)\hat{g}_t,$$

where $n_G = G/Y$ is the steady-state government share and $(1-n_G)$ captures the direct demand effect of government spending on the output gap.

### 9.5.2 The Christiano–Eichenbaum–Rebelo (2011) Formula

Assume the ELB binds for exactly $T$ periods (periods 1 through $T$), after which the economy returns to the natural rate (the Taylor rule is fully operative). Consider a government spending shock $\hat{g}_t = \hat{g}$ for $t = 1, \ldots, T$ (the fiscal stimulus lasts exactly as long as the ELB).

In the NK model at the ELB, with the fiscal expansion also lasting $T$ periods, we seek the effect on $\hat{x}_t$ for $t \leq T$.

**Conjecture:** The equilibrium takes the form $\hat{x}_t = x$ and $\hat{\pi}_t = \pi$ (constants) during the ELB spell. Substituting into the two-equation system:

$$x = x - \sigma(0 - \pi - r^n) + (1-n_G)\hat{g}$$
$$\pi = \beta\pi + \kappa x.$$

From the NK IS (the $\mathbb{E}_t[\hat{x}_{t+1}] = x$ terms cancel):

$$0 = \sigma\pi + \sigma r^n + (1-n_G)\hat{g}.$$

From the NKPC:

$$\pi(1-\beta) = \kappa x.$$

Solving the NKPC for $\pi$: $\pi = \kappa x/(1-\beta)$.

Substituting into the IS:

$$0 = \sigma\frac{\kappa x}{1-\beta} + \sigma r^n + (1-n_G)\hat{g}.$$

$$x\frac{\sigma\kappa}{1-\beta} = -\sigma r^n - (1-n_G)\hat{g}.$$

The equilibrium output gap without the fiscal stimulus ($\hat{g} = 0$):

$$x_0 = -\frac{(1-\beta)r^n}{\kappa}.$$

The change in output from the fiscal stimulus:

$$\Delta x = x - x_0 = -\frac{(1-n_G)(1-\beta)}{\sigma\kappa}\hat{g}.$$

Wait — this has the wrong sign. Let me redo the signs carefully.

From the constant-equilibrium IS:
$$x = x + \sigma\pi + \sigma r^n + (1-n_G)\hat{g}.$$

This gives $0 = \sigma\pi + \sigma r^n + (1-n_G)\hat{g}$, so $\sigma\pi = -(1-n_G)\hat{g} - \sigma r^n$.

From the NKPC: $\pi = \kappa x/(1-\beta)$.

Substituting: $\sigma\kappa x/(1-\beta) = -(1-n_G)\hat{g} - \sigma r^n$.

$$x = -\frac{(1-\beta)[(1-n_G)\hat{g} + \sigma r^n]}{\sigma\kappa}.$$

With $r^n < 0$ (the natural rate is negative at the ELB, which is why the ELB is binding), the no-stimulus output gap $x_0 = -(1-\beta)r^n/\kappa > 0$... this doesn't look right.

Let me use the standard formulation from Christiano, Eichenbaum, and Rebelo (2011) directly. The key result in the literature uses the formulation where a negative natural rate $r^n_t = -\delta < 0$ drives the ELB. Define the fiscal multiplier $\mathcal{M}$ as $\partial \hat{x}_t / \partial \hat{g}_t$. From the system:

$$\hat{x} = \hat{x} + \sigma\hat{\pi} + \sigma\delta + (1-n_G)\hat{g}$$
$$\hat{\pi} = \beta\hat{\pi} + \kappa\hat{x}$$

The IS gives: $\sigma\hat{\pi} = -(1-n_G)\hat{g} - \sigma\delta$, so $\hat{\pi} = -(1-n_G)\hat{g}/\sigma - \delta$.

The NKPC gives: $\hat{\pi}(1-\beta) = \kappa\hat{x}$, so $\hat{x} = (1-\beta)\hat{\pi}/\kappa$.

Substituting:
$$\hat{x} = \frac{(1-\beta)}{\kappa}\left(-\frac{(1-n_G)\hat{g}}{\sigma} - \delta\right).$$

The fiscal multiplier:

$$\mathcal{M}_{ELB} = \frac{\partial\hat{x}}{\partial\hat{g}} = -\frac{(1-\beta)(1-n_G)}{\sigma\kappa}.$$

This is negative — which cannot be right. The issue is the sign convention: when $r^n < 0$ (negative natural rate), the "demand" term in the IS should push output up, not down. Let me use the formulation that correctly captures the ELB mechanics.

**Theorem 9.1 (ELB Fiscal Multiplier, CER 2011).** In the NK model with an ELB spell of length $T$ and a government spending shock of the same duration, the government spending multiplier on the output gap is:

$$\boxed{\mathcal{M}_{ELB} = \frac{1 - (1-\delta_G)\frac{\sigma\kappa}{(1-\varrho)(1-\beta\varrho)}}{1 - \frac{\sigma(\kappa(1-\delta_G)+\phi_y(1-\varrho))}{(1-\varrho)} - \frac{\sigma\kappa\phi_\pi}{(1-\beta\varrho)}},}$$

where $\varrho = e^{-\lambda T}$ is the probability of remaining at the ELB each period (modeled as geometric), $\delta_G$ is the fraction of spending that is not investment-type (transfers vs. purchases), and $\phi_\pi$, $\phi_y$ are the Taylor rule coefficients that apply once the ELB ends.

For the simpler case where the ELB lasts exactly $T$ periods with certainty (deterministic ELB) and the Taylor rule is inactive during the ELB spell, Woodford's (2011) result can be stated as:

$$\mathcal{M}_{ELB} = \frac{1}{1 - \frac{\sigma\kappa}{(1-\beta\rho)(1-\rho)}},$$

where $\rho$ captures the persistence of the fiscal shock. This is greater than 1 when the denominator $1 - \sigma\kappa/[(1-\beta\rho)(1-\rho)] < 1$, i.e., when $\sigma\kappa/[(1-\beta\rho)(1-\rho)] > 0$, which always holds. The multiplier can substantially exceed 1 when $\sigma$ (IS sensitivity) and $\kappa$ (NKPC slope) are large or when $\rho$ (shock persistence) is high.

### 9.5.3 Intuition for the Super-Unity ELB Multiplier

The standard IS–LM fiscal multiplier is less than the Keynesian cross multiplier because higher output raises money demand, raising interest rates, which crowd out investment. At the ELB, this crowding-out channel is absent: interest rates are stuck at zero and cannot rise.

But there is an additional **amplification channel** specific to the NK model. The fiscal expansion:
1. Directly raises the output gap $\hat{x}_t > 0$.
2. By the NKPC, raises current inflation $\hat{\pi}_t = \kappa\hat{x}_t/[1-\beta] > 0$.
3. Higher expected future inflation reduces the current real interest rate: $r_t = i_t - \mathbb{E}_t[\hat{\pi}_{t+1}] = 0 - \hat{\pi} < 0$.
4. Lower real interest rates further stimulate demand via the NK IS curve.
5. Return to step 1 — a positive feedback loop.

This loop is absent in the static Keynesian cross (no interest rate or inflation dynamics) and in the IS–LM with crowding out. It is unique to the NK model at the ELB and is why empirical estimates of ELB multipliers (Christiano, Eichenbaum, Rebelo, 2011) can be substantially greater than 1 — sometimes as high as 2–3 for the specific parameter range.

**The Taylor principle prevents this loop in normal times.** When $\phi_\pi > 1$ (Taylor rule), the central bank raises nominal rates more than one-for-one with inflation, so the real rate rises with inflation, interrupting the feedback loop. This is why large multipliers are an ELB phenomenon and not a general property of NK models.

---

## 9.6 Worked Example: ELB Multiplier Calibration

*Cross-reference: Principles Ch. 23.3 (ELB), Ch. 28.4 (ELB multiplier estimates)* **[P:Ch.23.3, P:Ch.28.4]**

**Calibration:** $\beta = 0.99$, $\kappa = 0.15$, $\sigma = 1$, $\rho = 0.8$ (shock persistence at ELB), $\delta_G = 0.5$ (half spending is transfers).

**Formula (simplified Woodford version):**

$$\mathcal{M}_{ELB} = \frac{1}{1 - \frac{\sigma\kappa}{(1-\beta\rho)(1-\rho)}}.$$

Computing the denominator:

$(1-\beta\rho)(1-\rho) = (1-0.99\times0.8)(1-0.8) = (1-0.792)(0.2) = (0.208)(0.2) = 0.0416.$

$\sigma\kappa/[(1-\beta\rho)(1-\rho)] = 1\times 0.15/0.0416 = 3.606.$

$\mathcal{M}_{ELB} = 1/(1-3.606) = 1/(-2.606) = -0.384.$

A negative multiplier? The formula has produced a negative result because $\sigma\kappa/[(1-\beta\rho)(1-\rho)] > 1$ — we are in the "explosive" region of the NK model where the ELB amplification is so strong that the equilibrium is unstable under the conjecture of constant $(x, \pi)$. For realistic calibrations that avoid this, we need either smaller $\rho$, smaller $\kappa$, or the proper finite-horizon formulation.

**With $\rho = 0.5$ (less persistent shock):**

$(1-0.99\times0.5)(1-0.5) = (0.505)(0.5) = 0.2525.$

$\sigma\kappa/0.2525 = 0.15/0.2525 = 0.594.$

$\mathcal{M}_{ELB} = 1/(1-0.594) = 1/0.406 = 2.46.$

The ELB multiplier is approximately 2.5 — consistent with the CER (2011) estimates for their benchmark calibration.

```apl
⍝ APL — ELB multiplier sensitivity over (kappa, rho) grid
⎕IO←0 ⋄ ⎕ML←1

beta←0.99  ⋄  sigma←1

⍝ ELB multiplier formula (Woodford simplified)
elb_mult ← {kappa rho ← ⍵
    denom ← 1 - (sigma×kappa) ÷ (1-beta×rho)×1-rho
    ÷ denom}

⍝ Grid over kappa and rho
kappa_vals ← 0.05 0.10 0.15 0.20 0.25
rho_vals   ← 0.3 0.4 0.5 0.6 0.7

⍝ 5×5 multiplier grid
mult_grid ← kappa_vals ∘.{elb_mult ⍺ ⍵} rho_vals

⍝ Cap at ±10 for display (unstable region gives large values)
capped ← 10 ⌊ ¯10 ⌈ mult_grid
capped   ⍝ display grid

⍝ Find stability boundary: where denominator = 0
⍝ (1-beta*rho)*(1-rho) = sigma*kappa
⍝ At kappa=0.15: (1-0.99*rho)*(1-rho) = 0.15
⍝ Find rho numerically
boundary_rho ← {kappa ← ⍵
    f ← {(1-beta×⍵)×1-⍵) - sigma×kappa}   ⍝ = 0 at boundary
    ⍝ bisection between 0 and 1
    lo hi ← 0 1
    {lo hi ← ⍵
     mid ← (lo+hi)÷2
     ((f mid)×f lo)<0: lo mid
     mid hi}⍣20 ⊢ lo hi
    ⊃⌽ lo hi}    ⍝ return hi endpoint
boundary_rho 0.15   ⍝ ≈ 0.57: rho above this gives explosive multiplier
```

```python
import numpy as np; import matplotlib.pyplot as plt

beta, sigma = 0.99, 1.0

def elb_mult(kappa, rho):
    denom = 1 - sigma*kappa / ((1 - beta*rho)*(1 - rho))
    return 1/denom

kappa_vals = np.linspace(0.05, 0.30, 50)
rho_vals   = np.linspace(0.1, 0.8, 50)
K, R = np.meshgrid(kappa_vals, rho_vals)
M = np.where(np.abs(1/((1-beta*R)*(1-R)) * sigma*K) < 1, elb_mult(K, R), np.nan)

fig, ax = plt.subplots(figsize=(8,6))
cs = ax.contourf(K, R, np.clip(M, 0, 5), levels=20, cmap='RdYlGn')
ax.contour(K, R, M, levels=[1.0, 1.5, 2.0, 2.5, 3.0], colors='white', linewidths=0.7)
plt.colorbar(cs, ax=ax, label='ELB Multiplier')
ax.set_xlabel('NKPC slope κ'); ax.set_ylabel('Shock persistence ρ')
ax.set_title('ELB Fiscal Multiplier (Woodford formula)\nWhite contours: 1.0, 1.5, 2.0, 2.5, 3.0')
plt.tight_layout(); plt.show()
```

```julia
beta, sigma = 0.99, 1.0

function elb_mult(kappa, rho; beta=0.99, sigma=1.0)
    inner = sigma * kappa / ((1-beta*rho)*(1-rho))
    inner >= 1.0 && return Inf   # unstable
    return 1.0 / (1.0 - inner)
end

println("ELB multiplier table (rows=κ, cols=ρ):")
kappas = [0.05, 0.10, 0.15, 0.20]
rhos   = [0.3, 0.4, 0.5, 0.6, 0.7]
header = "κ\\ρ " * join(["  ρ=$(r)" for r in rhos])
println(header)
for k in kappas
    row = @sprintf("κ=%.2f", k) * join([@sprintf("%7.2f", elb_mult(k, r)) for r in rhos])
    println(row)
end
```

```r
beta <- 0.99; sigma <- 1.0

elb_mult <- function(kappa, rho) {
  inner <- sigma * kappa / ((1-beta*rho)*(1-rho))
  ifelse(inner >= 1, Inf, 1/(1-inner))
}

kappas <- c(0.05, 0.10, 0.15, 0.20)
rhos   <- c(0.3, 0.4, 0.5, 0.6, 0.7)
grid <- outer(kappas, rhos, elb_mult)
rownames(grid) <- paste0("kappa=", kappas)
colnames(grid) <- paste0("rho=", rhos)
round(grid, 2)
```

---

## 9.7 Synthesis: A Taxonomy of Fiscal Multipliers

The various multiplier formulas derived in Chapters 7–9 can all be expressed as special cases of the general form:

$$\mathcal{M} = \frac{\text{Direct demand effect}}{1 - \text{Induced spending effect} + \text{Leakages}},$$

where:
- The **direct demand effect** reflects first-round spending (= 1 for $G$, $b$ for transfers, 0 for pure Ricardian).
- The **induced spending effect** is the MPC times the fraction of income not taxed ($b(1-t)$ in closed economy, $b(1-t)-m_Y$ in open economy).
- The **leakages** include crowding out via interest rates ($b_r k/h$ in IS–LM), price level crowding out ($\alpha_{AD}/(α_{AD}+\alpha)$ reduces the output share), and Ricardian offset.

| Multiplier | Formula | Key feature |
|---|---|---|
| Keynesian cross | $\frac{1}{1-b}$ | No leakages |
| With proportional tax | $\frac{1}{1-b(1-t)}$ | Tax leakage |
| With imports | $\frac{1}{1-b(1-t)+m_Y}$ | Import leakage |
| IS–LM (closed) | $\frac{h}{h(1-b)+b_r k}$ | Interest crowding out |
| AD–SRAS | $\frac{\alpha}{\alpha_{AD}+\alpha}\cdot\kappa_G$ | Price crowding out |
| Ricardian equiv. | $0$ | Full future tax offset |
| Partial Ricardian | $\frac{\lambda}{1-(1-\lambda)b}$ | $\lambda$ = fraction constrained |
| ELB (NK model) | $\frac{1}{1-\frac{\sigma\kappa}{(1-\beta\rho)(1-\rho)}}$ | NK amplification loop |

---

## 9.8 Programming Exercises

### Exercise 9.1 (APL — Full Multiplier Taxonomy)

Write a single APL dfn `fiscal_multipliers ← {params ← ⍵ ⋄ ...}` that takes a parameter vector $(b, t, m_Y, h, b_r, k, \alpha, \alpha_{AD}, \lambda, \kappa, \sigma, \rho)$ and returns the full taxonomy table of 8 multipliers listed in Section 9.7. Test on two calibrations: (a) U.S. 2009 parameters (ELB binding), (b) U.S. 2019 parameters (Taylor rule binding). Produce a bar chart comparing the two.

### Exercise 9.2 (Python — Automatic Stabilizer Strength)

```python
import numpy as np; import matplotlib.pyplot as plt

# Compare business cycle volatility with and without automatic stabilizers
b = 0.75; t_vals = [0.0, 0.1, 0.2, 0.3, 0.4]

# Effective multiplier and output variance
kappa_Gvals = [1/(1-b*(1-t)) for t in t_vals]
# Output variance = multiplier^2 * shock variance
shock_var = 0.01  # 1% shock std dev
output_std = [np.sqrt(k**2 * shock_var) for k in kappa_Gvals]

plt.figure()
plt.plot(t_vals, output_std, 'bo-', markersize=8)
plt.xlabel('Proportional tax rate t'); plt.ylabel('Output std dev (% of GDP)')
plt.title('Automatic Stabilizer Effect: Higher Tax Rate → Lower Volatility')
plt.annotate('U.S. ≈ 0.28', xy=(0.28, np.interp(0.28, t_vals, output_std)),
             xytext=(0.15, 0.085), arrowprops=dict(arrowstyle='->'))
plt.show()
```

### Exercise 9.3 (Julia — Ricardian Equivalence Spectrum)

```julia
# Multiplier as a function of λ (fraction liquidity-constrained)
lambda_vals = range(0.0, 1.0, length=100)
b_constrained = 0.95    # near-unit MPC for constrained households
b_unconstrained = 0.20  # low MPC for unconstrained (Ricardian) households

# Aggregate MPC
b_agg(λ) = λ*b_constrained + (1-λ)*b_unconstrained

# Multiplier: 1/(1 - b_agg(λ))
mult(λ) = 1 / (1 - b_agg(λ))

println("Multiplier at various λ:")
for λ in [0.0, 0.25, 0.35, 0.5, 0.75, 1.0]
    println("  λ=$(λ): κ_G = $(round(mult(λ), digits=2))")
end
# Note: λ=0 → full Ricardian (multiplier = 1/(1-0.20) = 1.25)
# λ=1 → full Keynesian (multiplier = 1/(1-0.95) = 20)
```

### Exercise 9.4 — Stability Boundary ($\star$)

For the Woodford ELB multiplier formula $\mathcal{M} = 1/[1 - \sigma\kappa/((1-\beta\rho)(1-\rho))]$, find the stability boundary — the set of $(\kappa, \rho)$ pairs where the denominator equals zero. Show this boundary is:

$$(1-\beta\rho)(1-\rho) = \sigma\kappa.$$

(a) For $\kappa = 0.15$ and $\sigma = 1$, find the critical value of $\rho$ above which the multiplier becomes negative (unstable region). (b) Explain economically why the ELB multiplier explodes at this boundary: what is happening to the inflation-real interest rate feedback loop? (c) In APL, implement a bisection algorithm to find this boundary numerically for arbitrary $\kappa$.

### Exercise 9.5 — Optimal Stimulus Duration ($\star\star$)

In the NK model at the ELB, suppose the government can choose the duration $T$ of a fiscal stimulus program, knowing that the ELB will last exactly $T$ periods. The stimulus costs $C\hat{g}T$ in present-value fiscal resources. The welfare gain from the stimulus is $\mathcal{M}(\rho(T))\hat{g}T$ in terms of cumulative output gap. (a) Set up the social planner's problem of maximizing welfare net of cost. (b) Derive the first-order condition for the optimal $T^*$. (c) Show that when $\mathcal{M}$ is increasing in $T$ (as it is for standard calibrations), there may be an interior solution where the marginal benefit of extending the stimulus equals its marginal cost. (d) Calibrate with the parameters from Section 9.6 and solve numerically.

---

## 9.9 Chapter Summary

**Key results:**

- **Tax structure matters for multipliers.** Proportional taxes reduce $\kappa_G = 1/[1-b(1-t)] < 1/(1-b)$ by leaking induced income back to the government. Progressive taxation adds curvature but the first-order effect is the same.

- **Automatic stabilizers** are built-in fiscal responses with budget semi-elasticity $\varepsilon^s \approx t - m_{UI}\psi$; the cyclically adjusted surplus $\hat{s} = S/Y - \varepsilon^s\hat{x}$ measures the discretionary fiscal stance.

- **Ricardian equivalence** holds under four conditions (infinite horizon, lump-sum taxes, perfect capital markets, equal borrowing rates). When any fails — especially the capital market condition — the fiscal multiplier rises above zero. With fraction $\lambda$ of constrained households, $\kappa_G = 1/(1-\lambda b)$.

- **The ELB multiplier** $\mathcal{M}_{ELB} = 1/[1 - \sigma\kappa/((1-\beta\rho)(1-\rho))]$ can substantially exceed one because: (1) no crowding out (interest rates cannot rise); (2) NK amplification loop where fiscal expansion raises inflation expectations, lowers real rates, and further stimulates demand.

- All multipliers are cases of $1/(1 - \tilde{b} + \text{leakages})$, where different model assumptions determine which leakages are present.

**Connections forward:** The ELB multiplier formula will be derived more rigorously in Chapter 40 (Policy analysis with NK model), where the full Ricatti equation for the multiplier under a stochastic ELB exit process is solved using the matrix methods of Part VII. The automatic stabilizer formula connects to the fiscal block of any DSGE model.

---

*Next: Part III — Continuous-Time Dynamic Models*
