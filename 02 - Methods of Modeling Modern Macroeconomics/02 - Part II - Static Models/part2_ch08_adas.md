# Chapter 8: The AD–AS Model in Equation Form

*Solving for Equilibrium Price Level and Output*

> *"The AS–AD model is the macroeconomist's workhorse: simple enough to solve analytically, rich enough to address the central questions of stabilization policy."*

**Cross-reference:** *Principles* Ch. 7 (AS–AD derivation and shock analysis); Ch. 10 (Phillips curve, EAPC, NKPC); Ch. 30 (inflation and deflation) **[P:Ch.7, P:Ch.10, P:Ch.30]**

---

## 8.1 From Curves to Equations: The AD–AS System

*Principles* Chapter 7 presented the AS–AD framework graphically — a downward-sloping AD curve in $(Y, P)$ space and an upward-sloping SRAS — with qualitative analysis of demand and supply shocks. This chapter makes that analysis algebraic. We write explicit equations for AD, SRAS, and LRAS; solve for the short-run and long-run equilibria in closed form; and derive exact formulas for the effects of shocks on output and the price level.

The payoff of the algebraic treatment is twofold. First, it produces precise multipliers for the price level and output analogous to the income multipliers of Chapter 7. Second, it derives the **sacrifice ratio** — the output cost per unit of inflation reduction — algebraically from the structure of the model, providing the benchmark against which the empirical estimates of *Principles* Ch. 10 and Ch. 28 are evaluated.

---

## 8.2 The Aggregate Demand Equation

The AD equation is derived from the IS–LM system by eliminating the interest rate. From Chapter 6:

$$Y^*(P) = \frac{h\kappa_G\bar{A} + \beta_r(M/P)}{h + \beta_r k}.$$

This expresses equilibrium IS–LM output as a function of the price level $P$ (through the real money supply $M/P$). It is the AD equation.

**Definition 8.1 (Aggregate Demand Equation).** The AD equation expresses equilibrium output as a decreasing function of the price level:

$$Y = \frac{h\kappa_G\bar{A}}{h+\beta_r k} + \frac{\beta_r M}{(h+\beta_r k)P} \equiv A_0 + \frac{A_1}{P},$$

where $A_0 = h\kappa_G\bar{A}/(h+\beta_r k) > 0$ captures fiscal and autonomous demand components, and $A_1 = \beta_r M/(h+\beta_r k) > 0$ captures the monetary component. For a simpler linear approximation, log-linearize around a baseline price level $P_0$:

$$Y \approx A_0 + \frac{A_1}{P_0} - \frac{A_1}{P_0^2}(P - P_0) = \bar{Y}^{AD} - \alpha_{AD}(P - P_0),$$

where $\alpha_{AD} = A_1/P_0^2 = \beta_r M/(P_0^2(h+\beta_r k)) > 0$ is the slope of the linearized AD curve.

**What shifts AD?**

A change in any variable other than $P$ that affects IS or LM shifts the AD curve:
- Fiscal expansion ($\Delta G > 0$): raises $A_0$ by $h\kappa_G\Delta G/(h+\beta_r k) = \mu_G\Delta G$; AD shifts right.
- Monetary expansion ($\Delta M > 0$): raises $A_1/P$ at every $P$; AD shifts right.
- Consumer confidence shock ($\Delta a > 0$): raises $\bar{A}$, shifts AD right.

**The slope of the AD curve** in $(Y, P)$ space: $dP/dY|_{AD} = -P_0^2/(A_1) < 0$ (steep when $A_1$ is small — e.g., when $\beta_r$ is small or monetary policy is tight).

---

## 8.3 The Short-Run Aggregate Supply Equations

*Principles* Chapter 7 presented two microeconomic foundations for the upward-sloping SRAS: sticky wages and the Lucas imperfect-information model. We write both as explicit linear equations.

### 8.3.1 Sticky-Wage Model

When the nominal wage is fixed at $\bar{W}$ (by contracts), a higher price level $P$ reduces the real wage $\bar{W}/P$, inducing firms to hire more labour and produce more. From the firm's profit-maximization condition $F_N(K, N) = \bar{W}/P$ and inverting to get $N = N(\bar{W}/P)$, then substituting into $Y = F(K, N)$:

$$Y = Y^s(\bar{W}/P) \implies \frac{\partial Y^s}{\partial P} = -\frac{Y^s_N F_N}{P F_{NN}} > 0.$$

Linearizing around the flex-price equilibrium $(\bar{Y}, \bar{P})$ where $\bar{P} = P^e$ (expected price):

$$\boxed{Y = \bar{Y} + \alpha_{SW}(P - P^e), \quad \alpha_{SW} > 0.}$$

This is the **SRAS equation**: output exceeds potential $\bar{Y}$ when the actual price level $P$ exceeds the price level expected when wages were set $P^e$. Unexpected inflation erodes the real wage, stimulating production.

### 8.3.2 Lucas Imperfect-Information Supply Curve

*Principles* Chapter 16 derived the Lucas (1973) aggregate supply curve [P:Ch.16.1]:

$$Y = \bar{Y} + \alpha(P - P^e), \quad \alpha = \gamma\frac{\sigma_z^2}{\sigma_z^2 + \sigma_\eta^2},$$

where $\sigma_z^2$ is the variance of relative price shocks and $\sigma_\eta^2$ is the variance of aggregate nominal shocks. The parameter $\alpha$ falls when aggregate nominal volatility $\sigma_\eta^2$ rises: in high-inflation environments, firms attribute observed price increases to nominal shocks rather than real demand, reducing the output response. Both sticky-wage and Lucas supply functions take the same linear form $Y = \bar{Y} + \alpha(P - P^e)$.

**Definition 8.2 (SRAS Equation).** The **short-run aggregate supply (SRAS) equation** is:
$$Y = \bar{Y} + \alpha(P - P^e), \quad \alpha > 0,$$
where $\bar{Y}$ is potential output, $P$ is the actual price level, and $P^e$ is the price level expected when prices/wages were set.

**Definition 8.3 (LRAS Equation).** The **long-run aggregate supply (LRAS) equation** is:
$$Y = \bar{Y},$$
a vertical line at potential output. In the long run $P^e \to P$ (expectations adjust), and the SRAS collapses to the LRAS.

---

## 8.4 Short-Run Equilibrium: Solving the AD–SRAS System

Short-run equilibrium requires simultaneous clearing of goods and money markets (AD) and labour markets with sticky wages or imperfect information (SRAS). We have two equations in two unknowns $(Y^*, P^*)$:

$$Y = \bar{Y}^{AD} - \alpha_{AD}(P - P_0) \quad \text{(AD)}$$
$$Y = \bar{Y} + \alpha(P - P^e) \quad \text{(SRAS)}$$

Equating:
$$\bar{Y}^{AD} - \alpha_{AD}(P^* - P_0) = \bar{Y} + \alpha(P^* - P^e).$$

Solving for $P^*$:
$$P^*(\alpha_{AD} + \alpha) = \bar{Y}^{AD} - \bar{Y} + \alpha_{AD}P_0 + \alpha P^e.$$

$$\boxed{P^* = \frac{(\bar{Y}^{AD} - \bar{Y}) + \alpha_{AD}P_0 + \alpha P^e}{\alpha_{AD} + \alpha}.}$$

Substituting back to find $Y^*$:
$$Y^* = \bar{Y} + \alpha(P^* - P^e) = \bar{Y} + \frac{\alpha(\bar{Y}^{AD} - \bar{Y} + \alpha_{AD}(P_0 - P^e))}{\alpha_{AD} + \alpha}.$$

**Definition 8.4 (Output Multiplier in the AD–SRAS System).** The effect of a rightward AD shift $\Delta\bar{Y}^{AD}$ on equilibrium output:

$$\frac{\partial Y^*}{\partial\bar{Y}^{AD}} = \frac{\alpha}{\alpha_{AD} + \alpha}.$$

This is always strictly less than 1. When AD shifts right, the price level rises along the upward-sloping SRAS; the higher price level partially crowds out the demand stimulus through the real-balance effect (which is built into the AD slope $\alpha_{AD}$). The larger is $\alpha$ (steeper SRAS), the more of the AD shift goes into prices and the less into output.

**Two limiting cases:**

1. **Perfectly flexible prices** ($\alpha \to \infty$, vertical SRAS): $\partial Y^*/\partial\bar{Y}^{AD} \to 0$. AD shifts are entirely absorbed by price changes; output is unchanged at $\bar{Y}$. This is the classical result.

2. **Perfectly sticky prices** ($\alpha \to 0$, horizontal SRAS): $\partial Y^*/\partial\bar{Y}^{AD} \to 1$. AD shifts translate one-for-one into output changes; prices do not move. This is the extreme Keynesian result, coinciding with the Keynesian cross.

---

## 8.5 Shock Analysis: Analytical Comparative Statics

### 8.5.1 Demand Shock

A positive aggregate demand shock (e.g., fiscal expansion $\Delta G$) raises $\bar{Y}^{AD}$ by $\mu_G\Delta G$:

$$\Delta Y^* = \frac{\alpha}{\alpha_{AD}+\alpha}\mu_G\Delta G, \quad \Delta P^* = \frac{\alpha_{AD}}{\alpha_{AD}+\alpha}\frac{\mu_G\Delta G}{\alpha}.$$

Wait — let me re-express $\Delta P^*$ cleanly. From the $P^*$ formula:
$$\Delta P^* = \frac{\Delta\bar{Y}^{AD}}{\alpha_{AD}+\alpha} = \frac{\mu_G\Delta G}{\alpha_{AD}+\alpha}.$$

Both output and the price level rise. The fraction going to output vs. prices is determined by $\alpha/(\alpha_{AD}+\alpha)$ and $\alpha_{AD}/(\alpha_{AD}+\alpha)$ respectively.

### 8.5.2 Supply Shock

A negative supply shock (oil price increase) raises firms' marginal costs, shifting SRAS left: $\bar{Y} \to \bar{Y} - \Delta\bar{Y}^{supply}$ in the SRAS equation. The comparative statics:

$$\Delta Y^* = -\frac{\alpha_{AD}}{\alpha_{AD}+\alpha}\Delta\bar{Y}^{supply} < 0,$$
$$\Delta P^* = \frac{\alpha}{\alpha_{AD}+\alpha}\frac{\Delta\bar{Y}^{supply}}{\alpha} = \frac{\Delta\bar{Y}^{supply}}{\alpha_{AD}+\alpha} > 0.$$

A supply shock generates **stagflation**: output falls while prices rise. This is the algebraic confirmation of the *Principles* Chapter 7 graphical result [P:Ch.7.4]. Note the supply shock moves output and prices in *opposite* directions from a demand shock — a key identifying feature that empirical researchers exploit.

**The policy dilemma from a supply shock:** The central bank faces a trade-off. It can:
- **Accommodate** the supply shock (increase $M$, shifting AD right): stabilizes output at $Y = \bar{Y}$ but tolerates the higher price level.
- **Fight inflation** (decrease $M$, shifting AD left): stabilizes the price level but allows output to fall further.

This trade-off is captured by the **divine coincidence breakdown** of *Principles* Ch. 23 [P:Ch.23.2]: in the standard NK model without cost-push shocks, optimal policy achieves both stable inflation and a zero output gap simultaneously (divine coincidence). With supply shocks, this is impossible.

---

## 8.6 Long-Run Adjustment: The Price Adjustment Mechanism

In the short run, $P^e$ is predetermined. In the medium run, expectations adjust toward the realized price level. The adjustment mechanism:

$$\dot{P}^e = \pi(P - P^e), \quad \pi > 0 \text{ (speed of adjustment)}.$$

Substituting the SRAS equilibrium condition $P^* = P^e + (Y^* - \bar{Y})/\alpha$:

$$\dot{P}^e = \pi(P^* - P^e) = \pi\frac{Y^* - \bar{Y}}{\alpha}.$$

In long-run equilibrium, $\dot{P}^e = 0$, which requires $Y^* = \bar{Y}$ — output equals potential. The long-run price level $\bar{P}$ is then pinned by the AD equation:

$$\bar{Y} = \bar{Y}^{AD} - \alpha_{AD}(\bar{P} - P_0) \implies \bar{P} = P_0 + \frac{\bar{Y}^{AD} - \bar{Y}}{\alpha_{AD}}.$$

The adjustment path following a demand shock: output initially jumps above $\bar{Y}$, expectations gradually rise, SRAS shifts left, and output returns to $\bar{Y}$ at a permanently higher price level. This is the medium-run neutrality of aggregate demand shocks: real output returns to potential while the price level adjusts permanently. The speed of return is governed by $\pi$ — faster price adjustment means faster return to potential.

---

## 8.7 The Sacrifice Ratio: Algebraic Derivation

The **sacrifice ratio** measures the cumulative output loss required to reduce inflation by one percentage point. It is derived directly from the expectations-augmented Phillips curve [P:Ch.10.2].

The EAPC (in continuous time):

$$\dot{\pi} = -\alpha(u - u^*) + \varepsilon,$$

where $u$ is the unemployment rate, $u^*$ the natural rate, and $\varepsilon$ a supply shock term. Using Okun's Law $Y - \bar{Y} = -\psi(u - u^*)$ [P:Ch.3.4]:

$$\dot{\pi} = \frac{\alpha}{\psi}(Y - \bar{Y}) + \varepsilon.$$

**Definition 8.5 (Sacrifice Ratio).** The **sacrifice ratio** $SR$ is the cumulative output loss per unit of permanent inflation reduction:

$$SR = \frac{\int_0^\infty (\bar{Y} - Y_t)\,dt}{\Delta\pi_{permanent}}.$$

**Theorem 8.1 (Analytical Sacrifice Ratio from EAPC).** Under adaptive expectations and a disinflationary policy that holds $Y_t < \bar{Y}$ until the target inflation $\pi^{target}$ is reached, the sacrifice ratio is:

$$SR = \frac{\psi}{\alpha}.$$

*Proof.* The EAPC gives $d\pi = (\alpha/\psi)(Y - \bar{Y})dt$. Integrating over the disinflation:

$$\pi^{target} - \pi^{initial} = \frac{\alpha}{\psi}\int_0^T(Y_t - \bar{Y})dt.$$

The sacrifice ratio is the negative of this: $SR = -\int_0^T(Y_t-\bar{Y})dt/(\pi^{target}-\pi^{initial}) = \psi/\alpha$. $\square$

With $\alpha \approx 0.3$–$0.5$ (EAPC coefficient, U.S. estimates) and $\psi \approx 2$ (Okun coefficient):

$$SR = \frac{2}{0.3\text{ to }0.5} \approx 4 \text{ to } 7.$$

This is the output loss (in percent-years of GDP) per percentage point of inflation reduction — consistent with Ball (1994)'s empirical estimates of $SR \approx 1.4$–$2.8$ for fast disinflations and higher for slow ones [P:Ch.10.3].

**The New Keynesian modification:** Under rational expectations and the NKPC, the sacrifice ratio can be zero (costless credible disinflation) because forward-looking price-setters immediately revise inflation expectations downward upon a credible announcement. The hybrid NKPC (with backward-looking component $\omega$) gives a sacrifice ratio between 0 and $\psi/\alpha$.

---

## 8.8 The Dynamic Extension: Preview of the NK Three-Equation System

The static AD–SRAS framework is the limiting case of the dynamic New Keynesian model. This section makes the connection explicit, bridging to the derivations of Part VII.

In the static model, we have two equations in $(Y, P)$. The NK three-equation model [P:Ch.7.5] replaces these with three dynamic equations in $(\hat{x}, \hat{\pi}, i)$:

$$\hat{x}_t = \mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\hat{\pi}_{t+1}] - r_t^n) \quad \text{(NK IS)}$$

$$\hat{\pi}_t = \beta\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\hat{x}_t + u_t \quad \text{(NKPC)}$$

$$i_t = r^n + \pi^* + \phi_\pi(\hat{\pi}_t - \pi^*) + \phi_y\hat{x}_t \quad \text{(Taylor rule)}$$

**Mapping static to dynamic:** The static AD equation $Y = \bar{Y}^{AD} - \alpha_{AD}P$ is the static ($\mathbb{E}_t[\hat{x}_{t+1}] = 0$, $\sigma$ finite) special case of the NK IS curve: $\hat{x} = -\sigma(\phi_\pi\hat{\pi}+\phi_y\hat{x}) \Rightarrow \hat{x}(1+\sigma\phi_y) = -\sigma\phi_\pi\hat{\pi}$, giving $\hat{x} = -[\sigma\phi_\pi/(1+\sigma\phi_y)]\hat{\pi}$ — a negative relationship between output gap and inflation, analogous to AD.

The static SRAS $Y = \bar{Y} + \alpha(P - P^e)$ maps to the NKPC with $\beta \to 1$ and static expectations: $\hat{\pi} = \kappa\hat{x} + u$, where $\alpha \leftrightarrow 1/\kappa$.

**The key difference:** The NK model is forward-looking; the static model is static. In the NK model, the sacrifice ratio depends on the credibility of the disinflation (which determines whether $\mathbb{E}_t[\hat{\pi}_{t+1}]$ responds to the announcement), while in the static model it is simply $\psi/\alpha$.

---

## 8.9 Worked Example: Demand Shock with Full Price Adjustment

*Cross-reference: Principles Ch. 7.4 (demand and supply shock analysis)* **[P:Ch.7.4]**

**Setup:** An economy at potential output $\bar{Y} = 1000$ and equilibrium price level $P^e = P^* = 100$ experiences a positive demand shock from a fiscal expansion: $\Delta\bar{Y}^{AD} = 80$ (the AD curve shifts right by 80 units of output at any price level).

**Parameters:** $\alpha_{AD} = 0.5$ (AD slope), $\alpha = 2$ (SRAS slope).

**Short-run equilibrium:**

$$\Delta Y^* = \frac{\alpha}{\alpha_{AD}+\alpha}\Delta\bar{Y}^{AD} = \frac{2}{0.5+2}\times 80 = \frac{2}{2.5}\times 80 = 64.$$

$$\Delta P^* = \frac{\Delta\bar{Y}^{AD}}{\alpha_{AD}+\alpha} = \frac{80}{2.5} = 32.$$

So $Y^* = 1064$, $P^* = 132$. Output rises 64 (above potential by 64) and prices rise 32.

**Long-run adjustment:** In the long run, $P^e$ adjusts to $P^* = 132$. The SRAS shifts left until $Y = \bar{Y} = 1000$. The final price level is pinned by the new AD at $Y = 1000$:

$$\bar{P}^{LR} = P_0 + \frac{\Delta\bar{Y}^{AD}}{\alpha_{AD}} = 100 + \frac{80}{0.5} = 100 + 160 = 260.$$

Wait — let us recompute. The AD equation at $Y = 1000$: $1000 = (1000 + 80) - 0.5(P - 100) \Rightarrow 0.5(P-100) = 80 \Rightarrow P_{LR} = 100 + 160 = 260$.

So in the long run, output returns to 1000 but prices rise from 100 to 260 — a 160-unit increase in the price level, compared to the short-run impact of only 32. The fiscal expansion has no permanent output effect but generates substantial permanent inflation.

```apl
⍝ APL — AD-SRAS equilibrium solution
⎕IO←0 ⋄ ⎕ML←1

⍝ Parameters
alpha_AD ← 0.5   ⍝ AD slope
alpha    ← 2     ⍝ SRAS slope
Ybar     ← 1000  ⍝ potential output
Ybar_AD  ← 1000  ⍝ baseline AD intercept (Y at P=100)
Pe       ← 100   ⍝ expected price level
P0       ← 100   ⍝ baseline price level

⍝ Short-run equilibrium: solve AD = SRAS simultaneously
⍝ AD:   Y = Ybar_AD - alpha_AD*(P - P0)
⍝ SRAS: Y = Ybar + alpha*(P - Pe)
⍝ In matrix form: [1, alpha_AD; 1, -alpha] [Y; P] = [Ybar_AD + alpha_AD*P0; Ybar - alpha*Pe]

Amat ← 2 2 ⍴ 1 alpha_AD 1 (-alpha)
bvec ← (Ybar_AD + alpha_AD×P0) (Ybar - alpha×Pe)
eq   ← bvec ⌹ Amat
Y_sr ← eq[0]
P_sr ← eq[1]
Y_sr  ⍝ ≈ 1000 (baseline: at potential)
P_sr  ⍝ ≈ 100

⍝ Demand shock: Ybar_AD increases by 80
dYAD ← 80
Ybar_AD_new ← Ybar_AD + dYAD
bvec_new ← (Ybar_AD_new + alpha_AD×P0) (Ybar - alpha×Pe)
eq_new   ← bvec_new ⌹ Amat
Y_sr_new ← eq_new[0]
P_sr_new ← eq_new[1]
Y_sr_new  ⍝ ≈ 1064
P_sr_new  ⍝ ≈ 132

⍝ Long-run: Y = Ybar, solve AD for P
P_lr ← P0 + (Ybar_AD_new - Ybar) ÷ alpha_AD
P_lr  ⍝ ≈ 260

⍝ Multipliers
dY_dYAD ← alpha ÷ alpha_AD + alpha           ⍝ output multiplier
dP_dYAD ← 1 ÷ alpha_AD + alpha               ⍝ price multiplier (SR)
dP_dYAD_LR ← 1 ÷ alpha_AD                   ⍝ price multiplier (LR)
dY_dYAD  ⍝ 0.8
dP_dYAD  ⍝ 0.4  (SR)
dP_dYAD_LR ⍝ 2.0 (LR: all shock goes to prices)
```

```python
import numpy as np

def adas_equilibrium(alpha_AD, alpha, Ybar, Ybar_AD, Pe, P0):
    """Solve AD-SRAS for short-run (Y*, P*)."""
    A = np.array([[1, alpha_AD], [1, -alpha]])
    b = np.array([Ybar_AD + alpha_AD*P0, Ybar - alpha*Pe])
    return np.linalg.solve(A, b)

# Baseline
sol0 = adas_equilibrium(0.5, 2, 1000, 1000, 100, 100)
print(f"Baseline: Y* = {sol0[0]:.1f}, P* = {sol0[1]:.1f}")

# After demand shock: Ybar_AD -> 1080
sol1 = adas_equilibrium(0.5, 2, 1000, 1080, 100, 100)
print(f"After shock (SR): Y* = {sol1[0]:.1f}, P* = {sol1[1]:.1f}")

# Long-run: Y = 1000, solve for P from new AD
P_lr = 100 + (1080 - 1000)/0.5
print(f"Long-run: Y* = 1000, P_LR = {P_lr:.1f}")

# Multipliers
alpha_AD, alpha = 0.5, 2
print(f"\nOutput multiplier (SR): {alpha/(alpha_AD+alpha):.3f}")
print(f"Price multiplier (SR):  {1/(alpha_AD+alpha):.3f}")
print(f"Price multiplier (LR):  {1/alpha_AD:.3f}")
```

---

## 8.10 Programming Exercises

### Exercise 8.1 (APL — Supply Shock Stagflation)

Extend the APL code from Section 8.9 to compute the short-run and long-run effects of a negative supply shock $\Delta\bar{Y}^{supply} = -60$ (SRAS shifts left by 60). Verify: (a) output falls and prices rise (stagflation); (b) the effects on output and prices from a supply shock have opposite signs to those from a demand shock; (c) the long-run price level increases less from a supply shock than from an equivalent demand shock (because the policy accommodation is partial).

### Exercise 8.2 (Python — Phase Diagram)

```python
import numpy as np; import matplotlib.pyplot as plt

alpha_AD, alpha, Ybar, P0, Pe = 0.5, 2, 1000, 100, 100

# AD and SRAS curves as functions of P
P_range = np.linspace(50, 300, 500)
Y_AD    = lambda Ybar_AD, P: Ybar_AD - alpha_AD * (P - P0)
Y_SRAS  = lambda Pe_, P: Ybar + alpha * (P - Pe_)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: baseline and shock
ax = axes[0]
Ybar_AD_vals = [1000, 1080]
colors = ['blue', 'red']
labels = ['Baseline AD', 'AD after demand shock']
for Yad, c, lab in zip(Ybar_AD_vals, colors, labels):
    ax.plot(Y_AD(Yad, P_range), P_range, color=c, label=lab)
ax.plot(Y_SRAS(100, P_range), P_range, 'g-', label='SRAS (Pe=100)')
ax.plot(Y_SRAS(132, P_range), P_range, 'g--', label='SRAS after adjustment (Pe=132)')
ax.axvline(Ybar, color='k', linestyle=':', label='LRAS')
ax.set_xlim(900, 1150); ax.set_ylim(50, 300)
ax.set_xlabel('Y'); ax.set_ylabel('P')
ax.legend(fontsize=8); ax.set_title('Demand Shock: SR and LR')

# Right: sacrifice ratio as function of alpha (SRAS slope)
ax2 = axes[1]
alpha_vals = np.linspace(0.5, 5, 100)
psi = 2  # Okun coefficient
SR = psi / alpha_vals
ax2.plot(alpha_vals, SR)
ax2.axhline(2.0, linestyle='--', color='red', label='Empirical range bottom')
ax2.axhline(3.5, linestyle='--', color='orange', label='Empirical range top')
ax2.set_xlabel('α (SRAS slope / Phillips curve coefficient)')
ax2.set_ylabel('Sacrifice ratio SR = ψ/α')
ax2.set_title('Sacrifice Ratio vs. Price Flexibility')
ax2.legend()
plt.tight_layout(); plt.show()
```

### Exercise 8.3 (Julia — Sacrifice Ratio)

```julia
# Sacrifice ratio from EAPC with different expectations schemes
psi = 2.0  # Okun coefficient
alpha = 0.35  # EAPC slope

# Under adaptive expectations: SR = psi/alpha
SR_adaptive = psi / alpha

# Under hybrid NKPC with backward weight omega:
# SR = psi/alpha * (omega fraction)
# (forward-looking part reduces sacrifice ratio)
omega_vals = 0:0.1:1.0
SR_hybrid = psi ./ alpha .* omega_vals  # linear interpolation

println("SR under adaptive expectations: $(round(SR_adaptive, digits=2))")
println("\nSR under hybrid NKPC:")
for (ω, sr) in zip(omega_vals, SR_hybrid)
    println("  ω=$(round(ω,digits=1)): SR=$(round(sr,digits=2))")
end
println("\n(ω=0 means fully forward-looking → SR=0; ω=1 → fully adaptive → SR=$(round(SR_adaptive,digits=2)))")
```

### Exercise 8.4 (R — Policy Accommodation of Supply Shock)

```r
alpha_AD <- 0.5; alpha <- 2; Ybar <- 1000; P0 <- 100; Pe <- 100

# Supply shock: Ybar falls by 60
delta_supply <- 60
Ybar_shocked <- Ybar - delta_supply

# No accommodation: AD unchanged
A <- matrix(c(1, 1, alpha_AD, -alpha), 2, 2)
b_no_acc <- c(1000 + alpha_AD*P0, Ybar_shocked - alpha*Pe)
eq_no_acc <- solve(A, b_no_acc)

# Full accommodation: AD expanded to keep Y = Ybar
# Need Ybar_AD such that SR equilibrium gives Y=Ybar
# From SR formula: Ybar = Ybar_AD_new - alpha/(alpha_AD+alpha)*delta_supply
Ybar_AD_accommodate <- Ybar + (alpha_AD+alpha)/alpha * delta_supply
b_acc <- c(Ybar_AD_accommodate + alpha_AD*P0, Ybar_shocked - alpha*Pe)
eq_acc <- solve(A, b_acc)

cat("Supply shock: no accommodation\n")
cat(sprintf("  Y* = %.1f, P* = %.1f\n", eq_no_acc[1], eq_no_acc[2]))
cat("Supply shock: full accommodation (Y = Ybar)\n")
cat(sprintf("  Y* = %.1f, P* = %.1f\n", eq_acc[1], eq_acc[2]))
cat("Trade-off: accommodation prevents output loss but accepts higher price level\n")
```

### Exercise 8.5 — AS–AD Matrix Generalization ($\star$)

Generalize the AD–SRAS system to include a third equation: the price-adjustment mechanism $\dot{P}^e = \pi(P^* - P^e)$. Write the full three-dimensional system as an ODE $\dot{\mathbf{z}} = f(\mathbf{z})$ where $\mathbf{z} = (Y, P, P^e)'$. Linearize around the long-run equilibrium and compute the eigenvalues. Show that one eigenvalue is zero (reflecting the long-run neutrality of the price level), one is negative (stable convergence), and the Jacobian structure implies that the economy always returns to potential after demand shocks.

### Exercise 8.6 — Optimal Monetary Response to Supply Shock ($\star\star$)

Following a negative supply shock that reduces $\bar{Y}$ by $\delta$, the central bank must choose how much to accommodate (expand $M$). Model the bank's loss function as $\mathcal{L} = \lambda(Y-\bar{Y})^2 + (1-\lambda)(P-\bar{P})^2$ where $\bar{P}$ is the price target. (a) Derive the optimal accommodation level $\Delta M^*$ as a function of $\lambda$, $\delta$, and the model parameters. (b) Show that $\lambda = 0$ (pure inflation targeter) implies no accommodation ($\Delta M = 0$), while $\lambda = 1$ (pure output stabilizer) implies full accommodation. (c) For $\lambda = 0.5$ and the calibration from this chapter, compute $\Delta M^*$ numerically.

---

## 8.11 Chapter Summary

**Key results:**

- The **AD equation** $Y = A_0 + A_1/P$ is derived from IS–LM by substituting money-market equilibrium; linearized around $P_0$ it gives $Y \approx \bar{Y}^{AD} - \alpha_{AD}(P - P_0)$.
- The **SRAS equation** $Y = \bar{Y} + \alpha(P - P^e)$ arises from either sticky wages or Lucas imperfect information; the LRAS equation is $Y = \bar{Y}$.
- The **AD–SRAS short-run equilibrium** is a $2\times2$ linear system with unique solution: $\Delta Y^* = [\alpha/(\alpha_{AD}+\alpha)]\Delta\bar{Y}^{AD}$; $\Delta P^* = \Delta\bar{Y}^{AD}/(\alpha_{AD}+\alpha)$.
- **Demand shocks** raise both output and prices; **supply shocks** raise prices and lower output (stagflation). The sign of the output-price correlation distinguishes demand from supply shocks empirically.
- The **sacrifice ratio** $SR = \psi/\alpha$ is the output cost (in percent-years) per point of inflation reduction under adaptive expectations; under the NKPC it approaches zero for fully credible disinflation.
- In APL: the AD–SRAS system is the same `⌹` pattern as IS–LM, with a $2\times2$ matrix; the long-run price level requires only scalar arithmetic.

**Connections forward:** Chapter 9 derives the ELB analogue of these multipliers for the NK model. Part VII (Chapter 27–31) derives the AD–SRAS system from explicit household and firm optimization via log-linearization of the DSGE model.

---

*Next: Chapter 9 — Tax and Government Spending Multipliers: Analytical Solutions with Behavioral Foundations*
