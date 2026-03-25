# Part II: Static Macroeconomic Models and Comparative Statics

*Connects to: Principles Parts I–II*

---

The models of *Principles* Part II were presented qualitatively and graphically. The IS curve was a downward-sloping line; the LM curve an upward-sloping one; their intersection determined equilibrium — and shifts of those curves were analyzed by moving pictures in the $(Y, i)$ plane. That graphical approach builds invaluable intuition, but it has limits. It cannot tell you the *exact* size of the fiscal multiplier as a function of the model parameters. It cannot handle the open-economy extension (Mundell–Fleming) with three markets simultaneously without becoming unwieldy. And it provides no leverage for comparative statics — the precise derivative of equilibrium output with respect to, say, the money supply.

This part strips away the geometry and replaces it with algebra. Every model from *Principles* Parts I–II is translated into a system of equations, solved explicitly using the linear algebra tools of Chapter 2, and subjected to systematic comparative statics. The payoffs are:

- **Exact multiplier formulas.** Not "the output rises" but $\Delta Y^* = h/(h + b_r k)\cdot\Delta G$, from which we can immediately read off when fiscal policy is powerful (large $h$, small $b_r k$) and when it is impotent (small $h$, or $b_r k \to \infty$).
- **Policy counterfactuals.** Given a specific calibration (say, the U.S. economy in 2009), we can compute exactly what the model predicts about the ARRA multiplier, not merely what direction it points.
- **Transition to dynamic models.** The static multipliers derived here are the limiting cases of the dynamic impulse responses we will compute in Parts V–VII. Seeing the connection makes both clearer.

The four chapters cover: the IS–LM model (Chapter 6), the Keynesian cross and multiplier algebra (Chapter 7), the AS–AD model (Chapter 8), and the full taxonomy of fiscal multipliers including the ELB formula (Chapter 9).

---

# Chapter 6: Solving the IS–LM Model

*Matrix Inversion and Cramer's Rule for Policy Analysis*

> *"The IS–LM model is simple enough to solve in closed form and rich enough to generate most of the qualitative results of more sophisticated models."*

**Cross-reference:** *Principles* Ch. 9 (IS–LM derivation and policy analysis); Ch. 22 (fiscal policy); Ch. 23 (monetary policy); Ch. 21 (open economy, Mundell–Fleming) **[P:Ch.9, P:Ch.22, P:Ch.23, P:Ch.21]**

---

## 6.1 From Curves to Equations

*Principles* Chapter 9 derived the IS and LM curves graphically and used them to analyze fiscal and monetary policy. Here we write those curves as explicit linear equations, stack them into a matrix system, and use the tools of Chapter 2 to solve for equilibrium and compute all policy multipliers in closed form.

The strategy has three steps. First, write the IS and LM conditions as two equations in two unknowns $(Y, i)$. Second, solve the $2\times2$ system by matrix inversion or Cramer's rule. Third, differentiate the solution with respect to every policy variable to obtain the multipliers. Steps one and two require the algebra of Chapter 2; step three uses the implicit function theorem of Chapter 1.

### 6.1.1 The IS Equation

The IS curve is the goods-market equilibrium condition. In the closed economy, it requires that aggregate expenditure $\mathcal{E}$ equals aggregate output $Y$:

$$Y = C(Y - T) + I(r) + G.$$

Using a linear consumption function $C(y^d) = a + bY^d$ with marginal propensity to consume $b \in (0,1)$ and autonomous consumption $a > 0$, and the linear investment function $I(r) = \bar{I} - b_r r$ with sensitivity $b_r > 0$:

$$Y = a + b(Y - T) + \bar{I} - b_r r + G.$$

Rearranging, and using $r \approx i - \pi^e$ with $\pi^e$ exogenous (so that $r$ and $i$ differ by the constant $\pi^e$ which we absorb into parameters):

$$Y - bY = a - bT + \bar{I} - b_r i + G \implies Y(1-b) = \bar{A} - b_r i,$$

where $\bar{A} = a - bT + \bar{I} + G + b_r\pi^e$ collects all autonomous expenditure. The **IS equation** is thus:

$$Y + \frac{b_r}{1-b}i = \frac{\bar{A}}{1-b}, \quad \text{or equivalently} \quad Y = \frac{\bar{A}}{1-b} - \frac{b_r}{1-b}i.$$

Writing more compactly with $\kappa_G \equiv 1/(1-b)$ (the Keynesian multiplier) and letting $\beta_r \equiv b_r/(1-b)$:

$$\boxed{Y = \kappa_G\bar{A} - \beta_r i.}$$

This is the IS equation: a downward-sloping relationship in $(Y, i)$ space with slope $-1/\beta_r$ (or equivalently, when plotted with $i$ on the vertical axis, slope $-1/\beta_r$ in $i$ per unit of $Y$).

### 6.1.2 The LM Equation

The LM curve is the money-market equilibrium condition. Real money demand $L(Y, i) = kY - hi$ (where $k > 0$ is income elasticity and $h > 0$ is interest semi-elasticity) equals real money supply $M/P$:

$$kY - hi = \frac{M}{P}.$$

Solving for $i$:

$$\boxed{i = \frac{k}{h}Y - \frac{1}{h}\frac{M}{P}.}$$

This is the LM equation: an upward-sloping relationship in $(Y, i)$ space with slope $k/h$.

---

## 6.2 The IS–LM System in Matrix Form

Stacking the IS and LM equations in the form $A\mathbf{y} = \mathbf{b}$:

$$\underbrace{\begin{pmatrix} 1 & \beta_r \\ k & -h \end{pmatrix}}_{A}\underbrace{\begin{pmatrix} Y \\ i \end{pmatrix}}_{\mathbf{y}} = \underbrace{\begin{pmatrix} \kappa_G\bar{A} \\ -M/P \end{pmatrix}}_{\mathbf{b}}.$$

**Definition 6.1 (IS–LM Coefficient Matrix).** The matrix $A = \begin{pmatrix} 1 & \beta_r \\ k & -h \end{pmatrix}$ encodes the structural parameters of the IS–LM model. It is nonsingular (invertible) when $\det(A) \neq 0$.

Computing the determinant:
$$\det(A) = (1)(-h) - (\beta_r)(k) = -h - \beta_r k = -(h + \beta_r k) < 0.$$

Since $h > 0$, $\beta_r > 0$, and $k > 0$, the determinant is strictly negative, so $A$ is always invertible and the IS–LM system always has a unique solution. This is the algebraic counterpart to the graphical statement that the IS and LM curves always intersect in exactly one point (given upward-sloping LM and downward-sloping IS).

### 6.2.1 Direct Matrix Inversion

The inverse of the $2\times2$ matrix $A$:

$$A^{-1} = \frac{1}{\det(A)}\begin{pmatrix} -h & -\beta_r \\ -k & 1 \end{pmatrix} = \frac{1}{-(h+\beta_r k)}\begin{pmatrix} -h & -\beta_r \\ -k & 1 \end{pmatrix}.$$

The solution $\mathbf{y}^* = A^{-1}\mathbf{b}$:

$$\begin{pmatrix} Y^* \\ i^* \end{pmatrix} = \frac{1}{h + \beta_r k}\begin{pmatrix} h & \beta_r \\ k & -1 \end{pmatrix}\begin{pmatrix} \kappa_G\bar{A} \\ -M/P \end{pmatrix}.$$

Working out the products:

$$\boxed{Y^* = \frac{h\kappa_G\bar{A} + \beta_r(M/P)}{h + \beta_r k}}$$

$$\boxed{i^* = \frac{k\kappa_G\bar{A} - (M/P)}{h + \beta_r k}}$$

These are the **reduced-form equations** of the IS–LM model: they express the endogenous variables $(Y^*, i^*)$ entirely in terms of exogenous variables $(\bar{A}, M/P)$ and structural parameters $(b, b_r, k, h)$.

---

## 6.3 Policy Multipliers via Cramer's Rule

**Definition 6.2 (Policy Multiplier).** The **policy multiplier** with respect to instrument $z$ is the partial derivative $\partial Y^*/\partial z$ (or $\partial i^*/\partial z$), holding all other exogenous variables fixed.

We derive every policy multiplier from the reduced-form solution. Note that $\bar{A} = a - bT + \bar{I} + G + b_r\pi^e$, so $\partial\bar{A}/\partial G = 1$, $\partial\bar{A}/\partial T = -b$, and $\partial\bar{A}/\partial\bar{I} = 1$.

### 6.3.1 The Fiscal Multiplier

$$\mu_G = \frac{\partial Y^*}{\partial G} = \frac{h\kappa_G}{h + \beta_r k} = \frac{h/(1-b)}{h + b_r k/(1-b)} = \frac{h}{h(1-b) + b_r k}.$$

Writing $b_r = b_r$ (not divided by $1-b$) and $\kappa_G = 1/(1-b)$:

$$\boxed{\mu_G = \frac{h}{h + b_r k} \cdot \frac{1}{1-b} = \frac{h}{h(1-b) + b_r k}.}$$

This is always strictly less than $\kappa_G = 1/(1-b)$, the Keynesian cross multiplier. The ratio $\mu_G/\kappa_G = h(1-b)/(h(1-b)+b_r k) < 1$ measures how much of the Keynesian cross multiplier survives after accounting for the interest rate feedback — the **crowding-out factor**.

**Definition 6.3 (Crowding Out).** The reduction in private investment caused by the rise in the interest rate following a fiscal expansion is called **crowding out**. In the IS–LM model, a fiscal expansion of $\Delta G$ raises the interest rate by:

$$\Delta i^* = \frac{k}{h + \beta_r k}\cdot\kappa_G\Delta G = \frac{k}{h(1-b) + b_r k}\Delta G > 0,$$

which reduces investment by $b_r\Delta i^* > 0$. The gross fiscal effect on output ($\kappa_G\Delta G$) is partially offset by the investment crowding-out ($-b_r\Delta i^* = -b_r k\kappa_G\Delta G/(h + \beta_r k)$):

$$\Delta Y^* = \underbrace{\kappa_G\Delta G}_{\text{gross IS shift}} - \underbrace{\frac{b_r k\kappa_G}{h + \beta_r k}\Delta G}_{\text{crowding out}} = \frac{h}{h + \beta_r k}\kappa_G\Delta G.$$

**Theorem 6.1 (Crowding-Out Formula).** In the IS–LM model, the fraction of the Keynesian cross multiplier that survives crowding out is:

$$\text{Crowding-out factor} = 1 - \frac{b_r k}{h + b_r k} = \frac{h}{h + b_r k}.$$

*Proof.* Direct substitution of the two multiplier expressions above. $\square$

**Two limiting cases:**

1. **Liquidity trap:** $h \to \infty$. Then $\mu_G \to \kappa_G$ — the IS–LM multiplier equals the Keynesian cross multiplier. When money demand is perfectly elastic, the interest rate doesn't rise in response to fiscal expansion (the LM is horizontal), so there is zero crowding out.

2. **Classical case:** $h \to 0$ (money demand completely interest-inelastic). Then $\mu_G \to 0$ — fiscal policy is completely ineffective. The LM is vertical; any income increase raises money demand, which drives up the interest rate until investment falls by exactly $\Delta G$. Complete crowding out.

### 6.3.2 The Tax Multiplier

$$\mu_T = \frac{\partial Y^*}{\partial T} = \frac{h\kappa_G(-b)}{h + \beta_r k} = -\frac{bh}{h(1-b) + b_r k}.$$

$$\boxed{\mu_T = -\frac{bh}{h(1-b) + b_r k} = -b\cdot\mu_G/h \cdot h = -\frac{b}{1-b}\cdot\frac{h}{h+b_rk/(1-b)}}$$

More cleanly: $\mu_T = -b/(1-b) \cdot (h/(h+b_r k)) \cdot (1-b) = -bh/[h(1-b)+b_rk]$. Since $|\mu_T| = b\cdot|\mu_G/(1-b)| \cdot (1-b) < |\mu_G|$, we confirm $|\mu_T| < \mu_G$: tax cuts are less stimulative than equivalent spending increases (each dollar of spending creates a full dollar of first-round demand; each dollar of tax cut is only partially spent).

### 6.3.3 The Monetary Policy Multiplier

A monetary expansion increases the real money supply $M/P$:

$$\mu_M = \frac{\partial Y^*}{\partial(M/P)} = \frac{\beta_r}{h + \beta_r k} = \frac{b_r/(1-b)}{h + b_r k/(1-b)}.$$

$$\boxed{\mu_M = \frac{b_r}{h(1-b) + b_r k}.}$$

Monetary policy works by reducing the interest rate ($\partial i^*/\partial(M/P) = -1/(h+\beta_rk) < 0$), which stimulates investment, which raises income through the multiplier. The monetary multiplier is zero when $b_r = 0$ (investment insensitive to interest rates — the IS is vertical) or when $h \to \infty$ (liquidity trap).

**Summary table of IS–LM multipliers:**

| Instrument | Multiplier on $Y^*$ | Sign | Vanishes when |
|---|---|---|---|
| $G$ (spending) | $\frac{h}{h(1-b)+b_r k}$ | $+$ | $h\to 0$ (classical LM) |
| $T$ (taxes) | $\frac{-bh}{h(1-b)+b_r k}$ | $-$ | $h\to 0$ or $b=0$ |
| $M/P$ (money) | $\frac{b_r}{h(1-b)+b_r k}$ | $+$ | $b_r=0$ (vertical IS) or $h\to\infty$ (trap) |
| $\bar{I}$ (invest.) | $\frac{h}{h(1-b)+b_r k}$ | $+$ | Same as $G$ |

---

## 6.4 Comparative Statics via the Implicit Function Theorem

The multiplier formulas above were derived from the explicit reduced form. An alternative, more general method uses the IFT [M:Ch.1.5] directly on the equilibrium conditions without first solving for $Y^*$ explicitly.

The equilibrium conditions define an implicit function $F(Y, i; G, T, M/P) = \mathbf{0}$ where:

$$F_1 = Y - \kappa_G\bar{A} + \beta_r i = 0 \quad \text{(IS)}$$
$$F_2 = kY - hi - M/P = 0 \quad \text{(LM)}$$

The Jacobian with respect to endogenous variables:
$$\frac{\partial\mathbf{F}}{\partial(Y,i)} = \begin{pmatrix} 1 & \beta_r \\ k & -h \end{pmatrix} = A.$$

By the IFT (Theorem 1.8, multivariate version):

$$\frac{d(Y^*, i^*)}{d\mathbf{z}} = -A^{-1}\frac{\partial\mathbf{F}}{\partial\mathbf{z}},$$

where $\mathbf{z}$ is any exogenous variable vector. For $z = G$:

$$\frac{\partial\mathbf{F}}{\partial G} = \begin{pmatrix} -\kappa_G \\ 0 \end{pmatrix}.$$

Therefore:
$$\begin{pmatrix}\partial Y^*/\partial G \\ \partial i^*/\partial G\end{pmatrix} = -A^{-1}\begin{pmatrix}-\kappa_G \\ 0\end{pmatrix} = \frac{1}{h+\beta_r k}\begin{pmatrix}h\kappa_G \\ k\kappa_G\end{pmatrix},$$

recovering $\mu_G = h\kappa_G/(h+\beta_r k)$ and $\partial i^*/\partial G = k\kappa_G/(h+\beta_r k) > 0$. The IFT method is particularly convenient when the system has many variables and one does not want to invert the full matrix algebraically.

---

## 6.5 The Mundell–Fleming Model: A 3×3 System

The open economy adds the foreign exchange market (BP curve) to the IS–LM system, yielding a 3×3 linear system determining $(Y^*, i^*, e^*)$ where $e$ is the exchange rate (or the current-account balance, depending on the exchange rate regime).

### 6.5.1 The Three Equations

**IS (open economy):** Aggregate demand includes net exports $NX(Y, Y^*, e) = X(Y^*, e) - M^{imports}(Y, e)$:

$$Y = \kappa^{open}\bar{A}^{open} - \beta_r i + \beta_e e,$$

where $\kappa^{open} = 1/(1-b+m_Y)$ (with $m_Y$ the marginal propensity to import) and $\beta_e > 0$ captures the expenditure-switching effect: a higher $e$ (weaker domestic currency) improves net exports.

**LM (same as closed):**

$$kY - hi = M/P.$$

**BP (balance of payments):** For the capital account $KA = \kappa_i(i - i^*)$ where $i^*$ is the world interest rate, and trade balance $TB = TB_0 + \theta_e e - m_Y Y$, the BP equilibrium $TB + KA = 0$:

$$-m_Y Y + \theta_e e + \kappa_i(i - i^*) = 0 \implies -m_Y Y + \kappa_i i + \theta_e e = \kappa_i i^*.$$

The 3×3 system $A\mathbf{y} = \mathbf{b}$ with $\mathbf{y} = (Y, i, e)'$:

$$\underbrace{\begin{pmatrix} 1 & \beta_r & -\beta_e \\ k & -h & 0 \\ -m_Y & \kappa_i & \theta_e \end{pmatrix}}_{A}\begin{pmatrix} Y \\ i \\ e \end{pmatrix} = \underbrace{\begin{pmatrix} \kappa^{open}\bar{A}^{open} \\ M/P \\ \kappa_i i^* \end{pmatrix}}_{\mathbf{b}}.$$

The solution $\mathbf{y}^* = A^{-1}\mathbf{b}$ requires inverting this $3\times3$ matrix — straightforward with APL's `⌹` but algebraically complex by hand.

### 6.5.2 Exchange Rate Regimes and the Trilemma

The Mundell–Fleming model takes its most elegant form under the two polar exchange rate regimes, because each regime eliminates one endogenous variable.

**Fixed exchange rate ($e = \bar{e}$ fixed, $M/P$ endogenous):** The exchange rate is no longer a variable; the BP equation pins down the money supply required to maintain $\bar{e}$. The system reduces to a $2\times2$ problem in $(Y, i)$ with $e = \bar{e}$ substituted. The fiscal multiplier under fixed rates with perfect capital mobility ($\kappa_i \to\infty$):

$$\mu_G^{fixed,\, perfect} = \kappa^{open} = \frac{1}{1-b+m_Y}.$$

There is **no crowding out** under a fixed exchange rate with perfect capital mobility! The BP curve forces $i = i^*$; any tendency for the domestic interest rate to rise above $i^*$ immediately attracts capital inflows, expanding the money supply and keeping $i = i^*$. The fiscal expansion faces no interest-rate headwind.

**Flexible exchange rate ($e$ endogenous, $M/P$ fixed):** Under perfect capital mobility, the capital inflows triggered by a fiscal expansion appreciate the exchange rate ($e$ falls, domestic currency strengthens), crowding out net exports until $\Delta Y^* = 0$. **Fiscal policy is completely ineffective.** But monetary policy — which works through the exchange rate — is fully potent.

**Theorem 6.2 (Mundell–Fleming Trilemma, Algebraic Form).** Under perfect capital mobility ($\kappa_i \to\infty$), the fiscal multiplier satisfies:
$$\mu_G^{fixed} = \kappa^{open} > 0, \quad \mu_G^{flexible} = 0.$$

*Proof sketch.* Under perfect capital mobility, $i = i^*$ at all times. Under a fixed exchange rate, this constraint is satisfied by adjusting $M$, so the LM curve accommodates the fiscal expansion. Under a flexible exchange rate, $M$ is fixed, so the LM curve cannot shift — the only way to maintain $i = i^*$ is through exchange rate appreciation which eliminates the net-export stimulus exactly. $\square$

---

## 6.6 Worked Example: Full Solution Under Fixed and Flexible Rates

*Cross-reference: Principles Ch. 21 (exchange rates, UIP, trilemma)* **[P:Ch.21]**

**Calibration:**
- $b = 0.75$ (MPC), $m_Y = 0.15$ (marginal propensity to import), $b_r = 2$ (investment-interest sensitivity)
- $h = 4$ (interest semi-elasticity of money demand), $k = 0.5$ (income elasticity of money demand)
- $\beta_e = 1.5$ (exchange rate effect on net exports), $\theta_e = 0.8$, $\kappa_i = 10$ (capital mobility)
- $\bar{A}^{open} = 400$, $M/P = 500$, $i^* = 0.04$

**IS–LM solution (closed economy, for reference):**

$\kappa_G = 1/(1-0.75) = 4$, $\beta_r = 2/0.25 = 8$

$Y^* = (4 \times 400 \times h + 8 \times 500)/(h + 8 \times 0.5) = (6400 + 4000)/8 = 10400/8 = 1300$

$i^* = (0.5 \times 4 \times 400 - 500)/(4 + 8 \times 0.5) = (800 - 500)/8 = 37.5$

(Note: $i^*$ in these units is not a realistic interest rate — the calibration is illustrative.)

**Fiscal multipliers:**

$\mu_G^{closed} = h/(h(1-b)+b_r k) = 4/(4\times0.25 + 2\times0.5) = 4/(1+1) = 2.0$

$\mu_G^{open,\,fixed} \approx \kappa^{open} = 1/(1-0.75+0.15) = 1/0.40 = 2.5$

$\mu_G^{open,\,flexible} \approx 0$ (complete crowding out via exchange rate appreciation)

**Interpretation:** Moving from a closed economy to a fixed-exchange-rate open economy raises the fiscal multiplier from 2.0 to 2.5, because the fixed rate prevents the interest rate from rising (capital inflows finance the deficit), eliminating crowding out. Under flexible rates, the multiplier collapses to zero — consistent with *Principles* Ch. 21 [P:Ch.21.4].

```apl
⍝ APL — IS-LM solution and multipliers
⎕IO←0 ⋄ ⎕ML←1

⍝ Parameters
b←0.75  ⋄  b_r←2  ⋄  k←0.5  ⋄  h←4
kg ← ÷1-b                          ⍝ Keynesian multiplier
beta_r ← b_r × kg                  ⍝ β_r = b_r/(1-b)

⍝ IS-LM coefficient matrix and RHS
A ← 2 2 ⍴ 1 beta_r k (-h)
Abar ← 400  ⋄  MP ← 500
b_vec ← (kg×Abar) (-MP)

⍝ Solve: [Y*, i*] = A⁻¹ × b
eq ← b_vec ⌹ A
Y_star ← eq[0]
i_star ← eq[1]
Y_star  ⍝ equilibrium output
i_star  ⍝ equilibrium interest rate

⍝ All multipliers from the formula
denom ← h×(1-b) + b_r×k
mu_G  ← h     ÷ denom   ⍝ spending multiplier
mu_T  ← (-b×h)÷ denom   ⍝ tax multiplier
mu_M  ← b_r   ÷ denom   ⍝ monetary multiplier
mu_G  ⍝ ≈ 2.0
mu_T  ⍝ ≈ ¯1.5
mu_M  ⍝ ≈ 0.5

⍝ Parameter sensitivity: grid of mu_G over (h, b_r) values
H    ← 0.5 1 2 4 8 ∞
BR   ← 0.5 1 2 4
grid ← H ∘.{⍺÷(⍺×1-b)+⍵×k} BR    ⍝ mu_G for each (h,b_r) pair
grid   ⍝ 6×4 multiplier table
```

```python
# Python — IS-LM system: solve and compute all multipliers
import numpy as np
from sympy import symbols, Matrix, simplify, latex

# Numeric solution
b, b_r, k, h = 0.75, 2.0, 0.5, 4.0
kg = 1/(1-b); beta_r = b_r*kg
Abar, MP = 400, 500

A = np.array([[1, beta_r], [k, -h]])
rhs = np.array([kg*Abar, -MP])
Y_star, i_star = np.linalg.solve(A, rhs)
print(f"Y* = {Y_star:.2f},  i* = {i_star:.4f}")

denom = h*(1-b) + b_r*k
mu_G = h/denom; mu_T = -b*h/denom; mu_M = b_r/denom
print(f"μ_G = {mu_G:.3f},  μ_T = {mu_T:.3f},  μ_M = {mu_M:.3f}")

# Symbolic multipliers via sympy
b_s, br_s, k_s, h_s = symbols('b b_r k h', positive=True)
denom_s = h_s*(1-b_s) + br_s*k_s
print("Symbolic μ_G =", simplify(h_s/denom_s))
```

```julia
# Julia — IS-LM with parameter sweep
b, b_r, k, h = 0.75, 2.0, 0.5, 4.0
kg = 1/(1-b); beta_r = b_r*kg

A = [1 beta_r; k -h]
rhs = [kg*400.0, -500.0]
Y_star, i_star = A \ rhs
println("Y* = $(round(Y_star,digits=2)), i* = $(round(i_star,digits=4))")

# Sensitivity: multiplier grid over h and b_r
h_vals = [0.5, 1, 2, 4, 8, 100]
br_vals = [0.5, 1, 2, 4]
println("\nFiscal multiplier grid (rows=h, cols=b_r):")
grid = [h_v / (h_v*(1-b) + br_v*k) for h_v in h_vals, br_v in br_vals]
display(round.(grid, digits=3))
```

```r
# R — IS-LM with symbolic Cramer's rule via Ryacas
library(Ryacas)

# Define symbolic variables
b_s <- ysym("b"); br_s <- ysym("b_r"); k_s <- ysym("k"); h_s <- ysym("h")
Ag_s <- ysym("A_bar"); MP_s <- ysym("M_P")

# IS-LM matrix (symbolic)
A_sym <- ysym("{{1, b_r/(1-b)}, {k, -h}}")
b_sym <- ysym("{A_bar/(1-b), -M_P}")
sol   <- yac_str(paste0("Inverse(", as.character(A_sym), ") . ", as.character(b_sym)))
cat("Symbolic solution:\n", sol, "\n")

# Numeric
b <- 0.75; b_r <- 2; k <- 0.5; h <- 4
A_num <- matrix(c(1, k, b_r/(1-b), -h), 2, 2)
b_num <- c(400/(1-b), -500)
solution <- solve(A_num, b_num)
cat(sprintf("Y* = %.2f, i* = %.4f\n", solution[1], solution[2]))
```

---

## 6.7 Programming Exercises

### Exercise 6.1 (APL — Multiplier Sensitivity)

Write a dfn `islm_grid ← {b b_r k h ← ⍵ ⋄ ...}` that returns a $4\times1$ vector of multipliers $(\mu_G, \mu_T, \mu_M, \partial i^*/\partial G)$ for parameter vector input. Then generate the full $10\times10$ fiscal multiplier grid over $(h, b_r) \in \{0.5, 1, 2, 4, 8, 16, 32, 64, 128, \infty\} \times \{0.1, 0.5, 1, 2, 4, 8, 16, 32, 64, 128\}$ using `∘.f` outer product syntax. What is the fiscal multiplier when both $h \to\infty$ and $b_r \to 0$?

### Exercise 6.2 (Python — Mundell–Fleming 3×3)

Implement the full Mundell–Fleming system as a $3\times3$ linear system. For the calibration in Section 6.6, compute equilibrium $(Y^*, i^*, e^*)$ under both fixed and flexible exchange rates. Verify the theoretical prediction that $\mu_G^{flexible} \approx 0$ and $\mu_G^{fixed} \approx \kappa^{open}$ as $\kappa_i \to\infty$.

### Exercise 6.3 (Julia — Liquidity Trap Limit)

Demonstrate the liquidity trap numerically: compute $\mu_G$ and $\mu_M$ as functions of $h$ over $h \in [10^{-3}, 10^3]$ on a log scale. Plot both multipliers and verify: as $h \to\infty$, $\mu_G \to \kappa_G = 1/(1-b)$ and $\mu_M \to 0$. As $h \to 0$, $\mu_G \to 0$ and $\mu_M \to b_r/[b_r k] = 1/k$. Label the liquidity trap and the classical case on the plot.

### Exercise 6.4 — General Comparative Statics ($\star$)

Using the IFT approach from Section 6.4, derive the effect of a change in expected inflation $\pi^e$ on equilibrium output $Y^*$ and the interest rate $i^*$. Show that $\partial Y^*/\partial\pi^e = \mu_G b_r > 0$: higher expected inflation reduces the real interest rate for a given nominal rate, stimulating investment and output. Verify numerically using the calibration of Section 6.6.

### Exercise 6.5 — IS–LM Welfare Analysis ($\star\star$)

Define a social loss function $\mathcal{L} = \lambda_Y(Y - Y^*)^2 + \lambda_\pi(\pi - \pi^*)^2$ where $Y^*$ is potential output and $\pi^*$ the inflation target. In the static IS–LM framework, price level changes from the AD curve translate deviations of $Y$ from $\bar{Y}$ into inflation deviations. Derive the optimal fiscal and monetary policy mix — the $(G^*, M^*)$ combination that minimizes $\mathcal{L}$ subject to the IS–LM equilibrium conditions. Show that the optimal policy targets $Y = \bar{Y}$ when $\lambda_\pi = 0$, and compare with the AD–AS result in Chapter 8.

---

## 6.8 Chapter Summary

**Key results:**

- The IS–LM model is a $2\times2$ linear system $A\mathbf{y} = \mathbf{b}$ with coefficient matrix $A = \begin{pmatrix}1 & \beta_r \\ k & -h\end{pmatrix}$ and $\det(A) = -(h+\beta_r k) < 0$, guaranteeing a unique solution.
- The **fiscal multiplier** $\mu_G = h/[h(1-b)+b_r k]$ lies between 0 (classical LM) and $\kappa_G = 1/(1-b)$ (liquidity trap), with crowding out being the source of the gap.
- The **tax multiplier** $\mu_T = -bh/[h(1-b)+b_r k] = b\mu_T$ is smaller in absolute value than $\mu_G$: tax cuts are less potent than spending increases.
- The **monetary multiplier** $\mu_M = b_r/[h(1-b)+b_r k]$ vanishes in the liquidity trap ($h\to\infty$) and when investment is interest-insensitive ($b_r \to 0$).
- In the **Mundell–Fleming model**: under fixed exchange rates with perfect capital mobility, $\mu_G \to \kappa^{open}$ (no crowding out); under flexible rates, $\mu_G \to 0$ (complete crowding out via exchange rate appreciation).
- In APL: the entire IS–LM solution is `b_vec ⌹ A`; the multiplier grid is generated via `∘.f` outer product in a single expression.

**Connections forward:** Chapter 7 derives the Keynesian cross multiplier as the $h\to\infty$ limit of the IS–LM fiscal multiplier. Chapter 8 embeds IS–LM inside the AD–AS framework, adding price level determination. Chapter 9 derives the ELB multiplier — the $h\to\infty$ analogue for the dynamic NK model.

---

*Next: Chapter 7 — The Multiplier Effect in Closed and Open Economies*
