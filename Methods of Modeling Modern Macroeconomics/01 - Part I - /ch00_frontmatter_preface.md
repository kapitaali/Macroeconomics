# Methods of Modeling Modern Macroeconomics
## A Mathematical and Computational Approach

**Companion to:** *Principles of Modern Macroeconomics: From Theory to Everyday Application*

---

*First Edition*

---

> *"The purpose of models is not to fit the data but to sharpen the questions."*
> — Samuel Karlin

---

**Languages:** Dyalog APL · Python · Julia · R

**Software:** Dynare 6 · Dyalog APL 18.2+ · Python 3.11+ · Julia 1.10+ · R 4.3+

---

---

# Preface

This book was written to answer a question that every student of macroeconomics eventually asks, usually at the moment when the intuition finally clicks but the machinery still seems out of reach: *how do you actually compute this?*

The companion volume, *Principles of Modern Macroeconomics: From Theory to Everyday Application*, develops the ideas of modern macroeconomics with care and depth. It defines concepts precisely, derives core results, and traces the intellectual history of the discipline from Say's Law to HANK models. What it does not do — by design — is show you how to solve a Bellman equation numerically, how to implement the Kalman filter for estimating the output gap, how to write a Dynare `.mod` file for a medium-scale New Keynesian model, or how to simulate a DSGE model and compare its second moments to the data. This book does all of those things.

The relationship between the two books is not supplementary but constitutive. Every technique introduced in these pages is motivated by a question raised in *Principles*. Chapter 11 of *Principles* introduces the consumption Euler equation; Chapter 15 of this book shows you how to solve the full Bellman equation via value function iteration. Chapter 27 of *Principles* describes the RBC model and its calibration; Chapter 17 of this book walks through the complete computational solution. Chapter 23 of *Principles* derives optimal monetary policy; Chapter 40 of this book solves for the optimal Taylor rule coefficients numerically and computes the welfare cost of suboptimal policy. The cross-references are not decorative — they are the organizing principle of the entire volume.

### Who This Book Is For

The primary audience is advanced undergraduates and master's students in economics who have worked through *Principles* and now want to go further. A secondary audience is PhD students in the early years of their program who need a bridge between graduate coursework in mathematics and the practical quantitative skills expected in the research literature. A third audience consists of policy economists — at central banks, finance ministries, and international organizations — who use DSGE models or VAR analyses in their work and want to understand what the software is actually doing inside.

The book assumes college-level calculus (derivatives, integrals, partial derivatives, the chain rule, the implicit function theorem), one semester of linear algebra (vectors, matrices, eigenvalues, determinants), and basic probability theory (random variables, expectations, conditional distributions). It does not assume any prior programming experience, though readers who have written code in any language will find the material easier. Every algorithm is introduced from scratch with mathematical justification before code is shown.

### Four Languages

Every computational example in this book is implemented in four languages: **Dyalog APL**, **Python**, **Julia**, and **R**. This is not redundancy — it is a recognition that different computing communities have converged on different tools, and that seeing the same algorithm implemented differently illuminates the algorithm itself.

Dyalog APL is the primary matrix and time-series language throughout. APL's array-oriented design means that matrix operations, time-series transformations, and numerical algorithms that require dozens of lines in other languages can often be expressed in a single APL expression. The value function iteration inner loop, the Kalman filter update step, and the Tauchen discretisation procedure are each five lines or fewer in APL. More importantly, APL's notation maps naturally onto the mathematical objects it manipulates: `+.×` *is* matrix multiplication, `⌹` *is* matrix inversion (or least-squares solution), and `⍣≡` *is* iteration to a fixed point. Reading APL is, in many contexts, reading mathematics.

A crucial feature of APL that beginners must understand immediately: **APL evaluates from right to left**. The expression `f g h x` means $f(g(h(x)))$ — apply `h` first, then `g`, then `f`. This is not a quirk but a design choice that makes composition natural and mirrors the way mathematicians write composed functions. Throughout this book, all APL code uses the dfns (dynamic-function) style — `f←{body of f}` — with `⍺` for the left argument and `⍵` for the right argument. All code sets `⎕IO←0` (zero-based indexing) and `⎕ML←1` (Dyalog migration level 1) at the top of every workspace.

Python is used for its scientific ecosystem: NumPy and SciPy for numerical computation, statsmodels for econometrics, QuantEcon for dynamic programming, and matplotlib for visualization. Julia is used where performance matters at scale — large DSGE models, heterogeneous-agent computations, and simulation-heavy estimation. R is the standard for econometric work: the `vars`, `KFAS`, `gmm`, and `bvars` packages cover VAR estimation, state-space filtering, GMM, and Bayesian VAR respectively.

### On Dynare

Part VII assumes familiarity with the Dynare software platform, the standard tool for solving and estimating DSGE models at central banks worldwide. Dynare is not a programming language but a preprocessor: you write a `.mod` file describing your model, and Dynare's parser translates it into MATLAB or Octave code, solves the model, and generates impulse response functions, simulated moments, and — if requested — a posterior distribution via Metropolis–Hastings. Chapter 31 provides a complete guide to the Dynare workflow, including how to interpret Dynare's output and post-process it in APL, Python, Julia, or R.

### On Reproducibility

All code in this book is available in the companion repository described in Appendix N. The repository is organised by chapter, with one subdirectory per chapter containing the APL workspace, Python script, Julia script, and R script implementing that chapter's worked example. Every script produces exactly the figures and tables shown in the text when run with the specified software versions. Appendix N explains the version control strategy and how to use the repository.

### How This Book Is Structured

Part I builds the mathematical foundations: optimization, linear algebra, differential equations, difference equations, and stochastic processes. Part II solves the static models of *Principles* Part II analytically. Parts III and IV develop dynamic models in continuous and discrete time respectively. Part V introduces stochastic methods — rational expectations algebra, time series, the Kalman filter, and structural estimation. Part VI covers the numerical methods that underlie all modern macroeconomic computation. Part VII assembles the complete DSGE pipeline. Part VIII surveys the frontier: heterogeneous agents, network models, integrated assessment models, and agent-based computation. Part IX applies everything to real-world macroeconomic questions, with detailed replication exercises for the Great Recession and the COVID-19 pandemic.

---

---

# Acknowledgments

The intellectual debts accumulated in writing a comprehensive methods text in quantitative macroeconomics are many and deep. The theoretical frameworks described in these pages rest on foundations laid by Lars Peter Hansen, Thomas Sargent, Christopher Sims, Lawrence Christiano, Martin Eichenbaum, Charles Evans, Frank Smets, Rafael Wouters, Olivier Blanchard, Jordi Galí, Mark Gertler, and many others whose papers are cited throughout. The dynamic programming methods of Part IV owe everything to Richard Bellman and Nancy Stokey, Robert Lucas, and Edward Prescott. The numerical methods of Part VI are standard results from numerical analysis, presented here with macroeconomic motivation.

The Dyalog APL idioms developed throughout this book reflect the design philosophy of the APL language and the Dyalog implementation. The `⍣` power operator, the dfns namespace, and `⎕PY` for Python interoperability are features of Dyalog APL 18.2 and represent decades of language development.

The companion volume *Principles of Modern Macroeconomics* provides the economic content that gives every technique in this book its purpose. Any errors of exposition, omission, or commission in these pages are entirely the author's own.

---

---

# How to Use This Book

### Reading Paths

The book supports three independent reading paths, corresponding to different professional goals. The dependency diagram at the end of this section shows which chapters require which.

**Path A — Mathematical Macro Theory** develops the continuous-time and discrete-time dynamic models that underlie the DSGE literature, with emphasis on analytical results and phase diagrams. It covers Parts I–IV and Chapter 18. Estimated time: 30 hours. Recommended for students preparing for growth and business cycle theory courses.

**Path B — Computational DSGE** is the core track for students who want to build, solve, and estimate DSGE models. It covers Parts I (Chapters 1, 2, 4, 5), Part IV, Part V (Chapter 18), Part VI, Part VII, and Part IX (Chapters 40–42). Estimated time: 40 hours. This path ends with a complete working DSGE model estimated on real data.

**Path C — Empirical and Estimation** focuses on time-series methods, filtering, GMM, Bayesian estimation, and forecasting. It covers Part I (Chapters 1, 5), Part V (all), Chapter 26, Chapter 30, and Part IX (Chapters 39, 41). Estimated time: 30 hours. Recommended for students preparing for empirical macroeconomics research.

### Typographic Conventions

The book uses four named box types:

> **Definition.** A formal mathematical definition. Every technical term is defined before use.

> **Theorem / Proposition.** A formal result with either a full proof or an explicit proof sketch. Starred results ($\star$) have proofs in the appendices.

> **Algorithm.** Numbered pseudocode. Every algorithm is stated in language-neutral pseudocode before any implementation is shown.

> **Worked Example.** A fully worked numerical or analytical example. Never a sketch — always completed to the final numerical answer or closed-form expression.

Code appears in monospace font, tagged by language:

```apl
⍝ Dyalog APL
```
```python
# Python
```
```julia
# Julia
```
```r
# R
```

All APL code assumes `⎕IO←0` (zero-based indexing) and `⎕ML←1`. These two lines appear at the top of every APL workspace in the companion repository and should be set before running any APL example in this book.

### Cross-Reference Format

References to the companion volume use the format **[P:Ch.X]** for a chapter reference and **[P:Ch.X.Y]** for a section reference. For example, **[P:Ch.11.3]** refers to Section 11.3 of *Principles of Modern Macroeconomics* (the Euler equation approach to consumption). References to chapters within this volume use **[M:Ch.X]** format.

### Exercises

Each chapter contains exercises of three difficulty levels. Standard exercises (no marker) can be completed with the methods in that chapter. Starred exercises ($\star$) require combining techniques from multiple chapters. Double-starred exercises ($\star\star$) are research-level and are suitable for term papers or independent study. Hints and selected full solutions appear in Appendix M.

### A Note on Mathematical Notation

This book uses a consistent notation throughout. A summary table appears in the preliminary pages. Where a symbol carries different meanings in different contexts (as $\sigma$ does — it is the CRRA coefficient in Chapter 11 and the IS slope in Chapter 7 of *Principles*), the intended meaning is always stated explicitly at first use. Appendix K contains a complete glossary.

---

---

# Prerequisites

The following is a precise statement of what this book assumes. Readers who are comfortable with all items in each category should have no difficulty with the material.

### Calculus

- **Single-variable differentiation:** derivatives of polynomials, exponential, logarithmic, and trigonometric functions; product rule, quotient rule, chain rule; implicit differentiation
- **Single-variable integration:** antiderivatives; definite integrals; integration by substitution and by parts; improper integrals
- **Multivariate differentiation:** partial derivatives; the gradient vector $\nabla f$; the Hessian matrix $H_f$; directional derivatives; Young's theorem (symmetry of mixed partials)
- **The chain rule in vector form:** if $\mathbf{y} = g(\mathbf{x})$ and $z = f(\mathbf{y})$, then $\partial z/\partial x_i = \sum_j (\partial f/\partial y_j)(\partial g_j/\partial x_i)$
- **Taylor series:** first- and second-order Taylor expansions around a point; remainder term; multivariate Taylor expansion to first order (the linearization formula)
- **The implicit function theorem:** if $F(x, y) = 0$ defines $y$ implicitly as a function of $x$ near a point where $F_y \neq 0$, then $\mathrm{d}y/\mathrm{d}x = -F_x/F_y$

### Linear Algebra

- **Vectors and matrices:** addition, scalar multiplication, transposition
- **Matrix multiplication:** the product $AB$ requires the inner dimensions to match; $(AB)' = B'A'$
- **The determinant:** definition for $2\times 2$ and $3\times 3$ matrices; Sarrus's rule; the cofactor expansion
- **Matrix inversion:** definition; when it exists (nonsingular = full rank); the formula for $2\times 2$ inverses; Cramer's rule for small systems
- **Eigenvalues and eigenvectors:** definition ($Av = \lambda v$); the characteristic polynomial; the trace and determinant as sum and product of eigenvalues
- **Diagonalisation:** when a matrix is diagonalisable; $A = PDP^{-1}$ where $D$ is diagonal; powers $A^n = PD^nP^{-1}$

### Probability and Statistics

- **Random variables:** discrete and continuous; probability mass and density functions; the CDF
- **Expectations:** $\mathbb{E}[X]$, $\mathbb{E}[g(X)]$; linearity; iterated expectations $\mathbb{E}[\mathbb{E}[X|Y]] = \mathbb{E}[X]$
- **Variance and covariance:** definitions; $\mathrm{Var}(aX + b) = a^2\mathrm{Var}(X)$; the covariance matrix
- **Common distributions:** normal $\mathcal{N}(\mu, \sigma^2)$; standard normal $\Phi$; lognormal; uniform; the Bernoulli and binomial; the exponential
- **Conditional distributions:** $f(x|y) = f(x,y)/f(y)$; Bayes' theorem
- **Basic estimation:** sample mean and variance as estimators; the law of large numbers (informal); the central limit theorem (informal)

### Economics

- **Microeconomic optimization:** utility maximization subject to a budget constraint; profit maximization; the first-order conditions; the interpretation of Lagrange multipliers
- **The IS–LM and AS–AD models** at the level of *Principles* Chapters 7–10: graphical analysis, direction of shifts, policy multipliers (qualitative)
- **The Solow model** at the level of *Principles* Chapter 5: the fundamental ODE, the steady state, the Golden Rule (conceptual, not computational)
- **Rational expectations** at the level of *Principles* Chapter 15: the definition; the distinction from adaptive expectations; the forward solution of the NKPC as a present-value formula

### What Is Not Required

- No prior programming experience is assumed. Every algorithm is introduced from scratch with full pseudocode before any code is shown.
- Measure-theoretic probability is not required, though the concepts of a filtration $\{\mathcal{F}_t\}$ and a martingale are introduced carefully when needed.
- Optimal control theory (Hamiltonians, Pontryagin's maximum principle) is not assumed — it is developed in Chapter 11.
- Dynamic programming (Bellman equations, value function iteration) is not assumed — it is developed in Chapter 15.

---

---

# Reading Guide

### Estimated Time per Chapter

| Part | Chapter | Title (abbreviated) | Est. hours |
|---|---|---|---|
| I | 1 | Calculus and optimization | 2.5 |
| I | 2 | Linear algebra | 2.5 |
| I | 3 | Differential equations | 3.0 |
| I | 4 | Difference equations | 2.5 |
| I | 5 | Stochastic processes | 2.5 |
| II | 6 | Solving the IS–LM model | 2.0 |
| II | 7 | The multiplier effect | 2.0 |
| II | 8 | The AD–AS model | 2.0 |
| II | 9 | Tax and spending multipliers | 2.0 |
| III | 10 | Solow: solving the ODE | 2.5 |
| III | 11 | RCK: optimal control | 3.5 |
| III | 12 | Continuous-time OLG | 2.5 |
| III | 13 | Adjustment costs and Tobin's q | 2.5 |
| IV | 14 | Difference equations in macro | 2.5 |
| IV | 15 | Dynamic programming | 3.5 |
| IV | 16 | Discrete-time OLG | 2.5 |
| IV | 17 | Recursive RBC methods | 3.5 |
| V | 18 | Rational expectations | 3.0 |
| V | 19 | Time series methods | 3.5 |
| V | 20 | The Kalman filter | 3.0 |
| V | 21 | GMM and maximum likelihood | 3.0 |
| VI | 22 | Newton–Raphson | 2.0 |
| VI | 23 | Numerical integration | 2.0 |
| VI | 24 | Numerical optimization | 2.0 |
| VI | 25 | Systems of equations | 2.0 |
| VI | 26 | Monte Carlo methods | 2.5 |
| VII | 27 | Log-linearization | 2.5 |
| VII | 28 | Blanchard–Kahn and Sims | 3.5 |
| VII | 29 | Perturbation methods | 3.0 |
| VII | 30 | Bayesian estimation | 3.5 |
| VII | 31 | Dynare and software workflows | 3.0 |
| VIII | 32 | Heterogeneous agents | 4.0 |
| VIII | 33 | Inequality dynamics | 3.0 |
| VIII | 34 | Network models | 3.0 |
| VIII | 35 | Integrated assessment (DICE) | 3.0 |
| VIII | 36 | Agent-based models | 3.0 |
| IX | 37 | Replicating the Great Recession | 4.0 |
| IX | 38 | COVID-19 epidemic-economic models | 3.5 |
| IX | 39 | Forecasting GDP and inflation | 3.5 |
| IX | 40 | Policy analysis with NK model | 3.5 |
| IX | 41 | Model validation and sensitivity | 3.0 |
| IX | 42 | Capstone: DSGE from scratch | 4.5 |
| | | **Total** | **~120 hrs** |

*Note: hours include reading, working through derivations, and completing the programming exercises. The ~70-hour figure in the overview reflects a selective reading on one of the three paths.*

### Chapter Dependency Diagram

```
Part I: Foundations
  Ch.1 (Calculus/Optim)
  Ch.2 (Linear Algebra)
  Ch.3 (ODEs) ──────────────────────────────┐
  Ch.4 (Difference Eqs) ────────────────────┤
  Ch.5 (Stochastic Processes)               │
       │         │         │                │
       ▼         ▼         ▼                │
Part II (Static) Part III (Continuous)  Part IV (Discrete)
  Ch.6   Ch.7    Ch.10  Ch.11  Ch.12  Ch.13  Ch.14  Ch.15  Ch.16  Ch.17
  Ch.8   Ch.9                                
       │                   │                    │
       └───────────────────┴────────────────────┘
                           │
                    Part V (Stochastic)
                  Ch.18  Ch.19  Ch.20  Ch.21
                           │
                    Part VI (Numerical)
                  Ch.22  Ch.23  Ch.24  Ch.25  Ch.26
                           │
                    Part VII (DSGE Pipeline)
                  Ch.27  Ch.28  Ch.29  Ch.30  Ch.31
                     │                          │
              Part VIII (Advanced)      Part IX (Applications)
            Ch.32-36                   Ch.37-42
```

### Three Reading Paths in Detail

**Path A — Mathematical Macro Theory**
Chapters: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 18
Appendices: A, B, C, E, F
Hours: ~30

**Path B — Computational DSGE**
Chapters: 1, 2, 4, 5, 14, 15, 17, 18, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 37, 42
Appendices: H, I, N
Hours: ~40

**Path C — Empirical and Estimation**
Chapters: 1, 5, 18, 19, 20, 21, 22, 24, 26, 30, 39, 41
Appendices: D, G, J, L, M
Hours: ~30
