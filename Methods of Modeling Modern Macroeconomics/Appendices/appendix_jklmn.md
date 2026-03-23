# Appendix J: Data Sources for Macroeconomic Modeling

---

## J.1 FRED (Federal Reserve Economic Data)

**URL:** fred.stlouisfed.org | **API:** Free registration required.

**Key series for macroeconomic modeling:**

| Series ID | Description | Frequency |
|---|---|---|
| GDPC1 | Real GDP (chained 2017 $) | Quarterly |
| PCEPILFE | Core PCE price index | Monthly |
| CPIAUCSL | CPI, all items | Monthly |
| FEDFUNDS | Federal funds rate | Monthly |
| UNRATE | Unemployment rate | Monthly |
| INDPRO | Industrial production index | Monthly |
| HOUST | Housing starts | Monthly |
| GS10 | 10-year Treasury yield | Monthly |
| BAA | Moody's BAA corporate bond yield | Monthly |
| M2SL | M2 money supply | Monthly |

**Python (fredapi):**
```python
from fredapi import Fred
fred = Fred(api_key='YOUR_KEY')
gdp  = fred.get_series('GDPC1', frequency='q')
```

**R (fredr):**
```r
library(fredr)
fredr_set_key("YOUR_KEY")
gdp <- fredr(series_id="GDPC1", frequency="q")
```

**APL (HttpCommand):**
```apl
url ← 'https://api.stlouisfed.org/fred/series/observations?series_id=GDPC1&api_key=YOUR_KEY&file_type=json'
```

## J.2 Penn World Tables (PWT)

**URL:** ggdc.net/pwt | **Coverage:** 183 countries, 1950–present.

**Key variables for cross-country growth analysis:**

| Variable | Description |
|---|---|
| `rgdpna` | Real GDP (national accounts) |
| `ck` | Capital stock (2017 $ mn) |
| `labsh` | Labor share of income ($\alpha$) |
| `delta` | Depreciation rate |
| `rtfpna` | TFP relative to U.S. |
| `pop` | Population |
| `emp` | Employment |

**Use in Chapter 10:** Mankiw–Romer–Weil convergence regression uses PWT savings rates, population growth, and initial GDP.

## J.3 Other Sources

| Source | URL | Key data |
|---|---|---|
| BIS Statistical Warehouse | stats.bis.org | Credit gaps, global liquidity, banking |
| IMF WEO | imf.org/weo | Fiscal data, current accounts, EM |
| OECD Economic Outlook | stats.oecd.org | OECD quarterly data, forecasts |
| Philadelphia Fed RTDSM | philadelphiafed.org/rtdsm | Real-time GDP vintages |
| BEA Input-Output | bea.gov/resources/tables | 71-sector I-O tables |
| Compustat (via WRDS) | wrds-web.wharton.upenn.edu | Firm-level: Tobin's $q$, investment |

## J.4 Real-Time Data Best Practices

For forecasting exercises (Chapter 39):
- Use **real-time vintages** (data as it was available at each date), not the current vintage, to avoid look-ahead bias.
- The Philadelphia Fed's RTDSM provides real-time GDP and GNP data back to 1965.
- Document data vintage with ISO 8601 dates: `GDPC1_2019Q4_vintage_20200129`.

---

# Appendix K: Glossary of Modeling Terms

Selected terms; approximately 80 key entries.

**Approximate aggregation.** The property (Krusell–Smith) that the cross-sectional distribution of wealth can be well approximated by a small number of moments (usually just the mean) for forecasting aggregate prices. Ch. 32.

**Blanchard–Kahn conditions.** The requirement for a unique bounded solution to a rational expectations model: number of eigenvalues outside the unit circle must equal the number of jump (non-predetermined) variables. Ch. 28.

**Calvo pricing.** A price-setting model in which each firm independently adjusts its price with probability $1-\theta$ each period, generating price stickiness. Ch. 27.

**Certainty equivalence.** The property (of first-order approximations) that the policy function is independent of shock variances. Breaks down at second order. Ch. 29.

**Clearing payment vector.** The unique fixed point of the Eisenberg–Noe operator in financial network models; specifies the actual payments made under default. Ch. 34.

**Companion matrix.** The $np\times np$ matrix that converts a $p$-th order VAR into a first-order system. Eigenvalues = roots of the VAR characteristic polynomial. Ch. 14.

**Costate variable.** The shadow price in an optimal control problem; satisfies the adjoint (costate) equation. Denoted $\mu$ (continuous) or $q$ (discrete). Ch. 11.

**Certainty equivalent consumption.** The level of certain consumption that yields the same utility as the uncertain consumption distribution. Ch. 15.

**Contraction mapping.** A function $T$ satisfying $d(Tx,Ty) \leq \beta d(x,y)$ for $\beta < 1$. The Bellman operator is a contraction with modulus $\beta$. Ch. 15.

**Dynamic inefficiency.** The condition $r^* < n+g$ in growth models, implying the economy has over-accumulated capital; PAYG social security improves welfare. Ch. 16.

**EGM (Endogenous Grid Method).** A fast algorithm for solving household consumption problems by inverting the Euler equation. Ch. 32.

**Effective lower bound (ELB).** The lower bound on the nominal interest rate (approximately 0); monetary policy is constrained when the natural rate is negative. Ch. 40.

**Generalized eigenvalue.** The value $\lambda$ satisfying $A\mathbf{v} = \lambda B\mathbf{v}$; computed via the QZ decomposition. Ch. 28.

**Gini coefficient.** A measure of inequality: $G = 2\text{Cov}(X, F(X))/\mathbb{E}[X] \in [0,1]$; 0 = perfect equality, 1 = maximum inequality. Ch. 33.

**Golden Rule.** The capital stock that maximizes steady-state consumption per effective worker: $f'(k^{GR}) = n+g+\delta$. Ch. 10.

**HJB equation.** The Hamilton–Jacobi–Bellman equation — the continuous-time Bellman equation for value functions in dynamic programming. Ch. 32.

**HANK model.** Heterogeneous Agent New Keynesian model; generalizes the NK model to allow household wealth heterogeneity and incomplete markets. Ch. 32.

**HP filter.** The Hodrick–Prescott filter; minimizes a penalized sum of squared deviations to separate trend and cycle. Ch. 17.

**Identification.** A model is identified if distinct parameter values generate distinct observable implications. Local identification: Jacobian of moments w.r.t. parameters has full column rank. Ch. 41.

**Impulse response function (IRF).** The dynamic response of model variables to a structural shock. Computed as $H F^h D$ shock vector. Ch. 19.

**Inada conditions.** Conditions on the production function ensuring interior steady states: $f(0)=0$, $f'(0)=+\infty$, $f'(\infty)=0$. Ch. 10.

**Integrated Assessment Model (IAM).** A model combining an economic model with a climate model to determine optimal emissions trajectories and the social cost of carbon. Ch. 35.

**Kalman gain.** The matrix $K_t = P_{t|t-1}H'F_t^{-1}$ that optimally combines prediction and observation in the Kalman filter. Ch. 20.

**KFP equation.** The Kolmogorov–Fokker–Planck equation; governs the evolution of the density of a diffusion process; the dual of the HJB equation. Ch. 32.

**Log-linearization.** Approximating a nonlinear model by a linear one in log-deviations from the steady state $\hat{x}_t = \ln x_t - \ln x^*$. Ch. 27.

**MSV solution.** Minimum State Variable solution; the unique rational expectations solution that depends only on the minimum set of state variables. Ch. 18.

**Natural rate of interest.** The real interest rate consistent with output at its flexible-price level and stable inflation; denoted $r^n_t$. Ch. 27.

**Perron root.** The unique largest eigenvalue of a positive matrix; corresponds to an all-positive eigenvector (Perron–Frobenius theorem). App. B.

**Pontryagin's maximum principle.** Necessary conditions for optimal control: FOC, costate equation, TVC. Ch. 11.

**Pruning.** A technique for higher-order DSGE approximations that separates first- and second-order state components to prevent explosive paths. Ch. 29.

**QZ decomposition.** The generalized Schur form $Q\Gamma_0 Z = S$, $Q\Gamma_1 Z = T$; provides generalized eigenvalues as $T_{ii}/S_{ii}$. Used in gensys. Ch. 28.

**Rational expectations.** The assumption that agents form expectations consistent with the true model: $\mathbb{E}_t[x_{t+h}]$ equals the model's conditional expectation. Ch. 18.

**Saddle path.** The unique stable manifold in phase space along which the optimal trajectory approaches the steady state of a two-dimensional system. Ch. 11.

**Sims canonical form.** The standardized DSGE system $\Gamma_0\mathbf{y}_t = \Gamma_1\mathbf{y}_{t-1} + \Psi\mathbf{z}_t + \Pi\bm\eta_t$. Ch. 28.

**Social cost of carbon (SCC).** The negative of the shadow price of atmospheric CO₂ in the DICE model; the optimal carbon tax equals the SCC. Ch. 35.

**Sobol index.** A global sensitivity measure: $S_i = \text{Var}_{\theta_i}(\mathbb{E}[Y|\theta_i])/\text{Var}(Y)$; the fraction of output variance attributable to parameter $i$. Ch. 41.

**Sunspot equilibrium.** A rational expectations equilibrium in which the endogenous variable responds to extrinsic randomness; occurs when the Blanchard–Kahn conditions are not satisfied (too few unstable eigenvalues). Ch. 18.

**Taylor principle.** The condition $\phi_\pi > 1$ (more precisely $\phi_\pi + (1-\beta)\phi_y/\kappa > 1$) ensuring determinacy in the NK model. Ch. 28.

**Tauchen method.** A discretization of an AR(1) process into a Markov chain using the normal CDF and conditional probabilities. Ch. 26.

**Transversality condition (TVC).** The terminal condition $\lim_{t\to\infty}e^{-\rho t}\mu(t)x(t) = 0$ ruling out Ponzi schemes in infinite-horizon optimization. Ch. 11.

**URE (Unhedged Interest Rate Exposure).** Auclert's (2019) measure of each household's net exposure to interest rate changes; drives the redistribution channel of monetary policy. Ch. 33.

**Value function.** $V(s)$ = maximum discounted lifetime utility achievable from state $s$; satisfies the Bellman equation $V(s) = \max_c\{u(c,s) + \beta\mathbb{E}[V(s')]\}$. Ch. 15.

**Wold decomposition.** Any zero-mean stationary process has an MA($\infty$) representation: $y_t = \sum_{j=0}^\infty\psi_j\varepsilon_{t-j}$. App. G.

---

# Appendix L: Bibliography and Further Reading

## L.1 Mathematical Methods

**Chiang, A.C. (1992).** *Elements of Dynamic Optimization.* McGraw-Hill. [Optimal control theory; comprehensive and accessible; full proofs of Pontryagin conditions.]

**Stokey, N.L., Lucas, R.E., and Prescott, E.C. (1989).** *Recursive Methods in Economic Dynamics.* Harvard University Press. [The standard reference for dynamic programming in economics; full proofs of the contraction mapping theorem and its applications.]

**Stachurski, J. (2022).** *Economic Dynamics: Theory and Computation.* MIT Press (2nd ed.). [Modern treatment; Python-oriented; strong coverage of Markov chains and stability.]

## L.2 DSGE Models

**Woodford, M. (2003).** *Interest and Prices: Foundations of a Theory of Monetary Policy.* Princeton University Press. [The definitive treatment of the New Keynesian model; derivation of the welfare loss function; optimal commitment policy.]

**Gali, J. (2015).** *Monetary Policy, Inflation, and the Business Cycle.* Princeton University Press (2nd ed.). [Accessible NK model exposition; connects theory to data; covers HANK extensions.]

**Smets, F. and Wouters, R. (2007).** "Shocks and Frictions in US Business Cycles: A Bayesian DSGE Approach." *American Economic Review* 97(3): 586–606. [The canonical medium-scale DSGE; 41 parameters; estimated on U.S. data.]

**Sims, C.A. (2001).** "Solving Linear Rational Expectations Models." *Computational Economics* 20: 1–20. [The gensys algorithm; QZ decomposition for DSGE.]

## L.3 Heterogeneous-Agent Models

**Aiyagari, S.R. (1994).** "Uninsured Idiosyncratic Risk and Aggregate Saving." *Quarterly Journal of Economics* 109(3): 659–684.

**Krusell, P. and Smith, A.A. (1998).** "Income and Wealth Heterogeneity in the Macroeconomy." *Journal of Political Economy* 106(5): 867–896.

**Auclert, A. (2019).** "Monetary Policy and the Redistribution Channel." *American Economic Review* 109(6): 2333–2367.

**Carroll, C.D. (2006).** "The Method of Endogenous Gridpoints for Solving Dynamic Stochastic Optimization Problems." *Economics Letters* 91(3): 312–320.

## L.4 Bayesian Methods

**An, S. and Schorfheide, F. (2007).** "Bayesian Analysis of DSGE Models." *Econometric Reviews* 26(2-4): 113–172.

**Geweke, J. (2005).** *Contemporary Bayesian Econometrics and Statistics.* Wiley.

## L.5 Numerical Methods

**Judd, K. (1998).** *Numerical Methods in Economics.* MIT Press. [The standard reference; covers VFI, perturbation, projection methods, quadrature.]

**Miranda, M.J. and Fackler, P.L. (2002).** *Applied Computational Economics and Finance.* MIT Press. [Practical focus; MATLAB code (adaptable to Python/Julia).]

## L.6 Time Series

**Hamilton, J.D. (1994).** *Time Series Analysis.* Princeton University Press. [The definitive time series textbook; VAR, unit roots, cointegration, Kalman filter.]

**Lutkepohl, H. (2005).** *New Introduction to Multiple Time Series Analysis.* Springer.

---

# Appendix M: Exercise Solutions and Hints

## Selected Solutions

**Chapter 10, Exercise 10.1 (Solow Phase Diagram).**

The steady state satisfies $sf(k^*)^\alpha = \mu k^*$, giving $k^* = (s/\mu)^{1/(1-\alpha)}$. The convergence rate is $\lambda = (1-\alpha)\mu$. In APL:
```apl
kstar ← (s÷mu)*÷1-alpha
lambda ← (1-alpha)×mu
halflife ← (⍟2)÷lambda
```

**Chapter 14, Exercise 14.2 (Analytical Period).**

For the multiplier-accelerator with $b=0.75$, $v=1.0$: $r = \sqrt{bv} = \sqrt{0.75} \approx 0.866$, stable. Discriminant: $b^2(1+v)^2 = 0.5625 \times 4 = 2.25 > 4bv = 3$ → real roots, no oscillation for $v=1$.

For $v=2$: $r = \sqrt{1.5} > 1$ → unstable. For $v = 0.5, 1.5$: $r = 0.612, 1.06$, transition at $bv=1$, i.e., $v=1/b=4/3\approx1.33$.

**Chapter 15, Exercise 15.4 (EGM Hint).**

The EGM for CRRA utility: start with next-period asset grid $\{a'_j\}$; compute $c' = $ policy function at $a'_j$; compute $c = [\beta(1+r)c'^{-\sigma}]^{-1/\sigma}$; compute $a = c + a'/(1+r) - w$. Then interpolate $(a,c)$ back onto the original grid.

**Chapter 22, Exercise 22.2 (Newton Convergence).**

For $f(x) = x^3-2$, $f'(x) = 3x^2$, $f''(x) = 6x$. The convergence constant is $C = |f''(x^*)|/(2|f'(x^*)|) = 6\cdot 2^{1/3}/(2\cdot3\cdot2^{2/3}) = 1/2^{1/3} \approx 0.794$.

Starting from $x_0=2$: $|e_0| \approx 0.26$, $|e_1| \approx 0.053$, $|e_2| \approx 0.0022$. Ratio $|e_1|/|e_0|^2 \approx 0.053/0.068 \approx 0.78 \approx C$.

## Hints for Starred Exercises

**Ch. 11, Ex. 11.5 (Stochastic Euler, ★).** Log-linearize $c_t^{-\sigma} = \beta\mathbb{E}_t[(1+r_{t+1})c_{t+1}^{-\sigma}]$ around $c^* = C^*e^{\hat c_t}$. Use $e^{-\sigma\hat c} \approx 1 - \sigma\hat c$ and $\mathbb{E}_t[e^{-\sigma\hat c_{t+1}}] \approx 1-\sigma\mathbb{E}_t[\hat c_{t+1}]$. The log-linearized result: $\hat c_t = \mathbb{E}_t[\hat c_{t+1}] - (1/\sigma)\mathbb{E}_t[\hat r_{t+1}]$.

**Ch. 32, Ex. 32.4 (HJB, ★).** The upwind scheme: if $\dot a^* > 0$ (saving), use forward difference $V'_a \approx (V_{j+1}-V_j)/\Delta a$; if $\dot a^* < 0$ (borrowing against limit), use backward difference. This preserves the direction of information flow in the HJB.

**Ch. 35, Ex. 35.4 (Fat Tails, ★).** With $P(damage=90\%)=0.01$ and CRRA $\sigma=2$, the expected utility includes $0.01 \times u(0.10Y)^{-\sigma} = 0.01 \times (0.10)^{-2}u(Y) = u(Y)$. The marginal utility contribution diverges as $\sigma\to\infty$ (tail risk becomes infinite with infinite risk aversion), consistent with Weitzman's Dismal Theorem.

---

# Appendix N: Code Repository Guide and Version Control

---

## N.1 Repository Structure

The companion code repository for this book is organized as follows:

```
macroeconomics-methods/
│
├── README.md           # Book overview and quick-start
├── LICENSE             # MIT License
├── CITATION.cff        # How to cite the code
│
├── apl/                # Dyalog APL workspaces and source
│   ├── util.apl        # Shared utilities (RK4, bisection, VFI, MH)
│   ├── part1/          # Ch. 1–5: optimization, linear algebra
│   ├── part3/          # Ch. 10–13: Solow, RCK, OLG, q-model
│   ├── part7/          # Ch. 27–31: DSGE pipeline
│   └── ...
│
├── python/             # Python scripts and notebooks
│   ├── requirements.txt
│   ├── notebooks/      # Jupyter notebooks matching each chapter
│   └── src/
│       ├── dsge/       # DSGE model classes
│       ├── estimation/ # MH, Kalman filter
│       └── utils/      # VFI, quadrature, Tauchen
│
├── julia/              # Julia scripts
│   ├── Project.toml    # Package dependencies
│   ├── Manifest.toml
│   └── src/
│       └── ...
│
├── r/                  # R scripts
│   └── ...
│
└── data/               # Small reference datasets
    ├── us_quarterly_sample.csv
    └── README_data.md
```

## N.2 Git Workflow for Reproducible Research

**Initial setup:**
```bash
git init macroeconomics-methods
cd macroeconomics-methods
git checkout -b main
```

**Branching convention:**
- `main`: stable, paper-ready code.
- `dev`: development; merged into `main` after testing.
- `feature/chapter-32-hank`: feature branches for specific chapters.

**Commit message convention:**
```
[Ch.28] Add gensys QZ decomposition — verified against Dynare output
[Appendix I] Add Julia setup guide
[Bugfix] Fix sign error in net worth accumulation eq. Ch.37
```

**Tagging for paper submissions:**
```bash
git tag -a v1.0 -m "Code version as submitted to JME, March 2026"
git push origin v1.0
```

## N.3 File Naming Conventions

| File type | Convention | Example |
|---|---|---|
| APL source | `snake_case.apl` | `kalman_filter.apl` |
| Python script | `snake_case.py` | `bvar_minnesota.py` |
| Jupyter notebook | `chNN_topic.ipynb` | `ch28_gensys.ipynb` |
| Julia script | `snake_case.jl` | `rbc_vfi.jl` |
| Data file | `SERIES_VINTAGE.csv` | `GDPC1_2024Q4.csv` |
| Dynare .mod | `model_name.mod` | `capstone_nk.mod` |

## N.4 Reproducibility Checklist

Before sharing or submitting code:

- [ ] All random seeds fixed (`np.random.seed(42)`, `Random.seed!(42)`, `set.seed(42)`)
- [ ] Environment documented (`requirements.txt`, `Project.toml`, `sessionInfo()`)
- [ ] Data sources documented with download date and vintage
- [ ] All results reproducible from a single entry script (`run_all.sh` or `main.py`)
- [ ] Unit tests for key functions (steady-state residuals, Blanchard–Kahn check)
- [ ] README with system requirements and estimated runtime

## N.5 Licensing and Citation

All code in the companion repository is released under the **MIT License**: free to use, modify, and distribute with attribution.

**Recommended citation for code:**
```
Methods of Modeling Modern Macroeconomics: Companion Code Repository.
[Authors]. (2026). GitHub: https://github.com/[org]/macroeconomics-methods.
MIT License. Version 1.0.
```

---

*End of Appendices. Thank you for working through this book.*
