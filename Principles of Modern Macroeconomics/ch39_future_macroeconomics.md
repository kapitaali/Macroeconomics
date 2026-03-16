# Chapter 39 — The Future of Macroeconomics: New Theories and Methods

---

## 39.1 HANK Models

Heterogeneous-agent New Keynesian (HANK) models (Kaplan, Moll, and Violante, 2018) combine the Aiyagari–Bewley incomplete-markets problem with New Keynesian nominal rigidities. The key result: because high-MPC households are constrained, the aggregate consumption response to monetary policy depends primarily on the **indirect effect** (via employment and wages) rather than on the direct intertemporal substitution channel. This reverses the RANK conclusion.

The state of the economy must track the distribution of household wealth and income:

$$\hat{C}_t = \int \hat{c}(h,\, \hat{a}_t^h,\, \hat{y}_t^h)\,\mathrm{d}h.$$

Computationally tractable with Reiter's method, Winberry's perturbation, or Achdou et al.'s continuous-time approach.

---

## 39.2 Machine Learning and Causal Inference

Neural networks can represent value functions and policy functions of high-dimensional DSGE models more accurately than polynomial approximations (Duarte, 2018). The double-machine-learning estimator (Chernozhukov et al., 2018) allows high-dimensional nuisance controls in structural equation estimation, improving identification of policy multipliers in the presence of many potential confounders.

---

## 39.3 Rare Disasters and Fat Tails

Barro (2006): consumption disasters of 15% or more occurred in ~3.5% of country-years historically. The disaster risk model:

$$\ln c_{t+1} - \ln c_t = \mu_c + \epsilon_{t+1}^c + v_{t+1}^c\cdot\mathbf{1}\{J_{t+1}=1\}, \quad J_{t+1}\sim\mathrm{Bernoulli}(p),$$

with $p \approx 0.035$. Even low-probability disasters substantially raise the equity premium and resolve the Mehra–Prescott puzzle without requiring extreme risk aversion.

---

## 39.4 Intermediary Asset Pricing

He and Krishnamurthy (2013): financial intermediaries' balance sheets serve as the pricing kernel. When intermediaries are distressed, risk premia rise and investment contracts — connecting microstructure of financial markets to macroeconomic fluctuations. This macro-finance interface is the most active research frontier as of writing.

---

*Next: Chapter 40 — Case Study: The Great Recession of 2008*
