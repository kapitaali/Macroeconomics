# Chapter 13 — Labor Supply and Demand: Wages, Employment, and Unemployment

> *"The labor market is not like the market for widgets. It is a market in which human beings sell their time, and the terms of that sale shape not only their incomes but the structure of their days, their health, and their sense of dignity."*
> — Arthur Okun, *Prices and Quantities*, 1981

---

Labor markets sit at the intersection of virtually every major macroeconomic question. Business cycles manifest primarily as fluctuations in employment and hours; long-run growth translates into living standards through wages; inequality is shaped by the differential between high-skill and low-skill compensation; monetary policy reaches the real economy largely through its effects on the demand for labor. And yet the labor market resists the frictionless, price-clearing description that serves adequately for commodity markets. Workers are not homogeneous; jobs are not standardized; finding a suitable match takes time on both sides; and the wage is determined not only by supply and demand but by bargaining, norms, information asymmetries, and institutional constraints. This chapter develops a sequence of models — each capturing a different dimension of labor market behavior — that together provide the analytical foundation for the macroeconomic treatment of unemployment and wages in subsequent chapters.

---

## 13.1 Labor Supply: Intertemporal Substitution

The most important feature of labor supply for macroeconomics is not the static level of hours worked but how hours respond to temporary versus permanent changes in wages. In business cycle models, what matters is whether workers substitute labor across time — working harder in periods when wages are temporarily high and taking more leisure when wages are temporarily low. This intertemporal margin is what the real business cycle model relies on to generate employment fluctuations from technology shocks.

With separable utility over consumption and leisure, $U = \sum_t \beta^t [u(c_t) + v(1-n_t)]$, where $n_t \in [0,1]$ is the fraction of time spent working and $1-n_t$ is leisure, the intratemporal optimality condition equates the marginal rate of substitution between leisure and consumption to the real wage:

$$\frac{v'(1-n_t)}{u'(c_t)} = w_t.$$

A higher real wage raises the opportunity cost of leisure, inducing more work — this is the **substitution effect**. A higher wage also raises lifetime wealth, inducing more leisure — this is the **income effect**. For macroeconomic purposes, the crucial elasticity is not the Marshallian one (which conflates these two effects) but the Frisch elasticity, which holds marginal utility of wealth constant and thus isolates the pure substitution response.

**Definition (Frisch Elasticity of Labor Supply).** The **Frisch elasticity** $\varepsilon_F$ measures the percentage increase in hours worked in response to a one-percent temporary increase in the real wage, holding the marginal utility of wealth $\lambda$ constant:

$$\varepsilon_F = \frac{\partial \ln n_t}{\partial \ln w_t}\bigg|_{\lambda \text{ fixed}}.$$

With isoelastic disutility $v(1-n) = -v_0 n^{1+\eta}/(1+\eta)$, the Frisch elasticity equals $1/\eta$. With the alternative specification $v(1-n) = v_0(1-n)^{1-\eta}/(1-\eta)$, it equals $(1-n)/(\eta n)$ — decreasing in hours worked and in the curvature parameter $\eta$.

The calibration tension is severe. RBC models require $\varepsilon_F \geq 2$ to generate employment fluctuations comparable in magnitude to the data. With a smaller elasticity, a technology shock causes large wage movements but small hours movements — the model produces smooth employment and volatile real wages, the opposite of what is observed. But microeconometric estimates using household data consistently yield $\varepsilon_F \approx 0.1$–$0.3$ for prime-age males (Chetty et al., 2012), more than an order of magnitude below the macro requirement. This **labor supply puzzle** is one of the deepest unresolved tensions in quantitative macroeconomics.

Several resolutions have been proposed. **Indivisible labor** (Hansen, 1985): if the labor choice is binary — work full-time or not at all — the aggregate Frisch elasticity can be large even when individual elasticities are small, because aggregate hours adjust through the employment margin rather than hours-per-worker. **Extensive margin heterogeneity** (Rogerson and Wallenius, 2009): households with different productivities enter and exit the labor force in response to wage changes, generating a large aggregate elasticity from small individual margins. **Home production** (Benhabib, Rogerson, and Wright, 1991): some household time is devoted to productive non-market activities, changing the tradeoff between market work and leisure in ways that raise the effective aggregate elasticity.

---

## 13.2 Labor Demand

A competitive firm hires labor and capital to maximize profit $\Pi = F(K,N) - wN - rK$, taking factor prices as given. The first-order condition for labor:

$$F_N(K, N) = w.$$

The firm hires workers until the marginal product of labor equals the real wage. Since $F_{NN} < 0$ by diminishing marginal returns, this defines a downward-sloping labor demand schedule: at a higher real wage, fewer workers are hired.

With Cobb-Douglas technology $F(K,N) = AK^\alpha N^{1-\alpha}$, the marginal product of labor is $(1-\alpha)Y/N = (1-\alpha)AK^\alpha N^{-\alpha}$. Setting equal to $w$ and solving: $N^D = [(1-\alpha)AK^\alpha/w]^{1/\alpha}$. The labor demand elasticity with respect to the real wage is $-1/\alpha$, which exceeds one in absolute value when $\alpha < 1/2$, consistent with typical empirical estimates placing capital's share around $0.30$–$0.35$.

The **labor demand curve shifts** with any factor affecting the marginal product: an improvement in TFP $A$ raises MPN at every employment level; capital accumulation raises it further under capital-labor complementarity; and changes in product demand in imperfectly competitive markets affect the firm's effective demand for labor through the output price. In the medium-run framework developed in Chapter 19, the **price-setting curve** $W/P = F_N/(1+\mu)$ generalizes the competitive condition, where $\mu > 0$ is the product market markup. The wedge between the competitive and monopolistic labor demand conditions drives a gap between the marginal product and the real wage that generates equilibrium unemployment without any labor market frictions.

---

## 13.3 The Beveridge Curve and Labor Market Flows

Before turning to the micro-founded matching model, examining aggregate labor market flows is instructive. Workers move between three states: employed (E), unemployed (U), and out of the labor force (N). The steady-state unemployment rate is governed by two key transition rates: the **job-separation rate** $\delta$ (the fraction of employed workers who lose their jobs per period) and the **job-finding rate** $f$ (the fraction of unemployed workers who find jobs per period). Setting flows into and out of unemployment equal:

$$\delta E = f U \implies u^* = \frac{\delta}{\delta + f}.$$

Cyclical movements in unemployment can arise from changes in $\delta$, in $f$, or in both. U.S. evidence (Shimer, 2005) finds that fluctuations in the job-finding rate account for roughly three-quarters of cyclical unemployment variability, with separation rate variation more important during severe contractions like the Great Recession.

**Definition (Beveridge Curve).** The **Beveridge curve** is the empirically observed negative relationship between the unemployment rate and the vacancy rate over the business cycle. In expansions, vacancies are high and unemployment is low; in recessions, vacancies fall and unemployment rises. The curve reflects the underlying matching technology: each combination of vacancies and unemployment generates a given matching flow, and the curve traces out different levels of labor market tightness $\theta = V/U$.

Shifts of the Beveridge curve — outward movements that simultaneously raise both unemployment and vacancies — signal deterioration in matching efficiency: structural changes such as skills mismatch, geographic immobility, or barriers to labor reallocation that reduce the number of matches generated by a given level of search activity. The U.S. Beveridge curve shifted markedly outward following the 2007–09 recession and did not fully return to its pre-crisis position for several years, fueling debate about whether the financial crisis had caused lasting structural damage to the labor market's matching technology.

---

## 13.4 Search, Matching, and the Natural Rate

The Beveridge curve motivates a formal model of equilibrium unemployment based on costly search and matching. The **Mortensen–Pissarides model** (Mortensen, 1970; Pissarides, 1985) is the standard framework.

At the core is a **matching function** $m(U_t, V_t)$ that gives the flow of new employment relationships as a function of unemployed workers and open vacancies. Assuming constant returns to scale and the Cobb-Douglas form:

$$m(U, V) = s U^\alpha V^{1-\alpha}, \quad \alpha \in (0,1).$$

**Definition (Labor Market Tightness).** **Labor market tightness** $\theta \equiv V/U$ is the vacancy-to-unemployment ratio. The job-finding rate for workers is $f(\theta) = m/U = s\theta^{1-\alpha}$, increasing in $\theta$: more vacancies per unemployed worker means faster job-finding. The vacancy-filling rate for firms is $q(\theta) = m/V = s\theta^{-\alpha}$, decreasing in $\theta$: a tight market means each vacancy takes longer to fill.

The labor market equilibrium is determined by two conditions. The **job creation condition** requires that firms post vacancies until the value of a filled job equals the expected cost of posting: $J = c/q(\theta)$, where $c$ is the flow cost of maintaining a vacancy and $J$ is the present value of a filled job to the firm. This condition pins down equilibrium tightness $\theta^*$. Given $\theta^*$, the steady-state natural unemployment rate is:

$$u^* = \frac{\delta}{\delta + f(\theta^*)}.$$

Higher job destruction $\delta$ raises the natural rate; more efficient matching (higher $s$, which raises $f$) reduces it. The model provides microfoundations for the NAIRU that Chapter 3 introduced as a measurement concept and that Chapter 10 embedded in the Phillips curve.

---

## 13.5 Wage Determination: Nash Bargaining and the Hosios Condition

In the matching model, each matched worker-firm pair bargains bilaterally over the wage. The standard solution is the **Nash bargaining wage**, which maximizes the weighted product of the worker's and firm's surplus from the match:

$$w^* = \underset{w}{\arg\max}\;(W - U^{worker})^\eta \cdot (J - V)^{1-\eta},$$

where $W$ is the value of employment at wage $w$, $U^{worker}$ is the value of unemployment (the worker's outside option), $J$ is the value of a filled job to the firm, $V \approx 0$ is the value of a vacancy under free entry, and $\eta \in [0,1]$ is the worker's **bargaining power**.

Solving the Nash program and using the equilibrium value functions yields:

$$w^* = (1-\eta)b + \eta\big[y + c\theta\big],$$

where $b$ is the worker's flow utility when unemployed and $y$ is labor productivity. The term $c\theta$ reflects the labor market tightness premium: in a tight market ($\theta$ high), workers' outside options are better because the unemployed find jobs quickly, pushing up the bargained wage. Wages are procyclical — rising in booms as $\theta$ increases and falling in recessions as it falls — consistent with the empirical evidence.

**Definition (Hosios Condition).** The **Hosios condition** (Hosios, 1990) states that the decentralized matching equilibrium is socially efficient if and only if the worker's bargaining power equals the matching function's elasticity with respect to unemployment:

$$\eta = \alpha.$$

When $\eta > \alpha$, workers capture more than the socially optimal share of the match surplus, depressing firms' incentives to post vacancies, reducing tightness below its efficient level, and raising the unemployment rate above the social optimum. When $\eta < \alpha$, firms post excess vacancies. The Hosios condition is rarely satisfied in practice, providing a theoretical basis for policies such as hiring subsidies (which correct for $\eta < \alpha$) or unemployment benefits (which partially substitute for $\eta > \alpha$ scenarios by improving workers' outside options).

---

## 13.6 Efficiency Wages and Involuntary Unemployment

Frictional unemployment — workers between jobs while searching — is voluntary in a meaningful sense: workers have not yet found a match, but they are not excluded from employment at the prevailing wage. The deeper question is whether there exists **involuntary unemployment**: workers who would willingly accept a job at the current wage but cannot find one, not because of search frictions but because of deliberate exclusion.

The efficiency wage theory of Shapiro and Stiglitz (1984) generates this outcome. If employers cannot perfectly monitor worker effort, workers have an incentive to shirk — to supply minimum effort while collecting the full wage. Rational firms anticipate this and choose wages above the market-clearing level to make shirking unattractive.

The key mechanism: if a worker is caught shirking, she is fired and must search for a new job while unemployed. The value of a job therefore exceeds the value of unemployment by the entire present value of the employment premium. When unemployment is low, the expected duration of unemployment after firing is short — the cost of being caught is low — so firms must pay a larger premium to deter shirking. When unemployment is high, the threat of prolonged joblessness disciplines workers at a smaller wage premium.

The **no-shirking condition (NSC)** specifies the minimum wage at which a worker will not shirk:

$$w \geq b + e_H + \frac{e_H(r + q_f)}{q_f} \cdot \frac{1-u}{u} \cdot \frac{1}{r+\delta+f(u)},$$

where $b$ is the flow utility of unemployment, $e_H$ is the effort cost of working diligently, $q_f$ is the probability of being caught shirking per period, $r$ is the discount rate, and $u$ is the unemployment rate. The NSC is upward-sloping in $(u, w)$ space: higher unemployment tightens the discipline, allowing firms to deter shirking at a lower wage. The equilibrium is the intersection of the NSC with the profit-maximizing labor demand curve $w = F_N(K,N^{EW})$, which uniquely determines an efficiency wage $w^{EW}$ above the market-clearing level and an employment level $N^{EW}$ below full employment.

The workers who are unemployed in this equilibrium — $L - N^{EW}$ workers at the margin — would willingly accept a job at $w^{EW}$ if one were offered. They are involuntarily unemployed: excluded not by search frictions but by the deliberate design of firms' wage policy. This distinction matters for policy: efficiency wage unemployment does not respond to aggregate demand stimulus in the way that cyclical unemployment does, because hiring more workers at the current wage would reduce the unemployment premium and invite shirking. The efficiency wage equilibrium is not a market failure requiring demand management but a rational equilibrium response to moral hazard — though it may nonetheless justify interventions that improve monitoring technology or alter the structure of labor contracts.

---

## 13.7 Structural Unemployment and the Mismatch Problem

Neither frictional nor efficiency wage unemployment captures a third important category: unemployment that persists even when aggregate demand is sufficient to employ everyone, because workers' skills or locations do not correspond to available vacancies.

**Definition (Structural Unemployment).** **Structural unemployment** arises when there is a fundamental mismatch between the characteristics of unemployed workers — their skills, location, or industry experience — and the characteristics demanded by available vacancies. Unlike frictional unemployment, structural unemployment is not a transitional state in the job-matching process but a persistent gap between what workers supply and what the economy demands.

The Beveridge curve provides the natural diagnostic. Demand-driven recessions move the economy along a stable curve — unemployment up, vacancies down. Structural deterioration shifts the curve outward: both unemployment and vacancies rise simultaneously because the vacancies cannot be filled by the workers who are unemployed. The outward shift of the U.S. Beveridge curve after 2009, which persisted long after measured demand had recovered, generated significant debate about whether the financial crisis had induced lasting structural mismatches.

The policy implications diverge sharply from cyclical unemployment. Aggregate demand stimulus — monetary easing, fiscal expansion — can address cyclical unemployment by restoring the flow of job-creating demand. Structural unemployment is unresponsive to such tools: firms will not hire workers who lack the required skills or who live in distant locations, regardless of the level of aggregate demand. Addressing structural unemployment requires supply-side interventions: retraining and reskilling programs, geographic mobility subsidies, educational investment in forward-looking skills, and in some cases place-based policies that bring investment to depressed regions rather than requiring workers to move. The empirical challenge — separating structural from cyclical in real time — is one the hardest tasks in applied macroeconomics and has considerable stakes for policy.

---

## Chapter Summary

- Labor supply is characterized by the **Frisch elasticity** $\varepsilon_F$, measuring hours' response to a temporary wage change at fixed marginal utility of wealth. RBC models require $\varepsilon_F \geq 2$; microeconometric estimates yield $\varepsilon_F \approx 0.1$–$0.3$, a discrepancy partially resolved by extensive margin and home production arguments.

- Competitive labor demand equates $F_N(K,N) = w$; with product market power, the price-setting curve shifts this to $F_N/(1+\mu)$, creating a wedge that generates equilibrium unemployment even without frictions.

- The **steady-state unemployment rate** $u^* = \delta/(\delta + f)$ depends on separation and job-finding rates. The **Beveridge curve** traces vacancy-unemployment combinations; outward shifts signal structural deterioration in matching efficiency.

- The **Mortensen-Pissarides** model formalizes frictional unemployment. Equilibrium tightness $\theta^*$ is determined by free-entry vacancy posting; the **Hosios condition** $\eta = \alpha$ characterizes social efficiency of the equilibrium.

- **Efficiency wages** generate involuntary unemployment as a rational equilibrium: firms pay above the market-clearing wage to deter shirking, and the resulting unemployment pool disciplines workers by making job loss costly.

- **Structural unemployment** reflects skill and geographic mismatch, visible as Beveridge curve shifts, and responds to supply-side interventions rather than demand stimulus.

---

*Next: Chapter 14 — Money Demand and Supply: The Role of Central Banks*
