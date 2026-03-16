# Chapter 13 — Labor Supply and Demand: Wages, Employment, and Unemployment

---

Labor markets are at the intersection of several of the most important macroeconomic phenomena: business cycles, income distribution, technological change, and long-run growth all work primarily through their effects on employment and wages. Understanding how labor supply and demand jointly determine employment and wages — and understanding why labor markets do not always clear as smoothly as the commodity markets of intermediate microeconomics — is essential for the analysis of unemployment, the Phillips curve, and the transmission of monetary policy.

---

## 13.1 Labor Supply: Intertemporal Substitution

The most important feature of labor supply for macroeconomics is not the static level of hours worked but how hours respond to temporary versus permanent changes in wages. In business cycle models, what matters is whether workers substitute labor across time — working harder when wages are temporarily high and taking more leisure when wages are temporarily low.

With separable utility over consumption and leisure, $U = \sum_t \beta^t [u(c_t) + v(1-n_t)]$ where $n_t \in [0,1]$ is hours worked and $1-n_t$ is leisure, the static optimality condition is:

$$\frac{v'(1-n_t)}{u'(c_t)} = w_t,$$

equating the marginal rate of substitution between leisure and consumption to the real wage. Combined with the consumption Euler equation $u'(c_t) = \beta(1+r_{t+1})\mathbb{E}_t[u'(c_{t+1})]$, the intertemporal labor supply condition is derived: work more in periods when wages are high relative to lifetime wealth.

**Definition (Frisch Elasticity of Labor Supply).** The **Frisch elasticity** $\varepsilon_F$ measures the percentage increase in hours worked in response to a 1% temporary increase in wages, holding the marginal utility of wealth constant. It is the relevant elasticity for analyzing business cycle fluctuations, where wage changes are typically transitory. With $v(1-n) = v_0(1-n)^{1-\eta}/(1-\eta)$:

$$\varepsilon_F = \frac{1-n}{\eta\, n}.$$

RBC models require $\varepsilon_F \geq 2$ to match the large employment fluctuations observed in business cycles. Microeconometric estimates using household data typically yield $\varepsilon_F \approx 0.1$–$0.3$ for prime-age males (Chetty et al., 2012), more than an order of magnitude smaller. This **labor supply puzzle** is one of the central unresolved tensions between DSGE macroeconomic models and microeconomic evidence.

---

## 13.2 Labor Demand

A competitive firm chooses labor to maximize profits $\Pi = F(K,N) - wN - rK$. The first-order condition is simply:

$$F_N(K, N) = w.$$

The firm hires workers until the marginal product of labor equals the real wage. Since $F_{NN} < 0$ (diminishing marginal returns), this defines a downward-sloping labor demand curve: at a higher real wage, fewer workers are hired.

---

## 13.3 Search, Matching, and the Natural Rate

Competitive labor market models predict instantaneous wage and employment adjustment. But workers take time to find suitable jobs, and firms take time to fill vacancies. This matching friction gives rise to equilibrium unemployment even when wages are flexible — the **frictional unemployment** component of the natural rate.

The **Mortensen–Pissarides search-and-matching model** (Mortensen, 1970; Pissarides, 1985) provides the standard microeconomic foundation for the natural rate. The matching function $m(U_t, V_t)$ gives the flow of new employment relationships as a function of the stock of unemployed workers $U_t$ and open vacancies $V_t$. Under constant returns to scale: $m = m(U,V)$, with $m_U, m_V > 0$.

**Definition (Labor Market Tightness).** **Labor market tightness** $\theta_t \equiv V_t/U_t$ is the ratio of vacancies to unemployed workers. It measures how easy it is for workers to find jobs (high $\theta$ = many vacancies per unemployed worker = good for job-seekers) and how hard it is for firms to fill vacancies (high $\theta$ = tight market = costly vacancy-filling for firms). The job-finding rate for workers $f(\theta) = m(U,V)/U = m(1, \theta^{-1})^{-1}... = m(\theta^{-1}, 1)\theta$. With Cobb-Douglas matching $m = s U^\alpha V^{1-\alpha}$: $f(\theta) = s\theta^{1-\alpha}$ (increasing in $\theta$) and vacancy-filling rate $q(\theta) = m/V = s\theta^{-\alpha}$ (decreasing in $\theta$).

In steady state, with job-destruction rate $\delta$, the natural rate satisfies:

$$u^* = \frac{\delta}{\delta + f(\theta^*)}.$$

This says the natural rate is determined by the balance between job destruction (the flow into unemployment) and job finding (the flow out). Higher job destruction $\delta$ raises the natural rate; faster matching (higher $f(\theta^*)$) reduces it.

---

## 13.4 Wage Determination: Nash Bargaining

In the matching model, wage bargaining takes place bilaterally between each matched worker-firm pair. The Nash bargaining solution maximizes the product of the worker's and firm's net gains from the match:

$$w^* = \underset{w}{\arg\max}\;(w - b)^\eta\,(J - V)^{1-\eta},$$

where $b$ is the worker's outside option (unemployment benefit plus value of leisure), $J$ is the value of a filled job to the firm, $V$ is the value of a vacancy, and $\eta \in [0,1]$ is the worker's bargaining power — the weight on the worker's surplus in the Nash product. The solution gives:

$$w^* - b = \frac{\eta}{1-\eta}(J - V)\cdot q(\theta).$$

Wages increase with labor market tightness $\theta$ (more vacancies mean workers have better outside options) and with bargaining power $\eta$. The **Hosios condition** for efficient matching is $\eta = \alpha$: when workers' bargaining power equals the matching function's elasticity with respect to unemployment, the private equilibrium is socially optimal and unemployment equals the efficient level.

---

## 13.5 Efficiency Wages and Involuntary Unemployment

Search-and-matching models generate frictional unemployment — people between jobs. But what generates **involuntary unemployment** — workers who are unemployed but would willingly work at the current wage if only they could find a job? One prominent theory is efficiency wages.

Shapiro and Stiglitz (1984) argue that firms may intentionally pay wages above the market-clearing level to deter workers from shirking. If firms cannot perfectly monitor effort, workers will shirk (exert zero effort) whenever the expected cost of being caught is less than the effort cost $e_H$. The **no-shirking condition** requires:

$$w \geq b + \frac{e_H(r + q_f)}{q_f}\!\left(1 + \frac{r}{u/(1-u)}\right),$$

where $q_f$ is the probability of being caught shirking and $u$ is the unemployment rate. When unemployment is low, workers who are fired find new jobs quickly — the cost of shirking falls, and firms must pay higher wages to deter it. The NSC slopes upward in $(u, w)$ space.

The intersection of the NSC with the firm's labor demand curve $w = F_N(K,N)$ determines equilibrium employment and the efficiency wage. Crucially, the equilibrium unemployment rate exceeds the frictional natural rate: the economy is at $u > u^*_{frictional}$, and workers who are unemployed would accept jobs at the current wage if offered them — they are involuntarily unemployed. Unemployment is an equilibrium device for disciplining workers, not merely a transitional state between jobs.

---

*Next: Chapter 14 — Money Demand and Supply*
