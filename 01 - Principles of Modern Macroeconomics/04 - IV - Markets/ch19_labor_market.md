# Chapter 19 — The Labor Market: Wage Setting and Unemployment

> *"Labor markets clear only in economics textbooks."*
> — Alan Blinder

---

Of the three major markets examined in Part IV, the labor market is arguably the most important and the most controversial. It is where wages are determined and where the level of employment — and hence the unemployment rate — is established. It is also the market in which the most obvious departures from the competitive ideal are observed: wages do not instantly adjust to clear excess supply of labor (otherwise recessions would not cause unemployment), powerful institutions (unions, minimum wage laws, employment protection legislation) intervene between buyers and sellers, and informational asymmetries between employers and workers generate fundamental indeterminacies in wage determination.

This chapter develops the theory of labor market equilibrium in three frameworks of increasing complexity: the competitive benchmark, the search-and-matching model that provides the microeconomic foundations for the natural rate of unemployment, and the wage-setting/price-setting (WS-PS) framework that connects labor market outcomes to inflation dynamics. Together, these frameworks explain both the long-run level of unemployment (the natural rate) and the cyclical dynamics around it.

---

## 19.1 The Competitive Labor Market: A Benchmark and Its Failures

In the competitive benchmark, the real wage $w_t$ adjusts instantly to equate labor supply $N^S(w_t)$ and labor demand $N^D(w_t)$, giving a unique market-clearing employment level $N^* = N^S(w^*) = N^D(w^*)$ and zero involuntary unemployment. Every worker who wishes to work at the equilibrium wage is employed.

This benchmark is useful precisely because it is so clearly inconsistent with several important features of actual labor markets, and identifying the specific departures points toward better models.

**Failure 1: Wage rigidity and involuntary unemployment.** During recessions, employment falls sharply — by 5–10 million workers in the United States in 2008–09 — but wages fall very little. If markets were competitive and wages flexible, excess labor supply should drive wages down until employment recovered. Instead, the combination of persistent unemployment and sticky wages suggests that real wages are not the market-clearing price in the short run.

**Failure 2: Wage dispersion.** Even workers with similar observable skills earn very different wages across firms, industries, and regions. In a competitive market, wages should be equalized across equivalent workers by arbitrage. Persistent wage dispersion suggests frictions — the cost of searching, moving, or switching employers — that prevent workers from instantly accessing the best-paying job.

**Failure 3: Long-term relationships.** Labor contracts are typically long-term, explicit or implicit. Workers and firms invest in relationship-specific capital (skills, knowledge of the firm's processes, trust) that would be lost if the relationship ended. Long-term relationships create bilateral monopoly in the wage determination: neither the firm nor the worker can simply "take it or leave it."

---

## 19.2 The Search-and-Matching Model: Micro-Foundations of the Natural Rate

The search-and-matching model (Mortensen, 1970; Pissarides, 1985; Diamond, 1982) provides the most rigorous microeconomic foundation for equilibrium unemployment. The model has three key ingredients: a matching function, a job creation condition, and a wage-setting equation.

### The Matching Technology

Workers and firms find each other through a matching process that takes time and is subject to frictions. The aggregate matching function $m(U_t, V_t)$ gives the flow of new employment relationships per period as a function of the stock of unemployed workers $U_t$ and the stock of open vacancies $V_t$. The matching function is assumed to satisfy:

- $m_U > 0$, $m_V > 0$: more unemployed workers or more vacancies increases matches.
- Constant returns to scale: $m(\lambda U, \lambda V) = \lambda m(U,V)$ for all $\lambda > 0$.
- $m(0, V) = m(U, 0) = 0$: matches require both sides to be present.

The standard parametric form is **Cobb–Douglas matching**: $m(U,V) = sU^\alpha V^{1-\alpha}$, $\alpha \in (0,1)$.

**Definition (Labor Market Tightness).** **Labor market tightness** is $\theta_t \equiv V_t/U_t$, the ratio of vacancies to unemployed workers. With the Cobb–Douglas matching function, the **job-finding rate** for workers is:

$$f(\theta_t) = \frac{m(U_t, V_t)}{U_t} = s\theta_t^{1-\alpha}.$$

This is increasing in $\theta$: when vacancies are plentiful relative to unemployed workers, job-seekers find matches quickly. The **vacancy-filling rate** for firms is:

$$q(\theta_t) = \frac{m(U_t, V_t)}{V_t} = s\theta_t^{-\alpha}.$$

This is decreasing in $\theta$: when the market is tight, firms take longer to fill vacancies. The product $f(\theta)\cdot U = q(\theta)\cdot V = m(U,V)$ is consistent for any $\theta$.

### Job Creation Condition

A firm enters the labor market by posting a vacancy at flow cost $c > 0$ per period. When the vacancy is filled, the firm earns a match surplus — the value of the job to the firm. In steady state, the free-entry condition requires the expected value of posting a vacancy to equal zero (firms enter until the profit from posting is driven to zero):

$$\frac{c}{q(\theta^*)} = \frac{p - w}{r + \delta},$$

where $p$ is the flow productivity of a filled job (value of output per worker), $w$ is the wage, $r$ is the discount rate, and $\delta$ is the job-destruction rate. The left side is the expected cost of filling a vacancy (per-period cost divided by the vacancy-filling probability). The right side is the present value of the profit stream from a filled job. This **job creation (JC) condition** relates market tightness $\theta^*$ to the wage: a higher wage reduces the right side, requiring a lower $\theta$ (fewer vacancies relative to unemployed workers) to restore equilibrium.

### Wage Determination

Wages are determined by bilateral Nash bargaining between each matched worker-firm pair. The worker's surplus from employment is $W - U$, the difference between the value of employment and the value of unemployment. The firm's surplus is $J - V = J$ (since the free-entry condition drives the value of a vacancy to zero). The Nash bargaining solution:

$$w^* = \argmax_w\; (W - U)^\eta\, J^{1-\eta},$$

where $\eta \in [0,1]$ is the worker's **bargaining power** parameter. The solution (derived by taking the first-order condition with respect to $w$ and using the expressions for $W$, $U$, and $J$ in terms of the Bellman equations):

$$w^* = \eta(p + c\theta) + (1-\eta)(b + \rho),$$

where $b$ is the flow value of unemployment (benefits plus value of home production), $\rho$ is the worker's reservation utility, and $c\theta$ is the opportunity cost to the firm (since not hiring the worker means the vacancy continues to attract candidates). This **wage-setting (WS) condition** gives wages as a weighted average of the worker's outside option and the match surplus, weighted by bargaining power. Wages increase with $\theta$ (higher market tightness gives workers better outside options) and with $\eta$ (greater bargaining power).

### The Natural Rate

Combining the JC condition (which gives $\theta^*$ as a function of $w$) and the WS condition (which gives $w$ as a function of $\theta$) with the steady-state unemployment condition:

$$\delta(1 - u^*) = f(\theta^*)\, u^*,$$

(the flow of workers leaving employment equals the flow entering from unemployment), the system jointly determines the natural rate of unemployment $u^*$, market tightness $\theta^*$, and the equilibrium wage $w^*$:

$$u^* = \frac{\delta}{\delta + f(\theta^*)}.$$

The natural rate is increasing in $\delta$ (higher job destruction increases the flow into unemployment) and decreasing in $f(\theta^*)$ (faster matching reduces the stock of unemployed by shortening unemployment spells). The **Hosios (1990) condition** for efficiency is $\eta = \alpha$: when the worker's bargaining power equals the matching elasticity with respect to unemployment, the private equilibrium is socially optimal and the natural rate coincides with the efficient rate.

### The Beveridge Curve

The **Beveridge curve** plots the vacancy rate $v = V/L$ against the unemployment rate $u = U/L$ across time, tracing the locus of steady-state labor market conditions consistent with the matching technology. Using the steady-state condition and the Cobb–Douglas matching function:

$$\delta(1 - u) = su^\alpha v^{1-\alpha} \implies v = \left(\frac{\delta(1-u)}{su^\alpha}\right)^{1/(1-\alpha)}.$$

This defines a downward-sloping curve in $(u,v)$ space: recessions move the economy down and to the right (high unemployment, low vacancies); expansions move it up and to the left (low unemployment, high vacancies). Shifts in the curve itself — outward moves that correspond to higher unemployment for any given vacancy rate — indicate **structural deterioration** of the matching process: reduced matching efficiency, greater skill mismatch, increased geographic immobility, or higher search frictions. The dramatic outward shift of the U.S. Beveridge curve during 2021–22 — with both vacancies and unemployment high simultaneously — signaled significant sectoral reallocation frictions from the COVID pandemic.

---

## 19.3 The WS-PS Framework: Connecting Labor Market and Inflation

For the analysis of inflation dynamics, the wage-setting/price-setting (WS-PS) framework (Layard, Nickell, and Jackman, 1991; Blanchard and Katz, 1999) provides a tractable characterization of medium-run labor market equilibrium and its connection to the Phillips curve.

### The Wage-Setting Curve

The **wage-setting (WS) curve** gives the real wage that results from wage bargaining or efficiency-wage considerations as a function of unemployment. In the search model, the WS equation takes the form $w^{WS} = w(u, b, \sigma^2, \ldots)$ where wages decrease with unemployment (lower $u$ tightens the market, improving worker bargaining position). A simple log-linear approximation:

$$\ln w^{WS}_t = \ln w^0 - \gamma_u\ln u_t + \gamma_b\ln b_t + \epsilon_t^W,$$

where $b_t$ is the real value of unemployment benefits, $\gamma_u > 0$ captures the sensitivity of wages to unemployment, and $\gamma_b > 0$ captures the role of the replacement rate. This is an upward-sloping relationship in $(u, w)$ space: higher unemployment weakens worker bargaining and reduces wages.

### The Price-Setting Curve

Firms set prices as a markup $\mu$ over nominal marginal cost. With labor as the only variable input, nominal marginal cost is $W_t/\mathcal{A}_t$, where $\mathcal{A}_t$ is labor productivity. The profit-maximizing price:

$$P_t = (1 + \mu)\frac{W_t}{\mathcal{A}_t}.$$

Rearranging in terms of the real wage:

$$\frac{W_t}{P_t} = \frac{\mathcal{A}_t}{1+\mu} \equiv w^{PS}_t.$$

The **price-setting (PS) curve** gives the real wage consistent with profit maximization. It is horizontal in $(u, w)$ space (independent of unemployment in a simple model) at $w^{PS} = \mathcal{A}/(1+\mu)$: the real wage firms are willing to pay depends only on productivity and the markup, not on market conditions. More general formulations allow the markup to vary counter-cyclically (rising in recessions, falling in booms), making the PS curve slightly downward-sloping in unemployment.

### Medium-Run Equilibrium

Medium-run equilibrium in the labor market requires simultaneous satisfaction of both conditions:

$$w^{WS}(u^*) = w^{PS}.$$

This pins down the structural unemployment rate $u^*$ independently of the price level or aggregate demand — a crucial result. It says the natural rate of unemployment is determined by labor market institutions (bargaining power, union coverage), the generosity of unemployment insurance (replacement rate $b/w$), the degree of product market competition (markup $\mu$), and productivity $\mathcal{A}$. Monetary policy and fiscal policy can move $u$ temporarily away from $u^*$ (through the output gap and cyclical unemployment) but cannot change $u^*$ itself — except insofar as they affect the determinants of the WS and PS curves.

**Definition (Wage and Price Inflation in the Medium Run).** In the medium run, with full price and wage adjustment, the WS = PS equilibrium determines $u^*$. If actual unemployment $u_t \neq u^*$, there is persistent pressure on wages: $u_t < u^*$ generates accelerating wage inflation (workers have bargaining power above steady state); $u_t > u^*$ generates decelerating wage inflation. The resulting dynamic is the expectations-augmented Phillips curve:

$$\pi_t - \pi_t^e = -\alpha(u_t - u^*), \quad \alpha = \gamma_u/(1+\mu).$$

The coefficient $\alpha$ is jointly determined by labor market wage sensitivity ($\gamma_u$) and product market pricing behavior ($1+\mu$).

---

## 19.4 Labor Market Dynamics Over the Business Cycle

The adjustment of employment, wages, and unemployment over the business cycle is governed by the interaction of the matching process (governing flows in and out of employment) and the price- and wage-setting equations (governing the evolution of real wages relative to their equilibrium).

In the short run, aggregate demand shocks move the economy along the WS curve: a positive demand shock reduces unemployment below $u^*$, raising wages above $w^{PS}$. Firms respond to rising marginal costs by raising prices ($P$ rises), which — given a fixed nominal wage — reduces the real wage back toward $w^{PS}$. Eventually, workers renegotiate wages upward, restoring $w^{WS} = w^{PS}$ at $u = u^*$ — but at a higher price level. This medium-run adjustment mechanism is the labor-market foundation of the upward-sloping SRAS and the adjustment toward the vertical LRAS.

In the medium run, the key question is how fast wages adjust. Several empirical regularities bear on this. First, **nominal wage rigidity** — wages adjust more slowly downward than upward, a finding documented across many countries and decades. Second, **real wage rigidity** — the real wage responds sluggishly to unemployment, especially in European labor markets with high union coverage and employment protection legislation. Third, **efficiency wage premia** — many firms maintain wages above market-clearing for motivational reasons even when unemployment is high, generating persistent unemployment.

The interplay of these rigidities determines the sacrifice ratio: the greater the wage and price rigidity, the longer and more costly the disinflation when the central bank tightens policy to reduce inflation. This provides the connection between the labor market framework of this chapter and the Phillips curve analysis of Chapter 10.

---

*Next: Chapter 20 — Financial Markets*
