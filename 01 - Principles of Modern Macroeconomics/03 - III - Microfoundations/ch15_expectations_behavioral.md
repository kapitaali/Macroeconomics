# Chapter 15 — Expectations and Behavioral Macroeconomics: Psychology in Economics

> *"The expected never happens; it is the unexpected always."*
> — John Maynard Keynes, *A Tract on Monetary Reform*, 1923

---

Every macroeconomic model contains, somewhere, an assumption about how agents form expectations of the future. This assumption is not peripheral — it is structural. Whether fiscal stimulus raises output depends on whether households expect future taxes to offset current transfers. Whether monetary tightening reduces inflation quickly or slowly depends on whether wage setters expect the central bank to succeed. Whether a housing boom continues or collapses depends on whether buyers expect prices to keep rising. The expectations formation mechanism is therefore one of the most consequential modeling choices in macroeconomics, and it has been one of the most contested. This chapter develops the three main frameworks — static, adaptive, and rational expectations — examines the empirical and theoretical challenges to each, and surveys the behavioral economics research program that attempts to model the actual psychology of expectation formation with rigor rather than dismissing it as noise.

---

## 15.1 Three Benchmarks for Expectations Formation

Three assumptions span the spectrum from naive to fully rational, and understanding their properties — and the empirical evidence for each — is prerequisite to evaluating any macroeconomic model.

**Definition (Static Expectations).** An agent has **static expectations** if she expects each variable to remain at its current level: $\mathbb{E}_t^S[x_{t+1}] = x_t$. Static expectations are the implicit assumption of early IS–LM and original Phillips curve analyses. They are appropriate when variables change very slowly relative to the frequency of decisions, but in environments with persistent trends or systematic dynamics, they produce large and **predictable** forecast errors — errors that a rational agent would exploit to improve forecasts. Static expectations are therefore inconsistent with any forward-looking optimization.

**Definition (Adaptive Expectations).** An agent with **adaptive expectations** updates her forecast by a fraction $\lambda$ of the most recent forecast error:

$$\mathbb{E}_t^A[x_{t+1}] = \mathbb{E}_{t-1}^A[x_t] + \lambda\bigl(x_t - \mathbb{E}_{t-1}^A[x_t]\bigr), \quad \lambda \in (0,1).$$

Equivalently, adaptive expectations form a geometrically declining weighted average of all past observations:

$$\mathbb{E}_t^A[x_{t+1}] = \lambda\sum_{k=0}^\infty (1-\lambda)^k x_{t-k}.$$

Adaptive expectations substantially improve on static expectations when variables exhibit mean-reversion or slow trends: recent observations receive more weight, and the forecast catches up to reality. But they produce **systematic** forecast errors whenever variables follow persistent trends. In the context of the 1960s Phillips curve: when inflation rose persistently through 1965–1970, adaptive expectations systematically underpredicted inflation, implying that workers consistently accepted wage settlements lower than the correct inflation-adjusted rate. Eventually, as workers updated their expectations, the real wage moved toward equilibrium, and the unemployment gain dissipated. This is the Friedman-Phelps insight [Ch. 10]: adaptive expectations support only a temporary, expectations-distorted Phillips curve trade-off.

**Definition (Rational Expectations).** An agent has **rational expectations** if her forecast equals the mathematical conditional expectation of the variable given all available information:

$$\mathbb{E}_t^R[x_{t+1}] = \mathbb{E}[x_{t+1} \mid \mathcal{F}_t],$$

where $\mathcal{F}_t$ is the information set available at date $t$ and the expectation is computed using the model's true data-generating process. Rational expectations does not imply perfect foresight — agents can be surprised by genuinely random shocks — but it rules out **systematic, correctable forecast errors**: if there is a predictable pattern in your mistakes, you are not rational.

Rational expectations imposes internal consistency on a model: agents behave as if they know the model and the policy rule. This is a discipline on the analyst, preventing sleight of hand in which the public is assumed to be fooled by systematic policy actions. Its main practical virtue is that it prevents the assumption that predictable government policies can have indefinitely large real effects by exploiting an informationally disadvantaged public.

---

## 15.2 The Empirical Record on Expectations

The theoretical case for rational expectations rests on internal consistency and equilibrium discipline; the empirical case must rest on data. A substantial literature has tested whether actual expectations, as measured in surveys of households, firms, and professional forecasters, conform to the rational benchmark.

The evidence is decidedly mixed. **Professional forecasters** — economists at central banks, investment banks, and international organizations — have expectations that approximate rationality: their forecast errors are roughly uncorrelated with lagged information, consistent with the orthogonality test for rational expectations. The Federal Reserve's Greenbook forecasts, for example, are difficult to improve upon using publicly available information at the time they were produced.

**Household and firm expectations** are another matter. Survey evidence documents systematic departures from rationality. Households in the Michigan Survey of Consumer Sentiment consistently overpredict inflation compared to professional forecasters, with the overestimation particularly pronounced for lower-income and less-educated respondents. Business expectations of investment and sales exhibit extrapolation: Greenwood and Hanson (2015) show that firms over-invest in expanding sectors and under-invest in contracting ones, consistent with extrapolative rather than rational expectations about future profitability. Kohlhas and Broer (2022) document that households' inflation expectations respond more to salient personal experience — recently purchased goods, own housing costs — than to the aggregate data rational expectations requires agents to use.

These departures are not random noise. They are systematic in ways that make theory both harder (agents are not rational) and richer (the specific ways they deviate from rationality are themselves predictable and have macroeconomic implications).

---

## 15.3 Bounded Rationality: Simon, Satisficing, and Heuristics

Rational expectations imposes enormous cognitive and informational demands. Optimally processing all available information — including the correct model of the economy, the government's policy rule, other agents' strategies, and the complete probability distributions over future shocks — requires computational resources and informational access that real agents demonstrably do not possess. Herbert Simon (1955) proposed **bounded rationality** as a more realistic description: agents satisfice rather than optimize, using simplified decision rules that perform adequately in most states of the world but may generate systematic biases.

**Definition (Bounded Rationality).** An agent is **boundedly rational** if she uses simplified decision procedures — heuristics — that reduce the cognitive burden of optimization at the cost of some accuracy. A heuristic is itself rational in a deeper sense if the cost of acquiring and processing the additional information needed for full optimization exceeds the benefit from improved decisions. Bounded rationality does not mean irrational: it means rationally economizing on the costs of cognition.

Xavier Gabaix (2020) formalizes bounded rationality in a New Keynesian context through the **Behavioral New Keynesian Model**. Agents behave as if they apply a **cognitive discount factor** $M \in [0,1]$ to future macroeconomic variables: they underweight the future relative to a fully rational agent, treating distant future states as less relevant for current decisions. The behavioral IS curve and Phillips curve become:

$$\hat{x}_t = M\,\mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\pi_{t+1}] - r_t^n),$$

$$\hat{\pi}_t = \beta M_f\,\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\,\hat{x}_t,$$

where $M, M_f \in [0,1]$ are cognitive discount factors for output and inflation expectations respectively. When $M = M_f = 1$, the standard New Keynesian system is recovered. When $M < 1$, current output is less sensitive to expectations of future output and interest rates: the forward-looking IS relationship is attenuated.

The behavioral NK model directly addresses the **forward-guidance puzzle**: in the standard New Keynesian model, an announcement that interest rates will be held lower for five years has implausibly large effects on current output — the model's extreme sensitivity to distant future expectations defies both intuition and the measured modest effects of actual forward guidance. With cognitive discounting, the further away in time a policy promise is, the more steeply it is discounted by boundedly rational agents, generating responses to forward guidance that are both smaller and more realistic.

---

## 15.4 Extrapolative Expectations and Asset Prices

A particularly important and empirically well-documented departure from rational expectations in financial and housing markets is **extrapolation**: the tendency for agents to project recent trends forward, generating momentum in asset prices beyond what fundamentals justify.

**Definition (Extrapolative Expectations).** An agent has **extrapolative expectations** if she expects the future value of an asset to continue recent trends: $\mathbb{E}_t^E[p_{t+1}] = p_t + \mu(p_t - p_{t-1})$ with $\mu > 0$. Recent price growth creates expectations of further price growth — a self-reinforcing dynamic. When $\mu < 0$, expectations are mean-reverting and stabilizing; when $\mu > 0$, they are trend-following and potentially destabilizing.

Applied to housing markets: when $\mu > 0$, a period of rising house prices generates expectations of further increases, attracting additional buyers, driving prices up further, and attracting still more buyers in a feedback loop. This dynamic can sustain asset prices substantially above fundamental value for extended periods — until a triggering event reverses the trend and extrapolative expectations become destabilizing on the downside. Case, Shiller, and Thompson (2012) provide direct survey evidence of extrapolative housing expectations in the U.S. before and during the 2002–06 bubble: homebuyers in the bubble markets expected 10-year annual price appreciation of 10–15%, rates with no historical precedent, consistent with projection of the recent trend rather than mean-reversion toward long-run fundamentals.

Extrapolative expectations also feature prominently in asset pricing. Barberis et al. (2018) develop a model in which investors extrapolate past earnings growth into future earnings forecasts, generating the "excessive extrapolation" documented in equity markets: stocks with recent strong price performance attract capital that further bids up prices, creating momentum; eventual reversion generates the predictable negative returns to momentum strategies. Greenwood and Shleifer (2014) show that survey measures of stock market expectations are positively correlated with recent returns — directly measuring the extrapolation — and negatively correlated with subsequent returns — confirming the destabilizing effect on prices.

---

## 15.5 Animal Spirits and Coordination Failures

Keynes argued in the *General Theory* (1936) that investment is driven not only by expected returns but by "animal spirits" — the spontaneous confidence or pessimism of entrepreneurs facing genuine, irreducible uncertainty about the future. When long-run expectations are dominated by convention ("things will be much as they are") and sentiment ("the general opinion"), rather than by careful calculation of expected returns, aggregate investment becomes susceptible to large swings on thin evidential foundations. This was long treated as literary insight rather than formal theory, but modern equilibrium analysis has provided it with rigorous foundations.

**Definition (Sunspot Equilibrium).** A **sunspot equilibrium** is a rational expectations equilibrium in which economic outcomes depend on an extraneous, payoff-irrelevant signal — a "sunspot." The signal has no direct effect on technology, preferences, or endowments, but it influences outcomes because agents believe others are coordinating on it, and that belief is self-confirming. Sunspot equilibria exist only when the economy has multiple equilibria that the sunspot can select among.

In a model where aggregate demand depends on expectations of future aggregate demand — as in the New Keynesian IS curve, where current output depends on expected future output — multiple Nash equilibria can coexist: a high-output equilibrium in which optimistic expectations drive high investment, generating high output that validates the optimism; and a low-output equilibrium in which pessimism suppresses investment and demand, generating low output that confirms the pessimism. Farmer (1993, 2012) develops a full macroeconomic framework in which animal spirits — self-fulfilling beliefs about employment — can generate business cycle fluctuations entirely unrelated to changes in productivity, monetary policy, or other "fundamental" driving forces.

This framework provides a potential account of several empirical puzzles that fundamentals-driven models struggle to explain: why business cycle turning points are so difficult to predict from observable fundamentals; why recessions sometimes seem to feed on themselves through confidence channels long after any initial triggering event has passed; and why sentiment indicators (CEO confidence surveys, consumer confidence indices) have predictive power for investment and consumption beyond what is captured in financial variables.

### The Confidence Channel in Policy

The policy implications of animal spirits and sunspot multiplicity are both important and uncomfortable. If recessions can be triggered and sustained by self-fulfilling pessimism, then policies that restore confidence — announcements, institutional commitments, communication — may have real effects independent of the financial variables through which they operate. Mario Draghi's 2012 commitment to do "whatever it takes" to preserve the euro, which reversed the European sovereign debt crisis without any central bank actually purchasing sovereign debt, is perhaps the clearest modern example: the announcement shifted expectations from the bad to the good equilibrium, and the self-fulfilling improvement in bond markets validated the announcement.

Conversely, policy mistakes that trigger pessimism — premature fiscal consolidation, credibility-damaging reversals — may have costs that purely mechanistic models understate, because the negative-sentiment equilibrium they trigger can be self-reinforcing in ways the single-equilibrium model does not capture.

---

## Chapter Summary

- Three benchmarks span the expectations spectrum: **static expectations** ($\mathbb{E}_t^S[x_{t+1}] = x_t$, generating predictable errors in trending environments); **adaptive expectations** ($\mathbb{E}_t^A[x_{t+1}] = \lambda\sum(1-\lambda)^k x_{t-k}$, generating systematic errors during persistent trend changes); and **rational expectations** ($\mathbb{E}_t^R[x_{t+1}] = \mathbb{E}[x_{t+1}|\mathcal{F}_t]$, ruling out correctable errors).

- The empirical record is mixed: professional forecasters approximate rationality; household and firm expectations exhibit systematic biases including inflation overestimation, extrapolation, and anchoring to personal experience.

- **Bounded rationality** (Gabaix, 2020) formalizes cognitive limitations through a cognitive discount factor $M < 1$, generating behavioral IS and Phillips curves with attenuated forward-looking components. This resolves the forward-guidance puzzle by making the effect of distant future policy promises smaller and more realistic.

- **Extrapolative expectations** ($\mathbb{E}_t^E[p_{t+1}] = p_t + \mu(p_t - p_{t-1})$ with $\mu > 0$) generate self-reinforcing price dynamics in asset markets, documented empirically in housing (Case-Shiller-Thompson) and equities (Greenwood-Shleifer). They are a central mechanism in boom-bust cycles.

- **Animal spirits** and **sunspot equilibria** provide rigorous foundations for Keynes's insight about confidence: when models have multiple equilibria, payoff-irrelevant signals can coordinate expectations, generating fluctuations disconnected from fundamentals and giving communication policy direct real effects.

---

*Next: Chapter 16 — Rational Expectations and New Classical Economics*
