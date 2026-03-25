# Chapter 15 — Expectations and Behavioral Macroeconomics: Psychology in Economics

---

Every macroeconomic model contains, somewhere, an assumption about how agents form expectations of the future. This assumption is not peripheral — it is central. Whether fiscal stimulus raises output depends on whether households expect future taxes to offset current transfers. Whether monetary tightening reduces inflation quickly or slowly depends on whether wage setters expect inflation to fall. The choice of expectations formation mechanism is therefore one of the most consequential modeling choices in macroeconomics, and it has been one of the most contested.

---

## 15.1 Three Benchmarks for Expectations Formation

Three assumptions span the spectrum from naive to fully rational, and understanding their properties — and the empirical evidence for each — is prerequisite to evaluating any macroeconomic model.

**Definition (Static Expectations).** An agent has **static expectations** if she expects each variable to remain at its current level: $\mathbb{E}_t^S[x_{t+1}] = x_t$. This is the implicit assumption in early IS–LM and Phillips curve analyses. It is appropriate when variables change very slowly, but in environments with persistent trends or systematic dynamics, it implies large and predictable forecast errors — a violation of any reasonable notion of rationality.

**Definition (Adaptive Expectations).** An agent has **adaptive expectations** if she updates her forecast by a fraction $\lambda$ of the most recent forecast error:

$$\mathbb{E}_t^A[x_{t+1}] = \mathbb{E}_{t-1}^A[x_t] + \lambda\bigl(x_t - \mathbb{E}_{t-1}^A[x_t]\bigr), \quad \lambda \in (0,1).$$

Equivalently, this is a geometrically declining weighted average of past observations:

$$\mathbb{E}_t^A[x_{t+1}] = \lambda\sum_{k=0}^\infty (1-\lambda)^k x_{t-k}.$$

Adaptive expectations reduce but do not eliminate systematic forecast errors: when a variable is trending, adaptive expectations systematically underpredict (or overpredict) it, because past levels are poor guides to future levels. In the context of the Phillips curve, this is why the Friedman–Phelps model with adaptive expectations predicts only a temporary unemployment–inflation trade-off: eventually expectations adjust to match actual inflation, but the process of adjustment involves systematic forecast errors.

**Definition (Rational Expectations).** An agent has **rational expectations** if her forecasts are equal to the mathematical conditional expectation of the variable given all available information:

$$\mathbb{E}_t^R[x_{t+1}] = \mathbb{E}[x_{t+1} \mid \mathcal{F}_t],$$

where $\mathcal{F}_t$ is the information set at date $t$ and $\mathbb{E}[\cdot \mid \mathcal{F}_t]$ is the expectation computed using the model's true data-generating process. Rational expectations does not mean perfect foresight — agents can be surprised — but it rules out systematic, correctable forecast errors. If there is a pattern in your mistakes, you are not rational: you should exploit that pattern to improve your forecasts.

Rational expectations, as a methodological assumption, makes models internally consistent: agents behave as if they know the model. Its practical virtue is that it prevents the analyst from assuming that systematic policy rules can fool the public indefinitely.

---

## 15.2 Bounded Rationality: Simon, Satisficing, and Heuristics

Rational expectations imposes enormous cognitive demands. Optimally processing all available information — including the correct model of the economy, the government's policy rule, and all other agents' behavior — requires computational and informational resources that real agents demonstrably do not have. Herbert Simon (1955) proposed **bounded rationality** as a more realistic description: agents satisfice rather than optimize, using simplified rules of thumb that perform adequately in most states of the world but may generate systematic biases.

**Definition (Bounded Rationality).** An agent is **boundedly rational** if she uses simplified decision procedures — heuristics — that reduce the cognitive burden of optimization at the cost of some accuracy. Bounded rationality does not mean irrational: a heuristic is rational if the cost of acquiring and processing additional information exceeds the benefit of improved decisions.

Xavier Gabaix (2020) formalizes bounded rationality in a macroeconomic context via the **Behavioral New Keynesian Model**. Agents behave as if they apply a cognitive discount factor $M \in [0,1]$ to future variables: they discount future macroeconomic conditions more heavily than a fully rational agent would. The behavioral IS curve and Phillips curve are:

$$\hat{x}_t = M\,\mathbb{E}_t[\hat{x}_{t+1}] - \sigma(i_t - \mathbb{E}_t[\pi_{t+1}] - r_t^n),$$

$$\hat{\pi}_t = \beta M_f\,\mathbb{E}_t[\hat{\pi}_{t+1}] + \kappa\,\hat{x}_t,$$

where $M, M_f \in [0,1]$ are cognitive discount factors. When $M = M_f = 1$ the standard New Keynesian system is recovered; when $M < 1$, future variables are discounted more heavily, attenuating the forward-guidance effects of announced future policy changes. This directly addresses the **forward-guidance puzzle**: in the standard New Keynesian model, an announcement that interest rates will be lower in five years has implausibly large effects on current output; in the behavioral model, these effects are attenuated by cognitive discounting.

---

## 15.3 Extrapolative Expectations and Asset Prices

A particularly important departure from rational expectations in financial markets and housing is **extrapolation**: the tendency for agents to project recent trends into the future, generating momentum in asset prices that exceeds what fundamentals justify.

**Definition (Extrapolative Expectations).** An agent has **extrapolative expectations** if she expects future values to continue recent trends: $\mathbb{E}_t^E[x_{t+1}] = x_t + \mu(x_t - x_{t-1})$ with $\mu > 0$. When $\mu > 0$, positive recent growth generates positive expectations of future growth — a self-reinforcing dynamic. When $\mu < 0$, expectations are mean-reverting.

Applied to housing prices: if $\mu > 0$, a period of rising house prices generates expectations of further rises, attracting additional buyers, driving prices further up — a mechanism that can sustain house prices above fundamental value for extended periods. This extrapolative dynamic has been documented empirically by Case, Shiller, and Thompson (2012) using survey data on home buyers' expectations, and it is a central element of explanations for the U.S. housing bubble of 2002–06 (Chapter 40).

---

## 15.4 Animal Spirits and Coordination Failures

Keynes argued that investment is driven not only by expected returns but by "animal spirits" — the spontaneous confidence or pessimism of entrepreneurs in the face of genuine, irreducible uncertainty about the future. This concept was long dismissed as literary rather than scientific, but it has been given rigorous foundations in the theory of multiple equilibria and coordination games.

**Definition (Sunspot Equilibrium).** A **sunspot equilibrium** is a self-fulfilling equilibrium in which the economy coordinates on a particular outcome based on an extraneous, payoff-irrelevant signal — a "sunspot." The sunspot has no effect on fundamentals (technology, preferences, endowments) but nonetheless influences outcomes because agents believe others are paying attention to it, and their belief is self-confirming. Sunspot equilibria require that the economic model has multiple equilibria in the first place; the sunspot merely selects among them.

In a model where aggregate demand depends on expectations of aggregate demand — because higher expected output makes investment and consumption more attractive — multiple Nash equilibria can exist: one in which all agents are optimistic, investment is high, output is high, and expectations are confirmed; and one in which all are pessimistic, investment is low, output is low, and pessimism is confirmed. A sunspot — yesterday's stock market return, a newspaper headline, or a charismatic politician's speech — can shift the economy between these equilibria, generating aggregate fluctuations entirely unrelated to changes in fundamentals (Farmer, 1993). This provides a potential explanation for the observation that business cycles often seem driven by changes in sentiment that precede, rather than follow, changes in observable fundamentals.

---

*Next: Chapter 16 — Rational Expectations and New Classical Economics*
