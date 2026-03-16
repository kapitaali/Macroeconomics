# Chapter 2 — A Brief History of Macroeconomic Thought: From Classical to Keynesian to Modern

> *"Practical men, who believe themselves to be quite exempt from any intellectual influence, are usually the slaves of some defunct economist."*
> — John Maynard Keynes, *The General Theory*, 1936

---

The history of macroeconomic thought is not an antiquarian exercise. Every major theoretical development in the discipline was a response to an empirical anomaly or a policy failure that the prevailing framework could not accommodate. Understanding this dialectic — in which crises expose the limits of existing theory, driving innovations that in turn reveal new anomalies — is essential for evaluating the current state of the field and maintaining an appropriately humble attitude toward the models we use today. Each generation of economists has believed it had largely solved the core problems; each generation has been proved wrong in important ways. Humility about what we do not yet understand is not a counsel of despair but a precondition for progress.

---

## 2.1 Classical Political Economy and the Pre-Keynesian Consensus

The economists we now call "classical" — Adam Smith, David Ricardo, John Stuart Mill, and their contemporaries in the late eighteenth and nineteenth centuries — did not have a formal macroeconomic model. But their writings embody a coherent set of propositions about how aggregate economies work, propositions whose influence persisted for over a century and whose echoes can still be heard today.

The most important of these propositions is **Say's Law**, named after the French economist Jean-Baptiste Say. In its most common formulation, Say's Law holds that supply creates its own demand: the act of producing goods generates the income necessary to purchase them, so there can be no general excess supply — no situation in which producers collectively wish to sell more than buyers collectively wish to purchase. Excess supply in one market must, by arithmetic, be accompanied by excess demand in another. An economy may have sectoral mismatches, but it cannot be "depressed" in the aggregate.

The formal underpinning of Say's Law is **Walras' Law**, which states that the value of all excess demands sums to zero across all markets. Let there be $n$ markets with excess demands $z_i(p_1,\ldots,p_n)$ for $i = 1,\ldots,n$. Then at any price vector:

$$\sum_{i=1}^n p_i\, z_i(p_1,\ldots,p_n) = 0.$$

This is an accounting identity, not a statement that each individual market clears. It holds because every purchase is someone's sale: the budget constraints of all agents sum to zero in value. The classical misinterpretation was to read this identity as ruling out aggregate excess supply of produced goods — a logical error, since Walras' Law is consistent with any individual market being in disequilibrium, provided the values of all disequilibria sum to zero. The refutation of this misinterpretation, and the identification of the specific conditions under which general deficiency of demand is possible, was Keynes's central contribution.

Alongside Say's Law, the classical tradition held that prices and wages are fully flexible, adjusting rapidly to clear all markets. Under this assumption, unemployment is always voluntary — workers who want work at the market wage can find it — and recessions are either impossible (if Say's Law holds) or transitory (if they reflect temporary mismatches that price adjustment will eliminate). The policy implication was clear: government intervention to stimulate demand was unnecessary and potentially counterproductive, since markets would self-correct.

---

## 2.2 The Quantity Theory of Money

The classical treatment of the price level rests on what is called the **quantity theory of money**, a doctrine with roots reaching back to David Hume in the eighteenth century. The quantity theory asserts that the price level is determined, in the long run, by the supply of money.

**Definition (The Quantity Theory of Money).** The quantity theory holds that the price level $P$ is proportional to the money supply $M$, with the factor of proportionality determined by the velocity of money $V$ and the level of real output $Y$. This is captured by the **Fisher equation** (after Irving Fisher, who gave it its modern form):

$$MV = PY,$$

where $M$ is the nominal money supply, $V$ is the **income velocity of money** (the average number of times a unit of currency is used in transactions per period), $P$ is the price level, and $Y$ is real output. The velocity $V$ is defined by this equation itself; it becomes an independent empirical claim only when supplemented by a theory of what determines $V$.

Taking logarithms and differentiating with respect to time, the Fisher equation implies:

$$\hat{m} + \hat{v} = \pi + \hat{y},$$

where $\hat{x} \equiv \mathrm{d}\ln x / \mathrm{d}t$ denotes the growth rate of variable $x$. If velocity is approximately constant ($\hat{v} \approx 0$) and real output growth is determined by supply-side factors independently of money ($\hat{y}$ exogenous), then $\pi \approx \hat{m}$: the inflation rate equals the money growth rate. This is the essence of the quantity theory and it is one of the most robustly supported propositions in macroeconomics over long horizons (McCandless and Weber, 1995).

The quantity theory implies a **classical dichotomy**: real variables (output, employment, relative prices) are determined independently of nominal variables (the money supply, the price level). A doubling of $M$ doubles $P$ but leaves $Y$, $V$, wages, and all relative prices unchanged. Money is said to be **neutral**: it affects the scale of nominal transactions but not the real allocation of resources. This neutrality holds as a long-run proposition in modern macroeconomics, but fails in the short run precisely because of the nominal rigidities discussed in Chapter 1.

---

## 2.3 Wicksell and the Natural Rate of Interest

Before Keynes, the most important departure from the simple quantity theory was made by the Swedish economist Knut Wicksell. In *Geldzins und Güterpreise* (1898), Wicksell identified a tension that the static quantity theory could not resolve: if the price level is determined by $M$ and $V$, what determines the allocation of resources between consumption and investment? And what connects the financial sector — banks, interest rates, credit — to the real economy?

Wicksell's answer introduced the concept of the **natural rate of interest**.

**Definition (Natural Rate of Interest).** The **natural rate of interest** $r^n$ is the real interest rate at which saving equals investment in a hypothetical economy without money — the rate that would prevail in a barter economy and that reflects the marginal productivity of capital. The **market rate of interest** $r^m$ is the rate actually charged by banks on loans.

Wicksell argued that a divergence between the two rates generates a cumulative process of expansion or contraction. If $r^m < r^n$, borrowing is cheap relative to the productive opportunities available: firms borrow and invest aggressively, demand exceeds supply at current prices, and inflation accelerates. If $r^m > r^n$, the opposite obtains: investment is deterred, demand falls short of supply, and a deflationary spiral begins. This **Wicksellian mechanism** anticipates the modern concept of the output gap, the Taylor rule, and the neutral interest rate — topics that recur throughout Parts V and VI.

Irving Fisher built on Wicksell's insights to formalize the intertemporal consumption decision. Fisher's two-period model of a household with income endowments $y_1$ and $y_2$ is:

$$\max_{c_1, c_2}\; u(c_1) + \beta\,u(c_2) \quad \text{subject to} \quad c_1 + \frac{c_2}{1+r} = y_1 + \frac{y_2}{1+r},$$

where $\beta \in (0,1)$ is the **subjective discount factor** — a measure of the household's impatience, with lower $\beta$ reflecting greater preference for present consumption. The first-order condition:

$$u'(c_1) = \beta(1+r)\,u'(c_2),$$

is called the **Euler equation** for consumption. It says that at the optimum, the marginal utility of consumption today must equal the discounted marginal utility of consumption tomorrow, scaled by the gross return on saving. If $\beta(1+r) > 1$ — if the return on saving exceeds the rate of impatience — the household postpones consumption and saving is positive. This two-period Euler equation is the embryonic form of the intertemporal optimization that underlies all of Part III.

---

## 2.4 Keynes and the General Theory

The dominant influence on macroeconomics in the twentieth century is John Maynard Keynes, whose *General Theory of Employment, Interest and Money* (1936) was written in the shadow of the Great Depression. Between 1929 and 1933, U.S. output fell by nearly one-third and unemployment rose to 25%. The classical framework had no satisfactory explanation: it predicted either that depressions were impossible (Say's Law) or that they would be self-correcting through falling wages and prices. Wages and prices did fall, yet recovery did not come. Something was fundamentally wrong with the classical model.

Keynes attacked the classical framework on two fronts. The first was that money is not merely a medium of exchange but a store of value with unique liquidity properties. In a world of uncertainty, households and firms may prefer to hold money — which is certain to retain its face value — rather than lending at any positive interest rate, if they fear that the asset prices of bonds or equities may fall. This is the **liquidity trap**: a situation in which monetary expansion fails to reduce interest rates because additional money is simply hoarded rather than lent.

**Definition (Liquidity Trap).** A liquidity trap occurs when nominal interest rates reach their lower bound (approximately zero) and the demand for money becomes perfectly elastic with respect to the interest rate. In this situation, open-market operations that increase the money supply have no effect on interest rates or aggregate demand, because the public is willing to hold any amount of money at the prevailing interest rate. The liquidity trap makes monetary policy ineffective and creates a potential role for fiscal policy.

The second front was investment. Keynes argued that investment is driven not by the smooth marginal productivity calculations of the neoclassical model but by "animal spirits" — the spontaneous waves of optimism and pessimism that characterize entrepreneurial decisions in a world of genuine, irreducible uncertainty. Investment can collapse for reasons entirely unrelated to any change in fundamentals, and when it does, the economy can settle at an equilibrium with high unemployment indefinitely, because there is no automatic mechanism to restore full employment.

The core of Keynes's positive model can be stated concisely. Goods-market equilibrium requires that actual output equals desired expenditure:

$$Y = C(Y - T) + I(r) + G,$$

where $C(\cdot)$ is the consumption function — the relationship between disposable income $Y - T$ and desired consumption — $I(r)$ is investment, a decreasing function of the real interest rate $r$, $G$ is government spending, and $T$ is lump-sum taxes. Money-market equilibrium requires that real money demand $L(Y, i)$ — which increases with income and decreases with the nominal interest rate $i$ — equals the real money supply:

$$\frac{M}{P} = L(Y, i).$$

Combined with the Fisher relation $i = r + \pi^e$ linking the nominal rate, the real rate, and expected inflation $\pi^e$, these two equations in two unknowns $(Y, r)$ constitute the IS–LM model, which Hicks (1937) formalized and which is developed fully in Chapter 9. It is important to understand the IS–LM model not as a timeless truth but as a particular theoretical construct — one that captures important short-run dynamics but abstracts from many features (explicit intertemporal optimization, endogenous expectations, financial sector) that subsequent research has shown to be important.

---

## 2.5 The Neoclassical Synthesis and the Phillips Curve

The postwar synthesis, associated with Paul Samuelson, Franco Modigliani, and Robert Solow, attempted to reconcile Keynesian and classical ideas. The synthesis held that Keynesian demand management was necessary in the short run, when prices are sticky and markets do not clear, but that classical long-run results held once prices adjusted. This was a compromise position — arguably too comfortable, because it deferred rather than resolved the tension between sticky-price Keynesian models and the classical theory of long-run equilibrium.

The synthesis was operationalized through an empirical relationship identified by A.W. Phillips (1958). Using annual data on British wage inflation and unemployment from 1861 to 1957, Phillips documented a stable, negatively sloped curve relating wage inflation to unemployment:

$$\pi_t^w = f(u_t), \quad f' < 0.$$

Samuelson and Solow (1960) reinterpreted Phillips's result as a trade-off between price inflation and unemployment that policymakers could exploit. The implication seemed clear: by accepting modestly higher inflation, governments could permanently maintain lower unemployment. This apparent policy menu generated enormous optimism about the capacity of economic management to achieve full employment, and it shaped macroeconomic policy in many countries throughout the 1960s.

The problem with this interpretation — which became glaringly apparent during the stagflation of the 1970s, when both unemployment and inflation rose simultaneously — was that it confused a correlation that held under a specific monetary regime with a stable structural relationship. This confusion was the target of Friedman's and Phelps's contributions, discussed next.

---

## 2.6 Friedman, Phelps, and the Expectations Augmentation

In 1968, Milton Friedman and Edmund Phelps independently published papers that demolished the policy implications of the original Phillips curve. Their argument was both simple and devastating. The original Phillips curve relates nominal wage growth to unemployment. But workers and firms care about *real* wages — wages adjusted for inflation. A worker who receives a 5% wage increase when inflation is 5% has received no real improvement; if she expects 5% inflation, her wage demand already incorporates this expectation. The original Phillips curve, by omitting inflation expectations, was only valid if workers systematically failed to anticipate inflation — a condition that could not persist indefinitely.

Friedman and Phelps each proposed the **expectations-augmented Phillips curve** (EAPC):

$$\pi_t = \pi_t^e - \alpha(u_t - u^*) + \epsilon_t, \quad \alpha > 0,$$

where $\pi_t^e$ is expected inflation, $u^*$ is the **natural rate of unemployment** (defined in Chapter 1 and developed formally in Chapters 3 and 13), and $\epsilon_t$ captures supply-side disturbances. The critical implication is that the trade-off between inflation and unemployment exists only to the extent that actual inflation exceeds expected inflation — only when the economy is "fooled." In the long run, when $\pi_t^e$ adjusts to equal $\pi_t$, the EAPC implies $u_t = u^*$ regardless of $\pi_t$: the long-run Phillips curve is vertical, and there is no permanent trade-off between inflation and unemployment.

Under **adaptive expectations** — the assumption that $\pi_t^e = \pi_{t-1}$ (agents expect this period's inflation to equal last period's) — the EAPC becomes:

$$\pi_t - \pi_{t-1} = -\alpha(u_t - u^*) + \epsilon_t.$$

This is the **accelerationist Phillips curve**: when unemployment is held below its natural rate, inflation does not merely remain high — it continuously accelerates. The only inflation rate consistent with stable unemployment equal to $u^*$ is one that is neither accelerating nor decelerating. This is the origin of the term **NAIRU** (Non-Accelerating Inflation Rate of Unemployment), which will be developed formally in Chapters 3 and 10.

The stagflation of 1973–75, in which the first oil shock drove unemployment and inflation simultaneously upward, provided dramatic empirical confirmation of the Friedman–Phelps framework. The simple Phillips curve trade-off collapsed, and with it the intellectual foundations of the activist demand management that had characterized economic policy in the 1960s.

---

## 2.7 The Rational Expectations Revolution

The Friedman–Phelps critique assumed adaptive expectations — backward-looking, potentially improvable rules of thumb. Robert Lucas (1972, 1976) and Thomas Sargent and Neil Wallace (1975) pushed this critique to a stronger conclusion by replacing adaptive expectations with **rational expectations**.

**Definition (Rational Expectations).** An agent has rational expectations if her beliefs about future variables are equal to the mathematical conditional expectation of those variables given all available information:

$$\mathbb{E}_t^R[x_{t+1}] = \mathbb{E}[x_{t+1} \mid \mathcal{F}_t],$$

where $\mathcal{F}_t$ is the information set available at date $t$ and $\mathbb{E}[\cdot \mid \mathcal{F}_t]$ is the expectation formed using the model's true data-generating process. Rational expectations does not mean perfect foresight — agents can be surprised — but it rules out *systematic*, *correctable* forecast errors. If there is a pattern in your forecast errors, you are not rational, because you should exploit that pattern to improve your forecasts.

Rational expectations is a powerful assumption because it makes the model internally consistent: agents behave as if they know the model. It imposes discipline on theory-building by preventing the analyst from assuming that agents can be systematically fooled by a policy rule that is publicly known. Under rational expectations, only *unanticipated* changes in policy can have real effects — systematic, predictable policy moves are fully anticipated and already incorporated in prices and wages before they occur. This is the **policy ineffectiveness proposition** of Sargent and Wallace.

The rational expectations hypothesis also provided the foundation for Lucas's critique of econometric policy evaluation (defined formally in Chapter 1 and elaborated in Chapter 16). If estimated reduced-form relationships embed agents' beliefs about the current policy regime, then those relationships will shift when the regime changes — making historical estimates useless for predicting the effects of policy changes.

---

## 2.8 Real Business Cycles

Kydland and Prescott (1982) took the rational expectations, market-clearing framework to its logical extreme. In their **real business cycle (RBC) model**, business cycles are not failures of demand management or misallocations caused by monetary disturbances — they are the *optimal* responses of households and firms to exogenous fluctuations in total factor productivity (technology shocks). The definition is important enough to state explicitly.

**Definition (Real Business Cycle).** A real business cycle model is a fully microfounded dynamic general equilibrium model in which business cycle fluctuations are driven entirely by technology shocks (shocks to total factor productivity) and in which prices and wages are fully flexible, clearing all markets at every date. Because all agents optimize and all markets clear, the equilibrium is Pareto optimal: there is no role for stabilization policy, because fluctuations are efficient responses to genuine productivity changes.

The RBC approach was methodologically transformative: it established that business cycle models should be built from first principles, with explicit utility-maximizing households and profit-maximizing firms, and that the models' predictions should be evaluated quantitatively against the data rather than merely qualitatively. This methodological legacy — the DSGE approach — has become standard even among economists who reject RBC conclusions about the role of technology shocks.

---

## 2.9 The New Keynesian Synthesis

The New Keynesian program, developed during the 1980s and 1990s by Rotemberg, Mankiw, Woodford, Calvo, Galí, and Gertler, accepted the RBC methodological requirement of microfounded dynamic models with rational expectations, but added two ingredients that restore short-run non-neutrality of monetary policy: **monopolistic competition** (firms set prices rather than taking them from a Walrasian auctioneer) and **nominal rigidities** (firms cannot costlessly adjust prices every period).

**Definition (Nominal Rigidity — Calvo Pricing).** In Calvo's (1983) model of staggered price setting, in each period a randomly chosen fraction $1-\theta$ of firms can reset their price optimally. The remaining fraction $\theta$ must keep their price unchanged. The parameter $\theta \in [0,1]$ measures the degree of price stickiness: $\theta = 0$ means all firms reprice every period (no rigidity, classical case); $\theta \to 1$ means prices are completely fixed. The average duration of a price spell is $1/(1-\theta)$ periods.

From this micro structure, one can derive a forward-looking **New Keynesian Phillips Curve** in which current inflation depends on expected future inflation and current real marginal costs (proxied by the output gap), and a **New Keynesian IS curve** in which current output depends on expected future output and the real interest rate. Combined with a monetary policy rule (the Taylor rule), these three equations form the **three-equation New Keynesian model** — the dominant framework for monetary policy analysis in central banks worldwide. It is derived in full in Chapters 9, 10, and 23.

---

## 2.10 The Post-2008 Landscape and Open Questions

The Global Financial Crisis of 2007–09 exposed limitations in the pre-crisis consensus that were difficult to ignore. The dominant DSGE models of the era had no meaningful financial sector — they assumed complete asset markets and frictionless financial intermediation. They could not generate the kind of credit crunch, fire-sale dynamics, and balance-sheet recession that characterized the crisis and its aftermath. The slow recovery — output in many countries did not return to its pre-crisis trend for a decade — challenged models in which cyclical unemployment cannot permanently raise the natural rate.

Since 2010, the frontier has moved in several directions simultaneously. Models with **financial frictions** (Bernanke, Gertler, and Gilchrist, 1999; Gertler and Karadi, 2011) endogenize the cost of external finance and allow financial sector balance-sheet conditions to affect real investment. Models with **heterogeneous agents** (Aiyagari, 1994; Kaplan, Moll, and Violante, 2018) abandon the representative household and instead track the distribution of wealth across many households, allowing the aggregate consumption response to policy to depend on who is receiving income and who is liquidity-constrained. Models with **occasionally binding constraints** — in particular the zero lower bound on nominal interest rates — account for the large macroeconomic effects of the ELB that standard linearized models missed.

These developments are not the end of macroeconomic thought. They introduce new tractability challenges, new data requirements, and new debates about what the key mechanisms are. Understanding the current frontier requires first mastering the foundations that the frontier is extending, which is the purpose of the chapters that follow.

---

*Next: Chapter 3 — Key Concepts: GDP, Inflation, and Unemployment*
