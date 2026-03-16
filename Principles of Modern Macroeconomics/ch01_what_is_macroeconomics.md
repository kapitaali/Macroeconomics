# Chapter 1 — What Is Macroeconomics? The Big Picture and Why It Matters

> *"Economics is the science which studies human behaviour as a relationship between ends and scarce means which have alternative uses."*
> — Lionel Robbins, *An Essay on the Nature and Significance of Economic Science*, 1932

---

Why do some countries grow rich while others stagnate? What causes recessions, and can governments do anything to make them shorter? When a central bank raises interest rates, what exactly happens to employment, inflation, and the exchange rate — and over what time horizon? These questions matter not only to academic economists but to every person whose employment, savings, and standard of living depend on the answers. Macroeconomics is the branch of economics that attempts to answer them.

The word "macro" comes from the Greek for "large," and macroeconomics is distinguished from microeconomics precisely by its focus on the economy as a whole rather than on the decisions of individual agents. A microeconomist might ask why a particular firm chooses to hire fewer workers when the minimum wage rises, or why a specific consumer buys more wine when her income increases. A macroeconomist asks something different: what determines the *total* level of employment across all firms in an economy, or why household *aggregate* consumption rises in some years and falls in others. The shift from the individual to the aggregate is not merely a change in scale. It introduces entirely new phenomena — coordination failures, self-fulfilling expectations, monetary transmission, business cycles — that have no clean analogue at the level of a single household or firm.

---

## 1.1 The Scope and Method of Macroeconomic Analysis

The discipline organizes itself around three broad classes of question, each operating over a different time horizon.

The first is **business cycle analysis**: the study of short-run fluctuations in economic activity around its underlying trend. Economies do not grow smoothly. They expand rapidly in some periods, contract in others, and oscillate in a roughly recurrent pattern that economists call the business cycle. A business cycle, to be precise, is not a fixed-period oscillation like a pendulum; it is an irregular fluctuation in aggregate output, employment, and spending around a longer-run trend, driven by a mix of demand shocks, supply disturbances, and policy responses. Understanding what initiates these fluctuations, what propagates them across sectors and countries, and whether policy interventions can shorten recessions without creating other distortions is the central preoccupation of short-run macroeconomics.

The second is **long-run growth**: the study of what determines the level and trend rate of growth of a country's productive capacity over decades and generations. This is arguably the most consequential question in all of social science. A country that grows at 2% per year for a century ends up with output per person more than seven times larger than at the start; a country that grows at 1% ends up only 2.7 times richer. The difference between these trajectories — compounded over generations — dwarfs the effects of any plausible business cycle stabilization policy. Understanding why some countries sustain rapid growth for decades while others remain trapped at low income levels is therefore a question with profound implications for human welfare, and one to which economics has only partial answers.

The third is **stabilization policy**: the analysis of whether and how governments — through their control of public spending, taxation, and the money supply — can improve macroeconomic outcomes relative to a hypothetical unmanaged baseline. This question sits at the intersection of positive economics (what would actually happen if we changed policy?) and normative economics (what outcome do we want?), and it has been the most contested in the history of economic thought. Reasonable economists disagree sharply about the size of fiscal multipliers, the optimal inflation target, the appropriate response to a financial crisis, and the long-run effects of public debt. These disagreements reflect partly different readings of the evidence and partly different values, and it is important to keep both sources of disagreement visible rather than collapsing them into a single "the economists disagree" narrative.

### Formal Representation

Macroeconomic models are the tools through which these questions are addressed. A model is a deliberately simplified mathematical representation of the economy, specifying who the relevant agents are, what they are trying to achieve, what constraints they face, and how their individual decisions aggregate into economy-wide outcomes. The act of writing down a model forces the analyst to be explicit about every assumption. This makes those assumptions visible and therefore contestable — a virtue, not a weakness. A model that cannot be written down formally is too vague to generate testable predictions and too imprecise to communicate across researchers.

**Definition (Macroeconomic Model).** A macroeconomic model consists of: (i) a set of agents indexed by type — typically households, firms, a government, and a central bank; (ii) an objective function for each agent, such as a household's lifetime utility or a firm's present-discounted profits; (iii) a set of constraints for each agent, including budget constraints, production technology, and information sets; and (iv) an equilibrium concept that specifies how individual decisions are reconciled into consistent aggregate outcomes.

The standard equilibrium concept in modern macroeconomics is the **competitive equilibrium with rational expectations**.

**Definition (Macroeconomic Equilibrium).** A macroeconomic equilibrium is a tuple

$$\Bigl(\{c_{t}^{\ast} \, \ell_{t}^{\ast}\}_{t=0}^{\infty},\quad \{k_{t}^{*} \, n_{t}^{\ast}\}_{t=0}^{\infty},\quad \{w_{t}^{\ast} \, r_{t}^{\ast} \, p_{t}^{\ast}\}_{t=0}^{\infty}\Bigr)$$

of household plans for consumption $c_{t}^{*}$ and leisure $\ell_{t}^{*}$, firm plans for capital $k_{t}^{*}$ and labor demand $n_{t}^{*}$, and sequences of factor prices and the price level, such that:

1. Each household's plan maximizes lifetime utility subject to its intertemporal budget constraint, taking prices as given.
2. Each firm's plan maximizes the present discounted value of profits subject to its production technology and the capital accumulation equation, taking prices as given.
3. All markets clear at every date: goods ($c_{t}^{*} + i_{t}^{*} = y_{t}^{*}$), labor ($n_{t}^{s*} = n_{t}^{d*}$), and money (demand equals supply).
4. Expectations are **rational**: agents' beliefs about future variables are consistent with the distributions generated by the model itself.

Three features deserve comment. First, equilibrium is a *sequence* indexed by time, not a static vector of prices. Modern macroeconomics is inherently dynamic: today's consumption depends on expectations about tomorrow's income, and tomorrow's capital stock depends on today's investment. Second, equilibrium is *general*: every market clears simultaneously, so there is no partial equilibrium in which some prices are held fixed by assumption while others adjust. Third, the rational expectations requirement — condition (4) above — is the most substantive and contested element. It rules out systematic, predictable forecast errors: agents may be surprised, but they are not persistently wrong in the same direction. We examine this assumption critically in Chapters 15 and 16.

---

## 1.2 Positive Versus Normative Analysis

Before developing formal models it is worth distinguishing two types of economic claim, because confusion between them underlies many public debates about macroeconomic policy.

A **positive statement** is a claim about what is or what would happen under specific conditions. It is, in principle, testable against evidence. "When the central bank raises the nominal interest rate by one percentage point, output falls by approximately half a percent over the following year" is a positive claim. Whether it is true or false depends on the data, the estimation method, and the time period studied — but it is the kind of claim that evidence can inform. A **normative statement** is a claim about what ought to be. It expresses a value judgment that cannot be resolved by data alone. "The central bank should prioritize full employment over price stability" is normative: it depends on how one weighs the welfare costs of unemployment against the welfare costs of inflation, which is a choice, not a fact.

Positive macroeconomics attempts to explain observed outcomes: why did unemployment in the United States reach 25% in 1933 and 10% in 2009? Why did Japan experience two decades of near-zero inflation despite massive monetary expansion? Why do some countries sustain 7% annual growth for decades while others stagnate? Normative macroeconomics uses positive models to evaluate policy alternatives and recommend the best one, measured against an explicitly stated welfare criterion.

The distinction is important but porous in practice. Every positive model embeds normative choices: which variables to include, which welfare metric to maximize, which distributional weights to assign to rich and poor households, how to discount the welfare of future generations relative to the present. The standard welfare criterion — social welfare as the weighted sum of agents' lifetime utilities — is one choice among many, and its implicit value judgments are far from obvious. A model that maximizes aggregate consumption per capita treats a dollar gained by a billionaire as equivalent to a dollar gained by someone in poverty; whether this is the right criterion is a normative question that economics alone cannot answer. We return to these issues in Chapter 38 when discussing the distributional effects of macroeconomic policy.

---

## 1.3 Aggregate Variables and Their Measurement

A handful of economy-wide variables constitute the primary observables of macroeconomics. Understanding precisely what these variables measure — and what they do not — is a prerequisite for interpreting any empirical claim about the economy. We define each in turn and flag the most important measurement caveats.

**Definition (Real Gross Domestic Product).** Real gross domestic product $Y_t$ is the total market value of all *final* goods and services produced within a country's geographic borders during period $t$, measured at constant (base-year) prices. The qualification "final" excludes intermediate goods — inputs into further production — to avoid double-counting. The qualification "constant prices" removes the effect of changes in the general price level, so that changes in $Y_t$ reflect changes in physical quantities produced rather than mere inflation or deflation. Real GDP is the primary summary statistic for the level of economic activity.

**Definition (Inflation Rate).** Let $P_t$ denote the aggregate price level at date $t$ — a weighted average of the prices of all goods and services in the economy, constructed according to one of several index number formulas discussed in Chapter 3. The **inflation rate** is the proportionate change in the price level over one period:

$$\pi_t = \frac{P_t - P_{t-1}}{P_{t-1}}.$$

A positive $\pi_t$ indicates that the price level is rising (inflation); a negative $\pi_t$ indicates it is falling (deflation). The verbal definition of inflation is therefore: the rate at which the purchasing power of money is eroding over time.

**Definition (Labor Force, Employment, and the Unemployment Rate).** The **labor force** $L_t$ consists of all persons who are either employed or actively seeking employment. **Employment** $N_t$ is the number of persons in paid work. The number of **unemployed** persons is $U_t = L_t - N_t$ — those who are in the labor force but not currently working. The **unemployment rate** is the share of the labor force that is unemployed:

$$u_t = \frac{U_t}{L_t} = \frac{L_t - N_t}{L_t}.$$

One important subtlety: persons who have given up searching for work — so-called discouraged workers — are classified as out of the labor force rather than unemployed, which means $u_t$ understates the true extent of labor market slack. This point, and the broader question of which unemployment measure is most informative, is taken up in Chapter 3.

These three variables — $Y_t$, $\pi_t$, and $u_t$ — together with the nominal interest rate $i_t$, the exchange rate $e_t$, and the current account balance $CA_t$, constitute the core observables of macroeconomics. The study of how they are jointly determined, how they respond to shocks, and how policy can influence them is the subject of the remaining forty-one chapters.

| Variable | Symbol | Verbal Definition | Primary Data Source |
|---|---|---|---|
| Real GDP | $Y_t$ | Total real output produced domestically | BEA (US), ONS (UK), Eurostat |
| Employment | $N_t$ | Persons in paid work | BLS, Eurostat |
| Price level | $P_t$ | Weighted average of all prices | BLS, OECD |
| Inflation rate | $\pi_t$ | Rate of change of the price level | Derived: $(P_t - P_{t-1})/P_{t-1}$ |
| Unemployment rate | $u_t$ | Share of labor force not employed | BLS U-3, ILO |
| Nominal interest rate | $i_t$ | Cost of borrowing nominal funds | Federal Reserve, ECB |
| Exchange rate | $e_t$ | Price of foreign currency in domestic units | BIS, IMF |

---

## 1.4 Time Horizons: Why the Same Economy Has Different Models

One of the most initially confusing features of macroeconomics is that the same economy is analyzed using different models depending on the time horizon of interest. A student reading an introductory textbook learns that fiscal policy raises output; reading a graduate text she learns that in the long run money is neutral and output is determined entirely by supply-side factors. These statements are not contradictory — they apply to different time horizons — but the boundary between the horizons is not always crisp, and navigating between models requires care.

The fundamental source of the horizon dependence is the speed of **price adjustment**. Prices and wages in real economies do not jump instantaneously to market-clearing levels when demand conditions change. Firms have long-term customer relationships that make frequent price changes costly; workers and employers sign multi-year wage contracts; menu costs make it expensive for retailers to reprice constantly. These frictions mean that in the short run, when demand falls, firms cut output and employment rather than prices. Nominal rigidity — the stickiness of prices and wages in nominal (money) terms — is therefore the key distinction between short-run and long-run macroeconomics.

**Definition (Nominal Rigidity).** A nominal rigidity is any friction that prevents nominal prices or wages from adjusting instantaneously to changes in demand or money supply. In the presence of nominal rigidities, changes in the money supply or aggregate demand can have real effects — effects on employment, output, and real wages — in the short run. In the long run, as contracts expire and prices adjust, nominal rigidities disappear and money is neutral.

With this in mind, we can define the three time horizons precisely.

**Definition (Short Run).** The short run is the horizon over which nominal prices and wages are approximately fixed. On this horizon, fluctuations in aggregate demand translate primarily into fluctuations in real output and employment rather than in the price level. The short-run models of Parts II and IV are built around this assumption.

**Definition (Medium Run).** The medium run is the horizon over which price and wage adjustment is complete enough that the economy returns to its natural rate of output and unemployment. The **natural rate of output** $\bar{Y}_t$ is the level of real GDP consistent with a stable, non-accelerating inflation rate — the output the economy would produce if all prices and wages were fully flexible. It depends on factor endowments, technology, and institutions, but not on the price level or the money supply.

**Definition (Long Run).** The long run is the horizon over which the capital stock, the size and quality of the labor force, and the state of technology — the supply-side fundamentals — change substantially. On this horizon, the natural rate of output $\bar{Y}_t$ itself evolves, driven by capital accumulation and technological progress. Long-run growth theory, developed in Chapter 5, studies the forces governing this evolution.

The three time horizons are not separated by bright lines. A shock that generates a two-year recession has medium-run effects if it permanently reduces the labor force through hysteresis (a concept defined in Chapter 31) or long-run effects if it destroys productive capacity. Keeping track of which horizon a model applies to is one of the most important habits of mind a student of macroeconomics can develop.

---

## 1.5 The Interaction of Theory and Evidence

Macroeconomics aspires to be an empirical discipline. Its models generate predictions about the joint behavior of observable variables — how output, inflation, and employment should respond to identified shocks, and what statistical patterns (means, variances, autocorrelations) should appear in the data. These predictions can be confronted with evidence using econometric methods described in Appendix B.

The complication is that macroeconomic models are never tested in isolation. A test of whether fiscal multipliers exceed one is simultaneously a test of the consumption function, the investment function, the monetary policy rule, the expectations formation mechanism, and the degree of price flexibility. When the test fails, the data reveal that at least one of these components is misspecified, but they do not identify which one. Distinguishing between competing structural explanations for the same set of reduced-form correlations is the central inferential challenge of empirical macroeconomics.

This **joint hypothesis problem** is compounded by what is arguably the most important methodological insight of twentieth-century macroeconomics: the Lucas critique.

**Definition (The Lucas Critique).** A parameter is **structural** if it describes preferences or technology and is invariant to changes in the policy environment. A parameter is **reduced-form** if it is a function of both structural parameters and parameters describing the policy rule; it therefore changes whenever the policy rule changes. The Lucas critique (Lucas, 1976) is the observation that estimated reduced-form relationships between policy instruments and economic outcomes cannot be used to predict the effects of changes in policy rules, because such changes alter agents' decision rules and therefore the reduced-form parameters themselves. Only structural models — models built from preference and technology parameters — are immune to this critique.

To see why the Lucas critique matters concretely, consider the estimated relationship between money growth and inflation. Suppose a central bank observes a stable historical correlation $\hat{\pi} = a + b\hat{m}$ and uses this to predict that a 5-percentage-point increase in money growth will raise inflation by $5b$ percentage points. If the public believes this policy change is permanent and adjusts expectations accordingly, the actual inflation response will differ from $5b$ — because the inflation expectations term, omitted from the reduced-form, will have changed. The Lucas critique says that any model that does not explicitly model expectations formation is potentially vulnerable to this problem. This is why the DSGE methodology, which builds from structural first principles, is the standard for policy analysis in central banks and international financial institutions today.

---

## 1.6 A Map of the Book

Having established the scope of the discipline and its central methodological commitments, we can sketch the intellectual structure of what follows. Each part builds on the previous, and the cross-references are intentional: macroeconomics is a cumulative discipline in which later models do not replace earlier ones but rather expose their limitations and extend them.

Part I lays the foundations: intellectual history (Chapter 2), measurement of core variables (Chapter 3), national income accounting (Chapter 4), growth theory (Chapter 5), and the data infrastructure of empirical macroeconomics (Chapter 6).

Part II develops the core workhorse short-run models in roughly historical order: AS–AD (Chapter 7), the Keynesian cross (Chapter 8), IS–LM (Chapter 9), and the Phillips curve (Chapter 10). Each is presented with its original motivation and then evaluated against the limitations that prompted the next development.

Part III grounds the aggregate relationships of Part II in household and firm optimization: consumption (Chapter 11), investment (Chapter 12), labor supply and demand (Chapter 13), money demand (Chapter 14), expectations (Chapter 15), and rational expectations and the new classical school (Chapter 16). The consistent theme is that aggregate relationships are only as reliable as the microfoundations underlying them.

Part IV analyzes equilibrium in the major markets separately and together: goods (Chapter 17), money (Chapter 18), labor (Chapter 19), financial markets (Chapter 20), and the open economy (Chapter 21).

Part V disaggregates by institutional sector: the government (Chapter 22), the central bank (Chapter 23), the corporate sector (Chapter 24), households and demographics (Chapter 25), and the foreign sector (Chapter 26).

Part VI applies the framework to concrete policy questions: business cycles (Chapter 27), fiscal policy in practice (Chapter 28), monetary policy in practice (Chapter 29), inflation and deflation (Chapter 30), and unemployment (Chapter 31).

Part VII covers the domain-specific literatures: open-economy macroeconomics and global imbalances (Chapter 32), economic development (Chapter 33), financial crises and regulation (Chapter 34), and macroeconomic policy in developing countries (Chapter 35).

Part VIII surveys the frontier: the digital economy and productivity (Chapter 36), climate change (Chapter 37), inequality (Chapter 38), and the future of the discipline (Chapter 39).

Part IX concludes with two extended case studies — the Great Recession of 2008 (Chapter 40) and the COVID-19 pandemic (Chapter 41) — followed by a chapter on applying macroeconomic reasoning in professional practice (Chapter 42).

---

*Next: Chapter 2 — A Brief History of Macroeconomic Thought*
