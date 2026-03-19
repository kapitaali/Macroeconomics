# Chapter 5 — Economic Growth and Development: The Long-Run Perspective

> *"The consequences for human welfare involved in questions like these are simply staggering: Once one starts to think about them, it is hard to think about anything else."*
> — Robert E. Lucas Jr., *On the Mechanics of Economic Development*, 1988

---

The previous chapters addressed how we measure economic activity in the short and medium run. This chapter steps back and asks a more fundamental question: what determines the *level* of economic activity in the long run — the standard of living that a country can sustain not just this year but over decades and generations?

This is not merely an academic exercise. The gap in living standards between the world's richest and poorest countries is staggering: real income per person in the United States is roughly fifty times that of Malawi, sixty times that of Niger. A person born in Norway today can expect to live more than thirty years longer than one born in the Central African Republic. These differences are not explained by short-run fluctuations in aggregate demand; they reflect profound, persistent differences in productive capacity — in the amount of physical capital workers have at their disposal, the level of education and health they bring to their work, and the state of the technology they use. Growth theory attempts to explain why these productive capacities differ so dramatically and what, if anything, can change them.

---

## 5.1 What Growth Theory Is Explaining

Before presenting models, it is worth being precise about what we are trying to explain. The relevant outcome variable is **real output per worker**, or equivalently real output per capita in an economy where the employment rate is approximately constant. We care about output per worker rather than total output because a country is not richer simply by having more people — what matters for living standards is how much each person can produce.

Over the very long run, output per worker in currently advanced economies has grown at roughly 2% per year. This sounds modest, but compound growth at 2% per year means a doubling of living standards every 35 years. The standard of living of a typical American in 2000 was roughly eight times higher than in 1900, and roughly sixty times higher than in 1800. This sustained improvement is the central fact that any growth theory must explain.

Alongside this long-run trend, several robust regularities appear when comparing growth experiences across countries. These are the **Kaldor stylized facts** (Kaldor, 1961):

1. Output per worker grows at a roughly constant rate over long periods.
2. The capital-to-labor ratio grows at a roughly constant rate.
3. The real interest rate (return on capital) is roughly constant over long periods.
4. The capital-to-output ratio is roughly constant.
5. Labor and capital shares of national income are roughly constant.
6. Growth rates of output per worker differ substantially across countries.

Facts 1–5 together describe what growth theorists call a **balanced growth path** — a trajectory on which all quantities grow at constant rates and all ratios remain approximately stable. Fact 6 is the puzzle: if balanced growth is the long-run attractor, why do some countries grow along a high-income balanced path while others are stuck at a low-income one?

**Definition (Balanced Growth Path).** A balanced growth path is a trajectory of the economy on which output per effective worker $\tilde{y}_t \equiv Y_t/(A_t L_t)$, capital per effective worker $\tilde{k}_t \equiv K_t/(A_t L_t)$, and consumption per effective worker $\tilde{c}_t$ are all constant. Here $A_t$ is the level of labor-augmenting technology and $L_t$ the size of the labor force. Output per actual worker $y_t = Y_t/L_t$ grows along the balanced path at the rate $g$ of technological progress, consistent with Kaldor facts 1 and 2.

---

## 5.2 The Solow–Swan Growth Model

The first rigorous model of long-run growth is due independently to Robert Solow (1956) and Trevor Swan (1956). Its central insight is deceptively simple: capital accumulation alone cannot sustain indefinite growth, because of diminishing returns; sustained growth in output per worker requires sustained improvements in technology.

To see why, consider the production function. The economy produces a single good using capital $K_t$ and labor $L_t$ augmented by the level of technology $A_t$:

$$Y_t = F(K_t,\, A_t L_t).$$

The function $F$ is assumed to exhibit **constant returns to scale**: doubling both inputs doubles output, $F(\lambda K, \lambda AL) = \lambda F(K, AL)$ for any $\lambda > 0$. It also exhibits **diminishing marginal returns** to each input: $F_K > 0$, $F_{KK} < 0$, meaning each additional unit of capital adds positive but declining amounts to output. The **Inada conditions** — $\lim_{K \to 0} F_K = \infty$ and $\lim_{K \to \infty} F_K = 0$ — ensure that production cannot take place with zero capital and that capital is never so productive that the economy would want to accumulate it without limit.

The reason technology enters as $A_t L_t$ — multiplying labor rather than capital — is that this specification of **labor-augmenting** (or Harrod-neutral) technological progress is the only one consistent with a balanced growth path when preferences are of the constant-relative-risk-aversion form. Technology augmenting capital would violate the constancy of factor shares (Kaldor fact 5); technology entering in a Hicks-neutral form $A_t F(K_t, L_t)$ also creates difficulties. This is not a minor technical point: the form of technological progress determines whether a balanced growth path exists, and it is a constraint on model-building with empirical content.

**Definition (Effective Labor and Intensive Form).** Define **effective labor** as $A_t L_t$ — the product of the technology level and the raw number of workers. Define capital per effective worker $\tilde{k}_t \equiv K_t/(A_t L_t)$ and output per effective worker $\tilde{y}_t \equiv Y_t/(A_t L_t)$. Dividing both sides of the production function by $A_t L_t$, the **intensive form** of production is:

$$\tilde{y}_t = f(\tilde{k}_t), \quad \text{where } f(\tilde{k}) \equiv F(\tilde{k}, 1).$$

Working in intensive form collapses a model with two growing variables ($K_t$ and $A_t L_t$) into a model with one stationary variable $\tilde{k}_t$, which is much easier to analyze.

The economy's saving rate $s \in (0,1)$ is, in the basic Solow model, taken as exogenous — households save a fixed fraction $s$ of their income regardless of the interest rate or their age. This is the model's key simplifying assumption; it will be relaxed in the Ramsey–Cass–Koopmans model of Section 5.3. With exogenous saving, the capital accumulation equation in intensive form is:

$$\dot{\tilde{k}}_t = s\, f(\tilde{k}_t) - (n + g + \delta)\tilde{k}_t,$$

where $n \equiv \dot{L}_t/L_t$ is the population growth rate, $g \equiv \dot{A}_t/A_t$ is the rate of technological progress, and $\delta \in (0,1)$ is the depreciation rate of capital. The term $(n + g + \delta)\tilde{k}_t$ is the **effective depreciation**: even with zero physical depreciation, $\tilde{k}_t$ would shrink if population grows (each worker has less capital) or if technology improves (each unit of effective labor is endowed with less capital). The effective depreciation rate $\mu \equiv n + g + \delta$ plays an important role throughout.

The intuition for this equation is worth spelling out. The change in capital per effective worker equals investment per effective worker ($s\, f(\tilde{k}_t)$, since all saving is invested in a closed economy) minus the amount of investment needed just to keep $\tilde{k}_t$ constant as the labor force grows and technology improves. If investment exceeds break-even, $\tilde{k}_t$ rises; if it falls short, $\tilde{k}_t$ falls.

### The Steady State

The **steady state** of the Solow model is the value $\tilde{k}^*$ at which capital per effective worker is constant: $\dot{\tilde{k}} = 0$. Setting the right-hand side to zero:

$$s\, f(\tilde{k}^*) = \mu\, \tilde{k}^*.$$

This says that investment per effective worker exactly covers effective depreciation. The steady state exists and is unique under the Inada conditions, because $sf(\tilde{k})/\tilde{k}$ — the average product of capital, scaled by $s$ — is declining from infinity (as $\tilde{k} \to 0$) to zero (as $\tilde{k} \to \infty$), and it crosses the horizontal line $\mu$ exactly once.

With the **Cobb–Douglas production function** $f(\tilde{k}) = \tilde{k}^\alpha$, $\alpha \in (0,1)$, the steady state has a closed form:

$$\tilde{k}^* = \left(\frac{s}{\mu}\right)^{\frac{1}{1-\alpha}}, \qquad \tilde{y}^* = \left(\frac{s}{\mu}\right)^{\frac{\alpha}{1-\alpha}}.$$

These formulas have important implications. A higher saving rate $s$ raises $\tilde{k}^*$ and $\tilde{y}^*$ — countries that invest more are richer in the long run. A higher population growth rate $n$ or depreciation rate $\delta$ lowers $\tilde{k}^*$ and $\tilde{y}^*$ — these forces continuously dilute the capital stock. The parameter $\alpha$ (the capital share in income, approximately one-third in most countries) governs the sensitivity of income to saving: because $\alpha/(1-\alpha) \approx 0.5$, doubling the saving rate raises steady-state income by only about 41%, not 100%.

### Transitional Dynamics and Conditional Convergence

Away from the steady state, $\tilde{k}_t$ converges to $\tilde{k}^*$ monotonically. Linearizing around the steady state:

$$\frac{\mathrm{d}}{\mathrm{d}t}(\tilde{k}_t - \tilde{k}^*) \approx -\lambda(\tilde{k}_t - \tilde{k}^*), \qquad \lambda = (1-\alpha)\mu.$$

The solution is $\tilde{k}_t - \tilde{k}^* = (\tilde{k}_0 - \tilde{k}^*)\,e^{-\lambda t}$: the gap between current capital and the steady state shrinks exponentially at rate $\lambda$. With $\alpha = 1/3$ and $\mu \approx 0.06$, $\lambda \approx 0.04$, giving a half-life of $\ln 2 / \lambda \approx 17$ years. This means the Solow model predicts that economies converge to their steady states fairly slowly — a prediction consistent with observed growth patterns.

This convergence is **conditional**: countries converge to their own steady states, not to a common one. A poor country that also has a low saving rate and high population growth will have a low $\tilde{k}^*$ and will not catch up to a rich country with high $s$ and low $n$. The **conditional convergence hypothesis** — that countries converge in growth rates after conditioning on the determinants of the steady state — is strongly supported empirically (Mankiw, Romer, and Weil, 1992). **Unconditional convergence** — the idea that all poor countries simply grow faster than rich ones — is not supported in the broad cross-section of countries, though it holds within groups of similarly structured economies such as the OECD.

### The Golden Rule of Capital Accumulation

A natural question is: which saving rate $s$ maximizes household welfare? Consumption per effective worker on the balanced growth path is $\tilde{c}^* = (1-s)f(\tilde{k}^*)$. Since $\tilde{k}^*$ is an increasing function of $s$, raising $s$ simultaneously increases output ($f(\tilde{k}^*)$ rises) and decreases the fraction of output available for consumption ($1-s$ falls). The saving rate that maximizes steady-state consumption — the **Golden Rule saving rate** — balances these two forces:

$$f'(\tilde{k}^{GR}) = \mu = n + g + \delta.$$

Verbally: at the Golden Rule, the marginal product of capital equals the effective depreciation rate. Since the net marginal product of capital $r^{GR} = f'(\tilde{k}^{GR}) - \delta$, the Golden Rule condition is equivalent to $r^{GR} = n + g$.

An economy is **dynamically efficient** if $r > \delta + n + g$, meaning the return on capital exceeds what is needed to maintain $\tilde{k}^*$ at the Golden Rule level. Dynamically efficient economies have not overaccumulated capital — households are not saving too much relative to the welfare-maximizing amount. An economy with $r < \delta + n + g$ is **dynamically inefficient**: it has accumulated so much capital that reducing investment would allow more consumption at every future date, benefiting every generation. Most empirical evidence suggests that advanced economies are dynamically efficient (Abel et al., 1989).

---

## 5.3 The Ramsey–Cass–Koopmans Model

The Solow model's exogenous saving rate is its principal limitation. Households in the Solow model behave as if they follow a mechanical saving rule rather than choosing optimally. This is unsatisfying both intellectually — we should derive saving behavior from preferences — and practically, because it prevents us from analyzing how saving responds to changes in policy, taxation, or the return on capital.

The Ramsey–Cass–Koopmans (RCK) model endogenizes saving by replacing the mechanical saving rule with explicit intertemporal optimization by an infinitely-lived representative household. The "infinitely lived" assumption is a mathematical convenience: it can be justified either literally (people care about their descendants as much as themselves) or as an approximation to an overlapping-generations model where bequests are operative.

The representative household has $L_t = L_0 e^{nt}$ members and maximizes the present discounted value of utility:

$$U = \int_0^\infty e^{-(\rho - n)t}\, u(c_t)\, \mathrm{d}t, \qquad \rho > n,$$

where $\rho > 0$ is the **pure rate of time preference** — the rate at which the household discounts future utility relative to present utility, holding consumption constant. The requirement $\rho > n$ ensures that the utility integral converges. The per-period utility function $u(c_t)$ satisfies the standard conditions: $u' > 0$ (more consumption is better), $u'' < 0$ (diminishing marginal utility), and the Inada conditions $\lim_{c\to 0} u'(c) = \infty$, $\lim_{c\to\infty} u'(c) = 0$. The standard parametric form is the **constant relative risk aversion (CRRA) utility function**:

**Definition (CRRA Utility).** The CRRA utility function is:

$$u(c) = \frac{c^{1-\sigma} - 1}{1-\sigma}, \qquad \sigma > 0,\; \sigma \neq 1,$$

with the limiting case $u(c) = \ln c$ when $\sigma = 1$. The parameter $\sigma$ is simultaneously the **coefficient of relative risk aversion** (how much a household dislikes gambles over consumption) and the inverse of the **elasticity of intertemporal substitution** (EIS) $\varepsilon = 1/\sigma$. The EIS measures how responsive the household's consumption growth is to the real interest rate: a high $\sigma$ (low EIS) means households are reluctant to substitute consumption across time, so even a high interest rate induces little additional saving. The CRRA family is the only class of utility functions consistent with balanced growth (King, Plosser, and Rebelo, 1988).

Writing the household's problem in intensive form with $\tilde{c}_t = c_t/A_t$, the economy's resource constraint is:

$$\dot{\tilde{k}}_t = f(\tilde{k}_t) - \tilde{c}_t - (n + g + \delta)\tilde{k}_t,$$

the same as in the Solow model but with $\tilde{c}_t$ replacing $s\,f(\tilde{k}_t)$. The household now chooses the entire path $\{\tilde{c}_t\}_{t \geq 0}$ to maximize utility subject to this constraint. This is an **optimal control problem**: the state variable is $\tilde{k}_t$ (capital accumulated from past decisions) and the control variable is $\tilde{c}_t$ (the decision made at each instant).

The **current-value Hamiltonian** is:

$$\mathcal{H} = \frac{(A_t\tilde{c}_t)^{1-\sigma}-1}{1-\sigma} + \mu_t\bigl[f(\tilde{k}_t) - \tilde{c}_t - (n+g+\delta)\tilde{k}_t\bigr],$$

where $\mu_t$ is the **costate variable** — the shadow price of capital in utility units, measuring how much lifetime utility increases if the capital stock is marginally increased.

The **Pontryagin maximum principle** provides three conditions for optimality. The first-order condition for $\tilde{c}_t$:

$$\frac{\partial \mathcal{H}}{\partial \tilde{c}_t} = 0 \implies A_t^{1-\sigma}\tilde{c}_t^{-\sigma} = \mu_t.$$

The costate equation:

$$\dot{\mu}_t = (\rho - n)\mu_t - \frac{\partial \mathcal{H}}{\partial \tilde{k}_t} = (\rho - n)\mu_t - \mu_t f'(\tilde{k}_t).$$

The **transversality condition**:

$$\lim_{t \to \infty} e^{-(\rho-n)t}\,\mu_t\,\tilde{k}_t = 0.$$

The transversality condition rules out "Ponzi schemes" in which the household accumulates capital without limit, intending to consume it all at infinite time. It says that the present value of the capital stock must approach zero as the horizon extends to infinity — the household cannot plan to die with positive assets (that would mean it consumed too little) or negative assets (that would mean it consumed on borrowed time forever).

Differentiating the first-order condition with respect to time and substituting the costate equation yields the **Euler equation for consumption**:

$$\frac{\dot{\tilde{c}}_t}{\tilde{c}_t} = \frac{1}{\sigma}\bigl[f'(\tilde{k}_t) - \delta - \rho - \sigma g\bigr].$$

In per-capita terms, this simplifies to:

$$\frac{\dot{c}_t}{c_t} = \frac{r_t - \rho}{\sigma},$$

where $r_t = f'(\tilde{k}_t) - \delta$ is the **net marginal product of capital** — the real return on saving. This is the fundamental equation of consumption theory. It says: consumption grows when the real return on saving exceeds the rate of impatience. If $r_t > \rho$, it is worthwhile to sacrifice present consumption for future consumption; if $r_t < \rho$, the household prefers to consume now rather than wait.

### The Steady State and Phase-Plane Analysis

The steady state requires both $\dot{\tilde{c}} = 0$ and $\dot{\tilde{k}} = 0$. From the Euler equation, $\dot{\tilde{c}} = 0$ requires:

$$f'(\tilde{k}^*) = \delta + \rho + \sigma g.$$

This determines the steady-state capital stock independently of the saving rate — a sharp contrast with the Solow model, where the steady state depends on the exogenous $s$. In the RCK model, the steady-state capital stock is pinned down by preference parameters ($\rho$, $\sigma$) and technology ($g$, $\delta$). The corresponding steady-state consumption is determined by substituting $\tilde{k}^*$ into the resource constraint with $\dot{\tilde{k}} = 0$: $\tilde{c}^* = f(\tilde{k}^*) - \mu\tilde{k}^*$.

Comparing the RCK steady state with the Golden Rule: since $f'(\tilde{k}^*) = \delta + \rho + \sigma g > \delta + n + g$ (assuming $\rho > n + (1-\sigma)g$, which holds for standard parameter values), we have $\tilde{k}^* < \tilde{k}^{GR}$. The optimizing economy holds less capital than the Golden Rule would prescribe. Why? Because households are impatient ($\rho > 0$) and therefore discount future consumption. Accumulating more capital to reach the Golden Rule would require sacrificing present consumption, which households do not find worthwhile given their impatience. The RCK model thus endorses dynamic efficiency: rational households will never allow the economy to overaccumulate capital.

The dynamics of the two-dimensional system $(\tilde{k}_t, \tilde{c}_t)$ are analyzed using a phase diagram. The $\dot{\tilde{k}} = 0$ locus is a hump-shaped curve in $(\tilde{k}, \tilde{c})$ space (the resource constraint evaluated at $\dot{\tilde{k}} = 0$). The $\dot{\tilde{c}} = 0$ locus is a vertical line at $\tilde{k} = \tilde{k}^*$. The equilibrium $(\tilde{k}^*, \tilde{c}^*)$ is a **saddle point**: there are four regions of the phase plane corresponding to different combinations of rising/falling $\tilde{k}$ and $\tilde{c}$, and the stable manifold — called the **saddle path** — is the unique trajectory that converges to the steady state. Since $\tilde{c}_0$ is a free variable (the household can choose its initial consumption), the economy jumps immediately to the saddle path at $t=0$ and then travels along it to the steady state. This jump is not a discontinuity in behavior; it is the household choosing its initial consumption optimally, given that it knows the optimal path forward.

---

## 5.4 Exogenous Versus Endogenous Growth

In both the Solow and RCK models, the long-run growth rate of per-capita income equals $g$ — the rate of technological progress. This is **exogenous growth**: the engine of long-run improvement in living standards (technology) is taken as given by the model, falling like manna from heaven. The model explains the *level* of income but not the *source* of growth. This is a significant limitation. If we want to know whether subsidizing research, opening trade, or improving education will raise a country's long-run growth rate, a model with exogenous technology can only say "not if technology is exogenous" — which is circular.

**Definition (Exogenous Growth).** A growth model exhibits **exogenous growth** if the long-run growth rate of per-capita income is determined by variables — typically the rate of technological progress $g$ — that are taken as parameters outside the model. Policy can change the *level* of the balanced growth path but not the long-run *growth rate*.

**Definition (Endogenous Growth).** A growth model exhibits **endogenous growth** if the long-run growth rate of per-capita income is determined within the model — by the choices of households, firms, and the government regarding saving, education, research, and other investments. In endogenous growth models, policy can permanently raise the growth rate, not merely the level of income.

The distinction has profound policy implications. Under exogenous growth, a policy that permanently raises the saving rate shifts the balanced growth path upward but does not change its slope: the economy reaches a permanently higher level of income but eventually resumes growing at the same rate $g$. Under endogenous growth, a policy that permanently raises investment in knowledge or technology can permanently raise the growth rate — the slope of the income path changes. Whether this distinction corresponds to something empirically important is one of the central debates in growth economics.

### The AK Model

The simplest model of endogenous growth eliminates diminishing returns by assuming a linear production function:

$$Y_t = A K_t, \quad A > 0,$$

where $K_t$ is interpreted broadly to include both physical capital and human capital — knowledge, skills, and organizational capabilities that can be accumulated without limit. With saving rate $s$, capital accumulation is $\dot{K}_t = sY_t - \delta K_t = (sA - \delta)K_t$. If $sA > \delta$, output grows permanently at rate $g = sA - \delta$. Crucially, $g$ depends on the saving rate $s$: policies that raise saving raise the permanent growth rate. This is the defining feature of endogenous growth.

The AK model's linearity is, of course, a strong assumption. It requires that the external returns to capital (broadly defined) exactly offset the private diminishing returns, so that the social production function is linear. Whether this holds empirically — whether knowledge and human capital accumulation generate sufficient externalities to prevent long-run diminishing returns — is an empirical question on which the evidence is mixed.

### Romer's Knowledge Accumulation Model

A more satisfying approach to endogenous growth endogenizes the process of technological innovation itself. Romer (1990) models technology as the accumulated stock of designs — blueprints for new producer goods — produced by profit-seeking researchers. This makes technological progress the result of intentional investment decisions rather than an unexplained gift.

In the Romer model, the economy has three sectors. A **research sector** employs skilled labor $H_A$ to produce new designs; the flow of new designs is governed by the knowledge production function:

$$\dot{A}_t = \delta_A\, H_A\, A_t,$$

where $\delta_A > 0$ is a productivity parameter. Crucially, the stock of existing knowledge $A_t$ enters multiplicatively: standing on the shoulders of giants, today's researchers are more productive when more knowledge already exists. This generates the **standing-on-shoulders effect** that makes knowledge accumulation self-reinforcing.

A **goods sector** uses the existing stock of designs, along with human labor $H_Y$ and unskilled labor $L_Y$, to produce final output:

$$Y_t = H_Y^\alpha\, L_Y^{1-\alpha} \int_0^{A_t} x_i^\alpha\, \mathrm{d}i,$$

where $x_i$ is the quantity of intermediate input $i$ (embodied in design $i$). Each intermediate good is produced by a **monopolist** who holds the patent on the corresponding design — this monopoly profit is the reward that motivates researchers to invest in generating new designs.

In symmetric equilibrium, the long-run growth rate is $g = \delta_A H_A$, which is endogenous: it depends on the share of human capital devoted to research, $H_A/(H_A + H_Y)$, which is in turn affected by policies like R&D subsidies, intellectual property protection, and education spending. The model provides a rigorous foundation for policies aimed at stimulating innovation.

Jones (1995) identified a problem with the Romer model as stated: the **scale effect**. Because $g = \delta_A H_A$, doubling the population of researchers doubles the long-run growth rate. But over the past century, the number of scientists and engineers engaged in R&D in the United States has grown by more than tenfold, while the long-run growth rate of per-capita income has remained roughly constant at 2%. The scale effect prediction is decisively rejected. **Semi-endogenous growth models** (Jones, 1995b) resolve this by introducing diminishing returns in knowledge production: $\dot{A}_t = \delta_A H_A^\lambda A_t^\phi$ with $\phi < 1$. In the resulting model, the long-run growth rate is independent of the level of R&D investment (it depends only on population growth), but the transition path — and hence long-run income levels — depends on policy.

---

## 5.5 The Empirics of Convergence

The Solow model generates a specific, testable prediction about the growth rate of countries. Countries that are currently below their steady state should grow faster than countries that are at or above it. Converting the linearized dynamics to per-capita output growth rates gives the **convergence regression**:

$$\frac{\ln y_{iT} - \ln y_{i0}}{T} = \gamma - \frac{1 - e^{-\lambda T}}{T}\,\ln y_{i0} + \epsilon_{iT},$$

where the coefficient on initial income $\ln y_{i0}$ is negative: countries that start poor (low $y_{i0}$) grow faster, converging toward richer countries. The constant $\gamma$ absorbs the determinants of the steady state.

A crucial distinction governs the interpretation of this regression. **Unconditional convergence** — the hypothesis that poor countries simply grow faster regardless of their structural characteristics — would require the coefficient on $\ln y_{i0}$ to be negative in an unrestricted regression with no additional controls. This hypothesis is decisively rejected in cross-country data: poor countries are not, on average, growing faster than rich ones. Countries like Korea and Taiwan were poor in 1960 and have grown rapidly; countries like Zambia and Haiti were also poor in 1960 and have not grown at all. Unconditional convergence confuses the dynamics *within* an economy (toward its own steady state) with the dynamics *across* economies (toward a common steady state).

**Conditional convergence** — the hypothesis that countries converge to their own steady states, which differ across countries — is strongly supported. When the regression controls for the determinants of the steady state (saving rates, population growth, human capital), the coefficient on initial income becomes significantly negative. Mankiw, Romer, and Weil (1992) show that an augmented Solow model with human capital explains roughly 80% of the cross-country variation in per-capita income when these controls are included, and yields estimates of the physical and human capital shares ($\hat{\alpha} \approx 0.31$, $\hat{\beta} \approx 0.28$) consistent with factor income shares observed in national accounts.

The implication is that long-run differences in income are primarily explained by differences in steady-state characteristics — saving rates, human capital investment, population growth — rather than by idiosyncratic country-specific factors. Whether the determinants of these steady-state characteristics are themselves explainable (by institutions, geography, culture, historical accident) is the subject of Chapter 33.

---

*Next: Chapter 6 — Macroeconomic Data and Sources*
