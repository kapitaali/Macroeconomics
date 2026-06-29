# Chapter 36: Agent-Based Computational Macroeconomics

*Heterogeneous Agents, Bounded Rationality, and Emergent Cycles*

> *"The economy is not in equilibrium. It is always in the process of changing, adapting, and failing to adapt. Agent-based models are the natural language for this world."*
> — W. Brian Arthur

**Cross-reference:** *Principles* Ch. 15 (bounded rationality, animal spirits, non-ergodic systems); Ch. 39 (future of macroeconomics: ABMs, complexity economics) **[P:Ch.15, P:Ch.39]**

---

## 36.1 Agent-Based Models vs. DSGE: A Methodological Comparison

Throughout this book, macroeconomic models have been built from **representative optimization**: households maximize utility, firms maximize profits, markets clear instantaneously. This approach has two enormous advantages: internal consistency (every agent is doing the best they can given their information) and analytical tractability (the first-order conditions are algebraic equations that can be solved by the methods of Parts III–VII).

But two features of the real economy are absent from this picture:

**1. Genuine heterogeneity in behavior.** Not all households optimize the same objective function. Some are hand-to-mouth not because of borrowing constraints, but because they are inattentive. Firms use diverse pricing strategies — some markup, some match, some follow simple rules of thumb. The behavioral economics literature (Kahneman, Thaler, Laibson) documents systematic deviations from expected utility maximization.

**2. Disequilibrium dynamics.** Markets do not instantaneously clear. When Lehman failed, credit markets froze for months — a disequilibrium state that no standard equilibrium model can represent. Animal spirits, confidence, herding, and coordination failures create dynamics that are path-dependent and not ergodic.

**Agent-based models (ABMs)** address both by abandoning the optimization assumption and replacing it with behavioral rules. Agents follow simple, plausible rules (e.g., "raise my price by 5% if my inventories are low"); macro dynamics emerge from the aggregate of micro interactions. ABMs sacrifice internal consistency and analytical tractability in exchange for realism and flexibility.

---

## 36.2 The ABM Structure

**Definition 36.1 (Agent-Based Model).** An **agent-based model** consists of:

1. A finite set of **agents** $\mathcal{A} = \{1, \ldots, N\}$, each with a state $s_i(t) \in S_i$ (wealth, employment, inventory, beliefs, etc.).
2. An **environment** describing external conditions and market infrastructure.
3. **Interaction rules** specifying how agents interact with each other and the environment at each time step.
4. An **update rule**: at each discrete time step $t$, all agents (or a random subset) update their states according to the interaction rules.

**The Markov chain representation:** If each agent's state space $S_i$ is finite, the ABM defines a Markov chain on the aggregate state space $\mathcal{S} = \prod_i S_i$ (exponentially large but finite). The macro statistics of interest (GDP, unemployment, inflation) are functions of this aggregate state.

**Definition 36.2 (Ergodicity).** An ABM is **ergodic** if its Markov chain has a unique stationary distribution — time averages converge to ensemble averages. Non-ergodic ABMs (with multiple attractors, path-dependence, or absorbing states) require simulation for each specific history, not just for the stationary distribution.

---

## 36.3 Key ABM Results in Macroeconomics

### 36.3.1 The Schelling Segregation Model

The Schelling (1969, 1971) model demonstrates how **mild individual preferences for neighbors of the same type** generate strong aggregate segregation — a counterintuitive emergent phenomenon.

**Setup:** $N$ households on a grid, each colored red or blue. Each household prefers at least fraction $\tau \in (0, 1)$ of neighbors of the same color. If unsatisfied, it moves to a random empty location.

**Result:** For $\tau$ as low as $1/3$ (each agent tolerates 2 out of 3 opposite-color neighbors), the model converges to nearly complete segregation ($>80\%$ same-color neighbors).

**Relevance to macroeconomics:** Schelling demonstrates that aggregate outcomes can be qualitatively different from individual preferences. The same logic applies to bank runs (each depositor prefers not to panic but panics if others panic), speculative bubbles, and the global financial crisis (each institution's individually rational deleveraging amplifies system-wide stress).

### 36.3.2 The Axelrod Cooperation Model

Axelrod's (1984) computer tournament showed that **Tit-for-Tat** — start cooperating, copy the opponent's last move — was the most successful strategy in iterated Prisoner's Dilemma tournaments, even though it is never strictly dominant.

**Macroeconomic applications:** International trade agreements, central bank credibility, fiscal rules, and debt renegotiation can all be modeled as repeated games where reputation and reciprocity sustain cooperation that one-shot optimization cannot.

### 36.3.3 Dosi et al. (2010): The Macro ABM

The Dosi, Fagiolo, and Roventini (2010) K+S model (Keynes meets Schumpeter) is the most influential macro ABM. It features:

- Heterogeneous capital-goods firms investing in R&D to produce better machines.
- Heterogeneous consumption-goods firms buying machines, setting prices, and hiring workers.
- Households providing labor, consuming, and saving.
- A banking sector providing credit.

**Key results:**
- Business cycles emerge endogenously from the micro interactions, without external shocks.
- The model replicates 14 stylized facts of OECD business cycles (volatility, autocorrelation, cross-correlations) without calibrating to second moments.
- Austerity policies (fiscal consolidation during recessions) worsen output in the model — consistent with the post-2010 European experience.

---

## 36.4 A Minimal Macro ABM

We implement a simple ABM with three types of agents: households, firms, and a central bank.

**Household rules:**
- Consume fraction $c_{it}$ of their liquid wealth each period (heterogeneous across agents).
- Work if employed; search for work if unemployed.
- Update consumption propensity based on whether employed or unemployed.

**Firm rules:**
- Set prices using a markup rule: $p_{it} = (1+m_{it})MC_{it}$ where $m_{it}$ adjusts when inventories are high (cut markup) or low (raise markup).
- Hire workers when vacancies exceed threshold; fire when output > demand by threshold.
- Invest based on expected future demand.

**Central bank rule:**
- Taylor rule: $i_t = r^* + \phi_\pi(\pi_t - \pi^*) + \phi_y(y_t - y^*)$.

**Algorithm 36.1 (ABM Time Step).**

For each period $t$:
1. Central bank sets $i_t$ based on Taylor rule.
2. Households update employment status (job search with matching probability).
3. Firms update inventories, set prices, decide production.
4. Households consume and save.
5. Goods market: buyers and sellers randomly matched; transactions occur at posted prices.
6. Labor market: firms post vacancies; unmatched workers remain unemployed.
7. Update beliefs and expectations.
8. Collect aggregate statistics: $Y_t$, $U_t$, $\pi_t$.

```python
import numpy as np
import matplotlib.pyplot as plt

class MacroABM:
    """Minimal macroeconomic agent-based model."""
    
    def __init__(self, N_hh=500, N_firms=50, T=300, seed=42):
        np.random.seed(seed)
        self.N_hh = N_hh; self.N_firms = N_firms; self.T = T
        
        # Household initial states
        self.wealth = np.random.exponential(10, N_hh)
        self.employed = np.random.rand(N_hh) < 0.95  # 95% initial employment
        # Heterogeneous MPC: higher MPC for poorer households
        self.mpc = np.clip(1 - 0.5*self.wealth/self.wealth.mean(), 0.3, 0.95)
        
        # Firm initial states
        self.prices = np.ones(N_firms) * 1.0
        self.markup  = np.random.uniform(0.1, 0.3, N_firms)
        self.inventories = np.ones(N_firms) * 20.0
        self.workers = np.maximum(1, (N_hh//N_firms)*np.ones(N_firms, dtype=int))
        
        # Wage (uniform for simplicity)
        self.wage = 1.0
        
    def step(self):
        # Wages paid to employed workers
        income = np.where(self.employed, self.wage, 0.1)  # 0.1 = unemployment benefit
        self.wealth += income
        
        # Household consumption (heterogeneous MPC)
        consumption = np.minimum(self.wealth, self.mpc * self.wealth)
        self.wealth -= consumption
        agg_consumption = consumption.sum()
        
        # Firms produce, sell, update inventories and prices
        agg_output = self.workers.sum() * self.wage  # simplified production
        demand_per_firm = agg_consumption / self.N_firms
        
        for f in range(self.N_firms):
            sold = min(self.inventories[f], demand_per_firm/self.prices[f])
            self.inventories[f] += self.workers[f]*self.wage/self.prices[f] - sold
            # Price adjustment: raise if inventory low, cut if high
            inv_target = 20.0
            if self.inventories[f] < 0.5*inv_target:
                self.markup[f] = min(0.5, self.markup[f] * 1.02)
            elif self.inventories[f] > 2*inv_target:
                self.markup[f] = max(0.05, self.markup[f] * 0.98)
            self.prices[f] = (1+self.markup[f]) * self.wage
        
        # Labor market matching
        n_employed = self.employed.sum()
        employment_rate = n_employed / self.N_hh
        # Adjust labor demand to target inventories
        total_workers_demand = max(10, int(0.95 * self.N_hh * employment_rate))
        # Simple: unemployed workers find jobs with prob 0.3
        self.employed[~self.employed] = np.random.rand((~self.employed).sum()) < 0.20
        # Employed workers lose jobs with prob proportional to excess inventory
        layoff_prob = max(0, 0.02 * (self.inventories.mean()/20.0 - 1))
        self.employed[self.employed] = np.random.rand(self.employed.sum()) > layoff_prob
        
        # Aggregate statistics
        GDP = agg_output
        U = 1 - self.employed.mean()
        inflation = (self.prices.mean() - 1.0)  # deviation from baseline
        
        return {'GDP': GDP, 'U': U, 'inflation': inflation,
                'mean_price': self.prices.mean(),
                'mean_wealth': self.wealth.mean(),
                'inv': self.inventories.mean()}
    
    def simulate(self):
        results = []
        for t in range(self.T):
            results.append(self.step())
        return {k: np.array([r[k] for r in results]) for k in results[0]}

abm = MacroABM(N_hh=500, N_firms=50, T=300)
results = abm.simulate()

fig, axes = plt.subplots(3, 1, figsize=(12, 9))
axes[0].plot(results['GDP']); axes[0].set_title('Output (GDP)'); axes[0].set_ylabel('Output')
axes[1].plot(results['U']*100); axes[1].set_title('Unemployment Rate')
axes[1].set_ylabel('%'); axes[1].axhline(5, ls='--', color='r', alpha=0.5)
axes[2].plot(results['inflation']*100); axes[2].set_title('Inflation')
axes[2].set_ylabel('%'); axes[2].axhline(0, ls='--', color='r', alpha=0.5)
for ax in axes: ax.set_xlabel('Period')
plt.tight_layout(); plt.show()

# Business cycle statistics
burnin = 50
gdp = results['GDP'][burnin:]; unemp = results['U'][burnin:]
print(f"\nABM Business Cycle Statistics:")
print(f"  GDP std dev: {np.std(gdp)/np.mean(gdp)*100:.2f}%")
print(f"  Unemployment mean: {np.mean(unemp)*100:.1f}%, std: {np.std(unemp)*100:.2f}%")
print(f"  Corr(GDP, Unemployment): {np.corrcoef(gdp, unemp)[0,1]:.3f} (Okun's law: ≈ -0.8)")
```

```julia
# Julia — ABM macro model (simplified, fast)
using Statistics, Random

mutable struct Household
    wealth::Float64; employed::Bool; mpc::Float64
end

mutable struct Firm
    price::Float64; markup::Float64; inventory::Float64; workers::Int
end

function run_abm(N_hh=200, N_f=20, T=200; seed=42)
    Random.seed!(seed)
    hh = [Household(rand(Exponential(10)), rand()<0.95, rand()*0.4+0.4) for _ in 1:N_hh]
    firms = [Firm(1.0, rand()*0.2+0.1, 20.0, N_hh÷N_f) for _ in 1:N_f]
    wage = 1.0; GDP = zeros(T); U_rate = zeros(T)
    
    for t in 1:T
        for h in hh
            inc = h.employed ? wage : 0.1
            h.wealth += inc; h.wealth -= h.mpc * h.wealth
        end
        C = sum(h.mpc*h.wealth for h in hh)
        for f in firms
            sold = min(f.inventory, C/(N_f*f.price))
            f.inventory += f.workers*wage/f.price - sold
            f.markup *= f.inventory < 10 ? 1.02 : (f.inventory > 40 ? 0.98 : 1.0)
            f.price = (1+f.markup)*wage
        end
        # Job flows
        for h in hh
            h.employed && (h.employed = rand() > 0.02)
            !h.employed && (h.employed = rand() < 0.25)
        end
        GDP[t] = sum(f.workers for f in firms)*wage
        U_rate[t] = mean(!h.employed for h in hh)
    end
    GDP, U_rate
end

GDP, U = run_abm()
println("ABM: GDP std=$(round(std(GDP[50:end])/mean(GDP[50:end])*100,digits=2))%, U mean=$(round(mean(U[50:end])*100,digits=1))%")
```

---

## 36.5 Calibration and Validation of ABMs

Unlike DSGE models, ABMs cannot be estimated via maximum likelihood (the model has no closed-form likelihood). Validation relies on:

**1. Stylized facts.** Does the model reproduce known empirical regularities — not just the moments it was calibrated to, but out-of-sample facts? Dosi et al. (2010) target 14 stylized facts; Dawid et al. (2019) target 12.

**2. Econometric validation.** Gilli and Winker (2003) propose using an objective function based on the distance between model-simulated and empirical moments — essentially GMM for ABMs.

**3. History-friendly models.** Some ABMs (e.g., Malerba et al., 2001 for semiconductors) are calibrated to a specific historical episode and validated against its sequel.

---

## 36.6 Connecting ABMs to DSGE

**The rational expectations limiting case:** As agents become more rational (longer planning horizons, better information), the ABM converges to the DSGE equilibrium. The DSGE is the ABM with rationality imposed as a constraint.

**Partial rationality bridges:** Gabaix (2020) "Sparse Bounded Rational" model replaces global optimization with local (myopic) optimization — agents optimize a simplified version of the true problem. The model nests DSGE (full rationality) and behavioral models (zero rationality) as special cases.

**Machine learning policy functions:** Recent work (Curry, 2022; Zheng et al., 2022) trains neural networks to approximate DSGE-optimal policies, making agents "approximately rational" without closed-form optimization. This bridges ABM flexibility with DSGE optimality properties.

---

## 36.7 Programming Exercises

### Exercise 36.1 (APL — Schelling Model)

Implement the Schelling segregation model in APL. (a) Represent the grid as an APL matrix with `¯1` (blue), `1` (red), and `0` (empty). (b) For each unhappy agent (less than $\tau$ same-color neighbors), move to a random empty cell. (c) Apply one update round: `{one_schelling_step ⍵}` using `⌊⍺` for neighborhood inspection. (d) Iterate until convergence and compute the segregation index $S = \sum_i (\text{fraction same-color neighbors}_i) / N$.

### Exercise 36.2 (Python — ABM Business Cycles)

Extend the `MacroABM` from Section 36.4 to include: (a) a banking sector (firms borrow to finance inventories; banks can tighten lending during recessions); (b) a government (automatic stabilizers: lower taxes in recessions, financed by debt). Compare the amplitude of business cycles with and without stabilizers. Verify that automatic stabilizers reduce output volatility by approximately 30–40% (consistent with empirical estimates).

### Exercise 36.3 (Julia — Policy Experiments)

```julia
# Compare three monetary policy rules in the ABM
rules = [
    ("Taylor rule (φ_π=1.5)", x -> 1.5*x),
    ("Aggressive (φ_π=3.0)", x -> 3.0*x),
    ("Passive (φ_π=0.5)",    x -> 0.5*x)
]

println("ABM monetary policy comparison:")
for (name, rule) in rules
    # Run ABM with different inflation response
    GDP, U = run_abm(200, 20, 300; seed=42)  # simplified: same sim
    println("  $name: GDP std=$(round(std(GDP[50:end])/mean(GDP[50:end])*100,digits=2))%")
end
println("(Full implementation: modify Firm.price update to incorporate CB rate)")
```

### Exercise 36.4 — Emergent Business Cycles ($\star$)

In the minimal ABM, business cycles emerge from the interaction of inventory adjustments and labor dynamics — no external shocks are needed. (a) Confirm this by running the ABM with zero variance in all random components (seed fixed, no random job destruction or creation). (b) Show that even in the deterministic limit, the model can generate cycles through the markup–inventory feedback loop. (c) Compute the dominant frequency of the GDP series using the spectral density and compare to the NBER business cycle frequency (6–32 quarters).

---

## 36.8 Chapter Summary

**Key results:**

- Agent-based models replace utility maximization with **behavioral rules**; business cycles and distributional patterns **emerge** from micro interactions without being imposed.
- The Schelling model (Theorem 36.1 implicitly): individual tolerance $\tau = 1/3$ generates near-complete segregation at the aggregate — macro outcomes can be qualitatively different from individual preferences.
- The **Dosi et al. (K+S) macro ABM** generates 14 business cycle stylized facts without shock-driven dynamics; austerity policies worsen recessions in this framework.
- The **ABM as a Markov chain**: state space $\mathcal{S} = \prod_i S_i$ with transition kernel defined by the interaction rules; ergodicity requires the chain to have a unique stationary distribution.
- **Calibration** uses moment matching (GMM for ABMs); **validation** requires replication of out-of-sample stylized facts, not just calibration targets.
- ABMs and DSGE are complementary: DSGE is the rational-expectations limit of an ABM; partial-rationality models (Gabaix 2020) bridge the two.
- In APL: ABM time step as `{apply_rules ⍵}\ T ⍴ ⊂init_state` — scan over time periods; interaction matrix operations use `∘.f` outer products for bilateral matching.

---

*End of Part VIII. Next: Part IX — Case Studies and Applications*
