"""
Chapter 10: Agent-Based Models of Cooperation vs. Competition
=============================================================
Implements the BDI agent framework and the N=500 cooperation/competition
simulation described in the chapter.

Dependencies:
    pip install mesa numpy matplotlib pandas scipy

Mesa version: tested on mesa >= 2.0
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict

# ── try mesa import, give helpful message if missing ──────────────────────────
try:
    from mesa import Agent, Model
    from mesa.time import RandomActivation
    from mesa.space import NetworkGrid
    from mesa.datacollection import DataCollector
    import networkx as nx
except ImportError:
    raise ImportError(
        "Install mesa and networkx first:\n"
        "    pip install mesa networkx matplotlib numpy pandas scipy"
    )


# ═════════════════════════════════════════════════════════════════════════════
# 1.  BDI AGENT  (Definition 10.1)
# ═════════════════════════════════════════════════════════════════════════════

class BDIAgent(Agent):
    """
    Belief-Desire-Intention agent for the cooperation / competition simulation.

    Beliefs  : estimate of the fraction of cooperators in the neighbourhood
    Desires  : maximise own discounted payoff
    Intention: update strategy based on imitation of best-performing neighbour
               (bounded rationality — no global optimisation)

    Strategy space: {'cooperate', 'defect'}

    Payoff matrix (per interaction, before ecological penalty):
        CC : (R, R) = (3, 3)
        CD : (S, T) = (0, 5)
        DC : (T, S) = (5, 0)
        DD : (P, P) = (1, 1)

    Ecological penalty (Section 10.3):
        Each 'defect' action depletes the local natural capital stock by
        delta_NC per interaction; when NC < NC_crit the payoff from any
        action is multiplied by NC/NC_max (production function link).
    """

    # Payoff matrix parameters
    R, S, T, P = 3.0, 0.0, 5.0, 1.0

    def __init__(
        self,
        unique_id: int,
        model: "CooperationModel",
        initial_strategy: str = "cooperate",
        memory_length: int = 5,
        noise: float = 0.02,
    ):
        super().__init__(unique_id, model)
        self.strategy       = initial_strategy
        self.memory_length  = memory_length
        self.noise          = noise          # mutation / exploration rate

        # BDI components
        self.belief_coop_rate = 0.5          # belief: fraction cooperators nearby
        self.payoff           = 0.0          # current period payoff
        self.cumulative_payoff= 0.0
        self.natural_capital  = 1.0          # local NC stock (normalised)
        self.payoff_history   = []

    # ── BELIEFS: update from observed neighbourhood ───────────────────────────
    def update_beliefs(self):
        neighbours = self.model.grid.get_neighbors(self.pos, include_center=False)
        if not neighbours:
            return
        coop_count = sum(
            1 for nid in neighbours
            if self.model.schedule.agents[nid].strategy == "cooperate"
        )
        self.belief_coop_rate = coop_count / len(neighbours)

    # ── DESIRES + INTENTION: imitate best neighbour ───────────────────────────
    def update_strategy(self):
        neighbours = self.model.grid.get_neighbors(self.pos, include_center=False)
        if not neighbours:
            return

        # Find neighbour with highest recent payoff
        best_payoff   = self.payoff
        best_strategy = self.strategy

        for nid in neighbours:
            nb = self.model.schedule.agents[nid]
            if nb.payoff > best_payoff:
                best_payoff   = nb.payoff
                best_strategy = nb.strategy

        # Imitate with probability proportional to payoff difference
        if best_strategy != self.strategy:
            payoff_diff = best_payoff - self.payoff
            prob_switch = payoff_diff / (self.R - self.P + 1e-9)  # normalise
            prob_switch = min(max(prob_switch, 0), 0.9)
            if self.random.random() < prob_switch:
                self.strategy = best_strategy

        # Exploration / mutation (noise)
        if self.random.random() < self.noise:
            self.strategy = self.random.choice(["cooperate", "defect"])

    # ── PAYOFF: interact with each neighbour ──────────────────────────────────
    def compute_payoff(self):
        neighbours = self.model.grid.get_neighbors(self.pos, include_center=False)
        if not neighbours:
            self.payoff = 0.0
            return

        total = 0.0
        for nid in neighbours:
            nb = self.model.schedule.agents[nid]
            s1, s2 = self.strategy, nb.strategy

            if   s1 == "cooperate" and s2 == "cooperate": total += self.R
            elif s1 == "cooperate" and s2 == "defect":    total += self.S
            elif s1 == "defect"    and s2 == "cooperate": total += self.T
            else:                                          total += self.P

            # Ecological depletion: defection degrades local NC
            if s1 == "defect":
                self.natural_capital = max(
                    0.0,
                    self.natural_capital - self.model.nc_depletion_rate,
                )

        # NC penalty: production scales with NC stock
        nc_factor = self.natural_capital / 1.0  # NC_max = 1
        nc_factor = max(nc_factor, 0.1)         # floor to avoid zero payoff

        self.payoff = (total / len(neighbours)) * nc_factor
        self.cumulative_payoff += self.payoff

        self.payoff_history.append(self.payoff)
        if len(self.payoff_history) > self.memory_length:
            self.payoff_history.pop(0)

    # ── STEP ──────────────────────────────────────────────────────────────────
    def step(self):
        self.update_beliefs()
        self.compute_payoff()
        self.update_strategy()


# ═════════════════════════════════════════════════════════════════════════════
# 2.  MODEL
# ═════════════════════════════════════════════════════════════════════════════

class CooperationModel(Model):
    """
    N=500 agent cooperation-vs-competition simulation on a Watts-Strogatz
    small-world network.

    Parameters
    ----------
    n_agents         : number of agents (default 500 as in chapter)
    initial_coop_frac: fraction starting as cooperators (default 0.50)
    network_type     : 'watts_strogatz' | 'barabasi_albert' | 'random'
    k                : mean degree (WS) or edges per new node (BA)
    beta             : rewiring probability (WS only)
    nc_depletion_rate: natural capital lost per defection interaction
    noise            : mutation / exploration rate per agent per step
    seed             : random seed
    """

    def __init__(
        self,
        n_agents: int           = 500,
        initial_coop_frac: float= 0.50,
        network_type: str       = "watts_strogatz",
        k: int                  = 6,
        beta: float             = 0.10,
        nc_depletion_rate: float= 0.005,
        noise: float            = 0.02,
        seed: int               = 42,
    ):
        super().__init__()
        self.n_agents          = n_agents
        self.nc_depletion_rate = nc_depletion_rate
        self.schedule          = RandomActivation(self)
        self.random.seed(seed)

        # ── Build network ───────────────────────────────────────────────────
        if network_type == "watts_strogatz":
            G = nx.watts_strogatz_graph(n_agents, k, beta, seed=seed)
        elif network_type == "barabasi_albert":
            G = nx.barabasi_albert_graph(n_agents, k, seed=seed)
        else:
            G = nx.erdos_renyi_graph(n_agents, k / n_agents, seed=seed)

        self.G    = G
        self.grid = NetworkGrid(G)

        # ── Create and place agents ─────────────────────────────────────────
        n_coop = int(n_agents * initial_coop_frac)
        strategies = (["cooperate"] * n_coop +
                      ["defect"]    * (n_agents - n_coop))
        self.random.shuffle(strategies)

        for node_id, strategy in zip(G.nodes(), strategies):
            agent = BDIAgent(
                unique_id=node_id,
                model=self,
                initial_strategy=strategy,
                noise=noise,
            )
            self.schedule.add(agent)
            self.grid.place_agent(agent, node_id)

        # ── Data collection ─────────────────────────────────────────────────
        self.datacollector = DataCollector(
            model_reporters={
                "CoopFraction": lambda m: sum(
                    1 for a in m.schedule.agents
                    if a.strategy == "cooperate"
                ) / m.n_agents,
                "AvgPayoff": lambda m: np.mean(
                    [a.payoff for a in m.schedule.agents]
                ),
                "AvgNaturalCapital": lambda m: np.mean(
                    [a.natural_capital for a in m.schedule.agents]
                ),
                "GiniPayoff": lambda m: _gini(
                    np.array([a.cumulative_payoff for a in m.schedule.agents])
                ),
            },
            agent_reporters={
                "strategy":       "strategy",
                "payoff":         "payoff",
                "natural_capital":"natural_capital",
                "belief_coop":    "belief_coop_rate",
            },
        )

        self.datacollector.collect(self)

    def step(self):
        self.schedule.step()
        self.datacollector.collect(self)


# ═════════════════════════════════════════════════════════════════════════════
# 3.  UTILITIES
# ═════════════════════════════════════════════════════════════════════════════

def _gini(x: np.ndarray) -> float:
    """Gini coefficient — same formula as ch01."""
    x = np.sort(np.asarray(x, dtype=float))
    if x.sum() == 0:
        return 0.0
    n = len(x)
    i = np.arange(1, n + 1)
    return float((2 * (i * x).sum() - (n + 1) * x.sum()) / (n * x.sum()))


def run_simulation(
    steps: int = 200,
    network_type: str = "watts_strogatz",
    initial_coop_frac: float = 0.50,
    **model_kwargs,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run the CooperationModel and return (model_df, agent_df).
    """
    model = CooperationModel(
        network_type=network_type,
        initial_coop_frac=initial_coop_frac,
        **model_kwargs,
    )
    for _ in range(steps):
        model.step()

    model_df = model.datacollector.get_model_vars_dataframe()
    agent_df = model.datacollector.get_agent_vars_dataframe()
    return model_df, agent_df


# ═════════════════════════════════════════════════════════════════════════════
# 4.  MAIN EXPERIMENT — replicate Chapter 10 worked example
# ═════════════════════════════════════════════════════════════════════════════

def chapter10_experiment(steps: int = 200, seed: int = 42):
    """
    Compare cooperation dynamics on three network topologies:
      A. Watts-Strogatz small-world (cooperative baseline)
      B. Barabási-Albert scale-free (competitive topology)
      C. Erdős-Rényi random graph

    Reproduces the qualitative results of Table 10.1 and Figure 10.2.
    """
    configs = {
        "Watts-Strogatz\n(small-world, β=0.10)": dict(
            network_type="watts_strogatz", k=6, beta=0.10, seed=seed),
        "Barabási-Albert\n(scale-free)": dict(
            network_type="barabasi_albert", k=3, seed=seed),
        "Erdős-Rényi\n(random)": dict(
            network_type="random", k=6, seed=seed),
    }

    results = {}
    for label, kwargs in configs.items():
        print(f"  Running: {label.replace(chr(10), ' ')} …", end=" ", flush=True)
        mdf, _ = run_simulation(steps=steps, **kwargs)
        results[label] = mdf
        print(f"final coop fraction = {mdf['CoopFraction'].iloc[-1]:.3f}")

    # ── Plot ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes = axes.flatten()

    colors = ["steelblue", "firebrick", "seagreen"]
    metrics = ["CoopFraction", "AvgPayoff", "AvgNaturalCapital", "GiniPayoff"]
    ylabels = [
        "Cooperation fraction",
        "Average payoff per period",
        "Average natural capital",
        "Gini coefficient (cumulative payoff)",
    ]

    for ax, metric, ylabel in zip(axes, metrics, ylabels):
        for (label, mdf), color in zip(results.items(), colors):
            ax.plot(mdf.index, mdf[metric],
                    label=label.replace("\n", " "), color=color)
        ax.set_xlabel("Period")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)

    plt.suptitle(
        "Chapter 10: Cooperation vs. Competition ABM (N=500)\n"
        "Comparing Watts-Strogatz, Barabási-Albert, and Erdős-Rényi networks",
        fontsize=11,
    )
    plt.tight_layout()
    plt.savefig("ch10_cooperation_abm.png", dpi=150)
    plt.show()
    print("Saved: ch10_cooperation_abm.png")

    # ── Summary table ───────────────────────────────────────────────────────
    print("\n── Chapter 10, Table 10.1 (replicated) ──────────────────────────")
    print(f"{'Network':<30} {'Coop@t=200':>12} {'AvgPayoff@200':>15} "
          f"{'AvgNC@200':>12} {'Gini@200':>10}")
    print("-" * 80)
    for label, mdf in results.items():
        row = mdf.iloc[-1]
        print(f"{label.replace(chr(10),' '):<30} "
              f"{row['CoopFraction']:>12.3f} "
              f"{row['AvgPayoff']:>15.3f} "
              f"{row['AvgNaturalCapital']:>12.3f} "
              f"{row['GiniPayoff']:>10.3f}")
    print("-" * 80)


# ═════════════════════════════════════════════════════════════════════════════
# 5.  SOBOL SENSITIVITY ANALYSIS  (Chapter 10, Section 10.4)
# ═════════════════════════════════════════════════════════════════════════════

def sobol_sensitivity(n_samples: int = 64, steps: int = 100):
    """
    First-order Sobol sensitivity indices for three parameters:
        beta (WS rewiring), nc_depletion_rate, noise (mutation rate)
    Outcome: final cooperation fraction.

    Requires SALib:  pip install SALib
    """
    try:
        from SALib.sample import saltelli
        from SALib.analyze import sobol as sobol_analyze
    except ImportError:
        print("SALib not installed. Run:  pip install SALib")
        return

    problem = {
        "num_vars": 3,
        "names":    ["beta", "nc_depletion_rate", "noise"],
        "bounds":   [[0.01, 0.50], [0.001, 0.02], [0.005, 0.10]],
    }

    param_values = saltelli.sample(problem, n_samples, calc_second_order=False)

    Y = np.empty(len(param_values))
    for i, (beta, ncd, noise) in enumerate(param_values):
        mdf, _ = run_simulation(
            steps=steps,
            network_type="watts_strogatz",
            beta=float(beta),
            nc_depletion_rate=float(ncd),
            noise=float(noise),
            seed=0,
        )
        Y[i] = mdf["CoopFraction"].iloc[-1]
        if i % 20 == 0:
            print(f"  Sobol sample {i}/{len(param_values)}")

    Si = sobol_analyze.analyze(problem, Y, print_to_console=False)
    print("\n── Sobol First-Order Sensitivity Indices ─────────────")
    for name, s1, st in zip(problem["names"], Si["S1"], Si["ST"]):
        print(f"  {name:<22}: S1={s1:.3f}   ST={st:.3f}")


# ═════════════════════════════════════════════════════════════════════════════
# 6.  ENTRY POINT
# ═════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Chapter 10 — ABM: Cooperation vs. Competition")
    print("=" * 60)
    print("\nRunning three-topology comparison (N=500, 200 steps)…\n")
    chapter10_experiment(steps=200)

    print("\nRunning Sobol sensitivity (n_samples=64, 100 steps — may take ~2 min)…")
    sobol_sensitivity(n_samples=64, steps=100)
