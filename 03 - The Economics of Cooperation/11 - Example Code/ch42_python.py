"""
Chapter 42: Experimentation — Thompson Sampling for Policy Learning
===================================================================
Algorithm 42.1: Thompson Sampling for cooperative-regenerative policy arms.

Dependencies: numpy, matplotlib, pandas
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from dataclasses import dataclass, field


@dataclass
class PolicyArm:
    """A candidate cooperative-regenerative institutional design to be tested."""
    name: str
    true_welfare_effect: float       # unknown to planner (ground truth)
    prior_mean: float  = 0.0
    prior_var:  float  = 1.0

    # Running posterior (updated via Bayesian conjugate normal-normal)
    post_mean:  float  = field(init=False)
    post_var:   float  = field(init=False)
    n_trials:   int    = field(init=False, default=0)
    observations: list = field(init=False, default_factory=list)

    def __post_init__(self):
        self.post_mean = self.prior_mean
        self.post_var  = self.prior_var

    def observe(self, noise_var: float, rng: np.random.Generator) -> float:
        """Sample a welfare signal from this arm (noisy observation)."""
        obs = self.true_welfare_effect + rng.normal(0, noise_var ** 0.5)
        self.observations.append(obs)
        return obs

    def update(self, observation: float, noise_var: float) -> None:
        """Bayesian normal-normal conjugate update."""
        prior_prec  = 1.0 / self.post_var
        like_prec   = 1.0 / noise_var
        post_prec   = prior_prec + like_prec
        self.post_var  = 1.0 / post_prec
        self.post_mean = (prior_prec * self.post_mean + like_prec * observation) / post_prec
        self.n_trials += 1

    def sample_posterior(self, rng: np.random.Generator) -> float:
        return rng.normal(self.post_mean, self.post_var ** 0.5)


def thompson_sampling(
    arms: list[PolicyArm],
    n_periods: int,
    noise_var: float = 0.25,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Algorithm 42.1 — Thompson Sampling for Cooperative Policy Learning.

    Parameters
    ----------
    arms      : list of PolicyArm instances with unknown true effects
    n_periods : number of experimentation periods (e.g., electoral terms)
    noise_var : observation noise variance
    seed      : random seed

    Returns
    -------
    DataFrame with columns: period, chosen_arm, observed_welfare,
                             cumulative_welfare, [arm_post_mean for each arm]
    """
    rng     = np.random.default_rng(seed)
    records = []

    cumulative_welfare = 0.0
    best_true   = max(a.true_welfare_effect for a in arms)

    for t in range(1, n_periods + 1):
        # Thompson Sampling: draw from each arm's posterior
        samples    = [a.sample_posterior(rng) for a in arms]
        chosen_idx = int(np.argmax(samples))
        chosen_arm = arms[chosen_idx]

        # Observe welfare signal
        obs = chosen_arm.observe(noise_var, rng)
        chosen_arm.update(obs, noise_var)
        cumulative_welfare += obs

        row = {
            "period":             t,
            "chosen_arm":         chosen_arm.name,
            "observed_welfare":   obs,
            "cumulative_welfare": cumulative_welfare,
            "regret_cumulative":  t * best_true - cumulative_welfare,
        }
        for a in arms:
            row[f"post_mean_{a.name}"] = a.post_mean
        records.append(row)

    return pd.DataFrame(records)


def chapter42_demo():
    """
    Reproduce the worked example from Section 42.2.1:
    Four cooperative enterprise zone design variants as policy arms.

    Arm 4 is genuinely best (true effect = 6.0) but unknown ex ante.
    The planner starts with identical priors.
    """
    arms = [
        PolicyArm("CEZ-Design-A",   true_welfare_effect=3.0),
        PolicyArm("CEZ-Design-B",   true_welfare_effect=2.5),
        PolicyArm("CEZ-Design-C",   true_welfare_effect=3.8),
        PolicyArm("CEZ-Design-D",   true_welfare_effect=6.0),   # best arm
    ]

    print("\n── Thompson Sampling Policy Learning (Section 42.2) ────────")
    print(f"  Arms      : {[a.name for a in arms]}")
    print(f"  True effects: {[a.true_welfare_effect for a in arms]}")
    print(f"  Best arm  : {max(arms, key=lambda a: a.true_welfare_effect).name} (welfare = 6.0)")
    print(f"  Running {20} periods…\n")

    df = thompson_sampling(arms, n_periods=20, noise_var=0.5)

    # How quickly does sampling concentrate on arm D?
    best_arm_trials = (df.chosen_arm == "CEZ-Design-D").cumsum()
    best_arm_frac   = best_arm_trials / (df.period)

    print("  Period | Chosen arm   | Post mean D | Best arm %")
    print("  " + "-" * 52)
    for _, row in df.iterrows():
        p   = int(row.period)
        ca  = row.chosen_arm
        pmd = row[f"post_mean_CEZ-Design-D"]
        bf  = best_arm_frac.iloc[p-1]
        if p in [1, 3, 5, 8, 12, 16, 20]:
            print(f"  {p:>6} | {ca:<14}| {pmd:>11.2f} | {bf*100:>8.1f}%")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # Arm selection frequency over time
    for a in arms:
        mask = df.chosen_arm == a.name
        rolling = mask.rolling(5).mean()
        axes[0].plot(df.period, rolling, label=a.name)
    axes[0].set_title("Rolling arm selection frequency (5-period window)")
    axes[0].set_xlabel("Period")
    axes[0].set_ylabel("Selection frequency")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)

    # Posterior means over time
    for a in arms:
        col = f"post_mean_{a.name}"
        axes[1].plot(df.period, df[col], label=a.name)
        axes[1].axhline(a.true_welfare_effect, linestyle=":", color="grey", linewidth=0.5)
    axes[1].set_title("Posterior mean estimates (converging to truth)")
    axes[1].set_xlabel("Period")
    axes[1].set_ylabel("Posterior mean")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)

    # Cumulative regret
    axes[2].plot(df.period, df.regret_cumulative, color="firebrick")
    axes[2].set_title("Cumulative regret\n(welfare gap from optimal)")
    axes[2].set_xlabel("Period")
    axes[2].set_ylabel("Cumulative regret")
    axes[2].grid(alpha=0.3)

    plt.suptitle("Chapter 42: Thompson Sampling for Cooperative Policy Learning\n"
                 "(4 cooperative enterprise zone design variants)", fontsize=10)
    plt.tight_layout()
    plt.savefig("ch42_thompson_sampling.png", dpi=150)
    plt.show()
    print("\nSaved: ch42_thompson_sampling.png")


if __name__ == "__main__":
    print("Chapter 42 — Experimentation: Thompson Sampling")
    print("=" * 60)
    chapter42_demo()
