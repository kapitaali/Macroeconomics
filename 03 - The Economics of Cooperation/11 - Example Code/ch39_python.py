"""
Chapter 39: Digital Sovereignty and Data Cooperatives
======================================================
Algorithm 39.1: Monte Carlo Shapley Value for Data Contributions

Dependencies: numpy, matplotlib, pandas
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from typing import Callable


def compute_data_shapley(
    contributors: list,
    value_function: Callable[[frozenset], float],
    n_samples: int = 2_000,
    seed: int = 42,
) -> dict:
    """
    Algorithm 39.1 — Monte Carlo Shapley approximation for data cooperatives.

    Parameters
    ----------
    contributors   : list of contributor IDs
    value_function : v(S: frozenset) -> float
                     Returns the data value for coalition S.
    n_samples      : number of random permutations to sample
    seed           : random seed

    Returns
    -------
    dict with keys:
        shapley_values : {contributor_id: shapley_value}
        std_errors     : {contributor_id: std_error}  (across permutations)
        total_value    : v(grand coalition)
        sigma_v        : cooperative surplus fraction
    """
    rng   = np.random.default_rng(seed)
    n     = len(contributors)
    phi   = {c: 0.0 for c in contributors}
    phi2  = {c: 0.0 for c in contributors}   # for variance estimation

    grand_coalition = frozenset(contributors)
    total_value     = value_function(grand_coalition)
    singleton_sum   = sum(value_function(frozenset([c])) for c in contributors)
    sigma_v         = (total_value - singleton_sum) / total_value if total_value > 0 else 0.0

    for _ in range(n_samples):
        perm      = rng.permutation(contributors)
        coalition = frozenset()
        v_prev    = 0.0

        for contributor in perm:
            new_coalition = coalition | {contributor}
            v_new         = value_function(new_coalition)
            marginal      = v_new - v_prev
            phi[contributor]  += marginal
            phi2[contributor] += marginal ** 2
            coalition = new_coalition
            v_prev    = v_new

    # Average
    for c in contributors:
        phi[c]  /= n_samples
        phi2[c] /= n_samples

    # Standard error estimates
    std_errors = {
        c: ((phi2[c] - phi[c]**2) / n_samples) ** 0.5
        for c in contributors
    }

    return {
        "shapley_values": phi,
        "std_errors":     std_errors,
        "total_value":    total_value,
        "sigma_v":        sigma_v,
    }


# ─────────────────────────────────────────────────────────────────────────────
# HEALTH DATA VALUE FUNCTION  (power-law + diversity bonus)
# ─────────────────────────────────────────────────────────────────────────────

def health_data_value(
    coalition: frozenset,
    contributor_types: dict,   # {id: 'standard' | 'genomic' | 'longitudinal'}
    alpha: float = 1.6,
    beta: float  = 100.0,
    diversity_bonus: float = 0.25,
) -> float:
    """
    v(S) = β × |S|^α × (1 + diversity_bonus × D(S))

    where D(S) ∈ [0,1] is the normalised type diversity of S.
    """
    if not coalition:
        return 0.0

    s = len(coalition)
    types_in_S = set(contributor_types[c] for c in coalition)
    n_types    = len(types_in_S)
    max_types  = len(set(contributor_types.values()))
    D = (n_types - 1) / max(max_types - 1, 1)

    return beta * (s ** alpha) * (1 + diversity_bonus * D)


# ─────────────────────────────────────────────────────────────────────────────
# TAXI DRIVER DATA COOPERATIVE  (Section 39.6 worked example)
# ─────────────────────────────────────────────────────────────────────────────

def taxi_cooperative_shapley(
    n_drivers: int   = 1_000,
    n_samples: int   = 500,
    seed: int        = 42,
    annual_revenue:  float = 697_700.0,   # EUR
    driver_share:    float = 0.70,
) -> pd.DataFrame:
    """
    Compute Shapley dividends for a taxi driver data cooperative.

    Contribution score dimensions (weights from Section 39.6.2):
      - Trip volume   (weight 0.50)
      - Data quality  (weight 0.25)
      - Diversity     (weight 0.15)
      - Governance    (weight 0.10)
    """
    rng = np.random.default_rng(seed)

    # Generate synthetic driver profiles
    avg_trips = 8_000
    trips   = rng.lognormal(np.log(avg_trips), 0.4, n_drivers).astype(int)
    quality = rng.beta(8, 2, n_drivers)           # mostly high quality
    diversity = rng.beta(2, 3, n_drivers)         # skewed low
    governance= rng.beta(3, 5, n_drivers)

    # Contribution scores (OVA)
    scores = (0.50 * trips / avg_trips +
              0.25 * quality +
              0.15 * diversity +
              0.10 * governance)

    # Shapley dividend proportional to score
    distributable = annual_revenue * driver_share
    total_score   = scores.sum()
    dividends     = (scores / total_score) * distributable

    df = pd.DataFrame({
        "driver_id":   range(n_drivers),
        "trips":       trips,
        "quality":     quality.round(3),
        "diversity":   diversity.round(3),
        "governance":  governance.round(3),
        "ova_score":   scores.round(4),
        "dividend_EUR": dividends.round(2),
    })

    return df


def run_taxi_demo() -> None:
    df = taxi_cooperative_shapley()
    print("\n── Taxi Data Cooperative (Section 39.6) ───────────────────")
    print(f"  Drivers                   : {len(df):,}")
    print(f"  Annual revenue (EUR)      : 697,700")
    print(f"  Driver share (70%)        : 488,390")
    print(f"\n  Dividend distribution (EUR/year):")
    print(f"  Mean                      : {df.dividend_EUR.mean():.0f}")
    print(f"  Median                    : {df.dividend_EUR.median():.0f}")
    print(f"  Min (low contributor)     : {df.dividend_EUR.min():.0f}")
    print(f"  Max (high contributor)    : {df.dividend_EUR.max():.0f}")
    print(f"  Std deviation             : {df.dividend_EUR.std():.0f}")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].hist(df.dividend_EUR, bins=40, color="steelblue", edgecolor="white")
    axes[0].set_xlabel("Annual dividend (EUR)")
    axes[0].set_ylabel("Number of drivers")
    axes[0].set_title("Dividend distribution (Shapley-OVA allocation)")
    axes[0].axvline(df.dividend_EUR.mean(), color="red", linestyle="--",
                    label=f"Mean = €{df.dividend_EUR.mean():.0f}")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    axes[1].scatter(df.trips, df.dividend_EUR, alpha=0.2, s=4, color="steelblue")
    axes[1].set_xlabel("Annual trips")
    axes[1].set_ylabel("Annual dividend (EUR)")
    axes[1].set_title("Dividend vs. trip volume")
    axes[1].grid(alpha=0.3)

    plt.suptitle("Chapter 39: Taxi Driver Data Cooperative\n"
                 "Shapley-OVA allocation of EUR 488,390 among 1,000 drivers",
                 fontsize=10)
    plt.tight_layout()
    plt.savefig("ch39_taxi_cooperative.png", dpi=150)
    plt.show()
    print("Saved: ch39_taxi_cooperative.png")


def demo_health_data_shapley(n_contributors: int = 20, n_samples: int = 1_000):
    """Small health-data cooperative Shapley demonstration."""
    rng = np.random.default_rng(0)
    types_list = ["standard", "genomic", "longitudinal"]
    contributors = list(range(n_contributors))
    contributor_types = {
        c: types_list[rng.integers(0, 3)] for c in contributors
    }

    def vf(S):
        return health_data_value(S, contributor_types)

    result = compute_data_shapley(contributors, vf, n_samples=n_samples)
    phi    = result["shapley_values"]

    print("\n── Health Data Cooperative Shapley Values (n=20 contributors) ──")
    print(f"  Grand coalition value  : {result['total_value']:.1f}")
    print(f"  Cooperative surplus σ_v: {result['sigma_v']:.3f}")
    print(f"  Mean Shapley value     : {np.mean(list(phi.values())):.1f}")
    print(f"  Total check (= v(N)?)  : {sum(phi.values()):.1f}")


if __name__ == "__main__":
    print("Chapter 39 — Data Cooperatives: Shapley Value Allocation")
    print("=" * 60)
    demo_health_data_shapley()
    run_taxi_demo()
