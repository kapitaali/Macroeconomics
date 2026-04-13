"""
Chapter 36: Regenerative Agriculture and Landscape Restoration
==============================================================
Formal CBA for the UK Peak District conversion (10,000 ha)
and the soil capital dynamics model.

Dependencies: numpy, matplotlib, pandas
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ─────────────────────────────────────────────────────────────────────────────
# 1.  SOIL CAPITAL DYNAMICS MODEL  (Definition 36.1)
# ─────────────────────────────────────────────────────────────────────────────

def soil_capital_steady_state(
    alpha_M: float,   # SOC gain per unit organic input
    alpha_N: float,   # intrinsic regeneration rate
    K_s: float,       # carrying capacity (t C/ha)
    epsilon_s: float, # baseline respiration
    beta_T: float,    # tillage depletion rate
    T: float,         # tillage intensity
    M_s: float,       # organic input (t/ha/yr)
) -> float:
    """Steady-state SOC under a given management regime."""
    return (alpha_M * M_s / alpha_N - epsilon_s / alpha_N) / (
        1 + beta_T * T / alpha_N - 1 / K_s
    )


def soil_trajectory(
    N0: float, alpha_M, alpha_N, K_s, epsilon_s, beta_T,
    T_intensity, M_s, years=50
) -> np.ndarray:
    """Simulate SOC trajectory over time (discrete Euler)."""
    N = N0
    traj = [N]
    for _ in range(years):
        regen  = alpha_M * M_s + alpha_N * N * (1 - N / K_s)
        deplet = beta_T * T_intensity * N + epsilon_s
        N = max(0, N + regen - deplet)
        traj.append(N)
    return np.array(traj)


# ─────────────────────────────────────────────────────────────────────────────
# 2.  PEAK DISTRICT CBA  (Section 36.5)
# ─────────────────────────────────────────────────────────────────────────────

ECOSYSTEM_SERVICES = {
    "Food production (beef)":         {"baseline": 120,  "regenerative": 65},
    "Carbon sequestration":           {"baseline": -85,  "regenerative": 185},
    "Water regulation":               {"baseline": -40,  "regenerative": 95},
    "Biodiversity":                   {"baseline": 15,   "regenerative": 110},
    "Recreation/tourism":             {"baseline": 0,    "regenerative": 45},
}

TRANSITION_COSTS = {
    "Capital investment (fencing etc.)": 1_200,   # GBP/ha one-time
    "Foregone grouse shooting (5yr)":    120 * 5, # GBP/ha×yr × 5 years
}


def peak_district_cba(
    area_ha: float = 10_000,
    horizon: int   = 20,
    discount_rate: float = 0.035,
) -> dict:
    """
    Social CBA for converting 10,000 ha of Dark Peak grouse moor
    to regenerative upland farming.

    Returns
    -------
    dict with annual values, NPV, BCR
    """
    # Annual ecosystem service values (GBP/ha/year)
    baseline_total = sum(v["baseline"] for v in ECOSYSTEM_SERVICES.values())
    regen_total    = sum(v["regenerative"] for v in ECOSYSTEM_SERVICES.values())
    net_benefit_ha = regen_total - baseline_total

    # Transition cost (one-time, annualised over horizon)
    capital_cost = sum(TRANSITION_COSTS.values())  # GBP/ha

    # NPV calculation
    annual_benefits = net_benefit_ha * area_ha
    transition_cost = capital_cost  * area_ha

    discount_factors = np.array([(1 / (1 + discount_rate) ** t) for t in range(1, horizon + 1)])
    pv_benefits      = annual_benefits * discount_factors.sum()
    bcr              = pv_benefits / transition_cost

    return {
        "baseline_annual_ha":    baseline_total,
        "regen_annual_ha":       regen_total,
        "net_benefit_annual_ha": net_benefit_ha,
        "area_ha":               area_ha,
        "annual_benefit_total":  annual_benefits,
        "transition_cost_total": transition_cost,
        "pv_benefits":           pv_benefits,
        "npv":                   pv_benefits - transition_cost,
        "bcr":                   bcr,
        "horizon":               horizon,
        "discount_rate":         discount_rate,
    }


def plot_cba_results(result: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Waterfall chart of ecosystem services
    services = list(ECOSYSTEM_SERVICES.keys())
    baselines = [ECOSYSTEM_SERVICES[s]["baseline"] for s in services]
    regens    = [ECOSYSTEM_SERVICES[s]["regenerative"] for s in services]
    x = np.arange(len(services))
    width = 0.35
    axes[0].bar(x - width/2, baselines, width, label="Baseline (grouse moor)",
                color="firebrick", alpha=0.8)
    axes[0].bar(x + width/2, regens, width, label="Regenerative farming",
                color="steelblue", alpha=0.8)
    axes[0].axhline(0, color="black", linewidth=0.5)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([s[:20] for s in services], rotation=30, ha="right", fontsize=8)
    axes[0].set_ylabel("GBP / ha / year")
    axes[0].set_title("Ecosystem service values: baseline vs. regenerative")
    axes[0].legend()
    axes[0].grid(axis="y", alpha=0.3)

    # NPV accumulation over time
    discount_rate = result["discount_rate"]
    horizon       = result["horizon"]
    annual_benefit= result["annual_benefit_total"]
    trans_cost    = result["transition_cost_total"]

    cum_pv = np.array([
        annual_benefit * sum(1 / (1 + discount_rate) ** t for t in range(1, y + 1)) - trans_cost
        for y in range(1, horizon + 1)
    ])
    years = np.arange(1, horizon + 1)
    axes[1].plot(years, cum_pv / 1e6, color="steelblue")
    axes[1].axhline(0, color="grey", linestyle="--", linewidth=0.8)
    axes[1].fill_between(years, 0, cum_pv / 1e6,
                         where=(cum_pv >= 0), alpha=0.2, color="steelblue",
                         label="Positive NPV region")
    axes[1].set_xlabel("Year")
    axes[1].set_ylabel("Cumulative NPV (£M)")
    axes[1].set_title(f"NPV accumulation — BCR = {result['bcr']:.1f}:1")
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    plt.suptitle("Chapter 36: Peak District Regenerative Agriculture CBA\n"
                 f"(10,000 ha conversion, {horizon}-year horizon, {discount_rate*100:.1f}% discount)",
                 fontsize=10)
    plt.tight_layout()
    plt.savefig("ch36_peak_district_cba.png", dpi=150)
    plt.show()
    print("Saved: ch36_peak_district_cba.png")


if __name__ == "__main__":
    print("Chapter 36 — Regenerative Agriculture: CBA and Soil Dynamics")
    print("=" * 60)

    # Soil capital dynamics comparison
    print("\nSoil Carbon Steady States (UK upland agriculture):")
    params = dict(alpha_M=0.08, alpha_N=0.04, K_s=90, epsilon_s=0.8,
                  beta_T=0.02)
    N_conv = soil_capital_steady_state(**params, T=1.5, M_s=3)
    N_regen = soil_capital_steady_state(**params, T=0.2, M_s=8)
    print(f"  Conventional (T=1.5, M_s=3): N*_conv  = {N_conv:.1f} t C/ha")
    print(f"  Regenerative (T=0.2, M_s=8): N*_regen = {N_regen:.1f} t C/ha")

    # CBA
    print("\nPeak District CBA:")
    result = peak_district_cba()
    print(f"  Net annual benefit / ha  : GBP {result['net_benefit_annual_ha']:.0f}")
    print(f"  Total annual benefit     : GBP {result['annual_benefit_total']:,.0f}")
    print(f"  Transition cost          : GBP {result['transition_cost_total']:,.0f}")
    print(f"  PV of benefits (20yr)    : GBP {result['pv_benefits']:,.0f}")
    print(f"  NPV                      : GBP {result['npv']:,.0f}")
    print(f"  Benefit-Cost Ratio       : {result['bcr']:.1f}:1")
    plot_cba_results(result)
