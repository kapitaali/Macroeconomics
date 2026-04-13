"""
Chapter 28: SFC Comparison of Four Monetary Regimes
====================================================
Simulates debt-based, sovereign, mutual credit, and demurrage economies
over 50 years and compares stability, growth, inequality, and ecology.

Dependencies: numpy, matplotlib, pandas
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dataclasses import dataclass


@dataclass
class MonetaryRegimeParams:
    name:               str
    interest_rate:      float = 0.04    # i^L
    deposit_rate:       float = 0.01    # i^D
    money_growth_rule:  float = 0.03    # mu-bar (for sovereign/demurrage)
    demurrage_rate:     float = 0.00    # delta
    minsky_channel:     bool  = False   # whether Minsky amplification is active
    velocity_base:      float = 1.35    # V_0
    investment_rate:    float = 0.18    # s_K
    gini_shapley:       float = 0.000   # psi_Shapley (equalising force)
    seigniorage_div:    float = 0.000   # seigniorage dividend to households
    nc_growth_penalty:  float = 0.030   # D(Y) = coeff × Y (ecological depletion)


REGIMES = {
    "Debt-based": MonetaryRegimeParams(
        name="Debt-based",
        interest_rate=0.050, deposit_rate=0.010,
        minsky_channel=True, velocity_base=1.35, investment_rate=0.182,
        gini_shapley=0.000, seigniorage_div=0.000,
    ),
    "Sovereign money": MonetaryRegimeParams(
        name="Sovereign money",
        interest_rate=0.040, deposit_rate=0.000,
        minsky_channel=False, velocity_base=1.50, investment_rate=0.204,
        gini_shapley=0.000, seigniorage_div=0.002,
    ),
    "Mutual credit": MonetaryRegimeParams(
        name="Mutual credit",
        interest_rate=0.000, deposit_rate=0.000,
        minsky_channel=False, velocity_base=1.65, investment_rate=0.221,
        gini_shapley=0.000, seigniorage_div=0.000,
    ),
    "Demurrage": MonetaryRegimeParams(
        name="Demurrage",
        interest_rate=0.040, deposit_rate=0.000, demurrage_rate=0.025,
        minsky_channel=False, velocity_base=1.76, investment_rate=0.233,
        gini_shapley=0.000, seigniorage_div=0.002,
    ),
}


def simulate_regime(
    params: MonetaryRegimeParams,
    T: int              = 50,
    Y0: float           = 100.0,
    K0: float           = 320.0,
    N0: float           = 1.0,
    G0: float           = 0.68,
    alpha: float        = 0.35,
    gamma: float        = 0.07,
    delta_K: float      = 0.053,
    r_N: float          = 0.04,
    K_N: float          = 1.2,
    shock_period: int   = 25,
    shock_size: float   = 0.0,    # 0 = no shock; 0.25 = 25% productivity shock
    shock_duration: int = 5,
    seed: int           = 42,
) -> pd.DataFrame:
    """
    Simulate one monetary regime for T periods.
    Uses discrete-time approximation of the ODE system (Eqs 29.12–29.20).
    """
    rng = np.random.default_rng(seed)

    Y, K, N, G = Y0, K0, N0, G0
    M = Y * 0.85   # initial money stock ≈ 85% of GDP

    records = []

    for t in range(T + 1):
        # Productivity shock
        tfp = 1.0
        if shock_size > 0 and shock_period <= t < shock_period + shock_duration:
            tfp = 1.0 - shock_size

        # Production
        Y_prod = tfp * (K ** alpha) * (N ** gamma) * ((1 - alpha - gamma) ** 0)
        # Scale so Y0 = 100 at t=0
        if t == 0:
            scale = Y0 / Y_prod
        Y = Y_prod * scale

        # Investment
        I = params.investment_rate * Y

        # Capital accumulation
        K = K + I - delta_K * K
        K = max(K, 1.0)

        # Natural capital dynamics (logistic + depletion)
        regen  = r_N * N * (1 - N / K_N)
        deplet = params.nc_growth_penalty * (Y / Y0)
        N = max(0.01, N + regen - deplet)

        # Gini dynamics  (Theorem 32.1)
        r_eff = params.interest_rate - params.demurrage_rate
        g_eff = (K - K0) / (K0 * max(t, 1))    # approximate growth rate
        psi   = params.gini_shapley + params.seigniorage_div
        G = G + (r_eff - 0.02 - params.demurrage_rate) * G - psi
        G = min(max(G, 0.01), 0.99)

        # Minsky shock amplification (debt-based only)
        crisis_loss = 0.0
        if params.minsky_channel and shock_size > 0:
            if shock_period <= t < shock_period + shock_duration:
                crisis_loss = 0.13 * Y  # Minsky amplification ~13% additional GDP loss

        records.append({
            "period":    t,
            "Y":         Y - crisis_loss,
            "K":         K,
            "N":         N,
            "G":         G,
            "growth":    0.0 if t == 0 else (Y - records[-1]["Y"]) / records[-1]["Y"],
        })

    return pd.DataFrame(records)


def run_comparison(shock: bool = True, T: int = 50) -> dict[str, pd.DataFrame]:
    results = {}
    for name, params in REGIMES.items():
        df = simulate_regime(
            params, T=T,
            shock_size=0.25 if shock else 0.0,
        )
        results[name] = df
    return results


def plot_comparison(results: dict, shock: bool = True) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    colors = {"Debt-based": "firebrick", "Sovereign money": "steelblue",
              "Mutual credit": "seagreen", "Demurrage": "darkorange"}
    metrics = [
        ("Y",  "GDP (index, t=0 = 100)"),
        ("G",  "Wealth Gini coefficient"),
        ("N",  "Natural capital stock"),
        ("growth", "Annual GDP growth rate"),
    ]

    for ax, (col, ylabel) in zip(axes.flatten(), metrics):
        for name, df in results.items():
            ax.plot(df.period, df[col], label=name, color=colors[name])
        if shock:
            ax.axvline(25, color="grey", linestyle="--", linewidth=0.8,
                       label="Shock at t=25")
        ax.set_xlabel("Period (years)")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)

    title = ("Chapter 28: Four-Regime SFC Comparison (50 years)\n" +
             ("With 25% productivity shock at t=25" if shock else "Baseline (no shock)"))
    plt.suptitle(title, fontsize=10)
    plt.tight_layout()
    fname = "ch28_sfc_comparison_shock.png" if shock else "ch28_sfc_comparison_baseline.png"
    plt.savefig(fname, dpi=150)
    plt.show()
    print(f"Saved: {fname}")


def print_summary_table(results: dict) -> None:
    print("\n── Chapter 28, Table 28.4 (replicated) ─────────────────────────────")
    print(f"{'Regime':<18} {'GDP@50':>8} {'Gini@50':>9} {'NC@50':>8} {'CrisisLoss':>12}")
    print("-" * 62)
    baseline = {k: simulate_regime(v, shock_size=0.0) for k, v in REGIMES.items()}
    for name, df in results.items():
        base_y50 = baseline[name].iloc[-1]["Y"]
        y50  = df.iloc[-1]["Y"]
        perm_loss = (base_y50 - y50) / base_y50
        g50  = df.iloc[-1]["G"]
        n50  = df.iloc[-1]["N"]
        print(f"{name:<18} {y50:>8.1f} {g50:>9.3f} {n50:>8.3f} {perm_loss*100:>10.1f}%")


if __name__ == "__main__":
    print("Chapter 28 — SFC Comparison: Four Monetary Regimes")
    print("=" * 60)
    print("\nBaseline (no shock):")
    results_base = run_comparison(shock=False)
    plot_comparison(results_base, shock=False)

    print("\nWith 25% productivity shock at t=25:")
    results_shock = run_comparison(shock=True)
    plot_comparison(results_shock, shock=True)
    print_summary_table(results_shock)
