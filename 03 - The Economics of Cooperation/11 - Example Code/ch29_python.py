"""
Chapter 29: The Unified Model — Danish CRE Calibration
=======================================================
Simulates the unified cooperative-regenerative economy model
and computes the IPI welfare comparison (CE vs. CRE).

Dependencies: numpy, matplotlib, pandas
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def simulate_unified(
    cooperative: bool = True,
    T: int            = 50,
    Y0: float         = 100.0,
    K0: float         = 320.0,
    N0: float         = 0.71,
    G0: float         = 0.64,
    alpha: float      = 0.35,
    gamma: float      = 0.07,
    r_N: float        = 0.04,
    K_N: float        = 1.2,
    delta_K: float    = 0.053,
    sigma_v: float    = 0.18,      # cooperative surplus fraction (CRE)
    psi_shapley: float= 0.0028,    # Shapley equalising force
    lambda2_coop: float = 1.95,    # Fiedler value (CRE)
    lambda2_ce: float   = 1.28,    # Fiedler value (CE)
    beta: float       = 0.97,      # discount factor for IPI
    seed: int         = 42,
) -> pd.DataFrame:
    """
    Simulate the unified model (Eqs 29.12-29.20) for T periods.

    cooperative=True  →  CRE regime  (Shapley allocation, high Fiedler)
    cooperative=False →  CE regime   (competitive, no cooperative surplus)
    """
    rng = np.random.default_rng(seed)

    Y, K, N, G = Y0, K0, N0, G0
    records = []

    # Scale factor for production function
    Y_raw = (K ** alpha) * (N ** gamma)
    scale = Y0 / Y_raw

    for t in range(T + 1):
        Y = scale * (K ** alpha) * (N ** gamma)

        if cooperative:
            # Cooperative surplus boost to output
            Y = Y * (1 + sigma_v * 0.3)  # 30% of surplus passes through to Y
            s_K = 0.20
            psi = psi_shapley
            lambda2 = lambda2_coop
        else:
            s_K = 0.18
            psi = 0.0
            lambda2 = lambda2_ce

        # Investment
        I = s_K * Y

        # Capital
        K = max(K + I - delta_K * K, 1.0)

        # Natural capital — cooperative governance improves stewardship
        nc_bonus = 0.005 if cooperative else 0.0
        regen  = r_N * N * (1 - N / K_N) + nc_bonus
        deplet = 0.03 * (Y / Y0)
        N = max(0.01, N + regen - deplet)

        # Gini dynamics
        G = G + (0.04 - 0.02 - (0.025 if cooperative else 0)) * G - psi
        G = min(max(G, 0.05), 0.99)

        # Welfare: U(C) = log(C) approximately; C ≈ (1 - s_K) * Y
        C     = (1 - s_K) * Y
        U_per_capita = np.log(max(C, 0.01)) * (1 - G)   # inequality-adjusted

        records.append({
            "period":  t,
            "Y":       Y,
            "K":       K,
            "N":       N,
            "G":       G,
            "C":       C,
            "U":       U_per_capita,
            "IPI_disc": beta ** t * U_per_capita,
        })

    return pd.DataFrame(records)


def compute_ipi_gap(T: int = 50, beta: float = 0.97) -> dict:
    """Compute IPI(CRE) - IPI(CE) and welfare decomposition."""
    df_cre = simulate_unified(cooperative=True,  T=T, beta=beta)
    df_ce  = simulate_unified(cooperative=False, T=T, beta=beta)

    ipi_cre = df_cre["IPI_disc"].sum()
    ipi_ce  = df_ce["IPI_disc"].sum()

    return {
        "IPI_CRE":   ipi_cre,
        "IPI_CE":    ipi_ce,
        "IPI_gap":   ipi_cre - ipi_ce,
        "IPI_gap_pct": (ipi_cre - ipi_ce) / ipi_ce * 100,
        "Y_CRE_50":  df_cre.iloc[-1]["Y"],
        "Y_CE_50":   df_ce.iloc[-1]["Y"],
        "N_CRE_50":  df_cre.iloc[-1]["N"],
        "N_CE_50":   df_ce.iloc[-1]["N"],
        "G_CRE_50":  df_cre.iloc[-1]["G"],
        "G_CE_50":   df_ce.iloc[-1]["G"],
        "df_cre":    df_cre,
        "df_ce":     df_ce,
    }


def plot_cre_vs_ce(results: dict) -> None:
    df_cre = results["df_cre"]
    df_ce  = results["df_ce"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    for ax, col, ylabel in zip(
        axes.flatten(),
        ["Y", "G", "N", "IPI_disc"],
        ["GDP (index)", "Wealth Gini", "Natural capital", "Discounted welfare U"],
    ):
        ax.plot(df_cre.period, df_cre[col], label="CRE", color="steelblue")
        ax.plot(df_ce.period,  df_ce[col],  label="CE",  color="firebrick", linestyle="--")
        ax.set_xlabel("Period (years)")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)
        ax.legend()
        ax.grid(alpha=0.3)

    pct = results["IPI_gap_pct"]
    plt.suptitle(
        f"Chapter 29: CRE vs. CE Welfare Comparison (Danish calibration)\n"
        f"IPI(CRE) − IPI(CE) = +{pct:.1f}% over 50 years",
        fontsize=10,
    )
    plt.tight_layout()
    plt.savefig("ch29_cre_vs_ce.png", dpi=150)
    plt.show()
    print("Saved: ch29_cre_vs_ce.png")


if __name__ == "__main__":
    print("Chapter 29 — Unified Model: CRE vs. CE Comparison")
    print("=" * 60)
    results = compute_ipi_gap(T=50, beta=0.97)
    print(f"\nDanish calibration results:")
    print(f"  IPI(CRE)  = {results['IPI_CRE']:.2f}")
    print(f"  IPI(CE)   = {results['IPI_CE']:.2f}")
    print(f"  IPI gap   = +{results['IPI_gap_pct']:.1f}%  (CRE welfare advantage)")
    print(f"\n  GDP at year 50:    CRE={results['Y_CRE_50']:.1f}  CE={results['Y_CE_50']:.1f}")
    print(f"  Natural cap at 50: CRE={results['N_CRE_50']:.3f}  CE={results['N_CE_50']:.3f}")
    print(f"  Gini at year 50:   CRE={results['G_CRE_50']:.3f}  CE={results['G_CE_50']:.3f}")
    plot_cre_vs_ce(results)
