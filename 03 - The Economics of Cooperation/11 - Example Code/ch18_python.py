"""
Chapter 18: SFC-N — Stock-Flow Consistent Modeling with Natural Capital
========================================================================
Implements the Provisioning Balance Sheet (PBS) and SFC-N Transaction Flow
Matrix for the Darlington worked example from Chapter 44.

Dependencies: numpy, pandas, matplotlib
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dataclasses import dataclass, field


# ─────────────────────────────────────────────────────────────────────────────
# 1.  NATURAL CAPITAL STOCK
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class NaturalCapitalStock:
    """One natural capital stock in the SFC-N framework."""
    name:            str
    current_level:   float          # N_j (normalised, 1.0 = pre-industrial)
    carrying_capacity: float        # K_j
    regeneration_rate: float        # r_j
    critical_threshold: float       # N_j^critical
    depletion_coeff:  float         # d_j: depletion = d_j * Y
    shadow_price:     float         # p^N_j (£ or € per unit)

    @property
    def annual_regeneration(self) -> float:
        """Logistic regeneration: R(N) = r * N * (1 - N/K)"""
        return (self.regeneration_rate * self.current_level *
                (1 - self.current_level / self.carrying_capacity))

    def annual_depletion(self, gdp_index: float) -> float:
        return self.depletion_coeff * gdp_index

    def net_change(self, gdp_index: float) -> float:
        return self.annual_regeneration - self.annual_depletion(gdp_index)

    @property
    def value(self) -> float:
        """Market value: N_j × p^N_j"""
        return self.current_level * self.shadow_price

    @property
    def stewardship_met(self) -> bool:
        return self.net_change(1.0) >= 0   # at unit GDP index

    def step(self, gdp_index: float, investment: float = 0.0) -> None:
        """Advance one period: update N_j."""
        delta = self.net_change(gdp_index) + investment
        self.current_level = max(0.0, self.current_level + delta)


# ─────────────────────────────────────────────────────────────────────────────
# 2.  PROVISIONING BALANCE SHEET  (Definition 18.4 / Ch.44)
# ─────────────────────────────────────────────────────────────────────────────

class ProvisioningBalanceSheet:
    """
    PBS = Produced capital + Financial assets + Σ p^N_j N_j - Liabilities

    Stewardship Condition: PBS must be non-declining (dPBS/dt >= 0).
    """

    def __init__(
        self,
        produced_capital: float,
        financial_assets: float,
        liabilities: float,
        nc_stocks: list[NaturalCapitalStock],
    ):
        self.produced_capital = produced_capital
        self.financial_assets = financial_assets
        self.liabilities      = liabilities
        self.nc_stocks        = nc_stocks
        self.history: list[float] = [self.total]

    @property
    def nc_value(self) -> float:
        return sum(s.value for s in self.nc_stocks)

    @property
    def total(self) -> float:
        return self.produced_capital + self.financial_assets + self.nc_value - self.liabilities

    @property
    def stewardship_met(self) -> bool:
        if len(self.history) < 2:
            return True
        return self.history[-1] >= self.history[-2]

    def record(self) -> None:
        self.history.append(self.total)

    def report(self) -> pd.DataFrame:
        rows = []
        rows.append({"Component": "Produced capital", "Value (£M)": self.produced_capital})
        rows.append({"Component": "Financial assets",  "Value (£M)": self.financial_assets})
        for s in self.nc_stocks:
            rows.append({"Component": f"NC: {s.name}", "Value (£M)": s.value})
        rows.append({"Component": "Liabilities (−)",   "Value (£M)": -self.liabilities})
        rows.append({"Component": "TOTAL PBS",          "Value (£M)": self.total})
        df = pd.DataFrame(rows)
        return df


# ─────────────────────────────────────────────────────────────────────────────
# 3.  DARLINGTON CALIBRATION  (Chapter 44 worked example)
# ─────────────────────────────────────────────────────────────────────────────

def darlington_nc_stocks() -> list[NaturalCapitalStock]:
    """
    Five natural capital stocks for Darlington / Tees Valley,
    calibrated from Chapter 44 Section 2b.

    Shadow prices from the SFC-N framework (Chapter 18).
    Depletion coefficients calibrated so that at GDP index = 1.0
    the baseline depletion matches empirical estimates.
    """
    return [
        NaturalCapitalStock(
            name="Agricultural soil (Tees Valley)",
            current_level=0.78,           # 78% of pre-industrial SOC
            carrying_capacity=1.0,
            regeneration_rate=0.04,
            critical_threshold=0.50,
            depletion_coeff=0.036,        # net loss ~3.8M£/yr at Y=1
            shadow_price=285.0,           # £M total stock value
        ),
        NaturalCapitalStock(
            name="Urban green space",
            current_level=0.85,
            carrying_capacity=1.0,
            regeneration_rate=0.10,
            critical_threshold=0.40,
            depletion_coeff=0.012,
            shadow_price=45.0,
        ),
        NaturalCapitalStock(
            name="River Tees water quality",
            current_level=0.72,
            carrying_capacity=1.0,
            regeneration_rate=0.08,
            critical_threshold=0.50,
            depletion_coeff=0.028,
            shadow_price=62.0,
        ),
        NaturalCapitalStock(
            name="Local biodiversity",
            current_level=0.68,
            carrying_capacity=1.0,
            regeneration_rate=0.05,
            critical_threshold=0.30,
            depletion_coeff=0.022,
            shadow_price=38.0,
        ),
        NaturalCapitalStock(
            name="Carbon sequestration capacity",
            current_level=0.74,
            carrying_capacity=1.0,
            regeneration_rate=0.06,
            critical_threshold=0.40,
            depletion_coeff=0.030,
            shadow_price=118.0,
        ),
    ]


def run_pbs_simulation(
    stocks: list[NaturalCapitalStock],
    produced_capital: float = 850.0,    # £M
    financial_assets: float = 120.0,    # £M
    liabilities: float      = 380.0,    # £M
    gdp_index: float        = 1.0,
    investment_per_stock: float = 0.0,  # restoration investment (£M/yr each)
    T: int                  = 30,
) -> pd.DataFrame:
    """
    Simulate the PBS for T years.

    Returns a DataFrame with columns:
        year, PBS_total, NC_total, [stock_name for each NC stock], stewardship_met
    """
    pbs = ProvisioningBalanceSheet(produced_capital, financial_assets, liabilities, stocks)

    records = []
    for t in range(T + 1):
        row = {
            "year":           t,
            "PBS_total":      pbs.total,
            "NC_total":       pbs.nc_value,
            "stewardship":    pbs.stewardship_met,
        }
        for s in stocks:
            row[s.name] = s.current_level
        records.append(row)

        if t < T:
            # Step forward
            for s in stocks:
                s.step(gdp_index, investment=investment_per_stock)
            pbs.produced_capital *= (1 + 0.02)   # modest capital accumulation
            pbs.record()

    return pd.DataFrame(records)


def plot_pbs_comparison(T: int = 30) -> None:
    """Compare baseline (no restoration) vs. cooperative restoration investment."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    for (ax_idx, inv, label, color) in [
        (0, 0.0,  "Baseline (no restoration)",   "firebrick"),
        (1, 0.37, "Cooperative restoration (£11.1M/yr)", "steelblue"),
    ]:
        stocks = darlington_nc_stocks()
        df     = run_pbs_simulation(stocks, investment_per_stock=inv, T=T)

        axes[0, ax_idx].plot(df.year, df.PBS_total, color=color)
        axes[0, ax_idx].set_title(f"PBS Total — {label}")
        axes[0, ax_idx].set_xlabel("Year")
        axes[0, ax_idx].set_ylabel("£M")
        axes[0, ax_idx].grid(alpha=0.3)

        nc_cols = [c for c in df.columns if c not in ["year","PBS_total","NC_total","stewardship"]]
        for col in nc_cols:
            axes[1, ax_idx].plot(df.year, df[col], label=col[:25])
        axes[1, ax_idx].set_title(f"Natural Capital Levels — {label}")
        axes[1, ax_idx].set_xlabel("Year")
        axes[1, ax_idx].set_ylabel("Normalised stock level")
        axes[1, ax_idx].legend(fontsize=6)
        axes[1, ax_idx].axhline(0.5, color="grey", linestyle=":", linewidth=0.8,
                                label="Critical thresholds (approx)")
        axes[1, ax_idx].grid(alpha=0.3)

    plt.suptitle("Chapter 18: Provisioning Balance Sheet Simulation\n"
                 "Darlington / Tees Valley calibration (Chapter 44)", fontsize=11)
    plt.tight_layout()
    plt.savefig("ch18_pbs_darlington.png", dpi=150)
    plt.show()
    print("Saved: ch18_pbs_darlington.png")


if __name__ == "__main__":
    print("Chapter 18 — SFC-N: Provisioning Balance Sheet")
    print("=" * 60)

    stocks = darlington_nc_stocks()
    pbs    = ProvisioningBalanceSheet(850.0, 120.0, 380.0, stocks)

    print("\nInitial Provisioning Balance Sheet (Darlington):")
    print(pbs.report().to_string(index=False))

    print(f"\nStewardship audit (N̊_j ≥ 0 at baseline GDP?):")
    total_regen  = sum(s.annual_regeneration for s in stocks)
    total_deplet = sum(s.annual_depletion(1.0) for s in stocks)
    print(f"  Total annual regeneration : £{total_regen:.2f}M")
    print(f"  Total annual depletion    : £{total_deplet:.2f}M")
    print(f"  Net PBS change            : £{total_regen - total_deplet:.2f}M/yr")
    print(f"  Stewardship Condition met?: {total_regen >= total_deplet}")

    print("\nInvestment needed to meet Stewardship Condition:")
    gap = total_deplet - total_regen
    print(f"  Annual gap : £{gap:.2f}M/yr")

    print("\nRunning 30-year PBS simulation...")
    plot_pbs_comparison(T=30)
