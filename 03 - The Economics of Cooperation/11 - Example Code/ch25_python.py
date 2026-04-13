"""
Chapter 25: Mutual Credit and Clearing Systems
===============================================
Algorithm 25.1: Optimal Multilateral Clearing via NetworkX min-cost flow.
Includes the 50-firm B2B worked example.

Dependencies: numpy, networkx, matplotlib, pandas
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from typing import Optional


def optimal_clearing(
    members: list,
    obligations: dict,          # {(i, j): amount}
    balances: dict,             # {i: balance}  sum must be 0
    costs: Optional[dict] = None,
) -> dict:
    """
    Algorithm 25.1 — Optimal Multilateral Clearing.

    Solves the minimum-cost flow LP: clear as many obligations as
    possible while minimising total transaction cost.

    Parameters
    ----------
    members     : list of member IDs
    obligations : dict mapping (seller, buyer) -> gross obligation amount
    balances    : dict mapping member -> net balance (positive=creditor)
    costs       : optional dict mapping (i, j) -> unit clearing cost
                  (defaults to zero for all edges — free digital clearing)

    Returns
    -------
    dict with keys:
        flows              : {(i,j): flow} — clearing flows
        residual_balances  : {i: balance} after clearing
        clearing_efficiency: fraction of gross obligations cleared
        gross_cleared      : total volume cleared
        total_gross        : total gross obligations
    """
    if costs is None:
        costs = {}

    G = nx.DiGraph()

    # Add nodes with demand = -balance
    # NetworkX convention: positive demand = net consumer (needs inflow)
    # Creditors (balance > 0) need outflow  → negative demand
    # Debtors   (balance < 0) need inflow   → positive demand
    for m in members:
        G.add_node(m, demand=-balances.get(m, 0))

    # Add edges
    total_gross = 0.0
    for (i, j), amount in obligations.items():
        if amount > 1e-9:
            G.add_edge(
                i, j,
                capacity=amount,
                weight=costs.get((i, j), 0),
            )
            total_gross += amount

    # Solve
    try:
        flow_dict = nx.min_cost_flow(G)
    except nx.NetworkXUnfeasible:
        # Fall back to partial clearing
        return _partial_clearing(members, obligations, balances, total_gross)

    # Extract flows
    flows = {}
    cleared_volume = 0.0
    for i in G.nodes():
        for j, f in flow_dict.get(i, {}).items():
            if f > 1e-9:
                flows[(i, j)] = f
                cleared_volume += f

    # Residual balances
    residual = dict(balances)
    for (i, j), f in flows.items():
        residual[i] = residual.get(i, 0) + f      # debit settled
        residual[j] = residual.get(j, 0) - f      # credit settled

    total_abs_balance = sum(abs(b) for b in balances.values())
    residual_abs      = sum(abs(b) for b in residual.values())
    efficiency = 1 - residual_abs / (total_abs_balance + 1e-9)

    return {
        "flows":               flows,
        "residual_balances":   residual,
        "clearing_efficiency": efficiency,
        "gross_cleared":       cleared_volume,
        "total_gross":         total_gross,
    }


def _partial_clearing(members, obligations, balances, total_gross):
    """Fallback: greedy bilateral netting when full clearing is infeasible."""
    residual = dict(balances)
    flows = {}
    cleared = 0.0
    for (i, j), amount in sorted(obligations.items(), key=lambda x: -x[1]):
        flow = min(amount, max(0, -residual.get(i, 0)), max(0, residual.get(j, 0)))
        if flow > 1e-9:
            flows[(i, j)] = flow
            residual[i] = residual.get(i, 0) + flow
            residual[j] = residual.get(j, 0) - flow
            cleared += flow
    total_abs = sum(abs(b) for b in balances.values())
    residual_abs = sum(abs(b) for b in residual.values())
    return {
        "flows": flows,
        "residual_balances": residual,
        "clearing_efficiency": 1 - residual_abs / (total_abs + 1e-9),
        "gross_cleared": cleared,
        "total_gross": total_gross,
    }


# ─────────────────────────────────────────────────────────────────────────────
# 50-FIRM B2B WORKED EXAMPLE  (Section 25.8)
# ─────────────────────────────────────────────────────────────────────────────

def generate_b2b_network(n: int = 50, seed: int = 42) -> tuple:
    """
    Generate a synthetic 50-firm B2B mutual credit network.
    Uses a Watts-Strogatz small-world topology (k=8, β=0.10).

    Returns
    -------
    (members, obligations, balances, G)
    """
    rng = np.random.default_rng(seed)
    members = list(range(n))

    # Small-world trade network
    G = nx.watts_strogatz_graph(n, 8, 0.10, seed=seed)

    # Assign random gross trade obligations on each directed edge
    obligations = {}
    for (i, j) in G.edges():
        amt_ij = rng.uniform(5_000, 40_000)
        amt_ji = rng.uniform(5_000, 40_000)
        obligations[(i, j)] = round(amt_ij, 2)
        obligations[(j, i)] = round(amt_ji, 2)

    # Compute net balances (must sum to 0 in mutual credit)
    raw_balances = {}
    for m in members:
        credits = sum(v for (i, j), v in obligations.items() if j == m)
        debits  = sum(v for (i, j), v in obligations.items() if i == m)
        raw_balances[m] = round(credits - debits, 2)

    # Force exact zero-sum (floating point correction)
    total = sum(raw_balances.values())
    raw_balances[0] -= total
    balances = raw_balances

    return members, obligations, balances, G


def run_b2b_example() -> None:
    members, obligations, balances, G = generate_b2b_network(n=50)

    print("\n── 50-Firm B2B Mutual Credit Network (Section 25.8) ────────")
    total_gross = sum(obligations.values())
    print(f"  Firms              : {len(members)}")
    print(f"  Trade edges        : {len(obligations)} directed")
    print(f"  Total gross (€)    : {total_gross:,.0f}")
    print(f"  Max creditor (€)   : {max(balances.values()):,.0f}")
    print(f"  Max debtor   (€)   : {min(balances.values()):,.0f}")
    print(f"  System sum  (≡0?)  : {sum(balances.values()):.4f}")

    result = optimal_clearing(members, obligations, balances)

    print(f"\n  Gross obligations cleared: {result['gross_cleared']:,.0f} "
          f"({result['clearing_efficiency']*100:.1f}%)")
    print(f"  Remaining (uncleared)    : "
          f"{result['total_gross'] - result['gross_cleared']:,.0f}")

    residual = result["residual_balances"]
    max_residual = max(abs(v) for v in residual.values())
    print(f"  Max residual balance (€) : {max_residual:,.0f}")

    # Bilateral netting benchmark
    print("\n  Bilateral netting comparison:")
    bilateral_cleared = 0.0
    for (i, j), amount in obligations.items():
        reverse = obligations.get((j, i), 0)
        if amount > reverse:
            bilateral_cleared += reverse
        else:
            bilateral_cleared += amount
    bilateral_cleared /= 2   # avoid double counting
    bilateral_eff = bilateral_cleared / total_gross
    print(f"  Bilateral clearing efficiency: {bilateral_eff*100:.1f}%")
    print(f"  Multilateral additional gain : "
          f"{(result['clearing_efficiency'] - bilateral_eff)*100:.1f}pp")


def plot_clearing_network(n: int = 15, seed: int = 7) -> None:
    """Visualise the clearing flows on a small example network."""
    members, obligations, balances, G = generate_b2b_network(n=n, seed=seed)
    result = optimal_clearing(members, obligations, balances)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    pos = nx.spring_layout(G, seed=seed)

    # Pre-clearing
    colors = ["steelblue" if balances[m] >= 0 else "firebrick" for m in members]
    sizes  = [300 + abs(balances[m]) / 100 for m in members]
    nx.draw_networkx(G, pos, ax=axes[0], node_color=colors, node_size=sizes,
                     with_labels=True, font_size=7, edge_color="grey", alpha=0.7)
    axes[0].set_title("Pre-clearing network\n(blue=creditor, red=debtor)")

    # Post-clearing flows
    flow_edges = [(i, j) for (i, j) in result["flows"] if result["flows"][(i,j)] > 1e-3]
    flow_widths = [result["flows"][(i,j)] / 2000 for (i,j) in flow_edges]
    res_colors = ["steelblue" if result["residual_balances"][m] >= 0 else "firebrick"
                  for m in members]
    nx.draw_networkx(G, pos, ax=axes[1], node_color=res_colors, node_size=sizes,
                     with_labels=True, font_size=7, edge_color="grey", alpha=0.3)
    nx.draw_networkx_edges(G, pos, edgelist=flow_edges, width=flow_widths,
                           edge_color="darkorange", alpha=0.8, ax=axes[1],
                           arrows=True, arrowsize=15)
    axes[1].set_title(f"Post-clearing flows (efficiency={result['clearing_efficiency']*100:.0f}%)")

    plt.suptitle(f"Chapter 25: Multilateral Clearing (n={n} firms)", fontsize=11)
    plt.tight_layout()
    plt.savefig("ch25_mutual_credit_clearing.png", dpi=150)
    plt.show()
    print("Saved: ch25_mutual_credit_clearing.png")


if __name__ == "__main__":
    print("Chapter 25 — Mutual Credit and Clearing Systems")
    print("=" * 60)
    run_b2b_example()
    plot_clearing_network(n=15)
