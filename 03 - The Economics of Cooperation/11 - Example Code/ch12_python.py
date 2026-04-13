"""
Chapter 12: Network Structure and Economic Efficiency
======================================================
Implements:
  - Algebraic connectivity (Fiedler value) computation
  - Small-world vs. scale-free degree distribution comparison
  - Shock transmission eigenvalue analysis
  - The 2008 financial crisis network fragility model

Dependencies: numpy, scipy, networkx, matplotlib
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.linalg import eigh


def fiedler_value(G: nx.Graph) -> float:
    """Return the algebraic connectivity λ₂(L) of graph G."""
    L = nx.laplacian_matrix(G).toarray().astype(float)
    eigenvalues = np.sort(np.linalg.eigvalsh(L))
    return float(eigenvalues[1])


def spectral_gap(G: nx.Graph) -> float:
    """Return λ₂(L) / λ_n(L) — relative spectral gap."""
    L = nx.laplacian_matrix(G).toarray().astype(float)
    ev = np.sort(np.linalg.eigvalsh(L))
    return float(ev[1] / ev[-1]) if ev[-1] > 0 else 0.0


def network_stats(G: nx.Graph) -> dict:
    """Compute a standard battery of network statistics."""
    return {
        "n_nodes":          G.number_of_nodes(),
        "n_edges":          G.number_of_edges(),
        "density":          nx.density(G),
        "avg_degree":       np.mean([d for _, d in G.degree()]),
        "clustering":       nx.average_clustering(G),
        "avg_path_length":  nx.average_shortest_path_length(G) if nx.is_connected(G) else np.nan,
        "fiedler_value":    fiedler_value(G),
        "spectral_gap":     spectral_gap(G),
        "degree_cv":        np.std([d for _,d in G.degree()]) /
                            (np.mean([d for _,d in G.degree()]) + 1e-9),
    }


def compare_topologies(n: int = 200, k: int = 6, seed: int = 42) -> None:
    """
    Compare Watts-Strogatz, Barabási-Albert, and Erdős-Rényi topologies
    on economic network measures. Reproduces Table 12.1.
    """
    graphs = {
        "Watts-Strogatz (β=0.10)":  nx.watts_strogatz_graph(n, k, 0.10, seed=seed),
        "Watts-Strogatz (β=0.50)":  nx.watts_strogatz_graph(n, k, 0.50, seed=seed),
        "Barabási-Albert (m=3)":    nx.barabasi_albert_graph(n, 3, seed=seed),
        "Erdős-Rényi (p=k/n)":      nx.erdos_renyi_graph(n, k/n, seed=seed),
    }

    print(f"\n── Network Topology Comparison (n={n}) ─────────────────────")
    header = f"{'Topology':<28} {'λ₂':>8} {'Cluster':>9} {'AvgPath':>9} {'DegCV':>8}"
    print(header)
    print("-" * 70)
    for label, G in graphs.items():
        s = network_stats(G)
        print(f"{label:<28} {s['fiedler_value']:>8.3f} "
              f"{s['clustering']:>9.3f} {s['avg_path_length']:>9.3f} "
              f"{s['degree_cv']:>8.3f}")


def shock_transmission_analysis(n: int = 50, seed: int = 42) -> None:
    """
    Compare shock decay rates (max eigenvalue of A - γI) across topologies.
    Demonstrates the 2008 crisis fragility mechanism (Section 12.4).
    """
    gamma = 0.10  # natural shock absorption rate

    graphs = {
        "Scale-free (BA)":      nx.barabasi_albert_graph(n, 3, seed=seed),
        "Small-world (WS)":     nx.watts_strogatz_graph(n, 6, 0.10, seed=seed),
        "Cooperative (WS low)": nx.watts_strogatz_graph(n, 4, 0.05, seed=seed),
    }

    print("\n── Shock Transmission Analysis ─────────────────────────────")
    print(f"{'Topology':<25} {'ρ(A)':>8} {'ρ(A)-γ':>10} {'Stable?':>10}")
    print("-" * 58)
    for label, G in graphs.items():
        A = nx.to_numpy_array(G)
        A = A / A.sum(axis=1, keepdims=True)   # row-normalise
        rho = np.max(np.real(np.linalg.eigvals(A)))
        net_rate = rho - gamma
        stable = "YES" if net_rate < 0 else "NO (cascades)"
        print(f"{label:<25} {rho:>8.3f} {net_rate:>10.3f} {stable:>10}")


def plot_degree_distributions(n: int = 1000, seed: int = 42) -> None:
    """Plot degree distributions for scale-free vs. small-world networks."""
    G_ba = nx.barabasi_albert_graph(n, 3, seed=seed)
    G_ws = nx.watts_strogatz_graph(n, 6, 0.10, seed=seed)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    for ax, G, label, color in zip(
        axes,
        [G_ba, G_ws],
        ["Barabási-Albert (scale-free)", "Watts-Strogatz (small-world)"],
        ["firebrick", "steelblue"],
    ):
        degrees = [d for _, d in G.degree()]
        ax.hist(degrees, bins=30, color=color, edgecolor="white", alpha=0.8)
        ax.set_xlabel("Degree k")
        ax.set_ylabel("Count")
        ax.set_title(f"{label}\nλ₂ = {fiedler_value(G):.3f}")
        ax.grid(alpha=0.3)

    plt.suptitle("Degree distributions: scale-free vs. small-world\n"
                 "(Chapter 12, Figure 12.1)", fontsize=10)
    plt.tight_layout()
    plt.savefig("ch12_degree_distributions.png", dpi=150)
    plt.show()
    print("Saved: ch12_degree_distributions.png")


if __name__ == "__main__":
    print("Chapter 12 — Network Structure and Economic Efficiency")
    print("=" * 60)
    compare_topologies()
    shock_transmission_analysis()
    plot_degree_distributions()
