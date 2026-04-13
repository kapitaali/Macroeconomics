"""
Chapter 6: Cooperative Games and the Core
==========================================
Implements:
  - Characteristic function games
  - Core non-emptiness check (Bondareva-Shapley LP)
  - Shapley value (exact for n <= 10, Monte Carlo for larger n)
  - Nucleolus (via sequential LP)
  - Airport cost allocation worked example

Dependencies: numpy, scipy
"""

import numpy as np
from itertools import permutations, combinations
from scipy.optimize import linprog
from functools import lru_cache
from typing import Callable


# ─────────────────────────────────────────────────────────────────────────────
# 1.  CHARACTERISTIC FUNCTION HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def all_coalitions(n: int) -> list[frozenset]:
    """Return all non-empty subsets of {0, 1, ..., n-1}."""
    players = list(range(n))
    return [
        frozenset(s)
        for r in range(1, n + 1)
        for s in combinations(players, r)
    ]


def coalition_to_bitmask(S: frozenset, n: int) -> int:
    mask = 0
    for i in S:
        mask |= (1 << i)
    return mask


# ─────────────────────────────────────────────────────────────────────────────
# 2.  CORE NON-EMPTINESS  (Bondareva-Shapley LP, Theorem 6.1)
# ─────────────────────────────────────────────────────────────────────────────

def core_nonempty(n: int, v: Callable[[frozenset], float]) -> tuple[bool, np.ndarray | None]:
    """
    Test whether the core of (N, v) is non-empty using linear programming.

    The core is the set of payoff vectors x with:
        sum(x) = v(N)              (efficiency)
        sum(x_i for i in S) >= v(S)  for all S  (group rationality)

    We solve the LP:
        minimise   0  (feasibility check)
        s.t.       Ax >= b   (group rationality for every coalition)
                   1^T x = v(N)  (efficiency)

    Parameters
    ----------
    n : number of players
    v : characteristic function  v(S) -> float

    Returns
    -------
    (feasible, x) where x is a core imputation if feasible, else None
    """
    players = frozenset(range(n))
    grand   = v(players)
    coalitions = all_coalitions(n)
    m = len(coalitions)

    # Build constraint matrix A_ub x >= b_ub   (scipy uses <= so negate)
    A_ub = np.zeros((m, n))
    b_ub = np.zeros(m)
    for row, S in enumerate(coalitions):
        for i in S:
            A_ub[row, i] = -1.0     # negated for <=
        b_ub[row] = -v(S)

    # Efficiency equality constraint
    A_eq = np.ones((1, n))
    b_eq = np.array([grand])

    result = linprog(
        c       = np.zeros(n),
        A_ub    = A_ub,
        b_ub    = b_ub,
        A_eq    = A_eq,
        b_eq    = b_eq,
        method  = "highs",
    )

    if result.status == 0:
        return True, result.x
    return False, None


# ─────────────────────────────────────────────────────────────────────────────
# 3.  SHAPLEY VALUE  (exact)
# ─────────────────────────────────────────────────────────────────────────────

def shapley_exact(n: int, v: Callable[[frozenset], float]) -> np.ndarray:
    """
    Compute exact Shapley value by averaging marginal contributions
    over all n! permutations.

    Time complexity: O(n! * n)  — practical for n <= 10.

    Returns
    -------
    phi : array of shape (n,)
    """
    players = list(range(n))
    phi     = np.zeros(n)

    for perm in permutations(players):
        coalition = frozenset()
        for i, player in enumerate(perm):
            new_coalition = coalition | {player}
            marginal = v(new_coalition) - v(coalition)
            phi[player] += marginal
            coalition = new_coalition

    phi /= np.math.factorial(n)
    return phi


def shapley_monte_carlo(
    n: int,
    v: Callable[[frozenset], float],
    n_samples: int = 10_000,
    seed: int = 42,
) -> np.ndarray:
    """
    Monte Carlo Shapley approximation — O(n_samples * n).
    Use when n > 10.
    """
    rng     = np.random.default_rng(seed)
    players = np.arange(n)
    phi     = np.zeros(n)

    for _ in range(n_samples):
        perm      = rng.permutation(players)
        coalition = frozenset()
        for player in perm:
            new_coalition = coalition | {player}
            phi[player] += v(new_coalition) - v(coalition)
            coalition = new_coalition

    return phi / n_samples


# ─────────────────────────────────────────────────────────────────────────────
# 4.  NUCLEOLUS  (sequential LP method, small games)
# ─────────────────────────────────────────────────────────────────────────────

def nucleolus(n: int, v: Callable[[frozenset], float], tol: float = 1e-8) -> np.ndarray:
    """
    Compute the nucleolus via the sequential LP (Maschler) method.

    At each stage we find the minimum epsilon such that
        sum(x_i for i in S) >= v(S) - epsilon
    for all remaining (not yet tight) coalitions.

    This is a simplified implementation suitable for n <= 6.

    Returns
    -------
    x : nucleolus imputation, array of shape (n,)
    """
    players   = frozenset(range(n))
    grand     = v(players)
    coalitions = [S for S in all_coalitions(n) if S != players]
    m = len(coalitions)

    x       = np.zeros(n)
    fixed   = {}        # player index -> fixed value
    active  = set(range(m))

    for _iteration in range(m):
        if not active:
            break

        active_list = list(active)
        # Variables: x[0..n-1], epsilon
        n_vars = n + 1

        # Minimise epsilon (last variable)
        c = np.zeros(n_vars)
        c[-1] = 1.0

        A_ub = []
        b_ub_vals = []
        for row in active_list:
            S = coalitions[row]
            a = np.zeros(n_vars)
            for i in S:
                a[i] = -1.0
            a[-1] = -1.0          # - epsilon
            A_ub.append(a)
            b_ub_vals.append(-v(S))

        A_ub_arr = np.array(A_ub)
        b_ub_arr = np.array(b_ub_vals)

        A_eq = [np.zeros(n_vars)]
        A_eq[0][:n] = 1.0
        b_eq = [grand]

        # Fix already-determined players
        bounds = [(None, None)] * n + [(None, None)]
        for pi, val in fixed.items():
            bounds[pi] = (val, val)

        res = linprog(
            c, A_ub=A_ub_arr, b_ub=b_ub_arr,
            A_eq=np.array(A_eq), b_eq=np.array(b_eq),
            bounds=bounds, method="highs",
        )

        if res.status != 0:
            break

        x[:] = res.x[:n]
        eps_star = res.x[-1]

        # Identify tight coalitions (excess = eps_star)
        tight = []
        for row in active_list:
            S = coalitions[row]
            excess = sum(x[i] for i in S) - v(S)
            if abs(excess - eps_star) < tol:
                tight.append(row)
        for row in tight:
            active.discard(row)

    return x


# ─────────────────────────────────────────────────────────────────────────────
# 5.  AIRPORT COST ALLOCATION  (Worked Example, Section 6.4)
# ─────────────────────────────────────────────────────────────────────────────

def airport_game(costs: list[float]) -> Callable[[frozenset], float]:
    """
    Construct the airport cost-sharing game characteristic function.

    Players are aircraft types sorted by required runway length.
    The cost of serving coalition S is the cost of the longest runway needed.

    Parameters
    ----------
    costs : list of runway costs c_1 <= c_2 <= ... <= c_n
            (one per aircraft type / player)

    Returns
    -------
    Callable v(S) -> float
    """
    costs = sorted(costs)

    def v(S: frozenset) -> float:
        if not S:
            return 0.0
        return costs[max(S)]

    return v


def demo_airport():
    """
    Three-airport example from the worked example in Section 6.4.
    Aircraft types: small (c=100), medium (c=200), large (c=400)
    """
    print("\n── Airport Cost Allocation (Section 6.4) ──────────────────")
    costs = [100.0, 200.0, 400.0]   # runway costs, ascending
    n     = len(costs)
    v     = airport_game(costs)

    # Verify superadditivity
    print("\nCharacteristic function values:")
    for S in all_coalitions(n):
        print(f"  v({set(S)}) = {v(S):.0f}")

    # Core
    feasible, x_core = core_nonempty(n, v)
    print(f"\nCore non-empty: {feasible}")
    if feasible:
        print(f"Core imputation (LP): {np.round(x_core, 2)}")

    # Shapley value
    phi = shapley_exact(n, v)
    print(f"\nShapley value: {np.round(phi, 3)}")
    print(f"  Sum = {phi.sum():.1f}  (should equal v(N) = {v(frozenset(range(n)))})")

    # Analytical Shapley for airport game: phi_i = (c_i - c_{i-1}) / (n - i + 1)
    phi_analytic = np.zeros(n)
    c_prev = 0.0
    for i, c in enumerate(costs):
        phi_analytic[i] = (c - c_prev) / (n - i)
        c_prev = c
    print(f"Analytical Shapley: {np.round(phi_analytic, 3)}")

    # Nucleolus
    nucl = nucleolus(n, v)
    print(f"\nNucleolus: {np.round(nucl, 3)}")

    print("\nInterpretation:")
    labels = ["Small aircraft", "Medium aircraft", "Large aircraft"]
    for i, label in enumerate(labels):
        print(f"  {label}: Shapley = {phi[i]:.1f}  |  Nucleolus = {nucl[i]:.1f}")


# ─────────────────────────────────────────────────────────────────────────────
# 6.  COOPERATIVE vs. COMPETITIVE WELFARE COMPARISON
# ─────────────────────────────────────────────────────────────────────────────

def welfare_comparison(n: int, v: Callable[[frozenset], float]) -> dict:
    """
    Compare: grand coalition (cooperative) vs. singleton coalitions (competitive).

    Returns dict with:
      - coop_welfare    : v(grand coalition)
      - comp_welfare    : sum of v({i}) for all i
      - coop_surplus    : difference
      - sigma_v         : cooperative surplus fraction
      - shapley_phi     : Shapley allocation
    """
    players     = frozenset(range(n))
    coop_welfare = v(players)
    comp_welfare = sum(v(frozenset([i])) for i in range(n))
    surplus      = coop_welfare - comp_welfare
    sigma_v      = surplus / coop_welfare if coop_welfare > 0 else 0.0
    phi          = shapley_exact(n, v)

    return {
        "coop_welfare": coop_welfare,
        "comp_welfare": comp_welfare,
        "coop_surplus": surplus,
        "sigma_v":      sigma_v,
        "shapley_phi":  phi,
    }


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("Chapter 6 — Cooperative Games and the Core")
    print("=" * 60)

    demo_airport()

    print("\n── Simple superadditive game ───────────────────────────────")
    def v_simple(S):
        vals = {1: 10, 2: 15, 3: 12}
        sizes= {1: 10, 2: 15, 3: 12, 3: 40}  # grand coalition synergy
        if len(S) == 3: return 40.0
        if len(S) == 2: return sum(vals[i+1] for i in S) * 1.2
        return vals[list(S)[0] + 1]

    def v_game(S):
        if not S: return 0.0
        if S == frozenset([0]): return 10.0
        if S == frozenset([1]): return 15.0
        if S == frozenset([2]): return 12.0
        if S == frozenset([0,1]): return 30.0
        if S == frozenset([0,2]): return 26.0
        if S == frozenset([1,2]): return 32.0
        return 48.0   # grand coalition

    wc = welfare_comparison(3, v_game)
    print(f"  Cooperative welfare (grand coalition): {wc['coop_welfare']:.1f}")
    print(f"  Competitive welfare (singletons sum) : {wc['comp_welfare']:.1f}")
    print(f"  Cooperative surplus σ_v              : {wc['sigma_v']:.3f}")
    print(f"  Shapley allocation                   : {np.round(wc['shapley_phi'],2)}")

    feasible, x = core_nonempty(3, v_game)
    print(f"  Core non-empty: {feasible}  |  Core point: {np.round(x,2)}")
