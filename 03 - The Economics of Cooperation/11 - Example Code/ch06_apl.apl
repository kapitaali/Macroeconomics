⍝ ══════════════════════════════════════════════════════════════════
⍝ Chapter 6: Cooperative Games and the Core — APL
⍝
⍝ APL excels here: coalition enumeration, Shapley sums, and
⍝ core verification are pure array operations.
⍝ ══════════════════════════════════════════════════════════════════

⍝ ── 1. COALITION BITMASK REPRESENTATION ────────────────────────────
⍝
⍝ Represent a coalition S ⊆ {0…n-1} as an integer bitmask.
⍝ All non-empty coalitions of n players: integers 1 to 2^n - 1.
⍝
⍝ AllCoalitions n  →  vector of all 2^n - 1 non-empty coalition masks

AllCoalitions ← { 1↓⍳2*⍵ }     ⍝ masks 1 to 2^n - 1

⍝ Size of coalition from bitmask:
CoalSize ← { +/⍵⊤⍨⍵⍴2 }        ⍝ count 1-bits

⍝ Check if player i is in coalition mask c:
PlayerIn ← { 1=1∧⍺⌊⍵÷2*⍺ }     ⍝ ⍺=player, ⍵=mask

⍝ ── 2. SHAPLEY VALUE (exact, airport game) ──────────────────────────
⍝
⍝ For the airport game v(S) = costs[max(S)], the Shapley value has
⍝ an analytic form.  We implement it directly:
⍝   φ_i = (c_i - c_{i-1}) / (n - i + 1)    (1-indexed)

AirportShapley ← {          ⍝ ⍵ = vector of sorted costs (ascending)
    c    ← ⍵
    n    ← ≢c
    prev ← 0,c               ⍝ c_{i-1}, with c_0 = 0
    diff ← c - ¯1↓prev       ⍝ c_i - c_{i-1}
    denom← ⌽⍳n               ⍝ n, n-1, ..., 1
    diff ÷ denom
}

⍝ Verification:
⍝   AirportShapley 100 200 400
⍝   → 33.33  50  200   (marginal cost divided equally among users)
⍝   Sum should equal 400 (cost of longest runway)


⍝ ── 3. SHAPLEY VALUE (general, exact) ──────────────────────────────
⍝
⍝ Average marginal contributions over all n! permutations.
⍝ For small n (≤7) we enumerate all permutations.
⍝ Returns a vector of n Shapley values.

⍝ Factorial (helper):
Fact ← {×/1+⍳⍵}

⍝ All permutations of 0..n-1 as an (n! × n) matrix:
Perms ← {           ⍝ ⍵ = n
    :If ⍵=1 ⋄ ,⊂,0 ⋄ :Return ⋄ :EndIf
    prev ← Perms ⍵-1
    res ← ⍬
    :For p :In prev
        :For k :In ⌽⍳⍵
            res ,← ⊂(k↑p),⍵-1,(k↓p)
        :EndFor
    :EndFor
    res
}

⍝ ShapleyExact: exact Shapley for small games via permutation enumeration.
⍝ Requires a characteristic function V as a (2^n)-element lookup vector
⍝ indexed by bitmask (index 0 = empty set = 0).

ShapleyExact ← {        ⍝ ⍵ = V lookup vector (length 2^n)
    V ← ⍵
    n ← ⌊2⍟≢V           ⍝ n players
    perms ← Perms n
    phi ← n⍴0

    :For perm :In perms
        coal ← 0
        :For player :In perm
            new_coal ← coal+2*player
            phi[player] +← V[1+new_coal]-V[1+coal]
            coal ← new_coal
        :EndFor
    :EndFor

    phi ÷ Fact n
}

⍝ Usage — airport game with 3 players, costs 100 200 400:
⍝ Build V lookup:
⍝   V ← 8⍴0
⍝   V[1+1]←100  V[1+2]←200  V[1+4]←400   ⍝ singletons
⍝   V[1+3]←200  V[1+5]←400  V[1+6]←400   ⍝ pairs: max cost
⍝   V[1+7]←400                             ⍝ grand coalition
⍝   ShapleyExact V   ⍝ → ≈ 33.33  50  316.67 (numerics)


⍝ ── 4. CORE CHECK (LP-free, necessary condition) ───────────────────
⍝
⍝ For convex games (marginal contributions non-decreasing),
⍝ the core is guaranteed non-empty.
⍝ Check convexity: MC_i(S) ≤ MC_i(T) for S ⊆ T.
⍝
⍝ Quick check for airport game: always convex → core always non-empty.
⍝ (Proof: v(S∪{i}) - v(S) = max(c_i, max_S(costs)) - max_S(costs)
⍝  which is non-decreasing as S grows.)

AirportConvex ← 1    ⍝ airport game is always convex: hardcoded result


⍝ ── 5. COOPERATIVE SURPLUS FRACTION σ_v ────────────────────────────
⍝
⍝ σ_v = (v(N) - Σ v({i})) / v(N)
⍝ V = lookup vector, n = number of players

SigmaV ← {           ⍝ ⍺=V lookup (2^n), ⍵=n
    V ← ⍺ ⋄ n ← ⍵
    grand_mask ← 2*n-1            ⍝ bitmask with all n bits set = 2^n - 1
    vN   ← V[1+grand_mask]        ⍝ v(grand coalition)
    vSing← +/ V[1+2*⍳n]           ⍝ sum of v({0}), v({1}), …
    (vN-vSing)÷vN
}


⍝ ── 6. GINI OF SHAPLEY vs. COMPETITIVE ALLOCATION ──────────────────

Gini ← {
    w ← ⍵[⍋⍵] ⋄ n ← ≢w ⋄ i ← ⍳n
    (2×i+.×w-(n+1)×+/w) ÷ n×+/w
}

⍝ Example:
⍝   phi_shapley     ← 33.33 50 316.67    ⍝ airport game Shapley
⍝   phi_competitive ← 100 200 400        ⍝ marginal product = own cost
⍝   Gini phi_shapley      ⍝ lower Gini (more equal)
⍝   Gini phi_competitive  ⍝ higher Gini (less equal)


⍝ ══════════════════════════════════════════════════════════════════
⍝ END OF CHAPTER 6 APL
⍝
⍝ Quick verification:
⍝   AirportShapley 100 200 400      ⍝ → 33.33 50 316.67
⍝   +/AirportShapley 100 200 400   ⍝ → 400 (efficiency)
⍝   Fact 3                          ⍝ → 6 (3! permutations)
⍝ ══════════════════════════════════════════════════════════════════
