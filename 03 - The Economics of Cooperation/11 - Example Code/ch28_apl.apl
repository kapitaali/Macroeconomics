⍝ ══════════════════════════════════════════════════════════════════
⍝ Chapter 28: SFC Monetary Regimes — APL
⍝ SFC balance sheet row-sum verification and Minsky condition
⍝ ══════════════════════════════════════════════════════════════════

⍝ ── 1. SFC CONSISTENCY CHECK ────────────────────────────────────────
⍝ Every row of BSM must sum to zero (every asset is someone's liability)
⍝ Every column must equal sector net worth

SFCRowCheck ← { ∧/ 1e¯8 > |+/⍵ }     ⍝ 1 iff all rows sum to ~0
SFCColCheck ← { ∧/ 1e¯8 > |+⌿⍵ }     ⍝ 1 iff all cols sum to ~0 (net worth in last row)

⍝ ── 2. MINSKY STABILITY CONDITION ──────────────────────────────────
⍝ Minsky stable when i - π < g  (Theorem 23.1)
⍝ ⍵ = (i, pi, g)
MinskyStable ← {
    i pi g ← ⍵
    (i-pi) < g
}

⍝ ── 3. VELOCITY UNDER DEMURRAGE (Theorem 27.1) ──────────────────────
⍝ V(δ) = V₀ × √(1 + δ/δ₀)   where δ₀ is baseline demurrage level
⍝ ⍺=(V0 delta0), ⍵=delta
DemurrageVelocity ← {
    V0 d0 ← ⍺
    d ← ⍵
    V0 × (1+d÷d0)*0.5
}

⍝ ── 4. COMPARATIVE REGIME SUMMARY MATRIX ───────────────────────────
⍝ Rows: regimes (debt, sovereign, mutual, demurrage)
⍝ Cols: (velocity, investment_rate, minsky_stable, seigniorage_public)

RegimeMatrix ← 4 4 ⍴
    1.35 0.182 0 0     ⍝ Debt-based
    1.50 0.204 1 1     ⍝ Sovereign
    1.65 0.221 1 0     ⍝ Mutual credit
    1.76 0.233 1 1     ⍝ Demurrage

⍝ Which regimes are Minsky-stable?
⍝   RegimeMatrix[;3]     ⍝ → 0 1 1 1  (debt-based fails)

⍝ ── 5. OPTIMAL HYBRID ALLOCATION (Proposition 28.1) ─────────────────
⍝ Effective velocity of a 60/25/15 hybrid monetary system

HybridVelocity ← {     ⍝ ⍵ = (shares, velocities) — two n-vectors
    shares vels ← ⍵
    shares +.× vels
}

⍝ Chapter 28 calibration:
shares_hybrid ← 0.60 0.25 0.15
vels_hybrid   ← 1.50 1.65 1.76
⍝ HybridVelocity shares_hybrid vels_hybrid  ⍝ → 1.566

⍝ ══════════════════════════════════════════════════════════════════
⍝ END CH28 APL
⍝ Key checks:
⍝   MinskyStable 0.05 0.02 0.03   ⍝ → 0 (unstable: i-π=0.03 = g)
⍝   MinskyStable 0.04 0.02 0.035  ⍝ → 1 (stable)
⍝   HybridVelocity shares_hybrid vels_hybrid  ⍝ → 1.566
⍝ ══════════════════════════════════════════════════════════════════
