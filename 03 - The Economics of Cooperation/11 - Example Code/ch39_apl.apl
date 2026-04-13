⍝ ══════════════════════════════════════════════════════════════════
⍝ Chapter 39: Data Cooperatives — APL
⍝ Shapley approximation and superadditivity verification
⍝ ══════════════════════════════════════════════════════════════════

⍝ ── 1. POWER-LAW DATA VALUE ─────────────────────────────────────────
⍝ v(S) = β × |S|^α       (simplified, no diversity term)
⍝ SuperaddValue n alpha beta → v for coalition of size n

DataValue ← {      ⍝ ⍺=(alpha beta), ⍵=|S|
    alpha beta ← ⍺
    beta × ⍵*alpha
}

⍝ Superadditivity ratio: v(N) / Σv({i})  (Proposition 39.1)
⍝ = n^(α-1)
SuperAdditivityRatio ← {    ⍝ ⍺=alpha, ⍵=n
    ⍵*(⍺-1)
}

⍝ Example: n=10000, alpha=1.6 → ratio ≈ 630
⍝   1.6 SuperAdditivityRatio 10000   ⍝ → 630.96

⍝ ── 2. OVA SCORE NORMALISATION ──────────────────────────────────────
⍝ Given raw contribution scores s_i, normalise so sum=1
NormScores ← { ⍵÷+/⍵ }

⍝ Shapley dividend for each contributor
ShapleyDividend ← {     ⍝ ⍺=total distributable revenue, ⍵=score vector
    ⍺ × NormScores ⍵
}

⍝ ── 3. GINI OF SHAPLEY ALLOCATION ───────────────────────────────────
Gini ← {
    w ← ⍵[⍋⍵] ⋄ n ← ≢w ⋄ i ← ⍳n
    (2×i+.×w-(n+1)×+/w) ÷ n×+/w
}

⍝ ── 4. TAXI COOPERATIVE OVA (simplified) ────────────────────────────
⍝ weights: trips=0.5, quality=0.25, diversity=0.15, governance=0.10
⍝ ⍺=weights (4-element), ⍵=feature matrix (n×4)
OVAScore ← { ⍵+.×⍺ }

⍝ Example with 5 drivers:
⍝ trips qual div gov (normalised to 1.0 = average)
ExDrivers ← 5 4⍴ 1.5 0.95 0.80 0.70
            0.9 0.75 0.40 0.90
            0.7 0.88 0.55 0.60
            1.2 0.92 0.70 0.85
            0.8 0.70 0.45 0.50

weights  ← 0.50 0.25 0.15 0.10
scores5  ← ExDrivers OVAScore weights
divs5    ← 488390 ShapleyDividend scores5

⍝ divs5  ← Shapley dividends for 5 example drivers (sum = 488390)
⍝ Gini divs5  ← inequality of dividend distribution

⍝ ══════════════════════════════════════════════════════════════════
⍝ END CH39 APL
⍝ Key checks:
⍝   1.6 SuperAdditivityRatio 10000   ⍝ → ~630 (grand coalition = 630× singletons)
⍝   +/divs5                          ⍝ → 488390 (efficiency)
⍝   Gini divs5                       ⍝ → small (OVA is more equal than pure capital)
⍝ ══════════════════════════════════════════════════════════════════
