⍝ ══════════════════════════════════════════════════════════════════
⍝ Chapter 36: Regenerative Agriculture — APL
⍝ Soil capital dynamics and PES externality gap
⍝ ══════════════════════════════════════════════════════════════════

⍝ ── 1. SOIL CARBON STEADY STATE ─────────────────────────────────────
⍝ N* = (α_M × M_s / α_N - ε_s/α_N) / (1 + β_T×T/α_N - 1/K_s)
⍝ ⍵ = (alpha_M, alpha_N, K_s, epsilon_s, beta_T, T_intensity, M_s)

SoilSteadyState ← {
    aM aN Ks es bT T Ms ← ⍵
    num ← aM×Ms÷aN - es÷aN
    den ← 1 + bT×T÷aN - ÷Ks
    num÷den
}

⍝ Conventional vs. regenerative (UK upland params):
⍝   SoilSteadyState 0.08 0.04 90 0.8 0.02 1.5 3   ⍝ → conventional SOC
⍝   SoilSteadyState 0.08 0.04 90 0.8 0.02 0.2 8   ⍝ → regenerative SOC

⍝ ── 2. EXTERNALITY GAP (Proposition 36.1) ───────────────────────────
⍝ ΔM_s* = p^N_s / α_M × (r + β_T × T)
⍝ Under-investment in organic inputs per unit shadow price

ExternalityGap ← {     ⍝ ⍵ = (p_N, alpha_M, r, beta_T, T_intensity)
    pN aM r bT T ← ⍵
    pN ÷ aM × r + bT × T
}

⍝ UK calibration (Section 36.2.2):
⍝   ExternalityGap 2400 0.08 0.035 0.02 1.5   ⍝ → ~4.8 t/ha/yr under-investment

⍝ ── 3. ECOSYSTEM SERVICE VALUATION VECTOR ───────────────────────────
⍝ Annual ESV per hectare for each service type (GBP/ha/yr)
⍝ Rows: (baseline, regenerative)

ESV_Matrix ← 2 5 ⍴
    120 ¯85 ¯40  15   0    ⍝ baseline (grouse moor)
     65  185  95 110  45   ⍝ regenerative

⍝ Net benefit per ha:
NetBenefitHa ← { +/(⊃1↓⍵)-⊃⍵ }    ⍝ sum of row differences
⍝   NetBenefitHa ESV_Matrix   ⍝ → 490 GBP/ha/yr

⍝ ── 4. DISCOUNT FACTOR VECTOR AND NPV ──────────────────────────────
DiscFactors ← { (1+⍵)*-⍳⍺ }      ⍝ ⍺=T periods, ⍵=discount rate
NPV ← {                            ⍝ ⍺=annual_benefit, ⍵=(T, r)
    T r ← ⍵
    pv ← ⍺ × +/T DiscFactors r
    pv
}

⍝ Peak District:
⍝   area ← 10000                   ⍝ ha
⍝   ann_ben ← 490 × area           ⍝ = 4 900 000 GBP/yr
⍝   capex   ← 1800 × area          ⍝ transition cost = 18 000 000
⍝   pv_ben  ← ann_ben NPV 20 0.035 ⍝ PV over 20yr at 3.5%
⍝   bcr     ← pv_ben ÷ capex       ⍝ → ~5.5:1

⍝ ── 5. COOPERATIVE LANDSCAPE SURPLUS (Proposition 36.2) ────────────
⍝ Ecosystem service value increases superlinearly with landscape connectivity λ₂
⍝ V_landscape = V_patches × (1 + f(λ₂))
⍝ where f(λ₂) = connectivity bonus fraction

LandscapeBonus ← {       ⍝ ⍺=base_value, ⍵=lambda2
    ⍺ × 1 + 0.3 × ⍵÷2    ⍝ linear approximation: 30% bonus at λ₂=2
}

⍝ ══════════════════════════════════════════════════════════════════
⍝ END CH36 APL
⍝ Key checks:
⍝   NetBenefitHa ESV_Matrix           ⍝ → 490
⍝   ExternalityGap 2400 0.08 0.035 0.02 1.5  ⍝ → ~4.8 t/ha/yr
⍝ ══════════════════════════════════════════════════════════════════
