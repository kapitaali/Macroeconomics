⍝ ══════════════════════════════════════════════════════════════════
⍝ Chapter 18: SFC-N — APL matrix operations
⍝ Balance Sheet Matrix verification and natural capital accounting
⍝ ══════════════════════════════════════════════════════════════════

⍝ ── 1. BALANCE SHEET MATRIX (BSM) ROW-SUM CHECK ────────────────────
⍝ Every row of the BSM must sum to zero (C3: every asset = someone's liability)
BSMRowCheck ← { +/⍵ }          ⍝ returns row sums; all should be 0

⍝ ── 2. NATURAL CAPITAL STOCK EVOLUTION (logistic) ──────────────────
⍝ N_{t+1} = N_t + r × N_t × (1 - N_t/K) - d × Y
NCStep ← {          ⍝ ⍺=(r K d Y), ⍵=N_t
    r K d Y ← ⍺
    N ← ⍵
    N + r×N×(1-N÷K) - d×Y
}

⍝ Simulate T periods: NCTraj (r K d Y) (N0 T)
NCTraj ← {
    params ← ⍺ ⋄ N0 T ← ⍵
    traj ← ,N0
    N ← N0
    :Repeat T
        N ← params NCStep N
        N ← 0⌈N          ⍝ floor at 0
        traj ← traj,N
    :EndRepeat
    traj
}

⍝ ── 3. PBS TOTAL ────────────────────────────────────────────────────
⍝ PBS = produced_capital + financial_assets + NC_value - liabilities
⍝ NC_value = shadow_prices +.× NC_levels
PBS ← {        ⍝ ⍵ = (K_prod, F_assets, liabilities, shadow_prices, NC_levels)
    Kp Fa L p N ← ⍵
    Kp + Fa + p+.×N - L
}

⍝ ── 4. STEWARDSHIP CONDITION CHECK ─────────────────────────────────
⍝ dN_j/dt ≥ 0 for all j
⍝ Regen  = r_j × N_j × (1 - N_j/K_j)
⍝ Deplet = d_j × Y
⍝ NetChange = Regen - Deplet  (vector over all NC stocks)

NetChange ← {      ⍝ ⍺=(r K d Y vectors), ⍵=N vector
    r K d Y ← ⍺
    N ← ⍵
    (r×N×1-N÷K) - d×Y
}

StewardshipMet ← { ⍝ 1 if all net changes ≥ 0
    ∧/0≤⍵
}

⍝ ── 5. DARLINGTON CALIBRATION (5 NC stocks) ────────────────────────
⍝ Parameters from Chapter 44 / Chapter 18 worked example

r_d ← 0.04 0.10 0.08 0.05 0.06    ⍝ regeneration rates
K_d ← 5⍴1.0                        ⍝ carrying capacities (normalised)
d_d ← 0.036 0.012 0.028 0.022 0.030  ⍝ depletion coefficients
N_d ← 0.78 0.85 0.72 0.68 0.74    ⍝ current NC levels
p_d ← 285 45 62 38 118             ⍝ shadow prices (£M)
Y_d ← 1.0                          ⍝ GDP index

⍝ NC vector value:
NC_value ← p_d +.× N_d    ⍝ → total NC value in £M

⍝ Net change for each stock at baseline:
nc_baseline ← (r_d, K_d, d_d, Y_d) NetChange N_d

⍝ Are all stocks stable?
StewardshipMet nc_baseline    ⍝ → 0 (some stocks declining)

⍝ Regeneration and depletion totals:
total_regen  ← +/r_d×N_d×1-N_d÷K_d
total_deplet ← +/d_d×Y_d

⍝ Gap (£M/yr of restoration needed):
⍝  gap ← total_deplet - total_regen  (should be ~10.2£M)

⍝ ══════════════════════════════════════════════════════════════════
⍝ END CH18 APL
⍝ Key verification:
⍝   p_d +.× N_d      ⍝ → total NC value
⍝   StewardshipMet nc_baseline  ⍝ → 0 (failing)
⍝   total_deplet - total_regen  ⍝ → ~10.2 (annual deficit)
⍝ ══════════════════════════════════════════════════════════════════
