⍝ ══════════════════════════════════════════════════════════════════
⍝ Chapter 29: Unified Model — APL
⍝ Coupling terms, IPI computation, and the Shadow Price Identity
⍝ ══════════════════════════════════════════════════════════════════

⍝ ── 1. IPI (Intertemporal Provisioning Index) ───────────────────────
⍝ IPI = Σ_t β^t × Σ_i U_i(C_{i,t})
⍝ Discrete approximation: U_i(C) = ln(C) × (1 - G)  (inequality-adjusted log utility)

IPI ← {          ⍝ ⍺=discount factor β, ⍵=(C_vector, G_vector) over T+1 periods
    beta ← ⍺
    C G  ← ⍵
    T    ← ≢C
    disc ← beta * ⍳T    ⍝ β^0, β^1, …, β^(T-1)
    U    ← (⍟C) × 1-G   ⍝ inequality-adjusted log utility per period
    disc +.× U           ⍝ discounted sum
}

⍝ ── 2. COOPERATIVE SURPLUS EFFECT ON OUTPUT ─────────────────────────
⍝ CRE output premium: Y_CRE = Y_CE × (1 + σ_v × pass-through)
CoopOutputPremium ← {    ⍝ ⍺=Y_CE, ⍵=(sigma_v, passthrough)
    sigma pt ← ⍵
    ⍺ × 1 + sigma × pt
}

⍝ ── 3. SHADOW PRICE IDENTITY (Proposition 29.3) ─────────────────────
⍝ p^N_j(t) = λ_j(t) / λ_x(t)
⍝ Natural capital shadow price = co-state ratio
⍝ Given co-state vectors λ_N and λ_x at each period:

ShadowPrice ← { ⍺÷⍵ }    ⍝ ⍺=lambda_N, ⍵=lambda_x (element-wise)

⍝ ── 4. PLANETARY BOUNDARIES CONSTRAINT CHECK ───────────────────────
⍝ b_j(x) ≤ b_bar_j for all j  (Equation 29.5)
⍝ ⍺=current boundary values, ⍵=safe limits

PBcheck ← { ∧/⍺ ≤ ⍵ }    ⍝ 1 if all boundaries within safe limits

⍝ ── 5. KAKUTANI FIXED POINT EXISTENCE (numerical check) ─────────────
⍝ The CRE is guaranteed by Kakutani if the best-response correspondence
⍝ is non-empty, convex-valued, and upper-hemicontinuous.
⍝ We verify a necessary condition: feasible set is compact (bounded NC stocks).

FeasibleSetCompact ← {     ⍝ ⍵ = (N_vector, N_crit_vector, M, G)
    N Ncrit M G ← ⍵
    nc_ok  ← ∧/N ≥ Ncrit   ⍝ all NC above critical threshold
    money_ok ← M ≥ 0        ⍝ non-negative money stock
    gini_ok  ← (0 ≤ G) ∧ (G ≤ 1)
    nc_ok ∧ money_ok ∧ gini_ok
}

⍝ ── 6. COOPERATIVE STEWARDSHIP THEOREM — IPI COMPARISON ────────────
⍝ Compare discounted utility sums for CRE vs CE trajectories

IPIgap ← {     ⍝ ⍵=(ipi_cre, ipi_ce) — returns percentage gap
    ipi_cre ipi_ce ← ⍵
    100 × (ipi_cre - ipi_ce) ÷ ipi_ce
}

⍝ Danish calibration result:
⍝   ipi_cre ← 387.4   ipi_ce ← 281.5    (illustrative)
⍝   IPIgap ipi_cre ipi_ce    ⍝ → ~37.6% IPI advantage for CRE

⍝ ══════════════════════════════════════════════════════════════════
⍝ END CH29 APL
⍝ ══════════════════════════════════════════════════════════════════
