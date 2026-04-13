⍝ ══════════════════════════════════════════════════════════════════
⍝ Chapter 10: Agent-Based Models of Cooperation vs. Competition
⍝ APL — matrix-level operations on the ABM state
⍝
⍝ The Python file handles the full Mesa ABM.
⍝ APL handles the payoff matrix algebra and replicator dynamics —
⍝ operations that are one-liners in APL and ten lines in Python.
⍝ ══════════════════════════════════════════════════════════════════


⍝ ── 1. PAYOFF MATRIX ───────────────────────────────────────────────
⍝ Rows = own strategy (0=cooperate, 1=defect)
⍝ Cols = opponent strategy
⍝ PD payoffs: R=3, S=0, T=5, P=1

PayoffMatrix ← 2 2 ⍴ 3 0 5 1   ⍝ ┌ CC  CD ┐
                                 ⍝ └ DC  DD ┘

⍝ Expected payoff for a cooperator in a population with fraction p cooperating:
⍝   E[C] = p×R + (1-p)×S = p×3 + (1-p)×0 = 3p
ExpCoop ← { 3×⍵ }               ⍝ ⍵ = p (cooperation fraction)

⍝ Expected payoff for a defector:
⍝   E[D] = p×T + (1-p)×P = 5p + 1-p = 4p+1
ExpDef  ← { 1+4×⍵ }

⍝ Cooperation pays when E[C] > E[D]:  3p > 4p+1  =>  never in one-shot PD
⍝ (Confirms the theoretical result — cooperation requires repetition.)
CoopBetter ← { (ExpCoop ⍵) > ExpDef ⍵ }   ⍝ returns 0 for all p ∈ [0,1]


⍝ ── 2. REPLICATOR DYNAMICS  (Section 10.2) ─────────────────────────
⍝
⍝ ṗ = p(1-p)[E[C](p) - E[D](p)]
⍝
⍝ Discrete-time approximation with step h:
⍝   p_{t+1} = p_t + h × p_t(1-p_t)(E[C]-E[D])

ReplicatorStep ← {           ⍝ ⍺=h (step size), ⍵=p (current fraction)
    p   ← ⍵
    ec  ← ExpCoop p
    ed  ← ExpDef  p
    p + ⍺ × p × (1-p) × ec-ed
}

⍝ Simulate T steps of replicator dynamics
⍝ Usage: h ReplicatorSim (p0, T)
ReplicatorSim ← {                  ⍝ ⍺=h, ⍵=(p0 T)
    h ← ⍺
    p ← ⊃⍵ ⋄ T ← ⊃1↓⍵
    traj ← ,p
    :Repeat T
        p ← h ReplicatorStep p
        traj ← traj,p
    :EndRepeat
    traj
}

⍝ Example — one-shot PD always collapses to p=0 (defection ESS):
⍝   0.01 ReplicatorSim 0.8 100   ⍝ p decreases monotonically to 0


⍝ ── 3. NETWORK PAYOFF AGGREGATION ──────────────────────────────────
⍝
⍝ Given a strategy vector S (1=cooperate, 0=defect) and
⍝ adjacency matrix A, compute each agent's payoff.
⍝
⍝ Payoff_i = Σ_j A_ij × PayoffMatrix[S_i; S_j]
⍝
⍝ Vectorised: no loops needed in APL.

NetworkPayoffs ← {          ⍝ ⍺=adjacency matrix A (n×n), ⍵=strategy vector S (n,)
    A  ← ⍺
    S  ← ⍵                  ⍝ 1=cooperate, 0=defect (integer)
    n  ← ≢S

    ⍝ Build n×n payoff matrix where entry (i,j) = PayoffMatrix[S[i]; S[j]]
    ⍝ S_i repeated across columns, S_j repeated across rows
    Si ← S ∘.= 0 1          ⍝ n×2: one-hot encoding of own strategy
    Sj ← S ∘.= 0 1          ⍝ n×2: one-hot encoding of opponent strategy

    ⍝ Full payoff pair matrix: PM[i,j] = PayoffMatrix[S_i+1; S_j+1]
    ⍝ Exploit: S_i ∈ {0,1}, so index = S_i×2 + S_j + 1
    idx ← (2×S) ∘.+ S       ⍝ n×n matrix of flat indices into PayoffMatrix
    PM  ← (,PayoffMatrix)[1+idx]  ⍝ n×n payoff values

    ⍝ Each agent's payoff = row sum of (A ⊙ PM)
    +/A × PM                ⍝ n-vector of payoffs
}

⍝ Usage example (3 agents, fully connected):
⍝   A ← 3 3 ⍴ 0 1 1 1 0 1 1 1 0
⍝   S ← 1 1 0          ⍝ agents 1,2 cooperate; agent 3 defects
⍝   NetworkPayoffs A S  ⍝ → payoffs for each agent


⍝ ── 4. COOPERATION FRACTION DYNAMICS (population-level) ────────────
⍝
⍝ Given payoff vectors for cooperators and defectors, compute
⍝ the next-period cooperation fraction using the replicator rule.

PopulationUpdate ← {         ⍝ ⍵ = (p, avg_payoff_C, avg_payoff_D)
    p ← ⊃⍵
    fc← ⊃1↓⍵
    fd← ⊃2↓⍵
    fbar ← p×fc + (1-p)×fd   ⍝ population average payoff
    p×fc÷fbar                 ⍝ replicator: p' = p×f_C / f̄
}

⍝ Example — cooperative surplus from network reciprocity:
⍝ With network reciprocity, cooperators' avg payoff is boosted by
⍝ the clustering bonus (see Chapter 12 for Fiedler-value connection)
⍝   p  ← 0.5
⍝   fc ← 2.5   ⍝ boosted by network clustering
⍝   fd ← 1.8
⍝   PopulationUpdate p fc fd  ⍝ → > 0.5 (cooperation increases)


⍝ ── 5. SOBOL INDEX APPROXIMATION (simple, no SALib needed) ─────────
⍝
⍝ Estimates first-order Sobol index for parameter X_k
⍝ using the Saltelli estimator on pre-computed output vectors.
⍝
⍝ SobolS1 ← { ⍝ ⍺=f_A (n-vector), ⍵=f_AB_k (n-vector with k-th col swapped)
⍝     n ← ≢⍺
⍝     f0 ← (+/⍺)÷n          ⍝ grand mean
⍝     V  ← (+/⍺*2)÷n - f0*2 ⍝ total variance
⍝     ((+/⍵×⍺)÷n - f0*2) ÷ V
⍝ }

SobolS1 ← {
    n  ← ≢⍺
    f0 ← (+/⍺)÷n
    V  ← (+/⍺*2)÷n - f0*2
    ((+/⍵×⍺)÷n - f0*2) ÷ V
}


⍝ ══════════════════════════════════════════════════════════════════
⍝ END OF CHAPTER 10 APL
⍝
⍝ Key verifications in APL session:
⍝   CoopBetter 0.8       ⍝ → 0  (defection always dominates in one-shot PD)
⍝   CoopBetter 0.0       ⍝ → 0
⍝   PayoffMatrix         ⍝ → 2×2 matrix: R=3, S=0, T=5, P=1
⍝   0.01 ReplicatorSim 0.8 50   ⍝ → p decreasing from 0.8 toward 0
⍝
⍝ For the full ABM use ch10_python.py (Mesa implementation).
⍝ APL provides the analytical backbone: payoff algebra and
⍝ replicator dynamics that confirm the ABM's qualitative results.
⍝ ══════════════════════════════════════════════════════════════════
