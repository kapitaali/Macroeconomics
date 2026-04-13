⍝ ══════════════════════════════════════════════════════════════════
⍝ Chapter 12: Network Structure — APL
⍝ Graph Laplacian, Fiedler value, and spectral operations
⍝ ══════════════════════════════════════════════════════════════════

⍝ ── 1. GRAPH LAPLACIAN ─────────────────────────────────────────────
⍝ L = D - A  where D is the degree matrix (diagonal)
Laplacian ← { (⍤1⊢(+/⍵)) × =⍨⍳≢⍵ ) - ⍵ }
⍝ Cleaner version:
Laplacian ← {
    A ← ⍵
    D ← (+/A) × =/∘⍳¨⍨⍳≢A   ⍝ diagonal degree matrix: faster via outer product
    D - A
}

⍝ ── 2. ALGEBRAIC CONNECTIVITY (Fiedler value) ──────────────────────
⍝ Second-smallest eigenvalue of L.
⍝ In APL we rely on ⌹ (matrix divide) for eigenvalues — or call via ⎕NA.
⍝ Here we use the power-iteration approach for the 2nd eigenvalue:

⍝ Largest eigenvalue via power iteration:
PowerIter ← {           ⍝ ⍺=matrix, ⍵=n_iterations
    A ← ⍺ ⋄ n ← ≢A
    v ← n⍴1÷n*0.5       ⍝ initial unit vector
    :Repeat ⍵
        v ← A+.×v
        v ← v÷((+/v*2)*0.5)    ⍝ normalise
    :EndRepeat
    (+/v×A+.×v)÷+/v*2   ⍝ Rayleigh quotient
}

⍝ Fiedler value: λ₂ = λ_min(L) on the subspace orthogonal to 1.
⍝ Deflate L by removing the contribution of the zero eigenvector (1/√n):
FiedlerValue ← {        ⍝ ⍵ = adjacency matrix A (n×n)
    L ← Laplacian ⍵
    n ← ≢L
    ⍝ Project onto orthogonal complement of (1 1 … 1):
    P ← (=/∘⍳¨⍨⍳n) - ÷n     ⍝ projection matrix I - (1/n)11^T
    LP← P+.×L+.×P             ⍝ deflated Laplacian
    LP PowerIter 200           ⍝ largest eigenvalue of LP = λ₂ of L
}

⍝ ── 3. SHOCK TRANSMISSION SPECTRAL RADIUS ──────────────────────────
⍝ ρ(A) via power iteration on the adjacency matrix

SpectralRadius ← { ⍵ PowerIter 200 }

⍝ Stability condition: ρ(A) < γ (shock absorption rate)
ShockStable ← {    ⍝ ⍺=adjacency matrix, ⍵=gamma
    (SpectralRadius ⍺) < ⍵
}

⍝ ── 4. DEGREE DISTRIBUTION STATISTICS ──────────────────────────────
⍝ Degree vector, mean, variance, coefficient of variation

DegreeVec ← { +/⍵ }           ⍝ row sums of adjacency matrix = degrees
DegMean   ← { (+/DegreeVec ⍵) ÷ ≢⍵ }
DegVar    ← { (+/(d←DegreeVec ⍵-DegMean ⍵)*2) ÷ ≢⍵ }
DegCV     ← { (DegVar ⍵)*0.5 ÷ DegMean ⍵ }  ⍝ coefficient of variation

⍝ ── 5. SMALL EXAMPLE: 4-NODE RING ──────────────────────────────────
⍝ A_ring = [[0,1,0,1],[1,0,1,0],[0,1,0,1],[1,0,1,0]]
A_ring ← 4 4 ⍴ 0 1 0 1  1 0 1 0  0 1 0 1  1 0 1 0
⍝ FiedlerValue A_ring   ⍝ → 2 (ring: λ₂ = 2(1-cos(2π/4)) = 2)
⍝ DegCV A_ring          ⍝ → 0 (regular graph, all degrees equal)

⍝ ══════════════════════════════════════════════════════════════════
⍝ END CH12 APL — key checks:
⍝   DegMean A_ring      ⍝ → 2.0
⍝   DegCV A_ring        ⍝ → 0.0 (regular = zero inequality)
⍝ ══════════════════════════════════════════════════════════════════
