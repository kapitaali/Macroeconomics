⍝ ══════════════════════════════════════════════════════════════════
⍝ Chapter 25: Mutual Credit — APL
⍝ Balance vector operations and clearing efficiency
⍝ ══════════════════════════════════════════════════════════════════

⍝ ── 1. SYSTEM BALANCE CONSTRAINT ───────────────────────────────────
⍝ In any mutual credit system: sum of all balances = 0
BalanceOK ← { 1e¯9 > |+/⍵ }    ⍝ 1 if zero-sum (within tolerance)

⍝ ── 2. BILATERAL NETTING ────────────────────────────────────────────
⍝ Given gross obligation matrix G (n×n, G[i;j] = amount i owes j),
⍝ bilateral netting replaces each pair with net obligation:
⍝   G_net[i;j] = max(0, G[i;j] - G[j;i])

BilateralNet ← {
    G    ← ⍵
    Gt   ← ⍉G              ⍝ transpose = obligations in reverse direction
    0⌈G-Gt                 ⍝ element-wise max with 0
}

⍝ Clearing efficiency: gross cleared / gross total
⍝   Bilateral clears: sum(G) - sum(BilateralNet G)
⍝   Remaining after bilateral: sum(BilateralNet G)

BilateralEfficiency ← {
    G   ← ⍵
    GN  ← BilateralNet G
    gross ← +/+/G
    remain← +/+/GN
    (gross-remain)÷gross
}

⍝ ── 3. NET BALANCE VECTOR FROM OBLIGATION MATRIX ───────────────────
⍝ Balance_i = (sum of column i) - (sum of row i)
⍝ = credits received - debts owed

NetBalances ← {
    G ← ⍵
    (+/G) - +⌿G    ⍝ col sums - row sums
}

⍝ ── 4. LYAPUNOV IMBALANCE ENERGY ───────────────────────────────────
⍝ V(b) = (1/2) ||b||^2  — system stress measure (Section 25.5)

ImbalanceEnergy ← { 0.5 × +/⍵*2 }    ⍝ ⍵ = balance vector b

⍝ ── 5. SMALL EXAMPLE: 3-MEMBER NETWORK ─────────────────────────────
⍝ Baker (B), Carpenter (C), Farmer (F)
⍝ C owes B 100, B owes F 70, F owes C 50

G3 ← 3 3 ⍴ 0 0 70  100 0 0  0 50 0   ⍝ G[i;j] = i owes j

balances3 ← NetBalances G3    ⍝ → 30 ¯50 20  (B=+30, C=¯50, F=+20)
BalanceOK balances3            ⍝ → 1 (sums to 0)
ImbalanceEnergy balances3      ⍝ → 0.5×(900+2500+400) = 1900

GN3 ← BilateralNet G3         ⍝ after bilateral netting
BE3 ← BilateralEfficiency G3  ⍝ efficiency fraction

⍝ After full multilateral clearing (clearing the circular C→B→F→C):
⍝ flow = min(100, 70, 50) = 50 through the cycle
⍝ Residual gross obligations: (100-50)=50 (C→B), (70-50)=20 (B→F), 0 (F→C)
⍝ Cleared: 150 of original 220 = 68%

⍝ ══════════════════════════════════════════════════════════════════
⍝ END CH25 APL
⍝ Key checks:
⍝   NetBalances G3     ⍝ → 30 ¯50 20
⍝   BalanceOK balances3 ⍝ → 1
⍝   +/+/G3            ⍝ → 220 (total gross obligations)
⍝ ══════════════════════════════════════════════════════════════════
