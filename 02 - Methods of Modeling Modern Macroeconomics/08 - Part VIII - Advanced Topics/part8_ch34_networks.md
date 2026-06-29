# Chapter 34: Network Models for Financial Contagion and Systemic Risk

*Clearing Vectors, Cascade Algorithms, and Centrality*

> *"In a network, the failure of one institution is not the end of the story — it is the beginning."*

**Cross-reference:** *Principles* Ch. 20 (financial markets, SDF and interconnection); Ch. 34 (financial crises, systemic risk, Too-Big-To-Fail); Ch. 40 (Great Recession: Lehman failure and cascade) **[P:Ch.20, P:Ch.34, P:Ch.40]**

---

## 34.1 The Network Approach to Systemic Risk

The 2008 financial crisis demonstrated that standard macroeconomic models — which model the financial sector as a single representative bank — could not anticipate or explain the cascade of failures that followed Lehman Brothers' bankruptcy. The reason: financial networks. Banks are linked through bilateral exposures (interbank loans, derivatives, repo agreements); the failure of one institution reduces the assets of its creditors, potentially triggering secondary failures in a cascade.

Network models formalize this mechanism. They represent the financial system as a **graph** in which nodes are banks and edges are bilateral exposures, then characterize when shocks propagate through the network and when they are contained.

---

## 34.2 Graph Theory Basics

**Definition 34.1 (Directed Weighted Graph).** A **financial network** is a directed weighted graph $G = (N, E, W)$ where: $N = \{1, \ldots, n\}$ is the set of banks; $E \subseteq N\times N$ is the set of (directed) interbank exposures; and $W_{ij} \geq 0$ is the face value of bank $j$'s obligation to bank $i$.

**Definition 34.2 (Adjacency Matrix).** The **adjacency matrix** $A$ has $A_{ij} = 1$ if bank $j$ has an obligation to bank $i$ (edge from $j$ to $i$); $A_{ij} = 0$ otherwise.

**Definition 34.3 (Degree Distribution).** The **in-degree** of bank $i$: $d_i^{in} = \sum_j A_{ji}$ (how many banks owe money to $i$). The **out-degree**: $d_i^{out} = \sum_j A_{ij}$ (how many banks $i$ owes money to).

**Definition 34.4 (Reachability and Contagion Paths).** Bank $j$ is **reachable** from bank $i$ in $k$ steps if $(A^k)_{ij} > 0$. In APL: `A ⍣ k +.× I[i;]` gives the set of banks reachable from $i$ in exactly $k$ steps.

---

## 34.3 The Eisenberg–Noe Model

The Eisenberg–Noe (2001) model is the foundational framework for analyzing interbank clearing under default.

**Setup:** Each bank $i$ has:
- External assets $c_i \geq 0$ (loans to non-bank borrowers, securities).
- Interbank liabilities $\bar{p}_i = \sum_j W_{ij}$ (total face-value obligations to other banks).
- A relative liability matrix $\Pi_{ij} = W_{ij}/\bar{p}_i$ (fraction of bank $i$'s total obligations owed to bank $j$, with $\sum_j \Pi_{ij} = 1$ for banks with positive liabilities).

**Definition 34.5 (Clearing Payment Vector).** A **clearing payment vector** $\mathbf{p}^* = (p_1^*, \ldots, p_n^*)$ satisfies:

$$p_i^* = \min\!\left\{\bar{p}_i,\; c_i + \sum_j\Pi_{ji}p_j^*\right\}, \quad p_i^* \geq 0.$$

This says: bank $i$ pays the minimum of its face-value obligation $\bar{p}_i$ and its available assets (external assets $c_i$ plus interbank receipts $\sum_j\Pi_{ji}p_j^*$). If assets exceed obligations, the bank pays in full and is solvent. If not, it defaults and pays everything it can.

**Definition 34.6 (Eisenberg–Noe Operator).** Define $\Phi: [0,\bar{\mathbf{p}}] \to [0,\bar{\mathbf{p}}]$ by:

$$\Phi_i(\mathbf{p}) = \min\!\left\{\bar{p}_i,\; c_i + \sum_j\Pi_{ji}p_j\right\}.$$

The clearing vector is a fixed point: $\mathbf{p}^* = \Phi(\mathbf{p}^*)$.

**Theorem 34.1 (Eisenberg–Noe Existence and Uniqueness).** Under mild regularity conditions (the network is connected or all banks have sufficient external assets), the Eisenberg–Noe clearing payment vector $\mathbf{p}^*$ exists and is unique. Moreover, $\mathbf{p}^*$ is the **greatest fixed point** of $\Phi$ — the one that maximizes total payments.

*Proof.* The operator $\Phi$ maps the complete lattice $[0,\bar{\mathbf{p}}]$ to itself and is **monotone**: $\mathbf{p} \leq \mathbf{q} \Rightarrow \Phi(\mathbf{p}) \leq \Phi(\mathbf{q})$ (because $\Pi$ has non-negative entries). By Tarski's fixed-point theorem, the greatest fixed point exists. Uniqueness follows from the strict contractiveness of $\Phi$ when banks have sufficient external assets: $\|\Phi(\mathbf{p}) - \Phi(\mathbf{q})\| \leq \|\Pi\|\|\mathbf{p}-\mathbf{q}\| < \|\mathbf{p}-\mathbf{q}\|$ (since $\|\Pi\| \leq 1$). $\square$

**Algorithmic solution:** The greatest fixed point is found by iterating $\mathbf{p}^{k+1} = \Phi(\mathbf{p}^k)$ starting from $\mathbf{p}^0 = \bar{\mathbf{p}}$ (all banks pay in full). The iteration converges because $\Phi$ is a contraction.

In APL: the EN operator and iteration are concise:

```apl
⍝ APL — Eisenberg-Noe clearing vector
⎕IO←0 ⋄ ⎕ML←1

⍝ EN operator: Φ(p) = min(p_bar, c + Π' p)
⍝ Pi_T = transpose of relative liability matrix Π (n×n)
EN_op ← {p_bar c Pi_T p ← ⍵
    p_bar ⌊ c + Pi_T +.× p}

⍝ Find clearing vector: iterate from p = p_bar
⍝ Converge when two consecutive iterates agree to tolerance
clear_vec ← {p_bar c Pi ← ⍵
    step ← {EN_op p_bar c (⍉Pi) ⍵}
    step ⍣ (1e¯8∘>⌈/|⊢-step) ⊢ p_bar}

⍝ Example: 4-bank network
n ← 4
⍝ Face-value liabilities (p_bar) and external assets (c)
p_bar ← 100 80 60 90
c     ← 50 30 40 70

⍝ Relative liability matrix (rows sum to 1)
Pi ← 4 4 ⍴ 0 0.4 0.3 0.3
           0.5 0 0.3 0.2
           0.3 0.4 0 0.3
           0.2 0.3 0.5 0

⍝ Add shock: bank 1 loses 40 of external assets
c_shocked ← c - 0 40 0 0

p_baseline ← clear_vec p_bar c     Pi
p_shocked  ← clear_vec p_bar c_shocked Pi

p_baseline    ⍝ all pay in full
p_shocked     ⍝ bank 2 defaults; may cascade
```

---

## 34.4 Contagion Cascades and Default Propagation

When bank $j$ defaults and pays only $p_j^* < \bar{p}_j$, its creditors receive less than expected. If the shortfall is large enough, creditors may also default — a **contagion cascade**.

**Algorithm 34.1 (Default Cascade).**

Initialize: all banks alive; shocks reduce $c_i$ for some banks.

For each round $k = 0, 1, 2, \ldots$:
1. Compute each bank's assets: $A_i^{(k)} = c_i + \sum_j\Pi_{ji}p_j^{(k)}$.
2. Compute new payments: $p_i^{(k+1)} = \min(\bar{p}_i, A_i^{(k)})$.
3. If $\mathbf{p}^{(k+1)} = \mathbf{p}^{(k)}$: stop (fixed point reached).

Note: this is exactly the Eisenberg–Noe iteration. The algorithm converges to $\mathbf{p}^*$ in finite steps for finite networks.

**Percolation threshold for cascade vulnerability:** For large random networks (Erdős–Rényi), the cascade is large (infects a positive fraction of the network) iff the network exceeds a **percolation threshold**.

**Theorem 34.2 (Percolation Threshold for Erdős–Rényi Network).** For an Erdős–Rényi random graph $G(n, p)$ (each edge present independently with probability $p$), the giant connected component (GCC) appears when:

$$p > p_c = \frac{1}{n}.$$

For a financial contagion cascade, the cascade is large iff the average degree $\bar{d} = pn > 1$. Below the threshold: shocks are contained; above: they propagate to $O(n)$ banks.

*Proof sketch.* At the threshold, each infected bank has on average one neighbor. The cascade size follows a branching process with offspring distribution $\text{Poisson}(\bar{d})$. Survival probability of the branching process is positive iff $\bar{d} > 1$ (Galton–Watson criterion). $\square$

---

## 34.5 Network Centrality and Systemic Importance

Some banks are more systemically important than others — their failure would propagate more widely. Network centrality measures quantify this.

**Definition 34.7 (Eigenvector Centrality).** The **eigenvector centrality** vector $\mathbf{v}$ satisfies $A\mathbf{v} = \lambda_1\mathbf{v}$, where $\lambda_1$ is the largest eigenvalue of the adjacency matrix. Bank $i$'s centrality $v_i$ is proportional to the sum of centralities of its neighbors: a bank is central if it is connected to other central banks.

In APL: `{(A+.×⍵)÷+/A+.×⍵}⍣≡(N⍴1÷N)` — power iteration to find the leading eigenvector.

**Definition 34.8 (DebtRank).** The **DebtRank** of bank $i$ (Battiston et al., 2012):

$$DR_i = \sum_{j\neq i}h_j(i)\cdot\frac{e_j}{\sum_k e_k},$$

where $h_j(i)$ is the fractional impact on bank $j$ from bank $i$'s default ($h_j = \min(1, \Pi_{ij}\cdot\text{leverage}_i)$) and $e_j$ is bank $j$'s equity. DebtRank measures the fraction of total financial system equity that would be lost following bank $i$'s default.

---

## 34.6 Worked Example: 2008 Interbank Network

*Cross-reference: Principles Ch. 40 (Great Recession)* **[P:Ch.40]**

We simulate a stylized 20-bank network calibrated to the rough structure of the 2007 U.S. interbank market.

```python
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(42)
n = 20  # banks

# Generate a random financial network (Erdős–Rényi with p=0.3)
prob_link = 0.3
A_raw = (np.random.rand(n, n) < prob_link).astype(float)
np.fill_diagonal(A_raw, 0)

# Face-value liabilities (heterogeneous size)
p_bar = np.random.exponential(100, n) + 50
c_ext = np.random.exponential(60, n) + 30  # external assets

# Normalize to get relative liability matrix
row_sums = A_raw.sum(axis=1, keepdims=True)
row_sums[row_sums == 0] = 1  # avoid division by zero
W_raw = A_raw * p_bar[:, None]  # w_ij = p_bar_i * A_ij / degree
# Relative liability matrix (fraction of obligations to each creditor)
Pi = W_raw / np.maximum(row_sums * p_bar[:, None], 1)
# Normalize Pi so columns sum to <= 1
col_sum = Pi.sum(axis=0, keepdims=True); col_sum[col_sum == 0] = 1
Pi = Pi / col_sum  # re-normalize

def eisenberg_noe(p_bar, c, Pi, max_iter=1000, tol=1e-8):
    """Find Eisenberg-Noe clearing payment vector."""
    p = p_bar.copy()
    for _ in range(max_iter):
        p_new = np.minimum(p_bar, c + Pi.T @ p)
        if np.max(np.abs(p_new - p)) < tol:
            return p_new
        p = p_new
    return p

# Baseline (no shock)
p_star_baseline = eisenberg_noe(p_bar, c_ext, Pi)
n_default_baseline = np.sum(p_star_baseline < p_bar * 0.99)
print(f"Baseline defaults: {n_default_baseline}")

# Shock: biggest bank loses 70% of external assets
biggest = np.argmax(p_bar)
c_shocked = c_ext.copy()
c_shocked[biggest] *= 0.3  # 70% loss

p_star_shocked = eisenberg_noe(p_bar, c_shocked, Pi)
n_default_shocked = np.sum(p_star_shocked < p_bar * 0.99)
defaults = np.where(p_star_shocked < p_bar * 0.99)[0]
total_loss = np.sum(p_bar[defaults] - p_star_shocked[defaults])
print(f"Post-shock defaults: {n_default_shocked} banks")
print(f"Total loss: {total_loss:.1f}")

# Eigenvector centrality
def eigenvector_centrality(A, max_iter=100, tol=1e-10):
    n = A.shape[0]; v = np.ones(n)/n
    for _ in range(max_iter):
        v_new = A @ v
        norm = v_new.sum()
        if norm > 0: v_new /= norm
        if np.max(np.abs(v_new - v)) < tol: return v_new
        v = v_new
    return v

ec = eigenvector_centrality(A_raw)
print(f"\nMost central banks (top 3): {np.argsort(ec)[-3:][::-1]}")
print(f"Biggest bank centrality rank: {np.sum(ec > ec[biggest])}")

# DebtRank
equity = np.maximum(c_ext - p_bar * 0.3, 1)  # simplified equity
E_total = equity.sum()
debt_rank = np.zeros(n)
for i in range(n):
    leverage_i = p_bar[i] / np.maximum(equity[i], 1)
    h = np.minimum(1, Pi[i, :] * leverage_i)
    debt_rank[i] = np.sum(h * equity) / E_total

print(f"\nDebtRank — top 3 systemically important banks: {np.argsort(debt_rank)[-3:][::-1]}")
print(f"Max DebtRank: {debt_rank.max():.3f}")
```

```julia
using LinearAlgebra, SparseArrays

function eisenberg_noe(p_bar, c, Pi; tol=1e-8, maxiter=1000)
    p = copy(p_bar)
    for _ in 1:maxiter
        p_new = min.(p_bar, c .+ Pi' * p)
        maximum(abs.(p_new .- p)) < tol && return p_new
        p .= p_new
    end
    return p
end

# Eigenvector centrality via power iteration
function eigenvec_cent(A; maxiter=100)
    n = size(A,1); v = ones(n)/n
    for _ in 1:maxiter
        v_new = A * v; s = sum(v_new)
        s > 0 && (v_new ./= s)
        maximum(abs.(v_new .- v)) < 1e-10 && return v_new
        v .= v_new
    end; v
end

n=20; Random.seed!(42)
A = (rand(n,n) .< 0.3) .& .!I(n); A=Float64.(A)
p_bar = rand(Exponential(100),n).+50; c=rand(Exponential(60),n).+30
Pi = A .* p_bar; Pi ./= max.(sum(Pi,dims=1),1)

p_base = eisenberg_noe(p_bar,c,Pi)
println("Defaults (baseline): $(sum(p_base .< 0.99*p_bar))")

c_shock = copy(c); c_shock[argmax(p_bar)] *= 0.3
p_shock = eisenberg_noe(p_bar,c_shock,Pi)
println("Defaults (shocked): $(sum(p_shock .< 0.99*p_bar))")

ec = eigenvec_cent(A)
println("Top 3 central banks: $(sortperm(ec, rev=true)[1:3])")
```

---

## 34.7 Programming Exercises

### Exercise 34.1 (APL — Cascade Propagation)

Implement the default cascade in APL using Boolean matrix operations. (a) Represent solvency as a Boolean vector `solvent ← p_star ≥ 0.99 × p_bar`. (b) One cascade step: `cascade_step ← {solvent assets p ← ⍵ ⋄ assets_new ← c + Pi_T +.× p ⋄ p_new ← p_bar ⌊ assets_new ⋄ p_new solvent_new}`. (c) Show the cascade stops in $O(\log n)$ rounds using the Boolean reachability operator `A⍣k`.

### Exercise 34.2 (Python — Scale-Free Network)

A **scale-free** (Barabási–Albert) network has degree distribution $P(d) \propto d^{-3}$ — a few highly connected nodes (hubs). Generate a BA network: (a) start with 3 connected nodes; (b) add one node at a time, connecting to existing nodes with probability proportional to their degree; (c) when the network has 50 nodes, compare the default cascade size to an Erdős–Rényi network with the same average degree. Which architecture is more fragile? Which more robust? Connect to the Acemoglu et al. (2015) robustness-fragility result.

### Exercise 34.3 (Julia — Systemic Risk Measures)

```julia
# Compare three systemic risk measures on the same network
function systemic_risk_comparison(p_bar, c, Pi, equity, A)
    n = length(p_bar)
    
    # 1. SRISK: expected shortfall under stress scenario
    c_stress = c .* 0.5
    p_stress = eisenberg_noe(p_bar, c_stress, Pi)
    srisk = sum(max.(0, p_bar .- p_stress .- equity))
    
    # 2. Eigenvector centrality
    ec = eigenvec_cent(A)
    
    # 3. DebtRank
    E_tot = sum(equity)
    debt_rank = [sum(min.(1, Pi[i,:] .* (p_bar[i]/max(equity[i],1))) .* equity)/E_tot for i in 1:n]
    
    return srisk, ec, debt_rank
end

println("Systemic risk measures comparison complete")
```

### Exercise 34.4 — Eisenberg–Noe with Recovery Rates ($\star$)

In the basic EN model, defaulting banks pay everything available ($recovery = 100\%$). In practice, bankruptcy imposes costs: creditors receive only fraction $\alpha < 1$ of assets. Modify the EN operator: $p_i^* = \min(\bar{p}_i, c_i + \alpha\sum_j\Pi_{ji}p_j^*)$ for defaulting banks ($p_i^* < \bar{p}_i$) and $p_i^* = \bar{p}_i$ for solvent banks. (a) Show the modified operator is still monotone and the fixed point exists. (b) Compute how the default cascade size changes with $\alpha \in \{0.3, 0.5, 0.7, 1.0\}$ for the 20-bank example. (c) Interpret: why does lower recovery amplify contagion?

---

## 34.8 Chapter Summary

**Key results:**

- A financial network $G = (N, E, W)$ has the **adjacency matrix** $A_{ij} = 1$ if $j$ owes to $i$; **reachability** in $k$ steps is $(A^k)_{ij} > 0$; in APL: `A⍣n`.
- The **Eisenberg–Noe clearing vector** $\mathbf{p}^* = \Phi(\mathbf{p}^*)$ exists and is unique (Theorem 34.1 via Tarski's fixed-point theorem); computed by iterating $\mathbf{p}^{k+1} = \min(\bar{\mathbf{p}}, \mathbf{c}+\Pi'\mathbf{p}^k)$ from $\mathbf{p}^0 = \bar{\mathbf{p}}$; in APL: `EN_op ⍣ converged ⊢ p_bar`.
- The **percolation threshold** for cascade propagation in Erdős–Rényi networks: $\bar{d} = np > 1$ (Theorem 34.2); below threshold shocks are contained; above, $O(n)$ banks fail.
- **Eigenvector centrality** $A\mathbf{v} = \lambda_1\mathbf{v}$ identifies systemically important banks; APL power iteration: `{(A+.×⍵)÷+/A+.×⍵}⍣≡(N⍴1÷N)`.
- **DebtRank** quantifies the fraction of system equity at risk from each bank's failure; captures both direct and indirect contagion channels.

*Next: Chapter 35 — Integrated Assessment Models: Climate-Economy Modeling*
