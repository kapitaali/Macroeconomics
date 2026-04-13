⍝ ══════════════════════════════════════════════════════════════════
⍝ Chapter 42: Experimentation — APL
⍝ Multi-armed bandit analytics and policy diffusion
⍝ ══════════════════════════════════════════════════════════════════

⍝ ── 1. BAYESIAN POSTERIOR UPDATE (normal-normal) ────────────────────
⍝ Prior N(μ₀, σ₀²), Likelihood N(obs, σ_noise²)
⍝ Posterior: precision adds, mean is precision-weighted average
⍝ ⍺=(prior_mean prior_var noise_var), ⍵=observation

BayesUpdate ← {
    mu0 v0 vn ← ⍺
    obs ← ⍵
    prec0 ← ÷v0 ⋄ precn ← ÷vn
    post_prec  ← prec0+precn
    post_var   ← ÷post_prec
    post_mean  ← (prec0×mu0 + precn×obs) ÷ post_prec
    post_mean, post_var          ⍝ return (μ₁, σ₁²)
}

⍝ ── 2. REGRET CALCULATION ────────────────────────────────────────────
⍝ Cumulative regret = T × best_true_effect - Σ observed_welfare
CumRegret ← {           ⍝ ⍺=best_true_effect (scalar), ⍵=observed_welfare_vector
    (⍺ × ⍳≢⍵) - +\⍵
}

⍝ ── 3. OPTIMAL EXPERIMENT SCALE (Proposition 42.1) ──────────────────
⍝ n* = sqrt( 2σ²(1 + ρ×N_rep) / V'' )
⍝ ⍵ = (sigma2, rho, N_rep, V_double_prime)

OptExpScale ← {
    sigma2 rho Nrep Vpp ← ⍵
    (2×sigma2×1+rho×Nrep÷Vpp)*0.5
}

⍝ Calibration from Section 42.2.3:
⍝   sigma2=0.04, rho=0.4, N_rep=80, Vpp=2.5
⍝   OptExpScale 0.04 0.4 80 2.5   ⍝ → ~1.03 (one city-scale experiment optimal)

⍝ ── 4. POLICY DIFFUSION R₀ ──────────────────────────────────────────
⍝ R₀ = (A_bar × mu_demonstrated) / (theta × (1 - x_bar))
⍝ Policy goes viral when R₀ > 1

PolicyR0 ← {        ⍝ ⍵ = (A_bar, mu_demonstrated, theta, x_bar)
    A mu theta xbar ← ⍵
    (A×mu) ÷ theta×1-xbar
}

⍝ Preston Model calibration (Section 42.5):
⍝   A_bar   = 0.35 (UK local authority learning network)
⍝   mu      = 0.45 (demonstrated income premium ~4pp → normalised)
⍝   theta   = 0.30 (moderate political resistance)
⍝   x_bar   = 0.02 (2% of UK councils have adopted so far)
⍝ PolicyR0 0.35 0.45 0.30 0.02   ⍝ → ~0.54 (not yet viral — needs intervention)

⍝ ── 5. TIPPING THRESHOLD ─────────────────────────────────────────────
⍝ x_hat = (c - mu_beta) / v
TippingThreshold ← { (⊃⍵-1↓⍵) ÷ ⊃2↓⍵ }    ⍝ ⍵=(c, mu_beta, v)
⍝ Or more clearly:
TipThresh ← {           ⍝ ⍵=(switching_cost, avg_intrinsic_utility, network_coeff)
    c mu v ← ⍵
    (c-mu)÷v
}

⍝ ══════════════════════════════════════════════════════════════════
⍝ END CH42 APL
⍝ Key checks:
⍝   OptExpScale 0.04 0.4 80 2.5   ⍝ → ~1.03 (one experiment sufficient)
⍝   PolicyR0 0.35 0.45 0.30 0.02  ⍝ → ~0.54 (below viral threshold)
⍝   TipThresh 0.40 0.15 0.30      ⍝ → 0.833 (83% adoption needed)
⍝ ══════════════════════════════════════════════════════════════════
