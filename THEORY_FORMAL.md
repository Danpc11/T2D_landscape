# Formal framework: definitions, propositions and proofs

This document states precisely what the pipeline estimates and proves the
properties that the inferential claims rest on. `THEORY.md` gives the
conceptual framing; here we make it a theory in the mathematical sense:
objects, statements, proofs, and the assumptions each needs. Proofs are
complete where the result is ours and sketched where it is a direct
application of a published theorem (then cited).

Notation. T tissues, P common genes, stages s ∈ {0, 1, 2}. For tissue t and
stage s, W_ts ∈ ℝ^{P×P} is the weighted coexpression network (|bicor|^β,
zero diagonal), D_ts its strength matrix, L_ts = I − D^{-1/2} W D^{-1/2} the
normalised Laplacian, and V_ts ⊂ ℝ^P the r-dimensional eigenspace of the r
smallest non-trivial eigenvalues of L_ts. X_ts ∈ ℝ^{P×r} denotes any
orthonormal basis of V_ts (X_tsᵀX_ts = I_r). Gr(r, P) is the Grassmannian.

---

## Part I — Cross-tissue coherence as a cellular sheaf

### Definition 1 (tissue sheaf)

Fix a stage s. Let K_T be the complete graph on the T tissues. The
**coherence sheaf** F_s assigns
- to each vertex t the stalk F_s(t) = ℝ^{P×r};
- to each edge {t,u} the stalk F_s(tu) = ℝ^{P×r};
- restriction maps ρ_{t→tu}(X) = X Q_t, where Q_t ∈ O(r) are fixed
  orthogonal matrices (the alignment).

Its coboundary δ: ⊕_t F_s(t) → ⊕_{t<u} F_s(tu) is δ(X)_{tu} = X_t Q_t − X_u Q_u,
the sheaf Laplacian is L_F = δᵀδ (Hansen & Ghrist 2019), and the **sheaf
energy** of the assignment X = (X_ts)_t is

  E_s(X) = ‖δX‖²_F / Σ_t ‖X_ts‖²_F = ‖δX‖²_F / (T r).

The last equality uses orthonormality of each X_ts.

### Proposition 1 (global sections ⇔ zero energy)

H⁰(F_s) = ker δ, and for orthonormal stalk data, E_s(X) = 0 iff
X_ts Q_t = X_us Q_u for all t,u, iff V_ts = V_us for all t,u.

*Proof.* ker δ = H⁰ by definition of sheaf cohomology of a cellular sheaf on a
graph. E_s = 0 iff every summand ‖X_t Q_t − X_u Q_u‖² vanishes. Right
multiplication by an orthogonal matrix does not change the column span, so
X_t Q_t = X_u Q_u implies span(X_t) = span(X_u). Conversely if V_t = V_u
there exists Q ∈ O(r) with X_t Q = X_u, so choosing Q_t = Q, Q_u = I gives a
zero summand; the choice of a common reference in the pipeline realises
this simultaneously for all t when all subspaces coincide. ∎

### Proposition 2 (the aligned energy is a Grassmann distance)

For orthonormal X_t, X_u with principal angles θ_1,…,θ_r between V_t and V_u,

  min_{Q∈O(r)} ‖X_t Q − X_u‖²_F = 2 Σ_{i=1}^r (1 − cos θ_i) = 2r − 2 ‖X_tᵀX_u‖_*,

where ‖·‖_* is the nuclear norm. In particular the minimum depends only on
(V_t, V_u) ∈ Gr(r,P)², is invariant under any change of orthonormal basis
(eigenvector sign flips and rotations within degenerate eigenspaces), and
√(·) is a metric on Gr(r,P) (the chordal/Procrustes distance).

*Proof.* Expand ‖X_t Q − X_u‖² = ‖X_t Q‖² + ‖X_u‖² − 2 tr(QᵀX_tᵀX_u) =
2r − 2 tr(QᵀM) with M = X_tᵀX_u. Let M = U Σ Vᵀ be an SVD; the singular values
of M are cos θ_i (definition of principal angles). By von Neumann's trace
inequality, max_{Q∈O(r)} tr(QᵀM) = Σ σ_i(M) = ‖M‖_*, attained at Q = U Vᵀ
(the Procrustes solution). Substituting gives the formula. Basis changes
X_t ↦ X_t R, R ∈ O(r), leave the minimum unchanged because Q ↦ RᵀQ ranges
over O(r). That √(2Σ(1−cos θ_i)) is a metric on the Grassmannian is
classical (it is the Frobenius distance between the projectors
X_tX_tᵀ and X_uX_uᵀ up to a factor √2). ∎

**Remark (common reference vs pairwise).** The pipeline aligns every X_ts to
one reference R (the iteratively averaged healthy embeddings) rather than
solving each pair separately. Since the pairwise Procrustes minimum is a
lower bound for any common choice of Q_t, E_s(common) ≥ (1/Tr) Σ_{t<u} d²(V_t,V_u),
with equality iff a single set of Q_t is simultaneously optimal for all
pairs. The common-reference energy is therefore a conservative (upper)
estimate of pairwise incoherence, and it is the quantity that makes the
sheaf a sheaf (one restriction map per vertex, not per pair).

### Proposition 3 (sampling bias and its correction)

Let W_ts be estimated from n samples by a correlation estimator with
‖Ŵ_ts − W_ts‖_op = O_P(n^{-1/2}) (true for Pearson and biweight
midcorrelation under finite fourth moments), and let γ_ts > 0 be the
eigengap separating the r-dimensional eigenspace V_ts of L_ts from the rest
of the spectrum. Then

  ‖sin Θ(V̂_ts, V_ts)‖_F ≤ C ‖Ŵ_ts − W_ts‖_op / γ_ts = O_P(n^{-1/2}),

and consequently

  𝔼[Ê_s(n)] = E_s + b_s / n + o(1/n),  with b_s ≥ 0 depending on
  (γ_ts, noise level, T, r).

Hence (i) energies estimated at different n are not comparable, and (ii)
the Richardson-type estimator  Ẽ_s = 2 Ê_s(n) − Ê_s(n/2)  satisfies
𝔼[Ẽ_s] = E_s + o(1/n).

*Proof.* The first inequality is the Davis–Kahan sin Θ theorem (Yu, Wang &
Samworth 2015 version, which needs only the gap of the population matrix),
applied to L̂_ts vs L_ts; the normalised Laplacian is a Lipschitz function
of W on the set where strengths are bounded away from zero, which gives the
constant C. For the expansion, write each aligned stalk as
X̂_t = X_t + Δ_t with ‖Δ_t‖_F = O_P(n^{-1/2}) (Δ_t is the perturbation
after the optimal rotation, so its first-order part is orthogonal to V_t).
Then

  ‖X̂_t Q_t − X̂_u Q_u‖² = ‖X_tQ_t − X_uQ_u‖² + 2⟨X_tQ_t − X_uQ_u, Δ_tQ_t − Δ_uQ_u⟩ + ‖Δ_tQ_t − Δ_uQ_u‖².

The cross term has zero mean to first order because the leading term of Δ_t
is a linear function of the zero-mean estimation error Ŵ − W; the last term
has expectation of order 𝔼‖Δ‖² = O(1/n) and is non-negative, which gives
b_s ≥ 0. Summing over pairs and dividing by Tr yields the expansion. For
(ii), 2(E + b/n) − (E + 2b/n) = E. ∎

**Consequences implemented.** Equal-n subsampling within each tissue (removes
the between-stage difference in b/n), and the columns `E_n_extrapolated`
and `sampling_bias_hat` = Ê(n/2) − Ê(n) ≈ b/n in `sheaf_coherence.py`. On
synthetic data b/n ≈ 0.08 against between-stage differences ≈ 0.2, i.e. the
bias is not negligible and must be reported next to any effect.

### Proposition 4 (what the two nulls test)

Assume the latent decomposition V_ts = V_s^{sys} ⊕ V_ts^{spec} (a systemic
subspace common to tissues and a tissue-specific complement), and that
gene labels are exchangeable within tissue under the correspondence null.
Then:
(a) the correspondence-permutation null (rows of X_ts permuted independently
for t ≥ 1) has E_s^{null} equal in distribution to the energy of subspaces
with the same per-tissue spectra and no shared systemic component, so
E_s ≪ E_s^{null} is evidence that V_s^{sys} ≠ 0 (prediction P0);
(b) the stage-label permutation null preserves each tissue's sampling
distribution and n, so under H0 "stages are exchangeable" any statistic
that is a function of the per-stage energies has the same law as under
the permutation; rejecting it is evidence that E_s depends on s, at fixed n
(predictions P3, P5).

*Proof.* (a) Permuting gene rows of X_ts by a permutation matrix Π_t maps
V_ts to Π_t V_ts. If V_ts contained a component common to all tissues, it is
destroyed for t ≥ 1 with probability one (a generic permutation does not
preserve a fixed subspace), while the spectrum of Π_tL_tsΠ_tᵀ equals that of
L_ts, so per-tissue structure is preserved. (b) Standard permutation-test
argument: under exchangeability of samples across stages within a tissue,
the joint law of the data is invariant under relabelling, hence so is the
law of any statistic; the equal-n design guarantees that the statistic's
dependence on n does not differ between the observed and permuted
configurations. ∎

---

## Part II — The quasi-potential landscape

### Definition 2

Let p be the density of patient states in the r-dimensional embedding and
U = −ln p the quasi-potential. Û_h = −ln p̂_h is its kernel estimate with
bandwidth h. For a Morse function U, its **basins** are the stable
manifolds of its minima and its **passes** are the index-1 saddles.

### Proposition 5 (Morse correspondence for sublevel persistence)

For a Morse function U on a compact domain, the 0-dimensional persistence
diagram of the sublevel filtration {U ≤ c} has one point per minimum m that
is not the global minimum, with birth U(m) and death U(σ_m), where σ_m is
the index-1 saddle at which the basin of m merges with a deeper basin. The
persistence U(σ_m) − U(m) is the height of the barrier that has to be
crossed to leave the basin of m.

*Proof.* Classical Morse theory: the homotopy type of {U ≤ c} changes only at
critical values; at a minimum a new component is created (birth of an H₀
class), at an index-1 saddle two components merge (death of the younger,
by the elder rule), and higher-index critical points do not affect H₀
(Milnor 1963; Edelsbrunner & Harer 2010, ch. VII). ∎

### Proposition 6 (stability and a guaranteed basin count)

Let d_B denote the bottleneck distance. Then (Cohen-Steiner, Edelsbrunner &
Harer 2007)

  d_B(Dgm₀(Û_h), Dgm₀(U_h)) ≤ ‖Û_h − U_h‖_∞,

where U_h = −ln(p * K_h) is the smoothed population potential. Let ε_α be
such that ℙ(‖Û_h − U_h‖_∞ ≤ ε_α) ≥ α. Then, with probability ≥ α, every point
of Dgm₀(Û_h) with persistence > 2ε_α is matched to a point of Dgm₀(U_h) with
positive persistence, i.e. corresponds to a true basin of U_h.

*Proof.* A bottleneck matching of cost ≤ ε moves births and deaths by at most
ε each, so persistence changes by at most 2ε; a point with persistence > 2ε
cannot be matched to the diagonal. ∎

**Implementation.** ε_α is estimated by the bootstrap of Fasy et al. (2014)
restricted to the 80 % densest points (where minima and saddles live; in the
tails log p̂ is unstable and uninformative). `n_basins_guaranteed` uses the
threshold τ = max(ln 2, 2ε̂_{0.9}); the decision `two_attractors` uses the
nominal threshold ln 2 with scale stability, because Prop. 6 is a
sufficient condition and, at the present n, guarantees only barriers
above ≈ 1.4 nats (density ratio > 4).

### Theorem 1 (consistency of the two-attractor decision)

Let the embedding be fixed and suppose p is continuous with compact support
and bounded away from zero on it. Let h = h_n → 0 with n h_n^d → ∞ (d ≤ 3
here), so that ‖Û_{h_n} − U‖_∞ → 0 in probability on the core. Then:

(a) if U has exactly two minima whose merging saddle satisfies
U(σ) − max(U(m₁), U(m₂)) > ln 2 for all bandwidths in a neighbourhood of
h_n, the procedure returns `two_attractors` = 1 with probability → 1;
(b) if U has a single minimum, the procedure returns `two_attractors` = 0
with probability → 1.

*Proof.* Uniform consistency of the kernel density estimator on a compact
core under the stated bandwidth conditions is standard (Giné & Guillou
2002), and −ln is Lipschitz on densities bounded away from zero, so
‖Û − U‖_∞ → 0 on the core. By Prop. 6, the persistence of every point of
Dgm₀(Û) converges to that of Dgm₀(U). In case (a) the unique off-diagonal
point of Dgm₀(U) has persistence > ln 2, hence for large n so does its
estimate, at each of the three bandwidths (stability across scale), giving
n_basins = 2 stably; and the BIC of a two-component mixture exceeds that of
one component with probability → 1 because the true density is not in the
one-component family while a two-component Gaussian mixture reduces the
Kullback–Leibler divergence to it (BIC consistency for nested model
selection, Keribin 2000, under the usual regularity). Both conditions of the
decision rule hold, so it returns 1. In case (b) Dgm₀(U) is empty off the
diagonal; all estimated persistences → 0 < ln 2, so n_basins = 1 at every
scale in both planes; the rule returns 0 regardless of the BIC. ∎

**Remark on the finite-n regime.** Theorem 1 is asymptotic in n at fixed
embedding. The honest statement for the real data (n = 33–118) is the
power analysis in `python/simulation/`: the decision is reliable when the
latent separation exceeds ≈ 2.5× the Marchenko–Pastur noise of the top
principal components, and indeterminate below.

### Proposition 7 (Kramers asymmetry under gradient dynamics)

If the population dynamics are the overdamped Langevin diffusion
dx = −∇U(x) dt + √(2D) dB_t (gradient dynamics, J = 0), then for small D the
escape rates from the healthy and disease basins satisfy

  k(healthy→T2D) / k(T2D→healthy) ≍ exp[(ΔU_{T2D→h} − ΔU_{h→T2D}) / D]

up to a prefactor of order one (Kramers 1940; Freidlin & Wentzell). Hence
the barrier asymmetry reported by `landscape.py` is, under this assumption,
the log-ratio of return to onset rates; without the assumption it is a
purely geometric statement about U.

*Proof.* Direct application of Kramers' formula k ≍ exp(−ΔU/D) to each basin. ∎

### Proposition 8 (potential–flux decomposition of the inferred drift)

Let π be the entropic optimal transport plan between the healthy and T2D
empirical distributions (Schrödinger bridge with reference diffusion of
variance ε), b_i = Σ_j π_ij x_j / a_i − x_i the induced displacement at each
healthy sample, and f ∈ ℝ^E the projection of b onto the edges of the kNN
graph of all samples. Then f decomposes uniquely as f = B₁φ + f^⊥ with
φ ∈ ℝ^n (a node potential) and f^⊥ ∈ ker B₁ᵀ (a divergence-free flow), and

  ‖f‖² = ‖B₁φ‖² + ‖f^⊥‖².

If the true dynamics are gradient (J = 0) and the bridge recovers the drift,
then f^⊥ → 0 as the graph refines; conversely ‖f^⊥‖/‖B₁φ‖ bounded away from
zero, and larger than under label permutation, is evidence of a
non-conservative component J.

*Proof.* The decomposition is the discrete Hodge (Helmholtz) decomposition on
a graph: im B₁ ⊕ ker B₁ᵀ = ℝ^E, orthogonal because B₁ᵀ f^⊥ = 0. Least
squares gives φ. The statement about J is the definition of J as the
non-gradient part of the stationary drift (Wang, Xu & Wang 2008) combined
with the consistency of graph gradients for smooth fields on refining kNN
graphs. ∎

**Caveat (stated, not proved away).** The bridge identifies the most likely
diffusion between two *populations*; identifying it with the drift of an
individual's trajectory requires the ergodicity assumption of `THEORY.md`
§1. Under violation, ‖f^⊥‖ can be non-zero without any non-equilibrium
dynamics (e.g. if cohorts differ in a covariate not regressed out). This is
why the flux statistic is reported against the label-permutation null and
after within-stage covariate adjustment.

---

## Part III — What is new here, precisely

1. The identification of cross-tissue coherence with a cellular sheaf whose
   energy is a Grassmann distance between spectral subspaces (Def. 1,
   Props. 1–2), with the consequence that the "transition state" is
   characterised cohomologically (loss of global sections) rather than by
   any single-tissue statistic.
2. The finite-sample theory of that energy (Prop. 3): the 1/n bias, why
   equal n is necessary, and the extrapolated estimator that removes it.
3. The two-level evidence for basins (Prop. 6, Theorem 1): a decision rule
   with a consistency theorem, and a separately reported guaranteed count
   with a confidence level from the stability theorem.
4. The combination, at the population level, of a Morse-theoretic landscape,
   a Fisher–Rao geometry along a clinical control, and an OT-inferred
   potential–flux decomposition (Props. 7–8), each with its assumption
   made explicit and its null defined.

What is *not* claimed: mechanism, causality, individual dynamics, or
"irreversibility" in the thermodynamic sense.

## References

Cohen-Steiner, Edelsbrunner & Harer 2007, Discrete Comput Geom 37:103.
Edelsbrunner & Harer 2010, *Computational Topology*, AMS.
Fasy, Lecci, Rinaldo, Wasserman, Balakrishnan & Singh 2014, Ann Stat 42:2301.
Freidlin & Wentzell, *Random Perturbations of Dynamical Systems*, Springer.
Giné & Guillou 2002, Ann IHP Probab Stat 38:907.
Hansen & Ghrist 2019, J Appl Comput Topol 3:315.
Keribin 2000, Sankhyā A 62:49.
Kramers 1940, Physica 7:284.
Milnor 1963, *Morse Theory*, Princeton.
Wang, Xu & Wang 2008, PNAS 105:12271.
Yu, Wang & Samworth 2015, Biometrika 102:315.
