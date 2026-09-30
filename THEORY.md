# Theoretical framework

## 0. In one sentence

Healthy metabolic homeostasis and type 2 diabetes are two attractors of the
multi-organ transcriptional system; impaired glucose tolerance is the
transition state between them. The project asks whether this structure,
predicted by the bistable physiological models of T2D (Topp et al. 2000; Ha,
Satin & Sherman 2016), leaves a measurable footprint in the transcriptomes of
the organs, and whether the transition is a single-tissue event or a loss of
coordination between tissues.

## 1. One object: the quasi-potential landscape

Let x be the (multi-tissue) transcriptional state of a person, evolving on the
time scale of years as a stochastic process with stationary distribution P. We
define the **quasi-potential** U(x) = −ln P(x) (Freidlin–Wentzell; Jin Wang's
quantitative Waddington landscape). Each patient in a cohort is a sample from
P. Everything we measure is a property of U or of the dynamics that generate
it.

**Declared assumption.** Patients are independent samples from the same
dynamics (population ergodicity), not an individual trajectory. Covariates
that shape P without being part of the progression (BMI, age, sex) are
regressed out **centred within each stage**: regressing them raw erases the
between-group difference, because they are collinear with stage.

### 1.1 Topology of U: how many valleys?

Stable states are minima of U; transition states are index-1 saddle points. By
Morse theory, the persistent homology of the sublevel sets of U identifies
both: each H₀ bar is born at a minimum and dies at the saddle that merges its
basin with another; the **persistence** is the barrier height. The statement
"two equilibrium states and a point of no return" translates into: the
persistence diagram of U has two long H₀ bars, and the death of the shorter one
is the barrier.

Pre-specified decision rule for "two attractors": persistent basins (> ln 2
nats, i.e. the density at the saddle is less than half of that of the smaller
valley) **stable across bandwidths** ∧ ΔBIC(2 vs 1 components) > 0. No single
test suffices: a continuum can produce two basins at a single bandwidth.

### 1.2 Geometry of U: where does it change fastest?

If P is parametrised by a clinical control θ (HbA1c, glucose), the Fisher–Rao
metric g(θ) = E[(∂_θ ln P)²] measures how fast the distribution changes. It is
the unique natural metric on the space of distributions (Chentsov) and it is
maximised at phase transitions without requiring an order parameter.
Prediction: the peak of g falls in the IGT range. Without a continuous
covariate we use the critical-transition index I_c (Mojtahedi, Huang et al.
2016): the ratio of mean gene–gene correlation to mean patient–patient
correlation, which increases in the transition state.

### 1.3 Dynamics: potential and flux

The stationary dynamics decompose (Wang) as F = −D∇U + J, with J the
non-conservative probability flux. U describes the valleys; J breaks detailed
balance and makes the forward and return paths differ. With cross-sectional
data J is not directly observable, but with two ordered distributions (healthy,
T2D) the **Schrödinger bridge** (entropic optimal transport) gives the most
likely diffusion dynamics carrying one into the other, and from it a drift
b(x). The **Hodge decomposition** of b on the patient graph separates its
gradient part (dynamical −∇U) from its rotational + harmonic part (J). The
ratio ‖J‖/‖∇U‖ measures the distance from equilibrium; the barrier asymmetry
ΔU(T2D→healthy) − ΔU(healthy→T2D) measures, under gradient dynamics (Kramers),
why prevention is easier than reversal. Remission after bariatric surgery is
consistent with this: the disease valley is stable but can be left with a
large perturbation. We do not speak of "no return" but of a high return
barrier.

### 1.4 Gluing across tissues: one landscape or four?

Each tissue observes a projection P_t of the systemic state. Whether a single
systemic landscape exists of which the tissues are shadows is the question of
whether a **cellular sheaf** over the tissue graph has a global section. Stalk
of tissue t = position of each gene in the spectral embedding of its
coexpression network; restriction maps = O(r) rotations (Procrustes) onto a
common reference defined by the healthy state. The sheaf energy

E_s = Σ_{t<u} ‖X_ts − X_us‖² / Σ_t ‖X_ts‖²

is a normalised Grassmann distance between the spectral subspaces of the
tissues: zero if all share the same organisation, increasing with divergence.
E^Δ_s applies it to the changes relative to the healthy state. Its per-gene
decomposition localises genes whose reorganisation is tissue-specific. The
project's own hypothesis: **in the transition state the sheaf loses its global
sections**; the transition is a decoupling between organs, not a single-tissue
event.

### 1.5 Homeostasis

The healthy glycaemic plateau is an **infinitesimal homeostasis** in the sense
of Golubitsky and Stewart: a point where the derivative of the input–output
function (insulin resistance → glycaemia) vanishes, classifiable by singularity
theory. Its breakdown at a fold with hysteresis is the transition that the
bistable physiological models describe. This is the vocabulary for
"homeostatic state"; no thermodynamics is invoked.

## 2. Predictions (pre-specified)

| | Prediction | Statistic | Module |
|---|---|---|---|
| P0 | Cross-tissue structure exists in the healthy state | E_healthy ≪ correspondence null | `sheaf_coherence.py` |
| P1 | Two valleys in the tissues with power | `two_attractors` = 1 | `landscape.py` |
| P2 | IGT is the pass | I_c maximal in IGT; Fisher peak; IGT position on the route ≈ saddle | `landscape.py` |
| P3 | The pass is systemic | E_IGT, E^Δ_IGT maximal (`int_max`, `delta_int_gt_T2D`) | `sheaf_coherence.py` |
| P4 | T2D is an attractor, not noise | I_c and coherence recovered in T2D; modules differ from healthy | `02` + `sheaf` |
| P5 | High return barrier and flux | barrier asymmetry (bootstrap CI); ‖J‖/‖∇U‖ against null | `landscape.py` |
| S | Geometry of the valleys | EG, Ḡ, R̄, Q relative to the configuration null, at equal n | `02`–`05` |

Each prediction has three possible outcomes: support, refutation,
indeterminate (insufficient power). The power analysis with the real n
(`python/simulation/`) fixes in advance what is detectable: P3 requires
≥ 25–35 % of genes with tissue-specific reorganisation; P1 and P5 depend on
the separation between valleys relative to the Marchenko–Pastur noise of the
embedding.

## 3. Methodological safeguards

- n equalised across stages **within each tissue** (not globally).
- Network β unique per tissue (healthy state). In the sheaf, a single global β.
- Network metrics relative to the configuration null; permutation nulls for
  stage labels and for gene↔gene correspondence.
- Covariates centred by stage; sensitivity analysis without covariates.
- Liver (n = 5/4/9) excluded from P2–P3; only P1 (healthy vs T2D) and
  sensitivity.
- Sensitivity: leave-one-tissue-out, r ∈ {4, 8, 12}, β ∈ {4, 6, 8}, bandwidth
  ×{0.7, 1, 1.4}, bridge ε.

## 4. Limits

Cross-sectional design: no dynamics are observed; "transition" means
difference between stages, and J is inferred under a diffusion model.
Correlation networks: nothing is mechanistic; per-gene scores are topological
sensitivities. Different cohorts per tissue: the coherence we measure is the
one that survives that heterogeneity. And the sheaf energy also grows if each
tissue simply becomes noisier: the correspondence null and E^Δ control for
this only in part.

## Anchor references

Topp et al. 2000 J Theor Biol; Ha, Satin & Sherman 2016 Endocrinology;
Huang, Ernberg & Kauffman 2009; Mojtahedi et al. 2016 PLoS Biol; Chen, Liu &
Aihara 2012 Sci Rep; Wang, Xu & Wang 2008 PNAS; Freidlin & Wentzell;
Golubitsky & Stewart 2017 J Math Biol; Antoneli et al. 2025 Math Biosci;
Prokopenko et al. 2011 Phys Rev E; Schiebinger et al. 2019 Cell (Waddington-OT);
Hansen & Ghrist 2019 (sheaf Laplacians); Edelsbrunner & Harer (persistence).
