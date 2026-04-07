# Literature Grounding and Motivation

This document traces each module in `topo-sonify` to the publication, theorem, or
mathematical tradition that motivated its design. It is intended as a starting point
for a deeper lit review session with web access.

---

## 1. Simplicial Complexes (`src/topology/simplex.rs`, `complex.rs`)

### What we implement
- Abstract simplicial complexes as sorted vertex lists
- Face enumeration, closure, f-vector, Euler characteristic
- Star, link, join operations
- Vietoris-Rips construction from point clouds

### Grounding

The foundations are classical and appear in every algebraic topology textbook.
Our specific computational choices follow the tradition of:

**Edelsbrunner & Harer, _Computational Topology_ (2010), Ch. I-IV.**
The VR construction (Def. 1.21 in our docs) is defined as the flag complex
of the proximity graph. Our implementation uses the flag property explicitly:
we build the 1-skeleton first, then enumerate cliques via sorted adjacency
list intersection. This is the standard approach in computational topology
(cf. Zomorodian, _Topology for Computing_, 2005, §4.3).

**Missing depth:** We do not implement the **Čech complex** (which requires
computing circumradii, not just pairwise distances), nor **alpha complexes**
(Edelsbrunner 1993, which use the Delaunay triangulation for efficiency in
low ambient dimension). A proper lit review should cover:

- Edelsbrunner, E.P. Mücke, "Three-dimensional alpha shapes," _ACM Trans. Graphics_ 13(1), 1994.
- de Silva, V., Carlsson, G., "Topological estimation using witness complexes," _SPBG_ 2004.
  (Witness complexes: a landmark-based alternative to VR for large datasets.)
- Kerber, M., Sharathkumar, R., "Approximate Čech complex in low and high dimensions," _SoCG_ 2013.

**Open question for the review:** What is the state of the art for VR
construction in high ambient dimension? Ripser (Bauer 2021) uses implicit
representation — can we adopt a similar approach?

---

## 2. Homology (`src/topology/homology.rs`)

### What we implement
- Betti numbers β₀, β₁, β₂ via rank-nullity over Z/2Z
- β₀ via union-find (O(n·α(n)))
- Chain complexes, boundary matrices

### Grounding

**Poincaré, "Analysis Situs" (1895)** introduced homology as "connectivity numbers."
**Noether** (via Alexandroff, 1926) reformulated as group theory.
Our implementation follows **Munkres, _Elements of Algebraic Topology_ (1984), Ch. 5-8**.

The union-find optimization for β₀ is standard in computational geometry.
The Z/2Z coefficient choice avoids signs and Smith normal form, following
the practical recommendation of **Edelsbrunner & Harer (2010), §IV.2**.

**Missing depth:** We do not implement:

- **Smith normal form** for integer coefficients (detecting torsion).
  Key reference: Munkres (1984) Ch. 11; computational aspects in
  Dumas, Heckenbach, Saunders, Welker, "Computing simplicial homology
  based on efficient Smith normal form algorithms," _AAECC_ 2003.
- **Reduced homology** (H̃_k vs H_k) — differs only at β₀.
- **Relative homology** H_k(K, L) for a pair (K, L).
- **Cohomology** H^k — the dual theory, used by Ripser for efficiency.

**Open question:** The relationship between homology and cohomology
persistence (they produce the same diagram but different algorithms
have different complexity). See de Silva, Morozov, Vejdemo-Johansson,
"Dualities in persistent (co)homology," _Inverse Problems_ 2011.

---

## 3. Standard Persistence Algorithm (`src/topology/standard_persistence.rs`)

### What we implement
- The ELZ algorithm (Edelsbrunner, Letscher, Zomorodian 2002):
  column reduction on the filtration boundary matrix
- The twist optimization (Chen & Kerber 2011): clearing paired columns
- Sorted-vector column representation with merge-XOR for Z/2Z addition

### Grounding

**Edelsbrunner, Letscher, Zomorodian, "Topological persistence and
simplification," _DCG_ 28(4), 2002.** The foundational paper. Our
`compute_persistence()` is a direct implementation of Algorithm 1.

**Zomorodian, Carlsson, "Computing persistent homology," _DCG_ 33(2), 2005.**
The structure theorem for persistence modules: every persistence module
over a PID decomposes uniquely into interval modules. This is the algebraic
reason persistence diagrams are well-defined.

**Chen, Kerber, "Persistent homology computation with a twist," _EuroCG_ 2011.**
Our `compute_persistence_twist()` implements the clearing optimization.
When column j has pivot at row i, we clear column i since σ_i is now
known to be paired. This avoids redundant reductions.

**Missing depth:**

- **Cohomology approach** (used by Ripser): work with the coboundary matrix
  instead of the boundary matrix. The coboundary is sparser for VR complexes
  because each simplex has fewer cofaces than faces in the relevant range.
  Reference: Bauer, "Ripser: efficient computation of VR persistence barcodes,"
  _JACT_ 5, 2021.

- **Apparent pairs** (Bauer 2021): pairs that can be identified without
  any column reduction, just from the combinatorial structure.

- **Implicit matrix representation**: Ripser never materializes the boundary
  matrix — it computes columns on the fly from the combinatorics of the VR
  complex. This reduces memory from O(m²) to O(m).

- **Clearing + compression** (Bauer, Kerber, Reininghaus, Wagner 2017):
  combining clearing with column compression for further speedup.

**Critical gap:** Our current implementation materializes all simplices and
boundary columns. For a VR complex on n points up to dimension d, this is
O(n^{d+1}) simplices. Ripser avoids this by using the flag complex property
to compute boundary columns implicitly.

---

## 4. Filtered Complex (`src/topology/filtration.rs`)

### What we implement
- `FilteredComplex`: simplices ordered by (filtration_value, dimension, lexicographic)
- VR filtration from point clouds with diameter-based simplex weights
- Boundary-index lookup for the persistence algorithm

### Grounding

The concept of a filtration (nested sequence K₀ ⊆ K₁ ⊆ ... ⊆ K_n) is
from **Edelsbrunner & Harer (2010), §VI.1**. Our ordering convention
(lower dimension first at equal filtration values) ensures that faces
precede cofaces, which is required for the correctness of the column
reduction algorithm.

**Missing depth:**

- **Sublevel set filtrations** for functions on manifolds (not just distance-based).
- **Lower-star filtration** (used in our discrete Morse construction, but not
  exposed as a general-purpose filtration type).
- **Rips filtration as a 1-parameter family**: the theoretical framework of
  **persistence modules as functors** (Bubenik, Scott, "Categorification of
  persistent homology," _DCG_ 2014).
- **Multi-parameter persistence**: filtrations indexed by ℝ^d rather than ℝ.
  This is an active research frontier — no complete discrete invariant exists
  (Carlsson & Zomorodian 2009). Key reference: Lesnick, "The theory of the
  interleaving distance on multidimensional persistence modules," _Found.
  Comput. Math._ 2015.

---

## 5. Hodge Laplacian (`src/topology/hodge.rs`)

### What we implement
- L_k = ∂_{k+1} ∂_{k+1}^T + ∂_k^T ∂_k (the k-th combinatorial Laplacian)
- Signed boundary matrices over ℝ (not Z/2Z — signs matter for spectra)
- Eigenvalue computation via Jacobi iteration
- Spectral gap (algebraic connectivity)
- Kernel dimension = β_k (discrete Hodge theorem)

### Grounding

**Eckmann, "Harmonische Funktionen und Randwertaufgaben in einem Komplex,"
_Comment. Math. Helv._ 17, 1944-45.** The original discrete Hodge theory.

**Lim, "Hodge Laplacians on graphs," _SIAM Review_ 62(3), 2020.** The
modern treatment extending graph Laplacians to simplicial complexes. Our
implementation follows Lim's exposition.

**Chung, _Spectral Graph Theory_, CBMS 92, 1997.** For k=0, L_0 is the
graph Laplacian. The spectral gap λ_1 (Fiedler value) measures algebraic
connectivity. The Cheeger inequality relates it to edge expansion.

**Missing depth:**

- **Hodge decomposition**: every k-chain decomposes as
  c = ∂_{k+1} α + ∂_k^* β + h, where h is harmonic (in ker L_k).
  The three components are: gradient, curl, harmonic.
  Reference: Barbarossa & Sardellitti, "Topological signal processing over
  simplicial complexes," _IEEE TSP_ 68, 2020.

- **Heat kernel on simplicial complexes**: K_t(x,y) = Σ_i e^{-λ_i t} φ_i(x) φ_i(y).
  Controls diffusion on the complex. Relevant for reverb modeling.

- **Persistent Laplacians** (Wang, Nguyen, Wei, "Persistent spectral graph,"
  _IJNAM_ 17, 2020): track how the spectrum of L_k changes across the filtration.
  This is a spectral analogue of persistent homology — our sonification could
  use persistent eigenvalues instead of Betti numbers for richer audio control.

- **Normalized Laplacians** and their spectral properties.

---

## 6. Discrete Morse Theory (`src/topology/morse.rs`)

### What we implement
- Discrete Morse functions from vertex-function extensions (f(σ) = max vertex value)
- Gradient vector fields via the "unique unpaired face" rule
- Critical simplex enumeration
- Morse inequality verification (weak and strong)

### Grounding

**Forman, "Morse theory for cell complexes," _Adv. Math._ 134(1), 1998.**
The foundational paper. Defines discrete Morse functions, gradient vector
fields, and proves the discrete Morse inequalities.

**Forman, "A user's guide to discrete Morse theory," _Séminaire Lotharingien_ 48, 2002.**
An accessible exposition. Our "unique unpaired face" construction follows
the standard greedy algorithm described here.

**Missing depth:**

- **Optimal discrete Morse functions**: finding a discrete Morse function that
  minimizes the number of critical simplices is NP-hard (Joswig & Pfetsch 2006).
  But good heuristics exist: the **ProcessLowerStars** algorithm of
  Robins, Wood, Sheppard ("Theory and algorithms for constructing discrete
  Morse complexes from grayscale digital images," _IEEE PAMI_ 33, 2011).

- **Morse-Smale complexes**: the partition of the complex into cells based on
  the gradient flow. Each cell is associated with a pair of critical simplices.
  Reference: Edelsbrunner, Harer, Natarajan, Pascucci, "Morse-Smale complexes
  for piecewise linear 3-manifolds," _SoCG_ 2003.

- **Persistent homology via discrete Morse theory**: reducing the complex to
  its Morse complex before computing persistence. The Morse complex has far
  fewer cells (c_k instead of f_k), dramatically speeding up persistence.
  Reference: Mischaikow, Nanda, "Morse theory for filtrations and efficient
  computation of persistent homology," _DCG_ 50, 2013.

- **Connection to our MorseSynth instrument**: the current MorseSynth uses
  filtration value as a proxy for Morse function. A proper implementation
  would use actual critical simplices from the discrete Morse function to
  generate events, with the Morse index determining the voice character.

---

## 7. Diagram Distances (`src/topology/distances.rs`)

### What we implement
- Bottleneck distance via binary search + Kuhn bipartite matching
- Wasserstein distance (p-th order) via Hungarian algorithm
- Augmented bipartite graph with diagonal matching

### Grounding

**Cohen-Steiner, Edelsbrunner, Harer, "Stability of persistence diagrams,"
_DCG_ 37(1), 2007.** Proves d_B(Dgm(f), Dgm(g)) ≤ ||f-g||_∞.
This is the stability theorem that makes TDA work on noisy data.

**Cohen-Steiner, Edelsbrunner, Harer, Morozov, "Persistent homology for
kernels, images, and cokernels," _SoCG_ 2009.** Extended stability.

**Our implementation:** bottleneck via binary search on candidate thresholds
with Hopcroft-Karp-style augmenting path matching. Wasserstein via the
Kuhn-Munkres (Hungarian) algorithm. Both are standard.

**Missing depth:**

- **Sliced Wasserstein distance**: a fast approximation based on 1D projections.
  Reference: Carrière, Cuturi, Oudot, "Sliced Wasserstein kernel for
  persistence diagrams," _ICML_ 2017.

- **Persistence images** (Adams et al., _JMLR_ 18, 2017): a vectorized
  representation via Gaussian-weighted sums. Enables machine learning
  on persistence diagrams.

- **Persistence landscapes** (Bubenik, _JMLR_ 16, 2015): functional
  representation in a Banach space. We implement `landscape()` in
  persistence.rs but do not use it for distances.

---

## 8. Sonification Mapping (`src/mapping/`, `src/instruments/`, `src/audio/`)

### What we implement
- β_k → oscillator count, modulation depth, sub-oscillators
- f-vector → just intonation ratios
- χ → reverb character (via Gauss-Bonnet)
- Persistence pairs → delay taps
- Filtration sweep → Morse-style musical score
- Karplus-Strong string synthesis, biquad filters, Schroeder reverb

### Grounding

The sonification mapping is our original contribution. There is **very little
prior work** on topology-driven sonification specifically. The closest references:

**Hermann, Hunt, Neuhoff (eds.), _The Sonification Handbook_, Logos 2011.**
General sonification design principles. Our parameter mapping approach
follows the taxonomy in Ch. 15 (Grond & Berger).

**The ICAD proceedings** (International Conference on Auditory Display):
search for "topology" and "homology" yields almost no results as of 2024.
This is genuinely unexplored territory.

**Audio DSP references:**
- Karplus, Strong, "Digital synthesis of plucked-string and drum timbres,"
  _CMJ_ 7(2), 1983.
- Smith, J.O., _Physical Audio Signal Processing_, CCRMA online book.
- Välimäki et al., "Discrete-time modelling of musical instruments,"
  _Rep. Prog. Phys._ 69(1), 2006.

---

## Summary of Gaps for Lit Review

| Priority | Topic | Key Missing References |
|----------|-------|----------------------|
| HIGH | Ripser algorithm (implicit VR, apparent pairs, cohomology) | Bauer 2021 |
| HIGH | Persistent Laplacians | Wang, Nguyen, Wei 2020 |
| HIGH | Hodge decomposition for signal processing | Barbarossa & Sardellitti 2020 |
| HIGH | Multi-parameter persistence | Lesnick 2015, Carlsson & Zomorodian 2009 |
| MEDIUM | Morse-accelerated persistence | Mischaikow & Nanda 2013 |
| MEDIUM | Persistence images/landscapes for ML | Adams et al. 2017, Bubenik 2015 |
| MEDIUM | Alpha/witness complexes | Edelsbrunner 1994, de Silva & Carlsson 2004 |
| MEDIUM | Sliced Wasserstein | Carrière, Cuturi, Oudot 2017 |
| LOW | Sheaf Laplacians | Hansen & Ghrist 2019 |
| LOW | Zigzag persistence | Carlsson & de Silva 2010 |
| LOW | Cubical complexes | Kaczynski, Mischaikow, Mrozek 2004 |

## Recommended Lit Review Protocol

1. **Start with surveys**: Otter, Porter, Tillmann, Grindrod, Harrington,
   "A roadmap for the computation of persistent homology," _EPJ Data Sci._ 2017.
   This maps the entire field and identifies which algorithms apply where.

2. **Deep-dive Ripser**: Read Bauer (2021) in full. The key ideas (apparent
   pairs, emergent pairs, implicit matrix, cohomology) are the single most
   impactful upgrade path for our codebase.

3. **Spectral topology**: Read Lim (2020) and Barbarossa & Sardellitti (2020)
   back-to-back. The Hodge decomposition of signals on simplicial complexes
   is the theoretical bridge between our topology and audio.

4. **Sonification literature**: Search ICAD proceedings + Google Scholar for
   "topology sonification," "persistent homology audio," "TDA sound."
   Identify if anyone has published in this specific intersection.
