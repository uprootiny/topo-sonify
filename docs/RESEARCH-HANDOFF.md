# Research Handoff: Modern Computational Topology for Sonification

## Context

You are continuing work on `topo-sonify`, a Rust-based topology-first sonification engine at `/Users/uprootiny/Nissan/topo-sonify/`. The codebase currently implements (43 passing tests):

**Topology core:**
- **Simplicial complexes** — sorted-vec simplex representation, VR construction with squared-distance matrix and merge-intersection clique enumeration, union-find for β₀
- **Homology over Z/2Z** — boundary matrices, bitset column reduction, Betti numbers via rank-nullity with β₀ via union-find
- **Standard persistence algorithm** (ELZ 2002) — full column reduction on filtration boundary matrix, exact persistence pairs, twist optimization (Chen-Kerber 2011)
- **Filtered complexes** — diameter-weighted VR filtration with correct (value, dim, lex) ordering
- **Hodge Laplacian** — L_k = ∂_{k+1}∂_{k+1}^T + ∂_k^T∂_k over ℝ, signed boundary matrices, Jacobi eigenvalue algorithm, spectral gap, kernel dimension = β_k
- **Discrete Morse theory** — vertex-function extension, unique-unpaired-face gradient pairing, critical simplex enumeration, Morse inequality verification
- **Diagram distances** — bottleneck distance (binary search + Kuhn matching), Wasserstein distance (Hungarian algorithm)
- **Persistence analysis** — diagrams, landscapes, entropy, total/max persistence

**Instruments**: BettiDrone, SimplicialPluck, MorseSynth
- **Effects**: PersistenceDelay, CurvatureReverb, BettiDistortion
- **WASM bridge** for in-browser Web Audio playback
- **Literate HTML docs** (4,700+ lines) with MathJax, proofs, code stitches

All 24 tests pass. The math is correct for the basics we've implemented.

## What This Session Needs

A thorough literature review and gap analysis to bring the topology layer up to **publication-grade SOTA**. The current implementation covers the fundamentals but misses several important modern developments.

## Specific Research Tasks

### 1. Persistent Homology Algorithms (Priority: HIGH)

**Current state**: Naive filtration sweep — O(n⁴) for n points.

**Research needed**:
- **Standard algorithm** for persistence (Edelsbrunner, Letscher, Zomorodian 2002) — matrix reduction with lowest-one
- **Twist optimization** (Chen & Kerber 2011) — clearing columns
- **Chunk algorithm** for parallelism
- **Ripser** (Bauer 2021) — the current SOTA for Vietoris-Rips persistent homology. What makes it fast? (apparent pairs, implicit matrix, cohomology approach)
- **GUDHI** and **Dionysus** — what algorithms do these canonical libraries use?
- Find the exact algorithm descriptions and pseudocode sufficient to implement in Rust.

### 2. Coefficients Beyond Z/2Z (Priority: MEDIUM)

**Current state**: Everything is mod 2.

**Research needed**:
- Smith normal form over Z for integer homology
- Torsion in homology groups — what does it mean physically?
- When does working over Z vs Z/2Z vs Q change the Betti numbers?
- Are there sonification-relevant examples where Z/2Z misses structure?

### 3. Extended Persistent Homology (Priority: HIGH)

**Research needed**:
- **Extended persistence** (Cohen-Steiner, Edelsbrunner, Harer 2009) — what is the extended filtration?
- **Zigzag persistence** (Carlsson & de Silva 2010) — for time-varying data
- **Multi-parameter persistence** — is there a computable version?
- **Persistence images** and **persistence surfaces** as alternative representations
- **Persistence-weighted kernel methods**

### 4. Discrete Morse Theory (Priority: HIGH)

**Current state**: MorseSynth uses a simple filtration sweep as a "height function."

**Research needed**:
- **Forman's discrete Morse theory** (1998, 2002) — discrete Morse functions on CW-complexes
- Gradient vector fields and critical simplices
- How to find optimal discrete Morse functions
- The **Morse complex** and its relationship to homology
- Discrete Morse theory for simplicial complexes specifically
- Connection to Morse-Smale complexes in TDA

### 5. Spectral Methods (Priority: MEDIUM)

**Research needed**:
- **Hodge Laplacian** on simplicial complexes — the k-dimensional Laplacian
- Eigenvalues of the graph Laplacian → spectral graph theory → spectral clustering
- **Cheeger inequality** and its relationship to the spectral gap
- Connection to the **heat kernel** on manifolds
- How eigenvalues of the Hodge Laplacian relate to topological features
- The **spectral approach to persistent homology** (de Silva, Morozov, Vejdemo-Johansson)

### 6. Sheaves and Cosheaves (Priority: LOW but exciting)

**Research needed**:
- **Cellular sheaves** on simplicial complexes
- **Sheaf Laplacian** and its eigenvalues
- Ghrist & co.'s work on sheaf-theoretic signal processing
- Could sheaf cohomology provide richer sonification parameters?

### 7. Topological Signal Processing (Priority: HIGH)

**Research needed**:
- **Topological signal processing** (Barbarossa & Sardellitti 2020)
- Signal processing on simplicial complexes (not just graphs)
- Higher-order interactions in signal processing
- **Hodge decomposition** of signals: gradient + curl + harmonic components
- The harmonic component lives in the homology → direct connection to our Betti numbers

### 8. Mapper Algorithm (Priority: MEDIUM)

**Research needed**:
- The **Mapper** algorithm (Singh, Mémoli, Carlsson 2007)
- How it constructs a simplicial complex from data using a lens function
- Relationship to Reeb graphs
- Could Mapper provide a different sonification pathway?

### 9. Optimal Transport and Topology (Priority: LOW)

**Research needed**:
- Wasserstein distance between persistence diagrams
- Sliced Wasserstein kernels for persistence
- Stability theorems and their quantitative forms

### 10. Sonification-Specific Literature (Priority: HIGH)

**Research needed**:
- **Existing topology-sonification work** — what has been published?
- Search for: "sonification topology", "persistent homology sound", "topological data analysis audio", "Betti numbers music"
- The **auditory display** community (ICAD conference proceedings)
- Hermann, Hunt, Neuhoff: "The Sonification Handbook" — relevant chapters
- Any work on mapping TDA invariants to audio parameters
- Psychoacoustic considerations: which mappings are perceptually meaningful?

## Deliverables

For each research area, produce:

1. **Key references** (authors, year, title, DOI/arXiv ID where available)
2. **Core ideas** summarized in 2-3 paragraphs
3. **Algorithms/pseudocode** where applicable (especially for items 1, 4, 5, 7)
4. **Gap analysis**: what our current implementation lacks vs. SOTA
5. **Implementation recommendations**: what to add to `topo-sonify`, in priority order
6. **Open questions** worth investigating

## Output Format

Write results to `/Users/uprootiny/Nissan/topo-sonify/docs/RESEARCH-NOTES.md` as a structured document with the sections above.

Also update the literate HTML docs where appropriate — add "Further Reading" sections to each chapter pointing to the papers found.

## Constraints

- Prefer primary sources (original papers) over survey articles
- Distinguish between mathematically proven results and heuristic/empirical findings
- For algorithms, we need enough detail to implement in Rust (not just "use library X")
- Focus on what is relevant to sonification — not everything in TDA
- The codebase is Rust-only, targeting WASM. No Python dependencies.
