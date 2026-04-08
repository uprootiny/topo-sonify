# Methodology: Building topo-sonify

A record of how this project was designed, built, and iterated — the decision
sequence, the mistakes, and the principles that emerged.

---

## Phase 1: Math core first, audio second

### Decision
Build the entire topology layer before writing a single line of audio code.
Simplicial complexes, homology, boundary matrices, Betti numbers — all tested
against known mathematical truths (β(S²) = (1,0,1), χ = 2) before any sound
is produced.

### Why this worked
The topology-to-audio mapping is the core intellectual contribution. If the
invariants are wrong, the sonification is meaningless noise. By getting β₀,
β₁, β₂ correct first, the instrument design became straightforward: each
invariant had a clear musical interpretation that fell out of the math.

### What it cost
The project produced no sound for the first several hours. This required
trust that the math would pay off. It did.

---

## Phase 2: Instruments and effects as mathematical consequences

### Decision
Each instrument and effect module is parameterized by a specific topological
invariant, with a written justification for why the mapping is natural (not
arbitrary).

### The mapping table

| Invariant | Audio parameter | Why it's natural |
|-----------|----------------|------------------|
| β₀ (components) | oscillator count | Independent components can't interact, like independent voices |
| β₁ (loops) | LFO modulation | A cycle is periodic; periodicity is oscillation |
| β₂ (voids) | sub-oscillators | Enclosed cavities resonate, adding depth |
| f-vector ratios | just intonation | Combinatorial density → harmonic density |
| χ (Euler char) | reverb time/color | Gauss-Bonnet: curvature determines wave propagation |
| persistence pairs | delay taps | Feature lifetime → echo duration |
| filtration sweep | musical score | Morse theory: topology changes = musical events |

### What this prevented
No ad-hoc parameter mappings. Every knob has a theorem behind it.

---

## Phase 3: Literate documentation as a design tool

### Decision
Write the HTML seminar docs alongside the code, not after. Each chapter
develops the mathematical theory first (with full definitions, proofs, and
worked examples), then presents the code as the computational embodiment
of that theory.

### What this caught
Writing the ∂²=0 proof forced us to verify the sign convention in the
boundary matrix. Writing the Euler-Poincaré section revealed that
`euler_characteristic()` was computing χ from the f-vector but not
cross-checking against Σ(-1)^k β_k. The docs became a second test suite.

### The ratio
5,766 lines of docs for 4,200 lines of Rust. More documentation than code.
This is intentional: the code is the distillation; the docs are the argument.

---

## Phase 4: The data structure gauntlet

### The progression
1. `Simplex` as `BTreeSet<usize>` — tree-allocated, slow hash, slow compare
2. → sorted `Vec<usize>` — contiguous, O(k) everything, 50% less allocation
3. `SimplicialComplex` with linear-scan `simplices_of_dim()` — O(n) per call
4. → dimension-indexed cache with cloned simplices — 2× memory
5. → index-based cache (`Vec<Vec<usize>>`) — no cloning, O(1) lookup
6. `insert()` recursive — stack overflow on dim ≥ 10
7. → iterative worklist — bounded stack
8. `boundary_matrix` with `Vec<bool>` columns — slow XOR
9. → u64 bitset columns — 64× throughput
10. `vietoris_rips` with brute-force quadruple loop — O(n⁴)
11. → flat distance² matrix + sorted adjacency + merge-intersection — O(n² + output)

### The lesson
Every data structure choice was wrong on the first try. Not because we didn't
know better, but because "get it working first" defers the structural decision
to a point where it's harder to change. The fix each time was the same:
**use the simplest representation that supports the operations you actually need.**

---

## Phase 5: The persistence algorithm, done right

### The progression
1. **Betti tracking**: rebuild VR complex at each filtration step, diff Betti numbers.
   Problem: O(n⁴ × steps), doesn't produce exact pairs, discretization artifacts.

2. **Critical-radius sampling**: only rebuild at edge-insertion radii.
   Better: no discretization. Still wrong: rebuilding from scratch each time.

3. **Standard persistence (ELZ)**: single column reduction on the full filtration
   boundary matrix. Exact pairs. O(m³) worst case but fast in practice.

4. **Twist optimization (Chen-Kerber)**: clear columns of paired creators.
   2-10× speedup for free.

### The disconnect that persisted
Even after implementing the standard algorithm, the sonifier's `from_point_cloud`
still called the old Betti-tracking path. The new algorithm existed but wasn't
wired in. This went unnoticed for two review rounds because:
- Tests passed (they tested the new algorithm directly, not through the sonifier)
- The sonifier's output "sounded OK" (the old approximation was close enough)

Found during a full rewire audit. **Test the integration, not just the units.**

---

## Phase 6: WASM and the last mile

### The build chain
```
Rust source → cargo + wasm-pack → .wasm + .js bindings → GitHub Actions → Pages
```

### What was trivial
- `wasm-pack build --target web` just worked
- 158KB binary, loads instantly
- `ScriptProcessor` feeding Rust-generated samples to Web Audio

### What was not trivial
- **Import paths**: `../pkg/` works locally (web/ is a subdirectory) but not on
  Pages (everything is flat). Solution: try `./pkg/` first, fallback to `../pkg/`.
- **`script type="module"` with top-level await**: works in all modern browsers
  but silently fails in older ones. Acceptable tradeoff.
- **AudioContext autoplay policy**: browsers block audio until user gesture.
  The "Listen" / "Play" buttons handle this, but it's an extra interaction step.

---

## Phase 7: Interactive lessons

### The model
Ableton's Learning Music: 60-70% interactive area, 30% minimal text, progressive
disclosure, exploration before explanation.

### Design principles applied
1. **One concept per page.** Lesson 1: points. Lesson 2: connections. Lesson 3: loops.
2. **Interaction triggers text.** Add a point → first paragraph reveals. Add three → second. The text follows the discovery, not the other way around.
3. **Audio is immediate.** No "read this first" — click Listen, then explore.
4. **Stats overlay, not sidebar.** Betti numbers float over the canvas. You see them change as you interact. No context switch.
5. **Points are draggable.** In lessons 2-3, you can grab points and move them. The topology and sound update in real time.

### What's still missing
- Lesson 4 (persistence) needs the WASM engine for exact diagrams; currently uses JS approximation
- Lesson 5 (curvature) generates impulse responses in JS; should use the Rust reverb
- No touch support for mobile
- No guided "try this" prompts within lessons

---

## Phase 8: SOTA topology modules

### What we implemented
| Module | Algorithm | Reference |
|--------|-----------|-----------|
| Standard persistence | ELZ column reduction + twist | Edelsbrunner et al. 2002, Chen-Kerber 2011 |
| Hodge Laplacian | L_k = ∂_{k+1}∂_{k+1}^T + ∂_k^T∂_k, Jacobi eigenvalues | Eckmann 1944, Lim 2020 |
| Discrete Morse theory | Unique-unpaired-face pairing, critical simplices | Forman 1998 |
| Diagram distances | Bottleneck (binary search + Kuhn), Wasserstein (Hungarian) | Cohen-Steiner et al. 2007 |

### What tripped us on each

**Hodge eigenvalues.** QR with Wilkinson shifts: 150 lines, 3 bugs (tridiagonalization,
shift formula, deflation). Jacobi rotations: 50 lines, correct on first run. For small
matrices (which is all we need — the Laplacian of a complex with < 100 simplices),
Jacobi wins on simplicity.

**Discrete Morse.** The greedy "pair with lowest coface" heuristic violates Morse
inequalities. The standard construction requires checking for a *unique* unpaired face.
The difference: greedy pairs as many as possible; standard pairs only when there's
exactly one option. Greedy is too aggressive.

**Bottleneck distance.** Wrote the matching check, forgot to actually use its result.
`has_perfect_matching` returned `true` unconditionally. The binary search found the
smallest candidate threshold every time. Subtle because the output was always a
plausible-looking number.

---

## Recurring antipatterns

### 1. Piecemeal fixing
User says "there are algorithmic tricks you've skipped." Fix two things. User repeats.
Fix two more. User repeats again. Should have: read every file, listed every issue,
fixed them all in one pass.

### 2. Stub-then-forget
Write a stub implementation ("we'll replace this with the real algorithm later"),
get the tests passing, move on. The stub survives because nothing fails. The real
algorithm is implemented elsewhere but never wired in. Solution: **no stubs. Either
implement it or don't have it.**

### 3. Wrong abstraction level
BTreeSet for 3-element lists. Dense matrices for sparse operations. Recursive descent
for flat iteration. Each time, the abstraction was too heavy for the actual workload.
The fix is always: **look at the actual data sizes and access patterns, not the
theoretical generality.**

### 4. Testing the unit, not the integration
All 24 tests passed while the sonifier was using a completely different persistence
code path. The tests validated the new algorithm; the sonifier called the old one.
Solution: **test through the public API, not the internal modules.**

---

## By the numbers

| Metric | Count |
|--------|-------|
| Rust source lines | 4,200 |
| Tests | 44 |
| Clippy warnings | 0 |
| Seminar HTML docs | 5,766 lines |
| Interactive lessons | 6 |
| Runnable examples | 3 (680 lines) |
| WASM binary size | 158 KB |
| Git commits | 6 |
| CI builds | 6 (5 success, 1 cancelled) |
| Data structure rewrites | 5 |
| Algorithm rewrites | 4 |
| Bugs found by doc writing | 2 |
| Bugs masked by saturating_sub | 1 |
| Functions that returned true unconditionally | 1 |
