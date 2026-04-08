# Dev Notes: topo-sonify

What worked, what didn't, what kept tripping us up.

## What worked

**Math-first design.** Building the topology core before any audio code meant the instruments had solid invariants to map from. The Betti→voices, χ→reverb, persistence→delay mappings emerged naturally from the math rather than being forced.

**Sorted-vec simplices.** Switching from `BTreeSet<usize>` to sorted `Vec<usize>` for simplex representation was the single biggest performance win. Contiguous memory, O(k) hashing, O(n) face enumeration. Should have done this from the start.

**Union-find for β₀.** Computing connected components via matrix reduction is O(n³). Union-find is O(n·α(n)). Obvious in hindsight — β₀ is fundamentally a connectivity question, not a linear algebra question.

**Standard persistence algorithm.** The Betti-tracking approach (rebuild complex at each radius, compare Betti numbers) was a quick hack that gave approximate results. Replacing it with the real ELZ column reduction on the filtration boundary matrix gave exact persistence pairs. The twist optimization was trivial to add and measurably faster.

**WASM + Web Audio.** Compiling the full Rust engine to WASM and streaming audio via ScriptProcessor worked on the first try. 158KB binary, instant load.

**Literate style.** Writing the math docs forced us to think clearly about each algorithm before implementing it. Several bugs were caught during doc writing, not testing.

## What didn't work / kept tripping us

**Recursive face insertion.** The original `insert()` was recursive — for a k-simplex, it recursed k levels deep generating all faces. On high-dimensional simplices this blows the stack. Replaced with an iterative worklist. Should have been iterative from day one.

**Z/2Z column reduction — first attempt.** The initial column reduction only tried one elimination step per column (if the lowest entry collided with an existing pivot, it XOR'd once and moved on). Need to *loop* until either the column is zero or has a unique pivot. Classic off-by-one in the algorithm structure.

**BettiDrone frequency modulation.** First version mutated oscillator frequencies for LFO modulation (`set_frequency` every sample, then restore). This accumulates floating-point drift over millions of samples. Fixed by applying modulation as amplitude modulation on the mixed output instead.

**`saturating_sub` masking bugs.** When the Betti computation underflowed (negative result from wrong matrix rank), `saturating_sub` silently clamped to 0 instead of panicking. This hid the column reduction bug for several iterations. Should have used checked subtraction during development and only switched to saturating for production.

**Persistence diagram paths.** The `from_point_cloud` function used the old Betti-tracking hack even after we implemented the standard algorithm. The sonifier wired persistence into instruments via one path but computed it via another. This disconnect persisted across multiple "fix" rounds until a full rewire audit caught it.

**Dim cache cloning.** The dimension cache stored cloned simplices (`Vec<Vec<Simplex>>`) — doubling memory for no reason. Replaced with index-based caching (`Vec<Vec<usize>>` into a flat vec). The API mismatch between `simplices_of_dim(&self) -> Vec<&Simplex>` (immutable borrow, cache might be stale) and `simplices_of_dim_cached(&mut self) -> &[Simplex]` (mutable, always fresh) was confusing. Resolved by maintaining the cache incrementally during insert.

**Web asset paths.** `web/index.html` imported `../pkg/topo_sonify.js`. On GitHub Pages, the site is assembled flat (web/ contents copied to root), so `../pkg/` goes up one level too far. Needed a fallback: try `./pkg/` first, then `../pkg/`.

**Hodge eigenvalue algorithm.** First implementation used QR with Wilkinson shifts — a 150-line algorithm with subtle bugs in the tridiagonalization, shift computation, and deflation logic. It only returned 1 eigenvalue for a 3×3 matrix. Replaced with the Jacobi rotation algorithm (50 lines, always correct for symmetric matrices, converges reliably for small sizes).

**Discrete Morse pairing.** First attempt used a greedy "pair each face with its lowest coface" rule. This paired too many simplices — leaving fewer critical cells than Betti numbers, violating the Morse inequalities. Fixed with the standard "unique unpaired face" rule: only pair a simplex with a face if that face is its *unique* unpaired face.

**Bottleneck distance matching.** The `has_perfect_matching` function returned `true` unconditionally after checking diagonal reachability — never actually testing whether a valid bipartite matching existed at the given threshold. The binary search always converged to the smallest candidate distance regardless of the actual diagram structure. Fixed with a proper augmented bipartite graph and Kuhn's algorithm.

## Patterns

**Piecemeal fixes don't work.** Three rounds of "fix the algorithmic tricks" before doing a proper audit. Should have read all the code once, listed all issues, and fixed them in one pass.

**The simplest data structure wins.** BTreeSet → sorted Vec. HashMap → flat index. Recursive → iterative. Dense matrix → bitset columns. Every time, the simpler representation was both faster and more correct.

**Tests should encode mathematical truths.** The most valuable tests aren't "output is non-zero" — they're "β(S²) = (1,0,1)" and "χ = Σ(-1)^k β_k = Σ(-1)^k f_k". When these fail, the bug is always in the algorithm, not the test.

**The web frontend is the last mile.** Everything can be correct in Rust and still not work in the browser because of path issues, WASM loading, or AudioContext policies. Test the deployed artifact, not just `cargo test`.
