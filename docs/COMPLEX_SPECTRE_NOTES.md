# Complex Spectre & Dual-Plane 3-SAT Notes

> **Audit notes, added after an independent re-run.** The body below is the delivered text; bracketed markers such as [audit 1] point here. Details and numbers are in [experiments/04_complex_spectre/README.md](../experiments/04_complex_spectre/README.md).
>
> 1. **F_Φ is not holomorphic.** θ_k(z) is real-linear in (x, y), so e^{iθ_k(z)} also depends on z̄. Measured |∂F/∂z̄| is about as large as |∂F/∂z|, and roughly half the plane reverses orientation. The Stage 3 integral is therefore a topological degree: it counts vortices minus antivortices. The formula (1/2πi)∮F′/F dz and Cauchy's argument principle do not apply.
> 2. **Zeros of F_Φ are not satisfying assignments.** A zero needs the clause terms to cancel, not each term to vanish. None of 58 zeros across 40 n = 35 formulas rounds to a solution; zeros round to slightly worse than typical assignments. |F_Φ| is uncorrelated with the number of violated clauses (Spearman ≈ 0), while the nonnegative energy E is correlated (≈ 0.7–0.8). The clause weights e^{2πic/ρ} do not prevent cancellation.
> 3. **One z-plane shows few assignments.** The lines cos θ_k(z) = 0 cut the window into at most 1 + L + L(L−1)/2 cells. For `missed_sat_instance_35.cnf` (18 variables, not 35) the unit disk shows 128 of 262,144 assignments and none of its 9 solutions. The solver leaves the plane with its first flip, so Stages 3–4 only choose starting points. The Ω(2^{n/2}) frequency in Section 5 is right, and this counting bound proves it for the rounding map.
> 4. **The "Spectre mesh" is neither the Spectre nor a tiling.** The contour takes 13 unit steps and closes with an edge of length √7. Its copies sit on a golden-angle spiral and cover about 45% of the disk, with overlaps. The true Tile(1,1) is now in the code; swapping it in changes nothing measurable.
> 5. **Ablations.** Random starts in place of Stages 3–4 solve as many formulas, run 3–4× faster, and on 100 n = 35 formulas use fewer steps (p = 0.005). Replacing the logistic map with a PRNG makes no detectable difference. The SL(2, ℤ) fold count is computed but never used.
> 6. **Permutation result.** The reported formula is the most favorable of 20. Across 20 formulas, canonicalization lowers the mean on 13 and raises it on 7 (p = 0.26). Canonicalizing clause order as well makes every relabeling run identically.
> 7. **Critical "needles".** On 100 formulas, WalkSAT (Experiment 03's rule) matches the step counts (p = 0.23) and runs faster. The improvement over the deterministic golden sweep is the usual effect of adding noise.
> 8. **Chirality.** Random 3-SAT has no global spin-flip symmetry to break; NAE-SAT does. The solver also never uses tile chirality.
> 9. **Untested claims.** The Newton/Julia picture, the Beltrami vortices and "ghost basins" were not tested, and the Julia-set language needs a holomorphic map. RSB clustering near α_c is established physics; its link to F_Φ is not.

A mathematical notebook and post-mortem extending the `black-swan-jeweled-math` progression—from quasicrystal jewels (`01_initial_jewel`), Hat/Pinwheel adaptive tilings (`02_tiling_jewels`), and golden-ratio sweeps (`03_golden_sweep`)—into **complex dynamics, chiral monotiles, chaos universality, and dual 2D complex spaces ($\mathbb{C}_{\text{in}} \to \mathbb{C}_{\text{out}}$)**.

---

## 1. Constants That Alter Search Physics

| Constant | Value | Origin | Role in Search Physics |
| :--- | :--- | :--- | :--- |
| **Euler Identity** $e^{i\pi}+1=0$ | $-1$ | Complex unit circle | Promotes Boolean spins $\{-1,+1\}$ into smooth phasors $z_k = e^{i\theta_k}$; clauses become wave-interference products. |
| **Golden Ratio** $\phi$ | $1.6180339887\dots$ | $x^2 - x - 1 = 0$ | Optimal 2-term low-discrepancy angular packing ($\theta_g = 2\pi/\phi^2 \approx 137.508^\circ$). |
| **Plastic Number** $\rho$ | $1.3247179572\dots$ | $x^3 - x - 1 = 0$ | Smallest Pisot–Vijayaraghavan number; only other morphic number besides $\phi$ (Aarts et al., 2001). Matches 3-literal symmetry of 3-SAT; spaces complex clause weights $\omega_c = e^{i 2\pi c/\rho}$ to prevent false cancellation [audit 2]. |
| **Silver Ratio** $\delta_S$ | $2.4142135624\dots$ | $1 + \sqrt{2}$ | Governs 8-fold Ammann–Beenker quasicrystal tilings; $45^\circ$ corridors align cleanly with pairwise XOR/parity constraints. |
| **Feigenbaum** $\delta$ | $4.6692016091\dots$ | Period-doubling spacing | Governs bifurcation step-down $\mu_k = \mu_\infty + \delta^{-k}$ right across the logistic accumulation point $\mu_\infty \approx 3.56994567$. |
| **Feigenbaum** $|\alpha|$ | $2.5029078751\dots$ | Spatial rescaling | Governs step-size contraction $\alpha^{-k}$ during deterministic intermittent chaotic bursts. |
| **Pinwheel Angle** $\gamma$ | $\arctan(1/2) \approx 26.565^\circ$ | Conway–Radin Pinwheel tiling | Irrational multiple of $\pi$; generates an infinite dense stream of tile orientations that destroys periodic alignment blind spots. |

---

## 2. Next-Generation Tiling Geometries

1. **The Spectre Chiral Monotile (Smith, Myers, Kaplan & Goodman-Strauss, 2023/2024):**
   Unlike the Hat tile (which requires reflected left- and right-handed copies), the Spectre tiles the Euclidean plane aperiodically using **strictly rotations and translations—zero reflections**. In SAT search, eliminating reflections naturally breaks global spin-flip symmetry ($x \leftrightarrow \neg x$), discouraging oscillation between an assignment and its negated twin [audit 8].
2. **Hyperbolic $\{7, 3\}$ Tilings on the Poincaré Disk $\mathbb{D}$:**
   Euclidean tilings grow quadratically in tile count with radius ($O(r^2)$), whereas the Boolean hypercube $\{0,1\}^n$ grows exponentially ($2^n$). A hyperbolic heptagonal tiling has exponentially many tiles at radius $r$, mirroring DPLL/CDCL tree branching without Euclidean overlap.
3. **Wang Tiles (Turing-Complete Seam Matching):**
   Rather than using tiles solely to steer a sweep, finite Wang-tile matching is itself NP-complete (and infinite tiling is undecidable). Unsatisfied clauses compile directly into seam-color mismatches; local inflation/deflation acts as seam-healing rewrite rules.
4. **6D-to-3D Icosahedral Phason Flips:**
   Extending cut-and-project quasicrystals into physical ($\mathbb{R}^3_\parallel$) and perpendicular phason ($\mathbb{R}^3_\perp$) subspaces turns a single variable flip into a **phason flip**—an elastic tile rearrangement that propagates non-local strain across neighboring clauses.

---

## 3. Dual 2D Complex-Plane Transforms ($\mathbb{C}_{\text{in}} \to \mathbb{C}_{\text{out}}$)

Treating $\mathbb{C}_{\text{in}}$ as the tiled configuration space and $\mathbb{C}_{\text{out}}$ as the constraint-interference output plane yields four complementary geometric operators:

* **Möbius Lens Magnification:** $M_{z_0}(z) = e^{i\gamma}\frac{z - z_0}{1 - \bar{z}_0 z}$. Preserves angles and maps the Poincaré disk onto itself while blowing up a tiny high-conflict neighborhood around $z_0$ to fill the entire viewport.
* **$SL(2, \mathbb{Z})$ Modular Quotient Folding:** Mapping $\mathbb{D} \to \mathbb{H}$ via the Cayley transform $\tau = i\frac{1+z}{1-z}$ and reducing $\tau$ into the modular fundamental domain ($|\text{Re}(\tau)| \le \frac{1}{2}$, $|\tau| \ge 1$) models structural quotienting of variable-relabeling redundancies (paired in code with Laplacian Fiedler-vector spectral canonicalization) [audit 5].
* **SAT-Induced Rational/Newton Maps (Fatou & Julia Sets):** Iterating damped Newton steps $z_{t+1} = z_t - \eta F_\Phi(z_t) / F'_\Phi(z_t)$ partitions $\mathbb{C}_{\text{in}}$ into smooth **Fatou convergence basins** separated by a fractal **Julia boundary**. As clause density $\alpha = m/n$ approaches the 3-SAT phase transition ($\alpha_c \approx 4.267$), smooth basins fracture into dense Julia dust [audit 9].
* **Quasiconformal Beltrami Strain ($\partial_{\bar{z}} w = \mu(z)\,\partial_z w$):** Allowing the local clause-conflict gradient to act as a Beltrami coefficient $|\mu(z)| < 1$ warps circles into tilted ellipses; unsatisfied clauses emerge as topological vortices (disclinations) with nonzero winding numbers [audit 2, 9].

---

## 4. The Synthesis Algorithm: HOLOMORPHIC SPECTRE

An all-in attempt to fuse every constant and geometric idea into a single 6-stage solver (implemented in `experiments/04_complex_spectre/holomorphic_spectre_solver.py`):

1. **Stage 1 — Spectral Canonicalization & Hyperbolic Embedding:** Reorder variables using the graph-Laplacian Fiedler vector (plus degree and polarity tie-breaks) to defeat variable-permutation noise, then project onto the Poincaré disk $\mathbb{D}$ with $SL(2,\mathbb{Z})$ folding diagnostics.
2. **Stage 2 — Euler-Phasor Clause Field $F_\Phi(z)$:** Project $z \in \mathbb{C}_{\text{in}}$ into $n$ phasors $\theta_k(z)$ via multi-scale Golden/Plastic/Pinwheel wavevectors, then evaluate:
   $$F_\Phi(z) = \sum_{c=1}^{m} e^{i \cdot 2\pi c / \rho} \prod_{j=1}^{3} \frac{1 - s_{c,j} e^{i\theta_{v(c,j)}(z)}}{2}$$
   Every clause product vanishes identically if and only if at least one literal is satisfied ($s_{c,j} e^{i\theta} = +1$).
3. **Stage 3 — Chiral Spectre Mesh & Cauchy Winding-Number Pruning:** Tile $\mathbb{C}_{\text{in}}$ with rotation-only Spectre tiles and integrate $W(T) = \frac{1}{2\pi i}\oint_{\partial T} \frac{F'_\Phi(z)}{F_\Phi(z)}\,dz$ along each tile perimeter $\partial T$ to screen for interior zeros without enumerating interior states [audit 1–4].
4. **Stage 4 — Möbius Adaptive Zoom & Pinwheel Frame Kicks:** Re-center the Poincaré disk on top-ranked candidate tiles via $M_{z_0}(z)$ and rotate by irrational Pinwheel increments $\gamma = \arctan(1/2)$.
5. **Stage 5 — Continuous Phasor Descent + Feigenbaum $\mu_\infty$ Phason Bursts:** Interleave analytic gradient descent on $E(\theta) = \sum_c \prod_j \frac{1 - s_{c,j}\cos\theta_j}{2}$ with deterministic logistic-map perturbations right above $\mu_\infty \approx 3.56994567$, scaling step sizes by $\alpha^{-k}$ and triggering multi-clause phason flips during extreme intermittent tail excursions.
6. **Stage 6 — Discrete Rounding & Exact Witness Verification:** Snap phasors via $\text{sgn}(\cos \theta_k)$, run conflict-guided polish, and verify $0$ unsatisfied clauses.

---

## 5. Honest Verdict on $P \stackrel{?}{=} NP$: The Complexity Conservation Law

**HOLOMORPHIC SPECTRE does not prove $P = NP$.** More importantly, running it empirically exposes *precisely* how complexity theory defends itself when you try to compress an exponential discrete search space into 2D complex geometry:

### The Phase-Frequency Tax (Nyquist Aliasing Wall)
To make $z \mapsto \theta(z)$ injectively resolve $2^n$ Boolean assignments inside a bounded 2D disk, $F_\Phi(z)$ must oscillate at spatial frequency $\Omega(2^{n/2})$ [audit 3]. Consequently, accurately computing the Cauchy winding integral $\oint_{\partial T} d(\arg F_\Phi)$ without missing $2\pi$ phase wraps requires $\Omega(2^{n/2})$ boundary samples—**exponential work sneaks right back in through Nyquist quadrature!**

* In **Euclidean space**, exponentiality manifests as **state overlap collisions**.
* In **Hyperbolic Poincaré space**, tiles near $|z| \to 1$ shrink as $O(e^{-r})$, forcing **$\Omega(n)$-bit floating-point precision blowup**.
* In **Cauchy contour pruning**, exponentiality manifests as **boundary phase-aliasing**.
* Near **$\alpha_c \approx 4.267$**, replica symmetry breaking (RSB) shatters the continuous energy landscape into exponentially many frozen metastable clusters with $|F_\Phi| \approx 0$ ("ghost basins") [audit 9].

### What Genuinely Worked (Verified Empirically in Experiment 04) [audit 5–7]
Even though worst-case complexity remains exponential, the hybrid mechanisms delivered clear, measurable heuristic gains over fixed-stride deterministic sweeps:

* **Rescuing Deterministic Sweeps on Critical Needles ($n=35, \alpha=4.26$):** Pure Golden Sweep solved **58.3%** (7/12) with 739.0 median steps. **Holomorphic Spectre solved 100% (12/12)** with **185.0 median steps** (4x fewer median steps, also slightly beating Uniform WalkSAT's 239.5 median steps).
* **Spectral Canonicalization Crushes Permutation Drift ($n=30, \alpha=4.1$, 10 permutations):** Without canonicalization, mean steps = **123.7** (std **143.6**). With Fiedler-vector canonicalization, mean steps dropped to **23.3** (std **31.4**)—a **4.6x reduction in variance**.
* **Nyquist Aliasing Measured Directly:** Holding contour sampling fixed at 12 points/edge on a Spectre tile, mean boundary phase variation rose monotonically from **2.72 rad** ($n=6$) to **11.79 rad** ($n=50$), triggering Nyquist aliasing ($|\Delta \arg F_\Phi| > \pi/2$) starting at $n=35$ even under gentle $\sqrt{n}$ frequency scaling.
