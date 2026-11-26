# Fleurdelys: tiling jewels and adaptive weights

This playground follows the rotations, positivity identities, and spectral ideas in the author's Fleurdelys journey: https://shadowalzazel.github.io/fleurdelys-math/ . Its purpose is exploration. The construction and SAT experiments below do not require RH.

## Reproduce

From the extracted directory:

```sh
python tiling_adaptive.py
node tiling_research/generate_hat_patch.js
python tiling_research/validate_hat_patches.py
python tiling_research/generate_pinwheel.py
```

Only Python's standard library and Node.js are needed. The adaptive script embeds its two selected formulas and a snapshot of actual hat orientations; it can run independently. With the geometry JSON present, it verifies that the source orientations match the snapshot.

## Actual geometric constructions

The hat patches use the authors' H-metatile substitution code from https://github.com/isohedral/hatviz . They contain 4, 25, and 169 hats. The source code, commit hashes, and BSD-3-Clause license are retained in `tiling_research/`. Polygons are not manually packed or decorative imitations. Validation checks equal areas, connected edge-sharing adjacency, and pairwise non-overlap by independent polygon triangulation and intersection.

The pinwheel generator uses an exact five-triangle subdivision of the 1:2:sqrt(5) right triangle. Its levels contain 1, 5, 25, 125, and 625 congruent triangles within each level, with equal areas 5^-level. The number of orientation classes increases with refinement. See `tiling_research/README.md` for the explicit affine construction and checks.

Pinwheel substitution tilings have infinitely many orientations in the infinite construction. The undecorated triangle shape can also tile periodically; the substitution or matching rules matter. The hat itself forces aperiodicity but requires both handednesses. Tile(1,1) forces aperiodicity when reflections are forbidden; the edge-modified Spectre forces chiral aperiodicity even when reflections are allowed. Spectre geometry is discussed, but is not included in the generated patches.

## Positive adaptive weights

Let u_j(x) be 1 when clause j is violated and 0 otherwise. Weighted energy is E_w(x)=sum_j w_j u_j(x). Start every w_j at 1. At a strict local minimum that is not a SAT solution, increase the weights of currently violated clauses.

The orientation-driven 3-4-1 rule uses

g(a)=3+4 cos(a)+cos(2a)=2(1+cos(a))^2 >= 0,

beta_r=2*pi*r/phi^2,

w_j <- w_j + u_j(x) * [1/4 + g(theta_j-beta_r)].

Here theta_j is the orientation of tile j in the actual 25-hat patch, phi is the golden ratio, and r counts weight updates beginning at zero. Clause j is assigned to tile j. Reflection metadata are retained but do not affect this rule.

The bracket is at least 1/4, so weights stay strictly positive. Hence E_w(x)=0 exactly when x satisfies the original formula, at every update. At a weight update the current energy increases; only the intervening bit flips are required to decrease the fixed-weight energy. The changing scores do not form one globally decreasing energy function.

The tile assignment and phase schedule are experimental choices. Spatial proximity of tiles need not match the SAT interaction graph. No advantage from the tiling is assumed or established.

## Selected demonstrations

`tiling_adaptive_results.json` contains full clauses, tile assignments, weights, states, violations, and event histories. Integer bit j represents variable x_(j+1); labels list variables in that order.

| Case | Uniform greedy | Plain Breakout | Hat cosine | Hat 3-4-1 |
| --- | --- | --- | --- | --- |
| 14 clauses, 4 variables | Stalls | 1 update, 2 flips | 1 update, 2 flips | 1 update, 2 flips |
| 17 clauses, 4 variables | Stalls | 9 updates, 9 flips | 8 updates, 11 flips | 6 updates, 11 flips |

All adaptive runs listed reach a satisfying assignment. These deliberately selected examples illustrate the mechanism; their different boost scales prevent interpreting these counts as a fair performance ranking. They are not a representative benchmark.

Checks exhaustively compare all 16 assignments at every event to verify the exact SAT zero set. They also verify one-bit moves, strict descent between updates, boosts restricted to violated clauses, positive weights, final satisfaction, and the trigonometric identity.

## Basel and the integer lattice

Euler's sine product is sin(pi*z)/(pi*z)=product_(n>=1)(1-z^2/n^2). Its z^2 coefficient is -sum_(n>=1)1/n^2. Taylor expansion gives the same coefficient as -pi^2/6, proving zeta(2)=pi^2/6. Local uniform convergence justifies comparing the coefficients.

The density of positive integer pairs (a,b) with gcd(a,b)=1 is product_p(1-p^-2)=1/zeta(2)=6/pi^2. Such points are precisely the lattice points visible from the origin. One can justify the density by first excluding the common prime factors in a finite set using the Chinese remainder theorem, then bounding the remaining common-prime contribution by sum_(p>P)1/p^2. This is an unconditional bridge between the zeta function and planar lattice geometry.

## A next matrix experiment

Let A be the edge-sharing adjacency of a finite tile patch, and P its row-normalized neighbor-average matrix. A possible update is

w_next=(1-lambda)*w + lambda*P*w + eta*u(x),

with 0<=lambda<=1 and eta>=0. Since P has nonnegative entries and every row sums to one, positive weights remain positive. This connects penalty diffusion to the eigenmodes of the tiling graph. It is a proposed next experiment and is not implemented in this bundle.

## References and attribution

- Hat authors and source: https://cs.uwaterloo.ca/~csk/hat/ and https://github.com/isohedral/hatviz
- Spectre authors: https://cs.uwaterloo.ca/~csk/spectre/
- Radin, *The pinwheel tilings of the plane*: https://web.ma.utexas.edu/users/radin/papers/pinwheel.pdf
- Pinwheel construction overview: https://tilings.math.uni-bielefeld.de/substitution/pinwheel/
- NIST sine product: https://dlmf.nist.gov/4.22.E1
- NIST zeta Euler product: https://dlmf.nist.gov/25.2.E11
- Morris, *The Breakout Method for Escaping from Local Minima*: https://cse.unl.edu/~choueiry/Documents/BreakoutMorris1993.pdf

The adaptive schedules and pinwheel implementation are exploratory code produced for this conversation. Upstream hat source code retains its original BSD license and attribution.
