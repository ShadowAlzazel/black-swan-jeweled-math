# Complex tiling maps: two planes, one correspondence

This lab follows a tile from the complex input plane, z=x+iy, to the complex output plane, w=f(z). Shared colors and tile identities connect the two views. The full graph (x,y,Re f(z),Im f(z)) lies in four real dimensions; the paired planes show its correspondence without collapsing the input and output into one picture.

## Constructing the input patch

Parameter coordinates range over

u ∈ [-2 log κ, 2 log κ],  v ∈ [-π, π],

and enter the displayed input plane through the shear

z = u + i(v + τu),  -2 ≤ τ ≤ 2.

Available patterns include rectangular and triangular meshes, clipped hexagonal cells, and actual substitution patches containing 25 hats or 125 pinwheel triangles. The hat and pinwheel patches are uniformly fitted into the parameter rectangle and leave some of it blank. Their output patches consequently need not fill the entire image region.

The shear is an invertible, area-preserving real-linear transformation. For τ≠0 it is not holomorphic. This is separate from the complex map z→f(z), which is holomorphic wherever its formula is analytic.

## Three maps

**Exponential, w=e^z.** Its magnitude and argument are e^u and v+τu. At τ=0, horizontal parameter lines become rays and vertical lines become circles. With shear, a line v=constant becomes a logarithmic spiral. Euler's identity appears as the half-turn e^(z+iπ)=-e^z; a full 2πi shift returns to the same output.

The exponential is injective on the interior of this sheared strip: equal outputs require equal real parts u and imaginary parts differing by 2π, while two interior v values differ by less than 2π. Opposite v boundaries are identified. A mesh filling the parameter rectangle therefore maps to an annulus with radii κ^-2 and κ^2, with its boundary seam identified. The domain choice is what prevents interior overlap; a nonzero derivative alone would not ensure that.

**Squaring, w=z².** Radii square and angles double. Opposite points have the same image, so the full patch generally overlaps in the output. The origin is a critical point, where f′(z)=0. This is a branched covering pattern rather than a globally nonoverlapping tiling.

**Möbius, w=(z-0.3)/(1-0.3z).** This fractional linear map is injective away from its pole at z=10/3. The pole lies outside the selected parameter ranges, including the largest κ. The formula preserves the unit disk, although the entire input patch need not lie inside that disk. It provides a useful contrast with the many-to-one squaring map.

Whole tile edges are densely sampled before mapping. The resulting polylines approximate curved images instead of merely connecting transformed corners with straight chords. Away from critical points and poles, f′ determines local stretch |f′|, rotation arg f′, and area scale |f′|².

## Constants with different roles

κ can be φ, √2, the silver ratio 1+√2, e, |α|≈2.502907875, or δ≈4.669201609. Here κ sets logarithmic input width and exponential radial scale. Selecting a Feigenbaum constant as a scale does not by itself produce period-doubling universality.

In quadratic period-doubling dynamics, δ describes limiting ratios of parameter intervals; α describes spatial rescaling, with a sign that depends on convention. The accompanying finite computation through period 512 gives δ₉≈4.669171555380 and α₉≈-2.502906472219. These are numerical approximations approaching the standard constants, not exact identities.

## Next constructions

- **Elliptic functions:** compare Weierstrass ℘ on Gaussian and Eisenstein period lattices. Periodicity and evenness create repeated output images; poles require explicit treatment.
- **Hyperbolic tilings:** use disk-preserving Möbius maps to move polygon patterns in the Poincaré disk.
- **Renormalization:** compare composition-and-rescaling operators and their eigenmodes. Feigenbaum's δ arises as a relevant eigenvalue, connecting back to the project's spectral experiments.

## Primary references

- [NIST: exponential conformal maps](https://dlmf.nist.gov/4.3)
- [NIST: elliptic lattices and periodicity](https://dlmf.nist.gov/23.2)
- [Molteni: computing Feigenbaum constants and the renormalization operator](https://arxiv.org/html/1602.02357v1)
