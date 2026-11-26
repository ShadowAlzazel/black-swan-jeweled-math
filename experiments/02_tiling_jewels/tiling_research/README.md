# Hat and pinwheel construction data

## Hat

`geometry.js` and `hat.js` are unmodified source files from the primary authors' `isohedral/hatviz` repository. The BSD 3-Clause license is retained as `LICENSE` and must accompany redistribution of the code or adaptations.

Sources:

- https://github.com/isohedral/hatviz
- https://cs.uwaterloo.ca/~csk/hat/

File revisions: geometry.js `46e19f05545ed0152776c5fbeec74547932273f1`; hat.js `d8c1839bd4861c6f871eba83d6a42d010f1d9551`.

`generate_hat_patch.js` loads the authors' geometry and substitution code into a Node context and traverses actual H metatiles to obtain 4, 25, and 169 hats. It does not position copies heuristically. `validate_hat_patches.py` triangulates nonconvex polygons, checks equal areas and pairwise nonoverlap, and extracts adjacency from positive shared boundary length. Regenerate with:

```sh
node generate_hat_patch.js
python validate_hat_patches.py
```

Every tile includes its hierarchical address, full affine transform, polygon vertices, orientation in degrees, and a reflection flag. Orientation is the direction of the first affine column; all angles are multiples of 60 degrees. Reflection is the sign of the affine determinant. The 169-tile patch contains 22 reflected hats. This is precisely why rotation and reflection must not be conflated.

## Pinwheel

`generate_pinwheel.py` independently derives a five-triangle pinwheel-style substitution. Its mathematical basis is the Conway--Radin construction:

- https://web.ma.utexas.edu/users/radin/papers/pinwheel.pdf
- https://chaimgoodmanstrauss.com/the-pinwheel-tiling/

For the triangle T=[(0,0),(2,0),(0,1)], apply A=[[2,-1],[1,2]], a rotation and scale by sqrt(5). The resulting triangle has vertices (0,0),(4,2),(-1,2). The script gives an explicit dissection into five congruent triangles with sides 1,2,sqrt(5). Applying A^-1 gives five child affine maps inside T. Iterating those maps yields 1,5,25,125,625 triangles of equal area at each refinement. Each refinement fills the same parent triangle.

This explicit dissection and implementation are independently derived, rather than copied from the authors' software. The script checks the initial dissection, every tile's expected area, total area, and all pairwise overlaps through level 3. Level 4 inherits nonoverlap from recursive affine copies of the validated partition.

```sh
python generate_pinwheel.py
```

Every triangle is ordered right-angle vertex, long-leg endpoint, short-leg endpoint. `orientationDegrees` is the long-leg direction measured by the first affine column; `reflected` is the determinant sign. Retaining reflected child types is essential to the generation of increasing orientation classes under refinement. The undecorated right triangle alone also admits periodic tilings: the substitution hierarchy (or suitable matching rules) specifies the pinwheel family.

## Relation to the hat / Spectre distinction

The primary authors distinguish these carefully:

- Hat: every tiling requires left- and right-handed copies.
- Equilateral Tile(1,1): aperiodic when reflections are forbidden, but periodic tilings are possible if reflections are allowed.
- Spectre: modifies Tile(1,1)'s edges to force chiral aperiodic tilings by translations and rotations, even when reflections are permitted.

Source: https://cs.uwaterloo.ca/~csk/spectre/

Neither construction alone proves anything about SAT complexity. One may use tile positions, orientation, chirality, or hierarchy to propose clause orderings and positive weight schedules, while verifying the original Boolean clauses exactly.
