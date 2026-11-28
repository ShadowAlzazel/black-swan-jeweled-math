# Black Swan Jeweled Math: recovered additions review

Reviewed GitHub `main` at `626f565364cd3d9a6dd18c08349537f76a130276` against five saved ZIPs and the experiment left in this workspace. The other published branch pointed to the same commit. This package preserves candidate additions; it does not replace the current repository or its newer files.

## Archive comparison — observed

Payloads were compared using Git blob hashes against the complete, non-truncated remote tree (93 files). The five archives were retrieved again and match the earlier local copies by SHA-256. Git metadata and directory entries are excluded from payload counts.

| Archive | Identical files | Files present but changed | Absent files |
| --- | ---: | ---: | ---: |
| Black_Swan_Jewel_Lab.zip | 5 | 0 | 0 |
| Fleurdelys_Tiling_Jewels.zip | 12 | 0 | 0 |
| Golden_Jewel_Sweep.zip | 20 | 1 | 0 |
| black-swan-jeweled-math.zip | 46 | 4 | 0 |
| black-swan-jeweled-math-updated.zip | 46 | 4 | 5 |

The sweep script's change is the existing geometry-directory fallback. The four changed files in each full-repository ZIP are README.md, docs/GITHUB.md, provenance.json, and web/index.html. Keep the current GitHub versions; restoring these older files would discard later work.

## Useful additions

### 1. The hybrid and XOR-kernel experiment — observed

Twelve files under `recovered/experiments/05_black_swan_hybrid/` have not reached GitHub or the five saved ZIPs. The stored run contains 1,584 runs: 72 formulas, 11 configurations, two starts. Its independent audit reports 1,188 checked SAT witnesses, 390 independently checked UNSAT answers, six UNKNOWN answers, and no mismatches. These are existing saved results, inspected during this review; the experiments were not rerun here.

The addition beyond the current counting experiment is recognition of parity equations directly from suitable CNF clause groups, with source-clause tracking, a nullspace basis, and checkable contradiction certificates. Consistent mixed formulas retain their other constraints. The engine combines this exact recognition with bounded geometric guidance and DPLL fallback. It also retains golden/silver phases, hat/pinwheel/quasicrystal geometry, complex maps, Feigenbaum phases, finite zeta/prime and Collatz/Fourier callbacks, adaptive positive weights, and shuffled controls.

| Stored configuration | Determinate answers | DPLL nodes | Front-end assignment evaluations |
| --- | ---: | ---: | ---: |
| Occurrence-only DPLL | 142/144 | 55,522 | 0 |
| XOR recognition then DPLL | 144/144 | 2,326 | 0 |
| Combined guidance and XOR recognition | 144/144 | 1,651 | 15,227 |
| Shuffled combined guidance and XOR recognition | 144/144 | 1,672 | 15,720 |

The denominators are repeated runs, not 144 different formulas. Costs in the last two columns have different units and omit other setup/kernel costs; fewer DPLL nodes do not establish an end-to-end speedup. The close shuffled result does not establish a special constant or tiling advantage. XOR recognition exploits an established polynomial subclass; unrestricted general SAT retains potentially exponential branching. These findings supply no P=NP proof.

Recommended destination: `experiments/06_black_swan_hybrid/`. Before integration, add its own README, update paths in notes/provenance, and cross-reference the existing counting experiment. An interactive cycle-space/kernel lab could later use the next free web number, 08; no such lab is included here.

Two integration checks are still needed. The old geometric schedule generator divides by zero on zero-variable inputs or empty clauses; the saved dataset avoids those inputs, and the exact DPLL component already handles them. Preserve the delivered source before adding guards in a separate fix commit. The newer repository's conditional-expectation starting method should also be added as a control in any next comparison. This old run records total batch time, not per-run timings.

### 2. Complementary complex-map checks — observed

Five files from the updated ZIP are absent from GitHub: four under `experiments/04_complex_maps/` and `web/05_complex_maps.html`. They are preserved unchanged under `recovered/`.

The current `web/06_conformal_jewels.html` already covers the older lab's linked planes, shear, exponential/squaring maps, constants, hats, and pinwheels, and adds cubic, disk, and elliptic maps. The current Feigenbaum renormalization program goes beyond the older period-512 computation. A second copy of that visual lab and a second primary Feigenbaum implementation would largely duplicate existing work.

The pieces worth integrating selectively are:

* the fixed disk-preserving Möbius map `(z - 0.3)/(1 - 0.3*z)` as another option in the current conformal explorer;
* the numerical branch/derivative/area checks recorded in `map_validation.json` (540 states; explicitly not a browser-rendering test);
* the independently bracketed high-precision superstable-center calculation, with residual and proper-period checks, as a cross-check of the current Feigenbaum implementation.

The old validation file describes the old explorer's parameter choices. Adapt and rerun its checks before treating them as validation of the current explorer. The standalone old HTML is retained for reference, with no claim of current-browser validation.

The generator that produced `map_validation.json` was not retained among these files. The saved numbers are archived evidence; making this validation reproducible requires reconstructing that generator.

## Integration order

1. Keep the current GitHub checkout as the base.
2. Import recovered source bytes in a separate commit, retaining original paths in provenance or an archive record.
3. Add the hybrid as Experiment 06 and adapt the complementary complex-map checks/options in a following integration commit.
4. Preserve the current README, landing pages, experiment numbering 01–05, and newer audits. Regenerate documentation/provenance deliberately rather than copying older snapshots over them.

No remote commits or repository edits were made during this review. The `recovered/` folder is a review collection, not an installer: its experiment scripts expect the existing repository's geometry/data directories when integrated.

`audit/archive_vs_remote_audit.json` records file-level comparisons. `audit/manifest.json` records archive and recovered-file hashes. Third-party attribution and the retained hat-source license are included under `reference/`.
