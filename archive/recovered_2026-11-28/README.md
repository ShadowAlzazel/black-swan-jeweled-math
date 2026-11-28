# Recovered additions, 2026-11-28

These files came from `black-swan-recovered-additions.zip`, a review of five saved archives against `main` at `626f565`. Everything here is kept byte-for-byte as recovered (commit `b1f3e78`), under its original repository path.

* `REVIEW.md` and `audit/`: the review and its file-level comparison. Every older full-repository ZIP differs from GitHub only in README.md, docs/GITHUB.md, provenance.json and web/index.html, and the current versions are kept.
* `experiments/04_complex_maps/` and `web/05_complex_maps.html`: an older complex-maps lab and notes. Most of it is superseded by `web/06_conformal_jewels.html` and `experiments/04_complex_spectre/feigenbaum_renormalization.py`. Two pieces were integrated:
  * the disk-preserving Möbius map (z − 0.3)/(1 − 0.3z), now an option in `web/06`;
  * the superstable centers in `feigenbaum_centers.json`, computed independently by bisection at 80 digits. They agree with the Newton-method centers of Experiment 04 through period 512 to within 5 × 10⁻³⁰, the precision Experiment 04 stores (**observed**).
* `map_validation.json` records numerical branch, derivative and area checks of the old lab (540 states). The program that generated it was not recovered, so these numbers are archived evidence, not a validation of the current explorer.
* The old HTML lab is kept for reference only; it has not been re-tested in current browsers.

The recovered hybrid experiment itself was integrated as [Experiment 06](../../experiments/06_black_swan_hybrid/README.md).
