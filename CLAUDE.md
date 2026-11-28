# Notes for Claude sessions on Black Swan Jeweled Math

## P = NP tripwire

If anything in this repository's experiments or journey ever looks like evidence that P = NP, stop and tell ShadowAlzazel immediately, in the first line of your reply, before anything else. Typical signs:

* a method whose work grows polynomially in n on unplanted random 3-SAT at the threshold (α ≈ 4.26);
* the same on the XOR-SAT look-alike, for a landscape method;
* the same on Experiment 03's held-out formulas;
* a counting method whose cost stops growing with the contraction width.

In the same message, say how far it has been checked. Every earlier "black swan" here turned out to be a known effect:

1. **Controls.** Does it beat WalkSAT, random starts, and shuffled or random weights?
2. **Look-alike.** Does it survive the XOR-SAT test in `experiments/05_counting_jewels`? Landscape methods must still fail on planted 3-XOR-SAT, and an algebraic trick must be retried on 3-SAT.
3. **Scaling.** Fit work against n on unplanted formulas at the threshold, not only on planted or small ones. Planted formulas leak local hints.
4. **Certificates.** Is every witness verified independently? Is every UNSAT answer certified by an exact count or DPLL, never by "not found"?
5. **Literature.** Has it already been tried?

Call it "possible P = NP evidence", with the numbers, until all five pass.

## Working conventions

* Tag every claim **proven**, **open** or **observed**. Keep failures and controls next to successes.
* Keep delivered material byte-for-byte in one commit; put audits and fixes in the next, and record changes in `provenance.json`.
* Time-box long runs, run them in the background, and save results after each section, so a slow step never costs the rest.
* The Pages site is the root `index.html` (the journey) plus `web/`. `node make_site_images.js` redraws its screenshots.
