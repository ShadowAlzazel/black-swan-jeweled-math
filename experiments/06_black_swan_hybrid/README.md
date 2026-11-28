# Experiment 06 — The black swan hybrid and the parity island

This experiment was recovered on 2026-11-28 from an unfinished local run that never reached GitHub. It is kept byte-for-byte as delivered (commit `b1f3e78`; originally `experiments/05_black_swan_hybrid/`, renumbered because 05 is taken). The review that recovered it is in [archive/recovered_2026-11-28/REVIEW.md](../../archive/recovered_2026-11-28/REVIEW.md).

## What it does

A complete solver with three layers:

1. **Exact parity recognition** (`xor_kernel.py`). Groups of CNF clauses that forbid every assignment of one parity on a set of variables are recognized as XOR equations, with their source clauses tracked. Gaussian elimination over GF(2) then gives the rank, a particular solution, a nullspace basis, or a row combination reading 0 = 1, which is a contradiction certificate anyone can check. Partial groups are never treated as equations, and mixed formulas keep all their other clauses. See [XOR_KERNEL.md](XOR_KERNEL.md).
2. **Creative guidance** (`creative_schedules.py`). Bounded, positive-weight repair guided by every earlier idea: golden and silver phases, hat, pinwheel and quasicrystal geometry, complex maps, Feigenbaum phases, finite zeta and prime features, Collatz and Fourier schedules, and adaptive 3-4-1 boosts. A shuffled control keeps the same weights and graph with their attachment scrambled. See [CALLBACKS.md](CALLBACKS.md).
3. **Complete DPLL fallback** with a node cap. An interrupted search reports UNKNOWN, never UNSAT. See [PROOF_OBLIGATIONS.md](PROOF_OBLIGATIONS.md) for what a real P = NP claim would have to establish.

## Verification (2026-11-28)

* **Reproduction** (**observed**). `run_hybrid.py` regenerates `dataset.json`, `results.json`, `pilot_results.json` and `unresolved_cases.json` identically apart from timing fields, in 32 seconds. `audit_hybrid.py` reproduces `audit_results.json`: 1,584 runs, 1,188 SAT witnesses and 390 UNSAT answers checked independently, 6 UNKNOWN, no mismatches.
* **Kernel self-test** (**observed**). `xor_kernel.py` passes, including all 256 subsets of the eight three-variable clauses, random parity systems, mixed formulas and 80-variable systems.
* **The look-alike, solved** (**observed**). Run on Experiment 05's planted 3-XOR-SAT formulas, the kernel recognizes every equation from its four-clause CNF encoding and returns verified witnesses: 8/8 at every size up to n = 400, about 12 ms per formula there. WalkSAT already failed most of these at n = 100. On random 3-SAT (n = 60, α = 4.26) it returns UNKNOWN for all 20, as it should: there is no parity structure to find.
* **Fix** (integration commit). `make_schedule` divided by zero for zero-variable inputs or empty clauses. It now returns a neutral schedule for those, which exact DPLL already decides. Results on the stored dataset are unchanged.

## Results (stored run, 72 formulas × 2 starts = 144 runs per configuration)

| Configuration | SAT | UNSAT | UNKNOWN | DPLL nodes | Guidance evaluations |
| --- | ---: | ---: | ---: | ---: | ---: |
| occurrence-order DPLL | 108 | 34 | 2 | 55,522 | 0 |
| random-order DPLL | 108 | 36 | 0 | 50,177 | 0 |
| WalkSAT, then DPLL | 108 | 34 | 2 | 54,591 | 20,344 |
| **XOR recognition, then DPLL** | 108 | 36 | 0 | **2,326** | 0 |
| golden / hat | 108 | 36 | 0 | 45,963 | 17,085 |
| silver / pinwheel | 108 | 36 | 0 | 46,513 | 18,392 |
| zeta / quasicrystal | 108 | 36 | 0 | 47,052 | 18,504 |
| Feigenbaum / Möbius | 108 | 36 | 0 | 49,594 | 18,527 |
| Collatz / Fourier | 108 | 34 | 2 | 55,486 | 19,100 |
| all ideas, with XOR recognition | 108 | 36 | 0 | 1,651 | 15,227 |
| shuffled all ideas, with XOR recognition | 108 | 36 | 0 | 1,672 | 15,720 |

The formulas: 24 planted, 24 random, 8 XOR prisms, 6 XOR cubics, 6 pigeonhole, 4 split components.

What the table says (**observed**):

* **The XOR island accounts for nearly all the savings.** Of occurrence-order DPLL's 55,522 nodes, 51,468 are spent on the 12 XOR-cubic runs. XOR recognition removes those entirely and decides the two runs plain DPLL left UNKNOWN.
* **The geometry and constant schedules are equivalent to a random variable order.** They use 46,000–55,000 nodes, against 50,177 for random-order DPLL.
* **"All ideas" matches its shuffled control** (1,651 against 1,672 nodes), so the attachment of constants and tilings carries no measurable information. This agrees with Experiments 03–05.
* **Costs are in different units.** DPLL nodes and guidance evaluations cannot be added, and setup and kernel time are not included. Fewer nodes is not an end-to-end speedup.

Recognized parity systems are an established polynomial subclass (**proven**: Gaussian elimination). General 3-SAT is not a parity system. A single clause x ∨ y ∨ z has 7 satisfying assignments, and an affine set has a power of two, so no linear encoding can replace it (**proven**, [XOR_KERNEL.md](XOR_KERNEL.md)). Branching outside the island remains exponential in the worst case. No P = NP evidence; the tripwire in [CLAUDE.md](../../CLAUDE.md) is not triggered.

## How it connects

* [Experiment 05](../05_counting_jewels/README.md): the XOR look-alike showed landscape methods failing on parity. This kernel is the algebraic side of that test, and together they are the two required controls in CLAUDE.md.
* The contradiction certificates (rows XOR-ing to 0 = 1) are the first-cohomology obstruction in roadmap item R12, now computable for parity systems.
* Next comparison, per the recovery review: add Experiment 05's conditional-expectation starts as a control, and record per-run timings (this run logged only batch totals).

## Run

```sh
python experiments/06_black_swan_hybrid/xor_kernel.py     # kernel self-test, under a second
python experiments/06_black_swan_hybrid/run_hybrid.py     # full run, about 30 seconds; rewrites the result files
python experiments/06_black_swan_hybrid/audit_hybrid.py   # independent audit, about 1 second
```

`run_hybrid.py` overwrites the stored results with this machine's timings; `git checkout` restores the delivered files. The schedules read the hat, pinwheel and quasicrystal data from Experiments 01–03.
