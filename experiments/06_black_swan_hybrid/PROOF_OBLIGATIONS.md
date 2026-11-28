# What a Black Swan algorithm would have to establish

The experiment has three different kinds of mathematical result. They must remain separate.

1. **Representation identities.** Positive clause penalties preserve the satisfying assignments. Bijective variable relabelings and complemented variables preserve satisfiability. A recognized XOR constraint can be represented as a linear equation over GF(2).
2. **Finite numerical findings.** A bounded search can find solutions, encounter traps, or exhaust a budget on the formulas tested. These observations measure that implementation on those inputs.
3. **An algorithm for every input.** Establishing P=NP by an algorithm requires a correct decision procedure for all instances of an NP-complete problem, together with a worst-case polynomial bound in the encoded input length.

## Correct answers include the negative instances

A satisfying assignment is directly checkable. Failing to find one is not an UNSAT certificate. The hybrid must report UNKNOWN if a search cap interrupts the exact part. If it reports UNSAT, the complete case analysis or an independently checkable contradiction must justify that answer.

The complete DPLL fallback can establish correctness when allowed to finish. Its unrestricted search can still branch into both values of each Boolean variable. A generic recurrence of the form T(n) <= 2 T(n-1) + poly(input size) gives an exponential upper bound, not a polynomial one. A finite cap bounds the runtime of an incomplete experiment; it does not turn UNKNOWN into a correct decision.

## Geometry must retain logical information

The graph of a complex function has input and output coordinates, but a drawing may merge distinct inputs. In particular, z² maps z and -z to the same output, and exp(z+2πi)=exp(z). Neither operation is a lossless encoding of arbitrary Boolean states by itself.

All heuristic coordinates must therefore retain the Boolean variable or state labels. Coordinates may guide priorities, penalties, and proposals. They must not identify distinct logical states or discard assignments without an additional proof. A nonzero complex derivative establishes local conformality, not global injectivity or efficient inversion.

Cut-and-project windows select a quasicrystal point set. Applying such a window as a filter on arbitrary SAT assignments could delete the only solution. The window is a source of coordinates in these experiments, not an unproved satisfiability-preserving pruning rule.

## Positivity has a precise, limited consequence

For u_j(x) in {0,1} recording clause violations, E_w(x)=sum_j w_j u_j(x) has the exact SAT zero set whenever every w_j>0. Positive adaptive boosts and nonnegative stochastic averaging can preserve that invariant.

The invariant does not imply that a local search reaches zero, that there are no nonzero local minima, or that only polynomially many weight changes suffice. Weighted energy can rise when the weights change. Claims of descent apply only during moves made under a fixed set of weights.

## The polynomial XOR island

If a formula is exactly a conjunction of recognized parity equations, Gaussian elimination solves it in polynomial time. Rank, a particular solution, and a nullspace basis describe its affine solution set. A row combination yielding 0=1 certifies a contradiction.

This is an established tractable subclass. A partial parity clause group is not automatically an XOR equation. Ordinary three-literal OR clauses cannot all be replaced by linear parity equations while preserving their truth tables. Mixed formulas retain every constraint that has not been accounted for by a verified equivalence.

For example, one OR clause on three distinct variables has seven satisfying assignments. A consistent affine solution set over GF(2) has 2^d points. Projection of an affine set, including a system with linear auxiliary variables, is still affine. Thus a single affine system with existential linear auxiliaries cannot preserve exactly those seven assignments under ordinary coordinate projection. This statement concerns that particular representation; it does not rule out arbitrary nonlinear decision reductions, unions of affine pieces, or general SAT algorithms.

For parity constraints on the edges of a connected graph, the GF(2) incidence matrix has rank V-1. A consistent edge assignment therefore has a cycle-space nullity E-V+1. An odd total vertex charge gives a row-sum contradiction; an even total charge is consistent. This is an exact connection between cycles, kernels, and a tractable family, rather than an assumption that generic SAT clauses are linear.

## Constants and spectra do not provide an unproved search bound

Golden-angle rotations, Feigenbaum orbit phases, finite Collatz parity words, and zeta-zero frequencies are deterministic guidance choices. Random phases and shuffled controls help assess whether their structure contributes beyond generic perturbations.

The finite Collatz callback has an explicit step limit and assumes no Collatz conjecture. Finite computed critical-line zeta zeros can be used without assuming RH globally. Assuming RH supplies no automatic polynomial SAT bound. Feigenbaum's parameter and spatial scaling constants do not make arbitrary search dynamics universal merely by appearing in a coefficient.

Similarly, a useful spectral layout or a large graph spectral gap is not by itself a proof that all global Boolean constraints are simultaneously satisfiable or easy to search.

## Randomness and the meaning of P

Seeded random restarts and randomized priorities are experimental tools. Good success probabilities, expected runtimes, or high-probability runtime bounds are not the same as a deterministic worst-case polynomial algorithm. A P=NP argument based on randomized search would additionally need a valid deterministic construction or derandomization argument. Fixing a seed does not prove that the selected deterministic run succeeds on every input.

## Evidence that would trigger immediate re-examination

A prospective breakthrough needs an explicit algorithm, an invariant proving both SAT and UNSAT answers correct, and a worst-case operation/bit-complexity argument with every transformation and decoding step included. That argument must survive adversarial inputs and equivalent encodings, not only a favorable benchmark.

No finite success percentage, smooth picture, approximate zero, or easy family of inputs supplies those missing universal quantifiers. Restricted polynomial cases and useful heuristics remain worthwhile findings in their own right.
