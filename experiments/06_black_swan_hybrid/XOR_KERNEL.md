# The parity kernel: an exact special case

This component recognizes Boolean parity constraints inside a CNF formula. It
connects the cycle/kernel/Fourier ideas to an existing polynomial-time tool:
Gaussian elimination over the two-element field GF(2).

For a variable set S, the parity character is

chi_S(x) = (-1)^(sum_(i in S) x_i).

An XOR equation fixes that character to +1 or -1. A collection of such equations
is a linear system A*x=b over GF(2). Its satisfying assignments form either an
empty set or an affine translate of the nullspace of A.

## Exact CNF recognition

A clause with distinct variables is false at exactly one assignment to those
variables. For a three-variable XOR equation, all four assignments of the
opposite parity must be forbidden. Their four clauses encode exactly that XOR
equation. Width-one and width-two versions are handled as well.

The detector groups clauses by variable set and checks whether every forbidden
pattern of a given parity is present. Duplicate clauses are harmless. A partial
group is retained and is never treated as a complete XOR equation. If both
parities are fully forbidden, the two recognized equations are contradictory.
Tautological clauses are ignored; an empty clause is an immediate contradiction.
Clauses wider than three are retained without exponential pattern enumeration.

## Result semantics

`analyze_cnf(n, clauses)` returns:

- `SAT` only when the whole non-tautological CNF is covered by recognized XOR
  equations and the returned integer assignment satisfies the original CNF.
- `UNSAT` when Gaussian elimination derives a contradiction from recognized
  equations, including when other CNF clauses remain, or an empty clause exists.
- `UNKNOWN` for a consistent mixed or partial formula. The complete solver must
  retain the entire original CNF and continue its ordinary search.

An integer state uses bit j for variable x_(j+1). Equations are dictionaries with
`mask`, `rhs`, and original source `clauses` indices. The result also supplies the
rank, nullity, a particular solution when consistent, and a basis of integer
nullspace vectors.

For an inconsistent system, `contradiction_combination` selects original XOR
equations whose row-wise XOR has zero variable mask and right-hand side one:
0=1. This certificate can be checked without trusting the elimination process.
Each equation's source clauses can separately be truth-table checked.

## Scope and cost

For width at most three, recognizing each small group uses at most eight local
patterns. Gaussian elimination has a conventional cubic polynomial bound in the
system size (the implementation uses packed Python integers). This solves the
recognized parity subclass efficiently. General OR-3SAT is not a parity system,
so this does not give a polynomial-time algorithm for all 3SAT instances.

There is a useful obstruction to a direct solution-set-preserving linearization:
one clause x OR y OR z has seven satisfying assignments, while a nonempty affine
subspace of GF(2)^3 has 1, 2, 4, or 8 points. Adding existentially quantified
linear auxiliary variables does not fix this: a projection of an affine space
is still affine. Thus exact OR encoding needs some nonlinearity, disjunction,
or branching beyond a single linear system. This cardinality argument concerns
preserving the assignment set; it does not rule out arbitrary decision-problem
reductions.

In a connected graph's Tseitin parity system, edge variables meet vertex parity
equations. XORing every vertex equation cancels each edge twice. An odd total
charge therefore gives the kernel an explicit contradiction; consistent charges
leave an affine cycle space. That is a direct link between graph cycles,
nullspaces, and a family that a generic branching search can find difficult.

For example, the complete graph on four vertices has six edge variables and
four vertex parity equations. With all charges zero, the matrix has rank three
and nullity three: its eight solutions are a cycle space. Flipping one vertex
charge makes the XOR of all four equations read 0=1, an immediately checkable
UNSAT certificate.

## Verification

Run `python xor_kernel.py`. The self-test independently enumerates all 256
subsets of the eight full-width three-variable clauses, 800 random parity
systems through eight variables, 600 mixed formulas, edge cases, and consistent
and inconsistent 80-variable systems. It checks
SAT witnesses, contradiction row combinations, nullspace vectors, and equality
between the generated affine family and every satisfying assignment. These are
finite correctness checks in addition to the algebraic justification above.
