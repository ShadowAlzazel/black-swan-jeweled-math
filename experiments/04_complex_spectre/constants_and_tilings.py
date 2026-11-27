"""
constants_and_tilings.py
------------------------
Core mathematical constants, dual complex-plane transforms (C_in -> C_out),
chiral Spectre / Pinwheel / Hyperbolic tilings, Cauchy winding-number probes,
and Feigenbaum edge-of-chaos phason drivers for Experiment 04 (HOLOMORPHIC SPECTRE).
"""

from __future__ import annotations
import cmath
import math
from dataclasses import dataclass
from typing import List, Sequence, Tuple, Optional
import numpy as np

# ---------------------------------------------------------------------------
# High-precision mathematical constants
# ---------------------------------------------------------------------------
PHI = (1.0 + math.sqrt(5.0)) / 2.0                          # Golden ratio (~1.618033988749895)
SILVER_RATIO = 1.0 + math.sqrt(2.0)                         # Silver ratio (~2.414213562373095)
PLASTIC_RHO = 1.324717957244746025960908854                 # Real root of x^3 - x - 1 = 0 (smallest Pisot number)
PLASTIC_CONJ = complex(-0.6623589786223730, 0.5622795120623013)  # Complex conjugate root of x^3 - x - 1 = 0

FEIGENBAUM_DELTA = 4.669201609102990671853203820466         # Period-doubling spacing ratio
FEIGENBAUM_ALPHA = 2.502907875095892822283902873218         # Spatial rescaling factor (|alpha|)
LOGISTIC_MU_INF = 3.569945671870944901842005151386          # Accumulation point of logistic period doubling
MANDELBROT_C_INF = -1.4011551890920506000                   # Real-axis Feigenbaum point on Mandelbrot set

PINWHEEL_ANGLE = math.atan(0.5)                             # arctan(1/2) rad (~26.565 deg), irrational multiple of pi
GOLDEN_ANGLE = 2.0 * math.pi / (PHI * PHI)                  # ~2.3999632297 rad (~137.508 deg)


# ---------------------------------------------------------------------------
# 3-SAT Representation & DIMACS Parser
# ---------------------------------------------------------------------------
@dataclass
class CNFInstance:
    n_vars: int
    clauses: List[Tuple[int, ...]]
    name: str = "unnamed"
    planted_witness: Optional[List[bool]] = None

    @property
    def m_clauses(self) -> int:
        return len(self.clauses)

    @property
    def alpha(self) -> float:
        return self.m_clauses / max(1, self.n_vars)

    def evaluate(self, assignment: Sequence[bool]) -> int:
        """Return count of UNSATISFIED clauses (0 means SAT)."""
        unsat = 0
        for cl in self.clauses:
            ok = False
            for lit in cl:
                v = abs(lit) - 1
                val = assignment[v]
                if (lit > 0 and val) or (lit < 0 and not val):
                    ok = True
                    break
            if not ok:
                unsat += 1
        return unsat

    def verify_witness(self, assignment: Sequence[bool]) -> bool:
        return len(assignment) == self.n_vars and self.evaluate(assignment) == 0


def load_dimacs_cnf(path: str) -> CNFInstance:
    """Load a standard DIMACS .cnf file (e.g., ../03_golden_sweep/missed_sat_instance_35.cnf)."""
    n_vars = 0
    clauses: List[Tuple[int, ...]] = []
    current: List[int] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("c") or s.startswith("%") or s.startswith("0"):
                continue
            if s.startswith("p"):
                parts = s.split()
                if len(parts) >= 4:
                    n_vars = int(parts[2])
                continue
            for tok in s.split():
                val = int(tok)
                if val == 0:
                    if current:
                        clauses.append(tuple(current))
                        current = []
                else:
                    n_vars = max(n_vars, abs(val))
                    current.append(val)
    if current:
        clauses.append(tuple(current))
    return CNFInstance(n_vars=n_vars, clauses=clauses, name=path)


def check_clause_widths(cnf: CNFInstance) -> None:
    """
    The packed solver arrays hold three literals per clause. Shorter clauses are
    padded by repeating a literal (harmless for OR); longer clauses would be
    silently truncated, so they are rejected instead (audit fix).
    """
    for ci, cl in enumerate(cnf.clauses):
        if not 1 <= len(cl) <= 3:
            raise ValueError(
                f"clause {ci} has {len(cl)} literals; Experiment 04 handles 1-3 literals per clause"
            )


def generate_random_3sat(n_vars: int, alpha: float, seed: int = 0, planted: bool = False) -> CNFInstance:
    """Generate a random 3-SAT instance, optionally with a planted witness."""
    rng = np.random.default_rng(seed)
    m = max(1, int(round(n_vars * alpha)))
    witness = [bool(x) for x in rng.integers(0, 2, size=n_vars)] if planted else None
    clauses: List[Tuple[int, ...]] = []
    while len(clauses) < m:
        vars_pick = rng.choice(n_vars, size=3, replace=False) + 1
        signs = rng.choice([-1, 1], size=3)
        cl = tuple(int(v * s) for v, s in zip(vars_pick, signs))
        if witness is not None:
            # Ensure planted witness satisfies clause
            sat = any((lit > 0 and witness[abs(lit) - 1]) or (lit < 0 and not witness[abs(lit) - 1]) for lit in cl)
            if not sat:
                # Flip one literal sign at random
                idx = int(rng.integers(0, 3))
                cl_list = list(cl)
                cl_list[idx] = -cl_list[idx]
                cl = tuple(cl_list)
        clauses.append(cl)
    tag = f"planted_3sat_n{n_vars}_a{alpha:.2f}_s{seed}" if planted else f"rand_3sat_n{n_vars}_a{alpha:.2f}_s{seed}"
    return CNFInstance(n_vars=n_vars, clauses=clauses, name=tag, planted_witness=witness)


def permute_cnf(cnf: CNFInstance, seed: int = 42) -> Tuple[CNFInstance, List[int]]:
    """Randomly permute variable indices and clause order (for permutation stress audits)."""
    rng = np.random.default_rng(seed)
    perm = rng.permutation(cnf.n_vars).tolist()  # old_idx -> new_idx mapping via inverse
    inv_perm = [0] * cnf.n_vars
    for new_i, old_i in enumerate(perm):
        inv_perm[old_i] = new_i
    new_clauses = []
    for cl in cnf.clauses:
        mapped = tuple(int(math.copysign(inv_perm[abs(l) - 1] + 1, l)) for l in cl)
        new_clauses.append(mapped)
    rng.shuffle(new_clauses)
    new_witness = None
    if cnf.planted_witness is not None:
        new_witness = [cnf.planted_witness[old_i] for old_i in perm]
    return CNFInstance(
        n_vars=cnf.n_vars,
        clauses=new_clauses,
        name=f"{cnf.name}_perm{seed}",
        planted_witness=new_witness,
    ), perm


def spectral_canonical_order(cnf: CNFInstance) -> List[int]:
    """
    Build a permutation-invariant variable ordering using the Fiedler vector
    (2nd smallest eigenvector of the variable-co-occurrence graph Laplacian)
    with degree + polarity-balance tie-breaking.
    """
    n = cnf.n_vars
    adj = np.zeros((n, n), dtype=float)
    pos_deg = np.zeros(n, dtype=float)
    neg_deg = np.zeros(n, dtype=float)
    for cl in cnf.clauses:
        vs = [abs(l) - 1 for l in cl]
        for l in cl:
            if l > 0:
                pos_deg[abs(l) - 1] += 1.0
            else:
                neg_deg[abs(l) - 1] += 1.0
        for i in range(len(vs)):
            for j in range(i + 1, len(vs)):
                adj[vs[i], vs[j]] += 1.0
                adj[vs[j], vs[i]] += 1.0
    deg = adj.sum(axis=1)
    lap = np.diag(deg) - adj
    eigvals, eigvecs = np.linalg.eigh(lap)
    fiedler = eigvecs[:, 1] if n > 1 else np.zeros(n)
    # Fix sign ambiguity of Fiedler vector deterministically via skew / third moment
    skew = float(np.sum(fiedler ** 3))
    if skew < 0:
        fiedler = -fiedler
    keys = [(round(float(deg[i]), 6), round(float(abs(pos_deg[i] - neg_deg[i])), 6), round(float(fiedler[i]), 5), i)
            for i in range(n)]
    keys.sort(key=lambda t: (-t[0], -t[1], t[2]))
    return [k[3] for k in keys]


def canonicalize_cnf(cnf: CNFInstance) -> CNFInstance:
    """
    Audit addition: reorder clauses, and literals inside each clause, by the
    spectral canonical variable ranks. Variable names are unchanged, so a witness
    for the result is a witness for the input. Variable canonicalization alone
    leaves the solver sensitive to clause order (clause weights omega_c and the
    violated-clause picks are indexed by clause position); this removes that.
    """
    order = spectral_canonical_order(cnf) if cnf.n_vars >= 2 else list(range(cnf.n_vars))
    rank = {v: r for r, v in enumerate(order)}

    def lit_key(lit: int) -> Tuple[int, int]:
        return (rank[abs(lit) - 1], 0 if lit > 0 else 1)

    clauses = [tuple(sorted(cl, key=lit_key)) for cl in cnf.clauses]
    clauses.sort(key=lambda cl: [lit_key(l) for l in cl])
    return CNFInstance(n_vars=cnf.n_vars, clauses=clauses, name=f"{cnf.name}_canonical",
                       planted_witness=cnf.planted_witness)


# ---------------------------------------------------------------------------
# Dual Complex-Plane Transforms (C_in -> C_out)
# ---------------------------------------------------------------------------
def mobius_disk_automorphism(z: complex, z0: complex, rot: float = 0.0) -> complex:
    """
    Conformal automorphism of the unit Poincare disk D:
      M(z) = exp(i * rot) * (z - z0) / (1 - conj(z0) * z)
    Magnifies neighborhoods near z0 while preserving angles and circles.
    """
    if abs(z0) >= 0.999999:
        z0 = z0 / abs(z0) * 0.999999
    return cmath.exp(1j * rot) * (z - z0) / (1.0 - z0.conjugate() * z)


def poincare_distance(z1: complex, z2: complex) -> float:
    """Hyperbolic distance on the Poincare unit disk D."""
    r1 = min(abs(z1), 0.999999)
    r2 = min(abs(z2), 0.999999)
    num = 2.0 * abs(z1 - z2) ** 2
    den = max(1e-12, (1.0 - r1 * r1) * (1.0 - r2 * r2))
    return math.acosh(1.0 + num / den)


def reduce_sl2z_fundamental_domain(tau: complex, max_iter: int = 64) -> Tuple[complex, int]:
    """
    Fold a point tau in the upper half-plane H (Im(tau) > 0) into the canonical
    SL(2, Z) modular fundamental domain: |Re(tau)| <= 1/2, |tau| >= 1.
    Returns (reduced_tau, num_steps).
    """
    if tau.imag <= 0:
        tau = complex(tau.real, abs(tau.imag) + 1e-6)
    steps = 0
    for _ in range(max_iter):
        # T: shift Re(tau) into [-0.5, 0.5]
        shift = round(tau.real)
        if shift != 0:
            tau = complex(tau.real - shift, tau.imag)
            steps += 1
        # S: inversion -1 / tau if inside unit circle
        if abs(tau) < 1.0 - 1e-12:
            tau = -1.0 / tau
            steps += 1
        else:
            break
    return tau, steps


def disk_to_upper_half_plane(z: complex) -> complex:
    """Cayley transform mapping unit disk D to upper half-plane H."""
    denom = 1.0 - z
    if abs(denom) < 1e-12:
        denom = 1e-12 + 0j
    return 1j * (1.0 + z) / denom


# ---------------------------------------------------------------------------
# Chiral Spectre-Inspired Tile Patch & Low-Discrepancy Generators
# ---------------------------------------------------------------------------
def spectre_prototile_vertices(scale: float = 1.0, angle: float = 0.0, center: complex = 0j) -> List[complex]:
    """
    The 14-vertex contour used by the delivered Experiment 04 solver.

    AUDIT NOTE (audit_exp04.py): this is *not* the Spectre or Tile(1,1). Thirteen
    unit steps are taken and the polygon is then closed implicitly; the closing
    edge has length sqrt(7) ~= 2.646 at ~10.9 degrees instead of a unit edge at
    270 degrees. The result is a simple (non-self-intersecting) 14-gon. It is kept
    unchanged so the recorded benchmark reproduces exactly; use
    spectre_tile11_vertices() for the true polygonal Tile(1,1). Copies are only
    rotated and translated, never reflected.
    """
    # Delivered headings (multiples of 30 deg). They are not the Tile(1,1) edge sequence
    # (see TILE_AB_EDGES); only the first 13 are used and the polygon closes implicitly.
    headings_deg = [0, 60, 150, 90, 180, 120, 210, 150, 240, 300, 210, 270, 330, 270]
    pts = [0j]
    for h in headings_deg[:-1]:
        rad = math.radians(h)
        pts.append(pts[-1] + cmath.exp(1j * rad))
    # Center and normalize radius
    mean_pt = sum(pts) / len(pts)
    pts = [p - mean_pt for p in pts]
    max_r = max(abs(p) for p in pts) or 1.0
    rot = cmath.exp(1j * angle)
    return [center + scale * (p / max_r) * rot for p in pts]


# Edge sequence of the Smith-Myers-Kaplan-Goodman-Strauss family Tile(a, b):
# (edge type, heading in degrees). Tile(1, sqrt 3) is the hat (area 8*sqrt 3 with
# unit short edges), Tile(sqrt 3, 1) is the turtle, and Tile(1, 1) is the
# equilateral polygon whose curved-edge version is the chiral Spectre.
TILE_AB_EDGES: Tuple[Tuple[str, int], ...] = (
    ("a", 0), ("a", 300), ("b", 30), ("b", 90), ("a", 0), ("a", 60), ("b", 150),
    ("b", 210), ("a", 120), ("a", 180), ("a", 180), ("a", 240), ("b", 330), ("b", 270),
)


def tile_ab_vertices(a: float, b: float) -> List[complex]:
    """Raw vertices of Tile(a, b); the edge sequence closes exactly for all a, b."""
    pts = [0j]
    for kind, heading in TILE_AB_EDGES[:-1]:
        length = a if kind == "a" else b
        pts.append(pts[-1] + length * cmath.exp(1j * math.radians(heading)))
    return pts


def spectre_tile11_vertices(scale: float = 1.0, angle: float = 0.0, center: complex = 0j) -> List[complex]:
    """
    True polygonal Tile(1, 1), centered and normalized exactly like
    spectre_prototile_vertices() so the two contours can be swapped in audits.
    (The Spectre replaces each straight edge by a point-symmetric curve, which
    leaves the area unchanged.)
    """
    pts = tile_ab_vertices(1.0, 1.0)
    mean_pt = sum(pts) / len(pts)
    pts = [p - mean_pt for p in pts]
    max_r = max(abs(p) for p in pts) or 1.0
    rot = cmath.exp(1j * angle)
    return [center + scale * (p / max_r) * rot for p in pts]


def chiral_spectre_mesh_centers(count: int, radius: float = 0.92, twist: float = PINWHEEL_ANGLE) -> List[Tuple[complex, float, float]]:
    """
    Generate a chiral, aperiodic tiling center layout in the Poincare/complex disk
    combining Golden-angle phyllotaxis, Plastic-number radial spacing, and Pinwheel
    irrational orientations (strictly rotation-only, zero reflections).
    Returns list of (center_z, tile_radius, orientation_angle).
    """
    tiles = []
    for k in range(count):
        # Plastic-number modulated radial coordinate for 3-literal symmetry breaking
        u = (k + 0.5) / count
        r = radius * math.sqrt(u)
        theta = k * GOLDEN_ANGLE + (k % 3) * (2.0 * math.pi / (PLASTIC_RHO * 3.0))
        z = r * cmath.exp(1j * theta)
        tile_r = radius / math.sqrt(max(1, count)) * (1.0 + 0.15 * math.sin(k * PLASTIC_RHO))
        orient = (k * twist) % (2.0 * math.pi)  # Pinwheel irrational orientation
        tiles.append((z, tile_r, orient))
    return tiles


# ---------------------------------------------------------------------------
# Euler-Phasor Complex 3-SAT Potential F_Phi(z)
# ---------------------------------------------------------------------------
class EulerPhasorField:
    """
    Lifts a discrete 3-SAT formula Phi into a smooth complex potential
    F_Phi : C_in -> C_out using quasicrystal-projected variable phasors
    e^{i * theta_k(z)} and Plastic-number-spaced complex clause weights.
    """

    def __init__(self, cnf: CNFInstance, freq_scale: float = 1.0, canonicalize: bool = True, seed: int = 0):
        self.cnf = cnf
        self.n = cnf.n_vars
        self.m = cnf.m_clauses
        self.freq_scale = freq_scale

        # Variable ordering (spectral canonicalization counters variable-permutation drift)
        if canonicalize and self.n >= 2:
            self.order = spectral_canonical_order(cnf)
        else:
            self.order = list(range(self.n))
        self.inv_order = [0] * self.n
        for pos, orig in enumerate(self.order):
            self.inv_order[orig] = pos

        # Quasicrystal projection wavevectors k_j in C (using Golden, Plastic, and Pinwheel irrationalities)
        self.k_real = np.zeros(self.n, dtype=float)
        self.k_imag = np.zeros(self.n, dtype=float)
        self.phase_offset = np.zeros(self.n, dtype=float)
        rng = np.random.default_rng(seed)
        jitter = rng.uniform(0.0, 2.0 * math.pi, size=self.n) if seed != 0 else np.zeros(self.n)

        for rank in range(self.n):
            # Aperiodic multi-scale wavevectors inspired by 6D->2D cut-and-project jewels
            mag = freq_scale * (PLASTIC_RHO ** (1.0 + (rank % 7) * 0.35)) * (1.0 + 0.12 * (rank // 7))
            ang = rank * GOLDEN_ANGLE + (rank * PINWHEEL_ANGLE) * 0.5
            self.k_real[rank] = mag * math.cos(ang)
            self.k_imag[rank] = mag * math.sin(ang)
            self.phase_offset[rank] = (rank * math.pi / PHI + jitter[rank]) % (2.0 * math.pi)

        # Plastic-number complex weights omega_c = exp(i * 2*pi*c / rho) to suppress false cancellation
        self.clause_weights = np.array(
            [cmath.exp(1j * 2.0 * math.pi * (c + 1) / PLASTIC_RHO) for c in range(self.m)],
            dtype=np.complex128,
        )

        # Prepack clause indices (in canonical rank space) and literal signs
        check_clause_widths(cnf)
        self.c_vars = np.zeros((self.m, 3), dtype=np.int32)
        self.c_signs = np.zeros((self.m, 3), dtype=np.float64)
        for ci, cl in enumerate(cnf.clauses):
            padded = list(cl) + [cl[-1]] * max(0, 3 - len(cl))
            for j in range(3):
                lit = padded[j]
                self.c_vars[ci, j] = self.inv_order[abs(lit) - 1]
                self.c_signs[ci, j] = 1.0 if lit > 0 else -1.0

    def phases_at(self, z: complex) -> np.ndarray:
        """Compute variable phase vector theta(z) in canonical rank order."""
        return self.k_real * z.real + self.k_imag * z.imag + self.phase_offset

    def assignment_from_phases(self, theta: np.ndarray) -> List[bool]:
        """Round phasors e^{i*theta_k} to Booleans via sgn(cos(theta_k)) in original variable order."""
        bits_canon = np.cos(theta) >= 0.0
        out = [False] * self.n
        for rank, orig in enumerate(self.order):
            out[orig] = bool(bits_canon[rank])
        return out

    def eval_F(self, z: complex) -> complex:
        """
        Evaluate complex clause interference field F_Phi(z) in C_out.
        For a clause (l_1 v l_2 v l_3), literal factor is (1 - s_j * exp(i*theta_j))/2
        so that when s_j * exp(i*theta_j) = +1 (literal satisfied), the factor vanishes!
        Wait: True -> theta=0 -> e^{i0}=+1. If lit is positive (s=+1), satisfied => 1 - (+1)(+1) = 0.
        If lit is negative (s=-1) and var is False (e^{i*pi}=-1), satisfied => 1 - (-1)(-1) = 0.
        Thus every clause product vanishes identically iff at least one literal is satisfied!
        """
        theta = self.phases_at(z)
        phasors = np.exp(1j * theta)
        # Gather literals: shape (m, 3)
        lit_phasors = self.c_signs * phasors[self.c_vars]
        factors = 0.5 * (1.0 - lit_phasors)
        clause_vals = np.prod(factors, axis=1)
        return complex(np.dot(self.clause_weights, clause_vals))

    def eval_energy_and_grad_theta(self, theta: np.ndarray) -> Tuple[float, np.ndarray]:
        """
        Evaluate real nonnegative violation energy E(theta) = sum_c |prod_j (1 - s_j e^{i theta_j})/2|^2
        and its exact analytic gradient dE/dtheta for continuous phase polishing.
        Note: |(1 - s e^{i theta})/2|^2 = (1 - s cos(theta)) / 2 = sin^2((theta - pi*(1-s)/2)/2).
        """
        cos_t = np.cos(theta)
        sin_t = np.sin(theta)
        # per-literal violation v_{c,j} = (1 - s_{c,j} cos(theta)) / 2 in [0, 1]
        v = 0.5 * (1.0 - self.c_signs * cos_t[self.c_vars])
        # clause energy = v0 * v1 * v2
        e_c = v[:, 0] * v[:, 1] * v[:, 2]
        total_e = float(np.sum(e_c))

        # dv/dtheta = +0.5 * s * sin(theta)
        dv = 0.5 * self.c_signs * sin_t[self.c_vars]
        grad = np.zeros(self.n, dtype=float)
        # Product rule across the 3 literals
        np.add.at(grad, self.c_vars[:, 0], dv[:, 0] * v[:, 1] * v[:, 2])
        np.add.at(grad, self.c_vars[:, 1], v[:, 0] * dv[:, 1] * v[:, 2])
        np.add.at(grad, self.c_vars[:, 2], v[:, 0] * v[:, 1] * dv[:, 2])
        return total_e, grad


# ---------------------------------------------------------------------------
# Cauchy Winding-Number Contour Integrator & Nyquist Diagnostic
# ---------------------------------------------------------------------------
@dataclass
class WindingReport:
    winding_number: float
    rounded_winding: int
    total_abs_phase_variation: float
    max_step_phase_jump: float
    nyquist_aliased: bool
    min_abs_F: float
    best_z: complex


def cauchy_winding_number(
    field: EulerPhasorField,
    polygon: Sequence[complex],
    samples_per_edge: int = 16,
) -> WindingReport:
    """
    Approximate the Cauchy Argument Principle contour integral:
      W(T) = (1 / 2*pi) * oint_{dT} d(arg F_Phi(z))
    around polygonal tile boundary `polygon`.

    Also tracks total absolute phase variation and maximum per-step phase jump
    to empirically detect the Nyquist Phase-Aliasing Wall (|Delta arg| > pi/2).
    """
    pts: List[complex] = []
    n_verts = len(polygon)
    for i in range(n_verts):
        a = polygon[i]
        b = polygon[(i + 1) % n_verts]
        for s in range(samples_per_edge):
            t = s / samples_per_edge
            pts.append(a + t * (b - a))
    pts.append(polygon[0])

    vals = [field.eval_F(z) for z in pts]
    min_abs = float("inf")
    best_z = pts[0]
    for z, w in zip(pts, vals):
        aw = abs(w)
        if aw < min_abs:
            min_abs = aw
            best_z = z

    total_dtheta = 0.0
    total_abs_dtheta = 0.0
    max_jump = 0.0
    for i in range(len(vals) - 1):
        w0, w1 = vals[i], vals[i + 1]
        if abs(w0) < 1e-14 or abs(w1) < 1e-14:
            continue
        dtheta = cmath.phase(w1 / w0)
        total_dtheta += dtheta
        ad = abs(dtheta)
        total_abs_dtheta += ad
        if ad > max_jump:
            max_jump = ad

    w_num = total_dtheta / (2.0 * math.pi)
    return WindingReport(
        winding_number=w_num,
        rounded_winding=int(round(w_num)),
        total_abs_phase_variation=total_abs_dtheta,
        max_step_phase_jump=max_jump,
        nyquist_aliased=(max_jump > math.pi / 2.0),
        min_abs_F=min_abs,
        best_z=best_z,
    )


# ---------------------------------------------------------------------------
# Feigenbaum Edge-of-Chaos Phason Driver
# ---------------------------------------------------------------------------
class FeigenbaumPhasonDriver:
    """
    Deterministic edge-of-chaos perturbation engine tuned right at/above
    the Feigenbaum accumulation point mu_inf (~3.56994567), rescaling
    bifurcation offset by delta^{-k} and spatial kick radius by alpha^{-k}.
    Produces heavy-tailed escape bursts without a PRNG.
    """

    def __init__(self, r0: float = 1.0 / PHI, level_k: int = 1):
        self.r = r0 % 1.0
        if self.r <= 0.05 or self.r >= 0.95:
            self.r = 0.61803398875
        self.k = level_k
        self.step_count = 0

    @property
    def mu(self) -> float:
        # Slightly above accumulation point so intermittent chaotic bursts occur
        return min(3.99, LOGISTIC_MU_INF + (FEIGENBAUM_DELTA ** (-self.k)) * 0.85)

    def next_kick(self) -> Tuple[float, complex, bool]:
        """
        Advance logistic map and return (scalar_r, complex_kick, is_burst).
        """
        self.step_count += 1
        mu = self.mu
        self.r = mu * self.r * (1.0 - self.r)
        # Detect intermittent burst when trajectory visits extreme tails
        is_burst = (self.r < 0.22) or (self.r > 0.86)
        amp = (FEIGENBAUM_ALPHA ** (-max(0, self.k - 1))) * (2.2 if is_burst else 0.45) * self.r
        angle = self.step_count * PINWHEEL_ANGLE + (self.step_count % 3) * (2.0 * math.pi / PLASTIC_RHO)
        return self.r, amp * cmath.exp(1j * angle), is_burst

    def adapt_level(self, stagnant: bool) -> None:
        """Zoom out toward chaos on stagnation, zoom in along alpha^-k when descending."""
        if stagnant:
            self.k = max(0, self.k - 1)
        else:
            self.k = min(6, self.k + 1)


class PRNGPhasonDriver(FeigenbaumPhasonDriver):
    """
    Audit control: identical burst rule, amplitudes, angles, and level schedule,
    but r is drawn from a seeded uniform PRNG instead of the logistic map. If the
    solver behaves the same, the edge-of-chaos property is not what helps.
    """

    def __init__(self, seed: int = 0, level_k: int = 1):
        super().__init__(level_k=level_k)
        self._rng = np.random.default_rng(seed)

    def next_kick(self) -> Tuple[float, complex, bool]:
        self.step_count += 1
        self.r = float(self._rng.random())
        is_burst = (self.r < 0.22) or (self.r > 0.86)
        amp = (FEIGENBAUM_ALPHA ** (-max(0, self.k - 1))) * (2.2 if is_burst else 0.45) * self.r
        angle = self.step_count * PINWHEEL_ANGLE + (self.step_count % 3) * (2.0 * math.pi / PLASTIC_RHO)
        return self.r, amp * cmath.exp(1j * angle), is_burst
