"""
gelfand_graev_check.py -- numerical check of Lemma 14.33 and Proposition 14.34
(Gelfand-Graev weights) on GL_2(F_q) for q = 2, 3, 5, 7, plus a counting check of
the weight sum on GL_N(F_q) for N <= 4.

    python docs/gelfand_graev_check.py

Notation follows the shared text: Gamma_N = GL_N(F_q) (the text writes g for q),
U = upper unitriangular matrices, psi_0(x) = exp(2 pi i x / q) (q prime here),
e_U = |U|^-1 sum_u conj(psi_0(u_12)) u, and for a central idempotent E_c
    wt(e, c) = (e E_c)(1) / e(1).
Labels c are characteristic polynomials of the semisimple parameter s of each
Deligne-Lusztig (Lusztig) series. Standard library plus numpy only.
"""

from __future__ import annotations

import cmath
import itertools
import math
import time
from fractions import Fraction

import numpy as np

QS = (2, 3, 5, 7)


# ---------------------------------------------------------------------------
# GL_2(F_q) as tuples (a, b, c, d) = [[a, b], [c, d]]
# ---------------------------------------------------------------------------
def mul(x, y, q):
    a1, b1, c1, d1 = x
    a2, b2, c2, d2 = y
    return ((a1 * a2 + b1 * c2) % q, (a1 * b2 + b1 * d2) % q,
            (c1 * a2 + d1 * c2) % q, (c1 * b2 + d1 * d2) % q)


def inv(x, q):
    a, b, c, d = x
    det = (a * d - b * c) % q
    r = pow(det, q - 2, q)
    return ((d * r) % q, (-b * r) % q, (-c * r) % q, (a * r) % q)


def group(q):
    return [m for m in itertools.product(range(q), repeat=4) if (m[0] * m[3] - m[1] * m[2]) % q]


# ---------------------------------------------------------------------------
# F_{q^2} = F_q[x] / (x^2 - t0 x + d0), elements (u, v) = u + v x
# ---------------------------------------------------------------------------
class Fq2:
    def __init__(self, q):
        self.q = q
        self.t0, self.d0 = next((t, d) for t in range(q) for d in range(1, q)
                                if all((x * x - t * x + d) % q for x in range(q)))
        self.elems = [(u, v) for u in range(q) for v in range(q) if (u, v) != (0, 0)]
        self.gamma = next(g for g in self.elems if self.order(g) == q * q - 1)
        self.log = {}
        x = (1, 0)
        for j in range(q * q - 1):
            self.log[x] = j
            x = self.mul(x, self.gamma)
        g0 = self.pow(self.gamma, q + 1)
        assert g0[1] == 0, "the norm of a generator must lie in F_q"
        self.g0 = g0[0]

    def mul(self, a, b):
        q, t0, d0 = self.q, self.t0, self.d0
        u1, v1 = a
        u2, v2 = b
        return ((u1 * u2 - d0 * v1 * v2) % q, (u1 * v2 + u2 * v1 + t0 * v1 * v2) % q)

    def pow(self, a, e):
        r = (1, 0)
        for _ in range(e):
            r = self.mul(r, a)
        return r

    def order(self, a):
        x, k = a, 1
        while x != (1, 0):
            x = self.mul(x, a)
            k += 1
        return k

    def log_fq(self, a):
        """Discrete log of a in F_q^x to the base g0 = N(gamma)."""
        j = self.log[(a % self.q, 0)]
        assert j % (self.q + 1) == 0
        return j // (self.q + 1)


def zeta(n, k):
    return cmath.exp(2j * math.pi * (k % n) / n) if n > 1 else 1.0


# ---------------------------------------------------------------------------
# Classify elements and build the classical character table of GL_2(F_q)
# ---------------------------------------------------------------------------
def classify(g, q, F):
    a, b, c, d = g
    t, D = (a + d) % q, (a * d - b * c) % q
    roots = sorted({r for r in range(q) if (r * r - t * r + D) % q == 0})
    if len(roots) == 2:
        return ("split", tuple(roots), D)
    if len(roots) == 1:
        r = roots[0]
        return ("central" if g == (r, 0, 0, r) else "nonss", r, D)
    xi = next(x for x in F.elems
              if F.mul(x, x) == ((t * x[0] - D) % q, (t * x[1]) % q))     # x^2 = t x - D
    return ("elliptic", xi, D)


def character_table(q, F):
    """List of irreducibles: kind, parameters, degree, label (c1, c0) of X^2 + c1 X + c0, value function."""
    n1, n2 = q - 1, q * q - 1
    al = lambda k, a: zeta(n1, k * F.log_fq(a))
    ph = lambda k, x: zeta(n2, k * F.log[x])
    elt = lambda a: F.pow(F.gamma, 0) if a % q == 1 else (a % q, 0)
    irr = []
    for k in range(n1):
        a = F.pow(F.gamma, (q + 1) * k)[0]                      # a_k = g0^k
        lab = ((-2 * a) % q, (a * a) % q)
        irr.append(dict(kind="U", par=k, dim=1, label=lab,
                        f=lambda cl, k=k: al(k, cl[2])))
        def fV(cl, k=k):
            kind, data, D = cl
            if kind == "central":
                return q * al(k, data) ** 2
            if kind == "nonss":
                return 0.0
            if kind == "split":
                return al(k, D)
            return -al(k, D)
        irr.append(dict(kind="V", par=k, dim=q, label=lab, f=fV))
    for k, l in itertools.combinations(range(n1), 2):
        a, b = F.pow(F.gamma, (q + 1) * k)[0], F.pow(F.gamma, (q + 1) * l)[0]
        lab = ((-(a + b)) % q, (a * b) % q)
        def fW(cl, k=k, l=l):
            kind, data, D = cl
            if kind == "central":
                return (q + 1) * al(k, data) * al(l, data)
            if kind == "nonss":
                return al(k, data) * al(l, data)
            if kind == "split":
                r1, r2 = data
                return al(k, r1) * al(l, r2) + al(k, r2) * al(l, r1)
            return 0.0
        irr.append(dict(kind="W", par=(k, l), dim=q + 1, label=lab, f=fW))
    seen = set()
    for k in range(n2):
        if k % (q + 1) == 0 or k in seen:
            continue
        kq = (q * k) % n2
        seen |= {k, kq}
        xi = F.pow(F.gamma, k)
        xiq = F.pow(xi, q)
        s, nrm = F.mul(xi, (1, 0)), F.mul(xi, xiq)
        tr = ((xi[0] + xiq[0]) % q, (xi[1] + xiq[1]) % q)
        assert tr[1] == 0 and nrm[1] == 0
        lab = ((-tr[0]) % q, nrm[0])
        def fX(cl, k=k):
            kind, data, D = cl
            if kind == "central":
                return (q - 1) * ph(k, (data, 0))
            if kind == "nonss":
                return -ph(k, (data, 0))
            if kind == "split":
                return 0.0
            return -(ph(k, data) + ph(k, F.pow(data, q)))
        irr.append(dict(kind="X", par=k, dim=q - 1, label=lab, f=fX))
    return irr


def label_type(lab, q):
    """Factorization type of X^2 + c1 X + c0 over F_q, found independently by root search."""
    c1, c0 = lab
    roots = [r for r in range(q) if (r * r + c1 * r + c0) % q == 0]
    if len(roots) == 2:
        return "split", (1, 1)
    if len(roots) == 1:
        return "double", None
    return "irreducible", (2,)


# ---------------------------------------------------------------------------
# GL_2 checks
# ---------------------------------------------------------------------------
def check_gl2(q):
    t0 = time.time()
    F = Fq2(q)
    G = group(q)
    order = len(G)
    index = {g: i for i, g in enumerate(G)}
    U = [(1, x, 0, 1) for x in range(q)]
    psi = lambda u: zeta(q, u[1])
    cls = [classify(g, q, F) for g in G]
    irr = character_table(q, F)
    M = np.array([[ch["f"](c) for c in cls] for ch in irr], dtype=complex)
    dims = np.array([ch["dim"] for ch in irr])
    out = {"q": q, "|G|": order, "irreducibles": len(irr)}
    out["sum of squared degrees = |G|"] = int((dims ** 2).sum()) == order
    gram = M @ M.conj().T / order
    out["max |Gram - I|"] = float(np.abs(gram - np.eye(len(irr))).max())
    out["degrees match values at 1"] = bool(np.allclose(M[:, index[(1, 0, 0, 1)]], dims))

    # Gelfand-Graev character from the definition of induction, one evaluation per class
    key = lambda g: ((g[0] + g[3]) % q, (g[0] * g[3] - g[1] * g[2]) % q, g[1] == 0 and g[2] == 0 and g[0] == g[3])
    inverses = [inv(x, q) for x in G]
    phi_by_key = {}
    for g in G:
        kk = key(g)
        if kk in phi_by_key:
            continue
        total = 0j
        for x, xi in zip(G, inverses):
            y = mul(mul(x, g, q), xi, q)
            if y[0] == 1 and y[3] == 1 and y[2] == 0:
                total += psi(y)
        phi_by_key[kk] = total / len(U)
    Phi = np.array([phi_by_key[key(g)] for g in G])
    out["classes"] = len(phi_by_key)
    out["Phi(1) = [G:U]"] = abs(Phi[index[(1, 0, 0, 1)]] - order / q) < 1e-9
    mult = (M.conj() @ Phi / order).real
    mint = np.rint(mult).astype(int)
    out["multiplicities are integers"] = float(np.abs(mult - mint).max()) < 1e-9
    out["multiplicity one"] = set(mint.tolist()) <= {0, 1}
    out["Phi = sum of its constituents"] = bool(np.allclose(mint @ M, Phi))
    constituents = {ch["kind"] for ch, m in zip(irr, mint) if m}
    kinds_all = {ch["kind"] for ch in irr}
    out["constituent kinds"] = "".join(sorted(constituents))
    out["constituents = the generic (non-1-dimensional) irreducibles"] = constituents == kinds_all - {"U"}

    # central characters agree with the label's constant term under the fixed duality
    ok = True
    for ch, row in zip(irr, M):
        c0 = ch["label"][1]
        for a in range(1, q):
            val = row[index[(a, 0, 0, a)]] / ch["dim"]
            ok &= abs(val - zeta(q - 1, F.log_fq(c0) * F.log_fq(a))) < 1e-9
    out["central characters match labels"] = ok

    # weights, three ways
    GU = Fraction(order, q)                                     # [G:U]
    labels = sorted({ch["label"] for ch in irr})
    wt_exact, details = {}, []
    eU_one = 1 / len(U)
    perms = []                                                 # h -> index of u^{-1} h
    for u in U:
        ui = inv(u, q)
        perms.append(np.array([index[mul(ui, h, q)] for h in G]))
    eU = np.array([psi(u).conjugate() / len(U) for u in U])
    max_dev_norm = max_dev_direct = 0.0
    for lab in labels:
        members = [i for i, ch in enumerate(irr) if ch["label"] == lab]
        wt_exact[lab] = sum(Fraction(int(mint[i]) * irr[i]["dim"]) for i in members) / GU
        E = sum((irr[i]["dim"] / order) * M[i].conj() for i in members)       # E_c as a function on G
        # (e_U E_c)(1) = sum_u e_U(u) E_c(u^{-1})
        coef = sum(eU[j] * E[index[inv(u, q)]] for j, u in enumerate(U))
        conv = sum(eU[j] * E[perms[j]] for j in range(len(U)))                 # (e_U * E_c)(h)
        norm2 = float(np.vdot(conv, conv).real)
        max_dev_norm = max(max_dev_norm, abs(coef - norm2))
        max_dev_direct = max(max_dev_direct, abs(coef / eU_one - float(wt_exact[lab])))
        kind, degs = label_type(lab, q)
        gg = [irr[i] for i in members if mint[i]]
        details.append((lab, kind, degs, len(members), [c["kind"] + str(c["dim"]) for c in gg], wt_exact[lab]))
    out["(e_U E_c)(1) = ||lambda(e_U) lambda(E_c) delta_1||^2 (max dev)"] = max_dev_norm
    out["direct wt = Frobenius wt (max dev)"] = max_dev_direct
    out["Lemma 14.33: all weights >= 0"] = min(wt_exact.values()) >= 0
    out["Lemma 14.33: weights sum to 1"] = sum(wt_exact.values()) == 1
    out["e_U self-adjoint idempotent"] = bool(
        all(abs(eU[j] - eU[[U.index(inv(u, q)) for u in U][j]].conjugate()) < 1e-12 for j in range(len(U)))
        and all(abs(sum(eU[i] * eU[U.index(mul(inv(U[i], q), u, q))] for i in range(len(U))) - eU[k]) < 1e-12
                for k, u in enumerate(U)))

    # Proposition 14.34 on squarefree labels
    prop_ok, one_constituent = True, True
    for lab, kind, degs, n_members, gg, w in details:
        one_constituent &= len(gg) == 1
        if degs is not None:
            prop_ok &= w == Fraction(1, math.prod(q ** d - 1 for d in degs))
    out["each label has exactly one GG constituent"] = one_constituent
    out["Prop 14.34: wt = 1/prod(q^d_i - 1) on squarefree labels"] = prop_ok

    # the proof's degree claim chi_{s,reg}(1) = [G:C(s)]_{p'}, and the corrected version
    Gp, Gpp = q, (q - 1) ** 2 * (q + 1)                         # |G|_p, |G|_p'
    C = {"double": (Gp, Gpp), "split": (1, (q - 1) ** 2), "irreducible": (1, q * q - 1)}
    claim, corrected = {}, {}
    uncorrected_sum = Fraction(0)
    for lab, kind, degs, n_members, gg, w in details:
        deg = int(gg[0][1:])
        Cp, Cpp = C[kind]
        claim[kind] = claim.get(kind, True) and deg == Gpp // Cpp
        corrected[kind] = corrected.get(kind, True) and deg == (Gpp // Cpp) * Cp
        uncorrected_sum += Fraction(Gpp // Cpp) / GU
    out["claim chi_reg(1) = [G:C(s)]_p' holds, by label type"] = claim
    out["corrected chi_reg(1) = [G:C(s)]_p' |C(s)|_p holds, by label type"] = corrected
    out["sum of weights if the uncorrected degree were used"] = f"{uncorrected_sum} (q/(q+1) = {Fraction(q, q + 1)})"
    out["seconds"] = round(time.time() - t0, 2)
    return out, details, wt_exact


# ---------------------------------------------------------------------------
# GL_N weight sum by counting labels (no character table)
# ---------------------------------------------------------------------------
def poly_mod(a, b, q):
    """Remainder of a by monic b; polynomials are tuples of coefficients, constant term first."""
    a = list(a)
    while len(a) >= len(b):
        lead = a[-1]
        if lead:
            shift = len(a) - len(b)
            for i, coef in enumerate(b):
                a[shift + i] = (a[shift + i] - lead * coef) % q
        a.pop()
    while a and a[-1] == 0:
        a.pop()
    return a


def poly_div(a, b, q):
    a = list(a)
    quot = [0] * (len(a) - len(b) + 1)
    while len(a) >= len(b):
        lead = a[-1]
        shift = len(a) - len(b)
        quot[shift] = lead
        for i, coef in enumerate(b):
            a[shift + i] = (a[shift + i] - lead * coef) % q
        a.pop()
    return tuple(quot)


def irreducibles(q, maxdeg):
    irr = []
    for d in range(1, maxdeg + 1):
        for low in itertools.product(range(q), repeat=d):
            p = low + (1,)
            if not any(len(f) - 1 <= d // 2 and not poly_mod(p, f, q) for f in irr):
                irr.append(p)
    return irr


def factor(p, irr, q):
    out = []
    for f in irr:
        if len(f) > len(p):
            break
        m = 0
        while len(p) >= len(f) and not poly_mod(p, f, q):
            p = poly_div(p, f, q)
            m += 1
        if m:
            out.append((len(f) - 1, m))
        if len(p) == 1:
            break
    assert len(p) == 1, "factorization incomplete"
    return out


def weight_sums(q, N, irr):
    """Sum over labels (monic degree N, nonzero constant term) of the corrected and uncorrected weights."""
    corrected = uncorrected = Fraction(0)
    per_label = {}
    count = 0
    for low in itertools.product(range(q), repeat=N):
        if low[0] == 0:
            continue
        count += 1
        fac = factor(low + (1,), irr, q)
        w_c = w_u = Fraction(1)
        for d, m in fac:
            Q = q ** d
            gl = math.prod(Q ** k - 1 for k in range(1, m + 1))        # |GL_m(F_Q)|_p'
            w_c *= Fraction(Q ** (m * (m - 1) // 2), gl)               # |C(s)|_p / |C(s)|_p'
            w_u *= Fraction(1, gl)                                     # 1 / |C(s)|_p'
        corrected += w_c
        uncorrected += w_u
        per_label[low] = w_c
    return count, corrected, uncorrected, per_label


def mobius(n):
    result, k = 1, 2
    while k * k <= n:
        if n % k == 0:
            n //= k
            if n % k == 0:
                return 0
            result = -result
        k += 1
    return -result if n > 1 else result


def series_weight_sums(q, n_max):
    """Coefficients of t^0..t^n_max of prod over monic irreducible f != X of
    sum_m w(deg f, m) t^(m deg f), counting irreducibles by the necklace formula.
    Corrected w(d, m) = Q^(m(m-1)/2) / prod_k (Q^k - 1) with Q = q^d; uncorrected drops Q^(m(m-1)/2).
    The coefficient of t^N is the weight sum over all labels of GL_N(F_q)."""
    out = []
    for corrected in (True, False):
        series = [Fraction(0)] * (n_max + 1)
        series[0] = Fraction(1)
        for d in range(1, n_max + 1):
            count = sum(mobius(e) * q ** (d // e) for e in range(1, d + 1) if d % e == 0) // d
            if d == 1:
                count -= 1                                             # leave out f = X (s is invertible)
            Q = q ** d
            g = [Fraction(0)] * (n_max + 1)                           # the factor minus its constant term 1
            for m in range(1, n_max // d + 1):
                g[m * d] = Fraction(Q ** (m * (m - 1) // 2) if corrected else 1, math.prod(Q ** k - 1 for k in range(1, m + 1)))
            # (1 + g)^count = sum_k C(count, k) g^k; g^k starts at t^(k d), so k <= n_max // d
            power = [Fraction(1)] + [Fraction(0)] * n_max
            total = list(power)
            for k in range(1, n_max // d + 1):
                power = [sum(power[i] * g[n - i] for i in range(n + 1)) for n in range(n_max + 1)]
                total = [a + math.comb(count, k) * b for a, b in zip(total, power)]
            series = [sum(series[i] * total[n - i] for i in range(n + 1)) for n in range(n_max + 1)]
        out.append(series)
    return out


def main():
    print("Gelfand-Graev weights: Lemma 14.33 and Proposition 14.34 on GL_2(F_q)\n")
    for q in QS:
        out, details, wt_exact = check_gl2(q)
        print(f"q = {q}: |G| = {out['|G|']}, {out['irreducibles']} irreducibles, {out['classes']} classes ({out['seconds']} s)")
        for k, v in out.items():
            if k in ("q", "|G|", "irreducibles", "classes", "seconds"):
                continue
            print(f"   {k}: {v}")
        kinds = {}
        for lab, kind, degs, n, gg, w in details:
            kinds.setdefault((kind, tuple(gg_k[0] for gg_k in gg), w), 0)
            kinds[(kind, tuple(gg_k[0] for gg_k in gg), w)] += 1
        for (kind, gg, w), n in sorted(kinds.items(), key=lambda t: str(t)):
            print(f"   label type {kind:11s} x{n:2d}: GG constituent {gg}, wt = {w}")
        print()

    print("Weight sum over all labels of GL_N(F_q) by counting characteristic polynomials")
    print("(corrected degree |C(s)|_p [G:C(s)]_p'  versus the proof's [G:C(s)]_p')\n")
    for q in QS:
        irr = irreducibles(q, 4)
        for N in (1, 2, 3, 4):
            if q ** N > 3000:
                continue
            count, corr, unc, per_label = weight_sums(q, N, irr)
            line = f"   q = {q}, N = {N}: {count:5d} labels, corrected sum = {corr}, uncorrected sum = {unc}"
            if N == 2:
                _, _, wt_exact = check_gl2(q) if q <= 3 else (None, None, None)
                if wt_exact is not None:
                    same = all(wt_exact[((low[1]) % q, low[0])] == w for low, w in per_label.items())
                    line += f", matches GL_2 character-table weights label by label: {same}"
            print(line)

    print("\nThe same sums from the generating function, N = 1..16 (no enumeration of labels).")
    print("Proven for every N: the corrected series telescopes to 1/(1 - t); the uncorrected one is")
    print("H(t^2)/(1 - t), so its sums repeat in pairs and tend to 1/((p;p)(-p;p)^2) with p = 1/q.\n")
    for q in QS:
        corr, unc = series_weight_sums(q, 16)
        p, poch, pm = 1 / q, 1.0, 1.0
        for i in range(1, 400):
            poch *= 1 - p ** i
            pm *= 1 + p ** i
        print(f"   q = {q}: corrected sums all 1: {all(c == 1 for c in corr)}; uncorrected N = 1..6: "
              + ", ".join(str(u) for u in unc[1:7])
              + f"; equal in pairs (2K, 2K+1): {all(unc[2 * k] == unc[2 * k + 1] for k in range(1, 8))}"
              + f"; N = 16: {float(unc[16]):.12f}, limit {1 / (poch * pm * pm):.12f}")


if __name__ == "__main__":
    main()
