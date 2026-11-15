#!/usr/bin/env python3
"""
The ZK Hash Trapdoor: when "good diffusion" becomes deception
=============================================================

A self-contained, standard-library-only demo of ONE specific idea:

    A linear mixing matrix M can pass every standard check (it is MDS,
    invertible, "well diffused") and still hide an invariant direction that
    the nonlinear layer never touches during Poseidon-style PARTIAL rounds.

Run:    python3 zk_hash_trapdoor_demo.py
        python3 zk_hash_trapdoor_demo.py --pause      # press Enter between acts
        python3 zk_hash_trapdoor_demo.py --seed 7     # different random trapdoor

READ THE README FIRST if you want the 2-minute mental model.

IMPORTANT: this is a TOY. It is Poseidon-STYLE (same round shape, same field
as BN254, same x^5 S-box), but the weak matrix is engineered by us, and Act 3
uses a deliberately simplified permutation to isolate the mechanism. This
demo does NOT break real Poseidon. Act 4 shows exactly why not.
"""

from __future__ import annotations

import argparse
import hashlib
import random
import sys
import time
from itertools import combinations
from math import gcd
from typing import List, NamedTuple, Optional, Tuple

# --------------------------------------------------------------------------
# Public parameters (the "boring" part -- identical for honest and weak hash)
# --------------------------------------------------------------------------

# BN254 scalar field: the field used by many Groth16 / circom ZK circuits.
P = 0x30644E72E131A029B85045B68181585D2833E84879B9709143E1F593F0000001
ALPHA = 5          # S-box: x -> x^5 (a bijection because gcd(5, P-1) == 1)
WIDTH = 3          # state = (lane0, lane1, lane2)
RF_HALF = 4        # full rounds before and after the partial layer
RP = 57            # partial rounds (S-box on lane 0 ONLY)

assert gcd(ALPHA, P - 1) == 1, "S-box must be a permutation of the field"

Vec = List[int]
Mat = List[List[int]]


class Config(NamedTuple):
    name: str
    rf_first: int
    rp: int
    rf_last: int

    def schedule(self) -> List[str]:
        return ["F"] * self.rf_first + ["P"] * self.rp + ["F"] * self.rf_last


# Act 3 only: strips the permutation down to its partial layer.
SIMPLIFIED = Config("simplified: 57 partial rounds only", 0, RP, 0)
# Acts 4: the realistic Poseidon round structure.
FULL = Config("full structure: 4 full + 57 partial + 4 full", RF_HALF, RP, RF_HALF)


# --------------------------------------------------------------------------
# Tiny pretty-printing helpers
# --------------------------------------------------------------------------

USE_COLOR = sys.stdout.isatty()


def col(text: str, code: str) -> str:
    return f"\033[{code}m{text}\033[0m" if USE_COLOR else text


def good(t: str) -> str: return col(t, "32")
def bad(t: str) -> str: return col(t, "31")
def dim(t: str) -> str: return col(t, "2")
def bold(t: str) -> str: return col(t, "1")


def short(v: int) -> str:
    """First 8 hex digits of a field element -- enough to eyeball equality."""
    return "0x" + f"{v % P:064x}"[:8] + "…"


def banner(title: str) -> None:
    print()
    print(bold("=" * 78))
    print(bold(title))
    print(bold("=" * 78))


def say(text: str = "") -> None:
    print(text)


# --------------------------------------------------------------------------
# The Poseidon-style permutation   x_{i+1} = M * S(x_i) + c_i
# --------------------------------------------------------------------------

def round_constants(n: int) -> List[Vec]:
    """Nothing-up-my-sleeve constants: SHA-256 of a label. SAME for both matrices."""
    return [
        [
            int.from_bytes(hashlib.sha256(f"zk-trapdoor-demo:{r}:{j}".encode()).digest(), "big") % P
            for j in range(WIDTH)
        ]
        for r in range(n)
    ]


def mat_vec(M: Mat, v: Vec) -> Vec:
    return [sum(M[i][j] * v[j] for j in range(WIDTH)) % P for i in range(WIDTH)]


def permute(state: Vec, M: Mat, schedule: List[str], rc: List[Vec],
            trace: Optional[List[Tuple[str, Vec]]] = None) -> Vec:
    """Apply one round per entry of `schedule` ('F' = full, 'P' = partial)."""
    assert len(schedule) == len(rc)
    x = list(state)
    for kind, c in zip(schedule, rc):
        if kind == "F":                                   # S-box on every lane
            s = [pow(v, ALPHA, P) for v in x]
        else:                                             # S-box on lane 0 only
            s = [pow(x[0], ALPHA, P)] + x[1:]
        mixed = mat_vec(M, s)
        x = [(mixed[i] + c[i]) % P for i in range(WIDTH)]
        if trace is not None:
            trace.append((kind, list(x)))
    return x


def zk_hash(m1: int, m2: int, M: Mat, cfg: Config) -> int:
    """2-to-1 hash, circomlib style: state = (0, m1, m2); output = lane 0."""
    sched = cfg.schedule()
    return permute([0, m1 % P, m2 % P], M, sched, round_constants(len(sched)))[0]


# --------------------------------------------------------------------------
# Matrices: the honest one, the weak one, and the "standard audit"
# --------------------------------------------------------------------------

def det3(M: Mat) -> int:
    a, b, c = M[0]; d, e, f = M[1]; g, h, i = M[2]
    return (a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)) % P


def is_mds_3x3(M: Mat) -> bool:
    """MDS <=> every square submatrix is invertible (1x1, 2x2, 3x3)."""
    if any(v % P == 0 for row in M for v in row):
        return False
    for r in combinations(range(3), 2):
        for c in combinations(range(3), 2):
            minor = (M[r[0]][c[0]] * M[r[1]][c[1]] - M[r[0]][c[1]] * M[r[1]][c[0]]) % P
            if minor == 0:
                return False
    return det3(M) != 0


def honest_cauchy_matrix() -> Mat:
    """Cauchy matrix M[i][j] = 1/(x_i + y_j). Provably MDS; the textbook choice."""
    xs, ys = [0, 1, 2], [3, 4, 5]
    return [[pow(x + y, -1, P) for y in ys] for x in xs]


class Trapdoor(NamedTuple):
    matrix: Mat
    lam: int        # eigenvalue of the hidden direction
    direction: Tuple[int, int]   # (d1, d2): hidden eigenvector is (0, d1, d2)


def build_weak_matrix(rng: random.Random) -> Trapdoor:
    """
    The adversarial parameter generator.

    Goal: an MDS matrix M with an eigenvector v = (0, d1, d2) whose FIRST
    coordinate is zero. In terms of the three rows of M v = lam * v:

        row 0:  M01*d1 + M02*d2 = 0          <- the nonlinear lane never "sees" v
        row 1:  M11*d1 + M12*d2 = lam*d1
        row 2:  M21*d1 + M22*d2 = lam*d2

    Pick the first column and (M01, M11, M21) at random, then SOLVE for the
    last column. Retry until the result is MDS (so it passes the standard check).
    """
    while True:
        d1, d2 = rng.randrange(1, P), rng.randrange(1, P)
        lam = rng.randrange(2, P)
        m00, m10, m20, m01, m11, m21 = (rng.randrange(1, P) for _ in range(6))
        inv_d2 = pow(d2, -1, P)
        m02 = (-m01 * d1 * inv_d2) % P
        m12 = ((lam - m11) * d1 * inv_d2) % P
        m22 = ((lam * d2 - m21 * d1) * inv_d2) % P
        M = [[m00, m01, m02], [m10, m11, m12], [m20, m21, m22]]
        if is_mds_3x3(M):
            return Trapdoor(M, lam, (d1, d2))


def find_lane0_free_eigenvector(M: Mat) -> Optional[Tuple[int, Tuple[int, int]]]:
    """
    The "deeper audit": is there an eigenvector v = (0, d1, d2) ?

    Row 0 of M v must vanish, which forces (d1, d2) ∝ (M02, -M01): a UNIQUE
    candidate direction. So the test is just: is that candidate an eigenvector?
    Returns (lambda, (d1, d2)) if yes, else None.
    """
    d1, d2 = M[0][2] % P, (-M[0][1]) % P
    Mv = mat_vec(M, [0, d1, d2])
    lam = Mv[1] * pow(d1, -1, P) % P
    if Mv[0] == 0 and Mv[1] == lam * d1 % P and Mv[2] == lam * d2 % P:
        return lam, (d1, d2)
    return None


# --------------------------------------------------------------------------
# Demo acts
# --------------------------------------------------------------------------

def print_matrix(name: str, M: Mat) -> None:
    say(f"  {name}")
    for row in M:
        say("     [ " + "  ".join(short(v) for v in row) + " ]")


def act1(honest: Mat, weak: Mat) -> dict:
    banner("ACT 1 -- Two public matrices. Both pass the standard checks.")
    say("""
Same field, same x^5 S-box, same round counts, same round constants.
The ONLY thing that differs between the two hashes is the 3x3 matrix M.
""")
    print_matrix("Honest matrix (Cauchy construction):", honest)
    say()
    print_matrix("Weak matrix (from an adversarial parameter generator):", weak)
    say()
    say("  Standard parameter validation:")
    say(f"  {'check':<42}{'honest':>10}{'weak':>10}")
    rows = [
        ("S-box x^5 is a bijection (gcd(5,p-1)=1)", True, True),
        ("M is invertible", det3(honest) != 0, det3(weak) != 0),
        ("M is MDS (all square minors nonzero)", is_mds_3x3(honest), is_mds_3x3(weak)),
    ]
    for label, h, w in rows:
        cell = lambda ok: (good if ok else bad)(("PASS" if ok else "FAIL").rjust(10))
        say(f"  {label:<42}{cell(h)}{cell(w)}")
    say("\n  -> By these checks the two matrices are indistinguishable.")
    return {"both_mds": is_mds_3x3(honest) and is_mds_3x3(weak)}


def act2(honest: Mat, weak: Mat, trap: Trapdoor) -> dict:
    banner("ACT 2 -- The deeper audit: look for an eigenvector with lane-0 = 0")
    say("""
Why that specific shape?  In a PARTIAL round the S-box touches lane 0 only.
If a difference vector v = (0, d1, d2) is an eigenvector (M v = lam v), then:
   - lanes 1,2 carry no S-box, so v passes straight through the S-box layer,
   - M maps v to lam*v, which is still in the same direction,
   - and lane 0 of lam*v is still 0, so the S-box never sees the difference.
Row 0 of M v must be 0, which pins the candidate direction to (M02, -M01).
So the audit is one line of algebra: is that candidate an eigenvector?
""")
    h = find_lane0_free_eigenvector(honest)
    w = find_lane0_free_eigenvector(weak)
    say(f"  honest matrix : " + (bad("eigenvector FOUND") if h else good("none found (as expected)")))
    say(f"  weak matrix   : " + (bad("eigenvector FOUND") if w else good("none found")))
    if w:
        lam, (d1, d2) = w
        say(f"\n  hidden direction  v = (0, {short(d1)}, {short(d2)})")
        say(f"  eigenvalue        lam = {short(lam)}")
        same_dir = (d1 * trap.direction[1] - d2 * trap.direction[0]) % P == 0
        say(f"  matches the planted trapdoor (direction up to scale, same lam): {same_dir and lam == trap.lam}")
    say("""
Note: the matrix is PUBLIC. Anyone can run this audit. The danger is that
nobody was asked to, because the usual checklist (MDS, invertible) says "fine".""")
    return {"honest_clean": h is None, "weak_found": w is not None}


def table_row(label: str, a: Vec, b: Vec, prev_d1: Optional[int], lam: Optional[int]) -> int:
    d1 = (a[1] - b[1]) % P
    same0 = a[0] == b[0]
    ratio = ""
    if prev_d1 is not None and lam is not None:
        ratio = good("= lam") if (d1 * pow(prev_d1, -1, P)) % P == lam else bad("!= lam")
    say(f"  {label:>5} | {short(a[0])} | {short(b[0])} | "
        f"{good('EQUAL') if same0 else bad('differ'):>14} | {short(d1)} | {ratio}")
    return d1


def act3(honest: Mat, weak: Mat, trap: Trapdoor, rng: random.Random) -> dict:
    banner("ACT 3 -- The exploit: lane 0 never notices the difference")
    say(f"""
Setting: {SIMPLIFIED.name}.
(This is a deliberate simplification to isolate the mechanism. Act 4 puts
 the full round structure back.)

Hash:  H(m1, m2) = lane 0 of  Perm( 0, m1, m2 ).

Attack recipe  (needs only the PUBLIC matrix):
  1. Find the hidden eigenvector (0, d1, d2) -- Act 2.
  2. Pick any message (m1, m2) and any nonzero s.
  3. Second message:  (m1', m2') = (m1 + s*d1,  m2 + s*d2).
""")
    lam, (d1, d2) = find_lane0_free_eigenvector(weak)  # type: ignore[misc]
    m1, m2 = rng.randrange(P), rng.randrange(P)
    s = rng.randrange(1, P)
    n1, n2 = (m1 + s * d1) % P, (m2 + s * d2) % P
    say(f"  message A : ({short(m1)}, {short(m2)})")
    say(f"  message B : ({short(n1)}, {short(n2)})     <- differs from A by s*(d1,d2)")

    sched = SIMPLIFIED.schedule()
    rc = round_constants(len(sched))
    ta: List[Tuple[str, Vec]] = []
    tb: List[Tuple[str, Vec]] = []
    out_a = permute([0, m1, m2], weak, sched, rc, ta)[0]
    out_b = permute([0, n1, n2], weak, sched, rc, tb)[0]

    say("\n  Watching both states round by round (weak matrix):")
    say(f"  {'round':>5} | {'A: lane 0':<12} | {'B: lane 0':<12} | {'lane 0 A vs B':>14} | "
        f"{'diff in lane 1':<12} | diff(i)/diff(i-1)")
    say("  " + "-" * 82)
    prev = table_row("input", [0, m1, m2], [0, n1, n2], None, None)
    n = len(ta)
    show = set(range(1, 5)) | {n - 1, n}
    last_shown = 0
    d1_chain = prev
    ratios_ok = True
    for r in range(1, n + 1):
        a, b = ta[r - 1][1], tb[r - 1][1]
        d_now = (a[1] - b[1]) % P
        ratios_ok &= (d_now * pow(d1_chain, -1, P)) % P == lam
        if r in show:
            if r - last_shown > 1:
                say(f"  {'⋮':>5} |")
            table_row(str(r), a, b, d1_chain, lam)
            last_shown = r
        d1_chain = d_now
    say(f"\n  Every round: lane 0 identical, difference multiplied by exactly lam.  ({good('verified for all ' + str(n) + ' rounds') if ratios_ok else bad('FAILED')})")

    say(f"\n  H(A) = {short(out_a)}")
    say(f"  H(B) = {short(out_b)}")
    collision = out_a == out_b and (m1, m2) != (n1, n2)
    say("  " + (good(">>> COLLISION: two different inputs, same hash. <<<") if collision else bad("no collision")))

    # Same recipe against the honest matrix.
    h_lam_d = (honest[0][2] % P, (-honest[0][1]) % P)
    hd1, hd2 = h_lam_d
    hn1, hn2 = (m1 + s * hd1) % P, (m2 + s * hd2) % P
    h_a = zk_hash(m1, m2, honest, SIMPLIFIED)
    h_b = zk_hash(hn1, hn2, honest, SIMPLIFIED)
    say("\n  Same recipe, honest matrix (using its row-0-orthogonal direction):")
    say(f"    H(A) = {short(h_a)}\n    H(B) = {short(h_b)}   -> "
        + (bad("collision?!") if h_a == h_b else good("different (recipe fails, as it should)")))

    # Cost comparison
    t0 = time.perf_counter()
    reps = 300
    for k in range(reps):
        zk_hash(k, k + 1, weak, FULL)
    rate = reps / (time.perf_counter() - t0)
    birthday = 2 ** 127   # generic collision cost for a ~254-bit output
    years = birthday / rate / (365.25 * 86400)
    say(f"""
  Cost to find a collision:
     weak matrix   : 2 hash evaluations + one subtraction/multiplication
     honest matrix : generic birthday bound ~ 2^127 evaluations
                     (at this laptop's ~{rate:,.0f} hashes/s ≈ {years:.1e} years)""")
    return {"ratios_ok": ratios_ok, "collision": collision, "honest_resists": h_a != h_b}


def act4(honest: Mat, weak: Mat, trap: Trapdoor, rng: random.Random) -> dict:
    banner("ACT 4 -- Reality check: a weak matrix is NOT automatically a break")
    say(f"""
Now the realistic structure: {FULL.name}.
Same attack recipe, same weak matrix.
""")
    lam, (d1, d2) = find_lane0_free_eigenvector(weak)  # type: ignore[misc]
    m1, m2 = rng.randrange(P), rng.randrange(P)
    s = rng.randrange(1, P)
    n1, n2 = (m1 + s * d1) % P, (m2 + s * d2) % P
    h_a, h_b = zk_hash(m1, m2, weak, FULL), zk_hash(n1, n2, weak, FULL)
    say("  (a) Attack through the full hash:")
    say(f"      H(A) = {short(h_a)}\n      H(B) = {short(h_b)}   -> "
        + (bad("collision") if h_a == h_b else good("different: NO collision")))
    say("      The first 4 full rounds scramble the difference out of the hidden")
    say("      direction before it ever reaches the partial layer.")

    # Permutation-level: start INSIDE the partial layer.
    sched_all = FULL.schedule()
    rc_all = round_constants(len(sched_all))
    lo, hi = FULL.rf_first, FULL.rf_first + FULL.rp
    mid_a = [rng.randrange(P), rng.randrange(P), rng.randrange(P)]
    mid_b = [mid_a[0], (mid_a[1] + s * d1) % P, (mid_a[2] + s * d2) % P]
    after_a = permute(mid_a, weak, sched_all[lo:hi], rc_all[lo:hi])
    after_b = permute(mid_b, weak, sched_all[lo:hi], rc_all[lo:hi])
    lane0_same = after_a[0] == after_b[0]
    diff = [(after_b[i] - after_a[i]) % P for i in range(WIDTH)]   # B - A
    expected = [0, s * d1 * pow(lam, FULL.rp, P) % P, s * d2 * pow(lam, FULL.rp, P) % P]
    trail_ok = lane0_same and diff == expected
    say("\n  (b) Permutation-level view: start two states INSIDE the partial layer,")
    say("      differing by s*(0,d1,d2):")
    say(f"      lane 0 after 57 partial rounds equal? {good('yes') if lane0_same else bad('no')}")
    say(f"      difference == lam^57 * s*(0,d1,d2)  ? {good('yes') if trail_ok else bad('no')}")
    say("      -> the invariant trail is REAL and holds with probability 1.")

    out_a = permute(after_a, weak, sched_all[hi:], rc_all[hi:])
    out_b = permute(after_b, weak, sched_all[hi:], rc_all[hi:])
    say("\n  (c) ...then the last 4 full rounds:")
    say(f"      output A lane 0 = {short(out_a[0])}\n      output B lane 0 = {short(out_b[0])}   -> "
        + (bad("collision") if out_a[0] == out_b[0] else good("different: still no collision")))

    say("""
  What this shows
    * The weakness is real: 57 rounds of "nonlinear protection" collapse to
      linear scaling along one direction. A cryptanalyst now has a free,
      probability-1 differential across the entire partial layer.
    * It is NOT a collision on the full hash. The surrounding full rounds
      and the capacity=0 input constraint stand in the way.
    * Turning a trail into a real break needs a concrete attack on the
      COMPLETE construction. That is exactly the claim the README makes.
""")
    return {"no_full_collision": h_a != h_b, "trail_ok": trail_ok, "full_out_differs": out_a[0] != out_b[0]}


def scoreboard(results: dict) -> bool:
    banner("SCOREBOARD -- every claim above, machine-checked")
    checks = [
        ("Both matrices pass the standard MDS check", results["both_mds"]),
        ("Audit finds NO hidden eigenvector in honest matrix", results["honest_clean"]),
        ("Audit finds the hidden eigenvector in weak matrix", results["weak_found"]),
        ("Lane 0 identical + difference scales by lam every round", results["ratios_ok"]),
        ("Simplified hash: collision with weak matrix", results["collision"]),
        ("Simplified hash: same recipe fails on honest matrix", results["honest_resists"]),
        ("Full structure: invariant trail holds across partial layer", results["trail_ok"]),
        ("Full structure: NO collision (weakness != break)", results["no_full_collision"] and results["full_out_differs"]),
    ]
    for label, ok in checks:
        say(f"  [{good('✓') if ok else bad('✗')}] {label}")
    return all(ok for _, ok in checks)


def main() -> int:
    ap = argparse.ArgumentParser(description="ZK hash trapdoor demo")
    ap.add_argument("--seed", type=int, default=1337, help="RNG seed (default 1337)")
    ap.add_argument("--pause", action="store_true", help="wait for Enter between acts")
    # parse_known_args: tolerate extra flags injected by Jupyter/Colab (e.g. -f kernel.json)
    args, _ = ap.parse_known_args()
    rng = random.Random(args.seed)

    def pause() -> None:
        if args.pause:
            input(dim("\n[press Enter to continue]"))

    banner("THE ZK HASH TRAPDOOR  (toy, Poseidon-style, BN254 field)")
    say(
        "Question: can a PUBLIC, well-diffusing matrix still hide a weakness?\n"
        "Answer in 4 acts. Total runtime: a few seconds.")
    t0 = time.perf_counter()
    honest = honest_cauchy_matrix()
    trap = build_weak_matrix(rng)
    weak = trap.matrix

    results: dict = {}
    results.update(act1(honest, weak)); pause()
    results.update(act2(honest, weak, trap)); pause()
    results.update(act3(honest, weak, trap, rng)); pause()
    results.update(act4(honest, weak, trap, rng))
    ok = scoreboard(results)
    say(f"\n(done in {time.perf_counter() - t0:.1f}s)")
    return 0 if ok else 1


def _in_notebook() -> bool:
    return "ipykernel" in sys.modules or "google.colab" in sys.modules


if __name__ == "__main__":
    code = main()
    if not _in_notebook():      # in Colab/Jupyter, sys.exit() would raise SystemExit
        sys.exit(code)

"""

==============================================================================
THE ZK HASH TRAPDOOR  (toy, Poseidon-style, BN254 field)
==============================================================================
Question: can a PUBLIC, well-diffusing matrix still hide a weakness?
Answer in 4 acts. Total runtime: a few seconds.

==============================================================================
ACT 1 -- Two public matrices. Both pass the standard checks.
==============================================================================

Same field, same x^5 S-box, same round counts, same round constants.
The ONLY thing that differs between the two hashes is the 3x3 matrix M.

  Honest matrix (Cauchy construction):
     [ 0x2042def7…  0x244b3ad6…  0x135b5294… ]
     [ 0x244b3ad6…  0x135b5294…  0x285396b5… ]
     [ 0x135b5294…  0x285396b5…  0x06e9c210… ]

  Weak matrix (from an adversarial parameter generator):
     [ 0x01240e4a…  0x28dabb01…  0x2401c062… ]
     [ 0x0af539ea…  0x049f8b32…  0x06444368… ]
     [ 0x26b1a349…  0x0118d48b…  0x0ff24533… ]

  Standard parameter validation:
  check                                         honest      weak
  S-box x^5 is a bijection (gcd(5,p-1)=1)         PASS      PASS
  M is invertible                                 PASS      PASS
  M is MDS (all square minors nonzero)            PASS      PASS

  -> By these checks the two matrices are indistinguishable.

==============================================================================
ACT 2 -- The deeper audit: look for an eigenvector with lane-0 = 0
==============================================================================

Why that specific shape?  In a PARTIAL round the S-box touches lane 0 only.
If a difference vector v = (0, d1, d2) is an eigenvector (M v = lam v), then:
   - lanes 1,2 carry no S-box, so v passes straight through the S-box layer,
   - M maps v to lam*v, which is still in the same direction,
   - and lane 0 of lam*v is still 0, so the S-box never sees the difference.
Row 0 of M v must be 0, which pins the candidate direction to (M02, -M01).
So the audit is one line of algebra: is that candidate an eigenvector?

  honest matrix : none found (as expected)
  weak matrix   : eigenvector FOUND

  hidden direction  v = (0, 0x2401c062…, 0x07899371…)
  eigenvalue        lam = 0x28a7136c…
  matches the planted trapdoor (direction up to scale, same lam): True

Note: the matrix is PUBLIC. Anyone can run this audit. The danger is that
nobody was asked to, because the usual checklist (MDS, invertible) says "fine".

==============================================================================
ACT 3 -- The exploit: lane 0 never notices the difference
==============================================================================

Setting: simplified: 57 partial rounds only.
(This is a deliberate simplification to isolate the mechanism. Act 4 puts
 the full round structure back.)

Hash:  H(m1, m2) = lane 0 of  Perm( 0, m1, m2 ).

Attack recipe  (needs only the PUBLIC matrix):
  1. Find the hidden eigenvector (0, d1, d2) -- Act 2.
  2. Pick any message (m1, m2) and any nonzero s.
  3. Second message:  (m1', m2') = (m1 + s*d1,  m2 + s*d2).

  message A : (0x2ff170a8…, 0x02eeebb5…)
  message B : (0x07410774…, 0x2f9e7e72…)     <- differs from A by s*(d1,d2)

  Watching both states round by round (weak matrix):
  round | A: lane 0    | B: lane 0    |  lane 0 A vs B | diff in lane 1 | diff(i)/diff(i-1)
  ----------------------------------------------------------------------------------
  input | 0x00000000… | 0x00000000… |          EQUAL | 0x28b06933… | 
      1 | 0x204264de… | 0x204264de… |          EQUAL | 0x209a4c65… | = lam
      2 | 0x05c5956e… | 0x05c5956e… |          EQUAL | 0x0a9021c8… | = lam
      3 | 0x0d184749… | 0x0d184749… |          EQUAL | 0x2b0a6760… | = lam
      4 | 0x04c820cd… | 0x04c820cd… |          EQUAL | 0x14fe8c78… | = lam
      ⋮ |
     56 | 0x0f81b075… | 0x0f81b075… |          EQUAL | 0x06b860e1… | = lam
     57 | 0x1855f081… | 0x1855f081… |          EQUAL | 0x1d8a185b… | = lam

  Every round: lane 0 identical, difference multiplied by exactly lam.  (verified for all 57 rounds)

  H(A) = 0x1855f081…
  H(B) = 0x1855f081…
  >>> COLLISION: two different inputs, same hash. <<<

  Same recipe, honest matrix (using its row-0-orthogonal direction):
    H(A) = 0x1e199957…
    H(B) = 0x13ed23d4…   -> different (recipe fails, as it should)

  Cost to find a collision:
     weak matrix   : 2 hash evaluations + one subtraction/multiplication
     honest matrix : generic birthday bound ~ 2^127 evaluations
                     (at this laptop's ~300 hashes/s ≈ 1.8e+28 years)

==============================================================================
ACT 4 -- Reality check: a weak matrix is NOT automatically a break
==============================================================================

Now the realistic structure: full structure: 4 full + 57 partial + 4 full.
Same attack recipe, same weak matrix.

  (a) Attack through the full hash:
      H(A) = 0x03664ef8…
      H(B) = 0x0b76ff8c…   -> different: NO collision
      The first 4 full rounds scramble the difference out of the hidden
      direction before it ever reaches the partial layer.

  (b) Permutation-level view: start two states INSIDE the partial layer,
      differing by s*(0,d1,d2):
      lane 0 after 57 partial rounds equal? yes
      difference == lam^57 * s*(0,d1,d2)  ? yes
      -> the invariant trail is REAL and holds with probability 1.

  (c) ...then the last 4 full rounds:
      output A lane 0 = 0x275dc8df…
      output B lane 0 = 0x17deabce…   -> different: still no collision

  What this shows
    * The weakness is real: 57 rounds of "nonlinear protection" collapse to
      linear scaling along one direction. A cryptanalyst now has a free,
      probability-1 differential across the entire partial layer.
    * It is NOT a collision on the full hash. The surrounding full rounds
      and the capacity=0 input constraint stand in the way.
    * Turning a trail into a real break needs a concrete attack on the
      COMPLETE construction. That is exactly the claim the README makes.


==============================================================================
SCOREBOARD -- every claim above, machine-checked
==============================================================================
  [✓] Both matrices pass the standard MDS check
  [✓] Audit finds NO hidden eigenvector in honest matrix
  [✓] Audit finds the hidden eigenvector in weak matrix
  [✓] Lane 0 identical + difference scales by lam every round
  [✓] Simplified hash: collision with weak matrix
  [✓] Simplified hash: same recipe fails on honest matrix
  [✓] Full structure: invariant trail holds across partial layer
  [✓] Full structure: NO collision (weakness != break)

(done in 1.0s)

"""
