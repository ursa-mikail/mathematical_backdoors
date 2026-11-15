#!/usr/bin/env python3
"""
LLL / Coppersmith tutorial: recover a small secret hidden in a modular equation
================================================================================

SCENARIO  ("stereotyped message", RSA with e = 3)
  * Public key (N, e = 3).  Ciphertext  c = m^3 mod N.
  * The message is  m = M + secret,  where
        M       is KNOWN (a fixed header / template), low bits all zero,
        secret  is UNKNOWN, 0 <= secret < X   (X = 2^unknown_bits).
  * Brute force costs X tries.  LLL finds the secret in a fraction of a second.

THE IDEA IN ONE LINE
  Turn "find a root of a polynomial MODULO N" (hard) into
       "find a short vector in a lattice"       (LLL does this),
  and the short vector is a polynomial whose root is EXACT over the integers (easy).

Run:
    python3 lll_coppersmith_demo.py                       # narrated walk-through
    python3 lll_coppersmith_demo.py --unknown-bits 52     # watch it FAIL (too many bits)
    python3 lll_coppersmith_demo.py --sweep               # table: where does it stop working?
    python3 lll_coppersmith_demo.py --seed 7              # a different random instance

Needs:  Python 3.8+ and sympy (only for the very last, easy step; Colab has it).
READ README.md FIRST: every STEP below has the same number there.

TOY WARNING: 300-bit modulus, e = 3, no padding. Illustrates the method; not an
attack on modern RSA (random padding + e = 65537 + 2048-bit N).
"""

from __future__ import annotations

import argparse
import math
import random
import sys
import time
from fractions import Fraction
from typing import List, NamedTuple, Tuple

try:
    import sympy
except ImportError as exc:   # raise (not sys.exit) so notebooks show a clean message
    raise ImportError("This demo needs sympy for its last step:  pip install sympy") from exc

E = 3                      # the (historically common) small public exponent
DEFAULT_DELTA = Fraction(99, 100)   # LLL quality parameter


# ---------------------------------------------------------------------------
# Printing helpers
# ---------------------------------------------------------------------------

def banner(title: str) -> None:
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74)


def short(n: int) -> str:
    """Big number -> leading digits + size, so tables stay readable."""
    s = str(abs(n))
    lead = s if len(s) <= 14 else s[:12] + "..."
    return f"{'-' if n < 0 else ''}{lead} ({n.bit_length()} bits)"


# ---------------------------------------------------------------------------
# Toy RSA instance
# ---------------------------------------------------------------------------

def is_probable_prime(n: int, rng: random.Random, rounds: int = 16) -> bool:
    if n < 4:
        return n in (2, 3)
    if n % 2 == 0:
        return False
    r, s = 0, n - 1
    while s % 2 == 0:
        r += 1
        s //= 2
    for _ in range(rounds):
        a = rng.randrange(2, n - 1)
        x = pow(a, s, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def generate_prime(bits: int, rng: random.Random) -> int:
    """Prime with p mod 3 != 1, so that gcd(3, p-1) = 1 and e = 3 is a valid RSA exponent."""
    while True:
        cand = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if cand % 3 != 1 and is_probable_prime(cand, rng):
            return cand


class Instance(NamedTuple):
    N: int        # public modulus
    c: int        # ciphertext  m^3 mod N
    M: int        # known part of the message (low `unknown_bits` bits are zero)
    X: int        # bound on the unknown part: 0 <= secret < X
    secret: int   # the unknown part (printed only so you can check the answer)


def make_instance(rng: random.Random, bits_p: int, unknown_bits: int) -> Instance:
    while True:
        p, q = generate_prime(bits_p, rng), generate_prime(bits_p, rng)
        if p != q:
            break
    N = p * q
    X = 1 << unknown_bits
    secret = rng.randrange(X)
    M = rng.randrange(N >> unknown_bits) << unknown_bits
    m = M + secret
    assert m < N
    return Instance(N, pow(m, E, N), M, X, secret)


# ---------------------------------------------------------------------------
# STEP 1: the polynomial.   f(x) = (M + x)^3 - c   has the secret as a root mod N
# ---------------------------------------------------------------------------

def polynomial_coefficients(inst: Instance) -> Tuple[int, int, int]:
    """
    (M+x)^3 - c = x^3 + 3M x^2 + 3M^2 x + (M^3 - c).
    Returns (f2, f1, f0) = coefficients of x^2, x^1, x^0, each reduced mod N.
    Reducing mod N is allowed: it changes f by multiples of N, so f(secret) is
    still 0 mod N -- and it keeps the numbers (and the lattice) small.
    """
    N, M, c = inst.N, inst.M, inst.c
    return (3 * M) % N, (3 * M * M) % N, (M ** 3 - c) % N


# ---------------------------------------------------------------------------
# STEP 2: the lattice
# ---------------------------------------------------------------------------

def build_lattice(N: int, X: int, f2: int, f1: int, f0: int) -> List[List[int]]:
    """
    Rows = coefficient vectors [x^0, x^1, x^2, x^3] of four polynomials that are
    ALL 0 mod N at x = secret, with column j multiplied by X^j:

        row1 =  N            ->  [ N,   0,      0,       0   ]
        row2 =  N*x          ->  [ 0,   N*X,    0,       0   ]
        row3 =  N*x^2        ->  [ 0,   0,      N*X^2,   0   ]
        row4 =  f(x)         ->  [ f0,  f1*X,   f2*X^2,  X^3 ]
    """
    return [
        [N, 0, 0, 0],
        [0, N * X, 0, 0],
        [0, 0, N * X ** 2, 0],
        [f0, f1 * X, f2 * X ** 2, X ** 3],
    ]


# ---------------------------------------------------------------------------
# STEP 3: LLL (plain textbook version, exact rational arithmetic)
# ---------------------------------------------------------------------------

def _dot(u, v) -> Fraction:
    return sum(Fraction(a) * Fraction(b) for a, b in zip(u, v))


def _gram_schmidt(B):
    n = len(B)
    Bs = [[Fraction(0)] * len(B[0]) for _ in range(n)]
    mu = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        Bs[i] = [Fraction(x) for x in B[i]]
        for j in range(i):
            mu[i][j] = _dot(B[i], Bs[j]) / _dot(Bs[j], Bs[j])
            for k in range(len(B[i])):
                Bs[i][k] -= mu[i][j] * Bs[j][k]
    return Bs, mu


def lll_reduce(B, delta: Fraction = DEFAULT_DELTA) -> List[List[int]]:
    """
    Repeat until stable:
      (1) SIZE-REDUCE: subtract integer multiples of earlier rows so each row is
          as 'un-aligned' with the earlier ones as possible  (like Euclid's step).
      (2) LOVASZ TEST: if row k is much shorter than row k-1 (in the
          Gram-Schmidt sense), SWAP them and step back.
    Exact Fractions: no floating-point surprises (slow, but fine for 4x4).
    """
    B = [list(map(int, row)) for row in B]
    n = len(B)
    Bs, mu = _gram_schmidt(B)
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            if abs(mu[k][j]) > Fraction(1, 2):
                q = round(mu[k][j])
                B[k] = [a - q * b for a, b in zip(B[k], B[j])]
                Bs, mu = _gram_schmidt(B)
        if _dot(Bs[k], Bs[k]) >= (delta - mu[k][k - 1] ** 2) * _dot(Bs[k - 1], Bs[k - 1]):
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            Bs, mu = _gram_schmidt(B)
            k = max(k - 1, 1)
    return B


def norm_sq(v) -> int:
    return sum(int(c) ** 2 for c in v)


def bits_of_norm(v) -> float:
    return math.log2(norm_sq(v)) / 2


# ---------------------------------------------------------------------------
# STEP 4: from lattice vector back to a polynomial, and the size test
# ---------------------------------------------------------------------------

def vector_to_polynomial(vec, X: int) -> List[int]:
    """Undo the column scaling: [h0, h1, h2, h3] with h(x) = h0 + h1 x + h2 x^2 + h3 x^3."""
    out = []
    for j, comp in enumerate(vec):
        assert comp % (X ** j) == 0, "column j of any lattice vector is a multiple of X^j"
        out.append(comp // X ** j)
    return out


def howgrave_graham_holds(vec, N: int) -> bool:
    """
    Howgrave-Graham condition for a polynomial with 4 coefficients:
        || h(X*x) ||_2  <  N / sqrt(4)      <=>      4 * norm^2  <  N^2
    If true, then h(secret) = 0 EXACTLY (not just mod N).
    """
    return 4 * norm_sq(vec) < N * N


def howgrave_graham_margin_bits(vec, N: int) -> float:
    """log2( N/2 ) - log2( ||vec|| ).   Positive = condition holds, by that many bits."""
    return math.log2(N) - 1 - bits_of_norm(vec)


def provable_bound_bits(N: int, delta: Fraction = DEFAULT_DELTA) -> float:
    """
    Largest log2(X) for which LLL's *worst-case guarantee* proves success.
      ||b1|| <= alpha^(3/4) * det^(1/4),  det = N^3 X^6,  alpha = 1/(delta - 1/4)
      need   ||b1|| < N/2     =>   X < N^(1/6) / (2 * alpha^(3/4))^(2/3)
    In practice LLL usually does a little better than its guarantee.
    """
    alpha = 1 / (float(delta) - 0.25)
    return math.log2(N) / 6 - (2 / 3) * math.log2(2 * alpha ** 0.75)


# ---------------------------------------------------------------------------
# STEP 5: solve h(x) = 0 exactly  (ordinary algebra; sympy does the factoring)
# ---------------------------------------------------------------------------

def integer_roots_in_range(h: List[int], X: int) -> List[int]:
    x = sympy.symbols("x")
    poly = sympy.Poly(sum(int(c) * x ** i for i, c in enumerate(h)), x)
    if poly.is_zero:
        return []
    return sorted(int(r) for r in poly.ground_roots() if r.is_integer and 0 <= r < X)


# ---------------------------------------------------------------------------
# One-call version of steps 1-5 (used by --sweep)
# ---------------------------------------------------------------------------

class Outcome(NamedTuple):
    recovered: bool
    hg_holds: bool
    margin_bits: float


def solve_instance(inst: Instance, delta: Fraction = DEFAULT_DELTA) -> Outcome:
    f2, f1, f0 = polynomial_coefficients(inst)
    reduced = lll_reduce(build_lattice(inst.N, inst.X, f2, f1, f0), delta)
    best = min(reduced, key=norm_sq)
    h = vector_to_polynomial(best, inst.X)
    roots = integer_roots_in_range(h, inst.X)
    return Outcome(inst.secret in roots,
                   howgrave_graham_holds(best, inst.N),
                   howgrave_graham_margin_bits(best, inst.N))


# ---------------------------------------------------------------------------
# Narrated walk-through
# ---------------------------------------------------------------------------

def run_walkthrough(args) -> bool:
    rng = random.Random(args.seed)

    # ---- STEP 0 -------------------------------------------------------------
    banner("STEP 0 -- The scenario")
    inst = make_instance(rng, args.bits_p, args.unknown_bits)
    N, X = inst.N, inst.X
    print(f"  public N  = {N}   ({N.bit_length()} bits)")
    print(f"  public e  = {E}")
    print(f"  public c  = m^3 mod N = {short(inst.c)}")
    print(f"  KNOWN  M  = {short(inst.M)}      (the message's known high part)")
    print(f"  UNKNOWN secret: 0 <= secret < 2^{args.unknown_bits} = {X:,}")
    print(f"     (we print it only to check the final answer: secret = {inst.secret})")
    print(f"  message is  m = M + secret")

    t0 = time.perf_counter()
    tries = 20000
    for i in range(tries):
        _ = pow(inst.M + i, E, N) == inst.c
    rate = tries / (time.perf_counter() - t0)
    hours = X / rate / 3600
    print(f"\n  Brute force: up to {X:,} tries. Measured on this machine: ~{rate:,.0f} tries/s")
    print(f"  -> worst case ~{hours:,.1f} hours (~{hours / 24:,.1f} days). LLL below: a fraction of a second.")

    bound = provable_bound_bits(N)
    print(f"\n  How many unknown bits can THIS method handle for this N?")
    print(f"     provable limit ~ {bound:.1f} bits   (you chose {args.unknown_bits})")

    # ---- STEP 1 -------------------------------------------------------------
    banner("STEP 1 -- Turn the secret into a root of a polynomial modulo N")
    f2, f1, f0 = polynomial_coefficients(inst)
    print("""
  We know  (M + secret)^3 = c  (mod N).   Define  f(x) = (M + x)^3 - c.
  Expanding:   f(x) = x^3 + (3M) x^2 + (3M^2) x + (M^3 - c)
  Then  f(secret) = 0 (mod N).   Note f is MONIC (x^3 has coefficient 1): this
  matters in Step 2.  "= 0 mod N" is hard to solve directly: the true value of
  f(secret) could be 0, N, -7N, ... we cannot tell which.
""")
    print(f"  f(x) = x^3 + f2*x^2 + f1*x + f0   (coefficients reduced mod N)")
    print(f"     f2 = {short(f2)}")
    print(f"     f1 = {short(f1)}")
    print(f"     f0 = {short(f0)}")
    sec = inst.secret
    f_at_secret = (sec ** 3 + f2 * sec ** 2 + f1 * sec + f0) % N
    print(f"  check: f(secret) mod N == 0 ?  {f_at_secret == 0}")

    # ---- STEP 2 -------------------------------------------------------------
    banner("STEP 2 -- Build a lattice whose vectors are polynomials that vanish mod N")
    print("""
  Four polynomials, ALL equal to 0 mod N at x = secret:
      row1 = N,   row2 = N*x,   row3 = N*x^2,   row4 = f(x)
  Any integer combination of them also vanishes mod N at the secret.
  Write each as the vector of its coefficients [x^0, x^1, x^2, x^3], and
  multiply column j by X^j.  (Reason: the secret is at most X, so a term
  h_j * secret^j is at most h_j * X^j. The scaling makes the vector's length
  measure exactly "how big can this polynomial be at the secret".)

      [  N      0       0        0   ]
      [  0     N*X      0        0   ]
      [  0      0     N*X^2      0   ]
      [  f0   f1*X   f2*X^2     X^3  ]
""")
    B = build_lattice(N, X, f2, f1, f0)
    print("  size of each row's length, in bits (log2 of the norm):")
    for i, row in enumerate(B, 1):
        print(f"     row{i}: {bits_of_norm(row):7.1f} bits")
    det_bits = 3 * math.log2(N) + 6 * math.log2(X)
    print(f"  lattice determinant = N^3 * X^6  ->  {det_bits:.1f} bits   (this number decides success)")

    # ---- STEP 3 -------------------------------------------------------------
    banner("STEP 3 -- Run LLL: find a SHORT vector in that lattice")
    t0 = time.perf_counter()
    reduced = lll_reduce(B)
    print(f"  LLL finished in {time.perf_counter() - t0:.3f} s")
    print("  row lengths after LLL, in bits:")
    for i, row in enumerate(reduced, 1):
        print(f"     row{i}: {bits_of_norm(row):7.1f} bits")
    best = min(reduced, key=norm_sq)
    print(f"  shortest vector: {bits_of_norm(best):.1f} bits   (was >= {min(bits_of_norm(r) for r in B):.1f} bits)")

    # ---- STEP 4 -------------------------------------------------------------
    banner("STEP 4 -- Is the short vector small enough to be an EXACT equation?")
    h = vector_to_polynomial(best, X)
    print("  Undo the X^j scaling to read the vector as a polynomial h(x):")
    for j in range(3, -1, -1):
        print(f"     h{j} = {short(h[j])}")
    margin = howgrave_graham_margin_bits(best, N)
    holds = howgrave_graham_holds(best, N)
    print(f"""
  Howgrave-Graham test:  ||h(X*x)||  <  N / 2  ?
     ||h(X*x)||  = 2^{bits_of_norm(best):.1f}
     N / 2       = 2^{math.log2(N) - 1:.1f}
     margin      = {margin:+.1f} bits   -> {'HOLDS' if holds else 'FAILS'}
""")
    if holds:
        print("  Because it holds, h(secret) is not just 0 mod N: it is EXACTLY 0.")
        print("  (Proof sketch in README: |h(secret)| < N and N divides h(secret).)")
    else:
        print("  It does not hold: h(secret) may be a nonzero multiple of N. Expect failure.")
        if h[1] == h[2] == h[3] == 0:
            print("  Here the 'short vector' is just the constant polynomial h(x) = N (row1):")
            print("  LLL found nothing shorter, so there is no useful equation to solve.")
        print("  (The test is sufficient, not necessary, so success is still possible.)")

    # ---- STEP 5 -------------------------------------------------------------
    banner("STEP 5 -- Solve h(x) = 0 over the integers, read off the secret")
    roots = integer_roots_in_range(h, X)
    print(f"  integer roots of h(x) in [0, 2^{args.unknown_bits}): {roots}")
    ok = inst.secret in roots
    if ok:
        m_rec = inst.M + inst.secret
        print(f"\n  recovered secret  = {inst.secret}")
        print(f"  recovered message = M + secret = {short(m_rec)}")
        print(f"  pow(m, 3, N) == c ?  {pow(m_rec, E, N) == inst.c}")
        print("\n  SECRET RECOVERED. No brute force: one 4x4 lattice reduction.")
    else:
        print("\n  FAILED to recover the secret.")
        print(f"  This lattice proves success only up to ~{bound:.1f} unknown bits for this N;")
        print(f"  you asked for {args.unknown_bits}. See README, section 'Why it stops working'.")
    return ok


# ---------------------------------------------------------------------------
# Sweep: where does it stop working?
# ---------------------------------------------------------------------------

def run_sweep(args) -> bool:
    banner("SWEEP -- success rate vs. number of unknown bits")
    sample = make_instance(random.Random(0), args.bits_p, 8)
    print(f"  N is ~{sample.N.bit_length()} bits; provable limit ~ {provable_bound_bits(sample.N):.1f} bits")
    print(f"  {args.trials} random instances per row\n")
    print(f"  {'unknown bits':>12} | {'recovered':>9} | {'HG test holds':>13} | {'avg margin (bits)':>17}")
    print("  " + "-" * 62)
    for ub in (38, 42, 46, 48, 50, 52, 54):
        outs = [solve_instance(make_instance(random.Random(1000 * ub + t), args.bits_p, ub))
                for t in range(args.trials)]
        rec = sum(o.recovered for o in outs)
        hg = sum(o.hg_holds for o in outs)
        avg = sum(o.margin_bits for o in outs) / len(outs)
        print(f"  {ub:>12} | {rec:>4}/{args.trials:<4} | {hg:>8}/{args.trials:<4} | {avg:>+17.1f}")
    print("\n  Read it as: margin > 0 means the exactness test passes; it shrinks as the unknown grows.")
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description="LLL / Coppersmith small-root tutorial")
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--unknown-bits", type=int, default=40, help="bits of the message the attacker lacks")
    ap.add_argument("--bits-p", type=int, default=150, help="bits per RSA prime (N is about twice this)")
    ap.add_argument("--sweep", action="store_true", help="show success rate vs unknown bits")
    ap.add_argument("--trials", type=int, default=6, help="instances per row in --sweep")
    args, _ = ap.parse_known_args()      # ignore Colab/Jupyter's extra "-f kernel.json"

    ok = run_sweep(args) if args.sweep else run_walkthrough(args)
    return 0 if ok else 1


def _in_notebook() -> bool:
    return "ipykernel" in sys.modules or "google.colab" in sys.modules


if __name__ == "__main__":
    code = main()
    if not _in_notebook():
        sys.exit(code)

"""

==========================================================================
STEP 0 -- The scenario
==========================================================================
  public N  = 1831834123143190515947723801745646040780083037629283658263745698660031071344938158440493419   (300 bits)
  public e  = 3
  public c  = m^3 mod N = 792418144469... (299 bits)
  KNOWN  M  = 767174621134... (299 bits)      (the message's known high part)
  UNKNOWN secret: 0 <= secret < 2^40 = 1,099,511,627,776
     (we print it only to check the final answer: secret = 660408492882)
  message is  m = M + secret

  Brute force: up to 1,099,511,627,776 tries. Measured on this machine: ~409,652 tries/s
  -> worst case ~745.6 hours (~31.1 days). LLL below: a fraction of a second.

  How many unknown bits can THIS method handle for this N?
     provable limit ~ 49.1 bits   (you chose 40)

==========================================================================
STEP 1 -- Turn the secret into a root of a polynomial modulo N
==========================================================================

  We know  (M + secret)^3 = c  (mod N).   Define  f(x) = (M + x)^3 - c.
  Expanding:   f(x) = x^3 + (3M) x^2 + (3M^2) x + (M^3 - c)
  Then  f(secret) = 0 (mod N).   Note f is MONIC (x^3 has coefficient 1): this
  matters in Step 2.  "= 0 mod N" is hard to solve directly: the true value of
  f(secret) could be 0, N, -7N, ... we cannot tell which.

  f(x) = x^3 + f2*x^2 + f1*x + f0   (coefficients reduced mod N)
     f2 = 469689740259... (298 bits)
     f1 = 345203746513... (295 bits)
     f0 = 665451713135... (299 bits)
  check: f(secret) mod N == 0 ?  True

==========================================================================
STEP 2 -- Build a lattice whose vectors are polynomials that vanish mod N
==========================================================================

  Four polynomials, ALL equal to 0 mod N at x = secret:
      row1 = N,   row2 = N*x,   row3 = N*x^2,   row4 = f(x)
  Any integer combination of them also vanishes mod N at the secret.
  Write each as the vector of its coefficients [x^0, x^1, x^2, x^3], and
  multiply column j by X^j.  (Reason: the secret is at most X, so a term
  h_j * secret^j is at most h_j * X^j. The scaling makes the vector's length
  measure exactly "how big can this polynomial be at the secret".)

      [  N      0       0        0   ]
      [  0     N*X      0        0   ]
      [  0      0     N*X^2      0   ]
      [  f0   f1*X   f2*X^2     X^3  ]

  size of each row's length, in bits (log2 of the norm):
     row1:   299.8 bits
     row2:   339.8 bits
     row3:   379.8 bits
     row4:   377.9 bits
  lattice determinant = N^3 * X^6  ->  1139.5 bits   (this number decides success)

==========================================================================
STEP 3 -- Run LLL: find a SHORT vector in that lattice
==========================================================================
  LLL finished in 0.451 s
  row lengths after LLL, in bits:
     row1:   280.0 bits
     row2:   280.0 bits
     row3:   280.3 bits
     row4:   299.5 bits
  shortest vector: 280.0 bits   (was >= 299.8 bits)

==========================================================================
STEP 4 -- Is the short vector small enough to be an EXACT equation?
==========================================================================
  Undo the X^j scaling to read the vector as a polynomial h(x):
     h3 = -699790520238... (159 bits)
     h2 = 488454455040... (199 bits)
     h1 = 120783619334... (240 bits)
     h0 = -809138734254... (279 bits)

  Howgrave-Graham test:  ||h(X*x)||  <  N / 2  ?
     ||h(X*x)||  = 2^280.0
     N / 2       = 2^298.8
     margin      = +18.9 bits   -> HOLDS

  Because it holds, h(secret) is not just 0 mod N: it is EXACTLY 0.
  (Proof sketch in README: |h(secret)| < N and N divides h(secret).)

==========================================================================
STEP 5 -- Solve h(x) = 0 over the integers, read off the secret
==========================================================================
  integer roots of h(x) in [0, 2^40): [660408492882]

  recovered secret  = 660408492882
  recovered message = M + secret = 767174621134... (299 bits)
  pow(m, 3, N) == c ?  True

  SECRET RECOVERED. No brute force: one 4x4 lattice reduction.

"""