#!/usr/bin/env python3
"""
RSA-CRT: when partial knowledge becomes total compromise
========================================================

Two independent ways a CRT-RSA private key falls apart:

  Attack A  Leaked CRT exponent (almost all bits of d_p known)
  Attack B  One faulty CRT signature (Boneh-DeMillo-Lipton / "Bellcore")

...and one defense:

  Act 4     Verify the signature before releasing it.

Run:
    python3 rsa_crt_demo.py                  # default, ~1 second
    python3 rsa_crt_demo.py --seed 7         # different random key
    python3 rsa_crt_demo.py --missing-bits 6 # harder Attack A (16x slower)

Works in a terminal and in Colab/Jupyter. Standard library only (Python 3.8+).
READ README.md FIRST for the math behind every step.

TOY WARNING: 64-bit primes, unseeded-style PRNG, no padding. Real RSA uses
2048+ bit moduli, a CSPRNG and randomized padding. Do not reuse this code.
"""

from __future__ import annotations

import argparse
import math
import random
import sys
import time
from typing import NamedTuple, Optional

E = 65537   # the standard public exponent


# ---------------------------------------------------------------------------
# Printing helpers
# ---------------------------------------------------------------------------

def banner(title: str) -> None:
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74)


def tick(ok: bool) -> str:
    return "OK " if ok else "FAIL"


# ---------------------------------------------------------------------------
# Key generation (toy)
# ---------------------------------------------------------------------------

def is_probable_prime(n: int, rng: random.Random, rounds: int = 20) -> bool:
    """Miller-Rabin."""
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
    while True:
        cand = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_probable_prime(cand, rng):
            return cand


class Key(NamedTuple):
    p: int
    q: int
    N: int
    e: int
    d: int
    dp: int       # d mod (p-1)
    dq: int       # d mod (q-1)
    q_inv: int    # q^{-1} mod p


def make_key(bits: int, rng: random.Random, e: int = E) -> Key:
    """Retry until p != q and gcd(e, phi) == 1 (the old script could crash here)."""
    while True:
        p, q = generate_prime(bits, rng), generate_prime(bits, rng)
        if p == q:
            continue
        phi = (p - 1) * (q - 1)
        if math.gcd(e, phi) != 1:
            continue
        d = pow(e, -1, phi)
        return Key(p, q, p * q, e, d, d % (p - 1), d % (q - 1), pow(q, -1, p))


def d_from_factor(N: int, e: int, factor: int) -> int:
    """Once ONE prime factor is known, the whole private key follows."""
    other = N // factor
    assert factor * other == N
    return pow(e, -1, (factor - 1) * (other - 1))


# ---------------------------------------------------------------------------
# CRT signing (what a smart card / HSM / library does for speed)
# ---------------------------------------------------------------------------

def crt_sign(m: int, k: Key, fault_in_p_branch: bool = False) -> int:
    """
        s_p = m^dp mod p        s_q = m^dq mod q
        h   = q_inv * (s_p - s_q) mod p          (Garner recombination)
        s   = s_q + q*h
    """
    sp = pow(m, k.dp, k.p)
    sq = pow(m, k.dq, k.q)
    if fault_in_p_branch:
        sp = (sp ^ 1) % k.p          # flip one bit: a "glitch"
    h = (k.q_inv * (sp - sq)) % k.p
    return sq + k.q * h


def verified_crt_sign(m: int, k: Key, fault_in_p_branch: bool = False) -> Optional[int]:
    """DEFENSE: never release a signature that does not verify."""
    s = crt_sign(m, k, fault_in_p_branch)
    return s if pow(s, k.e, k.N) == m else None


# ---------------------------------------------------------------------------
# Attack A: leaked (almost-complete) dp
# ---------------------------------------------------------------------------

def attack_a(N: int, e: int, d0: int, missing_bits: int):
    """
    Facts:   e*dp - 1 = k*(p-1)   with 1 <= k < e,   and dp = d0 + x.
    So       k*p = e*(d0 + x) - 1 + k.
    Hence    gcd( e*(d0+x) - 1 + k , N ) = p   for the right (k, x).
    Returns (factor, k, x, number_of_guesses_tried).
    """
    tries = 0
    for k in range(1, e):
        for x in range(1 << missing_bits):
            tries += 1
            g = math.gcd(e * (d0 + x) - 1 + k, N)
            if 1 < g < N:
                return g, k, x, tries
    return None, None, None, tries


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description="RSA-CRT partial-leak and fault attacks")
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--bits", type=int, default=64, help="bits per prime (toy size)")
    ap.add_argument("--missing-bits", type=int, default=4,
                    help="low bits of dp the attacker does NOT know (default 4)")
    args, _ = ap.parse_known_args()      # ignore Colab/Jupyter's extra "-f kernel.json"
    rng = random.Random(args.seed)
    results = {}

    # ---- Act 1 ------------------------------------------------------------
    banner("ACT 1 -- A toy RSA-CRT key")
    k = make_key(args.bits, rng)
    print(f"  public   N = {k.N}")
    print(f"  public   e = {k.e}")
    print(f"  SECRET   p = {k.p}")
    print(f"  SECRET   q = {k.q}")
    print(f"  SECRET   dp = d mod (p-1) = {k.dp}")
    print("  (We print the secrets only so you can check the attacks. The attacker")
    print("   never sees p, q or d -- only the things each act says they get.)")

    # ---- Act 2 ------------------------------------------------------------
    banner("ACT 2 -- Attack A: dp leaks, except its lowest bits")
    mask = (1 << args.missing_bits) - 1
    d0 = k.dp & ~mask
    space = (E - 1) * (1 << args.missing_bits)
    print(f"  attacker knows : N, e, and d0 = dp with its lowest {args.missing_bits} bits zeroed")
    print(f"                   d0 = {d0}")
    print(f"  attacker guesses: k in [1, e-1]  and  x in [0, 2^{args.missing_bits}-1]")
    print(f"                   ({space:,} candidate pairs in the worst case)")
    t0 = time.perf_counter()
    g, kk, x, tries = attack_a(k.N, k.e, d0, args.missing_bits)
    dt = time.perf_counter() - t0
    ok_factor = g in (k.p, k.q)
    print(f"\n  found after {tries:,} guesses ({dt:.2f}s):  k = {kk},  x = {x}")
    print(f"  true missing bits were            x = {k.dp & mask}   -> {tick(x == (k.dp & mask))}")
    print(f"  recovered prime factor            {g}   -> {tick(ok_factor)}")
    d_rec = d_from_factor(k.N, k.e, g)
    print(f"  recovered private exponent d      {d_rec}")
    print(f"  equals the real d?                {tick(d_rec == k.d)}")
    results["A: factor recovered"] = ok_factor
    results["A: missing bits recovered exactly"] = x == (k.dp & mask)
    results["A: full private key d rebuilt"] = d_rec == k.d

    # ---- Act 3 ------------------------------------------------------------
    banner("ACT 3 -- Attack B: one faulty signature")
    m = 123456789
    assert math.gcd(m, k.N) == 1
    s_good = crt_sign(m, k)
    s_bad = crt_sign(m, k, fault_in_p_branch=True)
    print(f"  message           m = {m}")
    print(f"  correct signature s  = {s_good}   verifies? {pow(s_good, k.e, k.N) == m}")
    print(f"  faulty  signature s' = {s_bad}   verifies? {pow(s_bad, k.e, k.N) == m}")
    print("\n  The glitch hit the p-branch only, so s' is still correct modulo q.")
    val = pow(s_bad, k.e, k.N) - m
    print(f"  attacker computes  v = s'^e - m")
    print(f"     v mod q == 0 ?  {val % k.q == 0}   (s' is right modulo q)")
    print(f"     v mod p == 0 ?  {val % k.p == 0}   (s' is wrong modulo p)")
    f1 = math.gcd(val, k.N)
    print(f"  gcd(v, N) = {f1}   is q? {f1 == k.q}")
    f2 = math.gcd(s_good - s_bad, k.N)
    print(f"\n  Bonus: if the attacker holds BOTH signatures, no need to know m:")
    print(f"  gcd(s - s', N) = {f2}   is q? {f2 == k.q}")
    d_rec_b = d_from_factor(k.N, k.e, f1)
    print(f"\n  recovered d equals the real d?    {tick(d_rec_b == k.d)}")
    results["B: factor recovered from s'^e - m"] = f1 == k.q
    results["B: factor recovered from s - s' (no m)"] = f2 == k.q
    results["B: full private key d rebuilt"] = d_rec_b == k.d

    # ---- Act 4 ------------------------------------------------------------
    banner("ACT 4 -- Defense: verify before you release")
    ok_sig = verified_crt_sign(m, k)
    glitched = verified_crt_sign(m, k, fault_in_p_branch=True)
    print(f"  no fault : signature released      -> {ok_sig == s_good}")
    print(f"  with fault: output                 -> {glitched}   (suppressed)")
    print("  The attacker gets nothing to run gcd on.")
    results["Defense: honest signing still works"] = ok_sig == s_good
    results["Defense: faulty signature suppressed"] = glitched is None

    # ---- Scoreboard ---------------------------------------------------------
    banner("SCOREBOARD -- every claim, machine-checked")
    for label, ok in results.items():
        print(f"  [{'v' if ok else 'X'}] {label}")
    all_ok = all(results.values())
    print("\n  ALL CHECKS PASSED" if all_ok else "\n  SOMETHING FAILED")
    return 0 if all_ok else 1


def _in_notebook() -> bool:
    return "ipykernel" in sys.modules or "google.colab" in sys.modules


if __name__ == "__main__":
    code = main()
    if not _in_notebook():
        sys.exit(code)

"""
==========================================================================
ACT 1 -- A toy RSA-CRT key
==========================================================================
  public   N = 208949316736487877719575705392500277853
  public   e = 65537
  SECRET   p = 12916442222160946727
  SECRET   q = 16177002392964696539
  SECRET   dp = d mod (p-1) = 2867999255640113169
  (We print the secrets only so you can check the attacks. The attacker
   never sees p, q or d -- only the things each act says they get.)

==========================================================================
ACT 2 -- Attack A: dp leaks, except its lowest bits
==========================================================================
  attacker knows : N, e, and d0 = dp with its lowest 4 bits zeroed
                   d0 = 2867999255640113168
  attacker guesses: k in [1, e-1]  and  x in [0, 2^4-1]
                   (1,048,576 candidate pairs in the worst case)

  found after 232,818 guesses (1.30s):  k = 14552,  x = 1
  true missing bits were            x = 1   -> OK 
  recovered prime factor            12916442222160946727   -> OK 
  recovered private exponent d      169542362575281381081629845451553335421
  equals the real d?                OK 

==========================================================================
ACT 3 -- Attack B: one faulty signature
==========================================================================
  message           m = 123456789
  correct signature s  = 67994187799788993705083935259464190969   verifies? True
  faulty  signature s' = 169121646381512570962181718547268706795   verifies? False

  The glitch hit the p-branch only, so s' is still correct modulo q.
  attacker computes  v = s'^e - m
     v mod q == 0 ?  True   (s' is right modulo q)
     v mod p == 0 ?  False   (s' is wrong modulo p)
  gcd(v, N) = 16177002392964696539   is q? True

  Bonus: if the attacker holds BOTH signatures, no need to know m:
  gcd(s - s', N) = 16177002392964696539   is q? True

  recovered d equals the real d?    OK 

==========================================================================
ACT 4 -- Defense: verify before you release
==========================================================================
  no fault : signature released      -> True
  with fault: output                 -> None   (suppressed)
  The attacker gets nothing to run gcd on.

==========================================================================
SCOREBOARD -- every claim, machine-checked
==========================================================================
  [v] A: factor recovered
  [v] A: missing bits recovered exactly
  [v] A: full private key d rebuilt
  [v] B: factor recovered from s'^e - m
  [v] B: factor recovered from s - s' (no m)
  [v] B: full private key d rebuilt
  [v] Defense: honest signing still works
  [v] Defense: faulty signature suppressed

  ALL CHECKS PASSED

"""