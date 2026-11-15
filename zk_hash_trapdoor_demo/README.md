# The ZK Hash Trapdoor — when diffusion becomes deception

A single-file, dependency-free Python demo of **one idea**:

> A public mixing matrix can pass every standard check (MDS, invertible, "good diffusion")
> and still hide a direction that the nonlinear S-box **never touches** during Poseidon-style
> *partial* rounds. Along that direction, 57 rounds of "nonlinear protection" collapse into
> multiplication by a constant.

```bash
python3 zk_hash_trapdoor_demo.py            # runs in < 1 second
python3 zk_hash_trapdoor_demo.py --pause    # step through act by act
python3 zk_hash_trapdoor_demo.py --seed 7   # a different random weak matrix
```

Python 3.8+, standard library only. Exit code `0` means every claim the demo makes was machine-checked.

---

## Read this first: the 60-second mental model

A Poseidon-style round is `x_next = M · S(x) + c`.

In a **partial round**, the S-box `S` is applied to **lane 0 only**. Lanes 1 and 2 pass through untouched.

```
        lane 0   lane 1   lane 2
          │        │        │
        [ x^5 ]    │        │        <- S-box touches ONLY lane 0
          │        │        │
        ┌─┴────────┴────────┴─┐
        │      matrix  M      │      <- linear mixing
        └─┬────────┬────────┬─┘
          + c0     + c1     + c2     <- round constants
```

Now take two states that differ by a vector `v = (0, d1, d2)`: **lane 0 is identical**.

1. Lane 0 is identical → the S-box produces identical output → the S-box contributes **no difference**.
2. Lanes 1 and 2 had no S-box anyway → the difference `v` goes into `M` unchanged.
3. If `M·v = λ·v` (v is an **eigenvector**), the difference comes out as `λ·v`: same direction, still `0` in lane 0.

So the difference never reaches lane 0, so the S-box never sees it, so it repeats forever:

```
round 0:  difference = v
round 1:  difference = λ·v      lane 0 still identical
round 2:  difference = λ²·v     lane 0 still identical
   ⋮
round 57: difference = λ⁵⁷·v    lane 0 STILL identical
```

**The crux, in one sentence:** an eigenvector whose lane-0 coordinate is zero is invisible to the
only nonlinear component in the partial layer.

If the hash outputs lane 0 (as circomlib-style Poseidon does), two different inputs with that
difference produce the **same output**.

---

## Why would anyone choose such a matrix? (and why checks miss it)

An MDS matrix only promises that every square submatrix is invertible. That is a statement about
*branch numbers* and diffusion. It says nothing about whether `M` has an eigenvector with a zero
in a particular coordinate.

An adversarial parameter generator builds the matrix **backwards**:

1. Pick the hidden direction `v = (0, d1, d2)` and eigenvalue `λ`.
2. Pick the first column and `M01, M11, M21` at random.
3. **Solve** the last column so that `M·v = λ·v` holds. (Three linear equations.)
4. Retry until the matrix is MDS. It usually is on the first or second try.

Result: a matrix that looks random, passes the usual checklist, and carries a planted structure.

> **Note on the word "trapdoor":** there is no secret key here. The matrix is public and the
> eigenvector is recoverable by anyone who looks (Act 2 does it in one line). The adversary's
> advantage is *knowing where to look* while everyone else's checklist says "fine".

---

## What you will see, act by act

| Act | What happens | What to look for |
|-----|--------------|------------------|
| **1** | Honest (Cauchy) matrix vs. weak matrix | Both print `PASS` on every standard check. They are indistinguishable by that checklist. |
| **2** | The "deeper audit": is there an eigenvector `(0, d1, d2)`? | Honest: none. Weak: **found**, with its eigenvalue `λ`. |
| **3** | The exploit on a *simplified* hash (57 partial rounds only) | A round-by-round table: **lane 0 `EQUAL` in every row**, and the lane-1 difference grows by exactly `= lam` each round. Then `H(A) == H(B)` for different inputs. The identical recipe fails on the honest matrix. |
| **4** | The reality check on the *full* Poseidon round structure | The collision **disappears**. The invariant trail inside the partial layer still holds (probability 1), but the hash outputs differ. |
| **Scoreboard** | 8 machine-checked claims | All should show `✓`. |

### The table to stare at in Act 3

```
  round | A: lane 0    | B: lane 0    |  lane 0 A vs B | diff in lane 1 | diff(i)/diff(i-1)
  input | 0x00000000… | 0x00000000… |          EQUAL | 0x28b06933… |
      1 | 0x204264de… | 0x204264de… |          EQUAL | 0x209a4c65… | = lam
      2 | 0x05c5956e… | 0x05c5956e… |          EQUAL | 0x0a9021c8… | = lam
```

Two columns carry the whole story: lane 0 never differs, and the difference scales by the same
constant every round. That is the invariant subspace, made visible.

---

## What this demo is and is not

**Real (matches actual Poseidon-style constructions):**
- BN254 scalar field, `x^5` S-box, width 3, 8 full + 57 partial rounds
- Round formula `x_next = M·S(x) + c`; partial rounds apply S to lane 0 only
- Hash convention: state `(0, m1, m2)`, output lane 0 (circomlib style)
- The eigenvector condition is the genuine mechanism behind published invariant-subspace /
  subspace-trail analyses of partial layers

**Toy (be honest about these):**
- The weak matrix is **engineered by us**. This is not a flaw in any deployed parameter set.
- Act 3 **removes the full rounds** to isolate the mechanism. Real Poseidon surrounds the partial
  layer with full rounds, and that matters (Act 4).
- Round constants are SHA-256 derived, not the official Poseidon generation procedure.
- Width is 3 so the audit is one line. Larger widths need a subspace-trail algorithm, not a single check.

### Why Act 4 exists

A weak matrix does **not** automatically give a collision or forge a proof. In the full structure:

- The first 4 full rounds apply the S-box to every lane, which knocks a hash input out of the
  invariant direction before it reaches the partial layer.
- The last 4 full rounds re-expose any lane-0 equality to the S-box.
- The hash fixes the capacity lane to 0, which constrains what an attacker can feed in.

What the weakness **does** buy: a free, probability-1 differential across the whole partial layer,
meaning the 57 "expensive" nonlinear rounds contribute almost nothing along that direction.
Whether that can be extended to a collision, preimage, or Merkle forgery is a question for
analysis of the **complete construction**, not something the matrix alone answers.

---

## Why this matters for zero-knowledge systems

A ZK proof is only as sound as what's under it. Hash collisions or weak algebraic structure could
undermine Merkle-tree membership, commitments, or nullifiers, depending on how the hash is used.
The proof system can be flawless while the hash beneath it fails its assumptions.

The deeper lesson: **public parameters are not automatically trustworthy parameters.** Passing the
checklist shows the checklist was satisfied, not that the construction is free of structure.

### Defenses (what practitioners do)

- **Verifiable, reproducible parameter generation** ("nothing up my sleeve"): derive the matrix from
  a transparent seed so there is no freedom to engineer it backwards.
- **Run subspace-trail / invariant-subspace checks** on the matrix, not only the MDS check.
- **Independent third-party audits** of the *complete* parameter set, including round counts and constants.
- **Prefer widely reviewed parameter sets** over bespoke ones.

---

## File tour

| Function | Purpose |
|----------|---------|
| `permute()` | The Poseidon-style permutation, `x ← M·S(x) + c`, full or partial rounds |
| `zk_hash()` | 2-to-1 hash: `Perm(0, m1, m2)[0]` |
| `honest_cauchy_matrix()` | Provably-MDS textbook matrix |
| `build_weak_matrix()` | **The adversarial generator**: solves for a matrix with the planted eigenvector |
| `is_mds_3x3()` | The "standard check" both matrices pass |
| `find_lane0_free_eigenvector()` | The "deeper audit" that catches it |
| `act1()` … `act4()` | The four demo stages |

## FAQ

**Does this break Poseidon?**
No. It breaks a toy with a matrix we deliberately weakened, and Act 4 shows the full structure resists
the naive attack.

**Is the matrix secret?**
No. It's public. That is the point: public does not mean audited for *this* property.

**Why does the honest-matrix attack "fail"?** For the honest matrix, `(M02, −M01)` is not an eigenvector,
so the difference leaks into lane 0 on round 1, the S-box amplifies it, and outputs diverge.
Generic collision cost then falls back to the birthday bound, about 2^127 for a ~254-bit output.

**Can I try a different structure?** Edit `build_weak_matrix()` or change `RP` / `RF_HALF` at the top of the
script, then rerun. The scoreboard will tell you which claims still hold.

## Further reading

Search for: *Poseidon hash function (Grassi et al.)*, *invariant subspace attacks*, and
*subspace trails in partial SPN rounds / choosing the linear layer*. Verify citations against the
originals; this README is a teaching aid, not a cryptanalysis paper.

## Running in Colab / Jupyter

Paste the whole script into a cell and run it, or upload it and run `%run zk_hash_trapdoor_demo.py`.
The script ignores the extra `-f kernel.json` flag that notebooks inject, and does not call
`sys.exit()` there. To use options in a notebook: `%run zk_hash_trapdoor_demo.py --seed 7`.
