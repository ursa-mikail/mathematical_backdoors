# RSA-CRT: When Partial Knowledge Becomes Total Compromise

A single-file, standard-library-only Python demo of **two ways an RSA private key can be fully
recovered from partial information**, plus the one-line defense against the second.

```bash
python3 rsa_crt_demo.py                   # ~1 second
python3 rsa_crt_demo.py --seed 7          # a different random key
python3 rsa_crt_demo.py --missing-bits 6  # harder Attack A (about 4x more work; see FAQ)
```

Exit code `0` means every claim printed by the demo was machine-checked.
Math is written in LaTeX (renders on GitHub, VS Code, Jupyter; otherwise read the plain text in the code).

---

## 1. The one-minute picture

**The shared idea.** An RSA modulus is $N = p\,q$. Anyone who learns **one** of the primes can compute the other
($q = N/p$), then the whole private key. So an attack only needs to *find a single prime factor*.
Both attacks below are clever ways to make a factor fall out of a tiny piece of extra information.

| | **Attack A** | **Attack B** |
|---|---|---|
| Name | Leaked CRT exponent | Fault attack (Boneh–DeMillo–Lipton, "Bellcore") |
| What the attacker gets | $N$, $e$, and $d_p$ **except its lowest few bits** | $N$, $e$, the message $m$, and **one wrong signature** $s'$ |
| Key trick | $e\,d_p - 1$ is a multiple of $p-1$ | $s'$ is right mod $q$ but wrong mod $p$ |
| How a factor appears | $\gcd(\,e(d_0+x)-1+k,\ N\,) = p$ | $\gcd(\,s'^{\,e}-m,\ N\,) = q$ |
| Work | $\le (e-1)\cdot 2^{b}$ cheap gcds | **1** gcd |
| Defense | Never expose CRT exponents | Verify the signature before releasing it |

---

## 2. Notation (keep this open while reading)

| Symbol | Meaning | Public? |
|---|---|---|
| $p,\ q$ | the two secret primes | secret |
| $N = pq$ | the modulus | **public** |
| $e$ | public exponent (here $65537$) | **public** |
| $d$ | private exponent, $e\,d \equiv 1 \pmod{\varphi(N)}$ | secret |
| $\varphi(N) = (p-1)(q-1)$ | Euler's totient | secret |
| $d_p = d \bmod (p-1)$ | CRT exponent for the $p$ side | secret |
| $d_q = d \bmod (q-1)$ | CRT exponent for the $q$ side | secret |
| $q_{\text{inv}} = q^{-1} \bmod p$ | CRT recombination constant | secret |
| $m$ | message (a number $0<m<N$) | known to attacker in Attack B |
| $s$ | a correct signature $s = m^d \bmod N$ | public once released |
| $s'$ | a **faulty** signature | what the attacker captures |

---

## 3. Background: why CRT signing exists (and why it is fragile)

Plain RSA signing is $s = m^d \bmod N$, one huge exponentiation.
The **Chinese Remainder Theorem (CRT)** lets a device do two *half-size* exponentiations instead
(roughly 4x faster), which is why smart cards, HSMs, and libraries use it:

$$s_p = m^{d_p} \bmod p, \qquad s_q = m^{d_q} \bmod q$$

Then recombine so that $s \equiv s_p \pmod p$ and $s \equiv s_q \pmod q$ (Garner's formula):

$$h = q_{\text{inv}}\,(s_p - s_q) \bmod p, \qquad s = s_q + q\,h$$

*Why this works:* modulo $q$ the term $q\,h$ vanishes, so $s \equiv s_q$. Modulo $p$,
$s \equiv s_q + q\cdot q^{-1}(s_p - s_q) \equiv s_p$. ✔

**The fragility:** the key now contains extra secrets ($d_p,\ d_q,\ q_{\text{inv}}$), and the two
branches are computed *separately*. Leaking one of those extras, or corrupting one branch,
exposes a prime.

---

## 4. Attack A: a (nearly) leaked $d_p$

### 4.1 Scenario
A side channel, a buggy log file, or a memory dump reveals $d_p$ **except for its lowest $b$ bits**.
The attacker writes

$$d_p = d_0 + x, \qquad d_0 \text{ known},\quad 0 \le x < 2^{b}.$$

(Demo default: $b=4$, so there are only 16 possibilities for $x$.)

### 4.2 The identity everything rests on
Since $e\,d \equiv 1 \pmod{\varphi(N)}$ and $(p-1) \mid \varphi(N)$, we get $e\,d_p \equiv 1 \pmod{p-1}$. So there is an integer $k$ with

$$e\,d_p - 1 = k\,(p-1).$$

**$k$ is small.** Because $1 \le d_p < p-1$,

$$0 < e\,d_p - 1 < e\,(p-1) \quad\Longrightarrow\quad 1 \le k \le e-1.$$

For $e = 65537$ that is **at most 65,536 candidates**, regardless of how big $p$ is.

### 4.3 Making $p$ appear
Add $k$ to both sides:

$$k\,p \;=\; e\,d_p - 1 + k \;=\; e\,(d_0 + x) - 1 + k. \tag{$\ast$}$$

The right side uses only things the attacker can guess ($k$, $x$) or knows ($e$, $d_0$).
For the **correct** pair $(k, x)$ the value is $k\,p$, which is divisible by $p$ but **not** by $q$
(because $q$ is a huge prime and $k < e < q$). Hence

$$\gcd\bigl(e\,(d_0+x) - 1 + k,\ N\bigr) = p.$$

For a **wrong** pair the number is essentially random, so the gcd is $1$. The code accepts a result only
when $1 < g < N$.

### 4.4 The algorithm
```text
for k in 1 .. e-1:
    for x in 0 .. 2^b - 1:
        g = gcd( e*(d0 + x) - 1 + k , N )
        if 1 < g < N:  →  g is a prime factor of N, stop
```
Then: $q = N/p$, $\varphi = (p-1)(q-1)$, and $d = e^{-1} \bmod \varphi$. **The whole private key is rebuilt.**

### 4.5 Worked example (from the demo, seed 2024)
The demo prints $k = 14552$. You can check it by hand: $e\,d_p - 1 = 65537\cdot d_p - 1$ is exactly
$14552\cdot(p-1)$, with $p = 12916442222160946727$ and $d_p = 2867999255640113169$.
The hidden low bits were $x = 1$, and the attack finds the same.

### 4.6 Honest limits
- With only $b=4$ unknown bits this is **brute force over a tiny space**, not a deep lattice attack.
  The cost is $(e-1)\cdot 2^{b}$ gcds, so every extra missing bit doubles the work, and it becomes
  infeasible if most of $d_p$ is missing.
- Real research on partial key exposure (Coppersmith's method, Blömer–May and follow-ups) handles much larger
  missing fractions with lattice techniques. This demo shows the **core algebraic identity**, not those methods.

---

## 5. Attack B: a single faulty signature

### 5.1 Scenario
A voltage glitch, laser, or software bug corrupts **only the $p$-branch** of one signature:
$s_p$ is replaced by a wrong value $s_p'$, while $s_q$ stays correct. The device outputs the faulty $s'$.
(The demo flips one bit: $s_p' = s_p \oplus 1$.)

### 5.2 What the faulty signature looks like
Run Garner's recombination with the wrong $s_p'$:

$$h' = q_{\text{inv}}(s_p' - s_q) \bmod p, \qquad s' = s_q + q\,h'.$$

Reading $s'$ modulo each prime:

$$s' \equiv s_q \pmod q \quad(\text{the } qh' \text{ term vanishes}), \qquad s' \equiv s_p' \pmod p.$$

### 5.3 Right mod $q$, wrong mod $p$
**Modulo $q$:** since $e\,d_q \equiv 1 \pmod{q-1}$, Fermat's little theorem gives $s_q^{\,e} = m^{d_q e} \equiv m \pmod q$. So

$$s'^{\,e} - m \equiv 0 \pmod q.$$

**Modulo $p$:** $s'^{\,e} \equiv s_p'^{\,e}$. The map $x \mapsto x^e$ is a bijection mod $p$ (because $\gcd(e,p-1)=1$)
and the true $s_p$ satisfies $s_p^{\,e} \equiv m$. So $s_p'^{\,e} \equiv m$ would force $s_p' = s_p$, which is false. Therefore

$$s'^{\,e} - m \not\equiv 0 \pmod p.$$

### 5.4 The one-line attack
$s'^{\,e}-m$ is a multiple of $q$ and **not** of $p$, so

$$\boxed{\ \gcd\bigl(s'^{\,e} - m,\ N\bigr) = q\ }$$

One faulty signature, one gcd, full key. The demo prints the two residue checks ($v \bmod q = 0$, $v \bmod p \ne 0$)
so you can see this rather than trust it.

### 5.5 Bonus: no need to know $m$
If the attacker has the correct $s$ *and* the faulty $s'$ for the same message, then
$s - s' = q\,(h - h')$ is a multiple of $q$, while modulo $p$ it is $s_p - s_p' \ne 0$. So

$$\gcd(s - s',\ N) = q.$$

The demo checks this too.

### 5.6 What the attacker needs
- A signing operation using CRT with **no output check**.
- A way to cause **one** wrong branch (hardware glitching, Rowhammer-style memory faults, buggy code).
- Knowledge of the encoded message $m$ (deterministic padding, or raw RSA as in this demo). Randomized padding such as PSS
  hides the encoded value from the attacker and blocks this form, but the robust fix is still the output check in Section 6.

---

## 6. The defense (Act 4)

**Verify before you release.** After computing $s$, check

$$s^{\,e} \bmod N \overset{?}{=} m$$

and refuse to output anything if it fails. A faulty signature cannot pass this check (Section 5.3 shows it is wrong
mod $p$), so the attacker never sees $s'$ and has nothing to feed into a gcd. The demo shows honest signing still works
and the glitched run returns `None`.

Other practical defenses: never log or export $d_p$, $d_q$, $q_{\text{inv}}$; use constant-time, fault-hardened libraries;
use randomized padding (RSA-PSS).

---

## 7. How to run it, and what you should see

```bash
python3 rsa_crt_demo.py
```
Colab / Jupyter: paste the script into a cell, or `%run rsa_crt_demo.py --seed 7`. Extra notebook flags are ignored
and `sys.exit` is skipped there.

The demo has four acts:

| Act | What happens | Look for |
|---|---|---|
| 1 | Generates a toy key and prints it | Secrets are printed **only so you can verify**; the attacker never sees them |
| 2 | Attack A | `found after … guesses`, true $x$ matches recovered $x$, recovered $d$ equals real $d$ |
| 3 | Attack B | `v mod q == 0 True`, `v mod p == 0 False`, `gcd … is q? True` |
| 4 | Defense | Honest signature released, faulty one `None` |

Expected output for the default seed (abridged):

```text
ACT 2 -- Attack A: dp leaks, except its lowest bits
  found after 232,818 guesses (0.14s):  k = 14552,  x = 1
  true missing bits were            x = 1   -> OK
  recovered prime factor            12916442222160946727   -> OK
  equals the real d?                OK

ACT 3 -- Attack B: one faulty signature
  v mod q == 0 ?  True   (s' is right modulo q)
  v mod p == 0 ?  False  (s' is wrong modulo p)
  gcd(v, N) = 16177002392964696539   is q? True
  Bonus: gcd(s - s', N) = ...        is q? True

ACT 4 -- Defense: verify before you release
  with fault: output  -> None   (suppressed)

SCOREBOARD: all 8 checks [v]      ALL CHECKS PASSED
```
(Timings vary by machine.)

---

## 8. Real vs. toy

| Aspect | This demo | Real world |
|---|---|---|
| Prime size | 64 bits ($N$ is 128 bits) | 1024+ bits each ($N$ is 2048+ bits) |
| Public exponent | $65537$ | usually $65537$ (Attack A's cost depends on $e$, not on the key size) |
| Randomness | seeded `random` (reproducible) | cryptographic RNG |
| Padding | none (raw $m^d$) | PKCS#1 v1.5 / PSS / OAEP |
| Fault | bit flip we inject | physical glitch or software bug |

Attack B's one-gcd mechanism is **independent of key size**; it is why real libraries verify CRT signatures.
Attack A's cost grows with the number of *missing* bits, not with the size of $p$.

---

## 9. FAQ

**Does the demo "know" the secrets?** It prints them so you can compare. The attack code (`attack_a`, and the gcd lines
in Act 3) only uses public values plus the stated leak or faulty signature.

**Why is `--missing-bits 6` slower?** The search space is $(e-1)\cdot 2^{b}$ pairs, so each extra bit doubles it. Going from
$b=4$ to $b=6$ is about $4\times$ the worst-case work (on the default seed it took about $0.6$ s vs $0.14$ s).

**Could a wrong guess give a false factor?** Practically no. A wrong pair gives a random 128-bit number, and the chance it shares a
factor with $N$ is about $2^{-63}$. The code also requires $1<g<N$, and the demo double-checks the result by rebuilding $d$.

**Why is the fault "bit flip" fine?** Any wrong value of $s_p$ works, as long as $s_p' \not\equiv s_p \pmod p$.

**Is Attack A realistic with $b=4$?** It is a teaching setting. Real partial-exposure attacks use lattice methods to
handle far more missing bits; this demo isolates the algebraic identity they build on.

**Why did my earlier version fail on some keys?** An older estimate of $k$ from $\sqrt N$ assumed $p \approx q$ and could be off
by tens of thousands. This version searches the full range $1 \le k < e$, which is provably complete (Section 4.2).

---

## 10. Code map

| Function | Purpose |
|---|---|
| `make_key()` | Generates $p \neq q$ and retries until $\gcd(e,\varphi)=1$ (fixes a rare crash in the old script) |
| `crt_sign()` | CRT signing; `fault_in_p_branch=True` corrupts $s_p$ |
| `verified_crt_sign()` | **The defense**: checks $s^e \equiv m$ before releasing |
| `attack_a()` | The $k$/$x$ search with the gcd test |
| `d_from_factor()` | From one prime, rebuild $d$ (proves total compromise) |
| `main()` | Four acts plus the scoreboard |

## 11. What changed from the earlier version

- Added the missing step: the script now **actually rebuilds $d$** after each attack, so "total compromise" is demonstrated, not just claimed.
- Added the **defense act** and the no-$m$ bonus gcd.
- Safe key generation (distinct primes, $\gcd(e,\varphi)=1$), seeded RNG, `main()` guard, notebook-safe arguments.
- Removed unused variables and the dangling README reference; removed misleading "smuggled trapdoor" framing.
