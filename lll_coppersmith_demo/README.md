# LLL and Coppersmith, explained: recovering a small secret from a modular equation

A single-file Python demo that recovers a hidden **40-bit** number in a fraction of a second, where brute force
would take about two weeks. It does this by turning a hard problem ("find a root *modulo N*") into an easy one
("find a short vector in a lattice").

```bash
pip install sympy                                       # only needed for the very last step (Colab has it)
python3 lll_coppersmith_demo.py                         # narrated walk-through (about 1 second)
python3 lll_coppersmith_demo.py --unknown-bits 54       # watch it FAIL on purpose
python3 lll_coppersmith_demo.py --sweep                 # table: where does it stop working? (about 7 seconds)
python3 lll_coppersmith_demo.py --seed 7                # a different random instance
```

**How to read the math.** Formulas use a deliberately plain subset of LaTeX so they render everywhere.
A multi-letter name such as $\mathrm{f2}$ or $\mathrm{h3}$ is **one single quantity** (a coefficient), not a product of letters.
Each name is defined in the notation table (Section 2), and the code uses the same names in plain text (`f2`, `h3`).

---

## 1. The whole idea on one page

```text
   The problem:   m = M + secret,   we know c = m^3 mod N,   M is known,   secret is small but unknown

   STEP 1   write the problem as a polynomial f with   f(secret) = 0   (mod N)
              |
   STEP 2   build a LATTICE whose vectors are polynomials that ALL vanish mod N at the secret
              |
   STEP 3   run LLL  ->  a SHORT vector   =   a polynomial h with very small coefficients
              |
   STEP 4   small  =>  h(secret) is not just 0 mod N, it is EXACTLY 0 over the integers
              |        (this is the Howgrave-Graham lemma, the crux of the whole method)
              |
   STEP 5   solve h(x) = 0 with ordinary algebra  ->  the secret
```

**The crux, in two sentences.** A statement "this number is a multiple of $N$" only tells you something
if the number is *smaller than* $N$, because then the multiple must be zero. LLL is the tool that finds a polynomial small
enough for that to happen.

---

## 2. The scenario and notation

**Scenario (a "stereotyped message" attack).** Someone encrypts a message with RSA using the small exponent $e = 3$.
The message has a **known** part (a fixed header or template) and a short **unknown** part. The attacker sees the ciphertext
and wants the unknown part.

| Symbol | Meaning | Known to attacker? |
|---|---|---|
| $N$ | RSA modulus (about 300 bits in the demo) | yes |
| $e = 3$ | public exponent | yes |
| $c$ | ciphertext, $c = m^3 \bmod N$ | yes |
| $m$ | the full message | no |
| $M$ | known part of the message (its low bits are zero) | yes |
| $\hat{x}$ ("x-hat") | the **secret**: the unknown low bits, so $m = M + \hat{x}$ | **no** |
| $X$ | public bound on the secret: $0 \le \hat{x} \le X - 1$, here $X = 2^{40}$ | yes |
| $x$ | the polynomial variable (a placeholder, not the secret) | n/a |
| $f(x)$ | the polynomial $(M+x)^3 - c$ | yes |
| $\mathrm{f2},\ \mathrm{f1},\ \mathrm{f0}$ | coefficients of $f$ (of $x^2$, $x$, and $1$), reduced mod $N$ | yes |
| $h(x)$ | the small polynomial LLL finds | yes |
| $\mathrm{h3},\ \mathrm{h2},\ \mathrm{h1},\ \mathrm{h0}$ | coefficients of $h$ (of $x^3$, $x^2$, $x$, and $1$) | yes |
| $\delta$ | LLL quality parameter, the demo uses $0.99$ | n/a |

In the code: `inst.M`, `inst.X`, `inst.secret` (printed only so you can check the answer), `f2, f1, f0`, `h[3], h[2], h[1], h[0]`.

**Why not just brute-force?** There are $2^{40} \approx 1.1 \times 10^{12}$ candidates. On the test machine one candidate takes
about a microsecond, so the worst case is about 330 hours (roughly two weeks). The demo measures your machine's real rate and prints it.

---

## 3. Step by step (the script prints exactly these steps)

### STEP 0: the scenario
The script builds two 150-bit primes, so $N$ is about 300 bits, picks a random known $M$ and random secret $\hat{x}$, and computes
$c = (M + \hat{x})^3 \bmod N$. From here on the attacker's knowledge is only $N,\ e,\ c,\ M,\ X$.

### STEP 1: the secret is a root of a polynomial modulo N
We know $(M + \hat{x})^3 \equiv c \pmod N$. Define

$$f(x) = (M + x)^3 - c = x^3 + 3M \cdot x^2 + 3M^2 \cdot x + (M^3 - c)$$

so $f(\hat{x}) \equiv 0 \pmod N$. In the code the three lower coefficients are reduced mod $N$:

$$f(x) = x^3 + \mathrm{f2} \cdot x^2 + \mathrm{f1} \cdot x + \mathrm{f0} \pmod N$$

Two facts to notice:

* **Why this is hard.** The *true* integer $f(\hat{x})$ could be $0$, $N$, $-7N$, and so on. All we know is that it is a multiple of $N$.
  The usual formulas for solving cubics need an exact equation, not a "mod N" one.
* **Why reducing the coefficients mod N is allowed.** It changes $f$ only by multiples of $N$, so $f(\hat{x}) \equiv 0 \pmod N$ still holds.
  It keeps the numbers small, which matters in Step 2.
* **$f$ is monic** (its $x^3$ coefficient is exactly $1$). That makes the last lattice entry in Step 2 equal to $X^3$ rather than something huge.

### STEP 2: build a lattice of polynomials that vanish mod N

**What is a lattice?** Take some vectors (the rows of a matrix). The lattice is the set of *all integer combinations* of them.
The same lattice can be described by a long, skewed set of rows or by a short, tidy set of rows. LLL finds the tidy one (Step 3).

**Which polynomials?** Four of them, and the key property is that every one of them is $\equiv 0 \pmod N$ at $x = \hat{x}$:

| Row | Polynomial | Why it vanishes mod N at the secret |
|---|---|---|
| 1 | $N$ | it is a multiple of $N$ |
| 2 | $N \cdot x$ | a multiple of $N$ times anything |
| 3 | $N \cdot x^2$ | same |
| 4 | $f(x)$ | by Step 1 |

Any **integer combination** of these four therefore also vanishes mod $N$ at $\hat{x}$. So *every* vector in this lattice is a polynomial
with that property. We are free to pick the one we like best.

**Turning polynomials into vectors, and the scaling trick.** Write each polynomial as its coefficient vector $[x^0,\ x^1,\ x^2,\ x^3]$,
then multiply column $j$ by $X^j$. The matrix (a code block, because plain text is the clearest way to draw it):

```text
      [  N      0       0        0   ]      <- N
      [  0     N*X      0        0   ]      <- N * x
      [  0      0     N*X^2      0   ]      <- N * x^2
      [  f0   f1*X   f2*X^2     X^3  ]      <- f(x)
```

*Why scale by $X^j$?* The secret is at most $X$, so the term (coefficient of $x^j$) times $\hat{x}^j$ is at most the coefficient times $X^j$ in size.
After scaling, the **length of the vector measures how big the polynomial can possibly be at the secret.** That is exactly what Step 4 needs.

**The determinant decides everything.** The matrix is triangular, so its determinant is the product of the diagonal:

$$\det = N \cdot (N X) \cdot (N X^2) \cdot X^3 = N^3 \cdot X^6$$

A small determinant means short vectors must exist. In the demo run: $\log_2 \det = 1139.5$ bits.

### STEP 3: LLL finds a short vector

**What LLL does, by analogy.** You already know a one-dimensional version: **Euclid's algorithm** repeatedly subtracts multiples to
shrink numbers. LLL is the many-dimensional cousin: it subtracts integer multiples of some rows from others (to make rows shorter and less
"aligned") and swaps rows when that helps, until no more progress is possible.

**A two-dimensional picture.** Take the rows $(1, 0)$ and $(99, 1)$. They generate the same lattice as $(1, 0)$ and $(0, 1)$, because
$(99,1) - 99 \cdot (1,0) = (0,1)$. The first description is long and skewed, the second is short and tidy. LLL does this in 4 dimensions with enormous numbers.

**The guarantee.** For a reduction parameter $\delta$ and a lattice of dimension $n$, the first row $\mathbf{b}$ of the reduced basis satisfies

$$\lVert \mathbf{b} \rVert \le \alpha^{(n-1)/4} \cdot (\det)^{1/n} \qquad \text{where } \alpha = \frac{1}{\delta - 1/4}$$

With $\delta = 0.99$ we get $\alpha \approx 1.35$, and for $n = 4$ the factor $\alpha^{3/4} \approx 1.25$.
In the demo run: $(\det)^{1/4} = 2^{284.9}$, so the guarantee says the shortest vector has length at most about $2^{285.2}$.
LLL actually found $2^{280.0}$, slightly better than promised (this is typical).

**What the script shows.** Row lengths in bits (that is, $\log_2$ of the length):

| | row 1 | row 2 | row 3 | row 4 |
|---|---|---|---|---|
| before LLL | 299.8 | 339.8 | 379.8 | 377.9 |
| after LLL | 280.0 | 280.0 | 280.3 | 299.5 |

### STEP 4: small implies exact (the Howgrave-Graham lemma)

The shortest vector, un-scaled (divide column $j$ by $X^j$), is a polynomial

$$h(x) = \mathrm{h3} \cdot x^3 + \mathrm{h2} \cdot x^2 + \mathrm{h1} \cdot x + \mathrm{h0}$$

By construction it is a combination of the four rows, so $h(\hat{x}) \equiv 0 \pmod N$. The **lemma** says more:

> If the length of the scaled vector is smaller than $N/\sqrt{4} = N/2$ (the 4 is the number of coefficients), then $h(\hat{x}) = 0$ **exactly**, as integers.

**Proof (three lines).**

1. $N$ divides $h(\hat{x})$, because $h$ is a combination of polynomials that vanish mod $N$ at $\hat{x}$.
2. Since $0 \le \hat{x} \le X$, we have $\lvert h(\hat{x}) \rvert \le \lvert \mathrm{h0} \rvert + \lvert \mathrm{h1} \rvert X + \lvert \mathrm{h2} \rvert X^2 + \lvert \mathrm{h3} \rvert X^3$.
   That sum is the sum of the absolute values of the entries of our scaled vector. By Cauchy-Schwarz it is at most $\sqrt{4}$ times the vector's length.
   So $\lvert h(\hat{x}) \rvert \le 2 \cdot \text{length} \lt N$.
3. A multiple of $N$ whose absolute value is below $N$ must be $0$. $\blacksquare$

**What the script shows (default run).**

| Quantity | Value |
|---|---|
| length of shortest vector | $2^{280.0}$ |
| required: smaller than $N/2$ | $2^{298.8}$ |
| margin | **+18.9 bits, condition holds** |

So in the demo, $h(\hat{x}) = 0$ is a plain integer equation, with no "mod" anywhere.

### STEP 5: solve h(x) = 0 and read off the secret
Now it is ordinary algebra: find the integer roots of $h$ in $[0, X-1]$. The script uses `sympy` for this (the only place it is used).
It finds exactly one root, which equals the secret. It then checks $(M + \hat{x})^3 \bmod N = c$.

---

## 4. Why it works only up to a certain size, and then stops

Everything hinges on the Step 4 condition. Put in the LLL guarantee and the determinant $N^3 X^6$:

$$1.25 \cdot (N^3 X^6)^{1/4} \lt \frac{N}{2} \quad \Longrightarrow \quad X^{3/2} \lt \frac{N^{1/4}}{2.5} \quad \Longrightarrow \quad X \lt \frac{N^{1/6}}{1.84}$$

For a 299-bit $N$ this is about $2^{48.9}$, so **the method provably works up to about 49 unknown bits** (the script prints this number).

**Measured** (6 random instances per row, from `--sweep`):

| Unknown bits | Recovered | Exactness test passes (margin above 0) | Average margin (bits) |
|---|---|---|---|
| 38 | 6 / 6 | 6 / 6 | +23.2 |
| 42 | 6 / 6 | 6 / 6 | +15.5 |
| 46 | 6 / 6 | 6 / 6 | +7.2 |
| 48 | 6 / 6 | 6 / 6 | +3.2 |
| 50 | 2 / 6 | 1 / 6 | -0.6 |
| 52 | 0 / 6 | 0 / 6 | -1.0 |
| 54 | 0 / 6 | 0 / 6 | -1.0 |

**What failure looks like** (`--unknown-bits 54`). The margin drops by roughly 2 bits for every extra unknown bit and goes negative around 50 bits.
Past that, LLL can do no better than the trivial row $h(x) = N$.
That polynomial is just the constant $N$ (margin $-1.0$ bits, because $N$ is exactly twice as long as the threshold $N/2$), and it has no useful root. The script prints this explicitly.

**Why 50 bits sometimes still works.** The Step 4 test is *sufficient, not necessary*: LLL's guarantee is worst-case, so right at the edge
it sometimes finds a vector that is a little too long for the proof yet still has the secret as an exact root.

**What about Coppersmith's famous $N^{1/3}$?** Coppersmith's theorem promises roots up to about $N^{1/3}$ for a degree-3 polynomial
(for this $N$, about 100 bits). Reaching it needs a **bigger lattice** (extra rows built from powers of $f$ and shifts by $x$, and a larger dimension).
This demo uses the smallest lattice that shows the idea, which tops out near $N^{1/6}$, so it is a stripped-down teaching version of the theorem, not the full strength one.

---

## 5. Real world vs. this toy

| Aspect | This demo | Real world |
|---|---|---|
| Modulus | about 300 bits | 2048+ bits |
| Exponent | $e = 3$ | usually $65537$ (and then the attack above does not apply as stated) |
| Message | raw $m = M + \hat{x}$, no padding | randomized padding (OAEP), which destroys the "known prefix plus small unknown" structure |
| Lattice | 4 rows (minimal) | larger lattices to reach the $N^{1/3}$ bound |
| LLL | exact fractions, slow but clear | optimized libraries (fplll, Sage) |

**The practical lesson.** "Partially known plaintext" attacks with small $e$ are real (Håstad, Coppersmith, 1990s). They are why modern RSA uses
random padding and larger $e$.

---

## 6. FAQ

**Why do all four polynomials need to vanish mod N?** Because the lattice is built from *integer combinations*. If every building block vanishes mod $N$ at $\hat{x}$, every combination does,
including the short one LLL finds. That is what lets the short vector still "know" the secret.

**Why is the shortest vector not simply the row $N$?** Its length is $2^{299.8}$, while LLL finds vectors of length $2^{280}$. Small unknowns make the lattice's determinant small enough
for much shorter vectors to exist. When the unknown gets too big, $N$ really is the shortest, and the attack fails (Section 4).

**What if h has more than one integer root in range?** An attacker would simply test each candidate against $c$ (is $(M + \text{candidate})^3 \bmod N$ equal to $c$?). With these sizes there is exactly one root.

**Why does column scaling by $X^j$ not change which polynomial we get?** We divide it back out in Step 4. The scaling only affects what "short" means during LLL.

**Is the result random?** The instance is random (use `--seed`). The method is deterministic: the same instance always gives the same answer.

**Why is `sympy` needed at all?** Only for Step 5 (exact integer roots of a cubic). The idea of the demo is Steps 1 to 4. Sympy is preinstalled in Colab.

**Does it run in Colab / Jupyter?** Yes: paste the script into a cell, or `%run lll_coppersmith_demo.py --seed 7`. Extra notebook flags are ignored and `sys.exit` is skipped there.

---

## 7. Expected output (default seed, abridged)

```text
STEP 0 -- The scenario
  UNKNOWN secret: 0 <= secret < 2^40 = 1,099,511,627,776
  Brute force: ... ~931,220 tries/s  -> worst case ~328.0 hours (~13.7 days)
  provable limit ~ 49.1 bits   (you chose 40)

STEP 3 -- Run LLL
  LLL finished in 0.171 s
  shortest vector: 280.0 bits   (was >= 299.8 bits)

STEP 4 -- Is the short vector small enough to be an EXACT equation?
     margin      = +18.9 bits   -> HOLDS

STEP 5 -- Solve h(x) = 0 over the integers
  integer roots of h(x) in [0, 2^40): [660408492882]
  SECRET RECOVERED. No brute force: one 4x4 lattice reduction.
```
(Timings vary by machine. The tries per second figure is measured live.)

---

## 8. Code map

| Function | Step | Purpose |
|---|---|---|
| `make_instance()` | 0 | toy RSA key, known part $M$, secret, ciphertext |
| `polynomial_coefficients()` | 1 | $\mathrm{f2},\ \mathrm{f1},\ \mathrm{f0}$ |
| `build_lattice()` | 2 | the 4x4 matrix with $X^j$ column scaling |
| `lll_reduce()` | 3 | textbook LLL with exact fractions |
| `vector_to_polynomial()`, `howgrave_graham_holds()` | 4 | un-scale, then the size test |
| `integer_roots_in_range()` | 5 | exact integer roots via sympy |
| `provable_bound_bits()` | 4 (Section 4 here) | the $N^{1/6}$ limit |
| `solve_instance()`, `run_sweep()` | all | one-call version used by `--sweep` |

## 9. Glossary

* **Lattice:** all integer combinations of a set of vectors (rows).
* **Basis:** a particular set of rows generating a lattice. Many bases, one lattice.
* **LLL:** an algorithm that, given a basis, returns a nicer one with short, nearly perpendicular rows.
* **Determinant (of the lattice):** the volume of one cell; for a triangular matrix, the product of the diagonal.
* **Monic polynomial:** leading coefficient is $1$.
* **Howgrave-Graham lemma:** a polynomial that vanishes mod $N$ at a point and is small enough must vanish exactly there.
* **Stereotyped message:** a message with mostly known, fixed content and a small unknown part.
