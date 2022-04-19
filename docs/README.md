# Documentation

## BLS12-381

We refer to [BLS12-381 For The Rest Of Us](https://hackmd.io/@benjaminion/bls12-381) for details about the curve BLS12-381.

The curve `E` is

```
y^2 = x^3 + 4
```

over `F_q` for a `381`-bit prime `q`.

The twist `E2` of `E` is the curve

```
y^2 = x^3 + 4(1+u)
```

over `F_{q^2}` (see below).

For general implementations not particular to BLS12-381, we denote the prime by `p` below.

## Input format

### Proper BigInt representation

Almost every template we use contains two parameters `n` and `k` (in practice we use `n = 55, k = 7`), where `n` denotes the register bit size and `k` denotes the number of registers: for a large integer `A` we can represent it as an array `a` of length `k` such that

```
A = 2^0 * a[0] + 2^n * a[1] + ... + 2^{n(k-1)} * a[k-1]
```

where each **register** `0 <= a[i] < 2^n` for all `i`. In other words, this is base `2^n` little-endian notation. We call this a **proper** BigInt representation.

### Signed overflow representation

For the purposes of optimizing constraints, in many intermediate circuits we store a big integer in the form

```
A = 2^0 * a[0] + 2^n * a[1] + ... + 2^{n(k-1)} * a[k-1]
```

where each **register** `a[i]` is _signed_ and allowed to have absolute value possibly larger than `2^n`. Here _signed_ means: for a positive number `0 <= b < 2^252`, we can store `-b` as `BabyJubJubPrime - b > 2^252`, where `BabyJubJubPrime` is the inherent `254`-bit field modulus in circom.

Since we store numbers in this overflow representation, we need to make sure that as we perform other operations (such as multiplying two numbers), the absolute value of each array element stays `< 2^252`. This will involve keeping track of various bounds that we make note of both in this documentation and in the code.

### Fp element

An element in `F_p` uniquely corresponds to an integer `0 <= A < p` and we represent it as such.

While our witness generation always produces `F_p` elements `0 <= A < p`, to save constraints we usually **only constrain `A >= 0` to be a proper BigInt** and further constrain `A < p` only when necessary.

### Fp2 element

We construct `F_{p^2} = F_p[u]/(u^2+1)` (this requires `p = 3 (mod 4)`). An element in `F_{p^2}` has the form

```
A_0 + A_1 u
```

with `A_0, A_1` in `F_p`. We represent this element as a `2 x k` array `a` with `a[0], a[1]` in `F_p` (with the same caveats as above).

### Fp12 element

We construct `F_{p^12} = F_{p^2}[w]/(w^6 - (1+u))` (this requires `1+u` to not be a square or cube in `F_{p^2}`). An element in `F_{p^12}` has the form

```
\sum_{i=0}^5 A_i w^i
```

for `A_i` in `F_{p^2}`. We represent it as a `6 x 2 x k` array `a` where `a[i]` represents `A_i`.

### Point in E(Fp)

A point in `E(Fp)` is a pair `(x, y)` with `x,y` in `F_p` so we represent it as a `2 x k` array `P` where `P[0]` represents `x` and `P[1]` represents `y`.

### Point in E2(Fp2)

A point in `E2(Fp2)` is a pair `(x, y)` with `x,y` in `F_{p^2}` so we represent it as a `2 x 2 x k` array `Q` where `Q[0]` represents `x` and `Q[1]` represents `y`.

## BigInt templates

`bigint.circom` is a fork of [circom-ecdsa](https://github.com/0xPARC/circom-ecdsa/blob/master/circuits/bigint.circom), containing templates for operating on numbers represented in array form (usually because the number may exceed the size of the circom prime `BabyJubJubPrime`). See [Input format](#input-format) for a review of different ways to represent numbers.

Here are the main circuits we use from `bigint.circom`:

### `BigMultShortLong(n, k, m_out)`

Multiplication of two polynomials in one variable `X`.

Inputs:

- `a = a[0] + a[1] * X + ... + a[k-1] * X^{k-1}` a degree `k - 1` polynomial
- `b = b[0] + b[1] * X + ... + b[k-1] * X^{k-1}` a degree `k - 1` polynomial

Output:

- `out = out[0] + out[1] * X + ... + out[2 * k - 2] * X^{2*k - 2}`
- `out = a * b` as polynomials in `X`

Notes:

- Optimization due to xJsnark:
  - witness is calculated by normal polynomial multiplication
  - `out` is contrained by evaluating `out(X) === a(X) * b(X)` at `X = 0, ..., 2*k - 2`
- If all `a[i], b[i]` have absolute values `< B`, then `out[i]` has absolute value `< k * B^2`

For code readability/verification, `m_out` is passed as a parameter with the expected max number of bits in the output registers `out[i]`.

### `BigMultShortLongUnequal(n, ka, kb, m_out)`

This is the same as `BigMultShortLong` except `a` has degree `ka - 1` and `b` has degree `kb - 1`.

- If all `a[i], b[i]` have absolute values `< B`, then `out[i]` has absolute value `< min(ka, kb) * B^2`

### `BigLessThan(n, k)`

Inputs:

- `a,b` in [BigInt format](#bigint-representation)

Output:

- `out = (a < b) ? 1 : 0`

### `BigSub(n, k)`

Inputs:

- BigInts `a,b`
- Assume `a >= b`

Output:

- BigInt `out = a - b`
- `underflow` = how much is borrowed at the highest register of subtraction, only nonzero if `a < b`

### `CheckCarryToZero(n, m, k)`

Input:

- `in = in[0] + in[1] * X + ... + in[k-1] * X^{k-1}` as [signed overflow representation](#signed-overflow-representation)
- Assume each `in[i]` is in range `(-2^{m-1}, 2^{m-1})`

Implements:

- Constrain that `in` evaluated at `X = 2^n` as a big integer equals zero.
- Costs roughly `k * (m - n)` constraints.

### `PrimeReduce(n, k, m, p, m_out)`

This template compresses the array length of a [signed overflow BigInt](#signed-overflow-representation) at the cost of increasing the overflow size.

Input: Let `X = 2^n`.

- `in` is length `k + m` array in signed overflow representation
- `in = in[0] + in[1] * X + ... + in[k+m-1] * X^{k+m-1}`
- Assume each `in[i]` is a signed integer such that `abs(in[i] * 2^n) < 2^252`
- `p` is prime in [proper BigInt format](#proper-bigint-representation) passed as parameter.

Output:

- `out = out[0] + out[1] * X + ... + out[k-1] * X^{k-1}` is BigInt congruent to `in (mod p)`

Implementation:

- For `i >= k`, we precompute `X^i = r[i] mod p`, where `r[i][]` is a proper BigInt `< p`.
- `in[i] * X^i` is replaced by `sum_j in[i] * r[i][j] * X^j`
- If each `in[i]` has absolute value `<B`, then `out[i]` has absolute value `< (m+1) * 2^n * B`.

`m_out` is the expected max number of bits in the output registers.

### `BigMultShortLong2D(n, k, l)`

Polynomial multiplication in 2 variables.

Input:

- `a = sum_{i=0}^{l-1} sum_{j=0}^{k-1} a[i][j] * w^i * X^j`
- `b = sum_{i=0}^{l-1} sum_{j=0}^{k-1} b[i][j] * w^i * X^j`

Output:

- `out = sum_{i=0}^{2*l-2} sum_{j=0}^{2*k-1} out[i][j] * w^i * X^j`
- `out = a * b` as product of polynomials in two variables `w, X`

Implementation:

- Uses same xJsnark optimization as [`BigMultShortLong`](#bigmultshortlongn-k-mout)
- If `a[i][j], b[i][j]` have absolute value `< B`, then `out[i][j]` has absolute value `< l * k * B^2`.

Use case: one variable will end up being `2^n`; the other will be the field extension generator.

## BigInt functions

`bigint_func.circom` contains helper functions for witness generation on [proper](#proper-bigint-representation) BigInts represented in array form. The parameters `n,k` always mean the register bit size and number of registers. Unless otherwise specified, BigInts are always represented in proper format.

### `log_ceil(n)`

Parameters:

- `n`: an integer assumed `< 2^253`.

Outputs:

- Returns $\lceil \log_2 n\rceil$.

### `long_add(n, k, a, b)`

Parameters:

- `a`: BigInt with `k` registers
- `b`: BigInt with `k` registers

Outputs:

- Returns BigInt `out = a + b` with `k + 1` registers.

Assumptions:

- Assumes `k+1 <= 50` (because `var` array sizes are static).

### `long_add_unequal(n, k1, k2, a, b)`

Parameters:

- `a`: BigInt with `k1` registers
- `b`: BigInt with `k2` registers

Outputs:

- Returns BigInt `out = a + b` with `k1 + 1` registers.

Assumptions:

- Assumes `k1 >= k2` and `k1 + 1 <= 50`.

### `long_sub(n, k, a, b)`

Parameters:

- `a`: BigInt with `k` registers
- `b`: BigInt with `k` registers

Outputs:

- Returns BigInt `out = a - b`.

Assumptions:

- Assumes `a >= b` and `k <= 50`.

### `long_scalar_mult(n, k, a, b)`

Parameters:

- `a`: `n`-bit integer.
- `b`: BigInt with `k` registers

Outputs:

- Returns a BigInt `out = a * b` with `k+1` registers.

Assumptions:

- Assumes `k+1 <= 50`.

### `long_div2(n, k, m, a, b)`

Parameters:

- `a`: BigInt with `k + m` registers
- `b`: BigInt with `k` registers

Outputs:

- Returns BigInt `out[0]` with `m + 1` registers and `out[1]` with `k` registers such that `a = out[0] * b + out[1]`.

Assumptions:

- Assumes `k + m <= 50`.
- Assumes `b[k-1]` is nonzero.

### `signed_long_to_short(n, k, a)`

Parameters:

- `a`: BigInt in [signed overflow representation](#signed-overflow-representation) with `k` registers.

Outputs:

- Returns a length `50` array `out` representing `a` in proper BigInt form.
- `out[50] = 0` if `a` is positive, `1` if `a` is negative.

Assumptions:

- Assumes `a[i] < 2^{252}` for all `i` and that the proper representation of `a` uses at most `50` registers. (An `assert` will fail otherwise.)

### `prod(n, k, a, b)`

Parameters:

- `a`: BigInt with `k` registers
- `b`: BigInt with `k` registers

Outputs:

- Returns BigInt `out = a * b` with `2k - 1` registers.

Note that the difference between `prod` and `BigMultShortLong` is that the registers of `out` here are in `[0, 2^n)`. We are no longer doing polynomial multiplication and must "carry" using base `2^n`.

Assumptions:

- Assumes `k <= 25`.

### `prod2D(n, k, l, a, b)`

Parameters:

- `a`: length `l x k` array representing a degree `l - 1` polynomial `a = a[0] + a[1] * w + ... + a[l-1] * w^{l-1}` with `a[i]` BigInts with `k` registers.
- `b`: length `l x k` array representing a degree `l - 1` polynomial `b = b[0] + b[1] * w + ... + b[l-1] * w^{l-1}` with `b[i]` BigInts with `k` registers.

Outputs:

- Returns length `2l - 1 x 2k` array representing `out = a * b` as product of polynomials in `w`.
- `out = out[0] + out[1] * w + ... + out[2l-2] * w^{2l-2}`

Note that the difference between `prod2D` and `BigMultShortLong2D` is that the registers of `out` here are in `[0, 2^n)`. Despite the name, `prod2D` is multiplication of polynomials in **one** variable, whereas `BigMultShortLong2D` is multiplication of polynomials in **two** variables.

Assumptions:

- `k <= 25`
- `l <= 10`

### `long_add_mod(n, k, a, b, p)`

### `long_sub_mod(n, k, a, b, p)`

### `prod_mod(n, k, a, b, p)`

These are the same as `long_add, long_sub, prod` except we return `out % p` as BigInts. Note that `long_sub_mod` assumes that `a,b < p`.

### `mod_exp(n, k, a, p, e)`

Parameters:

- `a`: BigInt with `k` registers
- `p`: BigInt with `k` registers (note `p` doesn't need to be prime!)
- `e`: BigInt with `k` registers

Outputs:

- Returns BigInt `out = a^e (mod p)`.

Assumptions:

- Assumes `k <= 50` and `k * n <= 500`.

Implementation:

- A small optimization is to compute the bit length of `e` first, and then do the square-and-multiply method as a loop on the bit length.

### `mod_inv(n, k, a, p)`

Parameters:

- `a`: BigInt with `k` registers
- `p`: a **prime** BigInt with `k` registers

Outputs:

- Returns BigInt with `k` registers `out = a^{-1} (mod p)`. Returns `0` if `a = 0`.

Assumptions:

- Assumes `k <= 50` and `k * n <= 500`.

Implementation:

- Uses Fermat's Little Theorem and computes `out = a^{p-2} (mod p)`.
