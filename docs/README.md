# Documentation

## Table of Contents

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [BLS12-381](#bls12-381)
- [Input format](#input-format)
  - [Proper BigInt representation](#proper-bigint-representation)
  - [Signed overflow representation](#signed-overflow-representation)
  - [Fp element](#fp-element)
  - [Fp2 element](#fp2-element)
  - [Fp12 element](#fp12-element)
  - [Point in E(Fp)](#point-in-efp)
  - [Point in E2(Fp2)](#point-in-e2fp2)
  - [Point in E(Fp12)](#point-in-efp12)
- [BigInt templates](#bigint-templates)
  - [`BigMultShortLong(n, k, m_out)`](#bigmultshortlongn-k-m_out)
  - [`BigMultShortLongUnequal(n, ka, kb, m_out)`](#bigmultshortlongunequaln-ka-kb-m_out)
  - [`BigLessThan(n, k)`](#biglessthann-k)
  - [`BigSub(n, k)`](#bigsubn-k)
  - [`CheckCarryToZero(n, m, k)`](#checkcarrytozeron-m-k)
  - [`PrimeReduce(n, k, m, p, m_out)`](#primereducen-k-m-p-m_out)
  - [`BigMultShortLong2D(n, k, l)`](#bigmultshortlong2dn-k-l)
- [BigInt functions](#bigint-functions)
  - [`log_ceil(n)`](#log_ceiln)
  - [`long_add(n, k, a, b)`](#long_addn-k-a-b)
  - [`long_add_unequal(n, k1, k2, a, b)`](#long_add_unequaln-k1-k2-a-b)
  - [`long_sub(n, k, a, b)`](#long_subn-k-a-b)
  - [`long_scalar_mult(n, k, a, b)`](#long_scalar_multn-k-a-b)
  - [`long_div2(n, k, m, a, b)`](#long_div2n-k-m-a-b)
  - [`signed_long_to_short(n, k, a)`](#signed_long_to_shortn-k-a)
  - [`prod(n, k, a, b)`](#prodn-k-a-b)
  - [`prod2D(n, k, l, a, b)`](#prod2dn-k-l-a-b)
  - [`long_add_mod(n, k, a, b, p)`](#long_add_modn-k-a-b-p)
  - [`long_sub_mod(n, k, a, b, p)`](#long_sub_modn-k-a-b-p)
  - [`prod_mod(n, k, a, b, p)`](#prod_modn-k-a-b-p)
  - [`mod_exp(n, k, a, p, e)`](#mod_expn-k-a-p-e)
  - [`mod_inv(n, k, a, p)`](#mod_invn-k-a-p)
- [Fp templates](#fp-templates)
  - [`FpNegate(n, k, p)`](#fpnegaten-k-p)
  - [`CheckCarryModP(n, k, m, overflow, p)`](#checkcarrymodpn-k-m-overflow-p)
  - [`SignedFpCarryModP(n, k, overflow, p)`](#signedfpcarrymodpn-k-overflow-p)
  - [`FpMultiply(n, k, p)`](#fpmultiplyn-k-p)
  - [`SignedCheckCarryModToZero(n, k, overflow, p)`](#signedcheckcarrymodtozeron-k-overflow-p)
  - [`FpSgn0(n, k, p)`](#fpsgn0n-k-p)
  - [`FpIsZero(n, k, p)`](#fpiszeron-k-p)
  - [`FpIsEqual(n, k, p)`](#fpisequaln-k-p)
- [Fp2 templates](#fp2-templates)
  - [`SignedFp2MultiplyNoCarryUnequal(n, ka, kb, m_out)`](#signedfp2multiplynocarryunequaln-ka-kb-m_out)
  - [`SignedFp2MultiplyNoCarry(n, k, m_out)`](#signedfp2multiplynocarryn-k-m_out)
  - [`Fp2Compress(n, k, m, p, m_out)`](#fp2compressn-k-m-p-m_out)
  - [`SignedFp2MultiplyNoCarryCompress(n, k, p, m_in, m_out)`](#signedfp2multiplynocarrycompressn-k-p-m_in-m_out)
  - [`SignedFp2MultiplyNoCarryCompressThree(n, k, p, m_in, m_out)`](#signedfp2multiplynocarrycompressthreen-k-p-m_in-m_out)
  - [`RangeCheck2D(n, k)`](#rangecheck2dn-k)
  - [`SignedFp2CarryModP(n, k, overflow, p)`](#signedfp2carrymodpn-k-overflow-p)
  - [`Fp2Multiply(n, k, p)`](#fp2multiplyn-k-p)
  - [`Fp2MultiplyThree(n, k, p)`](#fp2multiplythreen-k-p)
  - [`Fp2Negate(n, k, p)`](#fp2negaten-k-p)
  - [`Fp2Invert(n, k, p)`](#fp2invertn-k-p)
  - [`SignedFp2Divide(n, k, overflowa, overflowb, p)`](#signedfp2dividen-k-overflowa-overflowb-p)
  - [`Fp2Conjugate(n, k, p)`](#fp2conjugaten-k-p)
  - [`Fp2FrobeniusMap(n, k, power, p)`](#fp2frobeniusmapn-k-power-p)
  - [`Fp2Sgn0(n, k, p)`](#fp2sgn0n-k-p)
  - [`Fp2IsZero(n, k, p)`](#fp2iszeron-k-p)
  - [`Fp2IsEqual(n, k, p)`](#fp2isequaln-k-p)
- [Fp12 templates](#fp12-templates)
  - [`Fp12FrobeniusMap(n, k, power)`](#fp12frobeniusmapn-k-power)
  - [`SignedFp12ScalarMultiplyNoCarryUnequal(n, ka, kb, m_out)`](#signedfp12scalarmultiplynocarryunequaln-ka-kb-m_out)
  - [`SignedFp12ScalarMultiplyNoCarry(n, k, m_out)`](#signedfp12scalarmultiplynocarryn-k-m_out)
  - [`SignedFp12Fp2MultiplyNoCarryUnequal(n, ka, kb, m_out)`](#signedfp12fp2multiplynocarryunequaln-ka-kb-m_out)
  - [`SignedFp12Fp2MultiplyNoCarry(n, k, m_out)`](#signedfp12fp2multiplynocarryn-k-m_out)
  - [`SignedFp12MultiplyNoCarryUnequal(n, ka, kb, m_out)`](#signedfp12multiplynocarryunequaln-ka-kb-m_out)
  - [`SignedFp12MultiplyNoCarry(n, k, m_out)`](#signedfp12multiplynocarryn-k-m_out)
  - [`Fp12Compress(n, k, m, p, m_out)`](#fp12compressn-k-m-p-m_out)
  - [`SignedFp12MultiplyNoCarryCompress(n, k, p, m_in, m_out)`](#signedfp12multiplynocarrycompressn-k-p-m_in-m_out)
  - [`SignedFp12CarryModP(n, k, overflow, p)`](#signedfp12carrymodpn-k-overflow-p)
  - [`Fp12Multiply(n, k, p)`](#fp12multiplyn-k-p)
  - [`Fp12Invert(n, k, p)`](#fp12invertn-k-p)
- [Finite field functions](#finite-field-functions)
  - [`get_fp_sng0(a)`](#get_fp_sng0a)
  - [`find_Fp_inverse(n, k, num, p)`](#find_fp_inversen-k-num-p)
  - [`get_signed_Fp_carry_witness(n, k, m, a, p)`](#get_signed_fp_carry_witnessn-k-m-a-p)
  - [`get_signed_Fp2_carry_witness(n, k, m, a, p)`](#get_signed_fp2_carry_witnessn-k-m-a-p)
  - [`get_fp2_sgn0(k, a)`](#get_fp2_sgn0k-a)
  - [`find_Fp2_product(n, k, a, b, p)`](#find_fp2_productn-k-a-b-p)
  - [`find_Fp2_sum(n, k, a, b, p)`](#find_fp2_sumn-k-a-b-p)
  - [`find_Fp2_diff(n, k, a, b, p)`](#find_fp2_diffn-k-a-b-p)
  - [`find_Fp2_exp(n, k, a, p, e)`](#find_fp2_expn-k-a-p-e)
  - [`is_equal_Fp2(n, k, a, b)`](#is_equal_fp2n-k-a-b)
  - [`signed_Fp2_mult_w6(k, a, XI0)`](#signed_fp2_mult_w6k-a-xi0)
  - [`find_Fp12_sum(n, k, a, b, p)`](#find_fp12_sumn-k-a-b-p)
  - [`find_Fp12_diff(n, k, a, b, p)`](#find_fp12_diffn-k-a-b-p)
  - [`find_Fp12_product(n, k, a, b, p)`](#find_fp12_productn-k-a-b-p)
  - [`find_Fp2_inverse(n, k, a, p)`](#find_fp2_inversen-k-a-p)
  - [`find_Fp6_inverse(n, k, p, a0, a1, a2)`](#find_fp6_inversen-k-p-a0-a1-a2)
  - [`find_Fp12_inverse(n, k, p, a)`](#find_fp12_inversen-k-p-a)
- [Elliptic curve templates](#elliptic-curve-templates)
  - [`PointOnLine(n, k, p)`](#pointonlinen-k-p)
  - [`PointOnCurve(n, k, a, b, p)`](#pointoncurven-k-a-b-p)
  - [`PointOnTangent(n, k, a, p)`](#pointontangentn-k-a-p)
  - [`EllipticCurveAddUnequal(n, k, p)`](#ellipticcurveaddunequaln-k-p)
  - [`EllipticCurveDouble(n, k, a, b, p)`](#ellipticcurvedoublen-k-a-b-p)
  - [`EllipticCurveAdd(n, k, p)`](#ellipticcurveaddn-k-p)
  - [`EllipticCurveScalarMultiply(n, k, b, x, p)`](#ellipticcurvescalarmultiplyn-k-b-x-p)
  - [`EllipticCurveScalarMultiplyUnequal(n, k, b, x, p)`](#ellipticcurvescalarmultiplyunequaln-k-b-x-p)
  - [`SignedLineFunctionUnequalNoCarry(n, k, m_out)`](#signedlinefunctionunequalnocarryn-k-m_out)
  - [`SignedLineFunctionEqualNoCarry(n, k, m_out)`](#signedlinefunctionequalnocarryn-k-m_out)
  - [`SignedLineFunctionUnequal(n, k, p)`](#signedlinefunctionunequaln-k-p)
  - [`SignedLineFunctionEqual(n, k, p)`](#signedlinefunctionequaln-k-p)
  - [`Fp12MultplyWithLineUnequal(n, k, kg, overflowg, p)`](#fp12multplywithlineunequaln-k-kg-overflowg-p)
  - [`MillerLoop(n, k, b, x, p)`](#millerloopn-k-b-x-p)
  - [`BLSMillerLoop(n, k, p)`](#blsmillerloopn-k-p)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

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

An element in `F_p` uniquely corresponds to an integer `0 <= A < p`.

While our witness generation always produces `F_p` elements `0 <= A < p`, to save constraints we usually **only constrain `A >= 0` to be a proper BigInt** and further constrain `A < p` only when necessary. In other words, unless otherwise specified, we constrain an `F_p` element representing `A` to be any proper BigInt (necessarily `< 2^{n*k}`) congruent to `A` modulo `p`.

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

A point in `E(F_p)` is a pair `(x, y)` with `x,y` in `F_p` so we represent it as a `2 x k` array `P` where `P[0]` represents `x` and `P[1]` represents `y`.

### Point in E2(Fp2)

A point in `E_2(F_{p^2})` is a pair `(x, y)` with `x,y` in `F_{p^2}` so we represent it as a `2 x 2 x k` array `Q` where `Q[0]` represents `x` and `Q[1]` represents `y`.

### Point in E(Fp12)

A point in `E(F_{p^12})` is a pair `(X, Y)` with `X,Y` in `F_{p^12}` so we represent it as a `2 x 6 x 2 x k` array `Q` where `Q[0]` represents `x` and `Q[1]` represents `y`.

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
- Assume each `in[i]` is in range `(-2^{m-1}, 2^{m-1})`.
- Assume `m < 252`.

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

## Fp templates

`fp.circom` contains templates for operating on elements of `F_p`. See [Fp element](#fp-element) for a reminder on how an element of `F_p` is represented. By default, an `F_p` element means a proper BigInt with `k` registers of bit size `n`. Unless explicitly stated, we **do not** constrain elements to be `< p`.

The templates all take a parameter `p`, which is a length `k` array representing `p` in proper BigInt format with bit size `n` registers.

### `FpNegate(n, k, p)`

Input:

- `in`: `F_p` element

Output:

- `out`: `-in (mod p) = p - in` if `in != 0`, else `0`.

Circuit constrains `in <= p `.

### `CheckCarryModP(n, k, m, overflow, p)`

Inputs:

- `in`: signed overflow representation with `k` registers with absolute value `< 2^overflow`
- `X`: signed overflow representation with `m` registers with absolute value `< 2^{overflow - n - log(min(k,m)) - 1}`
- `Y`: `F_p` element

Implements:

- Constrain `in = p * X + Y`. Uses `CheckCarryToZero`.

Assumptions:

- Assume `overflow < 251`
- Assume `overflow - 1 >= n`.

### `SignedFpCarryModP(n, k, overflow, p)`

Inputs:

- `in`: signed overflow representation with `k` registers with absolute value `< 2^overflow`

Outputs:

- `X`: signed overflow representation with `ceil(overflow / n)` registers in `[-2^n, 2^n)`
- `out`: `F_p` element such that `in = p * X + out` as integers.

Implements:

- Witnesses for `X, out` are precomputed. The witness for `out` satisfies `out < p`.
- Constrain `in === p * X + out` using `CheckCarryModP`. Does **not** constrain `out < p`.

Assumptions:

- Input registers have absolute value `< 2^overflow`
- `n + 1 <= overflow < 251`

### `FpMultiply(n, k, p)`

Input:

- `a,b`: `F_p` elements

Output:

- `out`: `F_p` element congruent to `a * b (mod p)`.

Implementation:

- We first prove a computation of `a * b` as an [overflow BigInt representation](#signed-overflow-representation) with `2k-1` registers using `BigMultShortLong`.
- We "compress" the above computation to a congruent representation `mod p` with `k` registers using `PrimeReduce`. This is an optimization for the next step.
- We use `SignedFpCarryModP` to generate a witness `out` and constrain that `out = a * b (mod p)`.

Assumption:

- `k^2 * 2^{3n} < 251`

### `SignedCheckCarryModToZero(n, k, overflow, p)`

Constrains that the input is congruent to `0 (mod p)`. The method is the same as [`SignedFpCarryModP`](#signedfpcarrymodpn-k-overflow-p) except we save constraints by not having to range check `out = 0`.

Input:

- `in`: [signed overflow BigInt representation](#signed-overflow-representation) with `k` registers
- Assume `in[i]` has absolute value `< 2^overflow`

Output:

- `X` is signed overflow BigInt with `ceil(overflow / n)` registers in `[-2^n, 2^n)` such that `in = p * X` as integers.

Assumption:

- `overflow < 251`

### `FpSgn0(n, k, p)`

This is a proof of computation of `sgn0` as described in the [IRTF draft](https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11#section-4.1).

Input:

- `in`: `F_p` element

Output:

- `out = sgn0(in)` is `0` or `1`.

Implementation:

- `sgn0(in) = in % 2`.
- Constrains `in < p` since the formula for `sgn0` assumes this.

### `FpIsZero(n, k, p)`

Input:

- `in`: `F_p` element

Output:

- `out = (in == 0 ? 1 : 0)`

Implementation:

- Constrains `in < p` and then checks if `in = 0` as BigInt.

### `FpIsEqual(n, k, p)`

Input:

- `in[0], in[1]`: two `F_p` elements

Output:

- `out = (in[0] == in[1])`

Implementation:

- Constrains `in[0] < p` and `in[1] < p` and then checks if they are equal as BigInts.

## Fp2 templates

`fp2.circom` contains templates for operating on elements of `F_{p^2}`. See [Fp2 element](#fp2-element) for a reminder on how an element of `F_{p^2}` is represented. Unless explicitly stated, we **do not** constrain the corresponding `F_p` elements to be `< p`.

### `SignedFp2MultiplyNoCarryUnequal(n, ka, kb, m_out)`

Inputs:

- `a = a[0] + a[1] u` where `a[i]` is [signed overflow representation](#signed-overflow-representation) with `ka` registers.
- `b = b[0] + b[1] u` where `b[i]` is signed overflow representation with `kb` registers

Output:

- `out = out[0] + out[1] u` where `out[i]` is signed overflow representation with `ka + kb - 1` registers such that `out = a * b` in `Z[u]/(u^2 + 1)` (`Z` denotes integers).

Implementation:

- Uses `BigMultShortLongUnequal` to multiply `a * b` without carries.
- Uses the identity `u^2 = -1`.
- Despite the naming, does not involve any prime `p`.
- If `a[i][j], b[i][j]` all have absolute value `< B`, then `out[i][j]` has absolute value `< 2k * B^2`.

Assumptions:

- Assume `2 * k * B^2 < 2^252` where `a[i][j], b[i][j] < B` for all `i,j`.

Parameter `m_out` is passed when calling template with the expected max number of bits in the output registers.

### `SignedFp2MultiplyNoCarry(n, k, m_out)`

Calls `SignedFp2MultiplyNoCarryUnequal(n, k, k, m_out)`.

### `Fp2Compress(n, k, m, p, m_out)`

Inputs:

- `in = in[0] + in[1] u` where `in[i]` is signed overflow representation with `k + m` registers.

Output:

- `out = out[0] + out[1] u` where `out[i]` is signed overflow representation with `k` registers
- `out[i] = in[i] (mod p)`

Implementation:

- Calls [`PrimeReduce`](#primereducen-k-m-p-mout) twice.
- If `in[i][j]` all have absolute value `< B`, then `out[i][j]` has absolute value `< (m+1) * 2^n * B`.

Assumptions:

- Assume `(m+1) * 2^n * B < 2^252`

Parameter `m_out` is passed when calling template with the expected max number of bits in the output registers.

### `SignedFp2MultiplyNoCarryCompress(n, k, p, m_in, m_out)`

Inputs:

- `a = a[0] + a[1] u` where `a[i]` is [signed overflow representation](#signed-overflow-representation) with `k` registers.
- `b = b[0] + b[1] u` where `b[i]` is signed overflow representation with `k` registers

Output:

- `out = out[0] + out[1] u` where `out[i]` is signed overflow representation with `k` registers
- `out` is equal to `a * b` as elements in `F_{p^2}`

Implementation:

- Combines `SignedFp2MultiplyNoCarry` and `Fp2Compress`.
- If `a[i][j], b[i][j]` all have absolute value `< B`, then `out[i][j]` has absolute value `< 2k^2 * 2^n * B^2`.

Parameters:

- `m_in`: expected max number of bits in the input registers (so `B = 2^{m_in}`)
- `m_out`: expected max number of bits in the output registers (so `2k^2 * 2^{n + 2*m_in} <= 2^{m_out}`)

### `SignedFp2MultiplyNoCarryCompressThree(n, k, p, m_in, m_out)`

Inputs:

- `a = a[0] + a[1] u` where `a[i]` is [signed overflow representation](#signed-overflow-representation) with `k` registers.
- `b = b[0] + b[1] u` where `b[i]` is signed overflow representation with `k` registers
- `c = c[0] + c[1] u` where `c[i]` is signed overflow representation with `k` registers

Output:

- `out = out[0] + out[1] u` where `out[i]` is signed overflow representation with `k` registers
- `out` is equal to `a * b * c` as elements in `F_{p^2}`

Implementation:

- Calls `SignedFp2MultiplyNoCarry` twice to compute `a * b` and `(a * b) * c` and then applies `Fp2Compress` to `a * b * c`.
- If `a[i][j], b[i][j], c[i][j]` all have absolute value `< B`, then `out[i][j]` has absolute value `< 4k^2 * (2k-1) * 2^n * B^3`.

Parameters:

- `m_in`: expected max number of bits in the input registers (so `B = 2^{m_in}`)
- `m_out`: expected max number of bits in the output registers (so `4k^2 * (2k-1) * 2^{n + 3*m_in} <= 2^{m_out}`)

### `RangeCheck2D(n, k)`

Inputs:

- `in`: length `2 x k` array of integers

Implements:

- Proves that `0 <= in[i][j] < 2^n` for all `i,j`. In other words, `in[0], in[1]` are proper BigInt representations.

### `SignedFp2CarryModP(n, k, overflow, p)`

Input:

- `in = in[0] + in[1] u` where `in[i]` is signed overflow representation with `k` registers

Outputs:

- `X = X[0] + X[1] u` where `X[i]` is signed overflow representation with `ceil(overflow / n)` registers in `[-2^n, 2^n)`
- `out`: `F_{p^2}` element
- `out = out[0] + out[1] u` such that `in[i] = p * X[i] + out[i]` as integers.

Implements:

- Calls [`SignedFpCarryModP`](#signedfpcarrymodpn-k-overflow-p) twice.
- We do not constrain `out[i] < p`.

Assumptions:

- Input registers have absolute value `< 2^overflow`
- `n + 1 <= overflow < 251`

### `Fp2Multiply(n, k, p)`

Inputs:

- `a,b`: `F_{p^2}` elements

Output:

- `out`: `F_{p^2}` element such that `out = a * b` in `F_{p^2}`

Implements:

- Combines `SignedFp2MultiplyNoCarryCompress` and `SignedFp2CarryModP`
- We do not constrain `out[i] < p`.

Assumptions:

- All BigInts are in proper format, i.e., registers in `[0, 2^n)`
- `2k^2 * 2^{3n} < 2^251`

### `Fp2MultiplyThree(n, k, p)`

Inputs:

- `a,b,c`: `F_{p^2}` elements

Output:

- `out`: `F_{p^2}` element such that `out = a * b * c` in `F_{p^2}`

Implements:

- Combines `SignedFp2MultiplyNoCarryCompressThree` and `SignedFp2CarryModP`
- We do not constrain `out[i] < p`.

Assumptions:

- All BigInts are in proper format, i.e., registers in `[0, 2^n)`
- `4k^2 * (2k-1) * 2^{4n} < 2^251`

### `Fp2Negate(n, k, p)`

Input:

- `in`: `F_{p^2}` element

Output:

- `out`: `F_{p^2}` element
- `out = (p - in[0]) + (p - in[1]) u` so `out = -in` in `F_{p^2}`

Implementation:

- Constrains `0 <= in[0], in[1] <= p`
- Calls `FpNegate` twice

### `Fp2Invert(n, k, p)`

Input:

- `in`: `F_{p^2}` element

Output:

- `out`: `F_{p^2}` element equal to `in^{-1}`

Implementation:

- Precompute witness for `in^{-1}` using `find_Fp2_inverse` helper function
- Constrain `out * in == 1` in `F_{p^2}` using `Fp2Multiply`

Assumptions:

- `2k^2 * 2^{3n} < 2^251`

### `SignedFp2Divide(n, k, overflowa, overflowb, p)`

Inputs:

- `a = a[0] + a[1] u` where `a[i]` are [signed overflow representations](#signed-overflow-representation) with `k` registers
- `b = b[0] + b[1] u` where `b[i]` are signed overflow representations with `k` registers

Output:

- `out`: `F_{p^2}` element equal to `a / b` in `F_{p^2}`.

Assumptions:

- Registers of `a` have absolute value `< 2^overflowa`
- Registers of `b` have absolute value `< 2^overflowb`
- `k + ceil(overflowa / n) <= 50`
- `2k + ceil(overflowb / n) <= 50`
- `overflowa < 250`
- `2k^2 * 2^{2n + overflowb} < 2^250`

### `Fp2Conjugate(n, k, p)`

Input:

- `in`: `F_{p^2}` element

Output:

- `out`: `F_{p^2}` element equal to conjugate of `in`
- `out = in[0] - in[1] u`

### `Fp2FrobeniusMap(n, k, power, p)`

Input:

- `in`: `F_{p^2}` element

Output:

- `out`: `F_{p^2}` element equal to `in^{p^power}`

Implements:

- Conjugates `in` if `power` is odd, `out = in` if `power` is even.

Assumption:

- `p = 3 (mod 4)` so that `a - b u = (a + b u)^p`

### `Fp2Sgn0(n, k, p)`

This is a proof of computation of `sgn0` as described in the [IRTF draft](https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11#section-4.1).

Input:

- `in`: `F_{p^2}` element

Output:

- `out = sgn0(in)` is either `0` or `1`

Implementation:

- `sgn0(in) = sgn0(in[0]) || ( in[0] == 0 && sgn0(in[1]) )`
- Calls `FpSgn0` twice, in the process constrains `in[0], in[1] < p`

### `Fp2IsZero(n, k, p)`

Input:

- `in`: `F_{p^2}` element

Output:

- `out = (in == 0)`

Implementation:

- Constrains `in[0], in[1] < p` and then checks if `in[0] == 0 && in[1] == 0`.

### `Fp2IsEqual(n, k, p)`

Input:

- `a, b`: `F_{p^2}` elements

Output:

- `out = (a == b)`

Implementation:

- Constrains `a[i], b[i] < p` and then checks if `a[i] == b[i]` as BigInts for `i=0,1`

## Fp12 templates

`fp12.circom` contains templates for operating on elements of `F_{p^12}`. See [Fp12 element](#fp12-element) for a reminder on how an element of F\_{p^12} is represented. Unless explicitly stated, we do not constrain the corresponding `F_p` elements to be `< p`.

### `Fp12FrobeniusMap(n, k, power)`

This circuit computes `power` applications of the Frobenius map (i.e., raising to `p`-th power). It is specialized to `p = q` the BLS12-381 base field modulus, due to the need to precompute Frobenius coefficients.

Input:

- `in`: `F_{p^12}` element

Output:

- `out` : `F_{p^12}` element equal to `in^{p^power}`

Implementation:

- We follow the methodology from this [paper](https://eprint.iacr.org/2010/354.pdf) and [@noble/bls12-381](https://github.com/paulmillr/noble-bls12-381/blob/main/math.ts)
- We load values for `FP12_FROBENIUS_COEFFICIENTS` precomputed in Python.

Assumption:

- Assume `p = 1 (mod 6)`

### `SignedFp12ScalarMultiplyNoCarryUnequal(n, ka, kb, m_out)`

Inputs:

- `a`: BigInt in [signed overflow representation](#signed-overflow-representation) with `ka` registers.
- `b = sum_{i=0}^5 sum_{j=0}^1 b[i][j] w^i u^j` where `b[i][j]` is signed overflow representation with `kb` registers

Output:

- `out = sum_{i=0}^5 sum_{j=0}^1 out[i][j] w^i u^j` where `out[i][j]` is signed overflow representation with `ka + kb - 1` registers
- `out = a * b` as `12`-dimensional vector with integer coefficients

Implementation:

- Uses `BigMultShortLongUnequal` to multiply `out[i][j] = a * b[i][j]` without carries.
- If input registers all have absolute value `< B`, then `out[i][j]` has absolute value `< min(ka, kb) * B^2`.

Assumption:

- Assume `min(ka, kb) * B^2 < 2^252`.

Parameter `m_out` is passed when calling template with the expected max number of bits in the output registers.

### `SignedFp12ScalarMultiplyNoCarry(n, k, m_out)`

Calls `SignedFp12ScalarMultiplyNoCarryUnequal(n, k, k, m_out)`.

### `SignedFp12Fp2MultiplyNoCarryUnequal(n, ka, kb, m_out)`

Inputs:

- `a = a[0] + a[1] u` where `a[i]` is signed overflow representation with `ka` registers.
- `b = sum_{i=0}^5 sum_{j=0}^1 b[i][j] w^i u^j` where `b[i][j]` is signed overflow representation with `kb` registers

Output:

- `out = sum_{i=0}^5 sum_{j=0}^1 out[i][j] w^i u^j` where `out[i][j]` is signed overflow representation with `ka + kb - 1` registers
- `out = a * b` in `F_{p^12}` as overflow representations

Implementation:

- Multiplies `out[i] = a * b[i]` as in `SignedFp2MultiplyNoCarry` using `u^2 = -1`.
- If input registers all have absolute value `< B`, then `out[i][j]` has absolute value `< 2*min(ka, kb) * B^2`.

Assumption:

- Assume `2*min(ka, kb) * B^2 < 2^252`.

Parameter `m_out` is passed when calling template with the expected max number of bits in the output registers.

### `SignedFp12Fp2MultiplyNoCarry(n, k, m_out)`

Calls `SignedFp12Fp2MultiplyNoCarryUnequal(n, k, k, m_out)`.

### `SignedFp12MultiplyNoCarryUnequal(n, ka, kb, m_out)`

Inputs:

- `a = sum_{i=0}^5 sum_{j=0}^1 a[i][j] w^i u^j` where `a[i][j]` is signed overflow representation with `ka` registers
- `b = sum_{i=0}^5 sum_{j=0}^1 b[i][j] w^i u^j` where `b[i][j]` is signed overflow representation with `kb` registers

Output:

- `out = sum_{i=0}^5 sum_{j=0}^1 out[i][j] w^i u^j` where `out[i][j]` is signed overflow representation with `ka + kb - 1` registers
- `out = a * b` in `F_{p^12}` as overflow representations

Implementation:

- Write `a = a[][0] + a[][1] u` and `b = b[][0] + b[][1] u`.
- Prove computation of `a * b = (a[][0] * b[][0] - a[][1] * b[][1]) + (a[][0] * b[][1] + a[][1] * b[][0]) u` using identity `u^2 = -1`.
  - Prove computation of `a[][i] * b[][j]` as polynomials in two variables `w, X` using [`BigMultShortLong2DUnequal`](#bigmultshortlongunequaln-ka-kb-mout) where `X` is variable for `2^n`.
- Result is polynomial of degree `10` in `w` and degree `ka + kb - 2` in `X`.
- Use identity `w^6 = 1 + u` to turn into degree `5` polynomial in `w`
- If input registers all have absolute value `< B`, then `out[i][j]` has absolute value `< 18*min(ka, kb) * B^2`.

Assumption:

- Assume `18*min(ka, kb) * B^2 < 2^252`.

Parameter `m_out` is passed when calling template with the expected max number of bits in the output registers.

### `SignedFp12MultiplyNoCarry(n, k, m_out)`

Calls `SignedFp12MultiplyNoCarryUnequal(n, k, k, m_out)`

### `Fp12Compress(n, k, m, p, m_out)`

Inputs:

- `in = sum_{i=0}^5 sum_{j=0}^1 in[i][j] w^i u^j` where `in[i][j]` is signed overflow representation with `k + m` registers

Outputs:

- `out = sum_{i=0}^5 sum_{j=0}^1 out[i][j] w^i u^j` where `out[i][j]` is signed overflow representation with `k` registers

Implements:

- Calls [`PrimeReduce`](#primereducen-k-m-p-mout) `12` times.
- If all input registers have absolute value `< B`, then output registers have absolute value `< (m+1) * 2^n * B`.

Assumption:

- `(m+1) * 2^n * B < 2^252`

### `SignedFp12MultiplyNoCarryCompress(n, k, p, m_in, m_out)`

Inputs:

- `a = sum_{i=0}^5 sum_{j=0}^1 a[i][j] w^i u^j` where `a[i][j]` is signed overflow representation with `k` registers
- `b = sum_{i=0}^5 sum_{j=0}^1 b[i][j] w^i u^j` where `b[i][j]` is signed overflow representation with `k` registers

Output:

- `out = sum_{i=0}^5 sum_{j=0}^1 out[i][j] w^i u^j` where `out[i][j]` is signed overflow representation with `k` registers
- `out = a * b` in `F_{p^12}` as overflow representations

Implementation:

- Combines `SignedFp12MultiplyNoCarry` and `Fp12Compress`.
- If all input registers have absolute value `< B`, then output registers have absolute value `< 18k^2 * 2^n * B^2`

Assumption:

- `18k^2 * 2^n * B^2 < 2^252`

### `SignedFp12CarryModP(n, k, overflow, p)`

Inputs:

- `in = sum_{i=0}^5 sum_{j=0}^1 in[i][j] w^i u^j` where `in[i][j]` is signed overflow representation with `k` registers

Outputs:

- `X = sum_{i=0}^5 sum_{j=0}^1 X[i][j] w^i u^j` where `X[i][j]` is signed overflow representation with `ceil(overflow / n)` registers in `[-2^n, 2^n)`
- `out`: proper `F_{p^12}` element
- `in = p * X + out` in `Z[u,w]/(u^2 + 1, w^6 - (1+u))`.

Implements:

- Calls [`SignedFpCarryModP`](#signedfpcarrymodpn-k-overflow-p) `12` times.
- We do not constrain `out[i][j] < p`.

Assumptions:

- Input registers have absolute value `< 2^overflow`
- `n + 1 <= overflow < 251`

### `Fp12Multiply(n, k, p)`

Inputs:

- `a,b`: `F_{p^12}` elements

Output:

- `out`: `F_{p^12}` element such that `out = a * b` in `F_{p^12}`

Implements:

- Combines `SignedFp12MultiplyNoCarryCompress` and `SignedFp12CarryModP`
- We do not constrain `out[i][j] < p`.

Assumptions:

- All BigInts are in proper format, i.e., registers in `[0, 2^n)`
- `18k^2 * 2^{3n} < 2^251`

### `Fp12Invert(n, k, p)`

Input:

- `in`: `F_{p^12}` element

Output:

- `out`: `F_{p^12}` element equal to `in^{-1}`

Implementation:

- Precompute witness for `in^{-1}` using `find_Fp12_inverse` helper function
- Constrain `out * in == 1` in `F_{p^12}` using `Fp12Multiply`

Assumptions:

- `18k^2 * 2^{3n} < 2^251`

## Finite field functions

`field_elements_func.circom` contains helper functions for witness generation involving field operations in `F_p, F_{p^2}, F_{p^12}`. Unless otherwise specified, all representations are in proper BigInt format with `k` registers of bit size `n`, and corresponding `F_p` elements are `< p`.

### `get_fp_sng0(a)`

Parameter:

- `a`: `F_p` element

Output:

- Returns `a % 2`

### `find_Fp_inverse(n, k, num, p)`

This is never used anywhere. Did not empirically lead to faster witness generation times versus Fermat's Little Theorem.

Parameter:

- `num`: `F_p` element

Output:

- Returns `num^{-1}` as `F_p` element using extended Euclidean algorithm.

### `get_signed_Fp_carry_witness(n, k, m, a, p)`

Parameter:

- `a`: BigInt in [signed overflow representation](#signed-overflow-representation) with `k` registers

Output:

- Return `out` as length `2 x 50` array
- `out[0]` has `m` registers in range `[-2^n, 2^n)`
- `out[1]` is `F_p` element with `k` registers, `0 <= out[1] < p`
- `a = p * out[0] + out[1]` so `out[1] = a (mod p)`

Assumptions:

- Assume actual integer value of `a` has absolute value `< 2^{n*(k+m)}`
- `k + m <= 50`

### `get_signed_Fp2_carry_witness(n, k, m, a, p)`

Parameter:

- `a`: length `2 x k` array with `a[i]` in [signed overflow representation](#signed-overflow-representation) with `k` registers

Output:

- Return `out` as length `2 x 2 x 50` array
- `out[i] = a[i] (mod p)` by calling `get_signed_Fp_carry_witness` twice

Assumptions:

- Assume actual integer value of `a[i]` has absolute value `< 2^{n*(k+m)}`
- `k + m <= 50`

### `get_fp2_sgn0(k, a)`

Parameter:

- `a`: `F_{p^2}` element

Output:

- Returns `sgn0(a[0]) || (a[0] == 0 && sgn0(a[1]))`

### `find_Fp2_product(n, k, a, b, p)`

Parameters:

- `a,b`: `F_{p^2}` elements

Output:

- Return `a * b` as proper `F_{p^2}` element

### `find_Fp2_sum(n, k, a, b, p)`

Parameters:

- `a,b`: `F_{p^2}` elements

Output:

- Return `a + b` as proper `F_{p^2}` element

### `find_Fp2_diff(n, k, a, b, p)`

Parameters:

- `a,b`: `F_{p^2}` elements

Output:

- Return `a - b` as proper `F_{p^2}` element

### `find_Fp2_exp(n, k, a, p, e)`

Parameters:

- `a`: `F_{p^2}` element
- `e`: proper BigInt representation with `2k` registers

Output:

- Return `a^e` as proper `F_{p^2}` element

Implementation:

- Small optimization to compute bit length of `e` first before running square-and-multiply loop on the bit length.

Assumption:

- Assume `k * n <= 400`

### `is_equal_Fp2(n, k, a, b)`

Parameters:

- `a,b`: `F_{p^2}` elements

Output:

- Returns `a == b`

### `signed_Fp2_mult_w6(k, a, XI0)`

This is a macro.

Parameter:

- `a = a[0] + a[1] u` with `a[i]` signed overflow representations.

Output:

- Returns `a * (XI0 + u)` in `Z[u]/(u^2 + 1)`.

### `find_Fp12_sum(n, k, a, b, p)`

Parameters:

- `a,b`: `F_{p^12}` elements

Output:

- Returns `a + b` as proper `F_{p^12}` element

### `find_Fp12_diff(n, k, a, b, p)`

Parameters:

- `a,b`: `F_{p^12}` elements

Output:

- Returns `a - b` as proper `F_{p^12}` element

### `find_Fp12_product(n, k, a, b, p)`

Parameters:

- `a,b`: `F_{p^12}` elements

Output:

- Returns `a * b` as proper `F_{p^12}` element

### `find_Fp2_inverse(n, k, a, p)`

Parameters:

- `a`: `F_{p^2}` element

Output:

- Returns `a^{-1}` as proper `F_{p^2}` element
- Uses fact that `(a + b u)^{-1} = (a^2 + b^2)^{-1} * (a - b u)` and reduces to compute an `F_p` inverse.

### `find_Fp6_inverse(n, k, p, a0, a1, a2)`

This is only used as an intermediary to find inverses in `F_{p^12}`.

Parameters:

- `a0, a1, a2`: `F_{p^2}` elements
- Consider `a = a0 + a1 v + a2 v^2` as `F_{p^6}` element, where `v = w^2` and `v^3 = 1 + u`.

Output:

- Returns `a^{-1}` in `F_{p^6}` as a proper `F_{p^12}` representation (length `6 x 2 x k` array)
- Computed by solving a system of equations with `F_{p^2}` coefficients.

### `find_Fp12_inverse(n, k, p, a)`

Parameters:

- `a`: `F_{p^12}` element

Output:

- Return `a^{-1}` as proper `F_{p^12}` element
- Write `a = a0 + a1 w` and use fact that `(a0 + a1 w) * (a0 - a1 w) = a0^2 - a1^2 v` is an element in `F_{p^6}`, where `v = w^2`, to reduce to finding inverses in `F_{p^6}`.

## Elliptic curve templates

`curve.circom` contains circuits for operations on `F_p` points of an elliptic curve. Most of the circuits work for an arbitrary elliptic curve `E` in short Weierstrass form `y^2 = x^3 + a x + b` where `a,b` are "short" integers, i.e., `< 2^n`. For optimization purposes and simplicity, a few circuits only work for curves of the form `y^2 = x^3 + b` (i.e., `a = 0`), but can be easily modified to work in the case `a != 0` if further need arises.

See [Point in `E(F_p)`](#point-in-efp) for a reminder on how we represent points in `E(F_p)`. Unless otherwise specifics, all BigInts use `k` registers of bit size `n`. Unless explicitly stated, we do not constrain a point in `E(F_p)` to actually lie on the curve.

### `PointOnLine(n, k, p)`

Prove that three points `(x_0, y_0), (x_1, y_1), (x_2, -y_2)` are co-linear.

Input:

- `in`: length `3 x 2 x k` array with `in[i][j]` element of `F_p`

Implementation:

- Let `in[i] = (x_i, y_i)`.
- Proves that `(y_0 + y_2) * (x_1 - x_0) - (y_1 - y_0) * (x_0 - x_2) == 0 (mod p)`.
- `(y_0 + y_2) * (x_1 - x_0) - (y_1 - y_0) * (x_0 - x_2)` is computed "without carry" using `BigMultShortLong` and then compressed using `PrimeReduce`.
- Use `SignedCheckCarryModToZero` to prove congruence to `0 (mod p)`.

Assumption:

- Requires `8k^2 * 2^{3n} < 2^251` to prevent overflow

### `PointOnCurve(n, k, a, b, p)`

Prove that a point `(x,y)` lies on an elliptic curve `y^2 = x^3 + a x + b`.

Input:

- `in`: length `2 x k` array representing two `F_p` elements

Implementation:

- Let `in = (x,y)`.
- Proves that `x^3 + a x + b - y^2 = 0 (mod p)`.
- Computes `x^3 + a x + b - y^2` without carries using `BigMultShortLong` and then compresses using `PrimeReduce`.
- Use `SignedCheckCarryModToZero` to prove congruence to `0 (mod p)`.

Assumptions:

- Curve parameters `0 <= a,b < 2^n`
- Requires `(k^2 + 1) * (2k-1) * 2^{4n} < 2^251` to prevent overflow.

### `PointOnTangent(n, k, a, p)`

Prove that a point `(x_1, y_1)` lies on the tangent line to curve `E : y^2 = x^3 + a x + b` at a point `(x_3, y_3)`. Does not check that `(x_3, y_3)` is on the curve. Note that the circuit does not require `b`.

Input:

- `in`: length `2 x 2 x k` array representing four `F_p` elements

Implementation:

- Let `in[0] = (x_1, y_1)` and `in[1] = (x_3, y_3)`.
- Proves that `2y_1 (y_1 + y_3) = (3 x_1^2 + a)(x_1 - x_3)`.

Assumptions:

- Curve parameter `0 <= a < 2^n`
- Require `(3k^2 + 1)*(2k - 1) * 2^{4n} < 2^251` to prevent overflow.
- `5k <= 2^n`

### `EllipticCurveAddUnequal(n, k, p)`

Verify the computation of addition of two points `P1 + P2` on an elliptic curve `E : y^2 = x^3 + a1 x + b1` assuming that `P1 != P2` and `P1 != -P2`. We do not verify these assumptions and do not check that `P1, P2` actually lie on the curve. Note that the circuit does not require knowledge of curve parameters `a1, b1`.

Inputs:

- `a`: point on `E(F_p)`
- `b`: point on `E(F_p)`

Output:

- `out`: point on `E(F_p)` equal to `a + b` using elliptic curve addition.

Implementation:

- Let `a = (x_1, y_1), b = (x_2, y_2), out = (x_3, y_3)`.
- Precomputes witness for `out` and then proves that:
- `(x_1 + x_2 + x_3)*(x_2 - x_1)^2 = (y_2 - y_1)^2`
- `(y_1 + y_3)*(x_2 - x_1) = (y_2 - y_1)*(x_1 - x_3)` using `PointOnLine`

Assumptions:

- `a, b` are on `E(F_p)`
- `x_1 != x_2`, i.e., `a != b` and `a != -b`.
- Require `(3k^2+1)(2k-1) * 2^{4n} < 2^251` to prevent overflow.
- `k <= 2^n`

### `EllipticCurveDouble(n, k, a, b, p)`

Verify the computation of doubling a point on the elliptic curve `E : y^2 = x^3 + a x + b`. We do not check that the point actually lies on the curve.

Inputs:

- `in`: point on `E(F_p)`

Output:

- `out`: point on `E(F_p)` equal to `2 * in` using elliptic curve doubling formula.

Implementation:

- Let `in = (x_1, y_1), out = (x_3, y_3)`.
- Precomputes witness for `out` and then proves that:
- `(x_3, y_3)` lies on `E` using `PointOnCurve`
- `(x_3, y_3)` is on the tangent line to `E` at `(x_1, y_1)` using `PointOnTangent`
- `x_1 != x_3`: the only points of intersection of the tangent line at `(x_1, y_1)` with curve `E` are `(x_1, y_1)` and `(x_3, y_3)`.

Assumptions:

- Curve parameters `0 <= a, b < 2^n`
- `y_1 != 0`, i.e., `in` is not a point of order `2`. (Fact: BLS12-381 has no points of order `2`.)
- Require `(3k^2 + 1)*(2k - 1) * 2^{4n} < 2^251` to prevent overflow.
- `5k <= 2^n`

### `EllipticCurveAdd(n, k, p)`

Verify the computation of addition of two points `P1 + P2` on an elliptic curve `E : y^2 = x^3 + a1 x + b1` without any assumptions on `P1, P2`. Moreover, we allow `P1, P2` to be the point at infinity.

We do not check that `P1, P2` actually lie on the curve. Due to the deterministic nature of `circom`, this circuit costs the sum of the constraints of `EllipticCurveAddUnequal` and `EllipticCurveDouble`.

Inputs:

- `a`: point on `E(F_p)`
- `aIsInfinity`: `1` if we should interpret `a` as point at infinity; `0` otherwise
- `b`: point on `E(F_p)`
- `bIsInfinity`: `1` if we should interpret `b` as point at infinity; `0` otherwise

Output:

- `out`: point on `E(F_p)` equal to `a + b` using elliptic curve addition.
- `isInfinity`: `1` if `a + b` is point at infinity; `0` otherwise
- If `a,b` are on `E(F_p)` (even if `aIsInfinity == 1`), then `out` is on `E(F_p)` (even if `isInfinity == 1`)

Implementation:

- Proves computations of both `EllipticCurveAddUnequal` and `EllipticCurveDouble`
- Case analysis of `a = b`, `a = -b`, `a, b` equal point at infinity

Assumptions:

- `a, b` are on `E(F_p)`
- `E(F_p)` has no points of order `2`
- Require `(3k^2 + 1)*(2k - 1) * 2^{4n} < 2^251` to prevent overflow.
- `5k <= 2^n`

### `EllipticCurveScalarMultiply(n, k, b, x, p)`

Proves computation of scalar multiplication `[x]P` on elliptic curve of form `E : y^2 = x^3 + b` using the double-and-add method.

Parameter:

- `x`: integer in `[0, 2^250)` to scalar multiply by

Inputs:

- `in`: point on `E(F_p)`
- `inIsInfinity`: `1` if `in` should be interpreted as point at infinity; `0` otherwise

Outputs:

- `out`: point on `E(F_p)` equal to `[x]P` where `P = in`
- `isInfinity`: `1` if `out` should be interpreted as point at infinity; `0` otherwise
- Assuming `in` is a point on `E(F_p)` (even when `inIsInfinity == 1`), `out` is a point on `E(F_p)` (even when `isInfinity == 1`).

Implementation:

- Double-and-add method. Uses `EllipticCurveAdd` in the "add" part since we cannot guarantee we are adding two unequal points.

Assumptions:

- `0 <= x < 2^250`
- `in` is point on `E(F_p)` even if `inIsInfinity == 1`
- `E(F_p)` has no points of order `2`

### `EllipticCurveScalarMultiplyUnequal(n, k, b, x, p)`

Proves computation of scalar multiplication `[x]P` on elliptic curve of form `E : y^2 = x^3 + b` using the double-and-add method, assuming that `P` has order `> x` so we never arrive at the point at infinity and we are always adding unequal points.

Parameter:

- `x`: integer in `[0, 2^250)` to scalar multiply by

Inputs:

- `in`: point on `E(F_p)`

Outputs:

- `out`: point on `E(F_p)` equal to `[x]P` where `P = in`.

Implementation:

- Double-and-add method. Uses `EllipticCurveAddUnequal` in the "add" part because if `P` has order `> x`, then `[i]P != [j]P` for any `0 < i,j <= x`.

Assumptions:

- `0 <= x < 2^250`
- `in` is point on `E(F_p)` of order `> x`
- `E(F_p)` has no points of order `2`

### `SignedLineFunctionUnequalNoCarry(n, k, m_out)`

Proves evaluation of a line through two **unequal** `E(F_p)` points at a `E(F_{p^12})` point, evaluated without carries.

Inputs:

- `P`: two `E(F_p)` points `P[0], P[1]`
- `Q`: a [point in `E(F_{p^12})`](#point-in-efp12)

Output:

- `out`: an element in `F_{p^12}` with `out[i][j]` in [signed overflow representation](#signed-overflow-representation) with `2k-1` registers

Implementation:

- Let `P[0] = (x_1, y_1), P[1] = (x_2, y_2)` and `Q = (X, Y)`.
- `out = (y_1 - y_2) X + (x_2 - x_1) Y + (x_1 y_2 - x_2 y_1)`
- We evaluate `out` using `SignedFp12ScalarMultiplyNoCarry`
- If all input registers are in `[0, B)`, then output registers have absolute value `< 3k * B^2`

Assumptions:

- `P[0] != P[1]`
- `3k * B^2 < 2^252`

Parameter `m_out` is the expected max number of bits in the output registers.

### `SignedLineFunctionEqualNoCarry(n, k, m_out)`

Proves evaluation of a tangent line to curve `E : y^2 = x^3 + b` through an `E(F_p)` point `P` at a `E(F_{p^12})` point `Q`, evaluated without carries.

Inputs:

- `P`: `E(F_p)` point
- `Q`: `E(F_{p^12})` point

Output:

- `out`: an element in `F_{p^12}` with `out[i][j]` in signed overflow representation with `3k-2` registers

Implementation:

- Let `P = (x, y)` and `Q = (X, Y)`.
- `out = 3 x^2 (-X + x) + 2 y (Y - y)`
- We evaluate `out` using `SignedFp12ScalarMultiplyNoCarry`
- If all input registers are in `[0, B)`, then output registers have absolute value `< (3k^2 + 2k/B) * B^3` assuming `2k <= B`

Assumptions:

- `(3k^2 + 2k/B) * B^3 < 2^252`

### `SignedLineFunctionUnequal(n, k, p)`

Proves evaluation of a line through two **unequal** `E(F_p)` points at a `E(F_{p^12})` point.

Inputs:

- `P`: two `E(F_p)` points `P[0], P[1]`
- `Q`: `E(F_{p^12})` point

Output:

- `out`: `F_{p^12}` element equal to `l_{P[0], P[1]}(Q)`

Implementation:

- Calls `LineFunctionUnequalNoCarry`, then compresses using `Fp12Compress`, and finally carries using `SignedFp12CarryModP`

Assumptions:

- `P[0] != P[1]`
- Requires `3k^2 * 2^{3n} < 2^251` to prevent overflow.

### `SignedLineFunctionEqual(n, k, p)`

Proves evaluation of a tangent line to curve `E : y^2 = x^3 + b` through an `E(F_p)` point `P` at a `E(F_{p^12})` point `Q`.

Inputs:

- `P`: `E(F_p)` point
- `Q`: `E(F_{p^12})` point

Output:

- `out`: `F_{p^12}` element

Implementation:

- Calls `SignedLineFunctionEqualNoCarry`, then compresses using `Fp12Compress`, and finally carries using `SignedFp12CarryModP`

Assumptions:

- Requires `(3k^2 + 1) * (2k-1) * 2^{4n} < 2^251` to prevent overflow.
- Assumes `2k <= 2^n`

### `Fp12MultplyWithLineUnequal(n, k, kg, overflowg, p)`

Inputs:

- `g`: element of `F_{p^12}` in signed overflow format with `kg` registers
- `P`: two `E(F_p)` points
- `Q`: `E(F_{p^12})` point

Output:

- `out`: `F_{p^12}` element equal to `g * l_{P[0], P[1]}(Q)`

Implementation:

- Calls `SignedLineFunctionUnequalNoCarry` to compute `l_{P[0], P[1]}(Q)`
- Calls `SignedFp12MultiplyNoCarryUnequal` to compute `g * l_{P[0], P[1]}(Q)` **all without carries**
- Compresses using `Fp12Compress` and finally carries **just once** using `SignedFp12CarryModP`

Assumptions:

- `P[0] != P[1]`
- Requires `3k * min(kg, 2k-1) * 18 * (k + kg - 1) * 2^{overflowg + 3n} < 2^251` to prevent overflow. (In practice `kg = k` and `overflowg = n`.)

The point is that `l_{P[0], P[1]}(Q)` is quadratic, so we can squeeze in another multiplication before we need to carry. This is used in `MillerLoop`.

### `MillerLoop(n, k, b, x, p)`

Proves computation of Miller's algorithm on a curve of the form `y^2 = x^3 + b`. Only used for the Tate pairing.

Parameters:

- `b`: curve parameter
- `x`: integer

Inputs:

- `in`: `F_{p^12}` element
- `P`: `E(F_p)` point
- `Q`: `E(F_{p^12})` point

Outputs:

- `out`: `F_{p^12}` element
- `xP`: `E(F_p)` point equal to `[x]P`

Implementation:

- `out = f_x(P, Q)` where `f_i(P,Q)` is defined recursively using Miller's algorithm:
  - `f_0(P,Q) = in`
  - `f_{i+j}(P,Q) = f_i(P,Q) * f_j(P,Q) * l_{[i]P, [j]P}(Q)`
- `f_x(P, Q)` and `[x]P` are computed using the double-and-add method.

Assumptions:

- `0 <= b < 2^n`
- `0 <= x < 2^250`
- `P` has order `> x` and `P` is not point at infinity
- Requires `(18k)^2 * (2k-1) * 2^{4n} < 2^251` to prevent overflow

### `BLSMillerLoop(n, k, p)`

Proves computation of Miller's algorithm on the curve BLS12-381 for the Tate pairing.

Inputs:

- `P`: `E(F_p)` point
- `Q`: `E(F_{p^12})` point

Output:

- `out`: `F_{p^12}` element equal to `f_r(P, Q)` where `r` is the prime order of torsion in Tate pairing

Implementation:

- Uses the specific form of `r = x^4 - x^2 + 1` where `x` is parameter of BLS12-381, 64-bit integer with low Hamming weight
- Calls `MillerLoop` four times

Assumptions:

- `P` has order `r` (but we do not constrain this) where `r` is prime
- `P` is not point at infinity
