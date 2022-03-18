### bigint
```bigint_func.circom``` contains functions for operations on BigInts in represented in array form. See representations.md (insert link) for a reminder on how points and integers are represented. 

```
log_ceil(n)
```
Parameters
- **n**: an integer $$n$$. 

Outputs
- Returns $$\lceil \log_2 n\rceil$$. 

Assumptions
- Assumes $$n < 2^{254}$$.

```
long_add(n, k, a, b)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[k]**: a length-$$k$$ array representing a BigInt $$a$$
- **b[k]**: a length-$$k$$ array representing a BigInt $$b$$

Outputs
- Returns a length-$$(k+1)$$ array **out** representing a BigInt $$out$$ with $$out = a + b$$. 

Assumptions
- Assumes $$k+1\le 50$$. 
```
long_add_unequal(n, k1, k2, a, b)
```
Parameters
- **n**: the register size. 
- **k1**: the number of registers in $$a$$.
- **k2**: the number of registers in $$b$$. 
- **a[k1]**: a length-$$k1$$ array representing a BigInt $$a$$
- **b[k2]**: a length-$$k2$$ array representing a BigInt $$b$$

Outputs
- Returns a length-$$(k1+1)$$ array **out** representing a BigInt $$out$$ with $$out = a + b$$. 

Assumptions
- Assumes $$k1 \ge k2$$ and $$k_1 + 1\le 50$$. 
```
long_sub(n, k, a, b)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[k]**: a length-$$k$$ array representing a BigInt $$a$$
- **b[k]**: a length-$$k$$ array representing a BigInt $$b$$

Outputs
- Returns a length-$$k$$ array **out** representing a BigInt $$out$$ with $$out = a - b$$. 

Assumptions
- Assumes $$a\ge b$$ and $$k\le 50$$. 

```
long_scalar_mult(n, k, a, b)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a**: an $$n$$-bit integer.
- **b[k]**: a length-$$k$$ array representing a BigInt $$b$$.

Outputs
- Returns a length-$$(k+1)$$ array **out** representing a BigInt $$out$$ with $$out = a\times b$$.

Assumptions
- Assumes $$k+1\le 50$$. 

```
long_div2(n, k, m, a, b)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers in $$b$$.
- **m**: an integer such that the number of registers in $$a$$ is $$k+m$$. 
- **a[k+m]**: a length-$$(k+m)$$ array representing a BigInt $$a$$
- **b[k]**: a length-$$k$$ array representing a BigInt $$b$$

Outputs
- Returns a $$2\times \max(m+1, k)$$ array **out** representing two BigInts $$out[0], out[1]$$ with $$a = out[0] \cdot b + out[1]$$.

Assumptions
- Assumes $$\max (m+1, k) \le 50$$. 

```
long_to_short(n, k, a)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[k]**: a length-$$k$$ array representing a BigInt $$a$$, except the entries of $$a$$ may exceed $$2^n$$.

Outputs
- Returns a length-$$50$$ array **out** representing $$a$$ in standard BigInt form (with entries less than $$2^n$$).

Assumptions
- Assumes $$a[i] < 2^{252}$$ for each $$i$$ and that the proper representation of $$a$$ uses at most $$50$$ registers. 

```
prod(n, k, a, b)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[k]**: a length-$$k$$ array representing a BigInt $$a$$
- **b[k]**: a length-$$k$$ array representing a BigInt $$b$$

Outputs
- Returns a length-$$(2k-1)$$ array **out** representing a BigInt $$out$$ with $$out = a\cdot b$$. 

Assumptions
- Assumes $$k\le 25$$. 
```
long_add_mod(n, k, a, b, p)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[k]**: a length-$$k$$ array representing a BigInt $$a$$
- **b[k]**: a length-$$k$$ array representing a BigInt $$b$$
- **p[k]**: a length-$$k$$ array representing a BigInt prime $$p$$

Outputs
- Returns a length-$$k$$ array **out** representing a BigInt $$out$$ with $$out \equiv a + b \pmod p$$. 

Assumptions
- Assumes $$k\le 50$$. 
```
long_sub_mod(n, k, a, b, p)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[k]**: a length-$$k$$ array representing a BigInt $$a$$
- **b[k]**: a length-$$k$$ array representing a BigInt $$b$$
- **p[k]**: a length-$$k$$ array representing a BigInt prime $$p$$

Outputs
- Returns a length-$$k$$ array **out** representing a BigInt $$out$$ with $$out \equiv a - b \pmod p$$. 

Assumptions
- Assumes $$k\le 50$$. 
```
prod_mod(n, k, a, b, p)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[k]**: a length-$$k$$ array representing a BigInt $$a$$
- **b[k]**: a length-$$k$$ array representing a BigInt $$b$$
- **p[k]**: a length-$$k$$ array representing a BigInt prime $$p$$

Outputs
- Returns a length-$$k$$ array **out** representing a BigInt $$out$$ with $$out \equiv a \times b \pmod p$$. 

Assumptions
- Assumes $$k\le 50$$. 
```
mod_exp(n, k, a, p, e)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[k]**: a length-$$k$$ array representing a BigInt $$a$$
- **p[k]**: a length-$$k$$ array representing a BigInt prime $$p$$
- **e**: a length-$$k$$ array representing a BigInt exponent $$e$$

Outputs
- Returns a length-$$k$$ array **out** representing a BigInt $$out$$ with $$out \equiv a^e\pmod p$$. 

Assumptions
- Assumes $$k\le 50$$ and $$k\cdot n\le 500$$. 
```
mod_inv(n, k, a, p)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[k]**: a length-$$k$$ array representing a BigInt $$a$$
- **p[k]**: a length-$$k$$ array representing a BigInt prime $$p$$

Outputs
- Returns a length-$$k$$ array **out** representing a BigInt $$out$$ with $$out \equiv a^{-1}\pmod p$$. Returns $$0$$ if $$a=0$$. 

Assumptions
- Assumes $$k\le 50$$ and $$k\cdot n\le 500$$. 
```
get_Fp_carry_witness(n, k, m, a, p)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **m**: an integer representing the number of overflow bits in $$a$$
- **a[2][k]**: a $$2\times k$$ array representing two BigInts $$a_0, a_1$$, except the registers can go up to $$2^{k+m}$$.
- **p[k]**: a length-$$k$$ array representing a BigInt prime $$p$$

Outputs
- Returns a $$2\times \max (k, m)$$ array $$out$$ representing two BigInts $$out_0, out_1$$ with $$a[0] - a[1] = p\cdot out_0 + out_1$$. $$out_0$$ may have negative entries. 

Assumptions
- Assumes $$k, m\le 50$$. 