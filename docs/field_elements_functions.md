### field_elements
```field_elements_func.circom``` contains functions for operations involving field elements. See representations.md (insert link) for a reminder on how points and integers are represented. 

```
find_Fp2_product(n, k, a, b, p)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[2][k]**: a $$2\times k$$ array representing an $$\mathbb F_{p^2}$$ element $$a$$
- **b[2][k]**: a $$2\times k$$ array representing an $$\mathbb F_{p^2}$$ element $$b$$
- **p[k]**: a length-$$k$$ array representing a BigInt prime $$p$$

Outputs
- Returns a $$2\times k$$ array **out** representing an $$\mathbb F_{p^2}$$ element $$out$$ with $$out = a \cdot b$$. 

```
find_Fp2_sum(n, k, a, b, p)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[2][k]**: a $$2\times k$$ array representing an $$\mathbb F_{p^2}$$ element $$a$$
- **b[2][k]**: a $$2\times k$$ array representing an $$\mathbb F_{p^2}$$ element $$b$$
- **p[k]**: a length-$$k$$ array representing a BigInt prime $$p$$

Outputs
- Returns a $$2\times k$$ array **out** representing an $$\mathbb F_{p^2}$$ element $$out$$ with $$out = a + b$$. 
```
find_Fp2_diff(n, k, a, b, p)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[2][k]**: a $$2\times k$$ array representing an $$\mathbb F_{p^2}$$ element $$a$$
- **b[2][k]**: a $$2\times k$$ array representing an $$\mathbb F_{p^2}$$ element $$b$$
- **p[k]**: a length-$$k$$ array representing a BigInt prime $$p$$

Outputs
- Returns a $$2\times k$$ array **out** representing an $$\mathbb F_{p^2}$$ element $$out$$ with $$out = a - b$$. 

```
find_Fp2_inverse(n, k, a, p)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[2][k]**: a $$2\times k$$ array representing an $$\mathbb F_{p^2}$$ element $$a$$
- **p[k]**: a length-$$k$$ array representing a BigInt prime $$p$$

Outputs
- Returns a $$2\times k$$ array **out** representing an $$\mathbb F_{p^2}$$ element $$out$$ with $$out = a^{-1}$$.

```
find_Fp12_sum(n, k, a, b, p)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[6][2][k]**: a $$6\times 2\times k$$ array representing an $$\mathbb F_{p^{12}}$$ element $$a$$
- **b[6][2][k]**: a $$6\times 2\times k$$ array representing an $$\mathbb F_{p^{12}}$$ element $$b$$
- **p[k]**: a length-$$k$$ array representing a BigInt prime $$p$$

Outputs
- Returns a $$6\times 2\times k$$ array **out** representing an $$\mathbb F_{p^{12}}$$ element $$out$$ with $$out = a + b$$. 

```
find_Fp2_diff(n, k, a, b, p)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[6][2][k]**: a $$6\times 2\times k$$ array representing an $$\mathbb F_{p^{12}}$$ element $$a$$
- **b[6][2][k]**: a $$6\times 2\times k$$ array representing an $$\mathbb F_{p^{12}}$$ element $$b$$
- **p[k]**: a length-$$k$$ array representing a BigInt prime $$p$$

Outputs
- Returns a $$6\times 2\times k$$ array **out** representing an $$\mathbb F_{p^{12}}$$ element $$out$$ with $$out = a - b$$. 

```
find_Fp12_product(n, k, a, b, p)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **a[6][2][k]**: a $$6\times 2\times k$$ array representing an $$\mathbb F_{p^{12}}$$ element $$a$$
- **b[6][2][k]**: a $$6\times 2\times k$$ array representing an $$\mathbb F_{p^{12}}$$ element $$b$$
- **p[k]**: a length-$$k$$ array representing a BigInt prime $$p$$

Outputs
- Returns a $$6\times 2\times k$$ array **out** representing an $$\mathbb F_{p^{12}}$$ element $$out$$ with $$out = a \cdot b$$. 

```
find_Fp12_inverse(n, k, p, a)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **p[k]**: a length-$$k$$ array representing a BigInt prime $$p$$
- **a[6][2][k]**: a $$6\times 2\times k$$ array representing an $$\mathbb F_{p^{12}}$$ element $$a$$

Outputs
- Returns a $$6\times 2\times k$$ array **out** representing an $$\mathbb F_{p^{12}}$$ element $$out$$ with $$out = a^{-1}$$.

```
get_Fp12_carry_witness(n, k, kX, p, a, b)
```
Parameters
- **n**: the register size. 
- **k**: the number of registers. 
- **kX**: a bound on the number of registers in the output. 
- **p[k]**: a length-$$k$$ array representing a BigInt prime
- **a[k]**: a length-$$k$$ array representing an $$\mathbb F_{p}$$ element $$a$$
- **b[k]**: a length-$$k$$ array representing an $$\mathbb F_{p}$$ element $$b$$

Outputs
- Returns a $$2\times kX$$ array **out[2][kX]** such that $$a - b = out[0] * p + out[1]$$. 

Assumptions
- Assumes $$\lfloor (a-b)/p\rfloor$$ has at most $$kX$$ registers. 