### bigint
```bigint.circom``` contains templates for operating on numbers represented in array form (usually because the number may exceed the size of the circom prime). See representations.md (insert link) for a reminder on how points and integers are represented. 

```
BigAddModP(n, k)
```
Inputs
- **a[k]**: a length-k array representing the BigInt $$a$$.
- **b[k]**: a length-k array representing the BigInt $$b$$.
- **p[k]**: a length-k array representing the BigInt prime $$p$$.

Outputs
- **out[k]**: a length-k array representing the BigInt $$out$$, such that $$out \equiv a+b\pmod p$$ and $$0 \le out < p$$

Assumptions
- $$n\le 252$$
```
BigSubModP(n, k)
```
Inputs
- **a[k]**: a length-k array representing the BigInt $$a$$.
- **b[k]**: a length-k array representing the BigInt $$b$$.
- **p[k]**: a length-k array representing the BigInt prime $$p$$.

Outputs
- **out[k]**: a length-k array representing the BigInt $$out$$, such that $$out \equiv a-b\pmod p$$ and $$0 \le out < p$$

Assumptions
- $$n\le 252$$
```
BigMultModP(n, k)
```
Inputs
- **a[k]**: a length-k array representing the BigInt $$a$$.
- **b[k]**: a length-k array representing the BigInt $$b$$.
- **p[k]**: a length-k array representing the BigInt prime $$p$$.

Outputs
- **out[k]**: a length-k array representing the BigInt $$out$$, such that $$out \equiv a\cdot b\pmod p$$ and $$0 \le out < p$$

Assumptions
- $$n\le 252$$
```
BigModInv(n, k)
```
Inputs
- **in[k]**: a length-k array representing the BigInt $$in$$.
- **p[k]**: a length-k array representing the BigInt prime $$p$$.

Outputs
- **out[k]**: a length-k array representing the BigInt $$out$$, such that $$out \equiv in^{-1}\pmod p$$ and $$0 \le out < p$$

Assumptions
- $$n\le 252$$
- Fails if $$in \equiv 0$$
