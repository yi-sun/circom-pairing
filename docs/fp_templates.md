### fp
```fp.circom``` contains templates for operating on elements of $$\mathbb F_p$$. See representations.md (insert link) for a reminder on how points and integers are represented. 

```
FpAdd(n, k, p)
```
Inputs
- **a[k]**: a length-k array representing the $$\mathbb F_p$$ element $$a$$.
- **b[k]**: a length-k array representing the $$\mathbb F_p$$ element $$b$$.

Outputs
- **out[k]**: a length-k array representing the $$\mathbb F_p$$ element $$out$$, such that $$out \equiv a+b\pmod p$$ and $$0 \le out < p$$

Assumptions
- $$n\le 252$$
```
FpSubtract(n, k)
```
Inputs
- **a[k]**: a length-k array representing the $$\mathbb F_p$$ element $$a$$.
- **b[k]**: a length-k array representing the $$\mathbb F_p$$ element $$b$$.

Outputs
- **out[k]**: a length-k array representing the $$\mathbb F_p$$ element $$out$$, such that $$out \equiv a-b\pmod p$$ and $$0 \le out < p$$

Assumptions
- $$n\le 252$$
```
FpMultiply(n, k)
```
Inputs
- **a[k]**: a length-k array representing the $$\mathbb F_p$$ element $$a$$.
- **b[k]**: a length-k array representing the $$\mathbb F_p$$ element $$b$$.

Outputs
- **out[k]**: a length-k array representing the $$\mathbb F_p$$ element $$out$$, such that $$out \equiv a\times b\pmod p$$ and $$0 \le out < p$$

Assumptions
- $$n\le 252$$
```
CheckCarryModToZero(n, k, overflow, p)
```
Inputs
- **in[2][k]**: a $$2\times k$$ array such that $$in[i]$$ represents the $$\mathbb F_p$$ element $$in_i$$

Outputs
- **X[m]**: a length-m array representing $$\frac{in_0 - in_1}{p}$$, where $$m = \lceil \text{overflow}/n\rceil$$. 

Assumptions
- $$in_0 \equiv in_1 \pmod p$$. 
- $$overflow < 252$$. 
