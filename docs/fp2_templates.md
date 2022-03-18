### fp2
```fp2.circom``` contains templates for operating on elements of $$\mathbb F_{p^2}$$. See representations.md (insert link) for a reminder on how points and integers are represented. 

```
Fp2Add(n, k, p)
```
Inputs
- **a[2][k]**: a $$2\times k$$ array representing the $$\mathbb F_{p^2}$$ element $$a$$.
- **b[2][k]**: a $$2\times k$$ array representing the $$\mathbb F_{p^2}$$ element $$b$$.

Outputs
- **out[2][k]**: a $$2\times k$$ array representing the $$\mathbb F_{p^2}$$ element $$out$$, such that $$out = a+b$$ and $$0 \le out[0], out[1] < p$$

Assumptions
- $$n\le 252$$
```
Fp2Subtract(n, k, p)
```
Inputs
- **a[2][k]**: a $$2\times k$$ array representing the $$\mathbb F_{p^2}$$ element $$a$$.
- **b[2][k]**: a $$2\times  k$$ array representing the $$\mathbb F_{p^2}$$ element $$b$$.

Outputs
- **out[2][k]**: a $$2\times k$$ array representing the $$\mathbb F_{p^2}$$ element $$out$$, such that $$out = a-b$$ and $$0 \le out[0], out[1] < p$$

Assumptions
- $$n\le 252$$
```
Fp2Multiply(n, k, p)
```
Inputs
- **a[2][k]**: a $$2\times k$$ array representing the $$\mathbb F_{p^2}$$ element $$a$$.
- **b[2][k]**: a $$2\times k$$ array representing the $$\mathbb F_{p^2}$$ element $$b$$.

Outputs
- **out[2][k]**: a $$2\times k$$ array representing the $$\mathbb F_{p^2}$$ element $$out$$, such that $$out= a\times b$$ and $$0 \le out[0], out[1] < p$$

Assumptions
- $$3n+2 + 2\log_2 k< 252$$
```
Fp2Invert(n, k, p)
```
Inputs
- **in[2][k]**: a $$2\times k$$ array representing the $$\mathbb F_{p^2}$$ element $$in$$.

Outputs
- **out[2][k]**: a $$2\times k$$ array representing the $$\mathbb F_{p^2}$$ element $$out$$, such that $$out= in^{-1}$$ and $$0 \le out[0], out[1] < p$$

Assumptions
- $$in\neq 0$$
