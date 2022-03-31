### fp12
```fp12.circom``` contains templates for operating on elements of $$\mathbb F_{p^{12}}$$. See representations.md (insert link) for a reminder on how points and integers are represented. 

```
Fp12FrobeniusMap(n, k, power)
```
Inputs
- **in[6][2][k]**: a $$6\times 2\times k$$ array representing the $$\mathbb F_{p^{12}}$$ element $$in$$.

Outputs
- **out[6][2][k]**: a $$6\times 2\times k$$ array representing the $$\mathbb F_{p^{12}}$$ element $$out$$, such that $$out = in^{p^{power}}$$ and $$0 < out[i][j] < p$$. Basically an optimized version of exponentiation using properties of the Frobenius map. 

```
Fp12Add(n, k, p)
```
Inputs
- **a[6][2][k]**: a $$6\times 2\times k$$ array representing the $$\mathbb F_{p^{12}}$$ element $$a$$.
- **b[6][2][k]**: a $$6\times 2\times k$$ array representing the $$\mathbb F_{p^{12}}$$ element $$b$$.

Outputs
- **out[6][2][k]**: a $$6\times 2\times k$$ array representing the $$\mathbb F_{p^{12}}$$ element $$out$$, such that $$out = a+b$$ and $$0 \le out[i][j] < p$$.

Assumptions
- $$n\le 252$$
```
Fp12Multiply(n, k, p)
```
Inputs
- **a[6][2][k]**: a $$6\times 2\times k$$ array representing the $$\mathbb F_{p^{12}}$$ element $$a$$.
- **b[6][2][k]**: a $$6\times 2\times k$$ array representing the $$\mathbb F_{p^{12}}$$ element $$b$$.

Outputs
- **out[6][2][k]**: a $$6\times 2\times k$$ array representing the $$\mathbb F_{p^{12}}$$ element $$out$$, such that $$out= a\times b$$ and $$0 \le out[i][j] < p$$.

Assumptions
- $$2n+4 + \log_2 k< 252$$
```
Fp12Invert(n, k, p)
```
Inputs
- **in[6][2][k]**: a $$6\times 2\times k$$ array representing the $$\mathbb F_{p^{12}}$$ element $$in$$.

Outputs
- **out[6][2][k]**: a $$6\times 2\times k$$ array representing the $$\mathbb F_{p^{12}}$$ element $$out$$, such that $$out= in^{-1}$$ and $$0 \le out[i][j] < p$$.
```
Fp12Exp(n, k, e, p)
```
Inputs
- **in[6][2][k]**: a $$6\times 2\times k$$ array representing the $$\mathbb F_{p^{12}}$$ element $$in$$.

Outputs
- **out[6][2][k]**: a $$6\times 2\times k$$ array representing the $$\mathbb F_{p^{12}}$$ element $$out$$, such that $$out = in^e$$ and $$0 < out[i][j] < p$$.

Assumptions
- $$0 < e < 2^{254}$$.