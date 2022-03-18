### final exp
```final_exp.circom``` contains templates for exponentiating elements in $$\mathbb F_{p^{12}}$$. See representations.md (insert link) for a reminder on how points and integers are represented. 

```
FinalExponentiate(n, k, p)
```
Inputs
- **in[6][2][k]**: a $$6\times 2\times k$$ array, representing the $$\mathbb F_{p^{12}}$$ element $$in$$.

Outputs
- **out[6][2][k]**: a $$6\times 2\times k$$ array, representing the $$\mathbb F_{p^{12}}$$ element $$in^{(p^{12}-1)/r}$$, where $$r$$ is the curve parameter