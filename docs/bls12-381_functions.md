### bls12-381
```bls12-381_func.circom``` contains functions for retrieving information about the BLS12-381 curve. See representations.md (insert link) for a reminder on how points and integers are represented. 

```
get_BLS12_381_parameter()
```
Parameters
- No parameters.

Outputs
- Returns the parameter $$x$$ of the BLS12-381 curve. 

```
get_BLS12_381_prime(n, k)
```
Parameters
- **n**: the register size.
- **k**: the number of registers. 

Outputs
- Returns a length-$$50$$ array of size-$$n$$ registers, such that the first $$k$$ registers together contain the BLS12-381 prime $$p$$ in BigInt representation. 

Assumptions
- $$(n,k) \in \{ (96, 4), (77, 5)\}$$. 
```
get_Fp12_frobenius(n, k)
```
Parameters
- **n**: the register size.
- **k**: the number of registers. 

Outputs
- Returns a $$12\times 6\times 2\times 20$$ array $$coeff$$ of size-$$n$$ registers, such that each $$coeff[j][i]$$ represents an element of $$\mathbb F_{p^2}$$ with $$(w^i)^{(p^j)} = coeff[j][i] \cdot w^i$$.

Assumptions
- $$(n,k) \in \{ (96, 4), (77, 5)\}$$. 
