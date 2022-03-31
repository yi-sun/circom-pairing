### pairing
```pairing.circom``` contains templates for elliptic curve pairings. See representations.md (insert link) for a reminder on how points and integers are represented. 

```
BLSTatePairing(n, k, q)
```
Inputs
- **P[2][k]**: a $$2\times k$$ array representing a $$\mathbb F_q^2$$ point on the BLS12-381 curve.
- **Q[2][6][2][k]**: a $$2\times 6\times 2 \times k$$ array representing a point in $$\mathbb F_{q^{12}}^2$$

Outputs
- **out[6][2][k]**: a $$6\times 2\times k$$ array representing an element of $$\mathbb F_{q^{12}}$$, which is the output of the Tate pairing between points $$P,Q$$ with respect to the BLS12-381 curve. 