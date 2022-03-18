### curve
```curve.circom``` contains templates for operations with points on an elliptic curve. See representations.md (insert link) for a reminder on how points and integers are represented. 

```
PointOnLine(n, k, p)
```
Inputs
- **in[3][2][k]**: a $$3\times 2\times k$$ array, such that $$in[i]$$ represents the $$\mathbb F_p^2$$ point $$(x_i, y_i)$$. 

Outputs
- No output signal, but the circuit will assert that the points $$(x_0,y_0), (x_1,y_1), (x_2,-y_2)$$ are collinear. 

Assumptions
- $$3n+2\log_2 k + 2< 253$$
```
PointOnCurve(n, k, a, b, p)
```
Inputs
- **in[2][k]**: a $$2\times k$$ array, representing the $$\mathbb F_p^2$$ point $$(x, y)$$. 

Outputs
- No output signal, but the circuit will assert that the point $$(x,y)$$ lies on the curve $$x^3 + ax + b \equiv y^2 \pmod p$$. 

Assumptions
- $$4n+3\log_2 k + 2< 253$$
```
PointOnTangent(n, k, a, p)
```
Inputs
- **in[2][2][k]**: a $$2\times 2\times k$$ array, such that $$in[i]$$ represents the $$\mathbb F_p^2$$ point $$(x_i, y_i)$$. 

Outputs
- No output signal, but the circuit will assert that the point $$(x_1, -y_1)$$ lies on the tangent to the elliptic curve $$x^3+ax+b = y^2$$ at $$(x_0, y_0)$$ (where $$b$$ is chosen so that $$(x_0, y_0)$$ lies on the curve). 

Assumptions
- $$4n+3\log_2 k + 3< 253$$
```
EllipticCurveAddUnequal(n, k, p)
```
Inputs
- **a[2][k]**: a $$2\times k$$ array, representing the $$\mathbb F_p^2$$ point $$(x_0, y_0)$$.
- **b[2][k]**: a $$2\times k$$ array, representing the $$\mathbb F_p^2$$ point $$(x_1, y_1)$$.

Outputs
- **out[2][k]**: a $$2\times k$$ array, representing the $$\mathbb F_p^2$$ point $$(x_2, y_2)$$, such that $$(x_0, y_0) + (x_1, y_1) = (x_2, y_2)$$, where $$+$$ denotes elliptic curve addition. 

Assumptions
- $$4n+3\log_2 k + 4< 253$$
```
EllipticCurveDouble(n, k, p)
```
Inputs
- **in[2][k]**: a $$2\times k$$ array, representing the $$\mathbb F_p^2$$ point $$(x_0, y_0)$$.

Outputs
- **out[2][k]**: a $$2\times k$$ array, representing the $$\mathbb F_p^2$$ point $$(x_2, y_2)$$, such that $$(x_0, y_0) + (x_1, y_1) = (x_2, y_2)$$, where $$+$$ denotes elliptic curve addition. 

Assumptions
- $$4n+3\log_2 k + 4< 253$$
```
LineFunctionUnequal(n, k, q)
```
Inputs
- **P[2][2][k]**: a $$2\times 2\times k$$ array, such that $$P[i]$$ represents the $$\mathbb F_p^2$$ point $$(x_i, y_i)$$.
- **Q[2][6][2][k]**: a $$2\times 6\times 2 \times k$$ array, representing the $$\mathbb F_{p^{12}}^2$$ point $$(X, Y)$$

Outputs
- **out[6][2][k]**: a $$6\times 2\times k$$ array representing the evaluation in $$\mathbb F_{p^{12}}$$ of the line function between $$(x_0,y_0), (x_1, y_1)$$ at the point $$(X,Y)$$
```
LineFunctionEqual(n, k, q)
```
Inputs
- **P[2][k]**: a $$2\times k$$ array, representing the $$\mathbb F_p^2$$ point $$(x_0, y_0)$$.
- **Q[2][6][2][k]**: a $$2\times 6\times 2 \times k$$ array, representing the $$\mathbb F_{p^{12}}^2$$ point $$(X, Y)$$

Outputs
- **out[6][2][k]**: a $$6\times 2\times k$$ array representing the evaluation in $$\mathbb F_{p^{12}}$$ of the line function between $$(x_0,y_0), (x_0, y_0)$$ at the point $$(X,Y)$$

