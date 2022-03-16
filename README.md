# circom-pairing

TODO: add table of contents with links

## Setup
First, install [yarn](https://classic.yarnpkg.com/en/) and [circom](https://docs.circom.io/getting-started/installation/). 
Then run ```yarn install``` in the root directory to install the dependencies in ```yarn.lock```. 

## Testing
See the ```/test``` directory for examples of tests. The circuits to be tested should be written in the ```/test/circuits``` folder, while the test execution code should be written in regular JavaScript files under ```/test```. A short description of each test can be passed in as the first parameter of the ```describe()``` function, and ```yarn --grep name``` will run all tests whose description contains ```name``` as a substring. 

## Useful templates
Remember that templates should be used to enforce and verify constraints on signals, and contribute to constraint costs. 
### bigint
```bigint.circom``` contains templates for operating on numbers represented in array form (usually because the number may exceed the size of the circom prime)

```
BigAddModP(n, k)
```
Inputs
- **a[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$a=2^0 \cdot a[0] + 2^n \cdot a[1] + \dots + 2^{n(k-1)} \cdot a[k-1]$$.
- **b[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$b=2^0 \cdot b[0] + 2^n \cdot b[1] + \dots + 2^{n(k-1)} \cdot b[k-1]$$.
- **p[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$p = 2^0\cdot p[0] + 2^n\cdot p[1] + \dots + 2^{n(k-1)} \cdot p[k-1]$$.

Outputs
- **out[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$out = 2^0\cdot p[0] + 2^n\cdot p[1] + \dots + 2^{n(k-1)} \cdot p[k-1]$$, such that $$out \equiv a+b\pmod p$$ and $$0 \le out < p$$

Assumptions
- $$n\le 252$$
```
BigSubModP(n, k)
```
Inputs
- **a[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$a=2^0 \cdot a[0] + 2^n \cdot a[1] + \dots + 2^{n(k-1)} \cdot a[k-1]$$.
- **b[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$b=2^0 \cdot b[0] + 2^n \cdot b[1] + \dots + 2^{n(k-1)} \cdot b[k-1]$$.
- **p[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$p = 2^0\cdot p[0] + 2^n\cdot p[1] + \dots + 2^{n(k-1)} \cdot p[k-1]$$.

Outputs
- **out[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$out = 2^0\cdot p[0] + 2^n\cdot p[1] + \dots + 2^{n(k-1)} \cdot p[k-1]$$, such that $$out \equiv a-b\pmod p$$ and $$0\le out < p$$

Assumptions
- $$n\le 252$$
```
BigMultModP(n, k)
```
Inputs
- **a[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$a=2^0 \cdot a[0] + 2^n \cdot a[1] + \dots + 2^{n(k-1)} \cdot a[k-1]$$.
- **b[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$b=2^0 \cdot b[0] + 2^n \cdot b[1] + \dots + 2^{n(k-1)} \cdot b[k-1]$$.
- **p[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$p = 2^0\cdot p[0] + 2^n\cdot p[1] + \dots + 2^{n(k-1)} \cdot p[k-1]$$.

Outputs
- **out[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$out = 2^0\cdot p[0] + 2^n\cdot p[1] + \dots + 2^{n(k-1)} \cdot p[k-1]$$, such that $$out \equiv a\cdot b\pmod p$$ and $$0\le out < p$$

Assumptions
- $$n\le 252$$
```
BigModInv(n, k)
```
Inputs
- **in[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$in=2^0 \cdot in[0] + 2^n \cdot in[1] + \dots + 2^{n(k-1)} \cdot in[k-1]$$.
- **p[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$p = 2^0\cdot p[0] + 2^n\cdot p[1] + \dots + 2^{n(k-1)} \cdot p[k-1]$$.

Outputs
- **out[k]**: a length-k array where each element is between $$0$$ and $$2^n-1$$, representing the number $$out = 2^0\cdot p[0] + 2^n\cdot p[1] + \dots + 2^{n(k-1)} \cdot p[k-1]$$, such that $$out \equiv in^{-1}\pmod p$$ and $$0\le out < p$$

Assumptions
- $$n\le 252$$
- Fails if $$in \equiv 0$$

## Useful functions
Remember that functions should be used to perform computations to generate signals, at no constraint cost. 