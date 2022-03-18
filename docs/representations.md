### BigInt Representation
Almost every template we use contains two parameters $$n$$ and $$k$$. $$n$$ denotes the register size and $$k$$ denotes the number of registers; for a large integer $$A$$ we can represent it as an array $$a$$ of length $$k$$ such that $$A = 2^0 \cdot a[0] + 2^n \cdot a[1] + \dots + 2^{n(k-1)} \cdot a[k-1]$$. 

### Fp Point Representation
Given a point $$A = (A_x, A_y)$$ in $$\mathbb F_p^2$$, we represent it as a $$2\times k$$ array $$a$$ such that $$a[0], a[1]$$ correspond to $$A_x,A_y$$ under the BigInt representation. 

### Fp2 Point Representation
Given a point $$A = (A_x, A_y)$$ in $$\mathbb F_{p^2}^2$$, we represent it as a $$2\times 2 \times k$$ array $$a$$ such that $$A_x$$ corresponds to $$a[0][0] + u a[0][1]$$ and $$A_y$$ corresponds to $$a[1][0] + u a[1][1]$$, where $$a[i][j]$$ are BigInts. Here, $$u$$ is an element of $$\mathbb F_{p^2}$$ with $$u^2=-1$$.

### Fp12 Point Representation
Given a point $$A = (A_x, A_y)$$ in $$\mathbb F_{p^{12}}^2$$, we represent it as a $$2\times 6\times 2 \times k$$ array $$a$$ such that $$A_x$$ corresponds to $$\sum_{i,j} a[0][i][j] w^i u^j$$ and $$A_y$$ corresponds to $$\sum_{i,j} a[1][i][j] \cdot w^i \cdot u^j$$, where $$a[i][j][k]$$ are BigInts. Here, $$u$$ is an element of $$\mathbb F_{p^2}$$ with $$u^2=-1$$, while $$w$$ is an element of $$\mathbb F_{p^{12}}$$ with $$w^6=u+1$$. 