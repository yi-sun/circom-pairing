from curve_field_elements import field_modulus, FQ, FQ2, FQ12, inv

def numberToBase(num, b):
    num = abs(num)
    # assume num >= 0
    if num==0:
        return [0]
    registers = []
    while num:
        registers.append(int(num % b))
        num //= b
    return registers

def numberToArray(num, n, k):
    num = abs(num)
    # assume num >= 0
    registers = []
    for i in range(k):
        registers.append(int(num % (2**n)))
        num //= 2**n
    return registers

def hamming_weight(x):
    num = abs(x)
    ones = 0;
    while num:
        if num & 1:
            ones = ones + 1
        num //= 2
    return ones

def printEllipticPoint(P, n, k):
    print(numberToArray(num=P[0].n, n=n, k=k), numberToArray(num=P[1].n, n=n, k=k))

def printFQ(x, n, k):
    print("[", end="")
    A = numberToArray(x.n, n, k)
    for idx in range(len(A)):
        print(f'"{A[idx]}"', end="")
        if idx != len(A)-1:
            print(",")
    print("]", end="")
    
def printFQ2(X, n, k):
    in22 = X.coeffs
    print("[")
    for i in range(len(in22)):
        printFQ(in22[i], n, k)
        if i != len(in22)-1:
            print(",")
    print("]")

def Fp12convert(X, n, k, xi=1):
    basis1 = X.coeffs
    ret = []
    for i in range(6):
        fq2elt = FQ2([basis1[i].n, 0]) + FQ2([basis1[i+6].n, 0]) * FQ2([xi,1])
        ret.append([ numberToArray(fq2elt.coeffs[0].n, n, k) , numberToArray(fq2elt.coeffs[1].n, n, k) ])
    return ret

def convert_out_to_Fp12(out, n, k):
    twelve = []
    for i in range(6):
        fp2 = []
        for j in range(2):
            num = 0;
            for idx in range(k):
                num = num + int(out[2*k * i + j*k + idx]) * (2 ** (n*idx))
            fp2.append(num)
        twelve.append( fp2 )
    coeff = [0] * 12
    for i in range(6):
        coeff[i] = twelve[i][0] - twelve[i][1]
        coeff[i+6] = twelve[i][1]
    return FQ12(coeff)

def printFQ12(X, n, k, xi=1):
    in62 = Fp12convert(X, n, k, xi)
    print("[")
    for i in range(len(in62)):
        print("[", end="")
        C = in62[i]
        for j in range(len(C)):
            print("[", end="")
            A = C[j]
            for idx in range(len(A)):
                print(f'"{A[idx]}"', end="")
                if idx != len(A)-1:
                    print(",")
            print("]", end="")
            if j != len(C)-1:
                print(",")
        print("]", end="")
        if i != len(in62)-1:
            print(",")
    print("]")

def print_fq12_frobenius_coeff(q, n, k, xi=1):
    gamma = [[0]*6]*12
    for j in range(12):
        gamma[j] = [ FQ2([xi,1]) ** ( (i*(q**j-1)//6) % (q**2-1) ) for i in range(6)]
    for j in range(12):
        for i in range(6):
            A, B = gamma[j][i].coeffs
            a = A.n
            b = B.n
            for r in range(k):
                print(f"coeff[{j}][{i}][0][{r}] = {a%(2**n)};")
                a //= 2**n
            print("")
            for r in range(k):
                print(f"coeff[{j}][{i}][1][{r}] = {b%(2**n)};")
                b //= 2**n
            print("")


