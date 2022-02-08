import random
q  = # a prime
r = # prime divisor of the size of E[F_q]. assume r doesn't divide q-1
bin_rep = [] # binary representation of r, without the leading 1
k = # minimal k with r | q^k - 1
phi = # distortion map sending P to a linearly independent phi(P)

def tate_pairing(P, Q):
    """
    evaluate the tate pairing at points P, Q
    returns an element of F_r
    """
    Q1 = phi(Q)
    return exp(e(P, Q1), (q**k - 1) // r)

def e(P, Q):
    """
    evaluate a bilinear pairing e with 
    e(aP, bQ) = e(P, Q)^{ab}
    returns an element of F_q^k
    """
    S = random()
    f_P1 = miller(P, Q+S)
    f_P2 = miller(P, S)
    return div(f_P1, f_P2)

def miller(P, Q):
    """
    use miller's algorithm to construct a rational function 
    with r zeros at P, r poles at 0, and then evaluate the 
    rational function at Q
    returns an element of F_q^k
    """
    f = 1
    point = P
    for i in range(len(bin_rep)):
        f = mul(mul(f, f), div(interpolate(point, point, Q), vertical(add(point, point), Q)))
        point = add(point, point)
        if bin_rep[i] == 1:
            f = div(mul(f, interpolate(point, P, Q)), vertical(add(point, P), Q))
            point = add(point, P)
    return f


def exp(a, n):
    """
    compute a^n in Fq^k
    """
    return

def mul(a, b):
    """
    multiply two elements of F_q^k
    """
    return 

def add(P, Q):
    """
    add two points on the elliptic curve. 
    done over F_q^k
    """
    return 

def div(a, b):
    """
    divide two elements of F_q^k
    """
    return 

def interpolate(P1, P2, Q):
    """
    construct a linear function through P1, P2; evaluate it at Q. 
    done over F_q^k
    """
    return 

def vertical(P, Q):
    """
    construct a vertical line through P; evaluate it at Q. 
    done over F_q^k
    """
    return 

# TODO: finite field arithmetic + elliptic curves