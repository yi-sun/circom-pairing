import random
import curve

def tate_pairing(P, Q, k, r, phi, S):
    """
    evaluate the optimized Tate pairing at points P, Q
    returns an element of F_r
    
    k is embedding degree
    r is prime divisor of size of E[F_q]. assume r doesn't divide q-1
    phi = distortion map sending P to a linearly independent phi(P)
    S is a random point not in [P, -P, Q, -Q, 0]
    """
    # it seems like for optimized Tate perhaps the random point S can just equal P?
    Q1 = phi(Q)
    return exp(e(P, Q1), (q**k - 1) // r)

# k is never explicitly needed below because we can use internal operators in FQk 

def e(P, Q, r, S):
    """
    evaluate the Tate pairing at points P, Q which are r-torsion
    S is a random point not in [P, -P, Q, -Q, 0]
    
    returns an element of F_q^k
    """ 
    f_P1 = miller(P, add(Q,S), r)
    f_P2 = miller(P, S, r)
    return div(f_P1, f_P2)

def miller(P, Q, r):
    """
    use miller's algorithm to construct a rational function 
    with r zeros at P, r poles at 0, and then evaluate the 
    rational function at Q
    returns an element of F_q^k
    """
    f = 1
    point = P
    while r:
        f = mul(mul(f, f), div(interpolate(point, point, Q), vertical(add(point, point), Q)))
        point = add(point, point)
        if 1 & r:
            f = div(mul(f, interpolate(point, P, Q)), vertical(add(point, P), Q))
            point = add(point, P)
        r //= 2
    return f


def exp(a, n):
    """
    compute a^n in Fq^k
    """
    a_pow = a
    res = None
    while n:
        if n % 2:
            res = res*a_pow if res else a_pow
        n //= 2
        a_pow *= a_pow        
    return res

def mul(a, b):
    """
    multiply two elements of F_q^k
    """
    return a * b 

def add(P, Q):
    """
    add two points on the elliptic curve. 
    done over F_q^k
    """
    return curve.add(P, Q)

def div(a, b):
    """
    divide two elements of F_q^k
    """
    return a / b

def interpolate(P1, P2, Q):
    """
    construct a linear function through P1, P2; evaluate it at Q. 
    done over F_q^k
    """
    lamb = (P2[1]-P1[1]) / (P2[0]-P1[0])
    return lamb * (P1[0]-Q[0]) - (P1[1]-Q[1]) 

def vertical(P, Q):
    """
    construct a vertical line through P; evaluate it at Q. 
    done over F_q^k
    """
    return Q[0] - P[0]

# TODO: finite field arithmetic + elliptic curves