pragma circom 2.0.3;

include "../curve.circom";
include "curve_fp2.circom";
include "bn254_func.circom";

// Untwist-Frobenius-twist endomorphism, introduced by Galbraith-Scott: https://eprint.iacr.org/2008/117.pdf

// Input: in = (x, y) in F_{p^2}^2 
// Output: out = psi(x,y) = ( (9 + u)^{(p-1)/3} x^p, (XI0 + u)^{(p-1)/2} y^p ) 
template EndomorphismPsi(n, k, p){
    signal input in[2][2][k];
    signal output out[2][2][k];
    
    // coeff[1][j] = (9+u)^{(p-1)/6 * j}
    var coeff[12][6][2][20] = get_Fp12_frobenius(n, k);
    component frob[2];
    component qx[2];
    for(var i=0; i<2; i++){
        frob[i] = Fp2FrobeniusMap(n, k, 1, p);
        for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
            frob[i].in[j][idx] <== in[i][j][idx];

        qx[i] = Fp2Multiply(n, k, p);
        for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            qx[i].a[j][idx] <== coeff[1][2+i][j][idx]; 
            qx[i].b[j][idx] <== frob[i].out[j][idx];
        } 
    }

    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        out[0][j][idx] <== qx[0].out[j][idx];
        out[1][j][idx] <== qx[1].out[j][idx];
    }
}

/* 
Subgroup check for G1: 
- G1 cofactor is 1 so just check point on curve 
*/

// `in` = P is 2 x k, pair of Fp elements
// check P is on curve E(Fp)
template SubgroupCheckG1(n, k){
    signal input in[2][k];

    var b = 3;
    var p[50] = get_bn254_prime(n, k);

    component is_on_curve = PointOnCurve(n, k, 0, b, p);
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        is_on_curve.in[i][idx] <== in[i][idx];
}

/*
Subgroup check for G2: 
use the latest method by El Housni, Guillevic, Piellard: https://eprint.iacr.org/2022/352.pdf
By Proposition 3, enough to check psi(P) = [lambda]P for lambda = t - 1 = 6 * x^2
*/

// `in` = P is 2 x 2 x k, pair of Fp2 elements 
// check P is on curve twist E2(Fp2)
// check psi(P) = [lambda]P where lambda = 6x^2, x is parameter for BN254
template SubgroupCheckG2(n, k){
    signal input in[2][2][k];
    
    var p[50] = get_bn254_prime(n, k);
    var x = get_bn254_parameter();
    var b[2][50] = get_bn254_b(n, k);

    component is_on_curve = PointOnCurveFp2(n, k, [0,0], b[0], b[1], p);

    for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        is_on_curve.in[i][j][idx] <== in[i][j][idx]; 

    component psiP = EndomorphismPsi(n, k, p); 
    var lambda = 6 * x * x;
    component lambdaP = EllipticCurveScalarMultiplyUnequalFp2(n, k, b[0], b[1], lambda, p); 

    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        psiP.in[0][j][idx] <== in[0][j][idx];
        psiP.in[1][j][idx] <== in[1][j][idx];
        lambdaP.in[0][j][idx] <== in[0][j][idx];
        lambdaP.in[1][j][idx] <== in[1][j][idx]; 
    }
    
    // psi(P) == [x]P
    component is_eq[2];
    for(var i=0; i<2; i++){
        is_eq[i] = Fp2IsEqual(n, k, p);
        for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            is_eq[i].a[j][idx] <== psiP.out[i][j][idx];
            is_eq[i].b[j][idx] <== lambdaP.out[i][j][idx];
        }
    }
    is_eq[0].out === 1;
    is_eq[1].out === 1;
}

