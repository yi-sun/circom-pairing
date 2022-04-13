pragma circom 2.0.3;

include "bigint.circom";
include "bigint_func.circom";
include "fp.circom";
include "fp2.circom";
include "curve.circom";
include "curve_fp2.circom";
include "bls12_381_func.circom";

/* 
implementation of optimized simplified SWU map to BLS12-381 G2
following Wahby-Boneh: Section 4.2 of https://eprint.iacr.org/2019/403.pdf
Python reference code: https://github.com/algorand/bls_sigs_ref/blob/master/python-impl/opt_swu_g2.py
Additional exposition: https://hackmd.io/@benjaminion/bls12-381#Simplified-SWU-map

E2 is y^2 = x^3 + 4(1+u) over Fp2
E2' is a curve of form y^2 = x^3 + a x + b that is 3-isogenous to E2
Constants are a = 240 u, b = 1012 + 1012 u where u = sqrt(-1)
*/

// Simplified SWU map, optimized and adapted to E2' 
// in = t: 2 x k array, element of Fp2 
// out: 2 x 2 x k array, point (out[0], out[1]) on curve E2' 
// 
// This is osswu2_help(t) in Python reference code
// See Section 4.2 of Wahby-Boneh: https://eprint.iacr.org/2019/403.pdf
// circom implementation is slightly different since sqrt and inversion are cheap
template OptSimpleSWU(n, k){
    signal input in[2][k];
    signal output out[2][2][k]; 
    signal output isInfinity;
    
    var p[50] = get_BLS12_381_prime(n, k);
    
    var a[2] = [0, 240];
    var b[2] = [1012, 1012];
    
    // distinguished non-square in Fp2 for SWU map: xi = -2 - u 
    var xi[2] = [-2, -1];

    var LOGK = log_ceil(k);
    // in = t, compute t^2, t^3
    component t_sq = SignedFp2MultiplyNoCarryCompress(n, k, p, n, 3*n + 2*LOGK + 1 );
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){
        t_sq.a[i][idx] <== in[i][idx];
        t_sq.b[i][idx] <== in[i][idx];
    }
    // compute xi * t^2 
    component xi_t_sq = SignedFp2CarryModP(n, k, 3*n + 2*LOGK + 3, p);
    for(var idx=0; idx<k; idx++){
        xi_t_sq.in[0][idx] <== xi[0] * t_sq.out[0][idx] - xi[1] * t_sq.out[1][idx];
        xi_t_sq.in[1][idx] <== xi[0] * t_sq.out[1][idx] + xi[1] * t_sq.out[0][idx];
    }

    // xi^2 * t^4
    component xi2t4 = SignedFp2MultiplyNoCarryCompress(n, k, p, n, 3*n + 2*LOGK + 1); 
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){
        xi2t4.a[i][idx] <== xi_t_sq.out[i][idx];
        xi2t4.b[i][idx] <== xi_t_sq.out[i][idx];
    }
    // xi^2 * t^4 + xi * t^2
    var num_den_common[2][k];
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        num_den_common[i][idx] = xi2t4.out[i][idx] + xi_t_sq.out[i][idx];
    
    // X0(t) = b (xi^2 * t^4 + xi * t^2 + 1) / (-a * (xi^2 * t^4 + xi * t^2)) 
    component X0_den = SignedFp2CarryModP(n, k, 3*n + 2*LOGK + 2 + 9, p);
    for(var idx=0; idx<k; idx++){
        X0_den.in[0][idx] <== -a[0] * num_den_common[0][idx] + a[1] * num_den_common[1][idx];
        X0_den.in[1][idx] <== -a[0] * num_den_common[1][idx] - a[1] * num_den_common[0][idx];
    }
    
    // Exception if X0_den = 0: 
    var den_zero_total = 2;
    component denIsZero[2];
    for(var i=0; i<2; i++){
        denIsZero[i] = BigIsZero(k);
        for(var idx=0; idx<k; idx++)
            denIsZero[i].in[idx] <== X0_den.out[i][idx];
        den_zero_total -= denIsZero[i].out;
    }
    component exception = IsZero();
    exception.in <== den_zero_total;
    isInfinity <== exception.out; 
    
    num_den_common[0][0]++;
    component X0_num = SignedFp2CarryModP(n, k, 3*n + 2*LOGK + 2 + 11, p);
    for(var idx=0; idx<k; idx++){
        X0_num.in[0][idx] <== b[0] * num_den_common[0][idx] - b[1] * num_den_common[1][idx];
        X0_num.in[1][idx] <== b[0] * num_den_common[1][idx] + b[1] * num_den_common[0][idx];
    }
    // division is same cost/constraints as multiplication, so we will compute X0 and avoid using projective coordinates 
    component X0 = SignedFp2Divide(n, k, n, n, p);
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){
        X0.a[i][idx] <== X0_num.out[i][idx];
        // denominator is X0_den if X0_den != 0, otherwise is 1 
        if(i==0 && idx==0)
            X0.b[i][idx] <== X0_den.out[i][idx] + isInfinity * (1 - X0_den.out[i][idx]);
        else
            X0.b[i][idx] <== X0_den.out[i][idx] - isInfinity * X0_den.out[i][idx];
    }
    
    // g(x) = x^3 + a x + b 
    // Compute g(X0(t)) 
    component gX0 = EllipticCurveFunction(n, k, a, b, p);
    // X1(t) = xi * t^2 * X0(t)
    component X1 = Fp2Multiply(n, k, p);

    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){
        gX0.in[i][idx] <== X0.out[i][idx];
        X1.a[i][idx] <== xi_t_sq.out[i][idx];
        X1.b[i][idx] <== X0.out[i][idx];
    }
    
    component xi3t6 = Fp2MultiplyThree(n, k, p); // shares a hidden component with xi2t4; I'll let compiler optimize that out for readability
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){ 
        xi3t6.a[i][idx] <== xi_t_sq.out[i][idx];
        xi3t6.b[i][idx] <== xi_t_sq.out[i][idx];
        xi3t6.c[i][idx] <== xi_t_sq.out[i][idx];
    }
    // g(X1(t)) = xi^3 * t^6 * g(X0(t)) 
    component gX1 = Fp2Multiply(n, k, p);
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){ 
        gX1.a[i][idx] <== xi3t6.out[i][idx];
        gX1.b[i][idx] <== gX0.out[i][idx];
    }
    /* 
    xi^3 is not a square, so one of gX0, gX1 must be a square 
    isSquare = 1 if gX0 is a square, = 0 if gX1 is a square
    sqrt = sqrt(gX0) if isSquare = 1, sqrt = sqrt(gX1) if isSquare = 0

    Implementation is special to p^2 = 9 mod 16
    References:
        p. 9 of https://eprint.iacr.org/2019/403.pdf
        F.2.1.1 for general version for any field: https://cfrg.github.io/draft-irtf-cfrg-hash-to-curve/draft-irtf-cfrg-hash-to-curve.html#straightline-sswu

    I do not use the trick for combining division and sqrt from Section 5 of 
    Bernstein, Duif, Lange, Schwabe, and Yang, "High-speed high-security signatures",
    since division is cheap in circom
    */
    
    signal isSquare;
    
    // Precompute sqrt_candidate = gX0^{(p^2 + 7) / 16} 
    // p^2 + 7
    var c1[50] = long_add_unequal(n, 2*k, 1, prod(n, k, p, p), [7]);
    // (p^2 + 7) // 16
    var c2[2][50] = long_div2(n, 1, 2*k-1, c1, [16]); 

    assert( c2[1][0] == 0 ); // assert p^2 + 7 is divisible by 16

    var sqrt_candidate[2][50] = find_Fp2_exp(n, k, gX0.out, p, c2[0]);
    // if gX0 is a square, square root must be sqrt_candidate * (8th-root of unity) 
    // -1 is a square in Fp2 (because p^2 - 1 is even) so we only need to check half of the 8th roots of unity
    var roots[4][2][50] = get_roots_of_unity(n, k);
    var sqrt_witness[2][2][50];
    for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        sqrt_witness[i][j][idx] = 0;

    var is_square = 0;
    for(var i=0; i<4; i++){
        var sqrt_tmp[2][50] = find_Fp2_product(n, k, sqrt_candidate, roots[i], p); 
        if(is_equal_Fp2(n, k, find_Fp2_product(n, k, sqrt_tmp, sqrt_tmp, p), gX0.out) == 1){
            is_square = 1;
            sqrt_witness[0] = sqrt_tmp;
        }
    }
    isSquare <-- is_square;
    isSquare * (1-isSquare) === 0; 
    
    var is_square1 = 0;
    var etas[4][2][50] = get_etas(n, k);
    // find square root of gX1 
    // square root of gX1 must be = sqrt_candidate * t^3 * eta 
    // for one of four precomputed values of eta
    // eta determined by eta^2 = xi^3 * (-1)^{-1/4}  
    var t_cu[2][50] = find_Fp2_product(n, k, find_Fp2_product(n, k, in, in, p), in, p);
    sqrt_candidate = find_Fp2_product(n, k, sqrt_candidate, t_cu, p);
    
    for(var i=0; i<4; i++){
        var sqrt_tmp[2][50] = find_Fp2_product(n, k, sqrt_candidate, etas[i], p); 
        if(is_equal_Fp2(n, k, find_Fp2_product(n, k, sqrt_tmp, sqrt_tmp, p), gX1.out) == 1){
            is_square1 = 1;
            sqrt_witness[1] = sqrt_tmp;
        }
    }
    assert(is_square == 1 || is_square1 == 1); // one of gX0 or gX1 must be a square!
        
     
    // X = out[0] = X0 if isSquare == 1, else X = X1    
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        out[0][i][idx] <== isSquare * (X0.out[i][idx] - X1.out[i][idx]) + X1.out[i][idx];  

    // sgn0(t) 
    component sgn_in = Fp2Sgn0(n, k, p);
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        sgn_in.in[i][idx] <== in[i][idx];

    var Y[2][50];
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        Y[i][idx] = is_square * sqrt_witness[0][i][idx] + (1-is_square) * sqrt_witness[1][i][idx];
    // Y = out[1] = +- sqrt_witness; sign determined by sgn0(Y) = sgn0(t) 
    if(get_fp2_sgn0(n, k, Y, p) != sgn_in.out){
        Y[0] = long_sub(n, k, p, Y[0]);
        Y[1] = long_sub(n, k, p, Y[1]);
    } 
    
    component Y_sq = Fp2Multiply(n, k, p);
    // Y^2 == g(X) 
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){
        out[1][i][idx] <-- Y[i][idx];
        Y_sq.a[i][idx] <== out[1][i][idx];
        Y_sq.b[i][idx] <== out[1][i][idx];
    }
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){
        Y_sq.out[i][idx] === isSquare * (gX0.out[i][idx] - gX1.out[i][idx]) + gX1.out[i][idx]; 
    }

    // sgn0(Y) == sgn0(t)
    component sgn_Y = Fp2Sgn0(n, k, p);
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        sgn_Y.in[i][idx] <== out[1][i][idx];

    sgn_Y.out === sgn_in.out;
    
}

