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
template OptSimpleSWU2(n, k){
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

/*
3-Isogeny from E2' to E2
References:
    Appendix E.3 of https://cfrg.github.io/draft-irtf-cfrg-hash-to-curve/draft-irtf-cfrg-hash-to-curve.html#name-3-isogeny-map-for-bls12-381
    Section 4.3 of Wahby-Boneh: https://eprint.iacr.org/2019/403.pdf
    iso3(P) in Python reference code: https://github.com/algorand/bls_sigs_ref/blob/master/python-impl/opt_swu_g2.py
*/ 

// Input:
//  in = (x', y') point on E2' 
//  inIsInfinity = 1 if input is point at infinity on E2' (in which case x', y' are arbitrary)
// Output:
//  out = (x, y) is point on E2 after applying 3-isogeny to in 
//  isInfinity = 1 if one of exceptional cases occurs and output should be point at infinity
// Exceptions:
//  inIsInfinity = 1
//  input is a pole of the isogeny, i.e., x_den or y_den = 0 
template Iso3Map(n, k){
    signal input in[2][2][k];
    signal input inIsInfinity;
    signal output out[2][2][k];
    signal output isInfinity;
    
    var p[50] = get_BLS12_381_prime(n, k);
    
    // load coefficients of the isogeny (precomputed)
    var coeffs[4][4][2][50] = get_iso3_coeffs(n, k);

    // x = x_num / x_den
    // y = y' * y_num / y_den
    // x_num = sum_{i=0}^3 coeffs[0][i] * x'^i
    // x_den = x'^2 + coeffs[1][1] * x' + coeffs[1][0] 
    // y_num = sum_{i=0}^3 coeffs[2][i] * x'^i
    // y_den = x'^3 + sum_{i=0}^2 coeffs[3][i] * x'^i
  
    var LOGK = log_ceil(k); 
    component xp2_nocarry = SignedFp2MultiplyNoCarry(n, k, 2*n + LOGK + 1); 
    component xp2 = SignedFp2CompressCarry(n, k, k-1, 2*n+LOGK+1, p);
    component xp3_nocarry = SignedFp2MultiplyNoCarryUnequal(n, 2*k-1, k, 3*n + 2*LOGK + 2); 
    component xp3 = SignedFp2CompressCarry(n, k, 2*k-2, 3*n+2*LOGK+2, p);

    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){
        xp2_nocarry.a[i][idx] <== in[0][i][idx];
        xp2_nocarry.b[i][idx] <== in[0][i][idx];
        xp3_nocarry.b[i][idx] <== in[0][i][idx];
    }
    for(var i=0; i<2; i++)for(var idx=0; idx<2*k-1; idx++){
        xp3_nocarry.a[i][idx] <== xp2_nocarry.out[i][idx];
        xp2.in[i][idx] <== xp2_nocarry.out[i][idx];
    }
    for(var i=0; i<2; i++)for(var idx=0; idx<3*k-2; idx++)
        xp3.in[i][idx] <== xp3_nocarry.out[i][idx]; 

    signal xp_pow[3][2][k]; 
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){
        xp_pow[0][i][idx] <== in[0][i][idx];
        xp_pow[1][i][idx] <== xp2.out[i][idx];
        xp_pow[2][i][idx] <== xp3.out[i][idx];
    }
     
    component coeffs_xp[4][3]; 
    var deg[4] = [3, 1, 3, 2];
    for(var i=0; i<4; i++)for(var j=0; j<deg[i]; j++){
        coeffs_xp[i][j] = SignedFp2MultiplyNoCarry(n, k, 2*n + LOGK + 1);
        for(var l=0; l<2; l++)for(var idx=0; idx<k; idx++){
            coeffs_xp[i][j].a[l][idx] <== coeffs[i][j+1][l][idx];
            coeffs_xp[i][j].b[l][idx] <== xp_pow[j][l][idx];
        }
    }
    var x_frac[4][2][50];  
    for(var i=0; i<4; i++){
        for(var l=0; l<2; l++)for(var idx=0; idx<2*k-1; idx++){
            if(idx<k)
                x_frac[i][l][idx] = coeffs[i][0][l][idx];
            else
                x_frac[i][l][idx] = 0;
        }
        for(var j=0; j<deg[i]; j++)for(var l=0; l<2; l++)for(var idx=0; idx<2*k-1; idx++)
            x_frac[i][l][idx] += coeffs_xp[i][j].out[l][idx];
    } 
    for(var l=0; l<2; l++)for(var idx=0; idx<k; idx++){
        x_frac[1][l][idx] += xp2.out[l][idx];
        x_frac[3][l][idx] += xp3.out[l][idx];
    }
    
    // carry the denominators since we need to check whether they are 0
    component den[2];
    component den_is_zero[2];
    for(var i=0; i<2; i++){
        den[i] = SignedFp2CompressCarry(n, k, k-1, 2*n + LOGK + 3, p); 
        for(var l=0; l<2; l++)for(var idx=0; idx<2*k-1; idx++)
            den[i].in[l][idx] <== x_frac[2*i+1][l][idx];

        den_is_zero[i] = Fp2IsZero(k);
        for(var l=0; l<2; l++)for(var idx=0; idx<k; idx++)
            den_is_zero[i].in[l][idx] <== den[i].out[l][idx];
    }

    component exception = IsZero();
    exception.in <== inIsInfinity + den_is_zero[0].out + den_is_zero[1].out; 
    isInfinity <== 1 - exception.out; 

    component num[2];
    for(var i=0; i<2; i++){
        num[i] = Fp2Compress(n, k, k-1, p, 3*n + 2*LOGK + 3); 
        for(var l=0; l<2; l++)for(var idx=0; idx<2*k-1; idx++)
            num[i].in[l][idx] <== x_frac[2*i][l][idx];
    }

    component x[2];
    // num / den if den != 0, else num / 1
    for(var i=0; i<2; i++){
        x[i] = SignedFp2Divide(n, k, 3*n + 2*LOGK + 3, n, p);
        for(var l=0; l<2; l++)for(var idx=0; idx<k; idx++){
            x[i].a[l][idx] <== num[i].out[l][idx];
            if(l==0 && idx==0)
                x[i].b[l][idx] <== isInfinity * (1 - den[i].out[l][idx]) + den[i].out[l][idx];
            else
                x[i].b[l][idx] <== -isInfinity * den[i].out[l][idx] + den[i].out[l][idx];
        } 
    }

    component y = Fp2Multiply(n, k, p);
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){
        y.a[i][idx] <== in[1][i][idx];
        y.b[i][idx] <== x[1].out[i][idx];
    } 

    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){
        out[0][i][idx] <== x[0].out[i][idx];
        out[1][i][idx] <== y.out[i][idx];
    }
}
