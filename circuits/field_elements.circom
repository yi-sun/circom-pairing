pragma circom 2.0.2;

include "bigint.circom";
// include "field_elements_func.circom";

// add two elements in Fp2
template Fp2Add(n, k) {
    signal input a[2][k];
    signal input b[2][k];
    signal input p[k];
    signal output c[2][k];

    component adders[2];
    for (var i = 0; i < 2; i++) {
        adders[i] = BigAddModP(n, k);
        for (var j = 0; j < k; j++) {
            adders[i].a[j] <== a[i][j];
            adders[i].b[j] <== b[i][j];
            adders[i].p[j] <== p[j];
        }   
        for (var j = 0; j < k; j ++) {
            c[i][j] <== adders[i].out[j];
        }
    }
}

// a[i][j], b[j][j] are short unsigned integers
// out[i][j] is a long unsigned integer
// basically multiply two-variable polynomials a, b
// use case: one variable will end up being 2**n; the other will be the field extension generator
template BigMultShortLong2D(n, k, l) {
    signal input a[l][k];
    signal input b[l][k];
    signal output out[2*l-1][2*k-1];

    var prod_val[2*l-1][2*k-1];
    for (var i = 0; i < 2*l-1; i++) {
        for (var j = 0; j < 2*k-1; j++) {
            prod_val[i][j] = 0;
        }
    }

    for (var i1 = 0; i1 < l; i1 ++) {
        for (var i2 = 0; i2 < l; i2 ++) {
            for (var j1 = 0; j1 < k; j1 ++) {
                for (var j2 = 0; j2 < k; j2 ++) {
                    var i = i1 + i2;
                    var j = j1 + j2;
                    prod_val[i][j] += a[i1][j1] * b[i2][j2];
                }
            }
        }
    }

    for (var i = 0; i < 2*l-1; i++) {
        for (var j = 0; j < 2*k-1; j++) {
            out[i][j] <-- prod_val[i][j];
        }
    }

    var a_poly[2*l-1][2*k-1];
    var b_poly[2*l-1][2*k-1];
    var out_poly[2*l-1][2*k-1];
    for (var i = 0; i < 2*l-1; i++) {
        for (var j = 0; j < 2*k-1; j++) {
            a_poly[i][j] = 0;
            b_poly[i][j] = 0;
            out_poly[i][j] = 0;
            for (var deg1 = 0; deg1 < l; deg1 ++) {
                for (var deg2 = 0; deg2 < k; deg2 ++) {
                    a_poly[i][j] = a_poly[i][j] + a[deg1][deg2] * (i ** deg1) * (j ** deg2);
                    b_poly[i][j] = b_poly[i][j] + b[deg1][deg2] * (i ** deg1) * (j ** deg2);
                }
            }
            for (var deg1 = 0; deg1 < 2*l-1; deg1 ++) {
                for (var deg2 = 0; deg2 < 2*k-1; deg2 ++) {
                    out_poly[i][j] = out_poly[i][j] + out[deg1][deg2] * (i ** deg1) * (j ** deg2);
                }
            }
        }
    }

    for (var i = 0; i < 2*l-1; i++) {
        for (var j = 0; j < 2*k-1; j++) {
            out_poly[i][j] === a_poly[i][j] * b_poly[i][j];
        }
    }
}

// multiplication specialized to Fp^2 
// (a0 + a1 u)*(b0 + b1 u) = (a0*b0 - a1*b1) + (a0*b1 + a1*b0)u
template Fp2multiply(n, k){
    signal input a[2][k];
    signal input b[2][k];
    signal input p[k];
    signal output out[2][k];

    var LOGK = 3;
    assert(k<7);
    assert(2*n + 1 + LOGK<254);

    // out[0] computation
    // solve for X and Y such that a0*b0 + (p-a1)*b1 = p*X + Y with Y in [0,p) 
    // -a1*b1 = (p-a1)*b1 mod p
    var a0b0_var[100] = prod(n, k, a[0], b[0]);
    var a1_neg[100] = long_sub(n, k, p, a[1]); 
    var a1b1_neg[100] = prod(n, k, a1_neg, b[1]);
    var diff[100] = long_add(n, 2*k, a0b0_var, a1b1_neg); // 2*k+1 registers
    var X_Y[2][100] = long_div2(n, k, k+1, diff, p); 
    // X = X_Y[0] has k+2 registers, Y = X_Y[1] has k registers 
    // out[0] = Y
    for(var i=0; i<k; i++)
        out[0][i] <-- X_Y[1][i];
    component range_checks[k];
    for(var i=0; i<k; i++){
        range_checks[i] = Num2Bits(n);
        range_checks[i].in <== out[0][i]; 
    }
    component lt = BigLessThan(n, k);
    for(var i=0; i<k; i++){
        lt.a[i] <== out[0][i];
        lt.b[i] <== p[i];
    }
    lt.out === 1;

    signal X[k+2];
    component X_range_checks[k+2];
    for(var i=0; i<k+2; i++){
        X[i] <-- X_Y[0][i];
        X_range_checks[i] = Num2Bits(n);
        X_range_checks[i].in <== X[i];
    }
    
    // constrain by Carry( a0 *' b0 +' p *' b1 -' a1 *' b1 - p *' X - Y ) = 0 
    // where all operations are performed without carry 
    // each register is an overflow representation in the range (-(k+1)*2^{2n+1}-2^n, (k+1)*2^{2n + 1} )
    //                                          which is inside (-2^{2n+1+LOGK}, 2^{2n+1+LOGK})
    
    component a0b0 = BigMultShortLong(n, k);
    component a1b1 = BigMultShortLong(n, k);
    component pb1 = BigMultShortLong(n, k); // 2*k-1 registers
    component pX = BigMultShortLong(n, k+2); // 2*k+3 registers
    for(var i=0; i<k; i++){
        a0b0.a[i] <== a[0][i];
        a0b0.b[i] <== b[0][i];

        a1b1.a[i] <== a[1][i];
        a1b1.b[i] <== b[1][i];

        pb1.a[i] <== p[i];
        pb1.b[i] <== b[1][i];

        pX.a[i] <== p[i];
        pX.b[i] <== X[i];
    }
    for(var i=k; i<k+2; i++){
        pX.a[i] <== 0;
        pX.b[i] <== X[i];
    }

    component carry_check = CheckCarryToZero(n, 2*n+2+LOGK, 2*k+3); 
    for(var i=0; i<k; i++)
        carry_check.in[i] <== a0b0.out[i] + pb1.out[i] - a1b1.out[i] - pX.out[i] - out[0][i]; 
    for(var i=k; i<2*k-1; i++)
        carry_check.in[i] <== a0b0.out[i] + pb1.out[i] - a1b1.out[i] - pX.out[i]; 
    for(var i=2*k-1; i<2*k+3; i++)
        carry_check.in[i] <== -pX.out[i];

    // now for out[1] computation
    // solve for Z and out[1] such that a0*b1 + a1*b0 = p*Z + out[1] with out[1] in [0,p) 
    var a0b1_var[100] = prod(n, k, a[0], b[1]);
    var a1b0_var[100] = prod(n, k, a[1], b[0]);
    var sum[100] = long_add(n, 2*k, a0b1_var, a1b0_var); // output 2*k+1 registers
    var sum_div[2][100] = long_div2(n, k, k+1, sum, p); 
    // Z = sum_div[0] has k+2 registers 
    for(var i=0; i<k; i++)
        out[1][i] <-- sum_div[1][i];
    component range_checks1[k];
    for(var i=0; i<k; i++){
        range_checks1[i] = Num2Bits(n);
        range_checks1[i].in <== out[1][i]; 
    }
    component lt1 = BigLessThan(n, k);
    for(var i=0; i<k; i++){
        lt1.a[i] <== out[1][i];
        lt1.b[i] <== p[i];
    }
    lt1.out === 1;

    signal Z[k+2];
    component Z_range_checks[k+2];
    for(var i=0; i<k+2; i++){
        Z[i] <-- sum_div[0][i];
        Z_range_checks[i] = Num2Bits(n);
        Z_range_checks[i].in <== Z[i];
    }

    // constrain by Carry( a0 *' b1 +' a1 *' b0 -' p *' Z - out[1]) = 0 
    // each register is an overflow representation in the range (-(k+1)*2^{2n}-2^n, (k+1)*2^{2n + 1} )
    //                                          which is inside (-2^{2n+1+LOGK}, 2^{2n+1+LOGK})

    component a0b1 = BigMultShortLong(n, k); // 2*k-1 registers
    component a1b0 = BigMultShortLong(n, k);
    component pZ = BigMultShortLong(n, k+2); // 2*k+3 registers
    for(var i=0; i<k; i++){
        a0b1.a[i] <== a[0][i];
        a0b1.b[i] <== b[1][i];

        a1b0.a[i] <== a[1][i];
        a1b0.b[i] <== b[0][i];
        
        pZ.a[i] <== p[i];
        pZ.b[i] <== Z[i];
    }
    for(var i=k; i<k+2; i++){
        pZ.a[i] <== 0;
        pZ.b[i] <== Z[i];
    }
    
    component carry_check1 = CheckCarryToZero(n, 2*n+2+LOGK, 2*k+3);
    for(var i=0; i<k; i++)
        carry_check1.in[i] <== a0b1.out[i] + a1b0.out[i] - pZ.out[i] - out[1][i]; 
    for(var i=k; i<2*k-1; i++)
        carry_check1.in[i] <== a0b1.out[i] + a1b0.out[i] - pZ.out[i]; 
    for(var i=2*k-1; i<2*k+3; i++)
        carry_check1.in[i] <== -pZ.out[i];
}

// input: in[0] + in[1] u
// output: (p-in[0]) + (p-in[1]) u
// assume 0 <= in < p
template Fp2negate(n, k){
    signal input in[2][k]; 
    signal input p[k];
    signal output out[2][k];
    
    component neg0 = BigSub(n, k);
    component neg1 = BigSub(n, k);
    for(var i=0; i<k; i++){
        neg0.a[i] <== p[i];
        neg1.a[i] <== p[i];
        neg0.b[i] <== in[0][i];
        neg1.b[i] <== in[1][i];
    }
    for(var i=0; i<k; i++){
        out[0][i] <== neg0.out[i];
        out[1][i] <== neg1.out[i];
    }
}

// input: a0 + a1 u, b0 + b1 u
// output: (a0-b0) + (a1-b1)u
template Fp2subtract(n, k){
    signal input a[2][k];
    signal input b[2][k];
    signal input p[k];
    signal output out[2][k];
    
    component sub0 = BigSubModP(n, k);
    component sub1 = BigSubModP(n, k);
    for(var i=0; i<k; i++){
        sub0.a[i] <== a[0][i];
        sub0.b[i] <== b[0][i];
        sub1.a[i] <== a[1][i];
        sub1.b[i] <== b[1][i];
    }
    for(var i=0; i<k; i++){
        out[0][i] <== sub0.out[i];
        out[1][i] <== sub1.out[i];
    }
}

// Src: https://github.com/paulmillr/noble-bls12-381/blob/23823d664b1767fb20c9c19c5800c66993b576a5/math.ts#L444
// We wish to find the multiplicative inverse of a nonzero
// element a + bu in Fp2. We leverage an identity
//
// (a + bu)(a - bu) = a² + b²
//
// which holds because u² = -1. This can be rewritten as
//
// (a + bu)(a - bu)/(a² + b²) = 1
//
// because a² + b² = 0 has no nonzero solutions for (a, b).
// This gives that (a - bu)/(a² + b²) is the inverse
// of (a + bu). Importantly, this can be computing using
// only a single inversion in Fp.
template Fp2invert(n, k){
    signal input in[2][k];
    signal input p[k];
    signal output out[2][k];

    // lambda = 1/(in0**2 + in1**2) % p
    
    var sq0[100] = prod(n, k, in[0], in[0]);
    var sq1[100] = prod(n, k, in[1], in[1]);
    var sq_sum[100] = long_add(n, 2*k, sq0, sq1);
    var sq_sum_div[2][100] = long_div2(n, k, k+1, sq_sum, p);
    // lambda = 1/(sq_sum)%p
    var lambda[100] = mod_inv(n, k, sq_sum_div[1], p);
    var out0[100] = prod(n, k, lambda, in[0]);
    var out0_div[2][100] = long_div(n, k, out0, p);
    for(var i=0; i<k; i++)
        out[0][i] <-- out0_div[1][i];
    
    var out1_pre[100] = long_sub(n, k, p, in[1]);
    var out1[100] = prod(n, k, lambda, out1_pre);
    var out1_div[2][100] = long_div(n, k, out1, p);
    for(var i=0; i<k; i++)
        out[1][i] <-- out1_div[1][i];

    //range checks
    component outRangeChecks[2][k];
    for(var i=0; i<2; i++) for(var j=0; j<k; j++){
        outRangeChecks[i][j] = Num2Bits(n);
        outRangeChecks[i][j].in <== out[i][j];
    }

    component in_out = Fp2multiply(n, k);
    for(var i=0; i<2; i++)for(var j=0; j<k; j++){
        in_out.a[i][j] <-- in[i][j];
        in_out.b[i][j] <-- out[i][j];
    }
    for(var i=0; i<k; i++) in_out.p[i] <-- p[i];

    for(var i=0; i<2; i++)for(var j=0; j<k; j++){
        if(i == 0 && j == 0)
            in_out.out[i][j] === 1;
        else
            in_out.out[i][j] === 0;
    }
}

// input: a+b u
// output: a-b u 
// IF p = 3 mod 4 THEN a - b u = (a+b u)^p <-- Frobenius map 
// aka Fp2frobeniusMap(n, k)
template Fp2conjugate(n, k){
    signal input in[2][k]; 
    signal input p[k];
    signal output out[2][k];
    
    component neg1 = BigSub(n, k);
    for(var i=0; i<k; i++){
        neg1.a[i] <== p[i];
        neg1.b[i] <== in[1][i];
    }
    for(var i=0; i<k; i++){
        out[0][i] <== in[0][i];
        out[1][i] <== neg1.out[i];
    }
}

// raises to q^power-th power 
template Fp2frobeniusMap(n, k, power){
    signal input in[2][k];
    signal input p[k];
    signal output out[2][k];
    
    var pow = power % 2;
    component neg1 = BigSub(n,k);
    if(pow == 0){
        for(var i=0; i<k; i++){
            out[0][i] <== in[0][i];
            out[1][i] <== in[1][i];
        }
    }else{
        for(var i=0; i<k; i++){
            neg1.a[i] <== p[i];
            neg1.b[i] <== in[1][i];
        }
        for(var i=0; i<k; i++){
            out[0][i] <== in[0][i];
            out[1][i] <== neg1.out[i];
        }
    }
}

// template Fp12frobeniusMap(n, k, power){
//     signal input in[6][2][k];
//     signal input p[k];
//     signal output out[6][2][k];

//     var FP12_FROBENIUS_COEFFICIENTS[12][6][2][k] = get_Fp12_frobenius(n, k);
//     var pow = power % 12;
    
//     // apply Frob to coefficients first
//     component in_frob[6]; 
//     for(var i=0; i<6; i++){
//         in_frob[i] = Fp2frobeniusMap(n, k, pow); 
//         for(var j=0; j<k; j++){
//             in_frob[i].in[0][j] <== in[i][0][j];
//             in_frob[i].in[1][j] <== in[i][1][j];
//             in_frob[i].p[j] <== p[j];
//         }
//     }
    
//     // multiply in_frob[i] by FP12_FROBENIUS_COEFFICIENTS[pow][i] 
//     // if pow is even, then FP12_FROBENIUS_COEFFICIENTS[pow][i] is in Fp instead of Fp2, so can optimize 
//     component mult_odd[6];
//     component mult_even[6][2];
//     if( (pow % 2) == 0 ){
//         for(var i=0; i<6; i++){
//             mult_even[i][0] = BigMultModP(n, k);
//             mult_even[i][1] = BigMultModP(n, k);
//             for(var j=0; j<k; j++){
//                 mult_even[i][0].a[j] <== in_frob[i].out[0][j];
//                 mult_even[i][1].a[j] <== in_frob[i].out[1][j];

//                 mult_even[i][0].b[j] <== FP12_FROBENIUS_COEFFICIENTS[pow][i][0][j];
//                 mult_even[i][1].b[j] <== FP12_FROBENIUS_COEFFICIENTS[pow][i][0][j];
                
//                 mult_even[i][0].p[j] <== p[j];
//                 mult_even[i][1].p[j] <== p[j];
//             }
//             for(var j=0; j<k; j++){
//                 out[i][0][j] <== mult_even[i][0].out[j];
//                 out[i][1][j] <== mult_even[i][1].out[j];
//             }
//         }
//     }else{
//         for(var i=0; i<6; i++){
//             mult_odd[i] = Fp2multiply(n, k);
//             for(var j=0; j<k; j++){
//                 for(var eps=0; eps<2; eps++){
//                     mult_odd[i].a[eps][j] <== in_frob[i].out[eps][j];
//                     mult_odd[i].b[eps][j] <== FP12_FROBENIUS_COEFFICIENTS[pow][i][eps][j];
//                 }
//                 mult_odd[i].p[j] <== p[j];
//             }
//             for(var j=0; j<k; j++){
//                 out[i][0][j] <== mult_odd[i].out[0][j];
//                 out[i][1][j] <== mult_odd[i].out[1][j];
//             }
//         }
//     }
// }

template Fp12Add(n, k) {
    signal input a[6][2][k];
    signal input b[6][2][k];
    signal input p[k];
    signal output c[6][2][k];
    component adders[6][2];
    for (var i = 0; i < 6; i ++) {
        for (var j = 0; j < 2; j ++) {
            adders[i][j] = BigAddModP(n,k);
            for (var m = 0; m < k; m ++) {
                adders[i][j].a[m] <== a[i][j][m];
                adders[i][j].b[m] <== b[i][j][m];
                adders[i][j].p[m] <== p[m];
            }
            for (var m = 0; m < k; m ++) {
                c[i][j][m] <== adders[i][j].out[m];
            }
        }
    }
}

// a = sum w^i u^j a_ij for w^6=u+1, u^2=-1. similarly for b
// we first write a = A + Bu, b = C + Du and compute ab = (AC + B(p-D)) + (AD+BC)u
// with deg(w) = 10, deg(u) = 1 and then simplify the representation
// first using w^6 = u + 1 to get deg(w) = 5, deg (u) = 2
// in addition to computing AD+BC we also compute A(p-D) + B(p-C) to avoid subtraction
// and then using u^2 = -1 to get deg(w) = 5, deg(u) = 1
template Fp12Multiply(n, k) {
    var l = 6;
    signal input a[l][2][k];
    signal input b[l][2][k];
    signal input p[k];
    signal output out[l][2][k];


    var LOGK = 3;
    var LOGL = 4;
    assert(l<15);
    assert(k<7);
    assert(2*n + 1 + LOGK + LOGL <254);

    var a0[l][k];
    var a1[l][k];
    var b0[l][k];
    var b1[l][k];
    var neg_b0[l][100];
    var neg_b1[l][100];
    for (var i = 0; i < l; i ++) { 
        for ( var j = 0; j < k; j ++) {
            a0[i][j] = a[i][0][j];
            a1[i][j] = a[i][1][j];
            b0[i][j] = b[i][0][j];
            b1[i][j] = b[i][1][j];
        }
    }
    for ( var i = 0; i < l; i ++) {
        neg_b0[i] = long_sub(n, k, p, b0[i]);
        neg_b1[i] = long_sub(n, k, p, b1[i]);
    }

    var real_init[2*l-1][100];
    var imag_init[2*l-1][100];
    var imag_init_neg[2*l-1][100];
    var real[l][2][100];
    var imag[l][2][100];
    // each product will be 2l-1 x 2k
    var a0b0_var[100][100] = prod2D(n, k, l, a0, b0);
    var a1b1_neg[100][100] = prod2D(n, k, l, a1, neg_b1);
    var a0b1_var[100][100] = prod2D(n, k, l, a0, b1);
    var a1b0_var[100][100] = prod2D(n, k, l, a1, b0);
    var a0b1_neg[100][100] = prod2D(n, k, l, a0, neg_b1);
    var a1b0_neg[100][100] = prod2D(n, k, l, a1, neg_b0);
    for (var i = 0; i < 2*l - 1; i ++) { // compute initial rep (deg w = 10)
        real_init[i] = long_add(n, 2*k, a0b0_var[i], a1b1_neg[i]); // 2*k+1 registers each
        imag_init[i] = long_add(n, 2*k, a0b1_var[i], a1b0_var[i]);
        imag_init_neg[i] = long_add(n, 2*k, a0b1_neg[i], a1b0_neg[i]);
    }
    var real_carry[l][100];
    var imag_carry[l][100];
    var real_final[l][100];
    var imag_final[l][100];
    var zeros[100]; // to balance register sizes
    for (var i = 0; i < 100; i ++) {
        zeros[i] = 0;
    }
    for (var i = 0; i < l; i ++) {
        if (i == l - 1) {
            real_carry[i] = long_add(n, 2*k+1, zeros, zeros);
            imag_carry[i] = long_add(n, 2*k+1, zeros, zeros);
        } else {
            real_carry[i] = long_add(n, 2*k+1, real_init[i+l], imag_init_neg[i+l]); // now 2*k+2 registers
            imag_carry[i] = long_add(n, 2*k+1, imag_init[i+l], real_init[i+l]);
        }
    }
    for (var i = 0; i < l; i ++) {
        real_final[i] = long_add_unequal(n, 2*k+2, 2*k+1, real_carry[i], real_init[i]); // now 2*k+3 registers
        imag_final[i] = long_add_unequal(n, 2*k+2, 2*k+1, imag_carry[i], imag_init[i]);
    }
    var XYreal_temp[l][2][100];
    var XYimag_temp[l][2][100];
    signal XYreal[l][2][k+4];
    signal XYimag[l][2][k+4];
    for (var i = 0; i < l; i ++) {
        XYreal_temp[i] = long_div2(n, k, k+3, real_final[i], p); // k+4 register quotient, k register remainder
        XYimag_temp[i] = long_div2(n, k, k+3, imag_final[i], p);
    }
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < k+4; j ++) {
            XYreal[i][0][j] <-- XYreal_temp[i][0][j];
            XYimag[i][0][j] <-- XYimag_temp[i][0][j];
            if (j < k) {
                XYreal[i][1][j] <-- XYreal_temp[i][1][j];
                XYimag[i][1][j] <-- XYimag_temp[i][1][j];
            } else {
                XYreal[i][1][j] <== 0;
                XYimag[i][1][j] <== 0;
            }
        }
    }
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < k; j ++) {
            out[i][0][j] <== XYreal[i][1][j];
            out[i][1][j] <== XYimag[i][1][j];
        }
    }

    component range_checks[l][2][k];
    for(var i=0; i<l; i++){
        for (var j = 0; j < 2; j ++) {
            for (var m = 0; m < k; m ++) {
                range_checks[i][j][m] = Num2Bits(n);
                range_checks[i][j][m].in <== out[i][j][m];
            }
        }
    }
    component lt[l][2];
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < 2; j ++) {
            lt[i][j] = BigLessThan(n, k);
            for (var m = 0; m < k; m ++) {
                lt[i][j].a[m] <== out[i][j][m];
                lt[i][j].b[m] <== p[m];
            }
            lt[i][j].out === 1;
        }
    }

    component X_range_checks[l][2][k+4];
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < k+4; j ++) {
            X_range_checks[i][0][j] = Num2Bits(n);
            X_range_checks[i][1][j] = Num2Bits(n);
            X_range_checks[i][0][j].in <== XYreal[i][0][j];
            X_range_checks[i][1][j].in <== XYimag[i][0][j];
        }
    }

    // constrain by: 
    // Carry( lower(a0 *' b0 +' p *' a1 -' a1 *' b1) + upper(a0 *' b0 +' p *' a1 -' a1 *' b1) 
    //       + upper(a0 *' p -' a0 *' b1 + a1 *' p -' a1 *' b0) - p *' Xreal - Yreal ) = 0 
    // Carry( lower(a0 *' p -' a0 *' b1 + a1 *' p -' a1 *' b0) + upper(a0 *' p -' a0 *' b1 + a1 *' p -' a1 *' b0)
    //      + upper(a0 *' b0 +' p *' a1 -' a1 *' b1) - p *' Ximag - Yimag) = 0
    // where all operations are performed without carry 
    // lower() refers to the coefficients with deg(w) <= 5
    // upper() refers to the coefficients with deg(w) >= 6
    // each register is an overflow representation in the range (-kl*2^{2n+4}, kl*2^{2n + 4} )
    
    component a0b0 = BigMultShortLong2D(n, k, l);
    component a1b1 = BigMultShortLong2D(n, k, l);
    component pa1 = BigMultShortLong2D(n, k, l);
    component pa0 = BigMultShortLong2D(n, k, l);
    component a1b0 = BigMultShortLong2D(n, k, l);
    component a0b1 = BigMultShortLong2D(n, k, l);
    component pXreal[l];
    component pXimag[l];
    for (var i = 0; i < l; i ++) {
        pXreal[i] = BigMultShortLong(n, k+4);
        pXimag[i] = BigMultShortLong(n, k+4);
        for (var j = 0; j < k; j ++) {
            a0b0.a[i][j] <== a[i][0][j];
            a0b0.b[i][j] <== b[i][0][j];

            a1b1.a[i][j] <== a[i][1][j];
            a1b1.b[i][j] <== b[i][1][j];

            if (i == 0) {
                pa1.a[i][j] <== p[j];
                pa0.a[i][j] <== p[j];
            } else {
                pa1.a[i][j] <== 0;
                pa0.a[i][j] <== 0;
            }

            pa1.b[i][j] <== a[i][1][j];
            pa0.b[i][j] <== a[i][0][j];

            a1b0.a[i][j] <== a[i][1][j];
            a1b0.b[i][j] <== b[i][0][j];

            a0b1.a[i][j] <== a[i][0][j];
            a0b1.b[i][j] <== b[i][1][j];

            pXreal[i].a[j] <== p[j];
            pXreal[i].b[j] <== XYreal[i][0][j];

            pXimag[i].a[j] <== p[j];
            pXimag[i].b[j] <== XYimag[i][0][j];
        }
        for (var j = k; j < k+4; j ++) {
            pXreal[i].a[j] <== 0;
            pXreal[i].b[j] <== XYreal[i][0][j];
            pXimag[i].a[j] <== 0;
            pXimag[i].b[j] <== XYimag[i][0][j];
        }
    }


    component carry_check[l][2];
    // Carry( lower(a0 *' b0 +' p *' a1 -' a1 *' b1) + upper(a0 *' b0 +' p *' a1 -' a1 *' b1) 
    //       + upper(a0 *' p -' a0 *' b1 + a1 *' p -' a1 *' b0) - p *' Xreal - Yreal ) = 0 
    // Carry( lower(a0 *' b1 + a1 *' b0) + upper(a0 *' b1 + a1 *' b0)
    //      + upper(a0 *' b0 +' p *' a1 -' a1 *' b1) - p *' Ximag - Yimag) = 0
    for (var i = 0; i < l-1; i ++) {
        carry_check[i][0] = CheckCarryToZero(n, 2*n+4+LOGK+LOGL, 2*k+4);
        carry_check[i][1] = CheckCarryToZero(n, 2*n+4+LOGK+LOGL, 2*k+4);
        for (var j = 0; j < k; j ++) {
            carry_check[i][0].in[j] <== (a0b0.out[i][j] + pa1.out[i][j] - a1b1.out[i][j]) + (a0b0.out[i+l][j] + pa1.out[i+l][j] - a1b1.out[i+l][j]) + (pa0.out[i+l][j] - a0b1.out[i+l][j] + pa1.out[i+l][j] - a1b0.out[i+l][j]) - pXreal[i].out[j] - out[i][0][j];
            carry_check[i][1].in[j] <== (a0b1.out[i][j] + a1b0.out[i][j]) + (a0b1.out[i+l][j] + a1b0.out[i+l][j]) + (a0b0.out[i+l][j] + pa1.out[i+l][j] - a1b1.out[i+l][j]) - pXimag[i].out[j] - out[i][1][j];
        }
        for (var j = k; j < 2*k-1; j ++) {
            carry_check[i][0].in[j] <== (a0b0.out[i][j] + pa1.out[i][j] - a1b1.out[i][j]) + (a0b0.out[i+l][j] + pa1.out[i+l][j] - a1b1.out[i+l][j]) + (pa0.out[i+l][j] - a0b1.out[i+l][j] + pa1.out[i+l][j] - a1b0.out[i+l][j]) - pXreal[i].out[j];
            carry_check[i][1].in[j] <== (a0b1.out[i][j] + a1b0.out[i][j]) + (a0b1.out[i+l][j] + a1b0.out[i+l][j]) + (a0b0.out[i+l][j] + pa1.out[i+l][j] - a1b1.out[i+l][j]) - pXimag[i].out[j];
        }
        for (var j = 2*k-1; j < 2*k+4; j ++) {
            carry_check[i][0].in[j] <== - pXreal[i].out[j];
            carry_check[i][1].in[j] <== - pXimag[i].out[j];
        }
    }

    carry_check[l-1][0] = CheckCarryToZero(n, 2*n+4+LOGK+LOGL, 2*k+4);
    carry_check[l-1][1] = CheckCarryToZero(n, 2*n+4+LOGK+LOGL, 2*k+4);
    // do the last part, less carries
    for (var j = 0; j < k; j ++) {
        carry_check[l-1][0].in[j] <== (a0b0.out[l-1][j] + pa1.out[l-1][j] - a1b1.out[l-1][j]) - pXreal[l-1].out[j] - out[l-1][0][j];
        carry_check[l-1][1].in[j] <== (a0b1.out[l-1][j] + a1b0.out[l-1][j]) - pXimag[l-1].out[j] - out[l-1][1][j];
    }
    for (var j = k; j < 2*k-1; j ++) {
        carry_check[l-1][0].in[j] <== (a0b0.out[l-1][j] + pa1.out[l-1][j] - a1b1.out[l-1][j]) - pXreal[l-1].out[j];
        carry_check[l-1][1].in[j] <== (a0b1.out[l-1][j] + a1b0.out[l-1][j]) - pXimag[l-1].out[j];
    }
    for (var j = 2*k-1; j < 2*k+4; j ++) {
        carry_check[l-1][0].in[j] <== - pXreal[l-1].out[j];
        carry_check[l-1][1].in[j] <== - pXimag[l-1].out[j];
    }
}
