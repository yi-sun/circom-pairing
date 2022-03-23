pragma circom 2.0.3;

include "bigint.circom";
include "field_elements_func.circom";
include "fp.circom";
include "fp2.circom";
include "bls12-381_func.circom";

template Fp12FrobeniusMap(n, k, power){
    signal input in[6][2][k];
    signal output out[6][2][k];

    var p[50] = get_BLS12_381_prime(n, k);
    var FP12_FROBENIUS_COEFFICIENTS[12][6][2][20] = get_Fp12_frobenius(n, k);
    var pow = power % 12;
 
    component in_frob[6]; 
 
    // multiply in_frob[i] by FP12_FROBENIUS_COEFFICIENTS[pow][i] 
    // if pow is even, then FP12_FROBENIUS_COEFFICIENTS[pow][i] is in Fp instead of Fp2, so can optimize 
    component mult_odd[6];
    component mult_even[6][2];
    if( (pow % 2) == 0 ){
        for(var j=0; j<k; j++){
            out[0][0][j] <== in[0][0][j];
            out[0][1][j] <== in[0][1][j];
        } 
        for(var i=1; i<6; i++){
            mult_even[i][0] = FpMultiply(n, k, p);
            mult_even[i][1] = FpMultiply(n, k, p);
            for(var j=0; j<k; j++){
                mult_even[i][0].a[j] <== in[i][0][j];
                mult_even[i][1].a[j] <== in[i][1][j];

                mult_even[i][0].b[j] <== FP12_FROBENIUS_COEFFICIENTS[pow][i][0][j];
                mult_even[i][1].b[j] <== FP12_FROBENIUS_COEFFICIENTS[pow][i][0][j];
            }
            for(var j=0; j<k; j++){
                out[i][0][j] <== mult_even[i][0].out[j];
                out[i][1][j] <== mult_even[i][1].out[j];
            }
        }
    }else{
        // apply Frob to coefficients first
        for(var i=0; i<6; i++){
            in_frob[i] = Fp2FrobeniusMap(n, k, pow, p); 
            for(var j=0; j<k; j++){
                in_frob[i].in[0][j] <== in[i][0][j];
                in_frob[i].in[1][j] <== in[i][1][j];
            }
        }
        for(var j=0; j<k; j++){
            out[0][0][j] <== in_frob[0].out[0][j];
            out[0][1][j] <== in_frob[0].out[1][j];
        } 
        for(var i=1; i<6; i++){
            mult_odd[i] = Fp2Multiply(n, k, p);
            for(var j=0; j<k; j++){
                for(var eps=0; eps<2; eps++){
                    mult_odd[i].a[eps][j] <== in_frob[i].out[eps][j];
                    mult_odd[i].b[eps][j] <== FP12_FROBENIUS_COEFFICIENTS[pow][i][eps][j];
                }
            }
            for(var j=0; j<k; j++){
                out[i][0][j] <== mult_odd[i].out[0][j];
                out[i][1][j] <== mult_odd[i].out[1][j];
            }
        }
    }
}

template Fp12Add(n, k, p) {
    signal input a[6][2][k];
    signal input b[6][2][k];
    signal output out[6][2][k];
    component adders[6][2];
    for (var i = 0; i < 6; i ++) {
        for (var j = 0; j < 2; j ++) {
            adders[i][j] = FpAdd(n,k, p);
            for (var m = 0; m < k; m ++) {
                adders[i][j].a[m] <== a[i][j][m];
                adders[i][j].b[m] <== b[i][j][m];
                adders[i][j].p[m] <== p[m];
            }
            for (var m = 0; m < k; m ++) {
                out[i][j][m] <== adders[i][j].out[m];
            }
        }
    }
}

// a is 2 x k array representing element a0 - a1 of Fp where we keep track of negatives
// b is 6 x 4 x k array representing element (b0 - b2) + (b1 - b3) u of Fp12 keeping track of negatives
//      where b_i = b[][i][] is 6 x k array
// out is a*b in Fp12 as 6 x 4 x (2k-1) array
template Fp12ScalarMultiplyNoCarry(n, k, m_out){
    signal input a[2][k];
    signal input b[6][4][k];
    signal output out[6][4][2*k-1];

    component ab[6][2]; 
    for(var i=0; i<6; i++)for(var j=0; j<2; j++){
        ab[i][j] = BigMultNoCarry(n, k, m_out); // 2k-1 registers 

        for(var eps=0; eps<2; eps++)for(var idx=0; idx<k; idx++){
            ab[i][j].a[eps][idx] <== a[eps][idx];
            ab[i][j].b[eps][idx] <== b[i][j+2*eps][idx]; 
        } 
    }
    
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)
        for(var eps=0; eps<2; eps++)for(var idx=0; idx<2*k-1; idx++)
            out[i][j+2*eps][idx] <== ab[i][j].out[eps][idx];
}


template Fp12ScalarMultiplyNoCarryUnequal(n, ka, kb, m_out){
    signal input a[2][ka];
    signal input b[6][4][kb];
    signal output out[6][4][ka+kb-1];

    component ab[6][2]; 
    for(var i=0; i<6; i++)for(var j=0; j<2; j++){
        ab[i][j] = BigMultNoCarryUnequal(n, ka, kb, m_out); // 2k-1 registers 

        for(var eps=0; eps<2; eps++)for(var idx=0; idx<ka; idx++)
            ab[i][j].a[eps][idx] <== a[eps][idx];
        for(var eps=0; eps<2; eps++)for(var idx=0; idx<kb; idx++)
            ab[i][j].b[eps][idx] <== b[i][j+2*eps][idx]; 
    }
    
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)
        for(var eps=0; eps<2; eps++)for(var idx=0; idx<ka+kb-1; idx++)
            out[i][j+2*eps][idx] <== ab[i][j].out[eps][idx];
}


// a = sum w^i u^j a_ij for w^6=u+1, u^2=-1. similarly for b
// we first write a = A + B u, b = C + D u and compute 
// ab = (AC + B(p-D)) + (AD+BC) u, and then simplify the representation
template Fp12Multiply(n, k, p) {
    var l = 6;
    signal input a[l][2][k];
    signal input b[l][2][k];
    signal output out[l][2][k];


    var LOGK = log_ceil(k);
    var LOGL = log_ceil(l);
    assert(2*n + 1 + LOGK + LOGL <254);

    var a0[l][k];
    var a1[l][k];
    var b0[l][k];
    var b1[l][k];
    var neg_b0[l][50];
    var neg_b1[l][50];
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

    var real_init[2*l-1][50];
    var imag_init[2*l-1][50];
    var imag_init_neg[2*l-1][50];
    var real[l][2][50];
    var imag[l][2][50];
    // each product will be 2l-1 x 2k
    var a0b0_var[20][50] = prod2D(n, k, l, a0, b0);
    var a1b1_neg[20][50] = prod2D(n, k, l, a1, neg_b1);
    var a0b1_var[20][50] = prod2D(n, k, l, a0, b1);
    var a1b0_var[20][50] = prod2D(n, k, l, a1, b0);
    var a0b1_neg[20][50] = prod2D(n, k, l, a0, neg_b1);
    var a1b0_neg[20][50] = prod2D(n, k, l, a1, neg_b0);
    for (var i = 0; i < 2*l - 1; i ++) { // compute initial rep (deg w = 10)
        real_init[i] = long_add(n, 2*k, a0b0_var[i], a1b1_neg[i]); // 2*k+1 registers each
        imag_init[i] = long_add(n, 2*k, a0b1_var[i], a1b0_var[i]);
        imag_init_neg[i] = long_add(n, 2*k, a0b1_neg[i], a1b0_neg[i]);
    }
    var real_carry[l][50];
    var imag_carry[l][50];
    var real_final[l][50];
    var imag_final[l][50];
    var zeros[50]; // to balance register sizes
    for (var i = 0; i < 50; i ++) {
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
    var XYreal_temp[l][2][50];
    var XYimag_temp[l][2][50];
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
        pXreal[i] = BigMultShortLong(n, k+4, 2*n+4+LOGK+LOGL);
        pXimag[i] = BigMultShortLong(n, k+4, 2*n+4+LOGK+LOGL);
        for (var j = 0; j < k; j ++) {
            a0b0.a[i][j] <== a[i][0][j];
            a0b0.b[i][j] <== b[i][0][j];

            a1b1.a[i][j] <== a[i][1][j];
            a1b1.b[i][j] <== b[i][1][j];

            pa1.a[i][j] <== p[j];
            pa0.a[i][j] <== p[j];

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

// we first write a = (a0 - a2) + (a1 - a3) u, b = (b0 - b2) + (b1 - b3) u for ai, bi being:
//     * length 6 vectors with k registers in [0, B_a) and [0, B_b)
// ab = (a0 b0 + a2 b2 + a1 b3 + a3 b1 - a0 b2 - a2 b0 - a1 b1 - a3 b3) + (a0 b1 + a2 b3 + a1 b0 + a3 b2 - a0 b3 - a2 b1 - a1 b2 - a3 b0) u
// set X = a0 b0 + a2 b2 + a1 b3 + a3 b1, Z = a0 b2 + a2 b0 + a1 b1 + a3 b3
//     Y = a0 b1 + a2 b3 + a1 b0 + a3 b2, W = a0 b3 + a2 b1 + a1 b2 + a3 b0 u
// Applying w^6 = u + 1 and splitting X, Y, Z, W into X_0, X_1 for w^0, ..., w^5 and w^6, ..., w^11 coeffs, get:
//     X_0 + X_1 + W_1 - Z_0 - Z_1 - Y_1 + (Y_0 + Y_1 + X_1 - W_0 - W_1 - Z_1) u
// The real and imaginary parts are
//     * length 6 vectors with 2k-1 registers in [0, B_a * B_b * 12 * k)
template Fp12MultiplyNoCarry(n, k, m_out){
    var l = 6;
    signal input a[l][4][k];
    signal input b[l][4][k];
    signal output out[l][4][2*k-1];

    component a0b0 = BigMultShortLong2D(n, k, l);
    component a0b1 = BigMultShortLong2D(n, k, l);
    component a0b2 = BigMultShortLong2D(n, k, l);
    component a0b3 = BigMultShortLong2D(n, k, l);
    component a1b0 = BigMultShortLong2D(n, k, l);
    component a1b1 = BigMultShortLong2D(n, k, l);
    component a1b2 = BigMultShortLong2D(n, k, l);
    component a1b3 = BigMultShortLong2D(n, k, l);
    component a2b0 = BigMultShortLong2D(n, k, l);
    component a2b1 = BigMultShortLong2D(n, k, l);
    component a2b2 = BigMultShortLong2D(n, k, l);
    component a2b3 = BigMultShortLong2D(n, k, l);
    component a3b0 = BigMultShortLong2D(n, k, l);
    component a3b1 = BigMultShortLong2D(n, k, l);
    component a3b2 = BigMultShortLong2D(n, k, l);
    component a3b3 = BigMultShortLong2D(n, k, l);
    
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < k; j ++) {
	    a0b0.a[i][j] <== a[i][0][j];
	    a0b1.a[i][j] <== a[i][0][j];
	    a0b2.a[i][j] <== a[i][0][j];
	    a0b3.a[i][j] <== a[i][0][j];

	    a1b0.a[i][j] <== a[i][1][j];
	    a1b1.a[i][j] <== a[i][1][j];
	    a1b2.a[i][j] <== a[i][1][j];
	    a1b3.a[i][j] <== a[i][1][j];

	    a2b0.a[i][j] <== a[i][2][j];
	    a2b1.a[i][j] <== a[i][2][j];
	    a2b2.a[i][j] <== a[i][2][j];
	    a2b3.a[i][j] <== a[i][2][j];

	    a3b0.a[i][j] <== a[i][3][j];
	    a3b1.a[i][j] <== a[i][3][j];
	    a3b2.a[i][j] <== a[i][3][j];
	    a3b3.a[i][j] <== a[i][3][j];

	    a0b0.b[i][j] <== b[i][0][j];
	    a1b0.b[i][j] <== b[i][0][j];
	    a2b0.b[i][j] <== b[i][0][j];
	    a3b0.b[i][j] <== b[i][0][j];

	    a0b1.b[i][j] <== b[i][1][j];
	    a1b1.b[i][j] <== b[i][1][j];
	    a2b1.b[i][j] <== b[i][1][j];
	    a3b1.b[i][j] <== b[i][1][j];

	    a0b2.b[i][j] <== b[i][2][j];
	    a1b2.b[i][j] <== b[i][2][j];
	    a2b2.b[i][j] <== b[i][2][j];
	    a3b2.b[i][j] <== b[i][2][j];

	    a0b3.b[i][j] <== b[i][3][j];
	    a1b3.b[i][j] <== b[i][3][j];
	    a2b3.b[i][j] <== b[i][3][j];
	    a3b3.b[i][j] <== b[i][3][j];
	}
    }

    signal X[2 * l - 1][2 * k - 1];
    signal Y[2 * l - 1][2 * k - 1];
    signal Z[2 * l - 1][2 * k - 1];
    signal W[2 * l - 1][2 * k - 1];
    for (var i = 0; i < 2 * l - 1; i++) {
	for (var j = 0; j < 2 * k - 1; j++) {
	    X[i][j] <== a0b0.out[i][j] + a2b2.out[i][j] + a1b3.out[i][j] + a3b1.out[i][j];
	    Z[i][j] <== a0b2.out[i][j] + a2b0.out[i][j] + a1b1.out[i][j] + a3b3.out[i][j];
	    Y[i][j] <== a0b1.out[i][j] + a2b3.out[i][j] + a1b0.out[i][j] + a3b2.out[i][j];
	    W[i][j] <== a0b3.out[i][j] + a2b1.out[i][j] + a1b2.out[i][j] + a3b0.out[i][j];	   
	}
    }

    for (var i = 0; i < l; i++) {
        for (var j = 0; j < 2 * k - 1; j++) {
            if (i < l - 1) {
                out[i][0][j] <== X[i][j] + X[l + i][j] + W[l + i][j];
                out[i][1][j] <== Y[i][j] + Y[l + i][j] + X[l + i][j];
                out[i][2][j] <== Z[i][j] + Z[l + i][j] + Y[l + i][j];
                out[i][3][j] <== W[i][j] + W[l + i][j] + Z[l + i][j];
            } else {
                out[i][0][j] <== X[i][j];
                out[i][1][j] <== Y[i][j];
                out[i][2][j] <== Z[i][j];
                out[i][3][j] <== W[i][j];
            }
        }
    }
    component range_checks[l][4][2*k-1];
    for (var outer = 0; outer < l; outer ++) {
        for (var i = 0; i < 4; i ++) {
            for (var j = 0; j < 2*k-1; j ++) {
                range_checks[outer][i][j] = Num2Bits(m_out);
                range_checks[outer][i][j].in <== out[outer][i][j];
            }
        }
    }
}

// The real and imaginary parts are
//     * length 6 vectors with 2k-1 registers in [0, B_a * B_b * 12 * min(ka, kb) )
template Fp12MultiplyNoCarryUnequal(n, ka, kb, m_out){
    var l = 6;
    signal input a[l][4][ka];
    signal input b[l][4][kb];
    signal output out[l][4][ka + kb -1];

    component a0b0 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a0b1 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a0b2 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a0b3 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a1b0 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a1b1 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a1b2 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a1b3 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a2b0 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a2b1 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a2b2 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a2b3 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a3b0 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a3b1 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a3b2 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    component a3b3 = BigMultShortLong2DUnequal(n, ka, kb, l, l);
    
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < ka; j ++) {
            a0b0.a[i][j] <== a[i][0][j];
            a0b1.a[i][j] <== a[i][0][j];
            a0b2.a[i][j] <== a[i][0][j];
            a0b3.a[i][j] <== a[i][0][j];

            a1b0.a[i][j] <== a[i][1][j];
            a1b1.a[i][j] <== a[i][1][j];
            a1b2.a[i][j] <== a[i][1][j];
            a1b3.a[i][j] <== a[i][1][j];

            a2b0.a[i][j] <== a[i][2][j];
            a2b1.a[i][j] <== a[i][2][j];
            a2b2.a[i][j] <== a[i][2][j];
            a2b3.a[i][j] <== a[i][2][j];

            a3b0.a[i][j] <== a[i][3][j];
            a3b1.a[i][j] <== a[i][3][j];
            a3b2.a[i][j] <== a[i][3][j];
            a3b3.a[i][j] <== a[i][3][j];
        }
        for (var j = 0; j < kb; j ++) {
            a0b0.b[i][j] <== b[i][0][j];
            a1b0.b[i][j] <== b[i][0][j];
            a2b0.b[i][j] <== b[i][0][j];
            a3b0.b[i][j] <== b[i][0][j];

            a0b1.b[i][j] <== b[i][1][j];
            a1b1.b[i][j] <== b[i][1][j];
            a2b1.b[i][j] <== b[i][1][j];
            a3b1.b[i][j] <== b[i][1][j];

            a0b2.b[i][j] <== b[i][2][j];
            a1b2.b[i][j] <== b[i][2][j];
            a2b2.b[i][j] <== b[i][2][j];
            a3b2.b[i][j] <== b[i][2][j];

            a0b3.b[i][j] <== b[i][3][j];
            a1b3.b[i][j] <== b[i][3][j];
            a2b3.b[i][j] <== b[i][3][j];
            a3b3.b[i][j] <== b[i][3][j];
        }
	}
    

    signal X[2 * l - 1][ka + kb - 1];
    signal Y[2 * l - 1][ka + kb - 1];
    signal Z[2 * l - 1][ka + kb - 1];
    signal W[2 * l - 1][ka + kb - 1];
    for (var i = 0; i < 2 * l - 1; i++) {
	for (var j = 0; j < ka + kb - 1; j++) {
	    X[i][j] <== a0b0.out[i][j] + a2b2.out[i][j] + a1b3.out[i][j] + a3b1.out[i][j];
	    Z[i][j] <== a0b2.out[i][j] + a2b0.out[i][j] + a1b1.out[i][j] + a3b3.out[i][j];
	    Y[i][j] <== a0b1.out[i][j] + a2b3.out[i][j] + a1b0.out[i][j] + a3b2.out[i][j];
	    W[i][j] <== a0b3.out[i][j] + a2b1.out[i][j] + a1b2.out[i][j] + a3b0.out[i][j];	   
	}
    }

    for (var i = 0; i < l; i++) {
        for (var j = 0; j < ka + kb - 1; j++) {
            if (i < l - 1) {
                out[i][0][j] <== X[i][j] + X[l + i][j] + W[l + i][j];
                out[i][1][j] <== Y[i][j] + Y[l + i][j] + X[l + i][j];
                out[i][2][j] <== Z[i][j] + Z[l + i][j] + Y[l + i][j];
                out[i][3][j] <== W[i][j] + W[l + i][j] + Z[l + i][j];
            } else {
                out[i][0][j] <== X[i][j];
                out[i][1][j] <== Y[i][j];
                out[i][2][j] <== Z[i][j];
                out[i][3][j] <== W[i][j];
            }
        }
    }
    component range_checks[l][4][ka+kb-1];
    for (var outer = 0; outer < l; outer ++) {
        for (var i = 0; i < 4; i ++) {
            for (var j = 0; j < ka+kb-1; j ++) {
                range_checks[outer][i][j] = Num2Bits(m_out);
                range_checks[outer][i][j].in <== out[outer][i][j];
            }
        }
    }
}

template Fp12Compress(n, k, m, p, m_out){
    var l = 6;
    signal input in[l][4][k+m];
    signal output out[l][4][k];

    component reduce[l][4];
    for (var i = 0; i < l; i++) {
        for (var j = 0; j < 4; j++){
            reduce[i][j] = PrimeReduce(n, k, m, p, m_out);

            for (var idx = 0; idx < k + m; idx++) 
                reduce[i][j].in[idx] <== in[i][j][idx];
        }
    }

    for (var i = 0; i < l; i++) 
        for (var j = 0; j < 4; j++) 
            for (var idx = 0; idx < k; idx++) 
                out[i][j][idx] <== reduce[i][j].out[idx];
}

// Input is same as for Fp12MultiplyNoCarry
// Our answer is the prime reduction of output of Fp12MultiplyNoCarry to
//     * length 6 vectors with k registers in [0, B_a * B_b * 2^n * 12 * k^2 )
// p is length k
template Fp12MultiplyNoCarryCompress(n, k, p, m_in, m_out) {
    var l = 6;
    signal input a[l][4][k];
    signal input b[l][4][k];
    signal output out[l][4][k];

    var LOGK = log_ceil(k);
    var LOGL = log_ceil(l);
    /*assert(2*n + 1 + LOGK + LOGL <254);*/

    component nocarry = Fp12MultiplyNoCarry(n, k, 2*m_in + 4 + LOGK);
    for (var i = 0; i < l; i ++)for(var j=0; j<4; j++)for(var idx=0; idx<k; idx++){ 
         nocarry.a[i][j][idx] <== a[i][j][idx];
         nocarry.b[i][j][idx] <== b[i][j][idx];
    }

    component reduce = Fp12Compress(n, k, k-1, p, m_out);
    for (var i = 0; i < l; i++)
        for (var j = 0; j < 4; j++)
            for (var idx = 0; idx < 2 * k - 1; idx++) 
                reduce.in[i][j][idx] <== nocarry.out[i][j][idx];

    for (var i = 0; i < l; i++) 
        for (var j = 0; j < 4; j++)
            for (var idx = 0; idx < k; idx++) 
                out[i][j][idx] <== reduce.out[i][j][idx];
}

// solve for: in0 - in2 = X * p + out
// X has Ceil( overflow / n ) registers, lying in [-2^n, 2^n)
// assume in has registers in [0, 2^overflow)
template Fp12CarryModP(n, k, overflow, p) {
    var l = 6;
    var kX = (overflow + n - 1) \ n;
    signal input in[l][4][k];
    signal output X[l][2][kX];
    signal output out[l][2][k];

    assert( overflow < 252 );

    // dimension [l, 2, k]
    var Xvar0[l][2][50];
    // dimension [l, 2, kX]
    var Xvar1[l][2][50];
    for (var i = 0; i < l; i++) {
        Xvar0[i] = get_Fp12_carry_witness(n, k, kX, p, in[i][0], in[i][2]);
        Xvar1[i] = get_Fp12_carry_witness(n, k, kX, p, in[i][1], in[i][3]);

        for (var idx = 0; idx < kX; idx++) {
            X[i][0][idx] <-- Xvar0[i][0][idx];
            X[i][1][idx] <-- Xvar1[i][0][idx];
        }
        for (var idx = 0; idx < k; idx++) {
            out[i][0][idx] <-- Xvar0[i][1][idx];
            out[i][1][idx] <-- Xvar1[i][1][idx];
        }
    }
    
    component carry_check[l][2]; 
    for(var i=0; i<l; i++)for(var j=0; j<2; j++){
        carry_check[i][j] = CheckCarryModP(n, k, kX, overflow, p); 
        for(var idx=0; idx<k; idx++){
            carry_check[i][j].in[idx] <== in[i][j][idx] - in[i][j+2][idx];
            carry_check[i][j].Y[idx] <== out[i][j][idx];
        }
        for(var idx=0; idx<kX; idx++)
            carry_check[i][j].X[idx] <== X[i][j][idx];
    }

}


// version of Fp12Multiply that uses the prime reduction trick
// takes longer to compile
// assumes p has k registers with kth register nonzero
template Fp12Multiply2(n, k, p) {
    var l = 6;
    signal input a[l][2][k];
    signal input b[l][2][k];
    
    signal output out[l][2][k];

    var LOGK = log_ceil(k); 
    // registers are in [0, 2^{3n} * 12 k )
    component no_carry = Fp12MultiplyNoCarryCompress(n, k, p, n, 3*n + 2*LOGK + 4);
    for (var i = 0; i < l; i++) {
	for (var idx = 0; idx < k; idx++) {
	    for (var j = 0; j < 2; j++) {
		no_carry.a[i][j][idx] <== a[i][j][idx];
		no_carry.b[i][j][idx] <== b[i][j][idx];
	    }
	    no_carry.a[i][2][idx] <== 0;
	    no_carry.a[i][3][idx] <== 0;
	    no_carry.b[i][2][idx] <== 0;
	    no_carry.b[i][3][idx] <== 0;
	}
    }
    // This is from old Fp12MultiplyNoCarryCompress: 
    // difference of registers of no_carry lie in (-2^{3n} * 12k, 2^{3n} * 12k)
    // |diff of no_carry| is bounded by 2^{n(k - 1)} * 2^{3n} * 12k
    // p is at least 2^{n (k - 1)} and is close to 2^{nk}
    // number of registers of X_0 in: no_carry[0] - no_carry[2] = X_0 * p + X_1
    // is bounded by log(2^{n(k - 1)} * 2^{3n} * 12k / p) / log(2^n) < 3 or 4 depending on k.
    // 
	// registers of in0 and in2 have bitlength at most (kX + 1) * n, so
	// registers of in0 - in2 - X * p - out have bitlength at most
	//     n + kX * n + LOGK + 2
	// we add 4 for safety for now...

    component carry_mod;
    /*if (12 * k < p[k - 1]) {
	carry_mod = Fp12CarryModP(n, k, 3, p);
    } else {
	carry_mod = Fp12CarryModP(n, k, 4, p);
    }*/
    carry_mod = Fp12CarryModP(n, k, 3*n + 2*LOGK + 4, p);
    for (var i = 0; i < l; i++) {
	for (var idx = 0; idx < k; idx++) {
	    for (var j = 0; j < 4; j++) {
		carry_mod.in[i][j][idx] <== no_carry.out[i][j][idx];
	    }
	}
    }
    
    for (var i = 0; i < l; i++) {
	for (var idx = 0; idx < k; idx++) {
	    for (var j = 0; j < 2; j++) {
		out[i][j][idx] <== carry_mod.out[i][j][idx];
	    }
	}
    }
}

// unoptimized squaring, just takes two elements of Fp12 and multiplies them
template Fp12Square(n, k, p) {
    signal input in[6][2][k];
    signal output out[6][2][k];

    // for now just use plain multiplication, this can be optimized later
    component square = Fp12Multiply(n, k, p);
    for(var i=0; i<6; i++)for(var j=0; j<k; j++){
        square.a[i][0][j] <== in[i][0][j];
        square.a[i][1][j] <== in[i][1][j];
    
        square.b[i][0][j] <== in[i][0][j];
        square.b[i][1][j] <== in[i][1][j];
    }

    for(var i=0; i<6; i++)for(var j=0; j<k; j++){
        out[i][0][j] <== square.out[i][0][j];
        out[i][1][j] <== square.out[i][1][j];
    }
}


// not actually a relevant circuit - this only exists to test find_Fp6_inverse
template Fp6Invert(n, k, p) {
    signal input a0[2][k];
    signal input a1[2][k];
    signal input a2[2][k];
    var out[6][2][50] = find_Fp6_inverse(n, k, p, a0, a1, a2);
    signal output real_out[6][2][k];
    for (var i = 0; i < 6; i++) {
        for (var j = 0; j < 2; j ++) {
            for (var idx = 0; idx < k; idx++) {
                real_out[i][j][idx] <-- out[i][j][idx];
            }
        }
    }
}

// Call find_Fp12_inverse to compute inverse
// Then check out * in = 1, out is an array of shorts
template Fp12Invert(n, k, p){
    signal input in[6][2][k];
    signal output out[6][2][k];

    var inverse[6][2][50] = find_Fp12_inverse(n, k, p, in); // 6 x 2 x 50, only 6 x 2 x k relevant
    for (var i = 0; i < 6; i ++) {
        for (var j = 0; j < 2; j ++) {
            for (var m = 0; m < k; m ++) {
                out[i][j][m] <-- inverse[i][j][m];
            }
        }
    }

    component outRangeChecks[6][2][k];
    for(var i=0; i<6; i++) for(var j=0; j<2; j++) for(var m=0; m<k; m++) {
        outRangeChecks[i][j][m] = Num2Bits(n);
        outRangeChecks[i][j][m].in <== out[i][j][m];
    }

    component in_out = Fp12Multiply(n, k, p);
    for(var i=0; i<6; i++) for(var j=0; j<2; j++) for(var m=0; m<k; m++) {
        in_out.a[i][j][m] <== in[i][j][m];
        in_out.b[i][j][m] <== out[i][j][m];
    }

    for(var i=0; i<6; i++)for(var j=0; j<2; j++) for(var m = 0; m < k; m ++) {
        if(i == 0 && j == 0 && m == 0)
            in_out.out[i][j][m] === 1;
        else
            in_out.out[i][j][m] === 0;
    }
}

// input is an element of Fp12 
// output is input raised to the e-th power
// use the square and multiply method
// assume 0 < e < 2^254
template Fp12Exp(n, k, e, p) {
    assert( e > 0 );

    signal input in[6][2][k];
    signal output out[6][2][k];

    var temp = e;
    var BITLENGTH;
    for(var i=0; i<254; i++){
        if( temp != 0 )
            BITLENGTH = i; 
        temp = temp>>1;
    }
    BITLENGTH++;
    component pow2[BITLENGTH]; // pow2[i] = in^{2^i} 
    component mult[BITLENGTH];

    signal first[6][2][k];
    var curid = 0;

    for(var i=0; i<BITLENGTH; i++){
        // compute pow2[i] = pow2[i-1]**2
        if( i > 0 ){ // pow2[0] is never defined since there is no squaring involved
            pow2[i] = Fp12Square(n, k, p);
            for(var j=0; j<k; j++) pow2[i].p[j] <== p[j];
            if( i == 1 ){
                for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                    pow2[i].in[id][eps][j] <== in[id][eps][j];
            }else{
                for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                    pow2[i].in[id][eps][j] <== pow2[i-1].out[id][eps][j];
            }
        }
        if( ((e >> i) & 1) == 1 ){
            if(curid == 0){ // this is the least significant bit
                if( i == 0 ){
                    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                        first[id][eps][j] <== in[id][eps][j];
                }else{
                    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                        first[id][eps][j] <== pow2[i].out[id][eps][j];
                }
            }else{
                // multiply what we already have with pow2[i]
                mult[curid] = Fp12Multiply(n, k, p); 
                for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                    mult[curid].a[id][eps][j] <== pow2[i].out[id][eps][j];
                if(curid == 1){
                    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                        mult[curid].b[id][eps][j] <== first[id][eps][j];
                }else{
                    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                        mult[curid].b[id][eps][j] <== mult[curid-1].out[id][eps][j];
                }
            } 
            curid++; 
        }
    }
    curid--;
    if(curid == 0){
        for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
            out[id][eps][j] <== first[id][eps][j];
    }else{
        for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
            out[id][eps][j] <== mult[curid].out[id][eps][j];
    }
}

