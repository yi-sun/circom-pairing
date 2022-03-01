pragma circom 2.0.2;

include "bigint.circom";
include "field_elements_func.circom";

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

// p has k registers 
// inputs: 
//  a[4][k] allow overflow
//  b[4][k] 
// outputs:
//  out[4][k] such that 
//      ( (a[0] - a[1]) + (a[2] - a[3])*u ) * ( (b[0] - b[1]) + (b[2] - b[3])*u ) 
//          is congruent to (out[0] - out[1]) + (out[2] - out[3])*u modulo p  
//      we keep track of "positive" and "negatives" since circom isn't able to 
//      out[4] has k registers because we use the "prime trick" to compress from 2*k-1 to k registers 
//      if each a[i][j] is in [0, B) then out[i][j] is in [0, 4*(k+1)*k*B^2 )
//          (k+1)B^2 from BigMultShortLong
//          *k from prime trick
//          *4 from adding 
template Fp2multiplyNoCarry(n, k, p){
    signal input a[4][k];
    signal input b[4][k];
    signal output out[4][k];
    
    component ab[4][4];
    for(var i=0; i<4; i++)for(var j=0; j<4; j++){
        ab[i][j] = BigMultShortLong(n, k); // output has 2*k-1 registers
        for(var l=0; l<k; l++){
            ab[i][j].a[l] <== a[i][l];
            ab[i][j].b[l] <== b[j][l];
        }
    }
    
    component c[4];
    for(var i=0; i<4; i++)
        c[i] = primeTrickCompression(n, k, k-1, p);
    for(var j=0; j<2*k-1; j++){
        c[0].in[j] <== ab[0][0].out[j] + ab[1][1].out[j] + ab[2][3].out[j] + ab[3][2].out[j];
        c[1].in[j] <== ab[0][1].out[j] + ab[1][0].out[j] + ab[2][2].out[j] + ab[3][3].out[j];
        c[2].in[j] <== ab[0][2].out[j] + ab[1][3].out[j] + ab[2][0].out[j] + ab[3][1].out[j];
        c[3].in[j] <== ab[0][3].out[j] + ab[1][2].out[j] + ab[2][1].out[j] + ab[3][0].out[j];
    }
 
    for(var i=0; i<4; i++)for(var j=0; j<k; j++)
        out[i][j] <== c[i].out[j];    
}

// constrain in[0] - in[1] = p * X + Y 
// in[i] in [0, 2^overflow) 
template BigModCheck(n, k, m, overflow, p){
    signal input in[2][k]; 
    signal input X[m];
    signal input Y[k];

    component pX;
    component carry_check;
    var maxkm;
    if(k < m) maxkm = m;
    else maxkm = k;

    pX = BigMultShortLong(n, maxkm); // p has k registers, X has m registers, so output really has k+m-1 registers 
    for(var i=0; i<maxkm; i++){
        if(i < k)
            pX.a[i] <== p[i];
        else
            pX.a[i] <== 0;
        if(i < m)
            pX.b[i] <== X[i];
        else 
            pX.b[i] <== 0;
    }
    carry_check = CheckCarryToZero(n, overflow+1, k+m-1 ); 
    for(var i=0; i<k; i++){
        carry_check.in[i] <== in[0][i] - in[1][i] - pX.out[i] - Y[i]; 
    }
    for(var i=k; i<k+m-1; i++)
        carry_check.in[i] <== -pX.out[i];
}

// multiplication specialized to Fp2
// outputs a*b in Fp2 
// (a0 + a1 u)*(b0 + b1 u) = (a0*b0 - a1*b1) + (a0*b1 + a1*b0)u 
// out[i] has k registers each in [0, 2^n)
// out[i] in [0, p)
template Fp2multiply(n, k, p){
    signal input a[2][k];
    signal input b[2][k];
    signal output out[2][k];

    assert(k<7);
    var LOGK = 6; // LOGK = ceil( log_2( (k+1)k ) )
    assert(3*n + 1 + LOGK<254);

    component c = Fp2multiplyNoCarry(n, k, p); 
    for(var i=0; i<k; i++){
        c.a[0][i] <== a[0][i];
        c.a[1][i] <== 0;
        c.a[2][i] <== a[1][i];
        c.a[3][i] <== 0;
        c.b[0][i] <== b[0][i];
        c.b[1][i] <== 0;
        c.b[2][i] <== b[1][i];
        c.b[3][i] <== 0;
    }
    // bounds below say X[eps] will only require 4 registers max 
    var m = 4;
    var Xvar[2][2][100] = Fp2prod(n, k, m, c.out, p); 
    component range_checks[2][k];
    component lt[2];
    signal X[2][m]; 
    component X_range_checks[2][m];
    
    for(var eps=0; eps<2; eps++){
        lt[eps] = BigLessThan(n, k);
        for(var i=0; i<k; i++){
            out[eps][i] <-- Xvar[eps][1][i];
            range_checks[eps][i] = Num2Bits(n);
            range_checks[eps][i].in <== out[eps][i];
            
            lt[eps].a[i] <== out[eps][i];
            lt[eps].b[i] <== p[i];
        }
        lt[eps].out === 1;
        
        for(var i=0; i<m; i++){
            X[eps][i] <-- Xvar[eps][0][i];
            X_range_checks[eps][i] = Num2Bits(n+1);
            X_range_checks[eps][i].in <== X[eps][i] + (1<<n); // X[eps][i] should be between [-2^n, 2^n)
        }
        
    }

    // out[0] constraint: X = X[0], Y = out[0] 
    // constrain by Carry( c0 -' c1 - p *' X - Y ) = 0 
    // where all operations are performed without carry 
    // each register is an overflow representation in the range 
    //      (-(k+1)*k*2^{3n}-2^{2n}-2^n, (k+1)*k*2^{3n}+2^{2n} )
    //      which is contained in (-2^{3n+LOGK+1}, 2^{3n+LOGK+1})

    // out[1] constraint: X = X[1], Y = out[1]
    // constrain by Carry( c2 -' c3 -' p *' X - Y) = 0 
    // each register is an overflow representation in the range 
    //      (-2^{2n}-2^n, (k+1)*k*2^{3n+1} + 2^{2n} )
    //      which is contained in (-2^{3n+LOGK+1}, 2^{3n+LOGK+1})
    
    component mod_check[2];
    for(var eps=0; eps<2; eps++){
        mod_check[eps] = BigModCheck(n, k, m, 3*n+1+LOGK, p);
        for(var i=0; i<k; i++){
            mod_check[eps].in[0][i] <== c.out[2*eps][i];
            mod_check[eps].in[1][i] <== c.out[2*eps+1][i];
            mod_check[eps].Y[i] <== out[eps][i];
        }
        for(var i=0; i<m; i++){
            mod_check[eps].X[i] <== X[eps][i];
        }
    }
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
template Fp2invert(n, k, p){
    signal input in[2][k];
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

    component in_out = Fp2multiply(n, k, p);
    for(var i=0; i<2; i++)for(var j=0; j<k; j++){
        in_out.a[i][j] <-- in[i][j];
        in_out.b[i][j] <-- out[i][j];
    }

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

template Fp12frobeniusMap(n, k, power){
    signal input in[6][2][k];
    signal output out[6][2][k];

    var p = get_BLS12_381_prime(n, k);
    var FP12_FROBENIUS_COEFFICIENTS[12][6][2][5] = get_Fp12_frobenius(n, k);
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
            mult_even[i][0] = BigMultModP(n, k);
            mult_even[i][1] = BigMultModP(n, k);
            for(var j=0; j<k; j++){
                mult_even[i][0].a[j] <== in[i][0][j];
                mult_even[i][1].a[j] <== in[i][1][j];

                mult_even[i][0].b[j] <== FP12_FROBENIUS_COEFFICIENTS[pow][i][0][j];
                mult_even[i][1].b[j] <== FP12_FROBENIUS_COEFFICIENTS[pow][i][0][j];
             
                mult_even[i][0].p[j] <== p[j];
                mult_even[i][1].p[j] <== p[j];
            }
            for(var j=0; j<k; j++){
                out[i][0][j] <== mult_even[i][0].out[j];
                out[i][1][j] <== mult_even[i][1].out[j];
            }
        }
    }else{
        // apply Frob to coefficients first
        for(var i=0; i<6; i++){
            in_frob[i] = Fp2frobeniusMap(n, k, pow); 
            for(var j=0; j<k; j++){
                in_frob[i].in[0][j] <== in[i][0][j];
                in_frob[i].in[1][j] <== in[i][1][j];
                in_frob[i].p[j] <== p[j];
            }
        }
        for(var j=0; j<k; j++){
            out[0][0][j] <== in_frob[0].out[0][j];
            out[0][1][j] <== in_frob[0].out[1][j];
        } 
        for(var i=1; i<6; i++){
            mult_odd[i] = Fp2multiply(n, k, p);
            for(var j=0; j<k; j++){
                for(var eps=0; eps<2; eps++){
                    mult_odd[i].a[eps][j] <== in_frob[i].out[eps][j];
                    mult_odd[i].b[eps][j] <== FP12_FROBENIUS_COEFFICIENTS[pow][i][eps][j];
                }
                mult_odd[i].p[j] <== p[j];
            }
            for(var j=0; j<k; j++){
                out[i][0][j] <== mult_odd[i].out[0][j];
                out[i][1][j] <== mult_odd[i].out[1][j];
            }
        }
    }
}

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
// we first write a = A + B u, b = C + D u and compute 
// ab = (AC + B(p-D)) + (AD+BC) u, and then simplify the representation
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

// unoptimized squaring, just takes two elements of Fp12 and multiplies them
template Fp12square(n, k) {
    signal input in[6][2][k];
    signal input p[k];
    signal output out[6][2][k];

    // for now just use plain multiplication, this can be optimized later
    component square = Fp12Multiply(n, k);
    for(var i=0; i<6; i++)for(var j=0; j<k; j++){
        square.a[i][0][j] <== in[i][0][j];
        square.a[i][1][j] <== in[i][1][j];
    
        square.b[i][0][j] <== in[i][0][j];
        square.b[i][1][j] <== in[i][1][j];
    }
    for(var i=0; i<k; i++) square.p[i] <== p[i];

    for(var i=0; i<6; i++)for(var j=0; j<k; j++){
        out[i][0][j] <== square.out[i][0][j];
        out[i][1][j] <== square.out[i][1][j];
    }
}


// assume input is an element of Fp12 in the cyclotomic subgroup GΦ₁₂
// A cyclotomic group is a subgroup of Fp^n defined by
//   GΦₙ(p) = {α ∈ Fpⁿ : α^{Φₙ(p)} = 1}

// below we implement compression and decompression for an element  GΦ₁₂ following Theorem 3.1 of https://www.ams.org/journals/mcom/2013-82-281/S0025-5718-2012-02625-1/S0025-5718-2012-02625-1.pdf
// Fp4 = Fp2(w^3) where (w^3)^2 = 1+u 
// Fp12 = Fp4(w) where w^3 = w^3 

// in = g0 + g2 w + g4 w^2 + g1 w^3 + g3 w^4 + g5 w^5 where g_i are elements of Fp2
// out = Compress(in) = [ g2, g3, g4, g5 ] 
template Fp12cyclotomicCompress(n, k) {
    signal input in[6][2][k];
    signal output out[4][2][k]; 

    for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        out[0][eps][j] <== in[1][eps][j];
        out[1][eps][j] <== in[4][eps][j];
        out[2][eps][j] <== in[2][eps][j];
        out[3][eps][j] <== in[5][eps][j];
    } 
}

// in = [g2, g3, g4, g5] where g_i are elements of Fp2
// out = Decompress(in) = g0 + g2 w + g4 w^2 + g1 w^3 + g3 w^4 + g5 w^5 where
// if g2 != 0:
//      g1 = (g5^2 * (1+u) + 3 g4^2 - 2 g3)/(4g2) 
//      g0 = (2 g1^2 + g2 * g5 - 3 g3*g4) * (1+u) + 1
// if g2 = 0:
//      g1 = (2 g4 * g5)/g3
//      g0 = (2 g1^2 - 3 g3 * g4) * (1+u)  + 1

// well this is going to be a slog... 
// not going to optimize too much since this isn't called often
template Fp12cyclotomicDecompress(n, k) {
    signal input in[4][2][k];
    signal input p[k];
    signal output out[6][2][k]; 

    // g2 = in[0], g3 = in[1], g4 = in[2], g5 = in[3]
    component g2Zero[2];
    for(var eps=0; eps<2; eps++){
        g2Zero[eps] = IsZero();
        for(var i=0; i<k; i++)
            g2Zero[eps].in[i] <== in[0][eps][i];
    }
    signal g2isZero;
    g2isZero <== g2Zero[0].out * g2Zero[1].out; 

    
}

// output is square of input 
template Fp12cyclotomicSquare(n, k) {
    signal input in[6][2][k];
    signal input p[k];
    signal output out[6][2][k];

    // for now just use plain multiplication, this can be optimized later
    component square = Fp12Multiply(n, k);
    for(var i=0; i<6; i++)for(var j=0; j<k; j++){
        square.a[i][0][j] <== in[i][0][j];
        square.a[i][1][j] <== in[i][1][j];
    
        square.b[i][0][j] <== in[i][0][j];
        square.b[i][1][j] <== in[i][1][j];
    }
    for(var i=0; i<k; i++) square.p[i] <== p[i];

    for(var i=0; i<6; i++)for(var j=0; j<k; j++){
        out[i][0][j] <== square.out[i][0][j];
        out[i][1][j] <== square.out[i][1][j];
    }
}

// assume input is an element of Fp12 in the cyclotomic subgroup GΦ12
// output is input raised to the e-th power
// use the square and multiply method
// assume 0 < e < 2^254
template Fp12cyclotomicExp(n, k, e) {
    assert( e > 0 );

    signal input in[6][2][k];
    signal input p[k];
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
            pow2[i] = Fp12cyclotomicSquare(n, k);
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
                mult[curid] = Fp12Multiply(n, k); 
                for(var j=0; j<k; j++) mult[curid].p[j] <== p[j];
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

// hard part of final exponentiation
// use equation at top of p.14 from https://eprint.iacr.org/2020/875.pdf
template hard_part(n, k){
    signal input in[6][2][k]; 
    signal input p[k];
    signal output out[6][2][k];

    var x = get_BLS12_381_parameter();  // absolute value of parameter for BLS12-381
    
    // in^{(x+1)/3} 
    component pow1 = Fp12cyclotomicExp(n, k, (x+1)\3 ); 
    for(var i=0; i<k; i++) pow1.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow1.in[id][eps][j] <== in[id][eps][j];
    
    // in^{(x+1)/3 * (x+1)}
    component pow2 = Fp12cyclotomicExp(n, k, x+1); 
    for(var i=0; i<k; i++) pow2.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow2.in[id][eps][j] <== pow1.out[id][eps][j];

    // in^{(x+1)^2/3 * -1} = pow2^-1  inverse = frob(6) in cyclotomic subgroup
    component pow3 = Fp12frobeniusMap(n, k, 6);
    for(var i=0; i<k; i++) pow3.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow3.in[id][eps][j] <== pow2.out[id][eps][j];

    // in^{(x+1)^2/3 * -x} = pow3^x 
    component pow4 = Fp12cyclotomicExp(n, k, x); 
    for(var i=0; i<k; i++) pow4.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow4.in[id][eps][j] <== pow3.out[id][eps][j];

    // in^{(x+1)^2/3 * p} = pow2^p 
    component pow5 = Fp12frobeniusMap(n, k, 1);
    for(var i=0; i<k; i++) pow5.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow5.in[id][eps][j] <== pow2.out[id][eps][j];

    // in^{(x+1)^2/3 * (-x+p)} = pow4 * pow5
    component pow6 = Fp12Multiply(n, k);
    for(var i=0; i<k; i++) pow6.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        pow6.a[id][eps][j] <== pow4.out[id][eps][j];
        pow6.b[id][eps][j] <== pow5.out[id][eps][j];
    }

    // in^{(x+1)^2/3 * (-x+p) * x}  = pow6^x
    component pow7 = Fp12cyclotomicExp(n, k, x);
    for(var i=0; i<k; i++) pow7.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow7.in[id][eps][j] <== pow6.out[id][eps][j];

    // in^{(x+1)^2/3 * (-x+p) * x^2}  = pow7^x
    component pow8 = Fp12cyclotomicExp(n, k, x);
    for(var i=0; i<k; i++) pow8.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow8.in[id][eps][j] <== pow7.out[id][eps][j];

    // in^{(x+1)^2/3 * (-x+p) * q^2} = pow6^{q^2}
    component pow9 = Fp12frobeniusMap(n, k, 2);
    for(var i=0; i<k; i++) pow9.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow9.in[id][eps][j] <== pow6.out[id][eps][j];
    
    // in^{(x+1)^2/3 * (-x+p) * -1} = pow6^{-1} = pow6^{q^6}
    component pow10 = Fp12frobeniusMap(n, k, 6);
    for(var i=0; i<k; i++) pow10.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow10.in[id][eps][j] <== pow6.out[id][eps][j];
    
    // in^{(x+1)^2/3 * (-x+p) * (x^2 + q^2)} = pow8 * pow9
    component pow11 = Fp12Multiply(n, k);
    for(var i=0; i<k; i++) pow11.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        pow11.a[id][eps][j] <== pow8.out[id][eps][j];
        pow11.b[id][eps][j] <== pow9.out[id][eps][j];
    }
    
    // in^{(x+1)^2/3 * (-x+p) * (x^2 + q^2 - 1)} = pow10 * pow11
    component pow12 = Fp12Multiply(n, k);
    for(var i=0; i<k; i++) pow12.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        pow12.a[id][eps][j] <== pow10.out[id][eps][j];
        pow12.b[id][eps][j] <== pow11.out[id][eps][j];
    }
    
    // final answer
    // in^{(x+1)^2/3 * (-x+p) * (x^2 + q^2 - 1) + 1} = pow12 * in 
    component pow13 = Fp12Multiply(n, k); 
    for(var i=0; i<k; i++) pow13.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        pow13.a[id][eps][j] <== pow12.out[id][eps][j];
        pow13.b[id][eps][j] <== in[id][eps][j];
    }
    
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        out[id][eps][j] <== pow13.out[id][eps][j];
}
