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

// a[i][j], b[j][j] are short unsigned integers
// out[i][j] is a long unsigned integer
// basically multiply two-variable polynomials a, b
// use case: one variable will end up being 2**n; the other will be the field extension generator
template BigMultShortLong2DUnequal(n, ka, kb, la, lb) {
    signal input a[la][ka];
    signal input b[lb][kb];
    signal output out[la + lb -1][ka + kb -1];

    var prod_val[la + lb -1][ka + kb -1];
    for (var i = 0; i < la + lb -1; i++) {
        for (var j = 0; j < ka + kb -1; j++) {
            prod_val[i][j] = 0;
        }
    }

    for (var i1 = 0; i1 < la; i1 ++) {
        for (var i2 = 0; i2 < lb; i2 ++) {
            for (var j1 = 0; j1 < ka; j1 ++) {
                for (var j2 = 0; j2 < kb; j2 ++) {
                    var i = i1 + i2;
                    var j = j1 + j2;
                    prod_val[i][j] += a[i1][j1] * b[i2][j2];
                }
            }
        }
    }

    for (var i = 0; i < la + lb -1; i++) {
        for (var j = 0; j < ka + kb -1; j++) {
            out[i][j] <-- prod_val[i][j];
        }
    }

    var a_poly[la + lb - 1][ka + kb -1];
    var b_poly[la + lb - 1][ka + kb -1];
    var out_poly[la + lb - 1][ka + kb -1];
    for (var i = 0; i < la + lb - 1; i++) {
        for (var j = 0; j < ka + kb - 1; j++) {
            a_poly[i][j] = 0;
            b_poly[i][j] = 0;
            out_poly[i][j] = 0;
            for (var deg1 = 0; deg1 < la + lb - 1; deg1 ++) {
		if (deg1 < la) {
                    for (var deg2 = 0; deg2 < ka; deg2 ++) {
			a_poly[i][j] = a_poly[i][j] + a[deg1][deg2] * (i ** deg1) * (j ** deg2);
                    }
		}
		if (deg1 < lb) {
		    for (var deg2 = 0; deg2 < kb; deg2 ++) {
			b_poly[i][j] = b_poly[i][j] + b[deg1][deg2] * (i ** deg1) * (j ** deg2);
		    }
		}
                for (var deg2 = 0; deg2 < ka + kb -1; deg2 ++) {
                    out_poly[i][j] = out_poly[i][j] + out[deg1][deg2] * (i ** deg1) * (j ** deg2);
                }
            }
        }
    }

    for (var i = 0; i < la + lb - 1; i++) {
        for (var j = 0; j < ka + kb - 1; j++) {
            out_poly[i][j] === a_poly[i][j] * b_poly[i][j];
        }
    }
}

// p has k registers 
// inputs: 
//  a[4][k] allow overflow
//  b[4][k] 
// outputs:
//  out[4][2*k-1] such that 
//      ( (a[0] - a[1]) + (a[2] - a[3])*u ) * ( (b[0] - b[1]) + (b[2] - b[3])*u ) 
//      = (out[0] - out[1]) + (out[2] - out[3])*u  
//      we keep track of "positive" and "negatives" since circom isn't able to 
//      if each a[i][j] is in [0, B) then out[i][j] is in [0, 4*(k+1)*B^2 )
//  out[i] has 2*k-1 registers since that's output of BigMultShortLong
template Fp2multiplyNoCarry(n, k, p){
    signal input a[4][k];
    signal input b[4][k];
    signal output out[4][2*k-1];

    component ab[4][4];
    for(var i=0; i<4; i++)for(var j=0; j<4; j++){
        ab[i][j] = BigMultShortLong(n, k); // output has 2*k-1 registers
        for(var l=0; l<k; l++){
            ab[i][j].a[l] <== a[i][l];
            ab[i][j].b[l] <== b[j][l];
        }
    }
    
    for(var j=0; j<2*k-1; j++){
        out[0].in[j] <== ab[0][0].out[j] + ab[1][1].out[j] + ab[2][3].out[j] + ab[3][2].out[j];
        out[1].in[j] <== ab[0][1].out[j] + ab[1][0].out[j] + ab[2][2].out[j] + ab[3][3].out[j];
        out[2].in[j] <== ab[0][2].out[j] + ab[1][3].out[j] + ab[2][0].out[j] + ab[3][1].out[j];
        out[3].in[j] <== ab[0][3].out[j] + ab[1][2].out[j] + ab[2][1].out[j] + ab[3][0].out[j];
    }
}

// same input as above
// outputs:
//  out[4][k] such that 
//      out[i] has k registers because we use the "prime trick" to compress from 2*k-1 to k registers 
//      if each a[i][j] is in [0, B) then out[i][j] is in [0, 4*(k+1)*k*B^2 )
//          (k+1)B^2 from BigMultShortLong
//          *4 from adding 
//          *k from prime trick
template Fp2multiplyNoCarryCompress(n, k, p){
    signal input a[4][k];
    signal input b[4][k];
    signal output out[4][k];
    
    component ab = Fp2multiplyNoCarry(n, k, p);
    for(var i=0; i<4; i++)for(var j=0; j<k; j++){
        ab.a[i][j] <== a[i][j];
        ab.b[i][j] <== b[i][j]; 
    }
    
    component c[4];
    for(var i=0; i<4; i++){
        c[i] = primeTrickCompression(n, k, k-1, p);
        for(var j=0; j<2*k-1; j++)
            c[i].in[j] <== ab.out[i][j]; 
    }
 
    for(var i=0; i<4; i++)for(var j=0; j<k; j++)
        out[i][j] <== c[i].out[j];
}

// constrain in = p * X + Y 
// in[i] in (-2^overflow, 2^overflow) 
template checkBigMod(n, k, m, overflow, p){
    signal input in[k]; 
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
        carry_check.in[i] <== in[i] - pX.out[i] - Y[i]; 
    }
    for(var i=k; i<k+m-1; i++)
        carry_check.in[i] <== -pX.out[i];
}

// check if in[0] + in[0]*u is a valid point of Fp2 with in[0],in[1] both with k registers in [0,2^n) and in[i] in [0,p)
template checkValidFp2(n, k, p){
    signal input in[2][k];
    component range_checks[2][k];
    component lt[2];
    
    for(var eps=0; eps<2; eps++){
        lt[eps] = BigLessThan(n, k);
        for(var i=0; i<k; i++){
            range_checks[eps][i] = Num2Bits(n);
            range_checks[eps][i].in <== in[eps][i];
            
            lt[eps].a[i] <== in[eps][i];
            lt[eps].b[i] <== p[i];
        }
        lt[eps].out === 1;
    }    
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

    component c = Fp2multiplyNoCarryCompress(n, k, p); 
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
    var Xvar[2][2][100] = Fp2_long_div(n, k, m, c.out, p); 
    component range_check = checkValidFp2(n, k, p);
    signal X[2][m]; 
    component X_range_checks[2][m];
    
    for(var eps=0; eps<2; eps++){
        for(var i=0; i<k; i++){
            out[eps][i] <-- Xvar[eps][1][i];
            range_check.in[eps][i] <== out[eps][i];
        }
        
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
    // c3 = 0 in this case
    // constrain by Carry( c2 -' c3 -' p *' X - Y) = 0 
    // each register is an overflow representation in the range 
    //      (-2^{2n}-2^n, (k+1)*k*2^{3n+1} + 2^{2n} )
    //      which is contained in (-2^{3n+LOGK+1}, 2^{3n+LOGK+1})
    
    component mod_check[2];
    for(var eps=0; eps<2; eps++){
        mod_check[eps] = checkBigMod(n, k, m, 3*n+1+LOGK, p);
        for(var i=0; i<k; i++){
            mod_check[eps].in[i] <== c.out[2*eps][i] - c.out[2*eps+1][i];
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

// Call Fp2invert_func to compute inverse
// Then check out * in = 1, out is an array of shorts
template Fp2invert(n, k, p){
    signal input in[2][k];
    signal output out[2][k];

    var inverse[2][100] = Fp2invert_func(n, k, p, in); // 2 x 100, only 2 x k relevant
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            out[i][j] <-- inverse[i][j];
        }
    }

    //range checks
    component outRangeChecks[2][k];
    for(var i=0; i<2; i++) for(var j=0; j<k; j++){
        outRangeChecks[i][j] = Num2Bits(n);
        outRangeChecks[i][j].in <== out[i][j];
    }

    component in_out = Fp2multiply(n, k, p);
    for(var i=0; i<2; i++)for(var j=0; j<k; j++){
        in_out.a[i][j] <== in[i][j];
        in_out.b[i][j] <== out[i][j];
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

// a = sum w^i u^j a_ij for w^6=u+1, u^2=-1. similarly for b, c
// we first write a = A + B u, b = C + D u, c = E + F u and compute 
// abc = (ACE - BDE - ADF - BCF) + (ADE + BCE + ACF - BCF) u, and then simplify the representation
// assumes n, k are chosen so that cubic carries are OK
template Fp12MultiplyThree(n, k) {
    var l = 6;
    signal input a[l][2][k];
    signal input b[l][2][k];
    signal input c[l][2][k];
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
    var c0[l][k];
    var c1[l][k];
    var neg_a0[l][100];
    var neg_a1[l][100];
   for (var i = 0; i < l; i ++) { 
        for ( var j = 0; j < k; j ++) {
            a0[i][j] = a[i][0][j];
            a1[i][j] = a[i][1][j];
            b0[i][j] = b[i][0][j];
            b1[i][j] = b[i][1][j];
            c0[i][j] = c[i][0][j];
            c1[i][j] = c[i][1][j];
        }
    }
    for ( var i = 0; i < l; i ++) {
        neg_a0[i] = long_sub(n, k, p, a0[i]);
        neg_a1[i] = long_sub(n, k, p, a1[i]);
    }

    var real_init[3*l-1][100];
    var imag_init[3*l-1][100];
    var imag_init_neg[3*l-1][100];
    var real[l][2][100];
    var imag[l][2][100];
    // each product will be 3l-1 x 3k
    var a0b0c0_var[100][100] = prod3D(n, k, l, a0, b0, c0);
    var a1b1c0_neg[100][100] = prod3D(n, k, l, neg_a1, b1, c0);
    var a1b0c1_neg[100][100] = prod3D(n, k, l, neg_a1, b0, c1);
    var a0b1c1_neg[100][100] = prod3D(n, k, l, neg_a0, b1, c1);

    var a1b0c0_var[100][100] = prod3D(n, k, l, a1, b0, c0);
    var a0b1c0_var[100][100] = prod3D(n, k, l, a0, b1, c0);
    var a0b0c1_var[100][100] = prod3D(n, k, l, a0, b0, c1);
    var a1b1c1_neg[100][100] = prod3D(n, k, l, neg_a0, b1, c1);

    var a1b0c0_neg[100][100] = prod3D(n, k, l, neg_a1, b0, c0);
    var a0b1c0_neg[100][100] = prod3D(n, k, l, neg_a0, b1, c0);
    var a0b0c1_neg[100][100] = prod3D(n, k, l, neg_a0, b0, c1);
    var a1b1c1_var[100][100] = prod3D(n, k, l, a0, b1, c1);

    for (var i = 0; i < 3 * l - 1; i ++) { // compute initial rep (deg w = 10)
        real_init[i] = long_add4(n, 3 * k, a0b0c0_var[i], a1b1c0_neg[i], a1b0c1_neg[i], a0b1c1_neg[i]); // 3 * k + 1 registers each
        imag_init[i] = long_add4(n, 3 * k, a1b0c0_var[i], a0b1c0_var[i], a0b0c1_var[i], a1b1c1_neg[i]);
	imag_init_neg[i] = long_add4(n, 3 * k, a1b0c0_neg[i], a0b1c0_neg[i], a0b0c1_neg[i], a1b1c1_var[i]);
    }

    // carries using w^6 = u + 1, w^12 = 2 u
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
            real_carry[i] = long_add4(n, 3*k+1, zeros, zeros, real_init[i + l], imag_init_neg[i + l]);
            imag_carry[i] = long_add4(n, 3*k+1, zeros, zeros, real_init[i + l], imag_init[i + l]);
        } else {
            real_carry[i] = long_add4(n, 3*k+1, real_init[i + l], imag_init_neg[i + l], imag_init_neg[i + 2 * l], imag_init_neg[i + 2 * l]); // now 3*k+2 registers
            imag_carry[i] = long_add4(n, 3*k+1, imag_init[i + l], real_init[i + l], real_init[i + 2 * l], real_init[i + 2 * l]);
        }
    }    
    for (var i = 0; i < l; i ++) {
        real_final[i] = long_add_unequal(n, 3*k+2, 3*k+1, real_carry[i], real_init[i]); // now 3*k+3 registers
        imag_final[i] = long_add_unequal(n, 3*k+2, 3*k+1, imag_carry[i], imag_init[i]);
    }

    // reduction mod p
    var prod_real_temp[l][2][100];
    var prod_imag_temp[l][2][100];

    // prod_real[*][0][2 * k + 4] * p + prod_real[*][1][k] = real_final[*]
    // prod_imag[*][0][2 * k + 4] * p + prod_imag[*][1][k] = imag_final[*]
    signal prod_real[l][2][2 * k + 4];
    signal prod_imag[l][2][2 * k + 4];
    for (var i = 0; i < l; i ++) {
        prod_real_temp[i] = long_div2(n, k, 2 * k + 3, real_final[i], p); // 2 * k + 4 register quotient, k register remainder
        prod_imag_temp[i] = long_div2(n, k, 2 * k + 3, imag_final[i], p);
    }
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < 2 * k + 4; j ++) {
            prod_real[i][0][j] <-- prod_real_temp[i][0][j];
            prod_imag[i][0][j] <-- prod_imag_temp[i][0][j];
            if (j < k) {
                prod_real[i][1][j] <-- prod_real_temp[i][1][j];
                prod_imag[i][1][j] <-- prod_imag_temp[i][1][j];
            } else {
                prod_real[i][1][j] <== 0;
                prod_imag[i][1][j] <== 0;
            }
        }
    }
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < k; j ++) {
            out[i][0][j] <== prod_real[i][1][j];
            out[i][1][j] <== prod_imag[i][1][j];
        }
    }

    component out_range_checks[l][2][k];
    for(var i=0; i<l; i++){
        for (var j = 0; j < 2; j ++) {
            for (var m = 0; m < k; m ++) {
                out_range_checks[i][j][m] = Num2Bits(n);
                out_range_checks[i][j][m].in <== out[i][j][m];
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

    component div_range_checks[l][2][2 * k + 4];
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < 2 * k + 4; j ++) {
            div_range_checks[i][0][j] = Num2Bits(n);
            div_range_checks[i][1][j] = Num2Bits(n);
            div_range_checks[i][0][j].in <== prod_real[i][0][j];
            div_range_checks[i][1][j].in <== prod_imag[i][0][j];
        }
    }

    // constrain by:
    // X = a0 *' b0 *' c0 +' (p -' a1) *' b1 *' c0 +' (p -' a1) *' b0 *' c1 +' (p -' a0) *' b1 *' c1
    // Y = a1 *' b0 *' c0 +' a0 *' b1 *' c0 +' a0 *' b0 *' c1 +' (p -' a1) *' b1 *' c1
    // Carry( X_0 +' X_1 -' Y_1 -' Y_2 -' Y_2 -' p *' prod_real[0] -' prod_real[1] ) = 0
    // Carry( Y_0 +' X_1 +' Y_1 +' X_2 +' X_2 -' p *' prod_imag[0] -' prod_imag[1] ) = 0
    // where all operations are performed without carry 
    // X_0 is the coeffs of w^0, ..., w^5
    // X_1 is the coeffs of w^6, ..., w^11
    // X_2 is the coeffs of w^12, ..., w^17
    // each register is an overflow representation in the range (-kl*2^{3n+4}, kl*2^{3n + 4} )    
    component b0c0 = BigMultShortLong2D(n, k, l);
    component b0c1 = BigMultShortLong2D(n, k, l);
    component b1c0 = BigMultShortLong2D(n, k, l);
    component b1c1 = BigMultShortLong2D(n, k, l);
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < k; j ++) {
            b0c0.a[i][j] <== b[i][0][j];
            b0c0.b[i][j] <== c[i][0][j];
            b0c1.a[i][j] <== b[i][0][j];
            b0c1.b[i][j] <== c[i][1][j];
            b1c0.a[i][j] <== b[i][1][j];
            b1c0.b[i][j] <== c[i][0][j];
            b1c1.a[i][j] <== b[i][1][j];
            b1c1.b[i][j] <== c[i][1][j];
	}
    }
    
    component a0b0c0 = BigMultShortLong2DUnequal(n, k, 2 * k - 1, l, 2 * l - 1);
    component a1b0c0 = BigMultShortLong2DUnequal(n, k, 2 * k - 1, l, 2 * l - 1);
    component a0b0c1 = BigMultShortLong2DUnequal(n, k, 2 * k - 1, l, 2 * l - 1);
    component a1b0c1 = BigMultShortLong2DUnequal(n, k, 2 * k - 1, l, 2 * l - 1);
    component a0b1c0 = BigMultShortLong2DUnequal(n, k, 2 * k - 1, l, 2 * l - 1);
    component a1b1c0 = BigMultShortLong2DUnequal(n, k, 2 * k - 1, l, 2 * l - 1);
    component a0b1c1 = BigMultShortLong2DUnequal(n, k, 2 * k - 1, l, 2 * l - 1);
    component a1b1c1 = BigMultShortLong2DUnequal(n, k, 2 * k - 1, l, 2 * l - 1);

    component pb0c1 = BigMultShortLong2DUnequal(n, k, 2 * k - 1, l, 2 * l - 1);
    component pb1c0 = BigMultShortLong2DUnequal(n, k, 2 * k - 1, l, 2 * l - 1);
    component pb1c1 = BigMultShortLong2DUnequal(n, k, 2 * k - 1, l, 2 * l - 1);
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < k; j ++) {
	    a0b0c0.a[i][j] <== a[i][0][j];
	    a0b0c1.a[i][j] <== a[i][0][j];
	    a0b1c0.a[i][j] <== a[i][0][j];
	    a0b1c1.a[i][j] <== a[i][0][j];
	    a1b0c0.a[i][j] <== a[i][1][j];
	    a1b0c1.a[i][j] <== a[i][1][j];
	    a1b1c0.a[i][j] <== a[i][1][j];
	    a1b1c1.a[i][j] <== a[i][1][j];

	    pb0c1.a[i][j] <== p[j];
	    pb1c0.a[i][j] <== p[j];
	    pb1c1.a[i][j] <== p[j];
	}
    }
    for (var i = 0; i < 2 * l - 1; i++) {
	for (var j = 0; j < 2 * k - 1; j++) {
	    a0b0c0.b[i][j] <== b0c0.out[i][j];
	    a1b0c0.b[i][j] <== b0c0.out[i][j];
	    a0b0c1.b[i][j] <== b0c1.out[i][j];
	    a1b0c1.b[i][j] <== b0c1.out[i][j];
	    a0b1c0.b[i][j] <== b1c0.out[i][j];
	    a1b1c0.b[i][j] <== b1c0.out[i][j];
	    a0b1c1.b[i][j] <== b1c1.out[i][j];
	    a1b1c1.b[i][j] <== b1c1.out[i][j];

	    pb0c1.b[i][j] <== b0c1.out[i][j];
	    pb1c0.b[i][j] <== b1c0.out[i][j];
	    pb1c1.b[i][j] <== b1c1.out[i][j];
	}
    }
    
    component p_prod_real0[l];
    component p_prod_imag0[l];
    for (var i = 0; i < l; i ++) {
        p_prod_real0[i] = BigMultShortLongUnequal(n, k, 2 * k + 4);
        p_prod_imag0[i] = BigMultShortLongUnequal(n, k, 2 * k + 4);

	for (var j = 0; j < k; j++) {
	    p_prod_real0[i].a[j] <== p[j];
	    p_prod_imag0[i].a[j] <== p[j];
	}
	for (var j = 0; j < 2 * k + 4; j++) {
	    p_prod_real0[i].b[j] <== prod_real[i][0][j];
	    p_prod_imag0[i].b[j] <== prod_imag[i][0][j];
	}
    }    

    var X0[l][3 * k - 2];
    var X1[l][3 * k - 2];
    var X2[l][3 * k - 2];
    var Y0[l][3 * k - 2];
    var Y1[l][3 * k - 2];
    var Y2[l][3 * k - 2];
    for (var i = 0; i < l; i++) {
	for (var j = 0; j < 3 * k - 2; j++) {
	    X0[i][j] = a0b0c0.out[i][j] + pb1c0.out[i][j] - a1b1c0.out[i][j] + pb0c1.out[i][j] - a1b0c1.out[i][j] + pb1c1.out[i][j] - a0b1c1.out[i][j];
	    X1[i][j] = a0b0c0.out[i + l][j] + pb1c0.out[i + l][j] - a1b1c0.out[i + l][j] + pb0c1.out[i + l][j] - a1b0c1.out[i + l][j] + pb1c1.out[i + l][j] - a0b1c1.out[i + l][j];
	    Y0[i][j] = a1b0c0.out[i][j] + a0b1c0.out[i][j] + a0b0c1.out[i][j] + pb1c1.out[i][j] - a1b1c1.out[i][j];
	    Y1[i][j] = a1b0c0.out[i + l][j] + a0b1c0.out[i + l][j] + a0b0c1.out[i + l][j] + pb1c1.out[i + l][j] - a1b1c1.out[i + l][j];
	    if (i < l - 2) {
		X2[i][j] = a0b0c0.out[i + 2 * l][j] + pb1c0.out[i + 2 * l][j] - a1b1c0.out[i + 2 * l][j] + pb0c1.out[i + 2 * l][j] - a1b0c1.out[i + 2 * l][j] + pb1c1.out[i + 2 * l][j] - a0b1c1.out[i + 2 * l][j];	    
		Y2[i][j] = a1b0c0.out[i + 2 * l][j] + a0b1c0.out[i + 2 * l][j] + a0b0c1.out[i + 2 * l][j] + pb1c1.out[i + 2 * l][j] - a1b1c1.out[i + 2 * l][j];
	    } else {
		X2[i][j] = 0;
		Y2[i][j] = 0;
	    }    
	}
    }
    
    component carry_check[l][2];
    for (var i = 0; i < l; i++) {
	if (3 * k - 2 < 2 * k + 4) {
            carry_check[i][0] = CheckCarryToZero(n, 3 * n + 4 + LOGK + LOGL, 2 * k + 4);
            carry_check[i][1] = CheckCarryToZero(n, 3 * n + 4 + LOGK + LOGL, 2 * k + 4);
	} else {
            carry_check[i][0] = CheckCarryToZero(n, 3 * n + 4 + LOGK + LOGL, 3 * k - 2);
            carry_check[i][1] = CheckCarryToZero(n, 3 * n + 4 + LOGK + LOGL, 3 * k - 2);
	}
	    
        for (var j = 0; j < k; j ++) {
            carry_check[i][0].in[j] <== X0[i][j] + X1[i][j] - Y1[i][j] - Y2[i][j] - Y2[i][j] - p_prod_real0[i].out[j] - prod_real[i][1][j];
	    carry_check[i][1].in[j] <== Y0[i][j] + X1[i][j] + Y1[i][j] + X2[i][j] + X2[i][j] - p_prod_imag0[i].out[j] - prod_imag[i][1][j];
        }
	if (3 * k - 2 < 2 * k + 4) {
            for (var j = k; j < 3 * k - 2; j ++) {
		carry_check[i][0].in[j] <== X0[i][j] + X1[i][j] - Y1[i][j] - Y2[i][j] - Y2[i][j] - p_prod_real0[i].out[j] - prod_real[i][1][j];
		carry_check[i][1].in[j] <== Y0[i][j] + X1[i][j] + Y1[i][j] + X2[i][j] + X2[i][j] - p_prod_imag0[i].out[j] - prod_imag[i][1][j];
            }
            for (var j = 3 * k - 2; j < 2 * k + 4; j++) {
		carry_check[i][0].in[j] <== - prod_real[i][1][j];
		carry_check[i][1].in[j] <== - prod_imag[i][1][j];
	    }
        } else {
	    for (var j = k; j < 2 * k + 4; j ++) {
		carry_check[i][0].in[j] <== X0[i][j] + X1[i][j] - Y1[i][j] - Y2[i][j] - Y2[i][j] - p_prod_real0[i].out[j] - prod_real[i][1][j];
		carry_check[i][1].in[j] <== Y0[i][j] + X1[i][j] + Y1[i][j] + X2[i][j] + X2[i][j] - p_prod_imag0[i].out[j] - prod_imag[i][1][j];
            }
            for (var j = 2 * k + 4; j < 3 * k - 2; j++) {
		carry_check[i][0].in[j] <== X0[i][j] + X1[i][j] - Y1[i][j] - Y2[i][j] - Y2[i][j] - p_prod_real0[i].out[j];
		carry_check[i][1].in[j] <== Y0[i][j] + X1[i][j] + Y1[i][j] + X2[i][j] + X2[i][j] - p_prod_imag0[i].out[j];
	    }

	}
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
// Our answer is the prime reduction of them to
//     * length 6 vectors with k registers in [0, B_a * B_b * 2^n * 12 * k)
// p is length k
template Fp12MultiplyNoCarry(n, k, p) {
    var l = 6;
    signal input a[l][4][k];
    signal input b[l][4][k];
    signal output out[l][4][k];

    var LOGK = 3;
    var LOGL = 4;
    assert(l<15);
    assert(k<7);
    assert(2*n + 1 + LOGK + LOGL <254);

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

    component reduce[l][4];
    for (var i = 0; i < l; i++) {
	for (var j = 0; j < 4; j++) {
	    reduce[i][j] = primeTrickCompression(n, k, k - 1, p);
	}

	for (var j = 0; j < 2 * k - 1; j++) {
	    if (i < l - 1) {
		reduce[i][0].in[j] <== X[i][j] + X[l + i][j] + W[l + i][j];
		reduce[i][1].in[j] <== Y[i][j] + Y[l + i][j] + X[l + i][j];
		reduce[i][2].in[j] <== Z[i][j] + Z[l + i][j] + Y[l + i][j];
		reduce[i][3].in[j] <== W[i][j] + W[l + i][j] + Z[l + i][j];
	    } else {
		reduce[i][0].in[j] <== X[i][j];
		reduce[i][1].in[j] <== Y[i][j];
		reduce[i][2].in[j] <== Z[i][j];
		reduce[i][3].in[j] <== W[i][j];
	    }
	}
    }

    for (var i = 0; i < l; i++) {
	for (var j = 0; j < 4; j++) {
	    for (var idx = 0; idx < k; idx++) {
		out[i][j][idx] <== reduce[i][j].out[idx];
	    }
	}
    }    
}

// solve for
// a - b = out0 * p + out1
// assumes out0 has at most kX registers
function get_fp12_carry_witness(n, k, kX, p, a, b) {
    var out[2][100];

    var a_short[100] = long_to_short(n, k, a);
    var b_short[100] = long_to_short(n, k, b);

    var X[2][2][100];
    X[0] = long_div2(n, k, kX, a_short, p);
    X[1] = long_div2(n, k, kX, b_short, p);
    
    var gt = long_gt(n, k, X[1][1], X[0][1]);
    if (gt == 0){
        out[1] = long_sub(n, k, X[0][1], X[1][1]); 
        for(var i = 0; i < kX; i++) {
            out[0][i] = X[0][0][i] - X[1][0][i];
	}
    } else{
        out[1] = long_add(n, k, X[0][1], long_sub(n, k, p, X[1][1]));
        out[0][0] = X[0][0][0] - X[1][0][0] - 1;
        for(var i = 1; i < kX; i++) {
            out[0][i] = X[0][0][i] - X[1][0][i];
	}
    }

    return out;
}

// solve for: in0 - in2 = X * p + out
// assumes X has at most kX registers, lying in (-2^n, 2^n)
template Fp12CarryModP(n, k, kX, p) {
    var l = 6;
    signal input in[l][4][k];

    signal output X[l][2][kX];
    signal output out[l][2][k];

    var LOGK = 3;
    assert(k < 7);

    // dimension [l, 2, k]
    var Xvar0[100][2][100];
    // dimension [l, 2, kX]
    var Xvar1[100][2][100];
    for (var i = 0; i < l; i++) {
	Xvar0[i] = get_fp12_carry_witness(n, k, kX, p, in[i][0], in[i][2]);
	Xvar1[i] = get_fp12_carry_witness(n, k, kX, p, in[i][1], in[i][3]);

	for (var idx = 0; idx < kX; idx++) {
	    X[i][0][idx] <-- Xvar0[i][0][idx];
	    X[i][1][idx] <-- Xvar1[i][0][idx];
	}
	for (var idx = 0; idx < k; idx++) {
	    out[i][0][idx] <-- Xvar0[i][1][idx];
	    out[i][1][idx] <-- Xvar1[i][1][idx];
	}
    }

    component X_range_checks[l][2][kX];
    component out_range_checks[l][2][k];
    component lt[l][2];
    for (var i = 0; i < l; i++) {
	for (var j = 0; j < 2; j++) {
	    for (var idx = 0; idx < k; idx++) {
		out_range_checks[i][j][idx] = Num2Bits(n);
		out_range_checks[i][j][idx].in <== out[i][j][idx];
	    }
	    for (var idx = 0; idx < kX; idx++) {	    
		X_range_checks[i][j][idx] = Num2Bits(n + 1);
		X_range_checks[i][j][idx].in <== X[i][j][idx] + (1 << n);
	    }

	    lt[i][j] = BigLessThan(n, k);
	    for (var idx = 0; idx < k; idx++) {
		lt[i][j].a[idx] <== out[i][j][idx];
		lt[i][j].b[idx] <== p[idx];
	    }
	    lt[i][j].out === 1;
	}
    }

    component pX[l][2];
    component carry_check[l][2];
    for (var i = 0; i < l; i++) {
	for (var j = 0; j < 2; j++) {
	    pX[i][j] = BigMultShortLongUnequal(n, k, kX);
	    for (var idx = 0; idx < k; idx++) {
		pX[i][j].a[idx] <== p[idx];
	    }
	    for (var idx = 0; idx < kX; idx++) {
		pX[i][j].b[idx] <== X[i][j][idx];
	    }

	    // registers of in0 and in2 have bitlength at most (kX + 1) * n, so
	    // registers of in0 - in2 - X * p - out have bitlength at most
	    //     n + kX * n + LOGK + 2
	    // we add 4 for safety for now...
	    carry_check[i][j] = CheckCarryToZero(n, n + kX * n + LOGK + 4, k + kX - 1);
	}
	  
	for (var idx = 0; idx < k; idx++) {
	    carry_check[i][0].in[idx] <== in[i][0][idx] - in[i][2][idx] - pX[i][0].out[idx] - out[i][0][idx];
	    carry_check[i][1].in[idx] <== in[i][1][idx] - in[i][3][idx] - pX[i][1].out[idx] - out[i][1][idx];
	}
	for (var idx = k; idx < k + kX - 1; idx++) {
	    carry_check[i][0].in[idx] <== - pX[i][0].out[idx];
	    carry_check[i][1].in[idx] <== - pX[i][1].out[idx];
	}	
    }
}

// assumes p has k registers with kth register nonzero
template Fp12Multiply2(n, k, p) {
    var l = 6;
    signal input a[l][2][k];
    signal input b[l][2][k];
    
    signal output out[l][2][k];

    // registers are in [0, 2^{3n} * 12 k)
    component no_carry = Fp12MultiplyNoCarry(n, k, p);
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

    // difference of registers of no_carry lie in (-2^{3n} * 12k, 2^{3n} * 12k)
    // |diff of no_carry| is bounded by 2^{n(k - 1)} * 2^{3n} * 12k
    // p is at least 2^{n (k - 1)} and is close to 2^{nk}
    // number of registers of X_0 in: no_carry[0] - no_carry[2] = X_0 * p + X_1
    // is bounded by log(2^{n(k - 1)} * 2^{3n} * 12k / p) / log(2^n) < 3 or 4 depending on k.
    component carry_mod;
    if (12 * k < p[k - 1]) {
	carry_mod = Fp12CarryModP(n, k, 3, p);
    } else {
	carry_mod = Fp12CarryModP(n, k, 4, p);
    }
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

template Fp12cyclotomicDecompress(n, k, p) {
    signal input in[4][2][k];
    signal output out[6][2][k]; 

    assert(k<7);
    var LOGK = 6; // LOGK = ceil( log_2( (k+1)k ) )
    assert(3*n + 5 + LOGK<254);

    // g2 = in[0], g3 = in[1], g4 = in[2], g5 = in[3]

    // detect if g2 is 0
    component g2Zero[2];
    for(var eps=0; eps<2; eps++){
        g2Zero[eps] = BigIsZero(k);
        for(var i=0; i<k; i++)
            g2Zero[eps].in[i] <== in[0][eps][i];
    }
    signal g2isZero;
    g2isZero <== g2Zero[0].out * g2Zero[1].out; 

    // COMPUTATION OF g1 when g2 != 0:
    component g5sq = Fp2multiplyNoCarryCompress(n, k, p); // overflow (k+1)*k * 2^{3n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        g5sq.a[2*eps][i] <== in[3][eps][i];
        g5sq.a[2*eps+1][i] <== 0;
        g5sq.b[2*eps][i] <== in[3][eps][i];
        g5sq.b[2*eps+1][i] <== 0;
    }
    // c = 1+u
    signal g5sqc[4][k]; // overflow 2*(k+1)*k * 2^{3n+1}
    for(var i=0; i<k; i++){
        g5sqc[0][i] <== g5sq.out[0][i] + g5sq.out[3][i];
        g5sqc[1][i] <== g5sq.out[1][i] + g5sq.out[2][i];
        g5sqc[2][i] <== g5sq.out[0][i] + g5sq.out[2][i];
        g5sqc[3][i] <== g5sq.out[1][i] + g5sq.out[3][i];
    }
    component g4sq3 = Fp2multiplyNoCarryCompress(n, k, p); // overflow 3*(k+1)*k * 2^{3n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        g4sq3.a[2*eps][i] <== 3*in[2][eps][i];
        g4sq3.a[2*eps+1][i] <== 0;
        g4sq3.b[2*eps][i] <== in[2][eps][i];
        g4sq3.b[2*eps+1][i] <== 0;
    }
    signal g1num[4][k];  // g5^2 * (1+u) + 3 g4^2 - 2 g3
                         // overflow 5*(k+1)*k * 2^{3n+1} + 2*2^n < 2^{4n}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        g1num[2*eps][i] <== g5sqc[2*eps][i] + g4sq3.out[2*eps][i];
        g1num[2*eps+1][i] <== g5sqc[2*eps+1][i] + g4sq3.out[2*eps+1][i] + 2*in[1][eps][i];
    }
    // precompute inverse of 4g2 
    // first compute 4g2 
    var fourg2[2][100]; 
    var four[100];
    four[0] = 4;
    for(var i=1; i<100; i++) four[i] = 0;
    for(var eps=0; eps<2; eps++)
        fourg2[eps] = prod_mod(n, k, four, in[0][eps], p);
    
    // 1/(4*g2)
    var fourg2inv[2][100] = Fp2invert_func(n, k, fourg2, p);
    // precompute g1_0 = g1num / (4*g2)  
    var g1num_mod[2][100]; 
    for(var eps=0; eps<2; eps++){
        var temp1[2][100] = long_div2(n,k,4,long_to_short(n, k, g1num[2*eps]),p);
        var temp2[2][100] = long_div2(n,k,4,long_to_short(n, k, g1num[2*eps+1]),p);
        g1num_mod[eps] = long_sub_mod(n,k,temp1[1],temp2[1],p);
    }

    var g1_1var[2][100] = find_Fp2_product(n, k, g1num_mod, fourg2inv, p);
    signal g1_1[2][k]; 
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++)
        g1_1[eps][i] <-- g1_1var[eps][i]; 
    
    // constraint is g1_1 * (4g2) = g1num + p*X'' for some X'' 
    // precompute g1_1 * (4g2) = p*X + Y,  g1num = p*X' + Y',  should have Y = Y' so X'' = X-X'
    // g1_1 * (4g2) 
    component multinv1 = Fp2multiplyNoCarryCompress(n, k, p);  // overflow 4*(k+1)*k * 2^{3n+1} 
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        multinv1.a[2*eps][i] <== g1_1[eps][i]; 
        multinv1.a[2*eps+1][i] <== 0;
        multinv1.b[2*eps][i] <== 4*in[0][eps][i]; 
        multinv1.b[2*eps+1][i] <== 0;
    }
    
    var m = 4;
    component check_g11 = checkValidFp2(n, k, p);
    for(var eps=0; eps<2; eps++)
        for(var i=0; i<k; i++)
            check_g11.in[eps][i] <== g1_1[eps][i];
        
    // get multinv1 = p*X + Y 
    var XY[2][2][100] = Fp2_long_div(n, k, m, multinv1.out, p); 
    // get g1num = p*X' + Y'
    var XY1[2][2][100] = Fp2_long_div(n, k, m, g1num, p); 

    signal X[4][2][m]; 
    component X_range_checks[4][2][m];
    for(var eps=0; eps<2; eps++){    
        for(var i=0; i<m; i++){
            // X'' = X-X'
            X[0][eps][i] <-- XY[eps][0][i] - XY1[eps][0][i];
            X_range_checks[0][eps][i] = Num2Bits(n+1);
            X_range_checks[0][eps][i].in <== X[0][eps][i] + (1<<n); // X[eps][i] should be between [-2^n, 2^n)
        }
    }
    // finally constrain multinv1 - g1num = p * X'' 
    component mod_check[2];  // overflow 9*(k+1)*k * 2^{3n+1} + 2*2^n < 2^{3n+LOGK+5} 
    for(var eps=0; eps<2; eps++){
        mod_check[eps] = checkBigMod(n, k, m, 3*n+5+LOGK, p);
        for(var i=0; i<k; i++){
            mod_check[eps].in[i] <== multinv1.out[2*eps][i] - multinv1.out[2*eps+1][i] - g1num[2*eps][i] + g1num[2*eps+1][i];
            mod_check[eps].Y[i] <== 0;
        }
        for(var i=0; i<m; i++){
            mod_check[eps].X[i] <== X[0][eps][i];
        }
    }
    // END OF COMPUTATION OF g1 when g2 != 0:

    // COMPUTATION OF g1 when g2 = 0:
    // g1 = 2*g4*g5 / g3
    component twog4g5 = Fp2multiplyNoCarryCompress(n, k, p); // overflow (k+1)*k * 2^{3n+2}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        twog4g5.a[2*eps][i] <== 2*in[2][eps][i];
        twog4g5.a[2*eps+1][i] <== 0;
        twog4g5.b[2*eps][i] <== in[3][eps][i];
        twog4g5.b[2*eps+1][i] <== 0;
    }
    // precompute inverse of g3 
    var g3inv[2][100] = Fp2invert_func(n, k, in[1], p);
    // precompute g1_0 = 2g4g5 / g3
    var twog4g5_mod[2][100]; 
    for(var eps=0; eps<2; eps++){
        var temp1[2][100] = long_div2(n,k,4,long_to_short(n, k, twog4g5.out[2*eps]),p);
        var temp2[2][100] = long_div2(n,k,4,long_to_short(n, k, twog4g5.out[2*eps+1]),p);
        twog4g5_mod[eps] = long_sub_mod(n,k,temp1[1],temp2[1],p);
    }

    var g1_0var[2][100] = find_Fp2_product(n, k, twog4g5_mod, g3inv, p);
    signal g1_0[2][k]; 
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++)
        g1_0[eps][i] <-- g1_0var[eps][i]; 

    // constraint is g1_0 * g3 = twog4g5 + p*X'' for some X'' 
    // precompute g1_0 * g3 = p*X + Y,  twog4g5 = p*X' + Y',  should have Y = Y' so X'' = X-X'
    // g1_0 * g3
    component multinv0 = Fp2multiplyNoCarryCompress(n, k, p);  // overflow (k+1)*k * 2^{3n+1} 
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        multinv0.a[2*eps][i] <== g1_0[eps][i]; 
        multinv0.a[2*eps+1][i] <== 0;
        multinv0.b[2*eps][i] <== in[1][eps][i]; 
        multinv0.b[2*eps+1][i] <== 0;
    }
    
    component check_g10 = checkValidFp2(n, k, p);
    for(var eps=0; eps<2; eps++)
        for(var i=0; i<k; i++)
            check_g10.in[eps][i] <== g1_0[eps][i];
        
    // get multinv0 = p*X + Y 
    XY = Fp2_long_div(n, k, m, multinv0.out, p); 
    // get twog4g5 = p*X' + Y'
    XY1 = Fp2_long_div(n, k, m, twog4g5.out, p); 

    for(var eps=0; eps<2; eps++){    
        for(var i=0; i<m; i++){
            // X'' = X-X'
            X[1][eps][i] <-- XY[eps][0][i] - XY1[eps][0][i];
            X_range_checks[1][eps][i] = Num2Bits(n+1);
            X_range_checks[1][eps][i].in <== X[1][eps][i] + (1<<n); // X[eps][i] should be between [-2^n, 2^n)
        }
    }
    // finally constrain multinv0 - twog4g5 = p * X'' 
    component mod_check1[2];  // overflow 3*(k+1)*k * 2^{3n+1} < 2^{3n+LOGK+3} 
    for(var eps=0; eps<2; eps++){
        mod_check1[eps] = checkBigMod(n, k, m, 3*n+3+LOGK, p);
        for(var i=0; i<k; i++){
            mod_check1[eps].in[i] <== multinv0.out[2*eps][i] - multinv0.out[2*eps+1][i] - twog4g5.out[2*eps][i] + twog4g5.out[2*eps+1][i];
            mod_check1[eps].Y[i] <== 0;
        }
        for(var i=0; i<m; i++){
            mod_check1[eps].X[i] <== X[1][eps][i];
        }
    }
    // END OF COMPUTATION OF g1 when g2 = 0.

    // COMPUTATIN OF g0 when g2 != 0:
    // g0 = (2 g1^2 + g2 g5 - 3 g3 g4 )(1+u) + 1
    component twog1sq= Fp2multiplyNoCarryCompress(n, k, p); // overflow 2*(k+1)*k * 2^{3n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        twog1sq.a[2*eps][i] <== 2*g1_1[eps][i];
        twog1sq.a[2*eps+1][i] <== 0;
        twog1sq.b[2*eps][i] <== g_1_1[eps][i];
        twog1sq.b[2*eps+1][i] <== 0;
    }
    
    component g2g5 = Fp2multiplyNoCarryCompress(n, k, p); // overflow (k+1)*k * 2^{3n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        g2g5.a[2*eps][i] <== in[0][eps][i];
        g2g5.a[2*eps+1][i] <== 0;
        g2g5.b[2*eps][i] <== in[3][eps][i];
        g2g5.b[2*eps+1][i] <== 0;
    }
    
    component threeg3g4 = Fp2multiplyNoCarryCompress(n, k, p); // overflow 3*(k+1)*k * 2^{3n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        threeg3g4.a[2*eps][i] <== 3*in[1][eps][i];
        threeg3g4.a[2*eps+1][i] <== 0;
        threeg3g4.b[2*eps][i] <== in[2][eps][i];
        threeg3g4.b[2*eps+1][i] <== 0;
    }
    // 2 g1^2 + g2 g5 - 3 g3 g4 
    var temp[4][100]; 
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        temp[2*eps][i] = twog1gsq.out[2*eps][i] + g2g5.out[2*eps][i] + threeg3g4.out[2*eps+1][i];
        temp[2*eps+1][i] = twog1gsq.out[2*eps+1][i] + g2g5.out[2*eps+1][i] + threeg3g4.out[2*eps][i];
    }
    // (2 g1^2 + g2 g5 - 3 g3 g4)(1+u)
    var tempc[4][100]; // overflow 2*(k+1)*k * 2^{3n+1}
    for(var i=0; i<k; i++){
        tempc[0][i] = temp[0][i] + temp[3][i];
        tempc[1][i] = temp[1][i] + temp[2][i];
        tempc[2][i] = temp[0][i] + temp[2][i];
        tempc[3][i] = temp[1][i] + temp[3][i];
    }
    // (2 g1^2 + g2 g5 - 3 g3 g4)(1+u) + 1
    tempc[0][0]++;
    
    
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
