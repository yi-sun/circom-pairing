pragma circom 2.0.2;

// a[4][k] registers can overflow - let's say in [0, B) 
//  assume actual value of each a[i] < 2^{k+m} 
// p[k] registers in [0, 2^n)
// out[2][2][k] solving
//      a[0] - a[1] = p * out[0][0] + out[0][1] with out[0][1] in [0,p) 
//      a[2] - a[3] = p * out[1][0] + out[1][1] with out[1][1] in [0,p) 
// out[i][0] has m registers in range [-2^n, 2^n)
// out[i][1] has k registers in range [0, 2^n)
function Fp2_long_div(n, k, m, a, p){
    var out[2][2][40];
    // solve for X and Y such that a0*b0 + (p-a1)*b1 = p*X + Y with Y in [0,p) 
    // -a1*b1 = (p-a1)*b1 mod p
    var a_short[4][40];
    for(var i=0; i<4; i++)
        a_short[i] = long_to_short(n, k, a[i]); 

    // let me make sure everything is in <= k+m registers
    for(var i=0; i<4; i++)
        for(var j=k+m; j<40; j++)
            assert( a_short[i][j] == 0 );

    var X[4][2][40];
    for(var i=0; i<4; i++)
        X[i] = long_div2(n, k, m, a_short[i], p);    

    for(var eps=0; eps<2; eps++){
        // compute X[2*eps][1] - X[2*eps+1][1] mod p 
        var gt = long_gt(n, k, X[2*eps+1][1], X[2*eps][1]);
        if(gt == 0){
            out[eps][1] = long_sub(n, k, X[2*eps][1], X[2*eps+1][1]); 
            for(var i=0; i<m; i++)
                out[eps][0][i] = X[2*eps][0][i] - X[2*eps+1][0][i];
        }else{
            // X[2*eps][1] - X[2*eps+1][1] + p 
            out[eps][1] = long_add(n, k, X[2*eps][1], long_sub(n, k, p, X[2*eps+1][1]) );
            out[eps][0][0] = X[2*eps][0][0] - X[2*eps+1][0][0] - 1;
            for(var i=1; i<m; i++)
                out[eps][0][i] = X[2*eps][0][i] - X[2*eps+1][0][i]; 
        }
    }
    return out;
}

// helper function to precompute the product of two elements a, b in Fp2
// a[2][k], b[2][k] all registers in [0, 2^n) 
// (a0 + a1 u)*(b0 + b1 u) = (a0*b0 - a1*b1) + (a0*b1 + a1*b0)u 
// this is a direct computation - totally distinct from the combo of Fp2multiplyNoCarry and Fp2_long_div
function find_Fp2_product(n, k, a, b, p){
    var out[2][40];
    var ab[2][2][40]; 
    for(var i=0; i<2; i++)for(var j=0; j<2; j++){
        ab[i][j] = prod_mod(n,k,a[i],b[j],p);
    }
    out[0] = long_sub_mod(n,k,ab[0][0],ab[1][1],p); 
    out[1] = long_add_mod(n,k,ab[0][1],ab[1][0],p);

    return out;
}

// helper function to precompute the sum of two elements a, b in Fp2
// a[2][k], b[2][k] all registers in [0, 2^n) 
// this is a direct computation
function find_Fp2_sum(n, k, a, b, p){
    var out[2][40];
    out[0] = long_add_mod(n,k,a[0],b[0],p); 
    out[1] = long_add_mod(n,k,a[1],b[1],p);
    return out;
}

// helper function to precompute the difference of two elements a, b in Fp2
// a[2][k], b[2][k] all registers in [0, 2^n) 
// this is a direct computation
function find_Fp2_diff(n, k, a, b, p){
    var out[2][40];
    out[0] = long_sub_mod(n,k,a[0],b[0],p); 
    out[1] = long_sub_mod(n,k,a[1],b[1],p);
    return out;
}

// adapted from mod_exp in bigint_func
// a is 2 x k
// p has k registers
// e has <=50 registers 
// p is a prime
// compute a^e in Fp2 
function find_Fp2_exp(n, k, a, p, e){
    var eBits[5000]; // assume n*50 <= 5000
    for (var i = 0; i < 50; i++) {
        for (var j = 0; j < n; j++) {
            eBits[j + n * i] = (e[i] >> j) & 1;
        }
    }

    var out[2][40]; 
    for (var i = 0; i < 40; i++) {
        out[0][i] = 0;
        out[1][i] = 0;
    }
    out[0][0] = 1;
    out[1][0] = 0;

    // repeated squaring
    for (var i = 50 * n - 1; i >= 0; i--) {
        // multiply by a if bit is 0
        if (eBits[i] == 1) {
            out = find_Fp2_product(n, k, out, a, p);
        }

        // square, unless we're at the end
        if (i > 0) {
            out = find_Fp2_product(n, k, out, out, p);
        }

    }
    return out;
}

// a[4][k] elt in Fp2 
// output multiplies by 1+u
function Fp2multc(k, a){
    var out[4][40];
    for(var i=0; i<k; i++){
        out[0][i] = a[0][i] + a[3][i];
        out[1][i] = a[1][i] + a[2][i];
        out[2][i] = a[0][i] + a[2][i];
        out[3][i] = a[1][i] + a[3][i];
    }
    return out;
}


function find_Fp12_sum(n, k, a, b, p) {
    var out[6][2][40];
    for(var i=0; i<6; i++)
        out[i] = find_Fp2_sum(n, k, a[i], b[i], p);
    return out;
}

function find_Fp12_diff(n, k, a, b, p) {
    var out[6][2][40];
    for(var i=0; i<6; i++)
        out[i] = find_Fp2_diff(n, k, a[i], b[i], p);
    return out;
}

function find_Fp12_product(n, k, a, b, p) {
    var l = 6;
    var a0[l][40];
    var a1[l][40];
    var b0[l][40];
    var b1[l][40];
    var neg_b0[l][40];
    var neg_b1[l][40];
    var out[l][2][40];
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

    var real_init[2*l-1][40];
    var imag_init[2*l-1][40];
    var imag_init_neg[2*l-1][40];
    // var real[l][2][40];
    // var imag[l][2][40];
    // each product will be 2l-1 x 2k
    var a0b0_var[40][40] = prod2D(n, k, l, a0, b0);
    var a1b1_neg[40][40] = prod2D(n, k, l, a1, neg_b1);
    var a0b1_var[40][40] = prod2D(n, k, l, a0, b1);
    var a1b0_var[40][40] = prod2D(n, k, l, a1, b0);
    var a0b1_neg[40][40] = prod2D(n, k, l, a0, neg_b1);
    var a1b0_neg[40][40] = prod2D(n, k, l, a1, neg_b0);
    for (var i = 0; i < 2*l - 1; i ++) { // compute initial rep (deg w = 10)
        real_init[i] = long_add(n, 2*k, a0b0_var[i], a1b1_neg[i]); // 2*k+1 registers each
        imag_init[i] = long_add(n, 2*k, a0b1_var[i], a1b0_var[i]);
        imag_init_neg[i] = long_add(n, 2*k, a0b1_neg[i], a1b0_neg[i]);
    }
    var real_carry[l][40];
    var imag_carry[l][40];
    var real_final[l][40];
    var imag_final[l][40];
    var zeros[40]; // to balance register sizes
    for (var i = 0; i < 40; i ++) {
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
    var XYreal_temp[l][2][40];
    var XYimag_temp[l][2][40];
    for (var i = 0; i < l; i ++) {
        XYreal_temp[i] = long_div2(n, k, k+3, real_final[i], p); // k+4 register quotient, k register remainder
        XYimag_temp[i] = long_div2(n, k, k+3, imag_final[i], p);
    }
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < k; j ++) {
            out[i][0][j] = XYreal_temp[i][1][j];
            out[i][1][j] = XYimag_temp[i][1][j];
        }
    }
    return out;
}

function get_BLS12_381_parameter(){
    return 15132376222941642752;
}

// a is 2 x k, represents element of Fp^2
// out is the inverse of a in Fp^2 (2 x k array of shorts)

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
// of (a + bu). 
function Fp2invert_func(n, k, a, p) {
    var sq0[40] = prod(n, k, a[0], a[0]);
    var sq1[40] = prod(n, k, a[1], a[1]);
    var sq_sum[40] = long_add(n, 2*k, sq0, sq1);
    var sq_sum_div[2][40] = long_div2(n, k, k+1, sq_sum, p);
    // lambda = 1/(sq_sum)%p
    var lambda[40] = mod_inv(n, k, sq_sum_div[1], p);
    var out0[40] = prod(n, k, lambda, a[0]);
    var out0_div[2][40] = long_div(n, k, out0, p);
    var out[2][40];
    out[0] = out0_div[1];
    
    var out1_pre[40] = long_sub(n, k, p, a[1]);
    var out1[40] = prod(n, k, lambda, out1_pre);
    var out1_div[2][40] = long_div(n, k, out1, p);
    out[1] = out1_div[1];
    return out;
}

// a is 6 x 2 x k element of Fp^12
// compute inverse. first multiply by conjugate a + bw (a,b in Fp^6, w^6=1+u, u^2=-1)
// then reduce to inverting in Fp^6
function Fp12invert_func(n, k, p, a) {
    var A[6][2][40];
    var B[6][2][40];
    var Bw[6][2][40];
    for (var i = 0; i < 3; i ++) {
        for (var j = 0; j < 2; j ++) {
            for (var m = 0; m < k; m ++) {
                A[2*i+1][j][m] = 0;
                B[2*i+1][j][m] = 0;
                A[2*i][j][m] = a[2*i][j][m];
                B[2*i][j][m] = a[2*i+1][j][m];
                Bw[2*i][j][m] = 0;
                Bw[2*i+1][j][m] = a[2*i+1][j][m];
            }
        }
    }
    var A2[6][2][40] = find_Fp12_product(n, k, A, A, p);
    var B2[6][2][40] = find_Fp12_product(n, k, B, B, p);
    var conj[6][2][40] = find_Fp12_diff(n, k, A, Bw, p);
    var w2[6][2][40];
    for (var i = 0; i < 6; i ++) {
        for (var j = 0; j < 2; j ++) {
            for (var m = 0; m < k; m ++) {
                if (i == 2 && j == 0 && m == 0) {
                    w2[i][j][m] = 1;
                } else {
                    w2[i][j][m] = 0;
                }
            }
        }
    }
    var B2w2[6][2][40] = find_Fp12_product(n, k, B2, w2, p);
    var conjProd[6][2][40] = find_Fp12_diff(n, k, A2, B2w2, p);
    var a0[2][40];
    var a1[2][40];
    var a2[2][40];
    for (var i = 0; i < 2; i ++) {
        for (var m = 0; m < k; m ++) {
            a0[i][m] = conjProd[0][i][m];
            a1[i][m] = conjProd[2][i][m];
            a2[i][m] = conjProd[4][i][m];
        }
    }
    var conjProdInv[6][2][40] = Fp6invert_func(n, k, p, a0, a1, a2);
    var out[6][2][40] = find_Fp12_product(n, k, conj, conjProdInv, p);
    return out;
}

// compute the inverse of a0 + a1v + a2v^2 in Fp6, where 
// v^3 = 1+u, u^2 = -1, a0 a1 a2 in Fp2 (2 x k)
// returns an element in standard Fp12 representation (6 x 2 x k)
function Fp6invert_func(n, k, p, a0, a1, a2) {
    var out[6][2][40];

    var a0_squared[2][40] = find_Fp2_product(n, k, a0, a0, p);
    var a1_squared[2][40] = find_Fp2_product(n, k, a1, a1, p);
    var a2_squared[2][40] = find_Fp2_product(n, k, a2, a2, p);
    var a0a1[2][40] = find_Fp2_product(n, k, a0, a1, p);
    var a0a2[2][40] = find_Fp2_product(n, k, a0, a2, p);
    var a1a2[2][40] = find_Fp2_product(n, k, a1, a2, p);
    var a0a1a2[2][40] = find_Fp2_product(n, k, a0a1, a2, p);

    var v3[2][40]; // v^3 = 1 + u
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            if (j == 0) {
                v3[i][j] = 1;
            } else {
                v3[i][j] = 0;
            }
        }
    }

    var three_v3[2][40]; // 3v^3 = 3 + 3u
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            if (j == 0) {
                three_v3[i][j] = 3;
            } else {
                three_v3[i][j] = 0;
            }
        }
    }

    var v6[2][40]; // v^6 = 2u
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            if (i == 1 && j == 0) {
                v6[i][j] = 2;
            } else {
                v6[i][j] = 0;
            }
        }
    }

    var v0_1[2][40] = find_Fp2_product(n, k, a1a2, v3, p);
    var v0_temp[2][40] = find_Fp2_diff(n, k, a0_squared, v0_1, p); // a0^2 - a1a2v^3
    var v1_1[2][40] = find_Fp2_product(n, k, a2_squared, v3, p);
    var v1_temp[2][40] = find_Fp2_diff(n, k, v1_1, a0a1, p); // v^3a2^2 - a0a1
    var v2_temp[2][40] = find_Fp2_diff(n, k, a1_squared, a0a2, p); // a1^2 - a0a2

    var a0_cubed[2][40] = find_Fp2_product(n, k, a0, a0_squared, p);
    var a1_cubed[2][40] = find_Fp2_product(n, k, a1, a1_squared, p);
    var a2_cubed[2][40] = find_Fp2_product(n, k, a2, a2_squared, p);
    var a13v3[2][40] = find_Fp2_product(n, k, a1_cubed, v3, p);
    var a23v6[2][40] = find_Fp2_product(n, k, a2_cubed, v6, p);
    var a0a1a23v3[2][40] = find_Fp2_product(n, k, a0a1a2, three_v3, p);

    var denom_1[2][40] = find_Fp2_sum(n, k, a0_cubed, a13v3, p);
    var denom_2[2][40] = find_Fp2_diff(n, k, a23v6, a0a1a23v3, p);
    var denom[2][40] = find_Fp2_sum(n, k, denom_1, denom_2, p); // a0^3 + a1^3v^3 + a2^3v^6 - 3a0a1a2v^3

    var denom_inv[2][40] = Fp2invert_func(n, k, denom, p);

    var v0_final[2][40] = find_Fp2_product(n, k, v0_temp, denom_inv, p);
    var v1_final[2][40] = find_Fp2_product(n, k, v1_temp, denom_inv, p);
    var v2_final[2][40] = find_Fp2_product(n, k, v2_temp, denom_inv, p);

    for (var i = 1; i < 6; i = i + 2) {
        for (var j = 0; j < 2; j ++) {
            for (var m = 0; m < 40; m ++) {
                if (i > 1)
                out[i][j][m] = 0;
                else 
                out[i][j][m] = 0;//v3[j][m];
            }
        }
    }
    out[0] = v0_final;
    out[2] = v1_final;
    out[4] = v2_final;
    return out;
}

// p[] is the prime represented with k registers using X = 2^n
// p = 1 (mod 6)
// outputs Frobenius coefficients:
//  coeff[j][i] represents an element in F_p^2
//  F_p^12 = F_p^2[w] / (w^6 - (u+1)) 
//  Apply Frobenius j times to w^i: (w^i)^(p^j) = coeff[j][i] * w^i 
//  coeff[j][i] = (1+u)^{ i * (p^j - 1) / 6 } as element of Fp2
function get_Fp12_frobenius(n, k, p){
    var coeff[12][6][2][40]; // 
/*  Tried to precompute Frobenius coefficients within circom, but might be too slow
    assert(k<=10); // need array size 2*11*k

    var one[40]; 
    for(var i=1; i<40; i++) one[i] = 0;
    one[0] = 1;

    var six[40]; 
    for(var i=1; i<40; i++) six[i] = 0;
    six[0] = 6;

    // 1 + u has order dividing p^2-1
    var modulus[300] = Prod(n, k, p, p);
    modulus[0]--;

    var ppow[12][300]; // ppow[i] stores p^i (mod p^2 - 1) 
    ppow[0] = one;
    ppow[1] = p;
    for(var i=k; i<300; i++){
        ppow[1][i] = 0;
        ppow[0][i] = 0;
    }

    for(var i=2; i<12; i++){
        ppow[i] = Prod(n, (i-1)*k, ppow[i], ppow[1]); // len i*k
    }
    
    var c[2][40]; // c = 1+u
    c[0] = one;
    c[1] = one;

    for(var j=0; j<12; j++){
        if(j == 0){
            for(var i=0; i<6; i++){
                coeff[j][i][0] = one; // coeff[j][i] = 1 + 0*u
                coeff[j][i][1] = one; 
                coeff[j][i][1][0] = 0;
            }
        }else{
            // get (p^j - 1)/6
            var temp[300] = Long_Sub(n, j*k, ppow[j], ppow[0]); 
            temp[j*k] = 0;
            var temp2[2][300] = Long_Div(n, 1, j*k-1, temp, six); // len j*k
            assert( temp2[1][0] == 0 );
            var temp_mod[300];
            if( j == 1 ){
                temp_mod = temp2[0];
                for(var id=k; id<2*k; id++)
                    temp_mod[id] = 0;
            }
            else{
                var temp3[2][300] = Long_Div(n, 2*k, (j-2)*k, temp2[0], modulus);
                temp_mod = temp3[1];
            }
            // temp_mod has len <= 2*k
            
            for(var i=0; i<6; i++){
                var longi[40] = one;
                longi[0] = i;
                // i * (q^j - 1) / 6
                var e[300] = Prod(n, 2*k, temp_mod, longi); // len 2*k+1
                 
                coeff[j][i] = find_Fp2_exp(n, k, c, p, e); 
            }
        }   
    }
*/
    assert( (n == 77 && k == 5) || (n == 43 && k == 9) );
// PRECOMPUTED FROM PYTHON
if( n == 43 && k == 9 ){
coeff[0][0][0][0] = 1;
coeff[0][0][0][1] = 0;
coeff[0][0][0][2] = 0;
coeff[0][0][0][3] = 0;
coeff[0][0][0][4] = 0;
coeff[0][0][0][5] = 0;
coeff[0][0][0][6] = 0;
coeff[0][0][0][7] = 0;
coeff[0][0][0][8] = 0;

coeff[0][0][1][0] = 0;
coeff[0][0][1][1] = 0;
coeff[0][0][1][2] = 0;
coeff[0][0][1][3] = 0;
coeff[0][0][1][4] = 0;
coeff[0][0][1][5] = 0;
coeff[0][0][1][6] = 0;
coeff[0][0][1][7] = 0;
coeff[0][0][1][8] = 0;

coeff[0][1][0][0] = 1;
coeff[0][1][0][1] = 0;
coeff[0][1][0][2] = 0;
coeff[0][1][0][3] = 0;
coeff[0][1][0][4] = 0;
coeff[0][1][0][5] = 0;
coeff[0][1][0][6] = 0;
coeff[0][1][0][7] = 0;
coeff[0][1][0][8] = 0;

coeff[0][1][1][0] = 0;
coeff[0][1][1][1] = 0;
coeff[0][1][1][2] = 0;
coeff[0][1][1][3] = 0;
coeff[0][1][1][4] = 0;
coeff[0][1][1][5] = 0;
coeff[0][1][1][6] = 0;
coeff[0][1][1][7] = 0;
coeff[0][1][1][8] = 0;

coeff[0][2][0][0] = 1;
coeff[0][2][0][1] = 0;
coeff[0][2][0][2] = 0;
coeff[0][2][0][3] = 0;
coeff[0][2][0][4] = 0;
coeff[0][2][0][5] = 0;
coeff[0][2][0][6] = 0;
coeff[0][2][0][7] = 0;
coeff[0][2][0][8] = 0;

coeff[0][2][1][0] = 0;
coeff[0][2][1][1] = 0;
coeff[0][2][1][2] = 0;
coeff[0][2][1][3] = 0;
coeff[0][2][1][4] = 0;
coeff[0][2][1][5] = 0;
coeff[0][2][1][6] = 0;
coeff[0][2][1][7] = 0;
coeff[0][2][1][8] = 0;

coeff[0][3][0][0] = 1;
coeff[0][3][0][1] = 0;
coeff[0][3][0][2] = 0;
coeff[0][3][0][3] = 0;
coeff[0][3][0][4] = 0;
coeff[0][3][0][5] = 0;
coeff[0][3][0][6] = 0;
coeff[0][3][0][7] = 0;
coeff[0][3][0][8] = 0;

coeff[0][3][1][0] = 0;
coeff[0][3][1][1] = 0;
coeff[0][3][1][2] = 0;
coeff[0][3][1][3] = 0;
coeff[0][3][1][4] = 0;
coeff[0][3][1][5] = 0;
coeff[0][3][1][6] = 0;
coeff[0][3][1][7] = 0;
coeff[0][3][1][8] = 0;

coeff[0][4][0][0] = 1;
coeff[0][4][0][1] = 0;
coeff[0][4][0][2] = 0;
coeff[0][4][0][3] = 0;
coeff[0][4][0][4] = 0;
coeff[0][4][0][5] = 0;
coeff[0][4][0][6] = 0;
coeff[0][4][0][7] = 0;
coeff[0][4][0][8] = 0;

coeff[0][4][1][0] = 0;
coeff[0][4][1][1] = 0;
coeff[0][4][1][2] = 0;
coeff[0][4][1][3] = 0;
coeff[0][4][1][4] = 0;
coeff[0][4][1][5] = 0;
coeff[0][4][1][6] = 0;
coeff[0][4][1][7] = 0;
coeff[0][4][1][8] = 0;

coeff[0][5][0][0] = 1;
coeff[0][5][0][1] = 0;
coeff[0][5][0][2] = 0;
coeff[0][5][0][3] = 0;
coeff[0][5][0][4] = 0;
coeff[0][5][0][5] = 0;
coeff[0][5][0][6] = 0;
coeff[0][5][0][7] = 0;
coeff[0][5][0][8] = 0;

coeff[0][5][1][0] = 0;
coeff[0][5][1][1] = 0;
coeff[0][5][1][2] = 0;
coeff[0][5][1][3] = 0;
coeff[0][5][1][4] = 0;
coeff[0][5][1][5] = 0;
coeff[0][5][1][6] = 0;
coeff[0][5][1][7] = 0;
coeff[0][5][1][8] = 0;

coeff[1][0][0][0] = 1;
coeff[1][0][0][1] = 0;
coeff[1][0][0][2] = 0;
coeff[1][0][0][3] = 0;
coeff[1][0][0][4] = 0;
coeff[1][0][0][5] = 0;
coeff[1][0][0][6] = 0;
coeff[1][0][0][7] = 0;
coeff[1][0][0][8] = 0;

coeff[1][0][1][0] = 0;
coeff[1][0][1][1] = 0;
coeff[1][0][1][2] = 0;
coeff[1][0][1][3] = 0;
coeff[1][0][1][4] = 0;
coeff[1][0][1][5] = 0;
coeff[1][0][1][6] = 0;
coeff[1][0][1][7] = 0;
coeff[1][0][1][8] = 0;

coeff[1][1][0][0] = 6517917179832;
coeff[1][1][0][1] = 5429504614638;
coeff[1][1][0][2] = 4234746852751;
coeff[1][1][0][3] = 2112089840866;
coeff[1][1][0][4] = 4217472791108;
coeff[1][1][0][5] = 6733099956857;
coeff[1][1][0][6] = 8439745676103;
coeff[1][1][0][7] = 6426130780557;
coeff[1][1][0][8] = 107455168258;

coeff[1][1][1][0] = 2278175820531;
coeff[1][1][1][1] = 6115366903537;
coeff[1][1][1][2] = 5088284968245;
coeff[1][1][1][3] = 8129181589039;
coeff[1][1][1][4] = 4927250153672;
coeff[1][1][1][5] = 726966144877;
coeff[1][1][1][6] = 2476195964910;
coeff[1][1][1][7] = 6764599985999;
coeff[1][1][1][8] = 4231932726;

coeff[1][2][0][0] = 0;
coeff[1][2][0][1] = 0;
coeff[1][2][0][2] = 0;
coeff[1][2][0][3] = 0;
coeff[1][2][0][4] = 0;
coeff[1][2][0][5] = 0;
coeff[1][2][0][6] = 0;
coeff[1][2][0][7] = 0;
coeff[1][2][0][8] = 0;

coeff[1][2][1][0] = 43692;
coeff[1][2][1][1] = 1374384390048;
coeff[1][2][1][2] = 5507500715325;
coeff[1][2][1][3] = 5164830044109;
coeff[1][2][1][4] = 8081740699602;
coeff[1][2][1][5] = 3659765840658;
coeff[1][2][1][6] = 144300128161;
coeff[1][2][1][7] = 4394637549586;
coeff[1][2][1][8] = 111687100985;

coeff[1][3][0][0] = 5480074431497;
coeff[1][3][0][1] = 6050511651344;
coeff[1][3][0][2] = 8493807746507;
coeff[1][3][0][3] = 7745978310882;
coeff[1][3][0][4] = 5721800802166;
coeff[1][3][0][5] = 5018446681989;
coeff[1][3][0][6] = 942123802623;
coeff[1][3][0][7] = 8770329198991;
coeff[1][3][0][8] = 28706735159;

coeff[1][3][1][0] = 5480074431497;
coeff[1][3][1][1] = 6050511651344;
coeff[1][3][1][2] = 8493807746507;
coeff[1][3][1][3] = 7745978310882;
coeff[1][3][1][4] = 5721800802166;
coeff[1][3][1][5] = 5018446681989;
coeff[1][3][1][6] = 942123802623;
coeff[1][3][1][7] = 8770329198991;
coeff[1][3][1][8] = 28706735159;

coeff[1][4][0][0] = 43693;
coeff[1][4][0][1] = 1374384390048;
coeff[1][4][0][2] = 5507500715325;
coeff[1][4][0][3] = 5164830044109;
coeff[1][4][0][4] = 8081740699602;
coeff[1][4][0][5] = 3659765840658;
coeff[1][4][0][6] = 144300128161;
coeff[1][4][0][7] = 4394637549586;
coeff[1][4][0][8] = 111687100985;

coeff[1][4][1][0] = 0;
coeff[1][4][1][1] = 0;
coeff[1][4][1][2] = 0;
coeff[1][4][1][3] = 0;
coeff[1][4][1][4] = 0;
coeff[1][4][1][5] = 0;
coeff[1][4][1][6] = 0;
coeff[1][4][1][7] = 0;
coeff[1][4][1][8] = 0;

coeff[1][5][0][0] = 3201898610966;
coeff[1][5][0][1] = 8731237770015;
coeff[1][5][0][2] = 3405522778261;
coeff[1][5][0][3] = 8412889744051;
coeff[1][5][0][4] = 794550648493;
coeff[1][5][0][5] = 4291480537112;
coeff[1][5][0][6] = 7262020859921;
coeff[1][5][0][7] = 2005729212991;
coeff[1][5][0][8] = 24474802433;

coeff[1][5][1][0] = 5594194389397;
coeff[1][5][1][1] = 2813633748160;
coeff[1][5][1][2] = 5917509042735;
coeff[1][5][1][3] = 1828381685854;
coeff[1][5][1][4] = 8350172296287;
coeff[1][5][1][5] = 3168585564622;
coeff[1][5][1][6] = 3653920781092;
coeff[1][5][1][7] = 2388908531357;
coeff[1][5][1][8] = 87212298552;

coeff[2][0][0][0] = 1;
coeff[2][0][0][1] = 0;
coeff[2][0][0][2] = 0;
coeff[2][0][0][3] = 0;
coeff[2][0][0][4] = 0;
coeff[2][0][0][5] = 0;
coeff[2][0][0][6] = 0;
coeff[2][0][0][7] = 0;
coeff[2][0][0][8] = 0;

coeff[2][0][1][0] = 0;
coeff[2][0][1][1] = 0;
coeff[2][0][1][2] = 0;
coeff[2][0][1][3] = 0;
coeff[2][0][1][4] = 0;
coeff[2][0][1][5] = 0;
coeff[2][0][1][6] = 0;
coeff[2][0][1][7] = 0;
coeff[2][0][1][8] = 0;

coeff[2][1][0][0] = 8796092956671;
coeff[2][1][0][1] = 1374394105919;
coeff[2][1][0][2] = 3815531105672;
coeff[2][1][0][3] = 5076441385796;
coeff[2][1][0][4] = 1062982245178;
coeff[2][1][0][5] = 3800300261076;
coeff[2][1][0][6] = 1975548490644;
coeff[2][1][0][7] = 194763;
coeff[2][1][0][8] = 0;

coeff[2][1][1][0] = 0;
coeff[2][1][1][1] = 0;
coeff[2][1][1][2] = 0;
coeff[2][1][1][3] = 0;
coeff[2][1][1][4] = 0;
coeff[2][1][1][5] = 0;
coeff[2][1][1][6] = 0;
coeff[2][1][1][7] = 0;
coeff[2][1][1][8] = 0;

coeff[2][2][0][0] = 8796092956670;
coeff[2][2][0][1] = 1374394105919;
coeff[2][2][0][2] = 3815531105672;
coeff[2][2][0][3] = 5076441385796;
coeff[2][2][0][4] = 1062982245178;
coeff[2][2][0][5] = 3800300261076;
coeff[2][2][0][6] = 1975548490644;
coeff[2][2][0][7] = 194763;
coeff[2][2][0][8] = 0;

coeff[2][2][1][0] = 0;
coeff[2][2][1][1] = 0;
coeff[2][2][1][2] = 0;
coeff[2][2][1][3] = 0;
coeff[2][2][1][4] = 0;
coeff[2][2][1][5] = 0;
coeff[2][2][1][6] = 0;
coeff[2][2][1][7] = 0;
coeff[2][2][1][8] = 0;

coeff[2][3][0][0] = 8796093000362;
coeff[2][3][0][1] = 2748778495967;
coeff[2][3][0][2] = 526938798789;
coeff[2][3][0][3] = 1445178407698;
coeff[2][3][0][4] = 348629922573;
coeff[2][3][0][5] = 7460066101735;
coeff[2][3][0][6] = 2119848618805;
coeff[2][3][0][7] = 4394637744349;
coeff[2][3][0][8] = 111687100985;

coeff[2][3][1][0] = 0;
coeff[2][3][1][1] = 0;
coeff[2][3][1][2] = 0;
coeff[2][3][1][3] = 0;
coeff[2][3][1][4] = 0;
coeff[2][3][1][5] = 0;
coeff[2][3][1][6] = 0;
coeff[2][3][1][7] = 0;
coeff[2][3][1][8] = 0;

coeff[2][4][0][0] = 43692;
coeff[2][4][0][1] = 1374384390048;
coeff[2][4][0][2] = 5507500715325;
coeff[2][4][0][3] = 5164830044109;
coeff[2][4][0][4] = 8081740699602;
coeff[2][4][0][5] = 3659765840658;
coeff[2][4][0][6] = 144300128161;
coeff[2][4][0][7] = 4394637549586;
coeff[2][4][0][8] = 111687100985;

coeff[2][4][1][0] = 0;
coeff[2][4][1][1] = 0;
coeff[2][4][1][2] = 0;
coeff[2][4][1][3] = 0;
coeff[2][4][1][4] = 0;
coeff[2][4][1][5] = 0;
coeff[2][4][1][6] = 0;
coeff[2][4][1][7] = 0;
coeff[2][4][1][8] = 0;

coeff[2][5][0][0] = 43693;
coeff[2][5][0][1] = 1374384390048;
coeff[2][5][0][2] = 5507500715325;
coeff[2][5][0][3] = 5164830044109;
coeff[2][5][0][4] = 8081740699602;
coeff[2][5][0][5] = 3659765840658;
coeff[2][5][0][6] = 144300128161;
coeff[2][5][0][7] = 4394637549586;
coeff[2][5][0][8] = 111687100985;

coeff[2][5][1][0] = 0;
coeff[2][5][1][1] = 0;
coeff[2][5][1][2] = 0;
coeff[2][5][1][3] = 0;
coeff[2][5][1][4] = 0;
coeff[2][5][1][5] = 0;
coeff[2][5][1][6] = 0;
coeff[2][5][1][7] = 0;
coeff[2][5][1][8] = 0;

coeff[3][0][0][0] = 1;
coeff[3][0][0][1] = 0;
coeff[3][0][0][2] = 0;
coeff[3][0][0][3] = 0;
coeff[3][0][0][4] = 0;
coeff[3][0][0][5] = 0;
coeff[3][0][0][6] = 0;
coeff[3][0][0][7] = 0;
coeff[3][0][0][8] = 0;

coeff[3][0][1][0] = 0;
coeff[3][0][1][1] = 0;
coeff[3][0][1][2] = 0;
coeff[3][0][1][3] = 0;
coeff[3][0][1][4] = 0;
coeff[3][0][1][5] = 0;
coeff[3][0][1][6] = 0;
coeff[3][0][1][7] = 0;
coeff[3][0][1][8] = 0;

coeff[3][1][0][0] = 3316018568866;
coeff[3][1][0][1] = 5494359866831;
coeff[3][1][0][2] = 829224074489;
coeff[3][1][0][3] = 2495293119023;
coeff[3][1][0][4] = 3422922142614;
coeff[3][1][0][5] = 2441619419745;
coeff[3][1][0][6] = 1177724816182;
coeff[3][1][0][7] = 4420401567566;
coeff[3][1][0][8] = 82980365825;

coeff[3][1][1][0] = 5480074431497;
coeff[3][1][1][1] = 6050511651344;
coeff[3][1][1][2] = 8493807746507;
coeff[3][1][1][3] = 7745978310882;
coeff[3][1][1][4] = 5721800802166;
coeff[3][1][1][5] = 5018446681989;
coeff[3][1][1][6] = 942123802623;
coeff[3][1][1][7] = 8770329198991;
coeff[3][1][1][8] = 28706735159;

coeff[3][2][0][0] = 0;
coeff[3][2][0][1] = 0;
coeff[3][2][0][2] = 0;
coeff[3][2][0][3] = 0;
coeff[3][2][0][4] = 0;
coeff[3][2][0][5] = 0;
coeff[3][2][0][6] = 0;
coeff[3][2][0][7] = 0;
coeff[3][2][0][8] = 0;

coeff[3][2][1][0] = 1;
coeff[3][2][1][1] = 0;
coeff[3][2][1][2] = 0;
coeff[3][2][1][3] = 0;
coeff[3][2][1][4] = 0;
coeff[3][2][1][5] = 0;
coeff[3][2][1][6] = 0;
coeff[3][2][1][7] = 0;
coeff[3][2][1][8] = 0;

coeff[3][3][0][0] = 3316018568866;
coeff[3][3][0][1] = 5494359866831;
coeff[3][3][0][2] = 829224074489;
coeff[3][3][0][3] = 2495293119023;
coeff[3][3][0][4] = 3422922142614;
coeff[3][3][0][5] = 2441619419745;
coeff[3][3][0][6] = 1177724816182;
coeff[3][3][0][7] = 4420401567566;
coeff[3][3][0][8] = 82980365825;

coeff[3][3][1][0] = 3316018568866;
coeff[3][3][1][1] = 5494359866831;
coeff[3][3][1][2] = 829224074489;
coeff[3][3][1][3] = 2495293119023;
coeff[3][3][1][4] = 3422922142614;
coeff[3][3][1][5] = 2441619419745;
coeff[3][3][1][6] = 1177724816182;
coeff[3][3][1][7] = 4420401567566;
coeff[3][3][1][8] = 82980365825;

coeff[3][4][0][0] = 8796093000362;
coeff[3][4][0][1] = 2748778495967;
coeff[3][4][0][2] = 526938798789;
coeff[3][4][0][3] = 1445178407698;
coeff[3][4][0][4] = 348629922573;
coeff[3][4][0][5] = 7460066101735;
coeff[3][4][0][6] = 2119848618805;
coeff[3][4][0][7] = 4394637744349;
coeff[3][4][0][8] = 111687100985;

coeff[3][4][1][0] = 0;
coeff[3][4][1][1] = 0;
coeff[3][4][1][2] = 0;
coeff[3][4][1][3] = 0;
coeff[3][4][1][4] = 0;
coeff[3][4][1][5] = 0;
coeff[3][4][1][6] = 0;
coeff[3][4][1][7] = 0;
coeff[3][4][1][8] = 0;

coeff[3][5][0][0] = 5480074431497;
coeff[3][5][0][1] = 6050511651344;
coeff[3][5][0][2] = 8493807746507;
coeff[3][5][0][3] = 7745978310882;
coeff[3][5][0][4] = 5721800802166;
coeff[3][5][0][5] = 5018446681989;
coeff[3][5][0][6] = 942123802623;
coeff[3][5][0][7] = 8770329198991;
coeff[3][5][0][8] = 28706735159;

coeff[3][5][1][0] = 3316018568866;
coeff[3][5][1][1] = 5494359866831;
coeff[3][5][1][2] = 829224074489;
coeff[3][5][1][3] = 2495293119023;
coeff[3][5][1][4] = 3422922142614;
coeff[3][5][1][5] = 2441619419745;
coeff[3][5][1][6] = 1177724816182;
coeff[3][5][1][7] = 4420401567566;
coeff[3][5][1][8] = 82980365825;

coeff[4][0][0][0] = 1;
coeff[4][0][0][1] = 0;
coeff[4][0][0][2] = 0;
coeff[4][0][0][3] = 0;
coeff[4][0][0][4] = 0;
coeff[4][0][0][5] = 0;
coeff[4][0][0][6] = 0;
coeff[4][0][0][7] = 0;
coeff[4][0][0][8] = 0;

coeff[4][0][1][0] = 0;
coeff[4][0][1][1] = 0;
coeff[4][0][1][2] = 0;
coeff[4][0][1][3] = 0;
coeff[4][0][1][4] = 0;
coeff[4][0][1][5] = 0;
coeff[4][0][1][6] = 0;
coeff[4][0][1][7] = 0;
coeff[4][0][1][8] = 0;

coeff[4][1][0][0] = 8796092956670;
coeff[4][1][0][1] = 1374394105919;
coeff[4][1][0][2] = 3815531105672;
coeff[4][1][0][3] = 5076441385796;
coeff[4][1][0][4] = 1062982245178;
coeff[4][1][0][5] = 3800300261076;
coeff[4][1][0][6] = 1975548490644;
coeff[4][1][0][7] = 194763;
coeff[4][1][0][8] = 0;

coeff[4][1][1][0] = 0;
coeff[4][1][1][1] = 0;
coeff[4][1][1][2] = 0;
coeff[4][1][1][3] = 0;
coeff[4][1][1][4] = 0;
coeff[4][1][1][5] = 0;
coeff[4][1][1][6] = 0;
coeff[4][1][1][7] = 0;
coeff[4][1][1][8] = 0;

coeff[4][2][0][0] = 43692;
coeff[4][2][0][1] = 1374384390048;
coeff[4][2][0][2] = 5507500715325;
coeff[4][2][0][3] = 5164830044109;
coeff[4][2][0][4] = 8081740699602;
coeff[4][2][0][5] = 3659765840658;
coeff[4][2][0][6] = 144300128161;
coeff[4][2][0][7] = 4394637549586;
coeff[4][2][0][8] = 111687100985;

coeff[4][2][1][0] = 0;
coeff[4][2][1][1] = 0;
coeff[4][2][1][2] = 0;
coeff[4][2][1][3] = 0;
coeff[4][2][1][4] = 0;
coeff[4][2][1][5] = 0;
coeff[4][2][1][6] = 0;
coeff[4][2][1][7] = 0;
coeff[4][2][1][8] = 0;

coeff[4][3][0][0] = 1;
coeff[4][3][0][1] = 0;
coeff[4][3][0][2] = 0;
coeff[4][3][0][3] = 0;
coeff[4][3][0][4] = 0;
coeff[4][3][0][5] = 0;
coeff[4][3][0][6] = 0;
coeff[4][3][0][7] = 0;
coeff[4][3][0][8] = 0;

coeff[4][3][1][0] = 0;
coeff[4][3][1][1] = 0;
coeff[4][3][1][2] = 0;
coeff[4][3][1][3] = 0;
coeff[4][3][1][4] = 0;
coeff[4][3][1][5] = 0;
coeff[4][3][1][6] = 0;
coeff[4][3][1][7] = 0;
coeff[4][3][1][8] = 0;

coeff[4][4][0][0] = 8796092956670;
coeff[4][4][0][1] = 1374394105919;
coeff[4][4][0][2] = 3815531105672;
coeff[4][4][0][3] = 5076441385796;
coeff[4][4][0][4] = 1062982245178;
coeff[4][4][0][5] = 3800300261076;
coeff[4][4][0][6] = 1975548490644;
coeff[4][4][0][7] = 194763;
coeff[4][4][0][8] = 0;

coeff[4][4][1][0] = 0;
coeff[4][4][1][1] = 0;
coeff[4][4][1][2] = 0;
coeff[4][4][1][3] = 0;
coeff[4][4][1][4] = 0;
coeff[4][4][1][5] = 0;
coeff[4][4][1][6] = 0;
coeff[4][4][1][7] = 0;
coeff[4][4][1][8] = 0;

coeff[4][5][0][0] = 43692;
coeff[4][5][0][1] = 1374384390048;
coeff[4][5][0][2] = 5507500715325;
coeff[4][5][0][3] = 5164830044109;
coeff[4][5][0][4] = 8081740699602;
coeff[4][5][0][5] = 3659765840658;
coeff[4][5][0][6] = 144300128161;
coeff[4][5][0][7] = 4394637549586;
coeff[4][5][0][8] = 111687100985;

coeff[4][5][1][0] = 0;
coeff[4][5][1][1] = 0;
coeff[4][5][1][2] = 0;
coeff[4][5][1][3] = 0;
coeff[4][5][1][4] = 0;
coeff[4][5][1][5] = 0;
coeff[4][5][1][6] = 0;
coeff[4][5][1][7] = 0;
coeff[4][5][1][8] = 0;

coeff[5][0][0][0] = 1;
coeff[5][0][0][1] = 0;
coeff[5][0][0][2] = 0;
coeff[5][0][0][3] = 0;
coeff[5][0][0][4] = 0;
coeff[5][0][0][5] = 0;
coeff[5][0][0][6] = 0;
coeff[5][0][0][7] = 0;
coeff[5][0][0][8] = 0;

coeff[5][0][1][0] = 0;
coeff[5][0][1][1] = 0;
coeff[5][0][1][2] = 0;
coeff[5][0][1][3] = 0;
coeff[5][0][1][4] = 0;
coeff[5][0][1][5] = 0;
coeff[5][0][1][6] = 0;
coeff[5][0][1][7] = 0;
coeff[5][0][1][8] = 0;

coeff[5][1][0][0] = 5594194389397;
coeff[5][1][0][1] = 2813633748160;
coeff[5][1][0][2] = 5917509042735;
coeff[5][1][0][3] = 1828381685854;
coeff[5][1][0][4] = 8350172296287;
coeff[5][1][0][5] = 3168585564622;
coeff[5][1][0][6] = 3653920781092;
coeff[5][1][0][7] = 2388908531357;
coeff[5][1][0][8] = 87212298552;

coeff[5][1][1][0] = 3201898610966;
coeff[5][1][1][1] = 8731237770015;
coeff[5][1][1][2] = 3405522778261;
coeff[5][1][1][3] = 8412889744051;
coeff[5][1][1][4] = 794550648493;
coeff[5][1][1][5] = 4291480537112;
coeff[5][1][1][6] = 7262020859921;
coeff[5][1][1][7] = 2005729212991;
coeff[5][1][1][8] = 24474802433;

coeff[5][2][0][0] = 0;
coeff[5][2][0][1] = 0;
coeff[5][2][0][2] = 0;
coeff[5][2][0][3] = 0;
coeff[5][2][0][4] = 0;
coeff[5][2][0][5] = 0;
coeff[5][2][0][6] = 0;
coeff[5][2][0][7] = 0;
coeff[5][2][0][8] = 0;

coeff[5][2][1][0] = 8796092956670;
coeff[5][2][1][1] = 1374394105919;
coeff[5][2][1][2] = 3815531105672;
coeff[5][2][1][3] = 5076441385796;
coeff[5][2][1][4] = 1062982245178;
coeff[5][2][1][5] = 3800300261076;
coeff[5][2][1][6] = 1975548490644;
coeff[5][2][1][7] = 194763;
coeff[5][2][1][8] = 0;

coeff[5][3][0][0] = 5480074431497;
coeff[5][3][0][1] = 6050511651344;
coeff[5][3][0][2] = 8493807746507;
coeff[5][3][0][3] = 7745978310882;
coeff[5][3][0][4] = 5721800802166;
coeff[5][3][0][5] = 5018446681989;
coeff[5][3][0][6] = 942123802623;
coeff[5][3][0][7] = 8770329198991;
coeff[5][3][0][8] = 28706735159;

coeff[5][3][1][0] = 5480074431497;
coeff[5][3][1][1] = 6050511651344;
coeff[5][3][1][2] = 8493807746507;
coeff[5][3][1][3] = 7745978310882;
coeff[5][3][1][4] = 5721800802166;
coeff[5][3][1][5] = 5018446681989;
coeff[5][3][1][6] = 942123802623;
coeff[5][3][1][7] = 8770329198991;
coeff[5][3][1][8] = 28706735159;

coeff[5][4][0][0] = 8796092956671;
coeff[5][4][0][1] = 1374394105919;
coeff[5][4][0][2] = 3815531105672;
coeff[5][4][0][3] = 5076441385796;
coeff[5][4][0][4] = 1062982245178;
coeff[5][4][0][5] = 3800300261076;
coeff[5][4][0][6] = 1975548490644;
coeff[5][4][0][7] = 194763;
coeff[5][4][0][8] = 0;

coeff[5][4][1][0] = 0;
coeff[5][4][1][1] = 0;
coeff[5][4][1][2] = 0;
coeff[5][4][1][3] = 0;
coeff[5][4][1][4] = 0;
coeff[5][4][1][5] = 0;
coeff[5][4][1][6] = 0;
coeff[5][4][1][7] = 0;
coeff[5][4][1][8] = 0;

coeff[5][5][0][0] = 2278175820531;
coeff[5][5][0][1] = 6115366903537;
coeff[5][5][0][2] = 5088284968245;
coeff[5][5][0][3] = 8129181589039;
coeff[5][5][0][4] = 4927250153672;
coeff[5][5][0][5] = 726966144877;
coeff[5][5][0][6] = 2476195964910;
coeff[5][5][0][7] = 6764599985999;
coeff[5][5][0][8] = 4231932726;

coeff[5][5][1][0] = 6517917179832;
coeff[5][5][1][1] = 5429504614638;
coeff[5][5][1][2] = 4234746852751;
coeff[5][5][1][3] = 2112089840866;
coeff[5][5][1][4] = 4217472791108;
coeff[5][5][1][5] = 6733099956857;
coeff[5][5][1][6] = 8439745676103;
coeff[5][5][1][7] = 6426130780557;
coeff[5][5][1][8] = 107455168258;

coeff[6][0][0][0] = 1;
coeff[6][0][0][1] = 0;
coeff[6][0][0][2] = 0;
coeff[6][0][0][3] = 0;
coeff[6][0][0][4] = 0;
coeff[6][0][0][5] = 0;
coeff[6][0][0][6] = 0;
coeff[6][0][0][7] = 0;
coeff[6][0][0][8] = 0;

coeff[6][0][1][0] = 0;
coeff[6][0][1][1] = 0;
coeff[6][0][1][2] = 0;
coeff[6][0][1][3] = 0;
coeff[6][0][1][4] = 0;
coeff[6][0][1][5] = 0;
coeff[6][0][1][6] = 0;
coeff[6][0][1][7] = 0;
coeff[6][0][1][8] = 0;

coeff[6][1][0][0] = 8796093000362;
coeff[6][1][0][1] = 2748778495967;
coeff[6][1][0][2] = 526938798789;
coeff[6][1][0][3] = 1445178407698;
coeff[6][1][0][4] = 348629922573;
coeff[6][1][0][5] = 7460066101735;
coeff[6][1][0][6] = 2119848618805;
coeff[6][1][0][7] = 4394637744349;
coeff[6][1][0][8] = 111687100985;

coeff[6][1][1][0] = 0;
coeff[6][1][1][1] = 0;
coeff[6][1][1][2] = 0;
coeff[6][1][1][3] = 0;
coeff[6][1][1][4] = 0;
coeff[6][1][1][5] = 0;
coeff[6][1][1][6] = 0;
coeff[6][1][1][7] = 0;
coeff[6][1][1][8] = 0;

coeff[6][2][0][0] = 1;
coeff[6][2][0][1] = 0;
coeff[6][2][0][2] = 0;
coeff[6][2][0][3] = 0;
coeff[6][2][0][4] = 0;
coeff[6][2][0][5] = 0;
coeff[6][2][0][6] = 0;
coeff[6][2][0][7] = 0;
coeff[6][2][0][8] = 0;

coeff[6][2][1][0] = 0;
coeff[6][2][1][1] = 0;
coeff[6][2][1][2] = 0;
coeff[6][2][1][3] = 0;
coeff[6][2][1][4] = 0;
coeff[6][2][1][5] = 0;
coeff[6][2][1][6] = 0;
coeff[6][2][1][7] = 0;
coeff[6][2][1][8] = 0;

coeff[6][3][0][0] = 8796093000362;
coeff[6][3][0][1] = 2748778495967;
coeff[6][3][0][2] = 526938798789;
coeff[6][3][0][3] = 1445178407698;
coeff[6][3][0][4] = 348629922573;
coeff[6][3][0][5] = 7460066101735;
coeff[6][3][0][6] = 2119848618805;
coeff[6][3][0][7] = 4394637744349;
coeff[6][3][0][8] = 111687100985;

coeff[6][3][1][0] = 0;
coeff[6][3][1][1] = 0;
coeff[6][3][1][2] = 0;
coeff[6][3][1][3] = 0;
coeff[6][3][1][4] = 0;
coeff[6][3][1][5] = 0;
coeff[6][3][1][6] = 0;
coeff[6][3][1][7] = 0;
coeff[6][3][1][8] = 0;

coeff[6][4][0][0] = 1;
coeff[6][4][0][1] = 0;
coeff[6][4][0][2] = 0;
coeff[6][4][0][3] = 0;
coeff[6][4][0][4] = 0;
coeff[6][4][0][5] = 0;
coeff[6][4][0][6] = 0;
coeff[6][4][0][7] = 0;
coeff[6][4][0][8] = 0;

coeff[6][4][1][0] = 0;
coeff[6][4][1][1] = 0;
coeff[6][4][1][2] = 0;
coeff[6][4][1][3] = 0;
coeff[6][4][1][4] = 0;
coeff[6][4][1][5] = 0;
coeff[6][4][1][6] = 0;
coeff[6][4][1][7] = 0;
coeff[6][4][1][8] = 0;

coeff[6][5][0][0] = 8796093000362;
coeff[6][5][0][1] = 2748778495967;
coeff[6][5][0][2] = 526938798789;
coeff[6][5][0][3] = 1445178407698;
coeff[6][5][0][4] = 348629922573;
coeff[6][5][0][5] = 7460066101735;
coeff[6][5][0][6] = 2119848618805;
coeff[6][5][0][7] = 4394637744349;
coeff[6][5][0][8] = 111687100985;

coeff[6][5][1][0] = 0;
coeff[6][5][1][1] = 0;
coeff[6][5][1][2] = 0;
coeff[6][5][1][3] = 0;
coeff[6][5][1][4] = 0;
coeff[6][5][1][5] = 0;
coeff[6][5][1][6] = 0;
coeff[6][5][1][7] = 0;
coeff[6][5][1][8] = 0;

coeff[7][0][0][0] = 1;
coeff[7][0][0][1] = 0;
coeff[7][0][0][2] = 0;
coeff[7][0][0][3] = 0;
coeff[7][0][0][4] = 0;
coeff[7][0][0][5] = 0;
coeff[7][0][0][6] = 0;
coeff[7][0][0][7] = 0;
coeff[7][0][0][8] = 0;

coeff[7][0][1][0] = 0;
coeff[7][0][1][1] = 0;
coeff[7][0][1][2] = 0;
coeff[7][0][1][3] = 0;
coeff[7][0][1][4] = 0;
coeff[7][0][1][5] = 0;
coeff[7][0][1][6] = 0;
coeff[7][0][1][7] = 0;
coeff[7][0][1][8] = 0;

coeff[7][1][0][0] = 2278175820531;
coeff[7][1][0][1] = 6115366903537;
coeff[7][1][0][2] = 5088284968245;
coeff[7][1][0][3] = 8129181589039;
coeff[7][1][0][4] = 4927250153672;
coeff[7][1][0][5] = 726966144877;
coeff[7][1][0][6] = 2476195964910;
coeff[7][1][0][7] = 6764599985999;
coeff[7][1][0][8] = 4231932726;

coeff[7][1][1][0] = 6517917179832;
coeff[7][1][1][1] = 5429504614638;
coeff[7][1][1][2] = 4234746852751;
coeff[7][1][1][3] = 2112089840866;
coeff[7][1][1][4] = 4217472791108;
coeff[7][1][1][5] = 6733099956857;
coeff[7][1][1][6] = 8439745676103;
coeff[7][1][1][7] = 6426130780557;
coeff[7][1][1][8] = 107455168258;

coeff[7][2][0][0] = 0;
coeff[7][2][0][1] = 0;
coeff[7][2][0][2] = 0;
coeff[7][2][0][3] = 0;
coeff[7][2][0][4] = 0;
coeff[7][2][0][5] = 0;
coeff[7][2][0][6] = 0;
coeff[7][2][0][7] = 0;
coeff[7][2][0][8] = 0;

coeff[7][2][1][0] = 43692;
coeff[7][2][1][1] = 1374384390048;
coeff[7][2][1][2] = 5507500715325;
coeff[7][2][1][3] = 5164830044109;
coeff[7][2][1][4] = 8081740699602;
coeff[7][2][1][5] = 3659765840658;
coeff[7][2][1][6] = 144300128161;
coeff[7][2][1][7] = 4394637549586;
coeff[7][2][1][8] = 111687100985;

coeff[7][3][0][0] = 3316018568866;
coeff[7][3][0][1] = 5494359866831;
coeff[7][3][0][2] = 829224074489;
coeff[7][3][0][3] = 2495293119023;
coeff[7][3][0][4] = 3422922142614;
coeff[7][3][0][5] = 2441619419745;
coeff[7][3][0][6] = 1177724816182;
coeff[7][3][0][7] = 4420401567566;
coeff[7][3][0][8] = 82980365825;

coeff[7][3][1][0] = 3316018568866;
coeff[7][3][1][1] = 5494359866831;
coeff[7][3][1][2] = 829224074489;
coeff[7][3][1][3] = 2495293119023;
coeff[7][3][1][4] = 3422922142614;
coeff[7][3][1][5] = 2441619419745;
coeff[7][3][1][6] = 1177724816182;
coeff[7][3][1][7] = 4420401567566;
coeff[7][3][1][8] = 82980365825;

coeff[7][4][0][0] = 43693;
coeff[7][4][0][1] = 1374384390048;
coeff[7][4][0][2] = 5507500715325;
coeff[7][4][0][3] = 5164830044109;
coeff[7][4][0][4] = 8081740699602;
coeff[7][4][0][5] = 3659765840658;
coeff[7][4][0][6] = 144300128161;
coeff[7][4][0][7] = 4394637549586;
coeff[7][4][0][8] = 111687100985;

coeff[7][4][1][0] = 0;
coeff[7][4][1][1] = 0;
coeff[7][4][1][2] = 0;
coeff[7][4][1][3] = 0;
coeff[7][4][1][4] = 0;
coeff[7][4][1][5] = 0;
coeff[7][4][1][6] = 0;
coeff[7][4][1][7] = 0;
coeff[7][4][1][8] = 0;

coeff[7][5][0][0] = 5594194389397;
coeff[7][5][0][1] = 2813633748160;
coeff[7][5][0][2] = 5917509042735;
coeff[7][5][0][3] = 1828381685854;
coeff[7][5][0][4] = 8350172296287;
coeff[7][5][0][5] = 3168585564622;
coeff[7][5][0][6] = 3653920781092;
coeff[7][5][0][7] = 2388908531357;
coeff[7][5][0][8] = 87212298552;

coeff[7][5][1][0] = 3201898610966;
coeff[7][5][1][1] = 8731237770015;
coeff[7][5][1][2] = 3405522778261;
coeff[7][5][1][3] = 8412889744051;
coeff[7][5][1][4] = 794550648493;
coeff[7][5][1][5] = 4291480537112;
coeff[7][5][1][6] = 7262020859921;
coeff[7][5][1][7] = 2005729212991;
coeff[7][5][1][8] = 24474802433;

coeff[8][0][0][0] = 1;
coeff[8][0][0][1] = 0;
coeff[8][0][0][2] = 0;
coeff[8][0][0][3] = 0;
coeff[8][0][0][4] = 0;
coeff[8][0][0][5] = 0;
coeff[8][0][0][6] = 0;
coeff[8][0][0][7] = 0;
coeff[8][0][0][8] = 0;

coeff[8][0][1][0] = 0;
coeff[8][0][1][1] = 0;
coeff[8][0][1][2] = 0;
coeff[8][0][1][3] = 0;
coeff[8][0][1][4] = 0;
coeff[8][0][1][5] = 0;
coeff[8][0][1][6] = 0;
coeff[8][0][1][7] = 0;
coeff[8][0][1][8] = 0;

coeff[8][1][0][0] = 43692;
coeff[8][1][0][1] = 1374384390048;
coeff[8][1][0][2] = 5507500715325;
coeff[8][1][0][3] = 5164830044109;
coeff[8][1][0][4] = 8081740699602;
coeff[8][1][0][5] = 3659765840658;
coeff[8][1][0][6] = 144300128161;
coeff[8][1][0][7] = 4394637549586;
coeff[8][1][0][8] = 111687100985;

coeff[8][1][1][0] = 0;
coeff[8][1][1][1] = 0;
coeff[8][1][1][2] = 0;
coeff[8][1][1][3] = 0;
coeff[8][1][1][4] = 0;
coeff[8][1][1][5] = 0;
coeff[8][1][1][6] = 0;
coeff[8][1][1][7] = 0;
coeff[8][1][1][8] = 0;

coeff[8][2][0][0] = 8796092956670;
coeff[8][2][0][1] = 1374394105919;
coeff[8][2][0][2] = 3815531105672;
coeff[8][2][0][3] = 5076441385796;
coeff[8][2][0][4] = 1062982245178;
coeff[8][2][0][5] = 3800300261076;
coeff[8][2][0][6] = 1975548490644;
coeff[8][2][0][7] = 194763;
coeff[8][2][0][8] = 0;

coeff[8][2][1][0] = 0;
coeff[8][2][1][1] = 0;
coeff[8][2][1][2] = 0;
coeff[8][2][1][3] = 0;
coeff[8][2][1][4] = 0;
coeff[8][2][1][5] = 0;
coeff[8][2][1][6] = 0;
coeff[8][2][1][7] = 0;
coeff[8][2][1][8] = 0;

coeff[8][3][0][0] = 1;
coeff[8][3][0][1] = 0;
coeff[8][3][0][2] = 0;
coeff[8][3][0][3] = 0;
coeff[8][3][0][4] = 0;
coeff[8][3][0][5] = 0;
coeff[8][3][0][6] = 0;
coeff[8][3][0][7] = 0;
coeff[8][3][0][8] = 0;

coeff[8][3][1][0] = 0;
coeff[8][3][1][1] = 0;
coeff[8][3][1][2] = 0;
coeff[8][3][1][3] = 0;
coeff[8][3][1][4] = 0;
coeff[8][3][1][5] = 0;
coeff[8][3][1][6] = 0;
coeff[8][3][1][7] = 0;
coeff[8][3][1][8] = 0;

coeff[8][4][0][0] = 43692;
coeff[8][4][0][1] = 1374384390048;
coeff[8][4][0][2] = 5507500715325;
coeff[8][4][0][3] = 5164830044109;
coeff[8][4][0][4] = 8081740699602;
coeff[8][4][0][5] = 3659765840658;
coeff[8][4][0][6] = 144300128161;
coeff[8][4][0][7] = 4394637549586;
coeff[8][4][0][8] = 111687100985;

coeff[8][4][1][0] = 0;
coeff[8][4][1][1] = 0;
coeff[8][4][1][2] = 0;
coeff[8][4][1][3] = 0;
coeff[8][4][1][4] = 0;
coeff[8][4][1][5] = 0;
coeff[8][4][1][6] = 0;
coeff[8][4][1][7] = 0;
coeff[8][4][1][8] = 0;

coeff[8][5][0][0] = 8796092956670;
coeff[8][5][0][1] = 1374394105919;
coeff[8][5][0][2] = 3815531105672;
coeff[8][5][0][3] = 5076441385796;
coeff[8][5][0][4] = 1062982245178;
coeff[8][5][0][5] = 3800300261076;
coeff[8][5][0][6] = 1975548490644;
coeff[8][5][0][7] = 194763;
coeff[8][5][0][8] = 0;

coeff[8][5][1][0] = 0;
coeff[8][5][1][1] = 0;
coeff[8][5][1][2] = 0;
coeff[8][5][1][3] = 0;
coeff[8][5][1][4] = 0;
coeff[8][5][1][5] = 0;
coeff[8][5][1][6] = 0;
coeff[8][5][1][7] = 0;
coeff[8][5][1][8] = 0;

coeff[9][0][0][0] = 1;
coeff[9][0][0][1] = 0;
coeff[9][0][0][2] = 0;
coeff[9][0][0][3] = 0;
coeff[9][0][0][4] = 0;
coeff[9][0][0][5] = 0;
coeff[9][0][0][6] = 0;
coeff[9][0][0][7] = 0;
coeff[9][0][0][8] = 0;

coeff[9][0][1][0] = 0;
coeff[9][0][1][1] = 0;
coeff[9][0][1][2] = 0;
coeff[9][0][1][3] = 0;
coeff[9][0][1][4] = 0;
coeff[9][0][1][5] = 0;
coeff[9][0][1][6] = 0;
coeff[9][0][1][7] = 0;
coeff[9][0][1][8] = 0;

coeff[9][1][0][0] = 5480074431497;
coeff[9][1][0][1] = 6050511651344;
coeff[9][1][0][2] = 8493807746507;
coeff[9][1][0][3] = 7745978310882;
coeff[9][1][0][4] = 5721800802166;
coeff[9][1][0][5] = 5018446681989;
coeff[9][1][0][6] = 942123802623;
coeff[9][1][0][7] = 8770329198991;
coeff[9][1][0][8] = 28706735159;

coeff[9][1][1][0] = 3316018568866;
coeff[9][1][1][1] = 5494359866831;
coeff[9][1][1][2] = 829224074489;
coeff[9][1][1][3] = 2495293119023;
coeff[9][1][1][4] = 3422922142614;
coeff[9][1][1][5] = 2441619419745;
coeff[9][1][1][6] = 1177724816182;
coeff[9][1][1][7] = 4420401567566;
coeff[9][1][1][8] = 82980365825;

coeff[9][2][0][0] = 0;
coeff[9][2][0][1] = 0;
coeff[9][2][0][2] = 0;
coeff[9][2][0][3] = 0;
coeff[9][2][0][4] = 0;
coeff[9][2][0][5] = 0;
coeff[9][2][0][6] = 0;
coeff[9][2][0][7] = 0;
coeff[9][2][0][8] = 0;

coeff[9][2][1][0] = 1;
coeff[9][2][1][1] = 0;
coeff[9][2][1][2] = 0;
coeff[9][2][1][3] = 0;
coeff[9][2][1][4] = 0;
coeff[9][2][1][5] = 0;
coeff[9][2][1][6] = 0;
coeff[9][2][1][7] = 0;
coeff[9][2][1][8] = 0;

coeff[9][3][0][0] = 5480074431497;
coeff[9][3][0][1] = 6050511651344;
coeff[9][3][0][2] = 8493807746507;
coeff[9][3][0][3] = 7745978310882;
coeff[9][3][0][4] = 5721800802166;
coeff[9][3][0][5] = 5018446681989;
coeff[9][3][0][6] = 942123802623;
coeff[9][3][0][7] = 8770329198991;
coeff[9][3][0][8] = 28706735159;

coeff[9][3][1][0] = 5480074431497;
coeff[9][3][1][1] = 6050511651344;
coeff[9][3][1][2] = 8493807746507;
coeff[9][3][1][3] = 7745978310882;
coeff[9][3][1][4] = 5721800802166;
coeff[9][3][1][5] = 5018446681989;
coeff[9][3][1][6] = 942123802623;
coeff[9][3][1][7] = 8770329198991;
coeff[9][3][1][8] = 28706735159;

coeff[9][4][0][0] = 8796093000362;
coeff[9][4][0][1] = 2748778495967;
coeff[9][4][0][2] = 526938798789;
coeff[9][4][0][3] = 1445178407698;
coeff[9][4][0][4] = 348629922573;
coeff[9][4][0][5] = 7460066101735;
coeff[9][4][0][6] = 2119848618805;
coeff[9][4][0][7] = 4394637744349;
coeff[9][4][0][8] = 111687100985;

coeff[9][4][1][0] = 0;
coeff[9][4][1][1] = 0;
coeff[9][4][1][2] = 0;
coeff[9][4][1][3] = 0;
coeff[9][4][1][4] = 0;
coeff[9][4][1][5] = 0;
coeff[9][4][1][6] = 0;
coeff[9][4][1][7] = 0;
coeff[9][4][1][8] = 0;

coeff[9][5][0][0] = 3316018568866;
coeff[9][5][0][1] = 5494359866831;
coeff[9][5][0][2] = 829224074489;
coeff[9][5][0][3] = 2495293119023;
coeff[9][5][0][4] = 3422922142614;
coeff[9][5][0][5] = 2441619419745;
coeff[9][5][0][6] = 1177724816182;
coeff[9][5][0][7] = 4420401567566;
coeff[9][5][0][8] = 82980365825;

coeff[9][5][1][0] = 5480074431497;
coeff[9][5][1][1] = 6050511651344;
coeff[9][5][1][2] = 8493807746507;
coeff[9][5][1][3] = 7745978310882;
coeff[9][5][1][4] = 5721800802166;
coeff[9][5][1][5] = 5018446681989;
coeff[9][5][1][6] = 942123802623;
coeff[9][5][1][7] = 8770329198991;
coeff[9][5][1][8] = 28706735159;

coeff[10][0][0][0] = 1;
coeff[10][0][0][1] = 0;
coeff[10][0][0][2] = 0;
coeff[10][0][0][3] = 0;
coeff[10][0][0][4] = 0;
coeff[10][0][0][5] = 0;
coeff[10][0][0][6] = 0;
coeff[10][0][0][7] = 0;
coeff[10][0][0][8] = 0;

coeff[10][0][1][0] = 0;
coeff[10][0][1][1] = 0;
coeff[10][0][1][2] = 0;
coeff[10][0][1][3] = 0;
coeff[10][0][1][4] = 0;
coeff[10][0][1][5] = 0;
coeff[10][0][1][6] = 0;
coeff[10][0][1][7] = 0;
coeff[10][0][1][8] = 0;

coeff[10][1][0][0] = 43693;
coeff[10][1][0][1] = 1374384390048;
coeff[10][1][0][2] = 5507500715325;
coeff[10][1][0][3] = 5164830044109;
coeff[10][1][0][4] = 8081740699602;
coeff[10][1][0][5] = 3659765840658;
coeff[10][1][0][6] = 144300128161;
coeff[10][1][0][7] = 4394637549586;
coeff[10][1][0][8] = 111687100985;

coeff[10][1][1][0] = 0;
coeff[10][1][1][1] = 0;
coeff[10][1][1][2] = 0;
coeff[10][1][1][3] = 0;
coeff[10][1][1][4] = 0;
coeff[10][1][1][5] = 0;
coeff[10][1][1][6] = 0;
coeff[10][1][1][7] = 0;
coeff[10][1][1][8] = 0;

coeff[10][2][0][0] = 43692;
coeff[10][2][0][1] = 1374384390048;
coeff[10][2][0][2] = 5507500715325;
coeff[10][2][0][3] = 5164830044109;
coeff[10][2][0][4] = 8081740699602;
coeff[10][2][0][5] = 3659765840658;
coeff[10][2][0][6] = 144300128161;
coeff[10][2][0][7] = 4394637549586;
coeff[10][2][0][8] = 111687100985;

coeff[10][2][1][0] = 0;
coeff[10][2][1][1] = 0;
coeff[10][2][1][2] = 0;
coeff[10][2][1][3] = 0;
coeff[10][2][1][4] = 0;
coeff[10][2][1][5] = 0;
coeff[10][2][1][6] = 0;
coeff[10][2][1][7] = 0;
coeff[10][2][1][8] = 0;

coeff[10][3][0][0] = 8796093000362;
coeff[10][3][0][1] = 2748778495967;
coeff[10][3][0][2] = 526938798789;
coeff[10][3][0][3] = 1445178407698;
coeff[10][3][0][4] = 348629922573;
coeff[10][3][0][5] = 7460066101735;
coeff[10][3][0][6] = 2119848618805;
coeff[10][3][0][7] = 4394637744349;
coeff[10][3][0][8] = 111687100985;

coeff[10][3][1][0] = 0;
coeff[10][3][1][1] = 0;
coeff[10][3][1][2] = 0;
coeff[10][3][1][3] = 0;
coeff[10][3][1][4] = 0;
coeff[10][3][1][5] = 0;
coeff[10][3][1][6] = 0;
coeff[10][3][1][7] = 0;
coeff[10][3][1][8] = 0;

coeff[10][4][0][0] = 8796092956670;
coeff[10][4][0][1] = 1374394105919;
coeff[10][4][0][2] = 3815531105672;
coeff[10][4][0][3] = 5076441385796;
coeff[10][4][0][4] = 1062982245178;
coeff[10][4][0][5] = 3800300261076;
coeff[10][4][0][6] = 1975548490644;
coeff[10][4][0][7] = 194763;
coeff[10][4][0][8] = 0;

coeff[10][4][1][0] = 0;
coeff[10][4][1][1] = 0;
coeff[10][4][1][2] = 0;
coeff[10][4][1][3] = 0;
coeff[10][4][1][4] = 0;
coeff[10][4][1][5] = 0;
coeff[10][4][1][6] = 0;
coeff[10][4][1][7] = 0;
coeff[10][4][1][8] = 0;

coeff[10][5][0][0] = 8796092956671;
coeff[10][5][0][1] = 1374394105919;
coeff[10][5][0][2] = 3815531105672;
coeff[10][5][0][3] = 5076441385796;
coeff[10][5][0][4] = 1062982245178;
coeff[10][5][0][5] = 3800300261076;
coeff[10][5][0][6] = 1975548490644;
coeff[10][5][0][7] = 194763;
coeff[10][5][0][8] = 0;

coeff[10][5][1][0] = 0;
coeff[10][5][1][1] = 0;
coeff[10][5][1][2] = 0;
coeff[10][5][1][3] = 0;
coeff[10][5][1][4] = 0;
coeff[10][5][1][5] = 0;
coeff[10][5][1][6] = 0;
coeff[10][5][1][7] = 0;
coeff[10][5][1][8] = 0;

coeff[11][0][0][0] = 1;
coeff[11][0][0][1] = 0;
coeff[11][0][0][2] = 0;
coeff[11][0][0][3] = 0;
coeff[11][0][0][4] = 0;
coeff[11][0][0][5] = 0;
coeff[11][0][0][6] = 0;
coeff[11][0][0][7] = 0;
coeff[11][0][0][8] = 0;

coeff[11][0][1][0] = 0;
coeff[11][0][1][1] = 0;
coeff[11][0][1][2] = 0;
coeff[11][0][1][3] = 0;
coeff[11][0][1][4] = 0;
coeff[11][0][1][5] = 0;
coeff[11][0][1][6] = 0;
coeff[11][0][1][7] = 0;
coeff[11][0][1][8] = 0;

coeff[11][1][0][0] = 3201898610966;
coeff[11][1][0][1] = 8731237770015;
coeff[11][1][0][2] = 3405522778261;
coeff[11][1][0][3] = 8412889744051;
coeff[11][1][0][4] = 794550648493;
coeff[11][1][0][5] = 4291480537112;
coeff[11][1][0][6] = 7262020859921;
coeff[11][1][0][7] = 2005729212991;
coeff[11][1][0][8] = 24474802433;

coeff[11][1][1][0] = 5594194389397;
coeff[11][1][1][1] = 2813633748160;
coeff[11][1][1][2] = 5917509042735;
coeff[11][1][1][3] = 1828381685854;
coeff[11][1][1][4] = 8350172296287;
coeff[11][1][1][5] = 3168585564622;
coeff[11][1][1][6] = 3653920781092;
coeff[11][1][1][7] = 2388908531357;
coeff[11][1][1][8] = 87212298552;

coeff[11][2][0][0] = 0;
coeff[11][2][0][1] = 0;
coeff[11][2][0][2] = 0;
coeff[11][2][0][3] = 0;
coeff[11][2][0][4] = 0;
coeff[11][2][0][5] = 0;
coeff[11][2][0][6] = 0;
coeff[11][2][0][7] = 0;
coeff[11][2][0][8] = 0;

coeff[11][2][1][0] = 8796092956670;
coeff[11][2][1][1] = 1374394105919;
coeff[11][2][1][2] = 3815531105672;
coeff[11][2][1][3] = 5076441385796;
coeff[11][2][1][4] = 1062982245178;
coeff[11][2][1][5] = 3800300261076;
coeff[11][2][1][6] = 1975548490644;
coeff[11][2][1][7] = 194763;
coeff[11][2][1][8] = 0;

coeff[11][3][0][0] = 3316018568866;
coeff[11][3][0][1] = 5494359866831;
coeff[11][3][0][2] = 829224074489;
coeff[11][3][0][3] = 2495293119023;
coeff[11][3][0][4] = 3422922142614;
coeff[11][3][0][5] = 2441619419745;
coeff[11][3][0][6] = 1177724816182;
coeff[11][3][0][7] = 4420401567566;
coeff[11][3][0][8] = 82980365825;

coeff[11][3][1][0] = 3316018568866;
coeff[11][3][1][1] = 5494359866831;
coeff[11][3][1][2] = 829224074489;
coeff[11][3][1][3] = 2495293119023;
coeff[11][3][1][4] = 3422922142614;
coeff[11][3][1][5] = 2441619419745;
coeff[11][3][1][6] = 1177724816182;
coeff[11][3][1][7] = 4420401567566;
coeff[11][3][1][8] = 82980365825;

coeff[11][4][0][0] = 8796092956671;
coeff[11][4][0][1] = 1374394105919;
coeff[11][4][0][2] = 3815531105672;
coeff[11][4][0][3] = 5076441385796;
coeff[11][4][0][4] = 1062982245178;
coeff[11][4][0][5] = 3800300261076;
coeff[11][4][0][6] = 1975548490644;
coeff[11][4][0][7] = 194763;
coeff[11][4][0][8] = 0;

coeff[11][4][1][0] = 0;
coeff[11][4][1][1] = 0;
coeff[11][4][1][2] = 0;
coeff[11][4][1][3] = 0;
coeff[11][4][1][4] = 0;
coeff[11][4][1][5] = 0;
coeff[11][4][1][6] = 0;
coeff[11][4][1][7] = 0;
coeff[11][4][1][8] = 0;

coeff[11][5][0][0] = 6517917179832;
coeff[11][5][0][1] = 5429504614638;
coeff[11][5][0][2] = 4234746852751;
coeff[11][5][0][3] = 2112089840866;
coeff[11][5][0][4] = 4217472791108;
coeff[11][5][0][5] = 6733099956857;
coeff[11][5][0][6] = 8439745676103;
coeff[11][5][0][7] = 6426130780557;
coeff[11][5][0][8] = 107455168258;

coeff[11][5][1][0] = 2278175820531;
coeff[11][5][1][1] = 6115366903537;
coeff[11][5][1][2] = 5088284968245;
coeff[11][5][1][3] = 8129181589039;
coeff[11][5][1][4] = 4927250153672;
coeff[11][5][1][5] = 726966144877;
coeff[11][5][1][6] = 2476195964910;
coeff[11][5][1][7] = 6764599985999;
coeff[11][5][1][8] = 4231932726;
}
// PRECOMPUTE FROM PYTHON
if( n == 77 && k == 5 ){
coeff[0][0][0][0] = 1;
coeff[0][0][0][1] = 0;
coeff[0][0][0][2] = 0;
coeff[0][0][0][3] = 0;
coeff[0][0][0][4] = 0;

coeff[0][0][1][0] = 0;
coeff[0][0][1][1] = 0;
coeff[0][0][1][2] = 0;
coeff[0][0][1][3] = 0;
coeff[0][0][1][4] = 0;

coeff[0][1][0][0] = 1;
coeff[0][1][0][1] = 0;
coeff[0][1][0][2] = 0;
coeff[0][1][0][3] = 0;
coeff[0][1][0][4] = 0;

coeff[0][1][1][0] = 0;
coeff[0][1][1][1] = 0;
coeff[0][1][1][2] = 0;
coeff[0][1][1][3] = 0;
coeff[0][1][1][4] = 0;

coeff[0][2][0][0] = 1;
coeff[0][2][0][1] = 0;
coeff[0][2][0][2] = 0;
coeff[0][2][0][3] = 0;
coeff[0][2][0][4] = 0;

coeff[0][2][1][0] = 0;
coeff[0][2][1][1] = 0;
coeff[0][2][1][2] = 0;
coeff[0][2][1][3] = 0;
coeff[0][2][1][4] = 0;

coeff[0][3][0][0] = 1;
coeff[0][3][0][1] = 0;
coeff[0][3][0][2] = 0;
coeff[0][3][0][3] = 0;
coeff[0][3][0][4] = 0;

coeff[0][3][1][0] = 0;
coeff[0][3][1][1] = 0;
coeff[0][3][1][2] = 0;
coeff[0][3][1][3] = 0;
coeff[0][3][1][4] = 0;

coeff[0][4][0][0] = 1;
coeff[0][4][0][1] = 0;
coeff[0][4][0][2] = 0;
coeff[0][4][0][3] = 0;
coeff[0][4][0][4] = 0;

coeff[0][4][1][0] = 0;
coeff[0][4][1][1] = 0;
coeff[0][4][1][2] = 0;
coeff[0][4][1][3] = 0;
coeff[0][4][1][4] = 0;

coeff[0][5][0][0] = 1;
coeff[0][5][0][1] = 0;
coeff[0][5][0][2] = 0;
coeff[0][5][0][3] = 0;
coeff[0][5][0][4] = 0;

coeff[0][5][1][0] = 0;
coeff[0][5][1][1] = 0;
coeff[0][5][1][2] = 0;
coeff[0][5][1][3] = 0;
coeff[0][5][1][4] = 0;

coeff[1][0][0][0] = 1;
coeff[1][0][0][1] = 0;
coeff[1][0][0][2] = 0;
coeff[1][0][0][3] = 0;
coeff[1][0][0][4] = 0;

coeff[1][0][1][0] = 0;
coeff[1][0][1][1] = 0;
coeff[1][0][1][2] = 0;
coeff[1][0][1][3] = 0;
coeff[1][0][1][4] = 0;

coeff[1][1][0][0] = 5857780092113332166584;
coeff[1][1][0][1] = 27555843941809817067324;
coeff[1][1][0][2] = 147854065178978622174689;
coeff[1][1][0][3] = 16480454558870818237447;
coeff[1][1][0][4] = 7384262935318800792611;

coeff[1][1][1][0] = 145252903046657682983667;
coeff[1][1][1][1] = 74116926119540154854243;
coeff[1][1][1][2] = 9107065692449161851484;
coeff[1][1][1][3] = 93599087433168491987599;
coeff[1][1][1][4] = 290816202565522499726;

coeff[1][2][0][0] = 0;
coeff[1][2][0][1] = 0;
coeff[1][2][0][2] = 0;
coeff[1][2][0][3] = 0;
coeff[1][2][0][4] = 0;

coeff[1][2][1][0] = 151070474438347898006188;
coeff[1][2][1][1] = 140545516233852444113487;
coeff[1][2][1][2] = 148207872154034921691459;
coeff[1][2][1][3] = 21270016808265337543434;
coeff[1][2][1][4] = 7675079137884323290816;

coeff[1][3][0][0] = 28127254136958132407305;
coeff[1][3][0][1] = 21583712728585527334752;
coeff[1][3][0][2] = 51663912557182296229312;
coeff[1][3][0][3] = 17835324027044025627323;
coeff[1][3][0][4] = 1972711818993931957891;

coeff[1][3][1][0] = 28127254136958132407305;
coeff[1][3][1][1] = 21583712728585527334752;
coeff[1][3][1][2] = 51663912557182296229312;
coeff[1][3][1][3] = 17835324027044025627323;
coeff[1][3][1][4] = 1972711818993931957891;

coeff[1][4][0][0] = 151070474438347898006189;
coeff[1][4][0][1] = 140545516233852444113487;
coeff[1][4][0][2] = 148207872154034921691459;
coeff[1][4][0][3] = 21270016808265337543434;
coeff[1][4][0][4] = 7675079137884323290816;

coeff[1][4][1][0] = 0;
coeff[1][4][1][1] = 0;
coeff[1][4][1][2] = 0;
coeff[1][4][1][3] = 0;
coeff[1][4][1][4] = 0;

coeff[1][5][0][0] = 33990078542129096261910;
coeff[1][5][0][1] = 98582514060874019318780;
coeff[1][5][0][2] = 42556846864733134377827;
coeff[1][5][0][3] = 75351964045704180477996;
coeff[1][5][0][4] = 1681895616428409458164;

coeff[1][5][1][0] = 117120604596641918888341;
coeff[1][5][1][1] = 3090256000475952602787;
coeff[1][5][1][2] = 114404284006694649648346;
coeff[1][5][1][3] = 34727577946335129747050;
coeff[1][5][1][4] = 5993183521455913834173;

coeff[2][0][0][0] = 1;
coeff[2][0][0][1] = 0;
coeff[2][0][0][2] = 0;
coeff[2][0][0][3] = 0;
coeff[2][0][0][4] = 0;

coeff[2][0][1][0] = 0;
coeff[2][0][1][1] = 0;
coeff[2][0][1][2] = 0;
coeff[2][0][1][3] = 0;
coeff[2][0][1][4] = 0;

coeff[2][1][0][0] = 40208700423117144063;
coeff[2][1][0][1] = 112242981279326174646352;
coeff[2][1][0][2] = 8753258717392862334713;
coeff[2][1][0][3] = 88809525183773972681612;
coeff[2][1][0][4] = 1521;

coeff[2][1][1][0] = 0;
coeff[2][1][1][1] = 0;
coeff[2][1][1][2] = 0;
coeff[2][1][1][3] = 0;
coeff[2][1][1][4] = 0;

coeff[2][2][0][0] = 40208700423117144062;
coeff[2][2][0][1] = 112242981279326174646352;
coeff[2][2][0][2] = 8753258717392862334713;
coeff[2][2][0][3] = 88809525183773972681612;
coeff[2][2][0][4] = 1521;

coeff[2][2][1][0] = 0;
coeff[2][2][1][1] = 0;
coeff[2][2][1][2] = 0;
coeff[2][2][1][3] = 0;
coeff[2][2][1][4] = 0;

coeff[2][3][0][0] = 151110683138771015150250;
coeff[2][3][0][1] = 101672770061349971921567;
coeff[2][3][0][2] = 5845403419599137187901;
coeff[2][3][0][3] = 110079541992039310225047;
coeff[2][3][0][4] = 7675079137884323292337;

coeff[2][3][1][0] = 0;
coeff[2][3][1][1] = 0;
coeff[2][3][1][2] = 0;
coeff[2][3][1][3] = 0;
coeff[2][3][1][4] = 0;

coeff[2][4][0][0] = 151070474438347898006188;
coeff[2][4][0][1] = 140545516233852444113487;
coeff[2][4][0][2] = 148207872154034921691459;
coeff[2][4][0][3] = 21270016808265337543434;
coeff[2][4][0][4] = 7675079137884323290816;

coeff[2][4][1][0] = 0;
coeff[2][4][1][1] = 0;
coeff[2][4][1][2] = 0;
coeff[2][4][1][3] = 0;
coeff[2][4][1][4] = 0;

coeff[2][5][0][0] = 151070474438347898006189;
coeff[2][5][0][1] = 140545516233852444113487;
coeff[2][5][0][2] = 148207872154034921691459;
coeff[2][5][0][3] = 21270016808265337543434;
coeff[2][5][0][4] = 7675079137884323290816;

coeff[2][5][1][0] = 0;
coeff[2][5][1][1] = 0;
coeff[2][5][1][2] = 0;
coeff[2][5][1][3] = 0;
coeff[2][5][1][4] = 0;

coeff[3][0][0][0] = 1;
coeff[3][0][0][1] = 0;
coeff[3][0][0][2] = 0;
coeff[3][0][0][3] = 0;
coeff[3][0][0][4] = 0;

coeff[3][0][1][0] = 0;
coeff[3][0][1][1] = 0;
coeff[3][0][1][2] = 0;
coeff[3][0][1][3] = 0;
coeff[3][0][1][4] = 0;

coeff[3][1][0][0] = 122983429001812882742946;
coeff[3][1][0][1] = 80089057332764444586815;
coeff[3][1][0][2] = 105297218314245487796861;
coeff[3][1][0][3] = 92244217964995284597723;
coeff[3][1][0][4] = 5702367318890391334446;

coeff[3][1][1][0] = 28127254136958132407305;
coeff[3][1][1][1] = 21583712728585527334752;
coeff[3][1][1][2] = 51663912557182296229312;
coeff[3][1][1][3] = 17835324027044025627323;
coeff[3][1][1][4] = 1972711818993931957891;

coeff[3][2][0][0] = 0;
coeff[3][2][0][1] = 0;
coeff[3][2][0][2] = 0;
coeff[3][2][0][3] = 0;
coeff[3][2][0][4] = 0;

coeff[3][2][1][0] = 1;
coeff[3][2][1][1] = 0;
coeff[3][2][1][2] = 0;
coeff[3][2][1][3] = 0;
coeff[3][2][1][4] = 0;

coeff[3][3][0][0] = 122983429001812882742946;
coeff[3][3][0][1] = 80089057332764444586815;
coeff[3][3][0][2] = 105297218314245487796861;
coeff[3][3][0][3] = 92244217964995284597723;
coeff[3][3][0][4] = 5702367318890391334446;

coeff[3][3][1][0] = 122983429001812882742946;
coeff[3][3][1][1] = 80089057332764444586815;
coeff[3][3][1][2] = 105297218314245487796861;
coeff[3][3][1][3] = 92244217964995284597723;
coeff[3][3][1][4] = 5702367318890391334446;

coeff[3][4][0][0] = 151110683138771015150250;
coeff[3][4][0][1] = 101672770061349971921567;
coeff[3][4][0][2] = 5845403419599137187901;
coeff[3][4][0][3] = 110079541992039310225047;
coeff[3][4][0][4] = 7675079137884323292337;

coeff[3][4][1][0] = 0;
coeff[3][4][1][1] = 0;
coeff[3][4][1][2] = 0;
coeff[3][4][1][3] = 0;
coeff[3][4][1][4] = 0;

coeff[3][5][0][0] = 28127254136958132407305;
coeff[3][5][0][1] = 21583712728585527334752;
coeff[3][5][0][2] = 51663912557182296229312;
coeff[3][5][0][3] = 17835324027044025627323;
coeff[3][5][0][4] = 1972711818993931957891;

coeff[3][5][1][0] = 122983429001812882742946;
coeff[3][5][1][1] = 80089057332764444586815;
coeff[3][5][1][2] = 105297218314245487796861;
coeff[3][5][1][3] = 92244217964995284597723;
coeff[3][5][1][4] = 5702367318890391334446;

coeff[4][0][0][0] = 1;
coeff[4][0][0][1] = 0;
coeff[4][0][0][2] = 0;
coeff[4][0][0][3] = 0;
coeff[4][0][0][4] = 0;

coeff[4][0][1][0] = 0;
coeff[4][0][1][1] = 0;
coeff[4][0][1][2] = 0;
coeff[4][0][1][3] = 0;
coeff[4][0][1][4] = 0;

coeff[4][1][0][0] = 40208700423117144062;
coeff[4][1][0][1] = 112242981279326174646352;
coeff[4][1][0][2] = 8753258717392862334713;
coeff[4][1][0][3] = 88809525183773972681612;
coeff[4][1][0][4] = 1521;

coeff[4][1][1][0] = 0;
coeff[4][1][1][1] = 0;
coeff[4][1][1][2] = 0;
coeff[4][1][1][3] = 0;
coeff[4][1][1][4] = 0;

coeff[4][2][0][0] = 151070474438347898006188;
coeff[4][2][0][1] = 140545516233852444113487;
coeff[4][2][0][2] = 148207872154034921691459;
coeff[4][2][0][3] = 21270016808265337543434;
coeff[4][2][0][4] = 7675079137884323290816;

coeff[4][2][1][0] = 0;
coeff[4][2][1][1] = 0;
coeff[4][2][1][2] = 0;
coeff[4][2][1][3] = 0;
coeff[4][2][1][4] = 0;

coeff[4][3][0][0] = 1;
coeff[4][3][0][1] = 0;
coeff[4][3][0][2] = 0;
coeff[4][3][0][3] = 0;
coeff[4][3][0][4] = 0;

coeff[4][3][1][0] = 0;
coeff[4][3][1][1] = 0;
coeff[4][3][1][2] = 0;
coeff[4][3][1][3] = 0;
coeff[4][3][1][4] = 0;

coeff[4][4][0][0] = 40208700423117144062;
coeff[4][4][0][1] = 112242981279326174646352;
coeff[4][4][0][2] = 8753258717392862334713;
coeff[4][4][0][3] = 88809525183773972681612;
coeff[4][4][0][4] = 1521;

coeff[4][4][1][0] = 0;
coeff[4][4][1][1] = 0;
coeff[4][4][1][2] = 0;
coeff[4][4][1][3] = 0;
coeff[4][4][1][4] = 0;

coeff[4][5][0][0] = 151070474438347898006188;
coeff[4][5][0][1] = 140545516233852444113487;
coeff[4][5][0][2] = 148207872154034921691459;
coeff[4][5][0][3] = 21270016808265337543434;
coeff[4][5][0][4] = 7675079137884323290816;

coeff[4][5][1][0] = 0;
coeff[4][5][1][1] = 0;
coeff[4][5][1][2] = 0;
coeff[4][5][1][3] = 0;
coeff[4][5][1][4] = 0;

coeff[5][0][0][0] = 1;
coeff[5][0][0][1] = 0;
coeff[5][0][0][2] = 0;
coeff[5][0][0][3] = 0;
coeff[5][0][0][4] = 0;

coeff[5][0][1][0] = 0;
coeff[5][0][1][1] = 0;
coeff[5][0][1][2] = 0;
coeff[5][0][1][3] = 0;
coeff[5][0][1][4] = 0;

coeff[5][1][0][0] = 117120604596641918888341;
coeff[5][1][0][1] = 3090256000475952602787;
coeff[5][1][0][2] = 114404284006694649648346;
coeff[5][1][0][3] = 34727577946335129747050;
coeff[5][1][0][4] = 5993183521455913834173;

coeff[5][1][1][0] = 33990078542129096261910;
coeff[5][1][1][1] = 98582514060874019318780;
coeff[5][1][1][2] = 42556846864733134377827;
coeff[5][1][1][3] = 75351964045704180477996;
coeff[5][1][1][4] = 1681895616428409458164;

coeff[5][2][0][0] = 0;
coeff[5][2][0][1] = 0;
coeff[5][2][0][2] = 0;
coeff[5][2][0][3] = 0;
coeff[5][2][0][4] = 0;

coeff[5][2][1][0] = 40208700423117144062;
coeff[5][2][1][1] = 112242981279326174646352;
coeff[5][2][1][2] = 8753258717392862334713;
coeff[5][2][1][3] = 88809525183773972681612;
coeff[5][2][1][4] = 1521;

coeff[5][3][0][0] = 28127254136958132407305;
coeff[5][3][0][1] = 21583712728585527334752;
coeff[5][3][0][2] = 51663912557182296229312;
coeff[5][3][0][3] = 17835324027044025627323;
coeff[5][3][0][4] = 1972711818993931957891;

coeff[5][3][1][0] = 28127254136958132407305;
coeff[5][3][1][1] = 21583712728585527334752;
coeff[5][3][1][2] = 51663912557182296229312;
coeff[5][3][1][3] = 17835324027044025627323;
coeff[5][3][1][4] = 1972711818993931957891;

coeff[5][4][0][0] = 40208700423117144063;
coeff[5][4][0][1] = 112242981279326174646352;
coeff[5][4][0][2] = 8753258717392862334713;
coeff[5][4][0][3] = 88809525183773972681612;
coeff[5][4][0][4] = 1521;

coeff[5][4][1][0] = 0;
coeff[5][4][1][1] = 0;
coeff[5][4][1][2] = 0;
coeff[5][4][1][3] = 0;
coeff[5][4][1][4] = 0;

coeff[5][5][0][0] = 145252903046657682983667;
coeff[5][5][0][1] = 74116926119540154854243;
coeff[5][5][0][2] = 9107065692449161851484;
coeff[5][5][0][3] = 93599087433168491987599;
coeff[5][5][0][4] = 290816202565522499726;

coeff[5][5][1][0] = 5857780092113332166584;
coeff[5][5][1][1] = 27555843941809817067324;
coeff[5][5][1][2] = 147854065178978622174689;
coeff[5][5][1][3] = 16480454558870818237447;
coeff[5][5][1][4] = 7384262935318800792611;

coeff[6][0][0][0] = 1;
coeff[6][0][0][1] = 0;
coeff[6][0][0][2] = 0;
coeff[6][0][0][3] = 0;
coeff[6][0][0][4] = 0;

coeff[6][0][1][0] = 0;
coeff[6][0][1][1] = 0;
coeff[6][0][1][2] = 0;
coeff[6][0][1][3] = 0;
coeff[6][0][1][4] = 0;

coeff[6][1][0][0] = 151110683138771015150250;
coeff[6][1][0][1] = 101672770061349971921567;
coeff[6][1][0][2] = 5845403419599137187901;
coeff[6][1][0][3] = 110079541992039310225047;
coeff[6][1][0][4] = 7675079137884323292337;

coeff[6][1][1][0] = 0;
coeff[6][1][1][1] = 0;
coeff[6][1][1][2] = 0;
coeff[6][1][1][3] = 0;
coeff[6][1][1][4] = 0;

coeff[6][2][0][0] = 1;
coeff[6][2][0][1] = 0;
coeff[6][2][0][2] = 0;
coeff[6][2][0][3] = 0;
coeff[6][2][0][4] = 0;

coeff[6][2][1][0] = 0;
coeff[6][2][1][1] = 0;
coeff[6][2][1][2] = 0;
coeff[6][2][1][3] = 0;
coeff[6][2][1][4] = 0;

coeff[6][3][0][0] = 151110683138771015150250;
coeff[6][3][0][1] = 101672770061349971921567;
coeff[6][3][0][2] = 5845403419599137187901;
coeff[6][3][0][3] = 110079541992039310225047;
coeff[6][3][0][4] = 7675079137884323292337;

coeff[6][3][1][0] = 0;
coeff[6][3][1][1] = 0;
coeff[6][3][1][2] = 0;
coeff[6][3][1][3] = 0;
coeff[6][3][1][4] = 0;

coeff[6][4][0][0] = 1;
coeff[6][4][0][1] = 0;
coeff[6][4][0][2] = 0;
coeff[6][4][0][3] = 0;
coeff[6][4][0][4] = 0;

coeff[6][4][1][0] = 0;
coeff[6][4][1][1] = 0;
coeff[6][4][1][2] = 0;
coeff[6][4][1][3] = 0;
coeff[6][4][1][4] = 0;

coeff[6][5][0][0] = 151110683138771015150250;
coeff[6][5][0][1] = 101672770061349971921567;
coeff[6][5][0][2] = 5845403419599137187901;
coeff[6][5][0][3] = 110079541992039310225047;
coeff[6][5][0][4] = 7675079137884323292337;

coeff[6][5][1][0] = 0;
coeff[6][5][1][1] = 0;
coeff[6][5][1][2] = 0;
coeff[6][5][1][3] = 0;
coeff[6][5][1][4] = 0;

coeff[7][0][0][0] = 1;
coeff[7][0][0][1] = 0;
coeff[7][0][0][2] = 0;
coeff[7][0][0][3] = 0;
coeff[7][0][0][4] = 0;

coeff[7][0][1][0] = 0;
coeff[7][0][1][1] = 0;
coeff[7][0][1][2] = 0;
coeff[7][0][1][3] = 0;
coeff[7][0][1][4] = 0;

coeff[7][1][0][0] = 145252903046657682983667;
coeff[7][1][0][1] = 74116926119540154854243;
coeff[7][1][0][2] = 9107065692449161851484;
coeff[7][1][0][3] = 93599087433168491987599;
coeff[7][1][0][4] = 290816202565522499726;

coeff[7][1][1][0] = 5857780092113332166584;
coeff[7][1][1][1] = 27555843941809817067324;
coeff[7][1][1][2] = 147854065178978622174689;
coeff[7][1][1][3] = 16480454558870818237447;
coeff[7][1][1][4] = 7384262935318800792611;

coeff[7][2][0][0] = 0;
coeff[7][2][0][1] = 0;
coeff[7][2][0][2] = 0;
coeff[7][2][0][3] = 0;
coeff[7][2][0][4] = 0;

coeff[7][2][1][0] = 151070474438347898006188;
coeff[7][2][1][1] = 140545516233852444113487;
coeff[7][2][1][2] = 148207872154034921691459;
coeff[7][2][1][3] = 21270016808265337543434;
coeff[7][2][1][4] = 7675079137884323290816;

coeff[7][3][0][0] = 122983429001812882742946;
coeff[7][3][0][1] = 80089057332764444586815;
coeff[7][3][0][2] = 105297218314245487796861;
coeff[7][3][0][3] = 92244217964995284597723;
coeff[7][3][0][4] = 5702367318890391334446;

coeff[7][3][1][0] = 122983429001812882742946;
coeff[7][3][1][1] = 80089057332764444586815;
coeff[7][3][1][2] = 105297218314245487796861;
coeff[7][3][1][3] = 92244217964995284597723;
coeff[7][3][1][4] = 5702367318890391334446;

coeff[7][4][0][0] = 151070474438347898006189;
coeff[7][4][0][1] = 140545516233852444113487;
coeff[7][4][0][2] = 148207872154034921691459;
coeff[7][4][0][3] = 21270016808265337543434;
coeff[7][4][0][4] = 7675079137884323290816;

coeff[7][4][1][0] = 0;
coeff[7][4][1][1] = 0;
coeff[7][4][1][2] = 0;
coeff[7][4][1][3] = 0;
coeff[7][4][1][4] = 0;

coeff[7][5][0][0] = 117120604596641918888341;
coeff[7][5][0][1] = 3090256000475952602787;
coeff[7][5][0][2] = 114404284006694649648346;
coeff[7][5][0][3] = 34727577946335129747050;
coeff[7][5][0][4] = 5993183521455913834173;

coeff[7][5][1][0] = 33990078542129096261910;
coeff[7][5][1][1] = 98582514060874019318780;
coeff[7][5][1][2] = 42556846864733134377827;
coeff[7][5][1][3] = 75351964045704180477996;
coeff[7][5][1][4] = 1681895616428409458164;

coeff[8][0][0][0] = 1;
coeff[8][0][0][1] = 0;
coeff[8][0][0][2] = 0;
coeff[8][0][0][3] = 0;
coeff[8][0][0][4] = 0;

coeff[8][0][1][0] = 0;
coeff[8][0][1][1] = 0;
coeff[8][0][1][2] = 0;
coeff[8][0][1][3] = 0;
coeff[8][0][1][4] = 0;

coeff[8][1][0][0] = 151070474438347898006188;
coeff[8][1][0][1] = 140545516233852444113487;
coeff[8][1][0][2] = 148207872154034921691459;
coeff[8][1][0][3] = 21270016808265337543434;
coeff[8][1][0][4] = 7675079137884323290816;

coeff[8][1][1][0] = 0;
coeff[8][1][1][1] = 0;
coeff[8][1][1][2] = 0;
coeff[8][1][1][3] = 0;
coeff[8][1][1][4] = 0;

coeff[8][2][0][0] = 40208700423117144062;
coeff[8][2][0][1] = 112242981279326174646352;
coeff[8][2][0][2] = 8753258717392862334713;
coeff[8][2][0][3] = 88809525183773972681612;
coeff[8][2][0][4] = 1521;

coeff[8][2][1][0] = 0;
coeff[8][2][1][1] = 0;
coeff[8][2][1][2] = 0;
coeff[8][2][1][3] = 0;
coeff[8][2][1][4] = 0;

coeff[8][3][0][0] = 1;
coeff[8][3][0][1] = 0;
coeff[8][3][0][2] = 0;
coeff[8][3][0][3] = 0;
coeff[8][3][0][4] = 0;

coeff[8][3][1][0] = 0;
coeff[8][3][1][1] = 0;
coeff[8][3][1][2] = 0;
coeff[8][3][1][3] = 0;
coeff[8][3][1][4] = 0;

coeff[8][4][0][0] = 151070474438347898006188;
coeff[8][4][0][1] = 140545516233852444113487;
coeff[8][4][0][2] = 148207872154034921691459;
coeff[8][4][0][3] = 21270016808265337543434;
coeff[8][4][0][4] = 7675079137884323290816;

coeff[8][4][1][0] = 0;
coeff[8][4][1][1] = 0;
coeff[8][4][1][2] = 0;
coeff[8][4][1][3] = 0;
coeff[8][4][1][4] = 0;

coeff[8][5][0][0] = 40208700423117144062;
coeff[8][5][0][1] = 112242981279326174646352;
coeff[8][5][0][2] = 8753258717392862334713;
coeff[8][5][0][3] = 88809525183773972681612;
coeff[8][5][0][4] = 1521;

coeff[8][5][1][0] = 0;
coeff[8][5][1][1] = 0;
coeff[8][5][1][2] = 0;
coeff[8][5][1][3] = 0;
coeff[8][5][1][4] = 0;

coeff[9][0][0][0] = 1;
coeff[9][0][0][1] = 0;
coeff[9][0][0][2] = 0;
coeff[9][0][0][3] = 0;
coeff[9][0][0][4] = 0;

coeff[9][0][1][0] = 0;
coeff[9][0][1][1] = 0;
coeff[9][0][1][2] = 0;
coeff[9][0][1][3] = 0;
coeff[9][0][1][4] = 0;

coeff[9][1][0][0] = 28127254136958132407305;
coeff[9][1][0][1] = 21583712728585527334752;
coeff[9][1][0][2] = 51663912557182296229312;
coeff[9][1][0][3] = 17835324027044025627323;
coeff[9][1][0][4] = 1972711818993931957891;

coeff[9][1][1][0] = 122983429001812882742946;
coeff[9][1][1][1] = 80089057332764444586815;
coeff[9][1][1][2] = 105297218314245487796861;
coeff[9][1][1][3] = 92244217964995284597723;
coeff[9][1][1][4] = 5702367318890391334446;

coeff[9][2][0][0] = 0;
coeff[9][2][0][1] = 0;
coeff[9][2][0][2] = 0;
coeff[9][2][0][3] = 0;
coeff[9][2][0][4] = 0;

coeff[9][2][1][0] = 1;
coeff[9][2][1][1] = 0;
coeff[9][2][1][2] = 0;
coeff[9][2][1][3] = 0;
coeff[9][2][1][4] = 0;

coeff[9][3][0][0] = 28127254136958132407305;
coeff[9][3][0][1] = 21583712728585527334752;
coeff[9][3][0][2] = 51663912557182296229312;
coeff[9][3][0][3] = 17835324027044025627323;
coeff[9][3][0][4] = 1972711818993931957891;

coeff[9][3][1][0] = 28127254136958132407305;
coeff[9][3][1][1] = 21583712728585527334752;
coeff[9][3][1][2] = 51663912557182296229312;
coeff[9][3][1][3] = 17835324027044025627323;
coeff[9][3][1][4] = 1972711818993931957891;

coeff[9][4][0][0] = 151110683138771015150250;
coeff[9][4][0][1] = 101672770061349971921567;
coeff[9][4][0][2] = 5845403419599137187901;
coeff[9][4][0][3] = 110079541992039310225047;
coeff[9][4][0][4] = 7675079137884323292337;

coeff[9][4][1][0] = 0;
coeff[9][4][1][1] = 0;
coeff[9][4][1][2] = 0;
coeff[9][4][1][3] = 0;
coeff[9][4][1][4] = 0;

coeff[9][5][0][0] = 122983429001812882742946;
coeff[9][5][0][1] = 80089057332764444586815;
coeff[9][5][0][2] = 105297218314245487796861;
coeff[9][5][0][3] = 92244217964995284597723;
coeff[9][5][0][4] = 5702367318890391334446;

coeff[9][5][1][0] = 28127254136958132407305;
coeff[9][5][1][1] = 21583712728585527334752;
coeff[9][5][1][2] = 51663912557182296229312;
coeff[9][5][1][3] = 17835324027044025627323;
coeff[9][5][1][4] = 1972711818993931957891;

coeff[10][0][0][0] = 1;
coeff[10][0][0][1] = 0;
coeff[10][0][0][2] = 0;
coeff[10][0][0][3] = 0;
coeff[10][0][0][4] = 0;

coeff[10][0][1][0] = 0;
coeff[10][0][1][1] = 0;
coeff[10][0][1][2] = 0;
coeff[10][0][1][3] = 0;
coeff[10][0][1][4] = 0;

coeff[10][1][0][0] = 151070474438347898006189;
coeff[10][1][0][1] = 140545516233852444113487;
coeff[10][1][0][2] = 148207872154034921691459;
coeff[10][1][0][3] = 21270016808265337543434;
coeff[10][1][0][4] = 7675079137884323290816;

coeff[10][1][1][0] = 0;
coeff[10][1][1][1] = 0;
coeff[10][1][1][2] = 0;
coeff[10][1][1][3] = 0;
coeff[10][1][1][4] = 0;

coeff[10][2][0][0] = 151070474438347898006188;
coeff[10][2][0][1] = 140545516233852444113487;
coeff[10][2][0][2] = 148207872154034921691459;
coeff[10][2][0][3] = 21270016808265337543434;
coeff[10][2][0][4] = 7675079137884323290816;

coeff[10][2][1][0] = 0;
coeff[10][2][1][1] = 0;
coeff[10][2][1][2] = 0;
coeff[10][2][1][3] = 0;
coeff[10][2][1][4] = 0;

coeff[10][3][0][0] = 151110683138771015150250;
coeff[10][3][0][1] = 101672770061349971921567;
coeff[10][3][0][2] = 5845403419599137187901;
coeff[10][3][0][3] = 110079541992039310225047;
coeff[10][3][0][4] = 7675079137884323292337;

coeff[10][3][1][0] = 0;
coeff[10][3][1][1] = 0;
coeff[10][3][1][2] = 0;
coeff[10][3][1][3] = 0;
coeff[10][3][1][4] = 0;

coeff[10][4][0][0] = 40208700423117144062;
coeff[10][4][0][1] = 112242981279326174646352;
coeff[10][4][0][2] = 8753258717392862334713;
coeff[10][4][0][3] = 88809525183773972681612;
coeff[10][4][0][4] = 1521;

coeff[10][4][1][0] = 0;
coeff[10][4][1][1] = 0;
coeff[10][4][1][2] = 0;
coeff[10][4][1][3] = 0;
coeff[10][4][1][4] = 0;

coeff[10][5][0][0] = 40208700423117144063;
coeff[10][5][0][1] = 112242981279326174646352;
coeff[10][5][0][2] = 8753258717392862334713;
coeff[10][5][0][3] = 88809525183773972681612;
coeff[10][5][0][4] = 1521;

coeff[10][5][1][0] = 0;
coeff[10][5][1][1] = 0;
coeff[10][5][1][2] = 0;
coeff[10][5][1][3] = 0;
coeff[10][5][1][4] = 0;

coeff[11][0][0][0] = 1;
coeff[11][0][0][1] = 0;
coeff[11][0][0][2] = 0;
coeff[11][0][0][3] = 0;
coeff[11][0][0][4] = 0;

coeff[11][0][1][0] = 0;
coeff[11][0][1][1] = 0;
coeff[11][0][1][2] = 0;
coeff[11][0][1][3] = 0;
coeff[11][0][1][4] = 0;

coeff[11][1][0][0] = 33990078542129096261910;
coeff[11][1][0][1] = 98582514060874019318780;
coeff[11][1][0][2] = 42556846864733134377827;
coeff[11][1][0][3] = 75351964045704180477996;
coeff[11][1][0][4] = 1681895616428409458164;

coeff[11][1][1][0] = 117120604596641918888341;
coeff[11][1][1][1] = 3090256000475952602787;
coeff[11][1][1][2] = 114404284006694649648346;
coeff[11][1][1][3] = 34727577946335129747050;
coeff[11][1][1][4] = 5993183521455913834173;

coeff[11][2][0][0] = 0;
coeff[11][2][0][1] = 0;
coeff[11][2][0][2] = 0;
coeff[11][2][0][3] = 0;
coeff[11][2][0][4] = 0;

coeff[11][2][1][0] = 40208700423117144062;
coeff[11][2][1][1] = 112242981279326174646352;
coeff[11][2][1][2] = 8753258717392862334713;
coeff[11][2][1][3] = 88809525183773972681612;
coeff[11][2][1][4] = 1521;

coeff[11][3][0][0] = 122983429001812882742946;
coeff[11][3][0][1] = 80089057332764444586815;
coeff[11][3][0][2] = 105297218314245487796861;
coeff[11][3][0][3] = 92244217964995284597723;
coeff[11][3][0][4] = 5702367318890391334446;

coeff[11][3][1][0] = 122983429001812882742946;
coeff[11][3][1][1] = 80089057332764444586815;
coeff[11][3][1][2] = 105297218314245487796861;
coeff[11][3][1][3] = 92244217964995284597723;
coeff[11][3][1][4] = 5702367318890391334446;

coeff[11][4][0][0] = 40208700423117144063;
coeff[11][4][0][1] = 112242981279326174646352;
coeff[11][4][0][2] = 8753258717392862334713;
coeff[11][4][0][3] = 88809525183773972681612;
coeff[11][4][0][4] = 1521;

coeff[11][4][1][0] = 0;
coeff[11][4][1][1] = 0;
coeff[11][4][1][2] = 0;
coeff[11][4][1][3] = 0;
coeff[11][4][1][4] = 0;

coeff[11][5][0][0] = 5857780092113332166584;
coeff[11][5][0][1] = 27555843941809817067324;
coeff[11][5][0][2] = 147854065178978622174689;
coeff[11][5][0][3] = 16480454558870818237447;
coeff[11][5][0][4] = 7384262935318800792611;

coeff[11][5][1][0] = 145252903046657682983667;
coeff[11][5][1][1] = 74116926119540154854243;
coeff[11][5][1][2] = 9107065692449161851484;
coeff[11][5][1][3] = 93599087433168491987599;
coeff[11][5][1][4] = 290816202565522499726;
}
     
    return coeff;
}
