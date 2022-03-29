pragma circom 2.0.3;

// n bits per register
// num has k registers
// p has k registers
// k * n <= 500
// p is a prime
// if num == 0 mod p, returns 0
// else computes inv = num^{-1} mod p using extended euclidean algorithm
// https://brilliant.org/wiki/extended-euclidean-algorithm/
function find_Fp_inverse(n, k, num, p) {
    var amodp[2][50] = long_div2(n, k, 0, num, p); 
    var a[50];
    var b[50]; 
    var x[50];
    var y[50];
    var u[50];
    var v[50];

    var ret[50];

    for(var i=0; i<k; i++){
        a[i] = amodp[1][i];
        b[i] = p[i];
        x[i] = 0;
        y[i] = 0;
        u[i] = 0;
        v[i] = 0;
    }
    y[0] = 1;
    u[0] = 1;
    // euclidean algorithm takes log_phi( min(a, p) ) iterations, where phi is golden ratio
    // should be less than 1000 for our cases...
    for(var l=0; l<1000; l++){
        var ka = 0;
        for (var i = 0; i < k; i++) {
            if (a[i] != 0) {
                ka = i + 1;
            }
        }
        if (ka == 0) {
            for (var i = 0; i < k; i++) {
                ret[i] = x[i];
            }
            return ret;
        }

        var r[2][50] = long_div2(n, ka, k - ka, b, a); 
        var q[50]; 
        for(var i = 0; i < k - ka + 1; i++)
            q[i] = r[0][i];
        for(var i = k - ka + 1; i < k; i++)
            q[i] = 0;
        
        var newu[50] = long_sub_mod(n, k, x, prod_mod(n, k, u, q, p), p); 
        var newv[50] = long_sub_mod(n, k, y, prod_mod(n, k, v, q, p), p); 
        
        for(var i = 0; i < k; i++){
            b[i] = a[i];
            if( i < ka )
                a[i] = r[1][i];
            else
                a[i] = 0;
            x[i] = u[i];
            y[i] = v[i];
            u[i] = newu[i];
            v[i] = newv[i];
        }
    }
    // should never reach here (loop should always return before now)
    assert(0 == 1);
    return ret;
}

// a[k] registers can overflow - let's say in [0, B) 
//  assume actual value of a in  2^{k+m} 
// p[k] registers in [0, 2^n)
// out[2][k] solving
//      a = p * out[0] + out[1] with out[1] in [0,p) 
// out[0] has m registers in range [-2^n, 2^n)
// out[1] has k registers in range [0, 2^n)
function get_signed_Fp_carry_witness(n, k, m, a, p){
    var out[2][50];
    var a_short[51] = signed_long_to_short(n, k, a); 

    // let me make sure everything is in <= k+m registers
    /* commenting out to improve speed
    for(var j=k+m; j<50; j++)
        assert( a_short[j] == 0 );
    */

    if(a_short[50] == 0){
        out = long_div2(n, k, m, a_short, p);    
    }else{
        var a_pos[50];
        for(var i=0; i<k+m; i++) 
            a_pos[i] = -a_short[i];

        var X[2][50] = long_div2(n, k, m, a_pos, p);
        // what if X[1] is 0? 
        var Y_is_zero = 1;
        for(var i=0; i<k; i++){
            if(X[1][i] != 0)
                Y_is_zero = 0;
        }
        if( Y_is_zero == 1 ){
            out[1] = X[1];
        }else{
            out[1] = long_sub(n, k, p, X[1]); 
            
            X[0][0]++;
            if(X[0][0] >= (1<<n)){
                for(var i=0; i<m-1; i++){
                    var carry = X[0][i] \ (1<<n); 
                    X[0][i+1] += carry;
                    X[0][i] -= carry * (1<<n);
                }
                assert( X[0][m-1] < (1<<n) ); 
            }
        }
        for(var i=0; i<m; i++)
            out[0][i] = -X[0][i]; 
    }

    return out;
}


// Implements: 
//      calls get_signed_Fp_carry_witness twice
// a[2][k] registers can overflow
//  assume actual value of each a[i] < 2^{k+m} 
// p[k] registers in [0, 2^n)
// out[2][2][k] solving
//      a[0] = p * out[0][0] + out[0][1] with out[0][1] in [0,p) 
//      a[1] = p * out[1][0] + out[1][1] with out[1][1] in [0,p) 
// out[i][0] has m registers in range [-2^n, 2^n)
// out[i][1] has k registers in range [0, 2^n)
function get_signed_Fp2_carry_witness(n, k, m, a, p){
    var out[2][2][50];

    for(var i=0; i<2; i++)
        out[i] = get_signed_Fp_carry_witness(n, k, m, a[i], p);

    return out;
}

// helper function to precompute the product of two elements a, b in Fp2
// a[2][k], b[2][k] all registers in [0, 2^n) 
// (a0 + a1 u)*(b0 + b1 u) = (a0*b0 - a1*b1) + (a0*b1 + a1*b0)u 
// this is a direct computation - totally distinct from the combo of Fp2multiplyNoCarry and get_Fp2_carry_witness
function find_Fp2_product(n, k, a, b, p){
    var out[2][50];
    var ab[2][2][50]; 
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
    var out[2][50];
    out[0] = long_add_mod(n,k,a[0],b[0],p); 
    out[1] = long_add_mod(n,k,a[1],b[1],p);
    return out;
}

// helper function to precompute the difference of two elements a, b in Fp2
// a[2][k], b[2][k] all registers in [0, 2^n) 
// this is a direct computation
function find_Fp2_diff(n, k, a, b, p){
    var out[2][50];
    out[0] = long_sub_mod(n,k,a[0],b[0],p); 
    out[1] = long_sub_mod(n,k,a[1],b[1],p);
    return out;
}

// a[2][k] elt in Fp2 
// output multiplies by XI0 +u
function signed_Fp2_mult_w6(k, a, XI0){
    var out[2][50];
    for(var i=0; i<k; i++){
        out[0][i] = a[0][i]*XI0 - a[1][i];
        out[1][i] = a[0][i] + a[1][i]*XI0;
    }
    return out;
}


function find_Fp12_sum(n, k, a, b, p) {
    var out[6][2][50];
    for(var i=0; i<6; i++)
        out[i] = find_Fp2_sum(n, k, a[i], b[i], p);
    return out;
}

function find_Fp12_diff(n, k, a, b, p) {
    var out[6][2][50];
    for(var i=0; i<6; i++)
        out[i] = find_Fp2_diff(n, k, a[i], b[i], p);
    return out;
}

function find_Fp12_product(n, k, a, b, p) {
    var l = 6;
    var a0[l][50];
    var a1[l][50];
    var b0[l][50];
    var b1[l][50];
    var neg_b0[l][50];
    var neg_b1[l][50];
    var out[l][2][50];
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

    var real_init[20][50];
    var imag_init[20][50];
    var imag_init_neg[20][50];
    // var real[l][2][50];
    // var imag[l][2][50];
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
function find_Fp2_inverse(n, k, a, p) {
    var sq0[50] = prod(n, k, a[0], a[0]);
    var sq1[50] = prod(n, k, a[1], a[1]);
    var sq_sum[50] = long_add(n, 2*k, sq0, sq1);
    var sq_sum_div[2][50] = long_div2(n, k, k+1, sq_sum, p);
    // lambda = 1/(sq_sum)%p
    var lambda[50] = mod_inv(n, k, sq_sum_div[1], p);
    var out0[50] = prod(n, k, lambda, a[0]);
    var out0_div[2][50] = long_div(n, k, out0, p);
    var out[2][50];
    out[0] = out0_div[1];
    
    var out1_pre[50] = long_sub(n, k, p, a[1]);
    var out1[50] = prod(n, k, lambda, out1_pre);
    var out1_div[2][50] = long_div(n, k, out1, p);
    out[1] = out1_div[1];
    return out;
}


// a is 6 x 2 x k element of Fp^12
// compute inverse. first multiply by conjugate a + bw (a,b in Fp^6, w^6=1+u, u^2=-1)
// then reduce to inverting in Fp^6
function find_Fp12_inverse(n, k, p, a) {
    var A[6][2][50];
    var B[6][2][50];
    var Bw[6][2][50];
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
    var A2[6][2][50] = find_Fp12_product(n, k, A, A, p);
    var B2[6][2][50] = find_Fp12_product(n, k, B, B, p);
    var conj[6][2][50] = find_Fp12_diff(n, k, A, Bw, p);
    var w2[6][2][50];
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
    var B2w2[6][2][50] = find_Fp12_product(n, k, B2, w2, p);
    var conjProd[6][2][50] = find_Fp12_diff(n, k, A2, B2w2, p);
    var a0[2][50];
    var a1[2][50];
    var a2[2][50];
    for (var i = 0; i < 2; i ++) {
        for (var m = 0; m < k; m ++) {
            a0[i][m] = conjProd[0][i][m];
            a1[i][m] = conjProd[2][i][m];
            a2[i][m] = conjProd[4][i][m];
        }
    }
    var conjProdInv[6][2][50] = find_Fp6_inverse(n, k, p, a0, a1, a2);
    var out[6][2][50] = find_Fp12_product(n, k, conj, conjProdInv, p);
    return out;
}

// compute the inverse of a0 + a1v + a2v^2 in Fp6, where 
// v^3 = 1+u, u^2 = -1, a0 a1 a2 in Fp2 (2 x k)
// returns an element in standard Fp12 representation (6 x 2 x k)
function find_Fp6_inverse(n, k, p, a0, a1, a2) {
    var out[6][2][50];

    var a0_squared[2][50] = find_Fp2_product(n, k, a0, a0, p);
    var a1_squared[2][50] = find_Fp2_product(n, k, a1, a1, p);
    var a2_squared[2][50] = find_Fp2_product(n, k, a2, a2, p);
    var a0a1[2][50] = find_Fp2_product(n, k, a0, a1, p);
    var a0a2[2][50] = find_Fp2_product(n, k, a0, a2, p);
    var a1a2[2][50] = find_Fp2_product(n, k, a1, a2, p);
    var a0a1a2[2][50] = find_Fp2_product(n, k, a0a1, a2, p);

    var v3[2][50]; // v^3 = 1 + u
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            if (j == 0) {
                v3[i][j] = 1;
            } else {
                v3[i][j] = 0;
            }
        }
    }

    var three_v3[2][50]; // 3v^3 = 3 + 3u
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            if (j == 0) {
                three_v3[i][j] = 3;
            } else {
                three_v3[i][j] = 0;
            }
        }
    }

    var v6[2][50]; // v^6 = 2u
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            if (i == 1 && j == 0) {
                v6[i][j] = 2;
            } else {
                v6[i][j] = 0;
            }
        }
    }

    var v0_1[2][50] = find_Fp2_product(n, k, a1a2, v3, p);
    var v0_temp[2][50] = find_Fp2_diff(n, k, a0_squared, v0_1, p); // a0^2 - a1a2v^3
    var v1_1[2][50] = find_Fp2_product(n, k, a2_squared, v3, p);
    var v1_temp[2][50] = find_Fp2_diff(n, k, v1_1, a0a1, p); // v^3a2^2 - a0a1
    var v2_temp[2][50] = find_Fp2_diff(n, k, a1_squared, a0a2, p); // a1^2 - a0a2

    var a0_cubed[2][50] = find_Fp2_product(n, k, a0, a0_squared, p);
    var a1_cubed[2][50] = find_Fp2_product(n, k, a1, a1_squared, p);
    var a2_cubed[2][50] = find_Fp2_product(n, k, a2, a2_squared, p);
    var a13v3[2][50] = find_Fp2_product(n, k, a1_cubed, v3, p);
    var a23v6[2][50] = find_Fp2_product(n, k, a2_cubed, v6, p);
    var a0a1a23v3[2][50] = find_Fp2_product(n, k, a0a1a2, three_v3, p);

    var denom_1[2][50] = find_Fp2_sum(n, k, a0_cubed, a13v3, p);
    var denom_2[2][50] = find_Fp2_diff(n, k, a23v6, a0a1a23v3, p);
    var denom[2][50] = find_Fp2_sum(n, k, denom_1, denom_2, p); // a0^3 + a1^3v^3 + a2^3v^6 - 3a0a1a2v^3

    var denom_inv[2][50] = find_Fp2_inverse(n, k, denom, p);

    var v0_final[2][50] = find_Fp2_product(n, k, v0_temp, denom_inv, p);
    var v1_final[2][50] = find_Fp2_product(n, k, v1_temp, denom_inv, p);
    var v2_final[2][50] = find_Fp2_product(n, k, v2_temp, denom_inv, p);

    for (var i = 1; i < 6; i = i + 2) {
        for (var j = 0; j < 2; j ++) {
            for (var m = 0; m < 50; m ++) {
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

