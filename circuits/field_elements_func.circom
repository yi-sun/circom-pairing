pragma circom 2.0.2;

// a[4][k] registers can overflow - let's say in [0, B) 
//  assume actual value of each a[i] < 2^{k+m} 
// p[k] registers in [0, 2^n)
// out[2][2][k] solving
//      a[0] - a[1] = p * out[0][0] + out[0][1] with out[0][1] in [0,p) 
//      a[2] - a[3] = p * out[1][0] + out[1][1] with out[1][1] in [0,p) 
// out[i][0] has m registers in range [-2^n, 2^n)
// out[i][1] has k registers in range [0, 2^n)
function get_Fp2_carry_witness(n, k, m, a, p){
    var out[2][2][20];
    // solve for X and Y such that a0*b0 + (p-a1)*b1 = p*X + Y with Y in [0,p) 
    // -a1*b1 = (p-a1)*b1 mod p
    var a_short[4][20];
    for(var i=0; i<4; i++)
        a_short[i] = long_to_short(n, k, a[i]); 

    // let me make sure everything is in <= k+m registers
    for(var i=0; i<4; i++)
        for(var j=k+m; j<20; j++)
            assert( a_short[i][j] == 0 );

    var X[4][2][20];
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
// this is a direct computation - totally distinct from the combo of Fp2multiplyNoCarry and get_Fp2_carry_witness
function find_Fp2_product(n, k, a, b, p){
    var out[2][20];
    var ab[2][2][20]; 
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
    var out[2][20];
    out[0] = long_add_mod(n,k,a[0],b[0],p); 
    out[1] = long_add_mod(n,k,a[1],b[1],p);
    return out;
}

// helper function to precompute the difference of two elements a, b in Fp2
// a[2][k], b[2][k] all registers in [0, 2^n) 
// this is a direct computation
function find_Fp2_diff(n, k, a, b, p){
    var out[2][20];
    out[0] = long_sub_mod(n,k,a[0],b[0],p); 
    out[1] = long_sub_mod(n,k,a[1],b[1],p);
    return out;
}

// a[4][k] elt in Fp2 
// output multiplies by 1+u
function Fp2multc(k, a){
    var out[4][20];
    for(var i=0; i<k; i++){
        out[0][i] = a[0][i] + a[3][i];
        out[1][i] = a[1][i] + a[2][i];
        out[2][i] = a[0][i] + a[2][i];
        out[3][i] = a[1][i] + a[3][i];
    }
    return out;
}


function find_Fp12_sum(n, k, a, b, p) {
    var out[6][2][20];
    for(var i=0; i<6; i++)
        out[i] = find_Fp2_sum(n, k, a[i], b[i], p);
    return out;
}

function find_Fp12_diff(n, k, a, b, p) {
    var out[6][2][20];
    for(var i=0; i<6; i++)
        out[i] = find_Fp2_diff(n, k, a[i], b[i], p);
    return out;
}

function find_Fp12_product(n, k, a, b, p) {
    var l = 6;
    var a0[l][20];
    var a1[l][20];
    var b0[l][20];
    var b1[l][20];
    var neg_b0[l][20];
    var neg_b1[l][20];
    var out[l][2][20];
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

    var real_init[2*l-1][20];
    var imag_init[2*l-1][20];
    var imag_init_neg[2*l-1][20];
    // var real[l][2][20];
    // var imag[l][2][20];
    // each product will be 2l-1 x 2k
    var a0b0_var[20][20] = prod2D(n, k, l, a0, b0);
    var a1b1_neg[20][20] = prod2D(n, k, l, a1, neg_b1);
    var a0b1_var[20][20] = prod2D(n, k, l, a0, b1);
    var a1b0_var[20][20] = prod2D(n, k, l, a1, b0);
    var a0b1_neg[20][20] = prod2D(n, k, l, a0, neg_b1);
    var a1b0_neg[20][20] = prod2D(n, k, l, a1, neg_b0);
    for (var i = 0; i < 2*l - 1; i ++) { // compute initial rep (deg w = 10)
        real_init[i] = long_add(n, 2*k, a0b0_var[i], a1b1_neg[i]); // 2*k+1 registers each
        imag_init[i] = long_add(n, 2*k, a0b1_var[i], a1b0_var[i]);
        imag_init_neg[i] = long_add(n, 2*k, a0b1_neg[i], a1b0_neg[i]);
    }
    var real_carry[l][20];
    var imag_carry[l][20];
    var real_final[l][20];
    var imag_final[l][20];
    var zeros[20]; // to balance register sizes
    for (var i = 0; i < 20; i ++) {
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
    var XYreal_temp[l][2][20];
    var XYimag_temp[l][2][20];
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
    var sq0[20] = prod(n, k, a[0], a[0]);
    var sq1[20] = prod(n, k, a[1], a[1]);
    var sq_sum[20] = long_add(n, 2*k, sq0, sq1);
    var sq_sum_div[2][20] = long_div2(n, k, k+1, sq_sum, p);
    // lambda = 1/(sq_sum)%p
    var lambda[20] = mod_inv(n, k, sq_sum_div[1], p);
    var out0[20] = prod(n, k, lambda, a[0]);
    var out0_div[2][20] = long_div(n, k, out0, p);
    var out[2][20];
    out[0] = out0_div[1];
    
    var out1_pre[20] = long_sub(n, k, p, a[1]);
    var out1[20] = prod(n, k, lambda, out1_pre);
    var out1_div[2][20] = long_div(n, k, out1, p);
    out[1] = out1_div[1];
    return out;
}


// solve for
// a - b = out0 * p + out1
// assumes out0 has at most kX registers
function get_Fp12_carry_witness(n, k, kX, p, a, b) {
    var out[2][20];

    var a_short[20] = long_to_short(n, k, a);
    var b_short[20] = long_to_short(n, k, b);

    var X[2][2][20];
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


// a is 6 x 2 x k element of Fp^12
// compute inverse. first multiply by conjugate a + bw (a,b in Fp^6, w^6=1+u, u^2=-1)
// then reduce to inverting in Fp^6
function find_Fp12_inverse(n, k, p, a) {
    var A[6][2][20];
    var B[6][2][20];
    var Bw[6][2][20];
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
    var A2[6][2][20] = find_Fp12_product(n, k, A, A, p);
    var B2[6][2][20] = find_Fp12_product(n, k, B, B, p);
    var conj[6][2][20] = find_Fp12_diff(n, k, A, Bw, p);
    var w2[6][2][20];
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
    var B2w2[6][2][20] = find_Fp12_product(n, k, B2, w2, p);
    var conjProd[6][2][20] = find_Fp12_diff(n, k, A2, B2w2, p);
    var a0[2][20];
    var a1[2][20];
    var a2[2][20];
    for (var i = 0; i < 2; i ++) {
        for (var m = 0; m < k; m ++) {
            a0[i][m] = conjProd[0][i][m];
            a1[i][m] = conjProd[2][i][m];
            a2[i][m] = conjProd[4][i][m];
        }
    }
    var conjProdInv[6][2][20] = find_Fp6_inverse(n, k, p, a0, a1, a2);
    var out[6][2][20] = find_Fp12_product(n, k, conj, conjProdInv, p);
    return out;
}

// compute the inverse of a0 + a1v + a2v^2 in Fp6, where 
// v^3 = 1+u, u^2 = -1, a0 a1 a2 in Fp2 (2 x k)
// returns an element in standard Fp12 representation (6 x 2 x k)
function find_Fp6_inverse(n, k, p, a0, a1, a2) {
    var out[6][2][20];

    var a0_squared[2][20] = find_Fp2_product(n, k, a0, a0, p);
    var a1_squared[2][20] = find_Fp2_product(n, k, a1, a1, p);
    var a2_squared[2][20] = find_Fp2_product(n, k, a2, a2, p);
    var a0a1[2][20] = find_Fp2_product(n, k, a0, a1, p);
    var a0a2[2][20] = find_Fp2_product(n, k, a0, a2, p);
    var a1a2[2][20] = find_Fp2_product(n, k, a1, a2, p);
    var a0a1a2[2][20] = find_Fp2_product(n, k, a0a1, a2, p);

    var v3[2][20]; // v^3 = 1 + u
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            if (j == 0) {
                v3[i][j] = 1;
            } else {
                v3[i][j] = 0;
            }
        }
    }

    var three_v3[2][20]; // 3v^3 = 3 + 3u
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            if (j == 0) {
                three_v3[i][j] = 3;
            } else {
                three_v3[i][j] = 0;
            }
        }
    }

    var v6[2][20]; // v^6 = 2u
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            if (i == 1 && j == 0) {
                v6[i][j] = 2;
            } else {
                v6[i][j] = 0;
            }
        }
    }

    var v0_1[2][20] = find_Fp2_product(n, k, a1a2, v3, p);
    var v0_temp[2][20] = find_Fp2_diff(n, k, a0_squared, v0_1, p); // a0^2 - a1a2v^3
    var v1_1[2][20] = find_Fp2_product(n, k, a2_squared, v3, p);
    var v1_temp[2][20] = find_Fp2_diff(n, k, v1_1, a0a1, p); // v^3a2^2 - a0a1
    var v2_temp[2][20] = find_Fp2_diff(n, k, a1_squared, a0a2, p); // a1^2 - a0a2

    var a0_cubed[2][20] = find_Fp2_product(n, k, a0, a0_squared, p);
    var a1_cubed[2][20] = find_Fp2_product(n, k, a1, a1_squared, p);
    var a2_cubed[2][20] = find_Fp2_product(n, k, a2, a2_squared, p);
    var a13v3[2][20] = find_Fp2_product(n, k, a1_cubed, v3, p);
    var a23v6[2][20] = find_Fp2_product(n, k, a2_cubed, v6, p);
    var a0a1a23v3[2][20] = find_Fp2_product(n, k, a0a1a2, three_v3, p);

    var denom_1[2][20] = find_Fp2_sum(n, k, a0_cubed, a13v3, p);
    var denom_2[2][20] = find_Fp2_diff(n, k, a23v6, a0a1a23v3, p);
    var denom[2][20] = find_Fp2_sum(n, k, denom_1, denom_2, p); // a0^3 + a1^3v^3 + a2^3v^6 - 3a0a1a2v^3

    var denom_inv[2][20] = find_Fp2_inverse(n, k, denom, p);

    var v0_final[2][20] = find_Fp2_product(n, k, v0_temp, denom_inv, p);
    var v1_final[2][20] = find_Fp2_product(n, k, v1_temp, denom_inv, p);
    var v2_final[2][20] = find_Fp2_product(n, k, v2_temp, denom_inv, p);

    for (var i = 1; i < 6; i = i + 2) {
        for (var j = 0; j < 2; j ++) {
            for (var m = 0; m < 20; m ++) {
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

