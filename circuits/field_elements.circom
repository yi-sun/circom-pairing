pragma circom 2.0.2;

include "bigint.circom";

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
    signal output c[2][k];

    var LOGK = 3;
    assert(k<7);
    assert(2*n + 1 + LOGK<254);

    // c[0] computation
    // solve for X and Y such that a0*b0 + (p-a1)*b1 = p*X + Y with Y in [0,p) 
    // -a1*b1 = (p-a1)*b1 mod p
    var a0b0_var[100] = prod(n, k, a[0], b[0]);
    var a1_neg[100] = long_sub(n, k, p, a[1]); 
    var a1b1_neg[100] = prod(n, k, a1_neg, b[1]);
    var diff[100] = long_add(n, 2*k, a0b0_var, a1b1_neg); // 2*k+1 registers
    var X_Y[2][100] = long_div2(n, k, k+1, diff, p); 
    // X = X_Y[0] has k+2 registers, Y = X_Y[1] has k registers 
    // c[0] = Y
    for(var i=0; i<k; i++)
        c[0][i] <-- X_Y[1][i];
    component range_checks[k];
    for(var i=0; i<k; i++){
        range_checks[i] = Num2Bits(n);
        range_checks[i].in <== c[0][i]; 
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
        a0b0.a[i] <-- a[0][i];
        a0b0.b[i] <-- b[0][i];

        a1b1.a[i] <-- a[1][i];
        a1b1.b[i] <-- b[1][i];

        pb1.a[i] <-- p[i];
        pb1.b[i] <-- b[1][i];

        pX.a[i] <-- p[i];
        pX.b[i] <-- X_Y[0][i];
    }
    for(var i=k; i<k+2; i++){
        pX.a[i] <-- 0;
        pX.b[i] <-- X_Y[0][i];
    }

    component carry_check = CheckCarryToZero(n, 2*n+2+LOGK, 2*k+3); 
    for(var i=0; i<k; i++)
        carry_check.in[i] <-- a0b0.out[i] + pb1.out[i] - a1b1.out[i] - pX.out[i] - c[0][i]; 
    for(var i=k; i<2*k-1; i++)
        carry_check.in[i] <-- a0b0.out[i] + pb1.out[i] - a1b1.out[i] - pX.out[i]; 
    for(var i=2*k-1; i<2*k+3; i++)
        carry_check.in[i] <-- -pX.out[i];

    // now for c[1] computation
    // solve for Z and c[1] such that a0*b1 + a1*b0 = p*Z + c[1] with c[1] in [0,p) 
    var a0b1_var[100] = prod(n, k, a[0], b[1]);
    var a1b0_var[100] = prod(n, k, a[1], b[0]);
    var sum[100] = long_add(n, 2*k, a0b1_var, a1b0_var); // output 2*k+1 registers
    var sum_div[2][100] = long_div2(n, k, k+1, sum, p); 
    // Z = sum_div[0] has k+2 registers 
    for(var i=0; i<k; i++)
        c[1][i] <-- sum_div[1][i];
    component range_checks1[k];
    for(var i=0; i<k; i++){
        range_checks1[i] = Num2Bits(n);
        range_checks1[i].in <== c[1][i]; 
    }

    // constrain by Carry( a0 *' b1 +' a1 *' b0 -' p *' Z - c[1]) = 0 
    // each register is an overflow representation in the range (-(k+1)*2^{2n}-2^n, (k+1)*2^{2n + 1} )
    //                                          which is inside (-2^{2n+1+LOGK}, 2^{2n+1+LOGK})

    component a0b1 = BigMultShortLong(n, k); // 2*k-1 registers
    component a1b0 = BigMultShortLong(n, k);
    component pZ = BigMultShortLong(n, k+2); // 2*k+3 registers
    for(var i=0; i<k; i++){
        a0b1.a[i] <-- a[0][i];
        a0b1.b[i] <-- b[1][i];

        a1b0.a[i] <-- a[1][i];
        a1b0.b[i] <-- b[0][i];
        
        pZ.a[i] <-- p[i];
        pZ.b[i] <-- sum_div[0][i];
    }
    for(var i=k; i<k+2; i++){
        pZ.a[i] <-- 0;
        pZ.b[i] <-- sum_div[0][i];
    }
    
    component carry_check1 = CheckCarryToZero(n, 2*n+2+LOGK, 2*k+3);
    for(var i=0; i<k; i++)
        carry_check1.in[i] <-- a0b1.out[i] + a1b0.out[i] - pZ.out[i] - c[1][i]; 
    for(var i=k; i<2*k-1; i++)
        carry_check1.in[i] <-- a0b1.out[i] + a1b0.out[i] - pZ.out[i]; 
    for(var i=2*k-1; i<2*k+3; i++)
        carry_check1.in[i] <-- -pZ.out[i];
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

// squaring can be optimized to save 2 multiplication
// (a**2-b**2) = (a+b)(a-b) 
// (a+b u)**2 = (a+b)(a-b) + (a*b+a*b)u
template Fp2square(n, k){
    signal input in[2][k];
    signal input p[k];
    signal output out[2][k];
    
    component sum = BigAdd(n, k);
    for(var i=0; i<k; i++){
        sum.a[i] <== in[0][i];
        sum.b[i] <== in[1][i];
    }
    component diff = BigSubModP(n, k);
    for(var i=0; i<k; i++){
        diff.a[i] <== in[0][i];
        diff.b[i] <== in[1][i];
        diff.p[i] <== p[i];
    }
    component prod = BigMult(n, k+1);
    for(var i=0; i<k; i++){
        prod.a[i] <== sum.out[i];
        prod.b[i] <== diff.out[i];
    }
    prod.a[k] <== sum.out[k];
    prod.b[k] <== 0;

    component prod_mod = BigMod2(n, k, 2*k+2);
    for(var i=0; i<2*k+2; i++){
        prod_mod.a[i] <== prod.out[i];
        if(i<k){
            prod_mod.b[i] <== p[i];
        }
    }
    for(var i=0; i<k; i++){
        out[0][i] <== prod_mod.mod[i];
    }
    
    component ab = BigMult(n, k);
    for(var i=0; i<k; i++){
        ab.a[i] <== in[0][i];
        ab.b[i] <== in[1][i];
    }
    component two_ab = BigAdd(n, 2*k); 
    for(var i=0; i<2*k; i++){
        two_ab.a[i] <== ab.out[i];
        two_ab.b[i] <== ab.out[i];
    }
    component two_ab_mod = BigMod2(n, k, 2*k+1);
    for(var i=0; i<2*k+1; i++){
        two_ab_mod.a[i] <== two_ab.out[i];
        if(i < k){
            two_ab_mod.b[i] <== p[i];
        }
    }
    for(var i=0; i<k; i++){
        out[1][i] <== two_ab_mod.mod[i];
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
            in_out.c[i][j] === 1;
        else
            in_out.c[i][j] === 0;
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
    signal input a[6][2][k];
    signal input b[6][2][k];
    signal input p[k];
    signal output c[6][2][k];
    component ac = BigMultShortLong2D(n, k, 6);
    component bd = BigMultShortLong2D(n, k, 6);
    component ad = BigMultShortLong2D(n, k, 6);
    component bc = BigMultShortLong2D(n, k, 6);
    component neg_bc = BigMultShortLong2D(n, k, 6);
    component neg_ad = BigMultShortLong2D(n, k, 6);

    signal p_minus_d[6][k];
    component p_subtract_d[6];
    for (var i = 0; i < 6; i ++) {
        p_subtract_d[i] = BigSub(n, k);
        for (var m = 0; m < k; m ++) {
            p_subtract_d[i].a[m] <== p[m];
            p_subtract_d[i].b[m] <== b[i][1][m];
        }
        for (var m = 0; m < k; m ++) {
            p_minus_d[i][m] <== p_subtract_d[i].out[m];
        }
    }

    signal p_minus_c[6][k];
    component p_subtract_c[6];
    for (var i = 0; i < 6; i ++) {
        p_subtract_c[i] = BigSub(n, k);
        for (var m = 0; m < k; m ++) {
            p_subtract_c[i].a[m] <== p[m];
            p_subtract_c[i].b[m] <== b[i][0][m];
        }
        for (var m = 0; m < k; m ++) {
            p_minus_c[i][m] <== p_subtract_c[i].out[m];
        }
    }

    for (var i = 0; i < 6; i ++) {
        for (var j = 0; j < k; j ++) {
            ac.a[i][j] <== a[i][0][j];
            ac.b[i][j] <== b[i][0][j];
            bd.a[i][j] <== a[i][1][j];
            bd.b[i][j] <== p_minus_d[i][j];
            ad.a[i][j] <== a[i][0][j];
            ad.b[i][j] <== b[i][1][j];
            bc.a[i][j] <== a[i][1][j];
            bc.b[i][j] <== b[i][0][j];
            neg_ad.a[i][j] <== a[i][0][j];
            neg_ad.b[i][j] <== p_minus_d[i][j];
            neg_bc.a[i][j] <== a[i][1][j];
            neg_bc.b[i][j] <== p_minus_c[i][j];
        }
    }
    // ac + bd, ad + bc would both be 11 x (2k-1)
    signal long_result[6][2][2*k-1];
    for (var i = 0; i < 5; i ++) {
        for (var j = 0; j < 2*k-1; j ++) {
            long_result[i][0][j] <== ac.out[i][j] + bd.out[i][j] + ac.out[i+6][j] + bd.out[i+6][j] + neg_ad.out[i+6][j] + neg_bc.out[i+6][j];
            long_result[i][1][j] <== ad.out[i][j] + bc.out[i][j] + ac.out[i+6][j] + bd.out[i+6][j] + ad.out[i+6][j] + bc.out[i+6][j];
        }
    }
    for (var j = 0; j < 2*k-1; j ++) {
        long_result[5][0][j] <== ac.out[5][j] + bd.out[5][j];
        long_result[5][1][j] <== ad.out[5][j] + bc.out[5][j];
    }
    component longshorts[6][2];
    for (var i = 0; i < 6; i ++) {
        for (var j = 0; j < 2; j ++) {
            longshorts[i][j] = LongToShortNoEndCarry(n, 2*k-1);
            for (var m = 0; m < 2*k-1; m ++) {
                longshorts[i][j].in[m] <== long_result[i][j][m];
            }
        }
    }
    component bigmods[6][2];
    for (var i = 0; i < 6; i ++) {
        for (var j = 0; j < 2; j ++) {
            bigmods[i][j] = BigMod(n, k);
            for (var m = 0; m < 2*k; m ++) {
                bigmods[i][j].a[m] <== longshorts[i][j].out[m];
            }
            for (var m = 0; m < k; m ++) {
                bigmods[i][j].b[m] <== p[m];
            }
        }
    }
    for (var i = 0; i < 6; i ++) {
        for (var j = 0; j < k; j ++) {
            c[i][0][j] <== bigmods[i][0].mod[j];
            c[i][1][j] <== bigmods[i][1].mod[j];
        }
    }
}
