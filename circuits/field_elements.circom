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

    // c[0] computation
    component a0b0 = BigMultShortLong(n, k);
    for(var i=0; i<k; i++){
        a0b0.a[i] <== a[0][i];
        a0b0.b[i] <== b[0][i];
    }
    // -a1*b1 = (p-a1)*b1 mod p
    component a1_neg = BigSub(n, k);
    for(var i=0; i<k; i++){
        a1_neg.a[i] <== p[i];
        a1_neg.b[i] <== a[1][i];
    }
    component a1b1_neg = BigMultShortLong(n, k);
    for(var i=0; i<k; i++){
        a1b1_neg.a[i] <== a1_neg.out[i];
        a1b1_neg.b[i] <== b[1][i];
    }
    component diff = LongToShortNoEndCarry(n, 2*k + 1);
    for(var i=0; i<2*k-1; i++){
        diff.in[i] <== a0b0.out[i] + a1b1_neg.out[i];
    }
    diff.in[2*k-1] <== 0;
    diff.in[2*k] <== 0;
    
    component diff_mod = BigMod2(n, k, 2*k+1);
    for(var i=0; i<2*k+1; i++){
        diff_mod.a[i] <== diff.out[i];
    }
    for(var i=0; i<k; i++){
        diff_mod.b[i] <== p[i];
    }
    for(var i=0; i<k; i++){
        c[0][i] <== diff_mod.mod[i];
    }


    // now for c[1] computation
    component a0b1 = BigMultShortLong(n, k);
    for(var i=0; i<k; i++){
        a0b1.a[i] <== a[0][i];
        a0b1.b[i] <== b[1][i];
    }
    component a1b0 = BigMultShortLong(n, k);
    for(var i=0; i<k; i++){
        a1b0.a[i] <== a[1][i];
        a1b0.b[i] <== b[0][i];
    }
    component sum = LongToShortNoEndCarry(n, 2*k + 1);
    for(var i=0; i<2*k-1; i++){
        sum.in[i] <== a0b1.out[i] + a1b0.out[i];
    }
    sum.in[2*k-1] <== 0;
    sum.in[2*k] <== 0;
    
    component sum_mod = BigMod2(n, k, 2*k+1);
    for(var i=0; i<2*k+1; i++){
        sum_mod.a[i] <== sum.out[i];
    }
    for(var i=0; i<k; i++){
        sum_mod.b[i] <== p[i];
    }
    for(var i=0; i<k; i++){
        c[1][i] <== sum_mod.mod[i];
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
    var sq_sum[100];
    var carry[100];
    carry[0] = 0;
    for(var i=0; i<2*k; i++){
        var sumAndCarry[2] = SplitFn(sq0[i] + sq1[i] + carry[i], n, n);
        sq_sum[i] = sumAndCarry[0];
        carry[i+1] = sumAndCarry[1];
    }
    sq_sum[2*k] = carry[2*k];
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
// we first write a = A + Bi, b = C + Di and compute ab = (AC + B(p-D)) + (AD+BC)i
// with deg(w) = 10, deg(u) = 1 and then simplify the representation
// first using w^6 = u + 1 to get deg(w) = 5, deg (u) = 2
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
        }
    }
    // ac + bd, ad + bc would both be 11 x (2k-1)
    signal long_result[6][3][2*k-1];
    for (var i = 0; i < 5; i ++) {
        for (var j = 0; j < 2*k-1; j ++) {
            long_result[i][0][j] <== ac.out[i][j] + bd.out[i][j] + ac.out[i+6][j] + bd.out[i+6][j];
            long_result[i][1][j] <== ad.out[i][j] + bc.out[i][j] + ac.out[i+6][j] + bd.out[i+6][j] + ad.out[i+6][j] + bc.out[i+6][j];
            long_result[i][2][j] <== ad.out[i+6][j] + bc.out[i+6][j];
        }
    }
    for (var j = 0; j < 2*k-1; j ++) {
        long_result[5][0][j] <== ac.out[5][j] + bd.out[5][j];
        long_result[5][1][j] <== ad.out[5][j] + bc.out[5][j];
        long_result[5][2][j] <== 0;
    }
    component longshorts[6][3];
    for (var i = 0; i < 6; i ++) {
        for (var j = 0; j < 3; j ++) {
            longshorts[i][j] = LongToShortNoEndCarry(n, 2*k-1);
            for (var m = 0; m < 2*k-1; m ++) {
                longshorts[i][j].in[m] <== long_result[i][j][m];
            }
        }
    }
    component bigmods[6][3];
    for (var i = 0; i < 6; i ++) {
        for (var j = 0; j < 3; j ++) {
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
            c[i][1][j] <== bigmods[i][1].mod[j];
        }
    }
    component bigsubs[6];
    for (var i = 0; i < 6; i ++) {
        bigsubs[i] = BigSubModP(n, k);
        for (var m = 0; m < k; m ++) {
            bigsubs[i].a[m] <== bigmods[i][0].mod[m];
            bigsubs[i].b[m] <== bigmods[i][2].mod[m];
            bigsubs[i].p[m] <== p[m];
        }
        for (var m = 0; m < k; m ++) {
            c[i][0][m] <== bigsubs[i].out[m];
        }
    }
}
