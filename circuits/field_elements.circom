pragma circom 2.0.2;

include "bigint.circom";

// a, b are elements of Fp^l
// a[i] represents a[i][0] + a[i][1] * 2**n + ... + a[i][l-1] * 2**(n*(k-1))
// compute a+b in Fp^l
template FieldAdd2D(n, k, l) {
    signal input a[l][k];
    signal input b[l][k];
    signal input p[k];
    signal output c[l][k];

    component adders[l];
    for (var i = 0; i < l; i++) {
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


// take a polynomial expression a[0] + omega^1 a[1] + ... + omega^(2l-2) a[2l-2]
// reduce it to degree l-1 using omega^l + omega^(l-1) poly[l-1] + ... + poly[0] = 0
// WARNING: can produce incorrectly handled negative coefficients. only here for reference; do not use
template PolynomialReduce(l) {
    signal input a[2*l-1];
    signal input poly[l];
    signal output out[l];

    var residue[2*l-1];
    signal quotient[l-1];
    for (var i = 0; i < 2*l-1; i++) {
        residue[i] = a[i];
    }
    for (var i = l-2; i >= 0; i --) {
        for (var j = 0; j < l; j ++) {
            residue[i + j] = residue[i + j] - residue[i + l] * poly[j];
        }
        quotient[i] <-- residue[i+l];
        residue[i+l] = 0;
    }
    component mult = BigMultShortLong(1, l+1);
    for (var i = 0; i < l-1; i ++) {
        mult.a[i] <== quotient[i];
        mult.b[i] <== poly[i];
    }
    mult.a[l-1] <== 0;
    mult.a[l] <== 0;
    mult.b[l-1] <== poly[l-1];
    mult.b[l] <== 1;
    signal a_out[2*l-1];
    for (var i = 0; i < 2*l-1; i++) {
        a_out[i] <== mult.out[i];
    }
    for (var i = 0; i < l; i ++ ) {
        out[i] <-- residue[i];
    }
    for (var i = 0; i < l; i ++) {
        a[i] === a_out[i] + out[i];
    }
    for (var i = l; i < 2*l-1; i ++) {
        a[i] === a_out[i];
    }
}

template Fp2PolynomialReduce(n, k) {
    var l = 2;
    signal input a[2*l-1][k];
    var poly[2] = [1, 0]; // x^2 + 1 = 0
    signal output out[l][k];
    signal input p[k];

    for (var i = 0; i < k; i ++) {
        out[1][i] <== a[1][i];
    }
    component sub = BigSubModP(n, k);
    for (var i = 0; i < k; i ++) {
        sub.a[i] <== a[0][i];
        sub.b[i] <== a[2][i];
        sub.p[i] <== p[i];
    }
    for (var i = 0; i < k; i ++) {
        out[0][i] <== sub.out[i];
    }
}

// A similar circuit can do multiplication in different fields. 
// The only difference is that Fp2PolynomialReduce (which reduces quadratics by x^2+1) 
// must be replaced with a different circuit specialized to the minimal polynomial
template Fp2Multiply(n, k) {
    // l is always 2. poly is always [1, 0]
    var l = 2;
    signal input a[l][k];
    signal input b[l][k];
    signal input p[k];
    signal output c[l][k];

    component mult = BigMultShortLong2D(n, k, l);
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < k; j ++) {
            mult.a[i][j] <== a[i][j];
            mult.b[i][j] <== b[i][j];
        }
    } // out: 2l-1 x 2k-1 array of longs
    component longshorts[2*l-1];
    for (var i = 0; i < 2*l-1; i++) {
        longshorts[i] = LongToShortNoEndCarry(n, 2*k-1);
        for (var j = 0; j < 2*k-1; j ++) {
            longshorts[i].in[j] <== mult.out[i][j];
        }
    } // out: 2l-1 x 2k array of shorts
    component bigmods[2*l-1];
    for (var i = 0; i < 2*l-1; i ++) {
        bigmods[i] = BigMod(n, k);
        for (var j = 0; j < 2*k; j ++) {
            bigmods[i].a[j] <== longshorts[i].out[j];
        }
        for (var j = 0; j < k; j ++) {
            bigmods[i].b[j] <== p[j];
        }
    } // out: 2l-1 x k array of shorts
    component reduce = Fp2PolynomialReduce(n, k);
    for (var i = 0; i < 2*l-1; i ++) {
        for (var j = 0; j < k; j ++) {
            reduce.a[i][j] <== bigmods[i].mod[j];
        }
    } // out: l x k array of shorts
    for (var j = 0; j < k; j ++) {
        reduce.p[j] <== p[j];
    }
    for (var i = 0; i < l; i++) {
        for (var j = 0; j < k; j++) {
            c[i][j] <== reduce.out[i][j];
        }
    } // out: l x k array of shorts
}

// more optimized multiplication specialized to Fp^2 
// (a0 + a1 u)*(b0 + b1 u) = (a0*b0 - a1*b1) + (a0*b1 + a1*b0)u
template Fp2multiply(n, k){
    signal input a[2][k];
    signal input b[2][k];
    signal input p[k];
    signal output c[2][k];

    component a0b0 = BigMult(n, k);
    for(var i=0; i<k; i++){
        a0b0.a[i] <== a[0][i];
        a0b0.b[i] <== b[0][i];
    }
    component a1b1 = BigMult(n, k);
    for(var i=0; i<k; i++){
        a1b1.a[i] <== a[1][i];
        a1b1.b[i] <== b[1][i];
    }
    component gt = BigLessThan(n, 2*k);
    for(var i=0; i<2*k; i++){
        gt.a[i] <== a1b1.out[i];
        gt.b[i] <== a0b0.out[i];
    }

    // following is only good if a0b0 >= a1b1
    var a0b0_minus_a1b1[100] = long_sub(n, 2*k, a0b0.out, a1b1.out); 
    // following is only good if a1b1 >= a0b0
    var a1b1_minus_a0b0[100] = long_sub(n, 2*k, a1b1.out, a0b0.out); 

    // absolute value of a0*b0 - a1*b1:
    signal abs_diff[100];
    for(var i=0; i<2*k; i++){
        abs_diff[i] <-- gt.out * a0b0_minus_a1b1[i] + (1-gt.out)*a1b1_minus_a0b0[i];
    }
    component range_checks[2*k]; 
    for(var i=0; i<2*k; i++){
        range_checks[i] = Num2Bits(n);
        range_checks[i].in <== abs_diff[i];
    }
    // constraint is a0*b0 + (1-gt)*abs_diff = a1*b1 + gt*abs_diff 
    component a0b0_check = BigAdd(n, 2*k);
    component a1b1_check = BigAdd(n, 2*k);
    for(var i=0; i<2*k; i++){
        a0b0_check.a[i] <== a0b0.out[i];
        a0b0_check.b[i] <== (1-gt.out) * abs_diff[i];
        a1b1_check.a[i] <== a1b1.out[i];
        a1b1_check.b[i] <== gt.out * abs_diff[i];
    }
    for(var i=0; i<2*k+1; i++){
        a0b0_check.out[i] === a1b1_check.out[i];
    }

    // gives abs_diff % p 
    component abs_diff_mod = BigMod(n, k); 
    for(var i=0; i<2*k; i++){
        abs_diff_mod.a[i] <== abs_diff[i];
    }
    for(var i=0; i<k; i++){
        abs_diff_mod.b[i] <== p[i];
    }
    
    // -abs_diff % p
    component neg = BigSubModP(n,k);
    for(var i=0; i<k; i++){
        neg.a[i] <== 0;
        neg.b[i] <== abs_diff_mod.mod[i];
        neg.p[i] <== p[i];
    }
    
    for(var i=0; i<k; i++){
        c[0][i] <== gt.out*(abs_diff_mod.mod[i]-neg.out[i]) + neg.out[i];    
    }

    
    component a0b1 = BigMult(n, k);
    for(var i=0; i<k; i++){
        a0b1.a[i] <== a[0][i];
        a0b1.b[i] <== b[1][i];
    }
    component a1b0 = BigMult(n, k);
    for(var i=0; i<k; i++){
        a1b0.a[i] <== a[1][i];
        a1b0.b[i] <== b[0][i];
    }
    component sum = BigAdd(n, 2*k);
    for(var i=0; i<2*k; i++){
        sum.a[i] <== a0b1.out[i];
        sum.b[i] <== a1b0.out[i];
    }
    component ans = BigMod2(n, k, 2*k+1);
    for(var i=0; i<2*k+1; i++){
        ans.a[i] <== sum.out[i];
    }
    for(var i=0; i<k; i++){
        ans.b[i] <== p[i];
    }
    for(var i=0; i<k; i++){
        c[1][i] <== ans.mod[i];
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
    
}
