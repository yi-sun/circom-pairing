pragma circom 2.0.2;

include "bigint.circom";

// Extra circuits with different approaches to field operations. 
// Kept here for reference; less efficient or specialized than the ones in field_elements.circom


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