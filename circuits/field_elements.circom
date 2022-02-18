pragma circom 2.0.2;

include "bigint.circom";

// a, b are elements of Fp^l
// a[i] represents a[i][0] + a[i][1] * 2**n + ... + a[i][l-1] * 2**(n*(k-1))
// compute a+b in Fp^l
template FieldAdd(n, k, l) {
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
// WARNING: can produce negative coefficients. unsure how to handle yet
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

// does not work yet. correct flow, but 
// need to fix signed outputs from PolynomialReduce first. 
template FieldMultiply(n, k, l) {
    signal input a[l][k];
    signal input b[l][k];
    signal input poly[l];
    signal input p[k];
    signal output c[l][k];

    component mult = BigMultShortLong2D(n, k, l);
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < k; j ++) {
            mult.a[i][j] <== a[i][j];
            mult.b[i][j] <== b[i][j];
        }
    }

    component polynomials[2*k-1];
    for (var i = 0; i < 2*k-1; i ++) {
        polynomials[i] = PolynomialReduce(l);
        for (var j = 0; j < 2*l-1; j ++) {
            polynomials[i].a[j] <== mult.out[j][i];
        }
        for (var j = 0; j < l; j ++) {
            polynomials[i].poly[j] <== poly[j];
        }
    }
    component longshorts[l];
    for (var i = 0; i < l; i++) {
        longshorts[i] = LongToShortNoEndCarry(n, 2*k-1);
        for (var j = 0; j < 2*k-1; j ++) {
            longshorts[i].in[j] <== polynomials[j].out[i];
        }
    }
    component bigmods[l];
    for (var i = 0; i < l; i ++) {
        bigmods[i] = BigMod(n, k);
        for (var j = 0; j < 2*k; j ++) {
            bigmods[i].a[j] <== longshorts[i].out[j];
        }
        for (var j = 0; j < k; j ++) {
            bigmods[i].b[j] <== p[j];
        }
    }
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < k; j ++) {
            c[i][j] <== bigmods[i].mod[j];
        }
    }
}