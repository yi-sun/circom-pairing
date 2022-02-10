pragma circom 2.0.2;

include "../node_modules/circomlib/circuits/bitify.circom";

include "./bigint.circom";
include "./bigint_func.circom";

// requires a[0] != b[0]
//
// Implements:
// lamb = (b[1] - a[1]) / (b[0] - a[0]) % q
// out[0] = lamb ** 2 - a[0] - b[0] % q
// out[1] = lamb * (a[0] - out[0]) - a[1] % q
template EllipticCurveAddUnequal3Reg(n, q0, q1, q2) {
    var k = 3;
    signal input a[2][k];
    signal input b[2][k];

    signal output out[2][k];

    var q[100];
    q[0] = q0;
    q[1] = q1;
    q[2] = q2;
    for (var idx = 3; idx < 100; idx++) {
	q[idx] = 0;
    }
    
    // b[1] - a[1]
    component sub1 = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        sub1.a[i] <== b[1][i];
        sub1.b[i] <== a[1][i];
        sub1.p[i] <== q[i];
    }

    // b[0] - a[0]
    component sub0 = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        sub0.a[i] <== b[0][i];
        sub0.b[i] <== a[0][i];
        sub0.p[i] <== q[i];
    }

    signal lambda[k];
    var sub0inv[100] = mod_inv(n, k, sub0.out, q);
    var sub1_sub0inv[100] = prod(n, k, sub1.out, sub0inv);
    var lamb_arr[2][100] = long_div(n, k, sub1_sub0inv, q);
    for (var i = 0; i < k; i++) {
        lambda[i] <-- lamb_arr[1][i];
    }
    component range_checks[k];
    for (var i = 0; i < k; i++) {
        range_checks[i] = Num2Bits(n);
        range_checks[i].in <== lambda[i];
    }
    component lt = BigLessThan(n, k);
    for (var i = 0; i < k; i++) {
        lt.a[i] <== lambda[i];
        lt.b[i] <== q[i];
    }
    lt.out === 1;

    component lambda_check = BigMultModP(n, k);
    for (var i = 0; i < k; i++) {
        lambda_check.a[i] <== sub0.out[i];
        lambda_check.b[i] <== lambda[i];
        lambda_check.p[i] <== q[i];
    }
    for (var i = 0; i < k; i++) {
        lambda_check.out[i] === sub1.out[i];
    }

    component lambdasq = BigMultModP(n, k);
    for (var i = 0; i < k; i++) {
        lambdasq.a[i] <== lambda[i];
        lambdasq.b[i] <== lambda[i];
        lambdasq.p[i] <== q[i];
    }
    component out0_pre = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        out0_pre.a[i] <== lambdasq.out[i];
        out0_pre.b[i] <== a[0][i];
        out0_pre.p[i] <== q[i];
    }
    component out0 = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        out0.a[i] <== out0_pre.out[i];
        out0.b[i] <== b[0][i];
        out0.p[i] <== q[i];
    }
    for (var i = 0; i < k; i++) {
        out[0][i] <== out0.out[i];
    }

    component out1_0 = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        out1_0.a[i] <== a[0][i];
        out1_0.b[i] <== out[0][i];
        out1_0.p[i] <== q[i];
    }
    component out1_1 = BigMultModP(n, k);
    for (var i = 0; i < k; i++) {
        out1_1.a[i] <== lambda[i];
        out1_1.b[i] <== out1_0.out[i];
        out1_1.p[i] <== q[i];
    }
    component out1 = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        out1.a[i] <== out1_1.out[i];
        out1.b[i] <== a[1][i];
        out1.p[i] <== q[i];
    }
    for (var i = 0; i < k; i++) {
        out[1][i] <== out1.out[i];
    }    
}

// requires a[0] != b[0]
//
// Implements:
// lamb = (b[1] - a[1]) / (b[0] - a[0]) % q
// out[0] = lamb ** 2 - a[0] - b[0] % q
// out[1] = lamb * (a[0] - out[0]) - a[1] % q
template EllipticCurveAddUnequal4Reg(n, q0, q1, q2, q3) {
    var k = 4;
    signal input a[2][k];
    signal input b[2][k];

    signal output out[2][k];

    var q[100];
    q[0] = q0;
    q[1] = q1;
    q[2] = q2;
    q[3] = q3;
    for (var idx = 4; idx < 100; idx++) {
	q[idx] = 0;
    }
    
    // b[1] - a[1]
    component sub1 = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        sub1.a[i] <== b[1][i];
        sub1.b[i] <== a[1][i];
        sub1.p[i] <== q[i];
    }

    // b[0] - a[0]
    component sub0 = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        sub0.a[i] <== b[0][i];
        sub0.b[i] <== a[0][i];
        sub0.p[i] <== q[i];
    }

    signal lambda[k];
    var sub0inv[100] = mod_inv(n, k, sub0.out, q);
    var sub1_sub0inv[100] = prod(n, k, sub1.out, sub0inv);
    var lamb_arr[2][100] = long_div(n, k, sub1_sub0inv, q);
    for (var i = 0; i < k; i++) {
        lambda[i] <-- lamb_arr[1][i];
    }
    component range_checks[k];
    for (var i = 0; i < k; i++) {
        range_checks[i] = Num2Bits(n);
        range_checks[i].in <== lambda[i];
    }
    component lt = BigLessThan(n, k);
    for (var i = 0; i < k; i++) {
        lt.a[i] <== lambda[i];
        lt.b[i] <== q[i];
    }
    lt.out === 1;

    component lambda_check = BigMultModP(n, k);
    for (var i = 0; i < k; i++) {
        lambda_check.a[i] <== sub0.out[i];
        lambda_check.b[i] <== lambda[i];
        lambda_check.p[i] <== q[i];
    }
    for (var i = 0; i < k; i++) {
        lambda_check.out[i] === sub1.out[i];
    }

    component lambdasq = BigMultModP(n, k);
    for (var i = 0; i < k; i++) {
        lambdasq.a[i] <== lambda[i];
        lambdasq.b[i] <== lambda[i];
        lambdasq.p[i] <== q[i];
    }
    component out0_pre = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        out0_pre.a[i] <== lambdasq.out[i];
        out0_pre.b[i] <== a[0][i];
        out0_pre.p[i] <== q[i];
    }
    component out0 = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        out0.a[i] <== out0_pre.out[i];
        out0.b[i] <== b[0][i];
        out0.p[i] <== q[i];
    }
    for (var i = 0; i < k; i++) {
        out[0][i] <== out0.out[i];
    }

    component out1_0 = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        out1_0.a[i] <== a[0][i];
        out1_0.b[i] <== out[0][i];
        out1_0.p[i] <== q[i];
    }
    component out1_1 = BigMultModP(n, k);
    for (var i = 0; i < k; i++) {
        out1_1.a[i] <== lambda[i];
        out1_1.b[i] <== out1_0.out[i];
        out1_1.p[i] <== q[i];
    }
    component out1 = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        out1.a[i] <== out1_1.out[i];
        out1.b[i] <== a[1][i];
        out1.p[i] <== q[i];
    }
    for (var i = 0; i < k; i++) {
        out[1][i] <== out1.out[i];
    }    
}

// Elliptic curve is E : y**2 = x**3 + a x + b
// Note that for BLS12-381, a = 0, b = 4

// Implements:
// computing 2P on elliptic curve E for P = (in[0], in[1])
// formula from https://crypto.stanford.edu/pbc/notes/elliptic/explicit.html

// lamb =  (3 * in[0] ** 2 + a) / (2 * in[1]) % q
// out[0] = lamb ** 2 - 2 * in[0] % q
// out[1] = lamb * (in[0] - out[0]) - in[1] % q
template EllipticCurveDouble(n, k, a, q0, q1, q2, q3) {
    signal input in[2][k];

    signal output out[2][k];

    // assuming q < 2**(4n) 
    // represent q = q0 + q1 * 2**n + q2 * 2**(2n) + q3 * 2**(3n)
    // not sure how I feel about this convention...
    var q[100];
    q[0] = q0;
    q[1] = q1;
    q[2] = q2;
    q[3] = q3;
    for (var idx = 4; idx < 100; idx++) {
	    q[idx] = 0;
    }

    // assuming a is small 
    var long_a[100];
    long_a[0] = a;
    for (var i = 1; i < 100; i++) {
        long_a[i] = 0;   
    }

    component in0_sq = BigMultModP(n, k);
    for (var i = 0; i < k; i++) {
        in0_sq.a[i] <== in[0][i];
        in0_sq.b[i] <== in[0][i];
        in0_sq.p[i] <== q[i];
    }

    var long_2[100];
    var long_3[100];
    long_2[0] = 2;
    long_3[0] = 3;
    for (var i = 1; i < k; i++) {
        long_a[i] = 0;
        long_2[i] = 0;
        long_3[i] = 0;
    }
    var inv_2[100] = mod_inv(n, k, long_2, q);
    var long_3_div_2[100] = prod(n, k, long_3, inv_2);
    var long_3_div_2_mod_q[2][100] = long_div(n, k, long_3_div_2, q);

    // numerator = 3/2 * in[0]**2 + a
    // numer1 = 3/2 * in[0]**2
    component numer1 = BigMultModP(n, k);
    for (var i = 0; i < k; i++) {
        numer1.a[i] <== long_3_div_2_mod_q[1][i];
        numer1.b[i] <== in0_sq.out[i];
        numer1.p[i] <== q[i];
    }
    component numer = BigAddModP(n, k);
    for (var i = 0; i < k; i++) {
        numer.a[i] <== numer1.out[i];
        numer.b[i] <== long_a[i];
        numer.p[i] <== q[i];
    }

    signal lambda[k];
    var denom_inv[100] = mod_inv(n, k, in[1], q);
    var product[100] = prod(n, k, numer.out, denom_inv);
    var lamb_arr[2][100] = long_div(n, k, product, q);
    for (var i = 0; i < k; i++) {
        lambda[i] <-- lamb_arr[1][i];
    }
    component lt = BigLessThan(n, k);
    for (var i = 0; i < k; i++) {
        lt.a[i] <== lambda[i];
        lt.b[i] <== q[i];
    }
    lt.out === 1;

    component lambda_range_checks[k];
    component lambda_check = BigMultModP(n, k);
    for (var i = 0; i < k; i++) {
        lambda_range_checks[i] = Num2Bits(n);
        lambda_range_checks[i].in <== lambda[i];

        lambda_check.a[i] <== in[1][i];
        lambda_check.b[i] <== lambda[i];
        lambda_check.p[i] <== q[i];
    }
    for (var i = 0; i < k; i++) {
        lambda_check.out[i] === numer.out[i];
    }

    component lambdasq = BigMultModP(n, k);
    for (var i = 0; i < k; i++) {
        lambdasq.a[i] <== lambda[i];
        lambdasq.b[i] <== lambda[i];
        lambdasq.p[i] <== q[i];
    }
    component out0_pre = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        out0_pre.a[i] <== lambdasq.out[i];
        out0_pre.b[i] <== in[0][i];
        out0_pre.p[i] <== q[i];
    }
    // out0 = lambda**2 - 2*in[0]
    component out0 = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        out0.a[i] <== out0_pre.out[i];
        out0.b[i] <== in[0][i];
        out0.p[i] <== q[i];
    }
    for (var i = 0; i < k; i++) {
        out[0][i] <== out0.out[i];
    }

    component out1_0 = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        out1_0.a[i] <== in[0][i];
        out1_0.b[i] <== out[0][i];
        out1_0.p[i] <== q[i];
    }
    component out1_1 = BigMultModP(n, k);
    for (var i = 0; i < k; i++) {
        out1_1.a[i] <== lambda[i];
        out1_1.b[i] <== out1_0.out[i];
        out1_1.p[i] <== q[i];
    }
    component out1 = BigSubModP(n, k);
    for (var i = 0; i < k; i++) {
        out1.a[i] <== out1_1.out[i];
        out1.b[i] <== in[1][i];
        out1.p[i] <== q[i];
    }
    for (var i = 0; i < k; i++) {
        out[1][i] <== out1.out[i];
    }
}

