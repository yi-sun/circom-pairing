pragma circom 2.0.3;

include "../node_modules/circomlib/circuits/bitify.circom";

include "./bigint.circom";
include "./bigint_func.circom";
include "./fp.circom";
include "./fp2.circom";
include "./fp12.circom";
include "./bls12_381_func.circom";

// in[i] = (x_i, y_i) 
// Implements constraint: (y_1 + y_3) * (x_2 - x_1) - (y_2 - y_1)*(x_1 - x_3) = 0 mod p
// used to show (x1, y1), (x2, y2), (x3, -y3) are co-linear
template PointOnLineFp2(n, k, p) {
    signal input in[3][2][2][k]; 

    var LOGK = log_ceil(k);
    var LOGK2 = log_ceil(16*k*k);
    assert(3*n + LOGK2 < 251);

    // AKA check point on line 
    component left = BigMultShortLong2D(n, k, 2); // 3 x 2k-1 registers in [0, 8k*2^{2n+1})
    for(var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            left.a[i][j] <== in[0][1][i][j] + in[2][1][i][j];
            left.b[i][j] <== in[1][0][i][j] - in[0][0][i][j];
        }
    }

    component right = BigMultShortLong2D(n, k, 2); // 3 x 2k-1 registers in [0, 8k*2^{2n+1})
    for(var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            right.a[i][j] <== in[1][1][i][j] - in[0][1][i][j];
            right.b[i][j] <== in[0][0][i][j] - in[2][0][i][j];
        }
    }
    
    component diff_red[2]; 
    diff_red[0] = PrimeReduce(n, k, k-1, p, 3*n + 2*LOGK + 4);
    diff_red[1] = PrimeReduce(n, k, k-1, p, 3*n + 2*LOGK + 4);
    for(var i=0; i<2*k-1; i++) {
        diff_red[0].in[i] <== left.out[0][i] - left.out[2][i] - right.out[0][i] + right.out[2][i];
        diff_red[1].in[i] <== left.out[1][i] - right.out[1][i]; 
    }
    // diff_red has k registers in [0, 16*k^2*2^{3n} )
    component diff_mod[2];
    for (var j = 0; j < 2; j ++) {
        diff_mod[j] = SignedCheckCarryModToZero(n, k, 3*n + LOGK2, p);
        for (var i = 0; i < k; i ++) {
            diff_mod[j].in[i] <== diff_red[j].out[i];
        }
    }
}

// in = (x, y)
// Implements:
// x^3 + ax + b - y^2 = 0 mod p
// Assume: a, b in [0, 2^n) 

// test: component main { public [in] } = PointOnCurveFp2(2, 2, 0, 3, [1,1]);
// 
// /* INPUT = {
//     "in": [[[2,0],[3,0]],[[1,0],[2,0]]]
// } */
template PointOnCurveFp2(n, k, a, b, p){
    signal input in[2][2][k]; 

    var LOGK = log_ceil(k);
    var LOGK2 = log_ceil( (2*k-1)*k*k*2 );
    assert(4*n + 3 * LOGK + 4 < 251);

    // compute x^3, y^2 
    component x_sq = SignedFp2MultiplyNoCarryUnequal(n, k, k, 2*n+1+LOGK); // 2k-1 registers in [0, 2*k*2^{2n}) 
    component y_sq = SignedFp2MultiplyNoCarryUnequal(n, k, k, 2*n+1+LOGK); // 2k-1 registers in [0, 2*k*2^{2n}) 
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k ; j ++) {
            x_sq.a[i][j] <== in[0][i][j];
            x_sq.b[i][j] <== in[0][i][j];
            y_sq.a[i][j] <== in[1][i][j];
            y_sq.b[i][j] <== in[1][i][j];
        }
    }
    component x_cu = SignedFp2MultiplyNoCarryUnequal(n, 2*k-1, k, 3*n+2*LOGK+2); // 3k-2 registers in [0, 4*k^2 * 2^{3n}) 
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < 2*k-1; j ++) {
            x_cu.a[i][j] <== x_sq.out[i][j];
        }
        for (var j = 0; j < k; j ++) {
            x_cu.b[i][j] <== in[0][i][j];
        }
    }

    // x_cu + a x + b has 3k-2 registers < 2^{3n + 2LOGK + 2} 
    component cu_red[2];
    for (var j = 0; j < 2; j ++) {
        cu_red[j] = PrimeReduce(n, k, 2*k-2, p, 4*n + 3*LOGK + 4);
        for(var i=0; i<3*k-2; i++){
            if(i == 0) {
                if (j == 0)
                    cu_red[j].in[i] <== x_cu.out[j][i] + a * in[0][j][i] + b;
                else
                    cu_red[j].in[i] <== x_cu.out[j][i] + a * in[0][j][i];
            }
            else{
                if(i < k)
                    cu_red[j].in[i] <== x_cu.out[j][i] + a * in[0][j][i]; 
                else
                    cu_red[j].in[i] <== x_cu.out[j][i];
            }
        }
    }
    // cu_red has k registers < (2k-1)*2^{4n + 2LOGK + 1} < 2^{4n + 3LOGK + 4}

    component y_sq_red[2];
    for (var i = 0; i < 2; i ++) {
        y_sq_red[i] = PrimeReduce(n, k, k-1, p, 4*n + 3*LOGK + 4);
        for(var j=0; j<2*k-1; j++){
            y_sq_red[i].in[j] <== y_sq.out[i][j];
        }
    }

    component constraint[2];
    constraint[0] = SignedCheckCarryModToZero(n, k, 4*n + 3*LOGK2+4, p);
    constraint[1] = SignedCheckCarryModToZero(n, k, 4*n + 3*LOGK2+4, p);
    for(var i=0; i<k; i++){
        constraint[0].in[i] <== cu_red[0].out[i] - y_sq_red[0].out[i]; 
        constraint[1].in[i] <== cu_red[1].out[i] - y_sq_red[1].out[i];
    }
}

// component main { public [in] } = PointOnTangentFp2(2, 2, 2, [1, 1]);

// /* INPUT = {
//     "in": [[[[1,0],[1,0]],[[0,0],[2,0]]], [[[1,0],[2,0]],[[2,0],[0,1]]]]
// } */

// in[0] = (x_1, y_1), in[1] = (x_3, y_3) 
// Checks that the line between (x_1, y_1) and (x_3, -y_3) is equal to the tangent line to the elliptic curve at the point (x_1, y_1)
// Implements: 
// (y_1 + y_3) = lambda * (x_1 - x_3)
// where lambda = (3 x_1^2 + a)/(2 y_1) 
// Actual constraint is 2y_1 (y_1 + y_3) = (3 x_1^2 + a ) ( x_1 - x_3 )
template PointOnTangentFp2(n, k, a, p){
    signal input in[2][2][2][k];
    
    var LOGK = log_ceil(k);
    var LOGK3 = log_ceil((2*k-1)*7*k*k);
    assert(4*n + LOGK3 < 251);
    component x_sq = SignedFp2MultiplyNoCarryUnequal(n, k, k, 2*n+1+LOGK); // 2k-1 registers in [0, 2*k*2^{2n}) 
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k ; j ++) {
            x_sq.a[i][j] <== in[0][0][i][j];
            x_sq.b[i][j] <== in[0][0][i][j];
        }
    }
    component right = SignedFp2MultiplyNoCarryUnequal(n, 2*k-1, k, 3*n + 2*LOGK + 3); // 3k-2 registers < 2*3*k^2*2^{3n} 
    for(var i=0; i<2*k-1; i++){
        if(i == 0) {
            right.a[0][i] <== 3 * x_sq.out[0][i] + a; // registers in [0, 3*k*2^{2n} + 2^n )  
            right.a[1][i] <== 3 * x_sq.out[1][i];
        }
        else {
            right.a[0][i] <== 3 * x_sq.out[0][i];
            right.a[1][i] <== 3 * x_sq.out[1][i];
        }
    }
    for(var i=0; i<k; i++){
        right.b[0][i] <== in[0][0][0][i] - in[1][0][0][i]; 
        right.b[1][i] <== in[0][0][1][i] - in[1][0][1][i];
    }
    
    component left = SignedFp2MultiplyNoCarryUnequal(n, k, k, 2*n + 3 + LOGK); // 2k-1 registers in [0, k * 2^{2n+3})
    for(var i=0; i<k; i++){
        left.a[0][i] <== 2*in[0][1][0][i];
        left.a[1][i] <== 2*in[0][1][1][i];
        left.b[0][i] <== in[0][1][0][i] + in[1][1][0][i];
        left.b[1][i] <== in[0][1][1][i] + in[1][1][1][i];
    }
    
    // prime reduce right - left 
    component diff_red[2];
    for (var i = 0; i < 2; i ++) {
        diff_red[i] = PrimeReduce(n, k, 2*k-2, p, 4*n + LOGK3);
        for (var j = 0; j < 3*k-2; j ++) {
            if (j < 2*k-1) {
                diff_red[i].in[j] <== right.out[i][j] - left.out[i][j];
            }
            else {
                diff_red[i].in[j] <== right.out[i][j];
            }
        }
    }
    // inputs of diff_red has registers < 7*k^2*2^{3n} 
    // diff_red.out has registers < (2k-1)*7*k^2 * 2^{4n}
    component constraint[2];
    for (var i = 0; i < 2; i ++) {
        constraint[i] = SignedCheckCarryModToZero(n, k, 4*n + LOGK3, p);
        for (var j = 0; j < k; j ++) {
            constraint[i].in[j] <== diff_red[i].out[j];
        }
    }
}

// requires x_1 != x_2
// assume p is size k array, the prime that curve lives over 
//
// Implements:
//  Given a = (x_1, y_1) and b = (x_2, y_2), 
//      assume x_1 != x_2 and a != -b, 
//  Find a + b = (x_3, y_3)
// By solving:
//  x_1 + x_2 + x_3 - lambda^2 = 0 mod p
//  y_3 = lambda (x_1 - x_3) - y_1 mod p
//  where lambda = (y_2-y_1)/(x_2-x_1) is the slope of the line between (x_1, y_1) and (x_2, y_2)
// these equations are equivalent to:
//  (x_1 + x_2 + x_3)*(x_2 - x_1)^2 = (y_2 - y_1)^2 mod p
//  (y_1 + y_3)*(x_2 - x_1) = (y_2 - y_1)*(x_1 - x_3) mod p
template EllipticCurveAddUnequalFp2(n, k, p) { // changing q's to p's for my sanity
    signal input a[2][2][k];
    signal input b[2][2][k];

    signal output out[2][2][k];

    var LOGK = log_ceil(k);
    var LOGK3 = log_ceil( (12*k+1)*k*(2*k-1)); 
    assert(4*n + LOGK3 + 2< 251);

    // precompute lambda and x_3 and then y_3
    var dy[2][50] = find_Fp2_diff(n, k, b[1], a[1], p);
    var dx[2][50] = find_Fp2_diff(n, k, b[0], a[0], p); 
    var dx_inv[2][50] = find_Fp2_inverse(n, k, dx, p);
    var lambda[2][50] = find_Fp2_product(n, k, dy, dx_inv, p);
    var lambda_sq[2][50] = find_Fp2_product(n, k, lambda, lambda, p);
    // out[0] = x_3 = lamb^2 - a[0] - b[0] % p
    // out[1] = y_3 = lamb * (a[0] - x_3) - a[1] % p
    var x3[2][50] = find_Fp2_diff(n, k, find_Fp2_diff(n, k, lambda_sq, a[0], p), b[0], p);
    var y3[2][50] = find_Fp2_diff(n, k, find_Fp2_product(n, k, lambda, find_Fp2_diff(n, k, a[0], x3, p), p), a[1], p);

    for(var i = 0; i < k; i++){
        out[0][0][i] <-- x3[0][i];
        out[0][1][i] <-- x3[1][i];
        out[1][0][i] <-- y3[0][i];
        out[1][1][i] <-- y3[1][i];
    }
    
    // constrain x_3 by CUBIC (x_1 + x_2 + x_3) * (x_2 - x_1)^2 - (y_2 - y_1)^2 = 0 mod p
    
    component dx_sq = BigMultShortLong2D(n, k, 2); // 2k-1 registers < 4k*2^{2n} 
    component dy_sq = BigMultShortLong2D(n, k, 2); // 2k-1 registers < 4k*2^{2n}
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            dx_sq.a[i][j] <== b[0][i][j] - a[0][i][j];
            dx_sq.b[i][j] <== b[0][i][j] - a[0][i][j];
            dy_sq.a[i][j] <== b[1][i][j] - a[1][i][j];
            dy_sq.b[i][j] <== b[1][i][j] - a[1][i][j];
        }
    }

    // x_1 + x_2 + x_3 has registers in [0, 3*2^n) 
    component cubic = BigMultShortLong2DUnequal(n, k, 2*k-1, 2, 2); // 3k-2 x 3 registers < 24 * k^2 * 2^{3n} ) 
    for(var i=0; i<k; i++) {
        cubic.a[0][i] <== a[0][0][i] + b[0][0][i] + out[0][0][i]; 
        cubic.a[1][i] <== a[0][1][i] + b[0][1][i] + out[0][1][i];
    }
    for(var i=0; i<2*k-1; i++){
        cubic.b[0][i] <== dx_sq.out[0][i] - dx_sq.out[2][i];
        cubic.b[1][i] <== dx_sq.out[1][i];
    }

    component cubic_red[2];
    cubic_red[0] = PrimeReduce(n, k, 2*k-2, p, 4*n + LOGK3 + 2);
    cubic_red[1] = PrimeReduce(n, k, 2*k-2, p, 4*n + LOGK3 + 2);
    for(var i=0; i<2*k-1; i++) {
        // get i^2 parts too!
        cubic_red[0].in[i] <== cubic.out[0][i] - cubic.out[2][i] - dy_sq.out[0][i] + dy_sq.out[2][i]; // registers in < 12*k^2*2^{3n} + 4k*2^{2n} < (12k+1)k * 2^{3n} )
        cubic_red[1].in[i] <== cubic.out[1][i] - dy_sq.out[1][i]; // registers in < 12*k^2*2^{3n} + 4k*2^{2n} < (12k+1)k * 2^{3n} )
    }
    for(var i=2*k-1; i<3*k-2; i++) {
        cubic_red[0].in[i] <== cubic.out[0][i] - cubic.out[2][i]; 
        cubic_red[1].in[i] <== cubic.out[1][i];
    }
    // cubic_red has k registers < (2k-1) (12k+1)k * 2^{4n}
    
    component cubic_mod[2];
    cubic_mod[0] = SignedCheckCarryModToZero(n, k, 4*n + LOGK3 + 2, p);
    cubic_mod[1] = SignedCheckCarryModToZero(n, k, 4*n + LOGK3 + 2, p);
    for(var i=0; i<k; i++) {
        cubic_mod[0].in[i] <== cubic_red[0].out[i];
        cubic_mod[1].in[i] <== cubic_red[1].out[i];
    }
    // END OF CONSTRAINING x3
    
    // constrain y_3 by (y_1 + y_3) * (x_2 - x_1) = (y_2 - y_1)*(x_1 - x_3) mod p
    component y_constraint = PointOnLineFp2(n, k, p); // 2k-1 registers in [0, k*2^{2n+1})
    for(var i = 0; i < k; i++)for(var j=0; j<2; j++){
        for(var ind = 0; ind < 2; ind ++) {
            y_constraint.in[0][j][ind][i] <== a[j][ind][i];
            y_constraint.in[1][j][ind][i] <== b[j][ind][i];
            y_constraint.in[2][j][ind][i] <== out[j][ind][i];
        }
    }
    // END OF CONSTRAINING y3

    // check if out[][] has registers in [0, 2^n) and each out[i] is in [0, p)
    // re-using Fp2 code by considering (x_3, y_3) as a 2d-vector over Fp
    component range_check[2];
    range_check[0] = CheckValidFp2(n, k, p);
    range_check[1] = CheckValidFp2(n, k, p);
    for(var j=0; j<2; j++)for(var i=0; i<k; i++) {
        range_check[0].in[j][i] <== out[0][j][i];
        range_check[1].in[j][i] <== out[1][j][i];
    }
}

// component main { public [a,b] } = EllipticCurveAddUnequalFp2(2, 2, [1,1]);

// /* INPUT = {
//     "a": [[[1,0],[1,0]],[[2,0],[3,0]]],
//     "b": [[[2,0],[1,0]],[[1,0],[2,0]]]
// } */
// requires x_1 != x_2
// assume p is size k array, the prime that curve lives over 
//
// Implements:
//  Given a = (x_1, y_1) and b = (x_2, y_2), 
//      assume x_1 != x_2 and a != -b, 
//  Find a + b = (x_3, y_3)
// By solving:
//  x_1 + x_2 + x_3 - lambda^2 = 0 mod p
//  y_3 = lambda (x_1 - x_3) - y_1 mod p
//  where lambda = (y_2-y_1)/(x_2-x_1) is the slope of the line between (x_1, y_1) and (x_2, y_2)
// these equations are equivalent to:
//  (x_1 + x_2 + x_3)*(x_2 - x_1)^2 = (y_2 - y_1)^2 mod p
//  (y_1 + y_3)*(x_2 - x_1) = (y_2 - y_1)*(x_1 - x_3) mod p
template EllipticCurveAddUnequalFp2(n, k, p) { // changing q's to p's for my sanity
    signal input a[2][2][k];
    signal input b[2][2][k];

    signal output out[2][2][k];

    var LOGK = log_ceil(k);
    var LOGK3 = log_ceil( (12*k+1)*k*(2*k-1)); 
    assert(4*n + LOGK3 + 2< 251);

    // precompute lambda and x_3 and then y_3
    var dy[2][50] = find_Fp2_diff(n, k, b[1], a[1], p);
    var dx[2][50] = find_Fp2_diff(n, k, b[0], a[0], p); 
    var dx_inv[2][50] = find_Fp2_inverse(n, k, dx, p);
    var lambda[2][50] = find_Fp2_product(n, k, dy, dx_inv, p);
    var lambda_sq[2][50] = find_Fp2_product(n, k, lambda, lambda, p);
    // out[0] = x_3 = lamb^2 - a[0] - b[0] % p
    // out[1] = y_3 = lamb * (a[0] - x_3) - a[1] % p
    var x3[2][50] = find_Fp2_diff(n, k, find_Fp2_diff(n, k, lambda_sq, a[0], p), b[0], p);
    var y3[2][50] = find_Fp2_diff(n, k, find_Fp2_product(n, k, lambda, find_Fp2_diff(n, k, a[0], x3, p), p), a[1], p);

    for(var i = 0; i < k; i++){
        out[0][0][i] <-- x3[0][i];
        out[0][1][i] <-- x3[1][i];
        out[1][0][i] <-- y3[0][i];
        out[1][1][i] <-- y3[1][i];
    }
    
    // constrain x_3 by CUBIC (x_1 + x_2 + x_3) * (x_2 - x_1)^2 - (y_2 - y_1)^2 = 0 mod p
    
    component dx_sq = BigMultShortLong2D(n, k, 2); // 2k-1 registers < 4k*2^{2n} 
    component dy_sq = BigMultShortLong2D(n, k, 2); // 2k-1 registers < 4k*2^{2n}
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            dx_sq.a[i][j] <== b[0][i][j] - a[0][i][j];
            dx_sq.b[i][j] <== b[0][i][j] - a[0][i][j];
            dy_sq.a[i][j] <== b[1][i][j] - a[1][i][j];
            dy_sq.b[i][j] <== b[1][i][j] - a[1][i][j];
        }
    }

    // x_1 + x_2 + x_3 has registers in [0, 3*2^n) 
    component cubic = BigMultShortLong2DUnequal(n, k, 2*k-1, 2, 2); // 3k-2 x 3 registers < 24 * k^2 * 2^{3n} ) 
    for(var i=0; i<k; i++) {
        cubic.a[0][i] <== a[0][0][i] + b[0][0][i] + out[0][0][i]; 
        cubic.a[1][i] <== a[0][1][i] + b[0][1][i] + out[0][1][i];
    }
    for(var i=0; i<2*k-1; i++){
        cubic.b[0][i] <== dx_sq.out[0][i] - dx_sq.out[2][i];
        cubic.b[1][i] <== dx_sq.out[1][i];
    }

    component cubic_red[2];
    cubic_red[0] = PrimeReduce(n, k, 2*k-2, p, 4*n + LOGK3 + 2);
    cubic_red[1] = PrimeReduce(n, k, 2*k-2, p, 4*n + LOGK3 + 2);
    for(var i=0; i<2*k-1; i++) {
        // get i^2 parts too!
        cubic_red[0].in[i] <== cubic.out[0][i] - cubic.out[2][i] - dy_sq.out[0][i] + dy_sq.out[2][i]; // registers in < 12*k^2*2^{3n} + 4k*2^{2n} < (12k+1)k * 2^{3n} )
        cubic_red[1].in[i] <== cubic.out[1][i] - dy_sq.out[1][i]; // registers in < 12*k^2*2^{3n} + 4k*2^{2n} < (12k+1)k * 2^{3n} )
    }
    for(var i=2*k-1; i<3*k-2; i++) {
        cubic_red[0].in[i] <== cubic.out[0][i] - cubic.out[2][i]; 
        cubic_red[1].in[i] <== cubic.out[1][i];
    }
    // cubic_red has k registers < (2k-1) (12k+1)k * 2^{4n}
    
    component cubic_mod[2];
    cubic_mod[0] = SignedCheckCarryModToZero(n, k, 4*n + LOGK3 + 2, p);
    cubic_mod[1] = SignedCheckCarryModToZero(n, k, 4*n + LOGK3 + 2, p);
    for(var i=0; i<k; i++) {
        cubic_mod[0].in[i] <== cubic_red[0].out[i];
        cubic_mod[1].in[i] <== cubic_red[1].out[i];
    }
    // END OF CONSTRAINING x3
    
    // constrain y_3 by (y_1 + y_3) * (x_2 - x_1) = (y_2 - y_1)*(x_1 - x_3) mod p
    component y_constraint = PointOnLineFp2(n, k, p); // 2k-1 registers in [0, k*2^{2n+1})
    for(var i = 0; i < k; i++)for(var j=0; j<2; j++){
        for(var ind = 0; ind < 2; ind ++) {
            y_constraint.in[0][j][ind][i] <== a[j][ind][i];
            y_constraint.in[1][j][ind][i] <== b[j][ind][i];
            y_constraint.in[2][j][ind][i] <== out[j][ind][i];
        }
    }
    // END OF CONSTRAINING y3

    // check if out[][] has registers in [0, 2^n) and each out[i] is in [0, p)
    // re-using Fp2 code by considering (x_3, y_3) as a 2d-vector over Fp
    component range_check[2];
    range_check[0] = CheckValidFp2(n, k, p);
    range_check[1] = CheckValidFp2(n, k, p);
    for(var j=0; j<2; j++)for(var i=0; i<k; i++) {
        range_check[0].in[j][i] <== out[0][j][i];
        range_check[1].in[j][i] <== out[1][j][i];
    }
}

// component main { public [in] } = EllipticCurveDoubleFp2(2, 2, 0, 2, [3,1]);

// /* INPUT = {
//     "in": [[[2,0],[1,0]],[[1,0],[2,0]]]
// } */

// Elliptic curve is E : y**2 = x**3 + ax + b
// assuming a < 2^n for now
// Note that for BLS12-381, a = 0, b = 4

// Implements:
// computing 2P on elliptic curve E for P = (x_1, y_1)
// formula from https://crypto.stanford.edu/pbc/notes/elliptic/explicit.html
// x_1 = in[0], y_1 = in[1]
// assume y_1 != 0 (otherwise 2P = O)

// lamb =  (3x_1^2 + a) / (2 y_1) % p
// x_3 = out[0] = lambda^2 - 2 x_1 % p
// y_3 = out[1] = lambda (x_1 - x_3) - y_1 % p

// We precompute (x_3, y_3) and then constrain by showing that:
// * (x_3, y_3) is a valid point on the curve 
// * the slope (y_3 - y_1)/(x_3 - x_1) equals 
// * x_1 != x_3 
template EllipticCurveDoubleFp2(n, k, a, b, p) {
    signal input in[2][2][k];
    signal output out[2][2][k];

    var long_a[2][k];
    var long_3[2][k];
    long_a[0][0] = a;
    long_3[0][0] = 3;
    long_a[1][0] = 0;
    long_3[1][0] = 0;
    for (var i = 1; i < k; i++) {
        long_a[0][i] = 0;
        long_3[0][i] = 0;
        long_a[1][i] = 0;
        long_3[1][i] = 0;
    }

    // precompute lambda 
    var lamb_num[2][50] = find_Fp2_sum(n, k, long_a, find_Fp2_product(n, k, long_3, find_Fp2_product(n, k, in[0], in[0], p), p), p);
    var lamb_denom[2][50] = find_Fp2_sum(n, k, in[1], in[1], p);
    var lamb[2][50] = find_Fp2_product(n, k, lamb_num, find_Fp2_inverse(n, k, lamb_denom, p), p);

    // precompute x_3, y_3
    var x3[2][50] = find_Fp2_diff(n, k, find_Fp2_product(n, k, lamb, lamb, p), find_Fp2_sum(n, k, in[0], in[0], p), p);
    var y3[2][50] = find_Fp2_diff(n, k, find_Fp2_product(n, k, lamb, find_Fp2_diff(n, k, in[0], x3, p), p), in[1], p);
    
    for(var i=0; i<k; i++){
        out[0][0][i] <-- x3[0][i];
        out[0][1][i] <-- x3[1][i];
        out[1][0][i] <-- y3[0][i];
        out[1][1][i] <-- y3[1][i];
    }
    // check if out[][] has registers in [0, 2^n) and each out[i] is in [0, p)
    // re-using Fp2 code by considering (x_3, y_3) as a 2d-vector over Fp
    component range_check[2];
    range_check[0] = CheckValidFp2(n, k, p);
    range_check[1] = CheckValidFp2(n, k, p);
    for(var j=0; j<2; j++)for(var i=0; i<k; i++) {
        range_check[0].in[j][i] <== out[0][j][i];
        range_check[1].in[j][i] <== out[1][j][i];
    }

    component point_on_tangent = PointOnTangentFp2(n, k, a, p);
    for(var j=0; j<2; j++)for(var i=0; i<k; i++){
        point_on_tangent.in[0][j][0][i] <== in[j][0][i];
        point_on_tangent.in[0][j][1][i] <== in[j][1][i];
        point_on_tangent.in[1][j][0][i] <== out[j][0][i];
        point_on_tangent.in[1][j][1][i] <== out[j][1][i];
    }
    
    component point_on_curve = PointOnCurveFp2(n, k, a, b, p);
    for(var j=0; j<2; j++)for(var i=0; i<k; i++) {
        point_on_curve.in[j][0][i] <== out[j][0][i];
        point_on_curve.in[j][1][i] <== out[j][1][i];
    }
    
    component x3_eq_x1 = IsArrayEqual(2*k);
    for(var i = 0; i < k; i++){
        x3_eq_x1.in[0][i] <== out[0][0][i];
        x3_eq_x1.in[1][i] <== in[0][0][i];
        x3_eq_x1.in[0][i+k] <== out[0][1][i];
        x3_eq_x1.in[1][i+k] <== in[0][1][i];
    }
    x3_eq_x1.out === 0;
}