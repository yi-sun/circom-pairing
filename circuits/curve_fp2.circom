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