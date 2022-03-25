pragma circom 2.0.3;

include "../node_modules/circomlib/circuits/bitify.circom";

include "./bigint.circom";
include "./bigint_func.circom";
include "./fp.circom";
include "./fp2.circom";
include "./fp12.circom";

// taken from https://zkrepl.dev/?gist=1e0a28ec3cc4967dc0994e30d316a8af
template IsArrayEqual(k){
    signal input in[2][k];
    signal output out;
    component isEqual[k+1];
    var sum = 0;
    for(var i = 0; i < k; i++){
        isEqual[i] = IsEqual();
        isEqual[i].in[0] <== in[0][i];
        isEqual[i].in[1] <== in[1][i];
        sum = sum + isEqual[i].out;
    }

    isEqual[k] = IsEqual();
    isEqual[k].in[0] <== sum;
    isEqual[k].in[1] <== k;
    out <== isEqual[k].out;
}

// in[i] = (x_i, y_i) 
// Implements constraint: (y_1 + y_3) * (x_2 - x_1) - (y_2 - y_1)*(x_1 - x_3) = 0 mod p
// used to show (x1, y1), (x2, y2), (x3, -y3) are co-linear
template PointOnLine(n, k, p) {
    signal input in[3][2][k]; 

    var LOGK = log_ceil(k);
    assert(3*n + 2*LOGK + 2 < 253);

    // AKA check point on line 
    component left = BigMultNoCarry(n, k, 2*n + 2*LOGK + 2); // 2k-1 registers in [0, k*2^{2n+1})
    for(var i = 0; i < k; i++){
        left.a[0][i] <== in[0][1][i] + in[2][1][i];
        left.a[1][i] <== 0;
        left.b[0][i] <== in[1][0][i];
        left.b[1][i] <== in[0][0][i]; 
    }

    component right = BigMultNoCarry(n, k, 2*n + 2*LOGK + 2); // 2k-1 registers in [0, k*2^{2n+1})
    for(var i = 0; i < k; i++){
        right.a[0][i] <== in[1][1][i];
        right.a[1][i] <== in[0][1][i];
        right.b[0][i] <== in[0][0][i];
        right.b[1][i] <== in[2][0][i];
    }
    
    component diff_red[2]; 
    for(var j=0; j<2; j++){
        diff_red[j] = PrimeReduce(n, k, k-1, p, 3*n + 2*LOGK + 2);
        for(var i=0; i<2*k-1; i++)
            diff_red[j].in[i] <== left.out[j][i] + right.out[j^1][i];  
    }
    // diff_red has 2 x k registers in [0, k^2*2^{3n+2} < 2^{3n + 2LOGK + 2} )
    component diff_mod = CheckCarryModToZero(n, k, 3*n + 2*LOGK + 2, p);
    for(var j=0; j<2; j++)for(var i=0; i<k; i++)
        diff_mod.in[j][i] <== diff_red[j].out[i]; 
}

// in = (x, y)
// Implements:
// x^3 + ax + b - y^2 = 0 mod p
// Assume: a, b in [0, 2^n) 
template PointOnCurve(n, k, a, b, p){
    signal input in[2][k]; 

    var LOGK = log_ceil(k);
    assert(4*n + 3*LOGK + 2 < 253);

    // compute x^3, y^2 
    component x_sq = BigMultShortLong(n, k, 2*n + LOGK); // 2k-1 registers in [0, k*2^{2n}) 
    component y_sq = BigMultShortLong(n, k, 2*n + LOGK); // 2k-1 registers in [0, k*2^{2n}) 
    for(var i=0; i<k; i++){
        x_sq.a[i] <== in[0][i];
        x_sq.b[i] <== in[0][i];

        y_sq.a[i] <== in[1][i];
        y_sq.b[i] <== in[1][i];
    }
    component x_cu = BigMultShortLongUnequal(n, 2*k-1, k, 3*n + 2*LOGK + 2); // 3k-2 registers in [0, k^2 * 2^{3n}) 
    for(var i=0; i<2*k-1; i++)
        x_cu.a[i] <== x_sq.out[i];
    for(var i=0; i<k; i++)
        x_cu.b[i] <== in[0][i];

    // x_cu + a x + b has 3k-2 registers < k^2 * 2^{3n} + 2^{2n} + 2^n < 2^{3n + 2LOGK + 1} 
    component cu_red = PrimeReduce(n, k, 2*k-2, p, 4*n + 3*LOGK + 2);
    for(var i=0; i<3*k-2; i++){
        if(i == 0)
            cu_red.in[i] <== x_cu.out[i] + a * in[0][i] + b; 
        else{
            if(i < k)
                cu_red.in[i] <== x_cu.out[i] + a * in[0][i]; 
            else
                cu_red.in[i] <== x_cu.out[i];
        }
    }
    // cu_red has k registers < (2k-1)*2^{4n + 2LOGK + 1} < 2^{4n + 3LOGK + 2}

    component y_sq_red = PrimeReduce(n, k, k-1, p, 4*n + 3*LOGK + 2);
    for(var i=0; i<2*k-1; i++)
        y_sq_red.in[i] <== y_sq.out[i]; 

    component constraint = CheckCarryModToZero(n, k, 4*n + 3*LOGK + 2, p);
    for(var i=0; i<k; i++){
        constraint.in[0][i] <== cu_red.out[i];
        constraint.in[1][i] <== y_sq_red.out[i]; 
    }
}

// in[0] = (x_1, y_1), in[1] = (x_3, y_3) 
// Checks that the line between (x_1, y_1) and (x_3, -y_3) is equal to the tangent line to the elliptic curve at the point (x_1, y_1)
// Implements: 
// (y_1 + y_3) = lambda * (x_1 - x_3)
// where lambda = (3 x_1^2 + a)/(2 y_1) 
// Actual constraint is 2y_1 (y_1 + y_3) = (3 x_1^2 + a ) ( x_1 - x_3 )
template PointOnTangent(n, k, a, p){
    signal input in[2][2][k];
    
    var LOGK = log_ceil(k);
    assert(4*n + 3*LOGK + 3 < 253);
    component x_sq = BigMultShortLong(n, k, 2*n + LOGK); // 2k-1 registers in [0, k*2^{2n}) 
    for(var i=0; i<k; i++){
        x_sq.a[i] <== in[0][0][i];
        x_sq.b[i] <== in[0][0][i];
    }
    component right = BigMultNoCarryUnequal(n, 2*k-1, k, 3*n + 2*LOGK + 2); // 3k-2 registers in [0, 3*k^2*2^{3n} + k*2^{2n} ) 
    for(var i=0; i<2*k-1; i++){
        if(i == 0)
            right.a[0][i] <== 3 * x_sq.out[i] + a; // registers in [0, 3*k*2^{2n} + 2^n )  
        else
            right.a[0][i] <== 3 * x_sq.out[i]; 
        right.a[1][i] <== 0;
    }
    for(var i=0; i<k; i++){
        right.b[0][i] <== in[0][0][i];
        right.b[1][i] <== in[1][0][i]; 
    }
    
    component left = BigMultShortLong(n, k, 2*n + 2 + LOGK); // 2k-1 registers in [0, k * 2^{2n+2})
    for(var i=0; i<k; i++){
        left.a[i] <== 2*in[0][1][i];
        left.b[i] <== in[0][1][i] + in[1][1][i];  
    }
    
    // prime reduce right - left 
    component diff_red[2]; 
    for(var i=0; i<2; i++)
        diff_red[i] = PrimeReduce(n, k, 2*k-2, p, 4*n + 3*LOGK + 3);
    for(var i=0; i<3*k-2; i++){
        diff_red[0].in[i] <== right.out[0][i]; 
        if(i < 2*k-1) 
            diff_red[1].in[i] <== right.out[1][i] + left.out[i]; 
        else
            diff_red[1].in[i] <== right.out[1][i];
    }
    // inputs of diff_red has registers in [0, 3*k^2*2^{3n} + k*2^{2n} + k*2^{2n+2} < 4*k^2*2^{3n} <= 2^{3n+2LOGK+2}) 
    // diff_red.out has registers < (2k-1)*2^{4n + 2LOGK + 2} <= 2^{4n + 3LOGK + 3}  
    component constraint = CheckCarryModToZero(n, k, 4*n + 3*LOGK + 3, p);
    for(var i=0; i<k; i++){
        constraint.in[0][i] <== diff_red[0].out[i];
        constraint.in[1][i] <== diff_red[1].out[i]; 
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
template EllipticCurveAddUnequal(n, k, p) { // changing q's to p's for my sanity
    signal input a[2][k];
    signal input b[2][k];

    signal output out[2][k];

    var LOGK = log_ceil(k);
    assert(4*n + 3*LOGK +4 < 253);

    // precompute lambda and x_3 and then y_3
    var dy[50] = long_sub_mod(n, k, b[1], a[1], p);
    var dx[50] = long_sub_mod(n, k, b[0], a[0], p); 
    var dx_inv[50] = mod_inv(n, k, dx, p);
    var lambda[50] = prod_mod(n, k, dy, dx_inv, p);
    var lambda_sq[50] = prod_mod(n, k, lambda, lambda, p);
    // out[0] = x_3 = lamb^2 - a[0] - b[0] % p
    // out[1] = y_3 = lamb * (a[0] - x_3) - a[1] % p
    var x3[50] = long_sub_mod(n, k, long_sub_mod(n, k, lambda_sq, a[0], p), b[0], p);
    var y3[50] = long_sub_mod(n, k, prod_mod(n, k, lambda, long_sub_mod(n, k, a[0], x3, p), p), a[1], p);

    for(var i = 0; i < k; i++){
        out[0][i] <-- x3[i];
        out[1][i] <-- y3[i];
    }
    
    // constrain x_3 by CUBIC (x_1 + x_2 + x_3) * (x_2 - x_1)^2 - (y_2 - y_1)^2 = 0 mod p
    
    component dx_sq = BigMultNoCarry(n, k, 2*n+1+LOGK); // 2k-1 registers in [0, k*2^{2n+1} )
    component dy_sq = BigMultNoCarry(n, k, 2*n+1+LOGK); // 2k-1 registers in [0, k*2^{2n+1} )
    for(var i = 0; i < k; i++){
        dx_sq.a[0][i] <== b[0][i];
        dx_sq.a[1][i] <== a[0][i];
        dx_sq.b[0][i] <== b[0][i];
        dx_sq.b[1][i] <== a[0][i];

        dy_sq.a[0][i] <== b[1][i];
        dy_sq.a[1][i] <== a[1][i];
        dy_sq.b[0][i] <== b[1][i];
        dy_sq.b[1][i] <== a[1][i];
    } 

    // x_1 + x_2 + x_3 has registers in [0, 3*2^n) 
    component cubic = BigMultNoCarry(n, 2*k-1, 3*n+4+2*LOGK); // 3k-2 registers in [0, 3 * k^2 * 2^{3n+1} ) 
    for(var i=0; i<2*k-1; i++){
        if(i < k)
            cubic.a[0][i] <== a[0][i] + b[0][i] + out[0][i]; 
        else
            cubic.a[0][i] <== 0;
        cubic.a[1][i] <== 0;
 
        cubic.b[0][i] <== dx_sq.out[0][i];
        cubic.b[1][i] <== dx_sq.out[1][i];
    }

    component cubic_red[2]; 
    for(var j=0; j<2; j++){
        cubic_red[j] = PrimeReduce(n, k, 2*k-2, p, 4*n + 3*LOGK + 4);
        for(var i=0; i<2*k-1; i++)
            cubic_red[j].in[i] <== cubic.out[j][i] + dy_sq.out[j ^ 1][i]; // registers in [0, 3*k^2*2^{3n+1} + k*2^{2n+1} < 2^{3n+2LOGK+3} )
        // j ^ 1 flips the bit so has the effect of subtracting  
        for(var i=2*k-1; i<3*k-2; i++)
            cubic_red[j].in[i] <== cubic.out[j][i]; 
    }
    // cubic_red has 2 x k registers in [0, (2k-1) 2^{4n+2LOGK+3} < 2^{4n+3LOGK+4} )
    
    component cubic_mod = CheckCarryModToZero(n, k, 4*n + 3*LOGK + 4, p);
    for(var j=0; j<2; j++)for(var i=0; i<k; i++)
        cubic_mod.in[j][i] <== cubic_red[j].out[i]; 
    // END OF CONSTRAINING x3
    
    // constrain y_3 by (y_1 + y_3) * (x_2 - x_1) = (y_2 - y_1)*(x_1 - x_3) mod p
    component y_constraint = PointOnLine(n, k, p); // 2k-1 registers in [0, k*2^{2n+1})
    for(var i = 0; i < k; i++)for(var j=0; j<2; j++){
        y_constraint.in[0][j][i] <== a[j][i];
        y_constraint.in[1][j][i] <== b[j][i];
        y_constraint.in[2][j][i] <== out[j][i];
    }
    // END OF CONSTRAINING y3

    // check if out[][] has registers in [0, 2^n) and each out[i] is in [0, p)
    // re-using Fp2 code by considering (x_3, y_3) as a 2d-vector over Fp
    component range_check = CheckValidFp2(n, k, p);
    for(var j=0; j<2; j++)for(var i=0; i<k; i++)
        range_check.in[j][i] <== out[j][i];
}


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
template EllipticCurveDouble(n, k, a, b, p) {
    signal input in[2][k];
    signal output out[2][k];

    /* 
    var LOGK = log_ceil(k);
    assert(4*n + 3*LOGK + 3 < 252);
    */

    var long_a[k];
    var long_3[k];
    long_a[0] = a;
    long_3[0] = 3;
    for (var i = 1; i < k; i++) {
        long_a[i] = 0;
        long_3[i] = 0;
    }

    // precompute lambda 
    var lamb_num[50] = long_add_mod(n, k, long_a, prod_mod(n, k, long_3, prod_mod(n, k, in[0], in[0], p), p), p);
    var lamb_denom[50] = long_add_mod(n, k, in[1], in[1], p);
    var lamb[50] = prod_mod(n, k, lamb_num, mod_inv(n, k, lamb_denom, p), p);

    // precompute x_3, y_3
    var x3[50] = long_sub_mod(n, k, prod_mod(n, k, lamb, lamb, p), long_add_mod(n, k, in[0], in[0], p), p);
    var y3[50] = long_sub_mod(n, k, prod_mod(n, k, lamb, long_sub_mod(n, k, in[0], x3, p), p), in[1], p);
    
    for(var i=0; i<k; i++){
        out[0][i] <-- x3[i];
        out[1][i] <-- y3[i];
    }
    // check if out[][] has registers in [0, 2^n) and each out[i] is in [0, p)
    // re-using Fp2 code by considering (x_3, y_3) as a 2d-vector over Fp
    component range_check = CheckValidFp2(n, k, p);
    for(var j=0; j<2; j++)for(var i=0; i<k; i++)
        range_check.in[j][i] <== out[j][i];

    component point_on_tangent = PointOnTangent(n, k, a, p);
    for(var j=0; j<2; j++)for(var i=0; i<k; i++){
        point_on_tangent.in[0][j][i] <== in[j][i];
        point_on_tangent.in[1][j][i] <== out[j][i];
    }
    
    component point_on_curve = PointOnCurve(n, k, a, b, p);
    for(var j=0; j<2; j++)for(var i=0; i<k; i++)
        point_on_curve.in[j][i] <== out[j][i];
    
    component x3_eq_x1 = IsArrayEqual(k);
    for(var i = 0; i < k; i++){
        x3_eq_x1.in[0][i] <== out[0][i];
        x3_eq_x1.in[1][i] <== in[0][i];
    }
    x3_eq_x1.out === 0;
}

// Inputs:
//  P is 2 x 2 x k array where P0 = (x_1, y_1) and P1 = (x_2, y_2) are points in E(Fp)
//  Q is 2 x 6 x 2 x k array representing point (X, Y) in E(Fp12)
// Assuming (x_1, y_1) != (x_2, y_2)
// Output:
//  out is 6 x 4 x (2k-1) array representing element of Fp12 (keeping track of negatives) equal to:
//  (y_1 - y_2) X + (x_2 - x_1) Y + (x_1 y_2 - x_2 y_1)
// We evaluate out without carries
// If all registers of P, Q are in [0, 2^n),
// Then all registers of out in [0, 2^{2n + log(k) + 2} )
//  more specifically: registers out[0][0 or 2][idx] are in [0, 3*2^{2n + log(k)} )
//                     all other registers are in [0, 2^{2n + log(k) + 1}) 
// m_out is the expected max number of bits in the output registers
template LineFunctionUnequalNoCarry(n, k, m_out){
    signal input P[2][2][k];
    signal input Q[2][6][2][k];
    signal output out[6][4][2*k-1];

    // (y_1 - y_2) X
    var LOGK = log_ceil(k);
    component Xmult = Fp12ScalarMultiplyNoCarry(n, k, 2*n + LOGK); // registers in [0, k*2^{2n} )
    // (x_2 - x_1) Y
    component Ymult = Fp12ScalarMultiplyNoCarry(n, k, 2*n + LOGK);
    for(var i=0; i<k; i++){
        Xmult.a[0][i] <== P[0][1][i];
        Xmult.a[1][i] <== P[1][1][i];
        
        Ymult.a[0][i] <== P[1][0][i];
        Ymult.a[1][i] <== P[0][0][i];
    }
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        Xmult.b[i][j][idx] <== Q[0][i][j][idx];
        Xmult.b[i][j+2][idx] <== 0; // expecting 0s to be optimized out by compiler
        
        Ymult.b[i][j][idx] <== Q[1][i][j][idx]; 
        Ymult.b[i][j+2][idx] <== 0;
    } 
    
    component x1y2 = BigMultShortLong(n, k, 2*n + LOGK); // registers in [0, k*2^{2n}) 
    component x2y1 = BigMultShortLong(n, k, 2*n + LOGK);
    for(var i=0; i<k; i++){
        x1y2.a[i] <== P[0][0][i]; 
        x1y2.b[i] <== P[1][1][i];
        
        x2y1.a[i] <== P[1][0][i]; 
        x2y1.b[i] <== P[0][1][i];
    }
    
    for(var i=0; i<6; i++)for(var j=0; j<4; j++)for(var idx=0; idx<2*k-1; idx++){
        if( i==0 && (j==0 || j==2) ){
            if(j==0)
                out[i][j][idx] <== Xmult.out[i][j][idx] + Ymult.out[i][j][idx] + x1y2.out[idx]; // register in [0, 3*k*2^{2n} ) 
            else
                out[i][j][idx] <== Xmult.out[i][j][idx] + Ymult.out[i][j][idx] + x2y1.out[idx];
        }else 
            out[i][j][idx] <== Xmult.out[i][j][idx] + Ymult.out[i][j][idx]; // register in [0, k*2^{2n+1} )
    }
    component range_checks[6][4][2*k-1];
    for (var outer = 0; outer < 6; outer ++) {
        for (var i = 0; i < 4; i ++) {
            for (var j = 0; j < 2*k-1; j ++) {
                range_checks[outer][i][j] = Num2Bits(m_out);
                range_checks[outer][i][j].in <== out[outer][i][j];
            }
        }
    }
}

// Assuming curve is of form Y^2 = X^3 + b for now (a = 0) for better register bounds 
// Inputs:
//  P is 2 x k array where P = (x, y) is a point in E(Fp) 
//  Q is 2 x 6 x 2 x k array representing point (X, Y) in E(Fp12) 
// Output: 
//  out is 6 x 4 x (3k-2) array representing element of Fp12 (keeping track of negatives) equal to:
//  3 x^2 (-X + x) + 2 y (Y - y)
// We evaluate out without carries 
// If P, Q have registers in [0, 2^n) 
// Then out has registers in [0, 3k^2*2^{3n} + 2k*2^{2n} < (3k + 2/2^n )*k*2^{3n})
// m_out is the expected max number of bits in the output registers
template LineFunctionEqualNoCarry(n, k, m_out){
    signal input P[2][k]; 
    signal input Q[2][6][2][k];
    signal output out[6][4][3*k-2];
    var LOGK = log_ceil(k);

    component x_sq3 = BigMultShortLong(n, k, 2*n + 2 + LOGK); // 2k-1 registers in [0, 3*k*2^{2n} )
    for(var i=0; i<k; i++){
        x_sq3.a[i] <== 3*P[0][i];
        x_sq3.b[i] <== P[0][i];
    } 
    
    // 3 x^2 (-X + x)
    component Xmult = Fp12ScalarMultiplyNoCarryUnequal(n, 2*k-1, k, 3*n + 2*LOGK + 2); // 3k-2 registers in [0, 3 * k^2 * 2^{3n})
    for(var idx=0; idx<2*k-1; idx++){
        Xmult.a[0][idx] <== x_sq3.out[idx];
        Xmult.a[1][idx] <== 0;
    }
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        Xmult.b[i][j+2][idx] <== Q[0][i][j][idx]; 
        if(i==0 && j==0)
            Xmult.b[i][j][idx] <== P[0][idx];
        else
            Xmult.b[i][j][idx] <== 0;
    }

    // 2 y (Y-y)
    component Ymult = Fp12ScalarMultiplyNoCarry(n, k, 2*n + LOGK + 1); // 2k-1 registers in [0, 2*k*2^{2n}) 
    for(var idx=0; idx < k; idx++){
        Ymult.a[0][idx] <== 2*P[1][idx];
        Ymult.a[1][idx] <== 0;
    }
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        Ymult.b[i][j][idx] <== Q[1][i][j][idx]; 
        if(i==0 && j==0)
            Ymult.b[i][j+2][idx] <== P[1][idx];
        else
            Ymult.b[i][j+2][idx] <== 0;
    }
    
    for(var i=0; i<6; i++)for(var j=0; j<4; j++)for(var idx=0; idx<3*k-2; idx++){
        if(idx < 2*k-1)
            out[i][j][idx] <== Xmult.out[i][j][idx] + Ymult.out[i][j][idx];
        else
            out[i][j][idx] <== Xmult.out[i][j][idx];
    }
    component range_checks[6][4][3*k-2];
    for (var outer = 0; outer < 6; outer ++) {
        for (var i = 0; i < 4; i ++) {
            for (var j = 0; j < 3*k-2; j ++) {
                range_checks[outer][i][j] = Num2Bits(m_out);
                range_checks[outer][i][j].in <== out[outer][i][j];
            }
        }
    }
}

// Inputs:
//  P is 2 x 2 x k array where P0 = (x_1, y_1) and P1 = (x_2, y_2) are points in E(Fp)
//  Q is 2 x 6 x 2 x k array representing point (X, Y) in E(Fp12)
// Assuming (x_1, y_1) != (x_2, y_2)
// Output:
//  Q is 6 x 2 x k array representing element of Fp12 equal to:
//  (y_1 - y_2) X + (x_2 - x_1) Y + (x_1 y_2 - x_2 y_1)
template LineFunctionUnequal(n, k, q) {
    signal input P[2][2][k];
    signal input Q[2][6][2][k];

    signal output out[6][2][k];
    var LOGK = log_ceil(k);

    component nocarry = LineFunctionUnequalNoCarry(n, k, 2 + 2 * n + LOGK);
    for (var i = 0; i < 2; i++) {
	for (var j = 0; j < 2; j++) {
	    for (var idx = 0; idx < k; idx++) {
		nocarry.P[i][j][idx] <== P[i][j][idx];
	    }
	}
    }

    for (var i = 0; i < 2; i++) {
	for (var j = 0; j < 6; j++) {
	    for (var l = 0; l < 2; l++) {
		for (var idx = 0; idx < k; idx++) {
		    nocarry.Q[i][j][l][idx] <== Q[i][j][l][idx];
		}
	    }
	}
    }
    component reduce[6][4];
    for (var i = 0; i < 6; i++) {
	for (var j = 0; j < 4; j++) {
	    reduce[i][j] = PrimeReduce(n, k, k - 1, q, 3 * n + 2 * LOGK + 2);
	}

	for (var j = 0; j < 4; j++) {
	    for (var idx = 0; idx < 2 * k - 1; idx++) {
		reduce[i][j].in[idx] <== nocarry.out[i][j][idx];
	    }
	}	
    }

    // max overflow register size is 3 * k * 2^{3n + log(k)}
    component carry = Fp12CarryModP(n, k, 3 * n + 2 * LOGK + 2, q);
    for (var i = 0; i < 6; i++) {
	for (var j = 0; j < 4; j++) {
	    for (var idx = 0; idx < k; idx++) {
		carry.in[i][j][idx] <== reduce[i][j].out[idx];
	    }
	}
    }
    
    for (var i = 0; i < 6; i++) {
	for (var j = 0; j < 2; j++) {
	    for (var idx = 0; idx < k; idx++) {
		out[i][j][idx] <== carry.out[i][j][idx];
	    }
	}
    }    
}


// Assuming curve is of form Y^2 = X^3 + b for now (a = 0) for better register bounds 
// Inputs:
//  P is 2 x k array where P = (x, y) is a point in E(Fp) 
//  Q is 2 x 6 x 2 x k array representing point (X, Y) in E(Fp12) 
// Output: 
//  out is 6 x 2 x k array representing element of Fp12 equal to:
//  3 x^2 (-X + x) + 2 y (Y - y)
template LineFunctionEqual(n, k, q) {
    signal input P[2][k];
    signal input Q[2][6][2][k];

    signal output out[6][2][k];

    var LOGK = log_ceil(k);
    component nocarry = LineFunctionEqualNoCarry(n, k, 3*n + 2*LOGK + 2);
    for (var i = 0; i < 2; i++) {
	for (var idx = 0; idx < k; idx++) {
	    nocarry.P[i][idx] <== P[i][idx];
	}
    }

    for (var i = 0; i < 2; i++) {
	for (var j = 0; j < 6; j++) {
	    for (var l = 0; l < 2; l++) {
		for (var idx = 0; idx < k; idx++) {
		    nocarry.Q[i][j][l][idx] <== Q[i][j][l][idx];
		}
	    }
	}
    }
    
    var logc = log_ceil((2*k-1)*(3*k+1)*k);
    component reduce[6][4]; 
    for (var i = 0; i < 6; i++) {
        for (var j = 0; j < 4; j++) {
            reduce[i][j] = PrimeReduce(n, k, 2 * k - 2, q, 4 * n + logc);
        }

        for (var j = 0; j < 4; j++) {
            for (var idx = 0; idx < 3 * k - 2; idx++) {
                reduce[i][j].in[idx] <== nocarry.out[i][j][idx];
            }
        }	
    }

    // max overflow register size is (2k - 1) * (3k+1)* k * 2^{4n}
    component carry = Fp12CarryModP(n, k, 4 * n + logc, q);
    for (var i = 0; i < 6; i++) {
	for (var j = 0; j < 4; j++) {
	    for (var idx = 0; idx < k; idx++) {
		carry.in[i][j][idx] <== reduce[i][j].out[idx];
	    }
	}
    }
    
    for (var i = 0; i < 6; i++) {
	for (var j = 0; j < 2; j++) {
	    for (var idx = 0; idx < k; idx++) {
		out[i][j][idx] <== carry.out[i][j][idx];
	    }
	}
    }    
}

// Input:
//  g is 6 x 4 x kg array representing element of Fp12, allowing overflow and negative
//  P0, P1, Q are as in inputs of LineFunctionUnequalNoCarry
// Assume:
//  all registers of g are in [0, 2^{overflowg}) 
//  all registers of P, Q are in [0, 2^n) 
// Output:
//  out = g * l_{P0, P1}(Q) as element of Fp12 with carry 
//  out is 6 x 2 x k
template Fp12MultiplyWithLineUnequal(n, k, kg, overflowg, q){
    signal input g[6][4][kg];
    signal input P[2][2][k];
    signal input Q[2][6][2][k];
    signal output out[6][2][k];

    var LOGK = log_ceil(k);
    component line = LineFunctionUnequalNoCarry(n, k, 2*n + LOGK + 2); // 6 x 4 x 2k - 1 registers in [0, 2^{2n + log(k) + 2})
    for(var l=0; l<2; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        line.P[l][j][idx] <== P[l][j][idx];
    for(var l=0; l<2; l++)for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        line.Q[l][i][j][idx] <== Q[l][i][j][idx];
    
    var mink;
    if(kg < 2*k - 1)
        mink = kg;
    else
        mink = 2*k - 1;
    var logc2 = log_ceil(12 * mink * k);
    component mult = Fp12MultiplyNoCarryUnequal(n, kg, 2*k - 1, overflowg + 2*n + logc2 + 2); // 6 x 4 x 2k + kg - 2 registers in [0, 12 * min(kg, 2k - 1) * 2^{overflowg + 2n + log(k) + 2} )
    
    for(var i=0; i<6; i++)for(var j=0; j<4; j++)for(var idx=0; idx<kg; idx++)
        mult.a[i][j][idx] <== g[i][j][idx];
    for(var i=0; i<6; i++)for(var j=0; j<4; j++)for(var idx=0; idx<2*k-1; idx++)
        mult.b[i][j][idx] <== line.out[i][j][idx];
    var logc = log_ceil(12 * (k + kg - 1) * mink * k);
    assert( overflowg + 3*n + logc + 2 < 252 );
    component reduce = Fp12Compress(n, k, k + kg - 2, q, overflowg + 3*n + logc + 2); // 6 x 4 x k registers in [0, 12 * (k + kg - 1) * min(kg, 2k - 1) * 2^{overflowg + 3n + log(k) + 2} )
    for(var i=0; i<6; i++)for(var j=0; j<4; j++)for(var idx=0; idx<2*k + kg - 2; idx++)
        reduce.in[i][j][idx] <== mult.out[i][j][idx];
    
    component carry = Fp12CarryModP(n, k, overflowg + 3*n + logc + 2, q);

    for(var i=0; i<6; i++)for(var j=0; j<4; j++)for(var idx=0; idx<k; idx++)
        carry.in[i][j][idx] <== reduce.out[i][j][idx];

    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[i][j][idx] <== carry.out[i][j][idx];
}

// Assuming curve is of form Y^2 = X^3 + b for now (a = 0) for better register bounds 
// Inputs:
//  in is 6 x 2 x k array representing element in Fq12
//  P is 2 x k array where P = (x, y) is a point in E[r](Fq) 
//  Q is 2 x 6 x 2 x k array representing point (X, Y) in E(Fq12) 
// Output:
//  out = f_x(P,Q) is 6 x 2 x k, where we start with f_0(P,Q) = in and use Miller's algorithm f_{i+j} = f_i * f_j * l_{i,j}(P,Q)
//  xP = [x]P is 2 x k array
// Assume:
//  r is prime (not a parameter in this template)
//  x in [0, 2^250) and x < r  (we will use this template when x has few significant bits in base 2)
//  q has k registers in [0, 2^n)
//  P != O so the order of P in E(Fq) is r, so [i]P != [j]P for i != j in Z/r 
template MillerLoop(n, k, b, x, q){
    signal input in[6][2][k];
    signal input P[2][k]; 
    signal input Q[2][6][2][k];

    signal output out[6][2][k];
    signal output xP[2][k];

    var Bits[250]; // length is k * n
    var BitLength;
    var SigBits=0;
    for (var i = 0; i < 250; i++) {
        Bits[i] = (x >> i) & 1;
        if(Bits[i] == 1){
            SigBits++;
            BitLength = i + 1;
        }
    }

    signal Pintermed[BitLength][2][k]; 
    signal f[BitLength][6][2][k];

    component Pdouble[BitLength];
    component fdouble[BitLength];
    component square[BitLength];
    component line[BitLength];
    component compress[BitLength];
    component nocarry[BitLength];
    component Padd[SigBits];
    component fadd[SigBits]; 
    var curid=0;

    var LOGK = log_ceil(k);
    var logc = log_ceil(144*k*k*k);
    var logc2 = log_ceil(144 * k * k);
    assert( 4*n + logc + 1 < 252 );
    
    for(var i=BitLength - 1; i>=0; i--){
        if( i == BitLength - 1 ){
            // f = 1 
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                f[i][l][j][idx] <== in[l][j][idx];
            }
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                Pintermed[i][j][idx] <== P[j][idx];
        }else{
            // compute fdouble[i] = f[i+1]^2 * l_{Pintermed[i+1], Pintermed[i+1]}(Q) 
            square[i] = Fp12MultiplyNoCarry(n, k, 2*n + 4 + LOGK); // 6 x 4 x 2k-1 registers in [0, 12 * k * 2^{2n} )
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                square[i].a[l][j][idx] <== f[i+1][l][j][idx];
                square[i].a[l][j+2][idx] <== 0;
                square[i].b[l][j][idx] <== f[i+1][l][j][idx];
                square[i].b[l][j+2][idx] <== 0;
            }

            line[i] = LineFunctionEqual(n, k, q); // 6 x 2 x k registers in [0, 2^n) 
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                line[i].P[j][idx] <== Pintermed[i+1][j][idx];            
            for(var eps=0; eps<2; eps++)for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                line[i].Q[eps][l][j][idx] <== Q[eps][l][j][idx];

            nocarry[i] = Fp12MultiplyNoCarryUnequal(n, 2*k-1, k, logc2 + 3*n); // 6 x 4 x 3k-2 registers in [0, 144 * k^2 * 2^{3n} ) 
            for(var l=0; l<6; l++)for(var j=0; j<4; j++)for(var idx=0; idx<2*k-1; idx++)
                nocarry[i].a[l][j][idx] <== square[i].out[l][j][idx];
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                nocarry[i].b[l][j][idx] <== line[i].out[l][j][idx];
                nocarry[i].b[l][j+2][idx] <== 0;
            }
            compress[i] = Fp12Compress(n, k, 2*k-2, q, logc + 1 + 4*n); // 6 x 4 x k registers in [0, 144 * k^2 * (2k-1) * 2^{4n} )
            for(var l=0; l<6; l++)for(var j=0; j<4; j++)for(var idx=0; idx<3*k-2; idx++)
                compress[i].in[l][j][idx] <== nocarry[i].out[l][j][idx];

            fdouble[i] = Fp12CarryModP(n, k, logc + 1 + 4*n, q);
            for(var l=0; l<6; l++)for(var j=0; j<4; j++)for(var idx=0; idx<k; idx++)
                fdouble[i].in[l][j][idx] <== compress[i].out[l][j][idx]; 

            Pdouble[i] = EllipticCurveDouble(n, k, 0, b, q);  
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                Pdouble[i].in[j][idx] <== Pintermed[i+1][j][idx]; 
            
            if(Bits[i] == 0){
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    f[i][l][j][idx] <== fdouble[i].out[l][j][idx]; 
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    Pintermed[i][j][idx] <== Pdouble[i].out[j][idx];
            }else{
                // fadd[curid] = fdouble * l_{Pdouble[i], P}(Q) 
                fadd[curid] = Fp12MultiplyWithLineUnequal(n, k, k, n, q); 
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                    fadd[curid].g[l][j][idx] <== fdouble[i].out[l][j][idx];
                    fadd[curid].g[l][j+2][idx] <== 0;
                }
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                    fadd[curid].P[0][j][idx] <== Pdouble[i].out[j][idx];            
                    fadd[curid].P[1][j][idx] <== P[j][idx];            
                }
                for(var eps=0; eps<2; eps++)for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    fadd[curid].Q[eps][l][j][idx] <== Q[eps][l][j][idx];

                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    f[i][l][j][idx] <== fadd[curid].out[l][j][idx]; 

                // Padd[curid] = Pdouble[i] + P 
                Padd[curid] = EllipticCurveAddUnequal(n, k, q); 
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                    Padd[curid].a[j][idx] <== Pdouble[i].out[j][idx];
                    Padd[curid].b[j][idx] <== P[j][idx];
                }

                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    Pintermed[i][j][idx] <== Padd[curid].out[j][idx];
                
                curid++;
            }
        }
    }
    for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[l][j][idx] <== f[0][l][j][idx];
    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        xP[j][idx] <== Pintermed[0][j][idx]; 
    
}


// Assuming curve is of form Y^2 = X^3 + b for now (a = 0) for better register bounds 
// Inputs:
//  P is 2 x k array where P = (x, y) is a point in E[r](Fq) 
//  Q is 2 x 6 x 2 x k array representing point (X, Y) in E(Fq12) 
// Output:
// f_r(Q) where <f_r> = [r]P - [r]O is computed using Miller's algorithm
// Assume:
//  r  = x^4 - x^2 + 1 where x is the parameter of the curve
//  q has k registers in [0, 2^n)
//  r is prime
//  P != O so the order of P in E(Fq) is r, so [i]P != [j]P for i != j in Z/r 
template MillerLoopOptimized(n, k, q){
    signal input P[2][k]; 
    signal input Q[2][6][2][k];
    signal output out[6][2][k];

    var b = 4; // Y^2 = X^3 + 4
    var x = get_BLS12_381_parameter();

    // fx[i] = f_{x^{i+1}} 
    component fx[4]; 
    for(var e=0; e<4; e++){
        fx[e] = MillerLoop(n, k, b, x, q);
        if( e == 0 ){
            for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                if(i == 0 && j == 0 && idx == 0)
                    fx[e].in[i][j][idx] <== 1;
                else    
                    fx[e].in[i][j][idx] <== 0;
            }
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                fx[e].P[j][idx] <== P[j][idx];
        }else{
            for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                fx[e].in[i][j][idx] <== fx[e - 1].out[i][j][idx];
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                fx[e].P[j][idx] <== fx[e - 1].xP[j][idx];            
        }
        for(var l=0; l<2; l++)for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
            fx[e].Q[l][i][j][idx] <== Q[l][i][j][idx];
    }
    
    // f_{x^4} * l_{x^4, 1}(P,Q) 
    component fx4l = Fp12MultiplyWithLineUnequal(n, k, k, n, q); 
    // assert( 4*n + log_ceil(12 * (2*k-1) *k * k) + 2 < 252 );  // need this to run MillerLoop anyways
    for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        fx4l.g[l][j][idx] <== fx[3].out[l][j][idx];
        fx4l.g[l][j+2][idx] <== 0;
    }
    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        fx4l.P[0][j][idx] <== fx[3].xP[j][idx];            
        fx4l.P[1][j][idx] <== P[j][idx];            
    }
    for(var l=0; l<2; l++)for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        fx4l.Q[l][i][j][idx] <== Q[l][i][j][idx];
    
    // f_{x^2} * l_{r,x^2}(P,Q) where l_{r,x^2}(P,Q) = Q.x - ([x^2]P).x 
    var LOGK = log_ceil(k);
    component fx2l = Fp12MultiplyNoCarryCompress(n, k, q, n, 3*n + 2*LOGK + 4); // registers in [0, 12*k^2*2^{3n} )
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        fx2l.a[i][j][idx] <== fx[1].out[i][j][idx];
        fx2l.a[i][j+2][idx] <== 0;
        
        fx2l.b[i][j][idx] <== Q[0][i][j][idx];
        if(i == 0 && j == 0)
            fx2l.b[i][j+2][idx] <== fx[1].xP[0][idx];
        else
            fx2l.b[i][j+2][idx] <== 0;
    }

    // find fx2l^{-1}. Not going to optimize this for now since it's just one call
    component carry = Fp12CarryModP(n, k, 3*n + 2*LOGK + 4, q);
    for(var i=0; i<6; i++)for(var j=0; j<4; j++)for(var idx=0; idx<k; idx++)
        carry.in[i][j][idx] <== fx2l.out[i][j][idx];
    
    component inv = Fp12Invert(n, k, q);
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        inv.in[i][j][idx] <== carry.out[i][j][idx];

    component fr = Fp12Multiply2(n, k, q);
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        fr.a[i][j][idx] <== fx4l.out[i][j][idx];
        fr.b[i][j][idx] <== inv.out[i][j][idx];
    }
    
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[i][j][idx] <== fr.out[i][j][idx];
}
