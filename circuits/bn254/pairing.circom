pragma circom 2.0.3;

include "../curve.circom";
include "curve_fp2.circom";
include "fp12.circom";
include "final_exp.circom";
include "bn254_func.circom";

// Inputs:
//  P is 2 x 2 x 2 x k array where P0 = (x_1, y_1) and P1 = (x_2, y_2) are points in E(Fp2)
//  Q is 2 x k array representing point (X, Y) in E(Fp)
// Assuming (x_1, y_1) != (x_2, y_2)
// Output:
//  out is 6 x 2 x (2k-1) array representing element of Fp12 equal to:
//  - line_{Psi(P0), Psi(P1)}(Q) where Psi(x,y) = (w^2 x, w^3 y)
//  - equals w^3 (y_1 - y_2) X + w^2 (x_2 - x_1) Y + w^5 (x_1 y_2 - x_2 y_1)
// We evaluate out without carries
// If all registers of P, Q are in [0, B),
// Then all registers of out have abs val < 2k * B^2 )
// m_out is the expected max number of bits in the output registers
template SignedLineFunctionUnequalNoCarryFp2(n, k, m_out){
    signal input P[2][2][2][k];
    signal input Q[2][k];
    signal output out[6][2][2*k-1];

    // (y_1 - y_2) X
    var LOGK = log_ceil(k);
    component Xmult = SignedFp2MultiplyNoCarry(n, k, 2*n + LOGK+1); // registers abs val < 2k*B^2
    // (x_2 - x_1) Y
    component Ymult = SignedFp2MultiplyNoCarry(n, k, 2*n + LOGK+1);
    for(var i = 0; i < 2; i ++) {
        for(var j=0; j<k; j++){
            Xmult.a[i][j] <== P[0][1][i][j] - P[1][1][i][j];
            
            Ymult.a[i][j] <== P[1][0][i][j] - P[0][0][i][j];
        }
    }
    for(var idx=0; idx<k; idx++){
        Xmult.b[0][idx] <== Q[0][idx];
        Xmult.b[1][idx] <== 0;

        Ymult.b[0][idx] <== Q[1][idx]; 
        Ymult.b[1][idx] <== 0;
    } 
    
    component x1y2 = BigMultShortLong2D(n, k, 2); // 2k-1 registers in [0, 2k*B^2) 
    component x2y1 = BigMultShortLong2D(n, k, 2);
    for(var i = 0; i < 2; i ++) {
        for(var j=0; j<k; j++){
            x1y2.a[i][j] <== P[0][0][i][j]; 
            x1y2.b[i][j] <== P[1][1][i][j];
            
            x2y1.a[i][j] <== P[1][0][i][j]; 
            x2y1.b[i][j] <== P[0][1][i][j];
        }
    }
    
    for(var idx=0; idx<2*k-1; idx++){
        out[0][0][idx] <== 0;
        out[0][1][idx] <== 0;
        out[1][0][idx] <== 0;
        out[1][1][idx] <== 0;
        out[4][0][idx] <== 0;
        out[4][1][idx] <== 0;

        out[5][0][idx] <== x1y2.out[0][idx] - x2y1.out[0][idx] - x1y2.out[2][idx] + x2y1.out[2][idx];
        out[5][1][idx] <== x1y2.out[1][idx] - x2y1.out[1][idx];
        out[3][0][idx] <== Xmult.out[0][idx];
        out[3][1][idx] <== Xmult.out[1][idx];
        out[2][0][idx] <== Ymult.out[0][idx];
        out[2][1][idx] <== Ymult.out[1][idx];
    }
}


// Assuming curve is of form Y^2 = X^3 + b for now (a = 0) for better register bounds 
// Inputs:
//  P is 2 x 2 x k array where P = (x, y) is a point in E(Fp2) 
//  Q is 2 x k array representing point (X, Y) in E(Fp) 
// Output: 
//  out is 6 x 2 x (3k-2) array representing element of Fp12 equal to:
//  - line_{Psi(P), Psi(P)}(Q) where Psi(x,y) = (w^2 x, w^3 y)
//  - equals (3x^3 - 2y^2)(XI0 + u) + w^4 (-3 x^2 X) + w^3 (2 y Y)
// We evaluate out without carries, with signs
// If P, Q have registers in [0, B) 
// Then out has registers with abs val < (12k^2 + 6k/B) * (XI0 + 1) * B^3
// m_out is the expected max number of bits in the output registers
template SignedLineFunctionEqualNoCarryFp2(n, k, m_out){
    signal input P[2][2][k]; 
    signal input Q[2][k];
    signal output out[6][2][3*k-2];
    var XI0 = 9;
    var LOGK = log_ceil(k);

    // 3 x^2
    component x_sq3 = BigMultShortLong2D(n, k, 2); // 2k-1 registers in [0, 6*k*2^{2n} )
    for(var i=0; i<2; i++){
        for(var j = 0; j < k; j ++) {
            x_sq3.a[i][j] <== 3*P[0][i][j];
            x_sq3.b[i][j] <== P[0][i][j];
        }
    }

    // 3 x^3
    component x_cu3 = BigMultShortLong2DUnequal(n, 2*k-1, k, 2, 2); // 3k-2 registers abs val < 12*k^2*2^{3n} 
    for (var i = 0; i < 2*k-1; i ++) {
        if (i < k) {
            x_cu3.b[0][i] <== P[0][0][i];
            x_cu3.b[1][i] <== P[0][1][i];
        }
        x_cu3.a[0][i] <== x_sq3.out[0][i] - x_sq3.out[2][i];
        x_cu3.a[1][i] <== x_sq3.out[1][i];
    }

    // 2 y^2
    component y_sq2 = BigMultShortLong2D(n, k, 2); // 2k-1 registers in [0, 6*k*2^{2n} )
    for(var i=0; i<2; i++){
        for(var j = 0; j < k; j ++) {
            y_sq2.a[i][j] <== 2*P[1][i][j];
            y_sq2.b[i][j] <== P[1][i][j];
        }
    } 
    
    // - 3 x^2 X 
    component Xmult = SignedFp2MultiplyNoCarryUnequal(n, 2*k-1, k, 3*n + 2*LOGK + 3); // 3k-2 registers abs val < 12 * k^2 * 2^{3n})
    for(var idx=0; idx<2*k-1; idx++){
        Xmult.a[0][idx] <== x_sq3.out[0][idx] - x_sq3.out[2][idx];
        Xmult.a[1][idx] <== x_sq3.out[1][idx];
    }
    for(var idx=0; idx<k; idx++){
        Xmult.b[0][idx] <== - Q[0][idx];
        Xmult.b[1][idx] <== 0;
    }

    // 2 y Y
    component Ymult = SignedFp2MultiplyNoCarryUnequal(n, k, k, 2*n + LOGK + 2); // 2k-1 registers abs val < 4k*2^{2n} 
    for(var idx=0; idx < k; idx++){
        Ymult.a[0][idx] <== 2*P[1][0][idx];
        Ymult.a[1][idx] <== 2*P[1][1][idx];
    }
    for(var idx=0; idx<k; idx++){
        Ymult.b[0][idx] <== Q[1][idx];
        Ymult.b[1][idx] <== 0;
    }
    
    // 3 x^3 - 2 y^2 
    var w6coeff[2][50]; // registers abs val < (12k^2 + 6k/2^n) * 2^{3n} 
    for(var idx=0; idx<3*k-2; idx++){
        if(idx < 2*k-1){
            w6coeff[0][idx] = x_cu3.out[0][idx] - x_cu3.out[2][idx] - y_sq2.out[0][idx] + y_sq2.out[2][idx];
            w6coeff[1][idx] = x_cu3.out[1][idx] - y_sq2.out[1][idx];
        }else{
            w6coeff[0][idx] = x_cu3.out[0][idx] - x_cu3.out[2][idx];
            w6coeff[1][idx] = x_cu3.out[1][idx];
        }
    }
    // (3x^3 - 2y^2)(XI0 + u)
    var out0[2][50] = signed_Fp2_mult_w6(3*k-2, w6coeff, XI0); // registers abs val < (12k^2 + 6k/2^n) * (XI0 + 1) * 2^{3n}

    for(var idx=0; idx<3*k-2; idx++){
        out[1][0][idx] <== 0;
        out[1][1][idx] <== 0;
        out[2][0][idx] <== 0;
        out[2][1][idx] <== 0;
        out[5][0][idx] <== 0;
        out[5][1][idx] <== 0;

        if (idx < 2*k-1) {
            out[3][0][idx] <== Ymult.out[0][idx];
            out[3][1][idx] <== Ymult.out[1][idx];
        }
        else {
            out[3][0][idx] <== 0;
            out[3][1][idx] <== 0;
        }
        out[0][0][idx] <== out0[0][idx];
        out[0][1][idx] <== out0[1][idx];
        out[4][0][idx] <== Xmult.out[0][idx];
        out[4][1][idx] <== Xmult.out[1][idx];
    }
}

// Inputs:
//  P is 2 x 2 x k array where P0 = (x_1, y_1) and P1 = (x_2, y_2) are points in E(Fp2)
//  Q is 2 x 6 x 2 x k array representing point (X, Y) in E(Fp12)
// Assuming (x_1, y_1) != (x_2, y_2)
// Output:
//  Q is 6 x 2 x k array representing element of Fp12 equal to:
//  - line_{Psi(P0), Psi(P1)}(Q) where Psi(x,y) = (w^2 x, w^3 y)
//  - equals w^3 (y_1 - y_2) X + w^2 (x_2 - x_1) Y + w^5 (x_1 y_2 - x_2 y_1)
template LineFunctionUnequalFp2(n, k, q) {
    signal input P[2][2][2][k];
    signal input Q[2][k];

    signal output out[6][2][k];
    var LOGK1 = log_ceil(2*k);
    var LOGK2 = log_ceil(2*k*k);

    component nocarry = SignedLineFunctionUnequalNoCarryFp2(n, k, 2 * n + LOGK1);
    for (var i = 0; i < 2; i++)for(var j = 0; j < 2; j++) {
	    for (var idx = 0; idx < k; idx++) {
            nocarry.P[i][j][0][idx] <== P[i][j][0][idx];
            nocarry.P[i][j][1][idx] <== P[i][j][1][idx];
	    }
    }

    for (var i = 0; i < 2; i++) {
        for (var idx = 0; idx < k; idx++) {
            nocarry.Q[i][idx] <== Q[i][idx];
        }
    }
    component reduce[6][2];
    for (var i = 0; i < 6; i++) {
        for (var j = 0; j < 2; j++) {
            reduce[i][j] = PrimeReduce(n, k, k - 1, q, 3 * n + LOGK2);
        }

        for (var j = 0; j < 2; j++) {
            for (var idx = 0; idx < 2 * k - 1; idx++) {
                reduce[i][j].in[idx] <== nocarry.out[i][j][idx];
            }
        }	
    }

    // max overflow register size is 2k^2 * 2^{3n}
    component carry = SignedFp12CarryModP(n, k, 3 * n + LOGK2, q);
    for (var i = 0; i < 6; i++) {
        for (var j = 0; j < 2; j++) {
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
//  P is 2 x 2 x k array where P = (x, y) is a point in E(Fp2) 
//  Q is 2 x k array representing point (X, Y) in E(Fp) 
// Output: 
//  out is 6 x 2 x k array representing element of Fp12 equal to:
//  - line_{Psi(P), Psi(P)}(Q) where Psi(x,y) = (w^2 x, w^3 y)
//  - equals (3x^3 - 2y^2)(XI0 + u) + w^4 (-3 x^2 X) + w^3 (2 y Y)
template LineFunctionEqualFp2(n, k, p) {
    signal input P[2][2][k];
    signal input Q[2][k];

    signal output out[6][2][k];

    var XI0 = 9;
    var LOGK2 = log_ceil( (12*k*k + 1) * (XI0 + 1) );
    var LOGK3 = log_ceil((2*k-1)*12*k*k * (XI0 + 1) + 1);
    assert( 4*n + LOGK3 < 251 );

    component nocarry = SignedLineFunctionEqualNoCarryFp2(n, k, 3*n + LOGK2);
    for (var i = 0; i < 2; i++) {
        for (var j = 0; j < 2; j ++) {
            for (var idx = 0; idx < k; idx++) {
                nocarry.P[i][j][idx] <== P[i][j][idx];
            }
        }
    }

    for (var i = 0; i < 2; i++) {
        for (var idx = 0; idx < k; idx++) {
            nocarry.Q[i][idx] <== Q[i][idx];
        }
    }
    
    component reduce[6][4]; 
    for (var i = 0; i < 6; i++) {
        for (var j = 0; j < 2; j++) {
            reduce[i][j] = PrimeReduce(n, k, 2 * k - 2, p, 4 * n + LOGK3);
        }

        for (var j = 0; j < 2; j++) {
            for (var idx = 0; idx < 3 * k - 2; idx++) {
                reduce[i][j].in[idx] <== nocarry.out[i][j][idx];
            }
        }	
    }

    // max overflow register size is (2k - 1) * 12k^2 * 2^{4n}
    component carry = SignedFp12CarryModP(n, k, 4 * n + LOGK3, p);
    for (var i = 0; i < 6; i++) {
        for (var j = 0; j < 2; j++) {
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
//  g is 6 x 2 x kg array representing element of Fp12, allowing overflow and negative
//  P0, P1, Q are as in inputs of SignedLineFunctionUnequalNoCarryFp2
// Assume:
//  all registers of g are in [0, 2^{overflowg}) 
//  all registers of P, Q are in [0, 2^n) 
// Output:
//  out = g * l_{P0, P1}(Q) as element of Fp12 with carry 
//  out is 6 x 2 x k
template Fp12MultiplyWithLineUnequalFp2(n, k, kg, overflowg, q){
    signal input g[6][2][kg];
    signal input P[2][2][2][k];
    signal input Q[2][k];
    signal output out[6][2][k];

    var XI0 = 9;
    var LOGK1 = log_ceil(12*k);
    var LOGK2 = log_ceil(12*k * min(kg, 2*k-1) * 6 * (2+XI0) );
    var LOGK3 = log_ceil(12*k * min(kg, 2*k-1) * 6 * (2+XI0) * (k + kg - 1) );
    assert( overflowg + 3*n + LOGK3 < 251 );

    component line = SignedLineFunctionUnequalNoCarryFp2(n, k, 2*n + LOGK1); // 6 x 2 x 2k - 1 registers in [0, 12k 2^{2n})
    for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var l=0; l<2; l++)for(var idx=0; idx<k; idx++)
        line.P[i][j][l][idx] <== P[i][j][l][idx];
    for(var l=0; l<2; l++)for(var idx=0; idx<k; idx++)
        line.Q[l][idx] <== Q[l][idx];
    
    component mult = SignedFp12MultiplyNoCarryUnequal(n, kg, 2*k - 1, overflowg + 2*n + LOGK2); // 6 x 2 x (2k + kg - 2) registers abs val < 12k * min(kg, 2k - 1) * 6 * (2+XI0)* 2^{overflowg + 2n} 
    
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<kg; idx++)
        mult.a[i][j][idx] <== g[i][j][idx];
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<2*k-1; idx++)
        mult.b[i][j][idx] <== line.out[i][j][idx];


    component reduce = Fp12Compress(n, k, k + kg - 2, q, overflowg + 3*n + LOGK3); // 6 x 2 x k registers abs val < 12 k * min(kg, 2k - 1) * 6*(2+XI0) * (k + kg - 1) *  2^{overflowg + 3n} 
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<2*k + kg - 2; idx++)
        reduce.in[i][j][idx] <== mult.out[i][j][idx];
    
    component carry = SignedFp12CarryModP(n, k, overflowg + 3*n + LOGK3, q);

    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        carry.in[i][j][idx] <== reduce.out[i][j][idx];

    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[i][j][idx] <== carry.out[i][j][idx];
}

// Input:
//  g is 6 x 2 x k array representing element of Fp12
//  P, Q are as in inputs of SignedLineFunctionEqualNoCarryFp2
// Assume:
//  all registers of g are in [0, 2^n) 
//  all registers of P, Q are in [0, 2^n) 
//  36k^2 (XI0 + 1) * (XI0 + 2) * (3k-2) <= 2^n
// Output:
//  out = g * l_{P, P}(Q) as element of Fp12 with carry 
//  out is 6 x 2 x k
template Fp12MultiplyWithLineEqualFp2(n, k, p){
    signal input g[6][2][k];
    signal input P[2][2][k];
    signal input Q[2][k];
    signal output out[6][2][k];

    var XI0 = 9;
    var LOGK2 = log_ceil(12*k*k * (XI0 + 1) + 1);
    var LOGK3 = log_ceil(12*k*k * (XI0 + 1) * 6 * k * (2+XI0) + 1);
    var LOGK4 = log_ceil(12*k*k * (XI0 + 1) * 6 * k * (2+XI0) * (3*k-2) + 1);
    assert( 5*n + LOGK4 < 251 );

    component line = SignedLineFunctionEqualNoCarryFp2(n, k, 3*n + LOGK2); // 6 x 2 x 3k - 2 registers abs val < (12k^2 * (XI0 + 1) + 1) 2^{3n}
    for(var j=0; j<2; j++)for(var l=0; l<2; l++)for(var idx=0; idx<k; idx++)
        line.P[j][l][idx] <== P[j][l][idx];
    for(var l=0; l<2; l++)for(var idx=0; idx<k; idx++)
        line.Q[l][idx] <== Q[l][idx];
    
    component mult = SignedFp12MultiplyNoCarryUnequal(n, k, 3*k - 2, 4*n + LOGK3); // 6 x 2 x (4k - 3) registers abs val < 12k^2 * (XI0 + 1) * 6 * (2+XI0) * 2^{4n} 
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        mult.a[i][j][idx] <== g[i][j][idx];
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<3*k-2; idx++)
        mult.b[i][j][idx] <== line.out[i][j][idx];


    component reduce = Fp12Compress(n, k, 3*k-3, p, 5*n + LOGK4); // 6 x 2 x k registers abs val < 12 k^2 (XI0+1) * 6*(2+XI0) * (3k - 2) *  2^{5n} 
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<4*k-3; idx++)
        reduce.in[i][j][idx] <== mult.out[i][j][idx];
    
    component carry = SignedFp12CarryModP(n, k, 5*n + LOGK4, p);
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        carry.in[i][j][idx] <== reduce.out[i][j][idx];

    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[i][j][idx] <== carry.out[i][j][idx];
}


// Assuming curve is of form Y^2 = X^3 + b for now (a = 0) for better register bounds 
// b = [b0, b1] is 2 x k array representing Fp2 element
// Inputs:
//  P is 2 x 2 x k array where P = (x, y) is a point in E[r](Fp2) 
//  Q is 2 x k array representing point (X, Y) in E(Fp) 
// Output:
//  out = f_loopCount(P,Q) * l_{[loopCount] P', Frob_p(P')}(Q) * l_{[loopCount] P' + Frob_p(P'), -Frob_p^2(P')}(Q) 
//  - is 6 x 2 x k, where we start with f_0(P,Q) = in and use Miller's algorithm f_{i+j} = f_i * f_j * l_{i,j}(P,Q)
//  - where P' = twist(P) in E(Fp12)
//  - Frob_p(x,y) = (x^p, y^p)
// Assume:
//  r is prime (not a parameter in this template)
//  loopCount in [0, 2^250) and loopCount < r and loopCount < p
//  p has k registers in [0, 2^n)
//  P != O so the order of P in E(Fp) is r, so [i]P != [j]P for i != j in Z/r 
//  X^3 + b = 0 has no solution in Fp2, i.e., the y-coordinate of P cannot be 0.
template MillerLoopFp2(n, k, b0, b1, loopCount, p){
    signal input P[2][2][k]; 
    signal input Q[2][k];

    signal output out[6][2][k];

    var LOGK = log_ceil(k);
    var XI0 = 9;
    var LOGK2 = log_ceil(36*(2+XI0)*(2+XI0) * k*k);
    var LOGK3 = log_ceil(36*(2+XI0)*(2+XI0) * k*k*(2*k-1));
    assert( 4*n + LOGK3 < 251 );
    
    var Bits[250]; 
    var BitLength;
    var SigBits=0;
    for (var i = 0; i < 250; i++) {
        Bits[i] = (loopCount >> i) & 1;
        if(Bits[i] == 1){
            SigBits++;
            BitLength = i + 1;
        }
    }

    signal R[BitLength][2][2][k]; 
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

    for(var i=BitLength - 1; i>=0; i--){
        if( i == BitLength - 1 ){
            // f = 1 
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                if(l==0 && j==0 && idx==0)
                    f[i][l][j][idx] <== 1;
                else    
                    f[i][l][j][idx] <== 0;
            }
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                R[i][j][l][idx] <== P[j][l][idx];
        }else{
            // compute fdouble[i] = f[i+1]^2 * l_{R[i+1], R[i+1]}(Q) 
            square[i] = SignedFp12MultiplyNoCarry(n, k, 2*n + 4 + LOGK); // 6 x 2 x 2k-1 registers in [0, 6 * k * (2+XI0) * 2^{2n} )
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                square[i].a[l][j][idx] <== f[i+1][l][j][idx];
                square[i].b[l][j][idx] <== f[i+1][l][j][idx];
            }
            
            line[i] = LineFunctionEqualFp2(n, k, p); // 6 x 2 x k registers in [0, 2^n) 
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                line[i].P[j][l][idx] <== R[i+1][j][l][idx];            
            for(var eps=0; eps<2; eps++)
                for(var idx=0; idx<k; idx++)
                    line[i].Q[eps][idx] <== Q[eps][idx];
            
            nocarry[i] = SignedFp12MultiplyNoCarryUnequal(n, 2*k-1, k, 3*n + LOGK2); // 6 x 2 x 3k-2 registers < (6 * (2+XI0))^2 * k^2 * 2^{3n} ) 
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<2*k-1; idx++)
                nocarry[i].a[l][j][idx] <== square[i].out[l][j][idx];
            
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                nocarry[i].b[l][j][idx] <== line[i].out[l][j][idx];
            
            compress[i] = Fp12Compress(n, k, 2*k-2, p, 4*n + LOGK3); // 6 x 2 x k registers < (6 * (2+ XI0))^2 * k^2 * (2k-1) * 2^{4n} )
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<3*k-2; idx++)
                compress[i].in[l][j][idx] <== nocarry[i].out[l][j][idx];
            
            fdouble[i] = SignedFp12CarryModP(n, k, 4*n + LOGK3, p);
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                fdouble[i].in[l][j][idx] <== compress[i].out[l][j][idx]; 
            
            Pdouble[i] = EllipticCurveDoubleFp2(n, k, [0,0], b0, b1, p);  
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                Pdouble[i].in[j][l][idx] <== R[i+1][j][l][idx]; 
            
            if(Bits[i] == 0){
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    f[i][l][j][idx] <== fdouble[i].out[l][j][idx]; 
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                    R[i][j][l][idx] <== Pdouble[i].out[j][l][idx];
            }else{
                fadd[curid] = Fp12MultiplyWithLineUnequalFp2(n, k, k, n, p); 
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    fadd[curid].g[l][j][idx] <== fdouble[i].out[l][j][idx];
                
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++){
                    fadd[curid].P[0][j][l][idx] <== Pdouble[i].out[j][l][idx];            
                    fadd[curid].P[1][j][l][idx] <== P[j][l][idx];            
                }
                for(var eps=0; eps<2; eps++)for(var idx=0; idx<k; idx++)
                    fadd[curid].Q[eps][idx] <== Q[eps][idx];

                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    f[i][l][j][idx] <== fadd[curid].out[l][j][idx]; 

                // Padd[curid] = Pdouble[i] + P 
                Padd[curid] = EllipticCurveAddUnequalFp2(n, k, p); 
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++){
                    Padd[curid].a[j][l][idx] <== Pdouble[i].out[j][l][idx];
                    Padd[curid].b[j][l][idx] <== P[j][l][idx];
                }

                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                    R[i][j][l][idx] <== Padd[curid].out[j][l][idx];
                
                curid++;
            }
        }
    }
    // R[0] = f_loopCount(P, Q) 

    // coeff[1][j] = (9+u)^{(p-1)/6 * j}
    var coeff[12][6][2][20] = get_Fp12_frobenius(n, k);
    // Frob_p( twist(P) ) = ( (w^2 x)^p, (w^3 y)^p ) = twist( coeff[1][2] * x^p, coeff[1][3] * y^p ) 
    component frobP[2];
    component twistedFrobP[2];
    for(var i=0; i<2; i++){
        frobP[i] = Fp2Conjugate(n, k, p); 
        for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
            frobP[i].in[j][idx] <== P[i][j][idx];
        twistedFrobP[i] = Fp2Multiply(n, k, p); 
        for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            twistedFrobP[i].a[j][idx] <== coeff[1][2+i][j][idx];
            twistedFrobP[i].b[j][idx] <== frobP[i].out[j][idx];
        } 
    }
    
    // l_{[loopCount] P', Frob_p(P')}(Q) 
    component line1 = LineFunctionUnequalFp2(n, k, p);
    for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        line1.P[0][i][j][idx] <== R[0][i][j][idx]; 
        line1.P[1][i][j][idx] <== twistedFrobP[i].out[j][idx];
    }
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        line1.Q[i][idx] <== Q[i][idx];

    // l_{[loopCount] P' + Frob_p(P'), -Frob_p^2(P')}(Q) 
    component R_add = EllipticCurveAddUnequalFp2(n, k, p);
    for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        R_add.a[i][j][idx] <== R[0][i][j][idx]; 
        R_add.b[i][j][idx] <== twistedFrobP[i].out[j][idx];
    }
    // -Frob_p^2(P') = twist( coeff[1][2] * x1^p, coeff[1][3] * -y1^p ) where twistedFrobP.out = (x1, y1)
    // x1^p is conjugation 
    component negFrob2[2]; 
    negFrob2[0] = FpNegate(n, k, p);
    negFrob2[1] = FpNegate(n, k, p);
    for(var idx=0; idx<k; idx++){
        negFrob2[0].in[idx] <== twistedFrobP[0].out[1][idx];
        negFrob2[1].in[idx] <== twistedFrobP[1].out[0][idx];
    }

    component negTwistedFrob2P[2];
    for(var i=0; i<2; i++){
        negTwistedFrob2P[i] = Fp2Multiply(n, k, p);
        for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            negTwistedFrob2P[i].a[j][idx] <== coeff[1][2+i][j][idx];
            if( (i==0 && j==1) || (i==1 && j==0) )
                negTwistedFrob2P[i].b[j][idx] <== negFrob2[i].out[idx];
            else
                negTwistedFrob2P[i].b[j][idx] <== twistedFrobP[i].out[j][idx];
        } 
    }
    component line2 = LineFunctionUnequalFp2(n, k, p);
    for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        line2.P[0][i][j][idx] <== R_add.out[i][j][idx]; 
        line2.P[1][i][j][idx] <== negTwistedFrob2P[i].out[j][idx];
    }
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        line2.Q[i][idx] <== Q[i][idx];
    
    component res = Fp12MultiplyThree(n, k, p);
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        res.a[i][j][idx] <== f[0][i][j][idx];
        res.b[i][j][idx] <== line1.out[i][j][idx];
        res.c[i][j][idx] <== line2.out[i][j][idx];
    }
    
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[i][j][idx] <== res.out[i][j][idx];
}

// Assuming curve is of form Y^2 = X^3 + b for now (a = 0) for better register bounds 
// Parameters:
// - b = [b0, b1] is 2 x k array representing Fp2 element
// - pseudoBinaryEncoding is length `loopBitLength` array consisting of {-1, 0, 1} entries such that `loopCount = sum pseudoBinaryEncoding[i] * 2^i for i in range(loopBitLength)`
// Input, Output, Assumptions same as for MillerLoopFp2(n, k, b0, b1, loopCount, p) 
// The only difference is that allowing -1 in the encoding leads to less nonzero entries 
template OptimizedMillerLoopFp2(n, k, b0, b1, pseudoBinaryEncoding, loopBitLength, p){
    signal input P[2][2][k]; 
    signal input Q[2][k];

    signal output out[6][2][k];

    var LOGK = log_ceil(k);
    var XI0 = 9;
    var LOGK2 = log_ceil(36*(2+XI0)*(2+XI0) * k*k);
    var LOGK3 = log_ceil(36*(2+XI0)*(2+XI0) * k*k*(2*k-1));
    assert( 4*n + LOGK3 < 251 );
    
    var sigBits=0;
    for (var i = 0; i < loopBitLength; i++) {
        assert( pseudoBinaryEncoding[i] >= -1 && pseudoBinaryEncoding[i] <= 1 );
        if(pseudoBinaryEncoding[i] != 0)
            sigBits++;
    }

    signal negP[2][2][k];
    component negPy = Fp2Negate(n, k, p);
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        negPy.in[i][idx] <== P[1][i][idx]; 
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){
        negP[0][i][idx] <== P[0][i][idx];
        negP[1][i][idx] <== negPy.out[i][idx];
    }

    signal R[loopBitLength][2][2][k]; 
    signal f[loopBitLength][6][2][k];

    component Pdouble[loopBitLength];
    component fdouble[loopBitLength];
    component square[loopBitLength];
    component line[loopBitLength];
    component compress[loopBitLength];
    component nocarry[loopBitLength];
    component Padd[sigBits];
    component fadd[sigBits]; 
    var curid=0;

    for(var i=loopBitLength - 1; i>=0; i--){
        if( i == loopBitLength - 1 ){
            // f = 1 
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                if(l==0 && j==0 && idx==0)
                    f[i][l][j][idx] <== 1;
                else    
                    f[i][l][j][idx] <== 0;
            }
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                R[i][j][l][idx] <== pseudoBinaryEncoding[i] == 1 ? P[j][l][idx] : negP[j][l][idx];
        }else{
            // compute fdouble[i] = f[i+1]^2 * l_{R[i+1], R[i+1]}(Q) 
            square[i] = SignedFp12MultiplyNoCarry(n, k, 2*n + 4 + LOGK); // 6 x 2 x 2k-1 registers in [0, 6 * k * (2+XI0) * 2^{2n} )
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                square[i].a[l][j][idx] <== f[i+1][l][j][idx];
                square[i].b[l][j][idx] <== f[i+1][l][j][idx];
            }
            
            line[i] = LineFunctionEqualFp2(n, k, p); // 6 x 2 x k registers in [0, 2^n) 
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                line[i].P[j][l][idx] <== R[i+1][j][l][idx];            
            for(var eps=0; eps<2; eps++)
                for(var idx=0; idx<k; idx++)
                    line[i].Q[eps][idx] <== Q[eps][idx];
            
            nocarry[i] = SignedFp12MultiplyNoCarryUnequal(n, 2*k-1, k, 3*n + LOGK2); // 6 x 2 x 3k-2 registers < (6 * (2+XI0))^2 * k^2 * 2^{3n} ) 
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<2*k-1; idx++)
                nocarry[i].a[l][j][idx] <== square[i].out[l][j][idx];
            
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                nocarry[i].b[l][j][idx] <== line[i].out[l][j][idx];
            
            compress[i] = Fp12Compress(n, k, 2*k-2, p, 4*n + LOGK3); // 6 x 2 x k registers < (6 * (2+ XI0))^2 * k^2 * (2k-1) * 2^{4n} )
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<3*k-2; idx++)
                compress[i].in[l][j][idx] <== nocarry[i].out[l][j][idx];
            
            fdouble[i] = SignedFp12CarryModP(n, k, 4*n + LOGK3, p);
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                fdouble[i].in[l][j][idx] <== compress[i].out[l][j][idx]; 
            
            Pdouble[i] = EllipticCurveDoubleFp2(n, k, [0,0], b0, b1, p);  
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                Pdouble[i].in[j][l][idx] <== R[i+1][j][l][idx]; 
            
            if(pseudoBinaryEncoding[i] == 0){
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    f[i][l][j][idx] <== fdouble[i].out[l][j][idx]; 
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                    R[i][j][l][idx] <== Pdouble[i].out[j][l][idx];
            }else{
                fadd[curid] = Fp12MultiplyWithLineUnequalFp2(n, k, k, n, p); 
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    fadd[curid].g[l][j][idx] <== fdouble[i].out[l][j][idx];
                
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++){
                    fadd[curid].P[0][j][l][idx] <== Pdouble[i].out[j][l][idx];            
                    fadd[curid].P[1][j][l][idx] <== pseudoBinaryEncoding[i] == 1 ? P[j][l][idx] : negP[j][l][idx];
                }
                for(var eps=0; eps<2; eps++)for(var idx=0; idx<k; idx++)
                    fadd[curid].Q[eps][idx] <== Q[eps][idx];

                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    f[i][l][j][idx] <== fadd[curid].out[l][j][idx]; 

                // Padd[curid] = Pdouble[i] + P 
                Padd[curid] = EllipticCurveAddUnequalFp2(n, k, p); 
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++){
                    Padd[curid].a[j][l][idx] <== Pdouble[i].out[j][l][idx];
                    Padd[curid].b[j][l][idx] <== pseudoBinaryEncoding[i] == 1 ? P[j][l][idx] : negP[j][l][idx];
                }

                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                    R[i][j][l][idx] <== Padd[curid].out[j][l][idx];
                
                curid++;
            }
        }
    }
    // R[0] = f_loopCount(P, Q) 

    // coeff[1][j] = (9+u)^{(p-1)/6 * j}
    var coeff[12][6][2][20] = get_Fp12_frobenius(n, k);
    // Frob_p( twist(P) ) = ( (w^2 x)^p, (w^3 y)^p ) = twist( coeff[1][2] * x^p, coeff[1][3] * y^p ) 
    component frobP[2];
    component P1[2];
    for(var i=0; i<2; i++){
        frobP[i] = Fp2Conjugate(n, k, p); 
        for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
            frobP[i].in[j][idx] <== P[i][j][idx];
        P1[i] = Fp2Multiply(n, k, p); 
        for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            P1[i].a[j][idx] <== coeff[1][2+i][j][idx];
            P1[i].b[j][idx] <== frobP[i].out[j][idx];
        } 
    }
    
    // l_{[loopCount] P', Frob_p(P')}(Q) 
    component line1 = LineFunctionUnequalFp2(n, k, p);
    for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        line1.P[0][i][j][idx] <== R[0][i][j][idx]; 
        line1.P[1][i][j][idx] <== P1[i].out[j][idx];
    }
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        line1.Q[i][idx] <== Q[i][idx];

    // l_{[loopCount] P' + Frob_p(P'), -Frob_p^2(P')}(Q) 
    component R_add = EllipticCurveAddUnequalFp2(n, k, p);
    for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        R_add.a[i][j][idx] <== R[0][i][j][idx]; 
        R_add.b[i][j][idx] <== P1[i].out[j][idx];
    }
    // -Frob_p^2(P') = twist( coeff[1][2] * x1^p, coeff[1][3] * -y1^p ) where P1.out = (x1, y1)
    // x1^p is conjugation 
    component negFrob2[2]; 
    negFrob2[0] = FpNegate(n, k, p);
    negFrob2[1] = FpNegate(n, k, p);
    for(var idx=0; idx<k; idx++){
        negFrob2[0].in[idx] <== P1[0].out[1][idx];
        negFrob2[1].in[idx] <== P1[1].out[0][idx];
    }

    component negP2[2];
    for(var i=0; i<2; i++){
        negP2[i] = Fp2Multiply(n, k, p);
        for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            negP2[i].a[j][idx] <== coeff[1][2+i][j][idx];
            if( (i==0 && j==1) || (i==1 && j==0) )
                negP2[i].b[j][idx] <== negFrob2[i].out[idx];
            else
                negP2[i].b[j][idx] <== P1[i].out[j][idx];
        } 
    }
    component line2 = LineFunctionUnequalFp2(n, k, p);
    for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        line2.P[0][i][j][idx] <== R_add.out[i][j][idx]; 
        line2.P[1][i][j][idx] <== negP2[i].out[j][idx];
    }
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        line2.Q[i][idx] <== Q[i][idx];
    
    component res = Fp12MultiplyThree(n, k, p);
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        res.a[i][j][idx] <== f[0][i][j][idx];
        res.b[i][j][idx] <== line1.out[i][j][idx];
        res.c[i][j][idx] <== line2.out[i][j][idx];
    }
    
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[i][j][idx] <== res.out[i][j][idx];
}


// Assuming curve is of form Y^2 = X^3 + b for now (a = 0) for better register bounds 
// b = [b0, b1] is 2 x k array representing Fp2 element
// Let numPairing = 3
// Inputs:
//  P is numPairing x 2 x 2 x k array where P[i] = (x_i, y_i) is a point in E[r](Fp2) 
//  Q is numPairing x 2 x k array where Q[i] = (X_i, Y_i) is a point in E(Fp) 
// Output:
//  out = e'(P[0], Q[0]) * ... * e'(P[numPairing-1], Q[numPairing-1]) 
//  - where e'(P,Q) is the output of `MillerLoopFp2`
// Assume:
//  r is prime (not a parameter in this template)
//  loopCount in [0, 2^250) and loopCount < r and loopCount < p
//  p has k registers in [0, 2^n)
//  P != O so the order of P in E(Fp) is r, so [i]P != [j]P for i != j in Z/r 
//  X^3 + b = 0 has no solution in Fp2, i.e., the y-coordinate of P cannot be 0.
template MillerLoopThreeFp2(n, k, b0, b1, loopCount, p){
    var numPairing = 3;
    signal input P[numPairing][2][2][k]; 
    signal input Q[numPairing][2][k];

    signal output out[6][2][k];

    var LOGK = log_ceil(k);
    var XI0 = 9;
    var LOGK2 = log_ceil(36*(2+XI0)*(2+XI0) * k*k);
    var LOGK3 = log_ceil(36*(2+XI0)*(2+XI0) * k*k*(2*k-1));
    assert( 4*n + LOGK3 < 251 );
    
    var Bits[250]; 
    var BitLength;
    var SigBits=0;
    for (var i = 0; i < 250; i++) {
        Bits[i] = (loopCount >> i) & 1;
        if(Bits[i] == 1){
            SigBits++;
            BitLength = i + 1;
        }
    }

    signal R[BitLength][numPairing][2][2][k]; 
    signal f[BitLength][6][2][k];

    component Pdouble[BitLength][numPairing];
    component fdouble[BitLength];
    component line[BitLength][numPairing];
    component lineProduct[BitLength];
    component Padd[SigBits][numPairing];
    component fadd[SigBits][numPairing]; 
    var curid=0;

    for(var i=BitLength - 1; i>=0; i--){
        if( i == BitLength - 1 ){
            // f = 1 
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                if(l==0 && j==0 && idx==0)
                    f[i][l][j][idx] <== 1;
                else    
                    f[i][l][j][idx] <== 0;
            }
            for(var idP=0; idP<numPairing; idP++)
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                    R[i][idP][j][l][idx] <== P[idP][j][l][idx];
        }else{
            // compute fdouble[i] = f[i+1]^2 * \prod_{idP=0}^{numPairing-1} l_{R[idP][i+1], R[idP][i+1]}(Q) 
            for(var idP=0; idP<numPairing; idP++){
                line[i][idP] = LineFunctionEqualFp2(n, k, p); // 6 x 2 x k registers in [0, 2^n) 
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                    line[i][idP].P[j][l][idx] <== R[i+1][idP][j][l][idx];            
                for(var l=0; l<2; l++)for(var idx=0; idx<k; idx++)
                    line[i][idP].Q[l][idx] <== Q[idP][l][idx];
                
                Pdouble[i][idP] = EllipticCurveDoubleFp2(n, k, [0,0], b0, b1, p);  
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                    Pdouble[i][idP].in[j][l][idx] <== R[i+1][idP][j][l][idx]; 
            }
            
            // Specializing to numPairing = 3: \prod_{idP=0}^{numPairing-1} l_{R[i+1][idP], R[i+1][idP]}(Q[idP]) 
            lineProduct[i] = Fp12MultiplyThree(n, k, p);
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                lineProduct[i].a[l][j][idx] <== line[i][0].out[l][j][idx];
                lineProduct[i].b[l][j][idx] <== line[i][1].out[l][j][idx];
                lineProduct[i].c[l][j][idx] <== line[i][2].out[l][j][idx];
            }
            
            fdouble[i] = Fp12MultiplyThree(n, k, p);
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                fdouble[i].a[l][j][idx] <== f[i+1][l][j][idx];
                fdouble[i].b[l][j][idx] <== f[i+1][l][j][idx];
                fdouble[i].c[l][j][idx] <== lineProduct[i].out[l][j][idx];
            }
            
            if(Bits[i] == 0){
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    f[i][l][j][idx] <== fdouble[i].out[l][j][idx]; 
                for(var idP=0; idP<numPairing; idP++)
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                        R[i][idP][j][l][idx] <== Pdouble[i][idP].out[j][l][idx];
            }else{
                // fdouble[i] * \prod_{idP=0}^{numPairing-1} l_{Pdouble[i][idP], P[idP]}(Q[idP])
                for(var idP=0; idP<numPairing; idP++){
                    fadd[curid][idP] = Fp12MultiplyWithLineUnequalFp2(n, k, k, n, p); 
                    for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                        fadd[curid][idP].g[l][j][idx] <== idP == 0 ? fdouble[i].out[l][j][idx] : fadd[curid][idP-1].out[l][j][idx];
                    }
                
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++){
                        fadd[curid][idP].P[0][j][l][idx] <== Pdouble[i][idP].out[j][l][idx]; 
                        fadd[curid][idP].P[1][j][l][idx] <== P[idP][j][l][idx];            
                    }
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                        fadd[curid][idP].Q[j][idx] <== Q[idP][j][idx];

                    // Padd[curid][idP] = Pdouble[i][idP] + P[idP] 
                    Padd[curid][idP] = EllipticCurveAddUnequalFp2(n, k, p); 
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++){
                        Padd[curid][idP].a[j][l][idx] <== Pdouble[i][idP].out[j][l][idx];
                        Padd[curid][idP].b[j][l][idx] <== P[idP][j][l][idx];
                    }

                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                        R[i][idP][j][l][idx] <== Padd[curid][idP].out[j][l][idx];
                
                }
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    f[i][l][j][idx] <== fadd[curid][numPairing-1].out[l][j][idx]; 
                
                curid++;
            }
        }
    }
    // R[0] = f_loopCount(P, Q) 

    // coeff[1][j] = (9+u)^{(p-1)/6 * j}
    var coeff[12][6][2][20] = get_Fp12_frobenius(n, k);
    // Frob_p( twist(P) ) = ( (w^2 x)^p, (w^3 y)^p ) = twist( coeff[1][2] * x^p, coeff[1][3] * y^p ) 
    component frobP[numPairing][2];
    component twistedFrobP[numPairing][2];
    component line1[numPairing];
    component line2[numPairing];
    component R_add[numPairing];
    component negFrob2[numPairing][2]; 
    component negTwistedFrob2P[numPairing][2];
    component outProd[numPairing];
    
    for(var idP=0; idP<numPairing; idP++){
        for(var i=0; i<2; i++){
            frobP[idP][i] = Fp2Conjugate(n, k, p); 
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                frobP[idP][i].in[j][idx] <== P[idP][i][j][idx];

            twistedFrobP[idP][i] = Fp2Multiply(n, k, p); 
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                twistedFrobP[idP][i].a[j][idx] <== coeff[1][2+i][j][idx];
                twistedFrobP[idP][i].b[j][idx] <== frobP[idP][i].out[j][idx];
            } 
        }
        
        // l_{[loopCount] P', Frob_p(P')}(Q) 
        line1[idP] = LineFunctionUnequalFp2(n, k, p);
        for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            line1[idP].P[0][i][j][idx] <== R[0][idP][i][j][idx]; 
            line1[idP].P[1][i][j][idx] <== twistedFrobP[idP][i].out[j][idx];
        }
        for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
            line1[idP].Q[i][idx] <== Q[idP][i][idx];

        // l_{[loopCount] P' + Frob_p(P'), -Frob_p^2(P')}(Q) 
        R_add[idP] = EllipticCurveAddUnequalFp2(n, k, p);
        for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            R_add[idP].a[i][j][idx] <== R[0][idP][i][j][idx]; 
            R_add[idP].b[i][j][idx] <== twistedFrobP[idP][i].out[j][idx];
        }
        // -Frob_p^2(P') = twist( coeff[1][2] * x1^p, coeff[1][3] * -y1^p ) where twistedFrobP.out = (x1, y1)
        // x1^p is conjugation 
        negFrob2[idP][0] = FpNegate(n, k, p);
        negFrob2[idP][1] = FpNegate(n, k, p);
        for(var idx=0; idx<k; idx++){
            negFrob2[idP][0].in[idx] <== twistedFrobP[idP][0].out[1][idx];
            negFrob2[idP][1].in[idx] <== twistedFrobP[idP][1].out[0][idx];
        }

        for(var i=0; i<2; i++){
            negTwistedFrob2P[idP][i] = Fp2Multiply(n, k, p);
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                negTwistedFrob2P[idP][i].a[j][idx] <== coeff[1][2+i][j][idx];
                if( (i==0 && j==1) || (i==1 && j==0) )
                    negTwistedFrob2P[idP][i].b[j][idx] <== negFrob2[idP][i].out[idx];
                else
                    negTwistedFrob2P[idP][i].b[j][idx] <== twistedFrobP[idP][i].out[j][idx];
            } 
        }
        line2[idP] = LineFunctionUnequalFp2(n, k, p);
        for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            line2[idP].P[0][i][j][idx] <== R_add[idP].out[i][j][idx]; 
            line2[idP].P[1][i][j][idx] <== negTwistedFrob2P[idP][i].out[j][idx];
        }
        for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
            line2[idP].Q[i][idx] <== Q[idP][i][idx];
        
        outProd[idP] = Fp12MultiplyThree(n, k, p);
        for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            outProd[idP].a[i][j][idx] <== idP == 0 ? f[0][i][j][idx] : outProd[idP-1].out[i][j][idx];
            outProd[idP].b[i][j][idx] <== line1[idP].out[i][j][idx];
            outProd[idP].c[i][j][idx] <== line2[idP].out[i][j][idx];
        }
    } 
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[i][j][idx] <== outProd[numPairing-1].out[i][j][idx];
}

// Assuming curve is of form Y^2 = X^3 + b for now (a = 0) for better register bounds 
// Parameters:
// - b = [b0, b1] is 2 x k array representing Fp2 element
// - pseudoBinaryEncoding is length `loopBitLength` array consisting of {-1, 0, 1} entries such that `loopCount = sum pseudoBinaryEncoding[i] * 2^i for i in range(loopBitLength)`
// Let numPairing = 3
// Input, Output, Assumptions same as for MillerLoopThreeFp2(n, k, b0, b1, loopCount, p) 
// The only difference is that allowing -1 in the encoding leads to less nonzero entries 
template OptimizedMillerLoopThreeFp2(n, k, b0, b1, pseudoBinaryEncoding, loopBitLength, p){
    var numPairing = 3;
    signal input P[numPairing][2][2][k]; 
    signal input Q[numPairing][2][k];

    signal output out[6][2][k];

    var XI0 = 9;
    var LOGK3 = log_ceil(36*(2+XI0)*(2+XI0) * k*k*(2*k-1));
    assert( 4*n + LOGK3 < 251 );
    
    var sigBits=0;
    for (var i = 0; i < loopBitLength; i++) {
        assert( pseudoBinaryEncoding[i] >= -1 && pseudoBinaryEncoding[i] <= 1 );
        if(pseudoBinaryEncoding[i] != 0)
            sigBits++;
    }

    signal negP[numPairing][2][2][k];
    component negPy[numPairing];
    for(var idP=0; idP<numPairing; idP++){
        negPy[idP] = Fp2Negate(n, k, p);
        for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
            negPy[idP].in[i][idx] <== P[idP][1][i][idx]; 
        for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){
            negP[idP][0][i][idx] <== P[idP][0][i][idx];
            negP[idP][1][i][idx] <== negPy[idP].out[i][idx];
        }
    }

    signal R[loopBitLength][numPairing][2][2][k]; 
    signal f[loopBitLength][6][2][k];

    component Pdouble[loopBitLength][numPairing];
    component fdouble[loopBitLength];
    component line[loopBitLength][numPairing];
    component lineProduct[loopBitLength];
    component Padd[sigBits][numPairing];
    component fadd[sigBits][numPairing]; 
    var curid=0;

    for(var i=loopBitLength - 1; i>=0; i--){
        if( i == loopBitLength - 1 ){
            // f = 1 
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                if(l==0 && j==0 && idx==0)
                    f[i][l][j][idx] <== 1;
                else    
                    f[i][l][j][idx] <== 0;
            }
            for(var idP=0; idP<numPairing; idP++)
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                    R[i][idP][j][l][idx] <== pseudoBinaryEncoding[i] == 1 ? P[idP][j][l][idx] : negP[idP][j][l][idx];
        }else{
            // compute fdouble[i] = f[i+1]^2 * \prod_{idP=0}^{numPairing-1} l_{R[idP][i+1], R[idP][i+1]}(Q) 
            for(var idP=0; idP<numPairing; idP++){
                line[i][idP] = LineFunctionEqualFp2(n, k, p); // 6 x 2 x k registers in [0, 2^n) 
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                    line[i][idP].P[j][l][idx] <== R[i+1][idP][j][l][idx];            
                for(var l=0; l<2; l++)for(var idx=0; idx<k; idx++)
                    line[i][idP].Q[l][idx] <== Q[idP][l][idx];
                
                Pdouble[i][idP] = EllipticCurveDoubleFp2(n, k, [0,0], b0, b1, p);  
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                    Pdouble[i][idP].in[j][l][idx] <== R[i+1][idP][j][l][idx]; 
            }
            
            // Specializing to numPairing = 3: \prod_{idP=0}^{numPairing-1} l_{R[i+1][idP], R[i+1][idP]}(Q[idP]) 
            lineProduct[i] = Fp12MultiplyThree(n, k, p);
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                lineProduct[i].a[l][j][idx] <== line[i][0].out[l][j][idx];
                lineProduct[i].b[l][j][idx] <== line[i][1].out[l][j][idx];
                lineProduct[i].c[l][j][idx] <== line[i][2].out[l][j][idx];
            }
            
            fdouble[i] = Fp12MultiplyThree(n, k, p);
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                fdouble[i].a[l][j][idx] <== f[i+1][l][j][idx];
                fdouble[i].b[l][j][idx] <== f[i+1][l][j][idx];
                fdouble[i].c[l][j][idx] <== lineProduct[i].out[l][j][idx];
            }
            
            if(pseudoBinaryEncoding[i] == 0){
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    f[i][l][j][idx] <== fdouble[i].out[l][j][idx]; 
                for(var idP=0; idP<numPairing; idP++)
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                        R[i][idP][j][l][idx] <== Pdouble[i][idP].out[j][l][idx];
            }else{
                // fdouble[i] * \prod_{idP=0}^{numPairing-1} l_{Pdouble[i][idP], P[idP]}(Q[idP])
                for(var idP=0; idP<numPairing; idP++){
                    fadd[curid][idP] = Fp12MultiplyWithLineUnequalFp2(n, k, k, n, p); 
                    for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                        fadd[curid][idP].g[l][j][idx] <== idP == 0 ? fdouble[i].out[l][j][idx] : fadd[curid][idP-1].out[l][j][idx];
                    }
                
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++){
                        fadd[curid][idP].P[0][j][l][idx] <== Pdouble[i][idP].out[j][l][idx]; 
                        fadd[curid][idP].P[1][j][l][idx] <== pseudoBinaryEncoding[i] == 1 ? P[idP][j][l][idx] : negP[idP][j][l][idx];
                    }
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                        fadd[curid][idP].Q[j][idx] <== Q[idP][j][idx];

                    // Padd[curid][idP] = Pdouble[i][idP] + P[idP] 
                    Padd[curid][idP] = EllipticCurveAddUnequalFp2(n, k, p); 
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++){
                        Padd[curid][idP].a[j][l][idx] <== Pdouble[i][idP].out[j][l][idx];
                        Padd[curid][idP].b[j][l][idx] <== pseudoBinaryEncoding[i] == 1 ? P[idP][j][l][idx] : negP[idP][j][l][idx];
                    }

                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                        R[i][idP][j][l][idx] <== Padd[curid][idP].out[j][l][idx];
                
                }
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    f[i][l][j][idx] <== fadd[curid][numPairing-1].out[l][j][idx]; 
                
                curid++;
            }
        }
    }
    // R[0] = f_loopCount(P, Q) 

    // coeff[1][j] = (9+u)^{(p-1)/6 * j}
    var coeff[12][6][2][20] = get_Fp12_frobenius(n, k);
    // Frob_p( twist(P) ) = ( (w^2 x)^p, (w^3 y)^p ) = twist( coeff[1][2] * x^p, coeff[1][3] * y^p ) 
    component frobP[numPairing][2];
    component P1[numPairing][2];
    component line1[numPairing];
    component line2[numPairing];
    component R_add[numPairing];
    component negFrob2[numPairing][2]; 
    component negP2[numPairing][2];
    component outProd[numPairing];
    
    for(var idP=0; idP<numPairing; idP++){
        for(var i=0; i<2; i++){
            frobP[idP][i] = Fp2Conjugate(n, k, p); 
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                frobP[idP][i].in[j][idx] <== P[idP][i][j][idx];

            P1[idP][i] = Fp2Multiply(n, k, p); 
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                P1[idP][i].a[j][idx] <== coeff[1][2+i][j][idx];
                P1[idP][i].b[j][idx] <== frobP[idP][i].out[j][idx];
            } 
        }
        
        // l_{[loopCount] P', Frob_p(P')}(Q) 
        line1[idP] = LineFunctionUnequalFp2(n, k, p);
        for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            line1[idP].P[0][i][j][idx] <== R[0][idP][i][j][idx]; 
            line1[idP].P[1][i][j][idx] <== P1[idP][i].out[j][idx];
        }
        for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
            line1[idP].Q[i][idx] <== Q[idP][i][idx];

        // l_{[loopCount] P' + Frob_p(P'), -Frob_p^2(P')}(Q) 
        R_add[idP] = EllipticCurveAddUnequalFp2(n, k, p);
        for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            R_add[idP].a[i][j][idx] <== R[0][idP][i][j][idx]; 
            R_add[idP].b[i][j][idx] <== P1[idP][i].out[j][idx];
        }
        // -Frob_p^2(P') = twist( coeff[1][2] * x1^p, coeff[1][3] * -y1^p ) where P1.out = (x1, y1)
        // x1^p is conjugation 
        negFrob2[idP][0] = FpNegate(n, k, p);
        negFrob2[idP][1] = FpNegate(n, k, p);
        for(var idx=0; idx<k; idx++){
            negFrob2[idP][0].in[idx] <== P1[idP][0].out[1][idx];
            negFrob2[idP][1].in[idx] <== P1[idP][1].out[0][idx];
        }

        for(var i=0; i<2; i++){
            negP2[idP][i] = Fp2Multiply(n, k, p);
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                negP2[idP][i].a[j][idx] <== coeff[1][2+i][j][idx];
                if( (i==0 && j==1) || (i==1 && j==0) )
                    negP2[idP][i].b[j][idx] <== negFrob2[idP][i].out[idx];
                else
                    negP2[idP][i].b[j][idx] <== P1[idP][i].out[j][idx];
            } 
        }
        line2[idP] = LineFunctionUnequalFp2(n, k, p);
        for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            line2[idP].P[0][i][j][idx] <== R_add[idP].out[i][j][idx]; 
            line2[idP].P[1][i][j][idx] <== negP2[idP][i].out[j][idx];
        }
        for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
            line2[idP].Q[i][idx] <== Q[idP][i][idx];
        
        outProd[idP] = Fp12MultiplyThree(n, k, p);
        for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            outProd[idP].a[i][j][idx] <== idP == 0 ? f[0][i][j][idx] : outProd[idP-1].out[i][j][idx];
            outProd[idP].b[i][j][idx] <== line1[idP].out[i][j][idx];
            outProd[idP].c[i][j][idx] <== line2[idP].out[i][j][idx];
        }
    } 
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[i][j][idx] <== outProd[numPairing-1].out[i][j][idx];
}

// Assuming curve is of form Y^2 = X^3 + b for now (a = 0) for better register bounds 
// Parameters:
// - b = [b0, b1] is 2 x k array representing Fp2 element
// - pseudoBinaryEncoding is length `loopBitLength` array consisting of {-1, 0, 1} entries such that `loopCount = sum pseudoBinaryEncoding[i] * 2^i for i in range(loopBitLength)`
// - numPairing is a small number 
// - p is prime, has k registers in [0, 2^n)
// Inputs:
//  P is numPairing x 2 x 2 x k array where P[i] = (x_i, y_i) is a point in E[r](Fp2) 
//  Q is numPairing x 2 x k array where Q[i] = (X_i, Y_i) is a point in E(Fp) 
// Output:
//  out = e'(P[0], Q[0]) * ... * e'(P[numPairing-1], Q[numPairing-1]) 
//  - where e'(P,Q) is the output of `MillerLoopFp2`
// Assume:
//  r is prime (not a parameter in this template)
//  loopCount in [0, 2^250) and loopCount < r and loopCount < p
//  P != O so the order of P in E(Fp) is r, so [i]P != [j]P for i != j in Z/r 
//  X^3 + b = 0 has no solution in Fp2, i.e., the y-coordinate of P cannot be 0.
template OptimizedMillerLoopProductFp2(n, k, b0, b1, numPairing, pseudoBinaryEncoding, loopBitLength, p){
    signal input P[numPairing][2][2][k]; 
    signal input Q[numPairing][2][k];

    signal output out[6][2][k];

    var XI0 = 9;
    var LOGK4 = log_ceil(12*k*k * (XI0 + 1) * 6 * k * (2+XI0) * (3*k-2) + 1);
    assert( 5*n + LOGK4 < 251 );
    
    var sigBits=0;
    for (var i = 0; i < loopBitLength; i++) {
        assert( pseudoBinaryEncoding[i] >= -1 && pseudoBinaryEncoding[i] <= 1 );
        if(pseudoBinaryEncoding[i] != 0)
            sigBits++;
    }

    signal negP[numPairing][2][2][k];
    component negPy[numPairing];
    for(var idP=0; idP<numPairing; idP++){
        negPy[idP] = Fp2Negate(n, k, p);
        for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
            negPy[idP].in[i][idx] <== P[idP][1][i][idx]; 
        for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){
            negP[idP][0][i][idx] <== P[idP][0][i][idx];
            negP[idP][1][i][idx] <== negPy[idP].out[i][idx];
        }
    }

    signal R[loopBitLength][numPairing][2][2][k]; 
    signal f[loopBitLength][6][2][k];

    component Pdouble[loopBitLength][numPairing];
    component fSquare[loopBitLength];
    component lineProduct[loopBitLength][numPairing];
    component Padd[sigBits][numPairing];
    component fadd[sigBits][numPairing]; 
    var curid=0;

    for(var i=loopBitLength - 1; i>=0; i--){
        if( i == loopBitLength - 1 ){
            // f = 1 
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                if(l==0 && j==0 && idx==0)
                    f[i][l][j][idx] <== 1;
                else    
                    f[i][l][j][idx] <== 0;
            }
            for(var idP=0; idP<numPairing; idP++)
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                    R[i][idP][j][l][idx] <== pseudoBinaryEncoding[i] == 1 ? P[idP][j][l][idx] : negP[idP][j][l][idx];
        }else{
            // compute fdouble[i] = f[i+1]^2 * \prod_{idP=0}^{numPairing-1} l_{R[idP][i+1], R[idP][i+1]}(Q) 
            fSquare[i] = Fp12Multiply(n, k, p);
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                fSquare[i].a[l][j][idx] <== f[i+1][l][j][idx];
                fSquare[i].b[l][j][idx] <== f[i+1][l][j][idx];
            }
            
            for(var idP=0; idP<numPairing; idP++){
                lineProduct[i][idP] = Fp12MultiplyWithLineEqualFp2(n, k, p); // 6 x 2 x k registers in [0, 2^n) 
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    lineProduct[i][idP].g[l][j][idx] <== idP == 0 ? fSquare[i].out[l][j][idx] : lineProduct[i][idP-1].out[l][j][idx];
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                    lineProduct[i][idP].P[j][l][idx] <== R[i+1][idP][j][l][idx];
                for(var l=0; l<2; l++)for(var idx=0; idx<k; idx++)
                    lineProduct[i][idP].Q[l][idx] <== Q[idP][l][idx];
                
                Pdouble[i][idP] = EllipticCurveDoubleFp2(n, k, [0,0], b0, b1, p);  
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                    Pdouble[i][idP].in[j][l][idx] <== R[i+1][idP][j][l][idx]; 
            }
            // fdouble[i] = lineProduct[i][numPairing-1]
            
            if(pseudoBinaryEncoding[i] == 0){
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    f[i][l][j][idx] <== lineProduct[i][numPairing-1].out[l][j][idx]; 
                for(var idP=0; idP<numPairing; idP++)
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                        R[i][idP][j][l][idx] <== Pdouble[i][idP].out[j][l][idx];
            }else{
                // fdouble[i] * \prod_{idP=0}^{numPairing-1} l_{Pdouble[i][idP], P[idP]}(Q[idP])
                for(var idP=0; idP<numPairing; idP++){
                    fadd[curid][idP] = Fp12MultiplyWithLineUnequalFp2(n, k, k, n, p); 
                    for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                        fadd[curid][idP].g[l][j][idx] <== idP == 0 ? lineProduct[i][numPairing-1].out[l][j][idx] : fadd[curid][idP-1].out[l][j][idx];
                    }
                
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++){
                        fadd[curid][idP].P[0][j][l][idx] <== Pdouble[i][idP].out[j][l][idx]; 
                        fadd[curid][idP].P[1][j][l][idx] <== pseudoBinaryEncoding[i] == 1 ? P[idP][j][l][idx] : negP[idP][j][l][idx];
                    }
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                        fadd[curid][idP].Q[j][idx] <== Q[idP][j][idx];

                    // Padd[curid][idP] = Pdouble[i][idP] + P[idP] 
                    Padd[curid][idP] = EllipticCurveAddUnequalFp2(n, k, p); 
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++){
                        Padd[curid][idP].a[j][l][idx] <== Pdouble[i][idP].out[j][l][idx];
                        Padd[curid][idP].b[j][l][idx] <== pseudoBinaryEncoding[i] == 1 ? P[idP][j][l][idx] : negP[idP][j][l][idx];
                    }

                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)for(var l=0; l<2; l++)
                        R[i][idP][j][l][idx] <== Padd[curid][idP].out[j][l][idx];
                
                }
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    f[i][l][j][idx] <== fadd[curid][numPairing-1].out[l][j][idx]; 
                
                curid++;
            }
        }
    }
    // R[0] = f_loopCount(P, Q) 

    // coeff[1][j] = (9+u)^{(p-1)/6 * j}
    var coeff[12][6][2][20] = get_Fp12_frobenius(n, k);
    // Frob_p( twist(P) ) = ( (w^2 x)^p, (w^3 y)^p ) = twist( coeff[1][2] * x^p, coeff[1][3] * y^p ) 
    component frobP[numPairing][2];
    component P1[numPairing][2];
    component line1[numPairing];
    component line2[numPairing];
    component R_add[numPairing];
    component negFrob2[numPairing][2]; 
    component negP2[numPairing][2];
    component outProd[numPairing];
    
    for(var idP=0; idP<numPairing; idP++){
        for(var i=0; i<2; i++){
            frobP[idP][i] = Fp2Conjugate(n, k, p); 
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                frobP[idP][i].in[j][idx] <== P[idP][i][j][idx];

            P1[idP][i] = Fp2Multiply(n, k, p); 
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                P1[idP][i].a[j][idx] <== coeff[1][2+i][j][idx];
                P1[idP][i].b[j][idx] <== frobP[idP][i].out[j][idx];
            } 
        }
        
        // l_{[loopCount] P', Frob_p(P')}(Q) 
        line1[idP] = LineFunctionUnequalFp2(n, k, p);
        for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            line1[idP].P[0][i][j][idx] <== R[0][idP][i][j][idx]; 
            line1[idP].P[1][i][j][idx] <== P1[idP][i].out[j][idx];
        }
        for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
            line1[idP].Q[i][idx] <== Q[idP][i][idx];

        // l_{[loopCount] P' + Frob_p(P'), -Frob_p^2(P')}(Q) 
        R_add[idP] = EllipticCurveAddUnequalFp2(n, k, p);
        for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            R_add[idP].a[i][j][idx] <== R[0][idP][i][j][idx]; 
            R_add[idP].b[i][j][idx] <== P1[idP][i].out[j][idx];
        }
        // -Frob_p^2(P') = twist( coeff[1][2] * x1^p, coeff[1][3] * -y1^p ) where P1.out = (x1, y1)
        // x1^p is conjugation 
        negFrob2[idP][0] = FpNegate(n, k, p);
        negFrob2[idP][1] = FpNegate(n, k, p);
        for(var idx=0; idx<k; idx++){
            negFrob2[idP][0].in[idx] <== P1[idP][0].out[1][idx];
            negFrob2[idP][1].in[idx] <== P1[idP][1].out[0][idx];
        }

        for(var i=0; i<2; i++){
            negP2[idP][i] = Fp2Multiply(n, k, p);
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                negP2[idP][i].a[j][idx] <== coeff[1][2+i][j][idx];
                if( (i==0 && j==1) || (i==1 && j==0) )
                    negP2[idP][i].b[j][idx] <== negFrob2[idP][i].out[idx];
                else
                    negP2[idP][i].b[j][idx] <== P1[idP][i].out[j][idx];
            } 
        }
        line2[idP] = LineFunctionUnequalFp2(n, k, p);
        for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            line2[idP].P[0][i][j][idx] <== R_add[idP].out[i][j][idx]; 
            line2[idP].P[1][i][j][idx] <== negP2[idP][i].out[j][idx];
        }
        for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
            line2[idP].Q[i][idx] <== Q[idP][i][idx];
        
        outProd[idP] = Fp12MultiplyThree(n, k, p);
        for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
            outProd[idP].a[i][j][idx] <== idP == 0 ? f[0][i][j][idx] : outProd[idP-1].out[i][j][idx];
            outProd[idP].b[i][j][idx] <== line1[idP].out[i][j][idx];
            outProd[idP].c[i][j][idx] <== line2[idP].out[i][j][idx];
        }
    } 
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[i][j][idx] <== outProd[numPairing-1].out[i][j][idx];
}


template OptimalAtePairing(n, k, p){
    signal input P[2][2][k]; 
    signal input Q[2][k];
    signal output out[6][2][k];

    //var loop_count = get_bn254_ate_loop_count();
    var bn254_b[2][50] = get_bn254_b(n, k);

    component miller = OptimizedMillerLoopFp2(n, k, bn254_b[0], bn254_b[1], get_bn254_pseudo_binary_encoding(), 65, p);
    for(var i=0; i<2; i++)for(var j=0; j <2; j++)for(var idx=0; idx<k; idx++)
        miller.P[i][j][idx] <== P[i][j][idx];
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        miller.Q[i][idx] <== Q[i][idx];
    
    component finalexp = FinalExponentiate(n, k, p);
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        finalexp.in[i][j][idx] <== miller.out[i][j][idx];

    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[i][j][idx] <== finalexp.out[i][j][idx]; 

    /*
    // constrain out[i][j] < p since we never did that previously
    component lt[6][2];
    for(var i=0; i<6; i++)for(var j=0; j<2; j++){
        lt[i][j] = BigLessThan(n, k);
        for(var idx=0; idx<k; idx++){
            lt[i][j].a[idx] <== out[i][j][idx];
            lt[i][j].b[idx] <== q[idx];
        }
        lt[i][j].out === 1;
    }*/
}

