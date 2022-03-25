pragma circom 2.0.2;

include "../node_modules/circomlib/circuits/bitify.circom";

include "./bigint.circom";
include "./bigint_func.circom";
include "./fp.circom";
include "./fp2.circom";
include "./fp12.circom";
include "./curve.circom";

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

    var q[20];
    q[0] = q0;
    q[1] = q1;
    q[2] = q2;
    for (var idx = 3; idx < 20; idx++) {
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
    var sub0inv[20] = mod_inv(n, k, sub0.out, q);
    var sub1_sub0inv[20] = prod(n, k, sub1.out, sub0inv);
    var lamb_arr[2][20] = long_div(n, k, sub1_sub0inv, q);
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

    var q[20];
    q[0] = q0;
    q[1] = q1;
    q[2] = q2;
    q[3] = q3;
    for (var idx = 4; idx < 20; idx++) {
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
    var sub0inv[20] = mod_inv(n, k, sub0.out, q);
    var sub1_sub0inv[20] = prod(n, k, sub1.out, sub0inv);
    var lamb_arr[2][20] = long_div(n, k, sub1_sub0inv, q);
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
    var q[20];
    q[0] = q0;
    q[1] = q1;
    q[2] = q2;
    q[3] = q3;
    for (var idx = 4; idx < 20; idx++) {
	    q[idx] = 0;
    }

    // assuming a is small 
    var long_a[20];
    long_a[0] = a;
    for (var i = 1; i < 20; i++) {
        long_a[i] = 0;   
    }

    component in0_sq = BigMultModP(n, k);
    for (var i = 0; i < k; i++) {
        in0_sq.a[i] <== in[0][i];
        in0_sq.b[i] <== in[0][i];
        in0_sq.p[i] <== q[i];
    }

    var long_2[20];
    var long_3[20];
    long_2[0] = 2;
    long_3[0] = 3;
    for (var i = 1; i < k; i++) {
        long_a[i] = 0;
        long_2[i] = 0;
        long_3[i] = 0;
    }
    var inv_2[20] = mod_inv(n, k, long_2, q);
    var long_3_div_2[20] = prod(n, k, long_3, inv_2);
    var long_3_div_2_mod_q[2][20] = long_div(n, k, long_3_div_2, q);

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
    var denom_inv[20] = mod_inv(n, k, in[1], q);
    var product[20] = prod(n, k, numer.out, denom_inv);
    var lamb_arr[2][20] = long_div(n, k, product, q);
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

// Assuming curve is of form Y^2 = X^3 + b for now (a = 0) for better register bounds 
// Inputs:
//  P is 2 x k array where P = (x, y) is a point in E[r](Fq) 
//  Q is 2 x 6 x 2 x k array representing point (X, Y) in E(Fq12) 
// Output:
// f_r(Q) where <f_r> = [r]P - [r]O is computed using Miller's algorithm
// Assume:
//  r has k registers in [0, 2^n)
//  q has k registers in [0, 2^n)
//  r is prime
//  P != O so the order of P in E(Fq) is r, so [i]P != [j]P for i != j in Z/r 
template MillerLoop1(n, k, b, r, q){
    signal input P[2][k]; 
    signal input Q[2][6][2][k];
    signal output out[6][2][k];

    var rBits[500]; // length is k * n
    var rBitLength;
    var rSigBits=0;
    for (var i = 0; i < k; i++) {
        for (var j = 0; j < n; j++) {
            rBits[j + n * i] = (r[i] >> j) & 1;
            if(rBits[j + n * i] == 1){
                rSigBits++;
                rBitLength = j + n * i + 1;
            }
        }
    }

    signal Pintermed[rBitLength][2][k]; 
    signal f[rBitLength][6][2][k];

    component Pdouble[rBitLength];
    component fdouble[rBitLength];
    component square[rBitLength];
    component line[rBitLength];
    component compress[rBitLength];
    component nocarry[rBitLength];
    component Padd[rSigBits];
    component fadd[rSigBits]; 
    var curid=0;

    var LOGK = log_ceil(k);
    var logc = log_ceil(144*k*k*k);
    assert( 4*n + logc + 1 < 252 );
    
    for(var i=rBitLength - 1; i>=0; i--){
        if( i == rBitLength - 1 ){
            // f = 1 
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                if(l==0 && j==0 && idx==0)
                    f[i][l][j][idx] <== 1;
                else
                    f[i][l][j][idx] <== 0;
            }
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                Pintermed[i][j][idx] <== P[j][idx];
        }else{
            // compute fdouble[i] = f[i+1]^2 * l_{Pintermed[i+1], Pintermed[i+1]}(Q) 
            square[i] = Fp12MultiplyNoCarry(n, k); // 6 x 4 x 2k-1 registers in [0, 12 * k * 2^{2n} )
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

            nocarry[i] = Fp12MultiplyNoCarryUnequal(n, 2*k-1, k); // 6 x 4 x 3k-2 registers in [0, 144 * k^2 * 2^{3n} ) 
            for(var l=0; l<6; l++)for(var j=0; j<4; j++)for(var idx=0; idx<2*k-1; idx++)
                nocarry[i].a[l][j][idx] <== square[i].out[l][j][idx];
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                nocarry[i].b[l][j][idx] <== line[i].out[l][j][idx];
                nocarry[i].b[l][j+2][idx] <== 0;
            }
            compress[i] = Fp12Compress(n, k, 2*k-2, q); // 6 x 4 x k registers in [0, 144 * k^2 * (2k-1) * 2^{4n} )
            for(var l=0; l<6; l++)for(var j=0; j<4; j++)for(var idx=0; idx<3*k-2; idx++)
                compress[i].in[l][j][idx] <== nocarry[i].out[l][j][idx];

            fdouble[i] = Fp12CarryModP(n, k, logc + 1 + 4*n, q);
            for(var l=0; l<6; l++)for(var j=0; j<4; j++)for(var idx=0; idx<k; idx++)
                fdouble[i].in[l][j][idx] <== compress[i].out[l][j][idx]; 

            if( i != 0 || (i == 0 && rBits[i] == 1) ){
                Pdouble[i] = EllipticCurveDouble(n, k, 0, b, q);  
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    Pdouble[i].in[j][idx] <== Pintermed[i+1][j][idx]; 
            }
            
            if(rBits[i] == 0){
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    f[i][l][j][idx] <== fdouble[i].out[l][j][idx]; 
                if( i != 0 ){
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                        Pintermed[i][j][idx] <== Pdouble[i].out[j][idx];
                }
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

                if(i != 0){
                    // Padd[curid] = Pdouble[i] + P 
                    Padd[curid] = EllipticCurveAddUnequal(n, k, q); 
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                        Padd[curid].a[j][idx] <== Pdouble[i].out[j][idx];
                        Padd[curid].b[j][idx] <== P[j][idx];
                    }

                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                        Pintermed[i][j][idx] <== Padd[curid].out[j][idx];
                }
                
                curid++;
            }
        }
    }
    for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[l][j][idx] <== f[0][l][j][idx];
    
}


// Input:
//  g is 6 x 4 x kg array representing element of Fp12, allowing overflow and negative
//  P, Q are as in inputs of LineFunctionEqualNoCarry
// Assume:
//  all registers of g are in [0, 2^{overflowg}) 
//  all registers of P, Q are in [0, 2^n) 
// Output:
//  out = g * l_{P, P}(Q) as element of Fp12 with carry 
//  out is 6 x 2 x k
template Fp12MultiplyWithLineEqual(n, k, kg, overflowg, q){
    signal input g[6][4][kg];
    signal input P[2][k];
    signal input Q[2][6][2][k];
    signal output out[6][2][k];

    component line = LineFunctionEqualNoCarry(n, k); // 6 x 4 x (3k - 2) registers in [0, 2^{3n + 2log(k) + 2})
    for(var l=0; l<2; l++)for(var idx=0; idx<k; idx++)
        line.P[l][idx] <== P[l][idx];
    for(var l=0; l<2; l++)for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        line.Q[l][i][j][idx] <== Q[l][i][j][idx];
    
    component mult = Fp12MultiplyNoCarryUnequal(n, kg, 3*k - 2); // 3k + kg - 3 registers in [0, 12 * min(kg, 3k - 2) * 2^{overflowg + 3n + 2log(k) + 2} )
    
    for(var i=0; i<6; i++)for(var j=0; j<4; j++)for(var idx=0; idx<kg; idx++)
        mult.a[i][j][idx] <== g[i][j][idx];
    for(var i=0; i<6; i++)for(var j=0; j<4; j++)for(var idx=0; idx<3*k-2; idx++)
        mult.b[i][j][idx] <== line.out[i][j][idx];

    component reduce = Fp12Compress(n, k, 2*k + kg - 3, q); // k registers in [0, 12 * (2k + kg - 2) * min(kg, 3k - 2) * 2^{overflowg + 4n + 2log(k) + 2} )
    for(var i=0; i<6; i++)for(var j=0; j<4; j++)for(var idx=0; idx<3*k + kg - 3; idx++)
        reduce.in[i][j][idx] <== mult.out[i][j][idx];

    var mink;
    if(kg < 3*k - 2)
        mink = kg;
    else
        mink = 3*k - 2;
    var logc = log_ceil(12 * (2*k + kg - 2) * mink * k * k);
    
    assert( overflowg + 4*n + logc + 2 < 252 );
    component carry = Fp12CarryModP(n, k, overflowg + 4*n + logc + 2, q);

    for(var i=0; i<6; i++)for(var j=0; j<4; j++)for(var idx=0; idx<k; idx++)
        carry.in[i][j][idx] <== reduce.out[i][j][idx];

    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[i][j][idx] <== carry.out[i][j][idx];
}


// version with one less carry per loop that requires 6n + ... overflow 
// doesn't actually reduce constraints for some reason
template MillerLoop2(n, k, b, r, q){
    signal input P[2][k]; 
    signal input Q[2][6][2][k];
    signal output out[6][2][k];

    var rBits[500]; // length is k * n
    var rBitLength;
    var rSigBits=0;
    for (var i = 0; i < k; i++) {
        for (var j = 0; j < n; j++) {
            rBits[j + n * i] = (r[i] >> j) & 1;
            if(rBits[j + n * i] == 1){
                rSigBits++;
                rBitLength = j + n * i + 1;
            }
        }
    }

    signal Pintermed[rBitLength][2][k]; 
    signal f[rBitLength][6][2][k];

    component Pdouble[rBitLength];
    component Padd[rSigBits];
    component fdouble[rBitLength];
    component square[rBitLength];
    component fadd[rSigBits]; 
    var curid=0;

    var LOGK = log_ceil(k);
    assert( 6*n + LOGK + 6 + log_ceil( 12*(4*k-3) * (2*k-1) * k * k ) < 252 );
    
    for(var i=rBitLength - 1; i>=0; i--){
        if( i == rBitLength - 1 ){
            // f = 1 
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                if(l==0 && j==0 && idx==0)
                    f[i][l][j][idx] <== 1;
                else
                    f[i][l][j][idx] <== 0;
            }
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                Pintermed[i][j][idx] <== P[j][idx];
        }else{
            // compute fdouble[i] = f[i+1]^2 * l_{Pintermed[i+1], Pintermed[i+1]}(Q) 
            square[i] = Fp12MultiplyNoCarry(n, k); // 2k-1 registers in [0, 12 * k * 2^{2n} )
            for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                square[i].a[l][j][idx] <== f[i+1][l][j][idx];
                square[i].a[l][j+2][idx] <== 0;
                square[i].b[l][j][idx] <== f[i+1][l][j][idx];
                square[i].b[l][j+2][idx] <== 0;
            }

            fdouble[i] = Fp12MultiplyWithLineEqual(n, k, 2*k-1, 2*n + LOGK + 4, q);
            // assert ( 6n + log(k) + 6 + log(12 * (4k-3) * (2k-1) * k * k ) ) < 252
            for(var l=0; l<6; l++)for(var j=0; j<4; j++)for(var idx=0; idx<2*k-1; idx++)
                fdouble[i].g[l][j][idx] <== square[i].out[l][j][idx];
            for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                fdouble[i].P[j][idx] <== Pintermed[i+1][j][idx];            
            for(var eps=0; eps<2; eps++)for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                fdouble[i].Q[eps][l][j][idx] <== Q[eps][l][j][idx];

            if( i != 0 || (i == 0 && rBits[i] == 1) ){
                Pdouble[i] = EllipticCurveDouble(n, k, 0, b, q);  
                for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    Pdouble[i].in[j][idx] <== Pintermed[i+1][j][idx]; 
            }
            
            if(rBits[i] == 0){
                for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                    f[i][l][j][idx] <== fdouble[i].out[l][j][idx]; 
                if( i != 0 ){
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                        Pintermed[i][j][idx] <== Pdouble[i].out[j][idx];
                }
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

                if(i != 0){
                    // Padd[curid] = Pdouble[i] + P 
                    Padd[curid] = EllipticCurveAddUnequal(n, k, q); 
                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
                        Padd[curid].a[j][idx] <== Pdouble[i].out[j][idx];
                        Padd[curid].b[j][idx] <== P[j][idx];
                    }

                    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
                        Pintermed[i][j][idx] <== Padd[curid].out[j][idx];
                }
                
                curid++;
            }
        }
    }
    for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[l][j][idx] <== f[0][l][j][idx];
    
}
