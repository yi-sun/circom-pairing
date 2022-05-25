pragma circom 2.0.3;

include "../curve.circom";

// Curve E : y^2 = x^3 + b
// Inputs:
//  in is 2 x k array where P = (x, y) is a point in E(Fp) 
//  inIsInfinity = 1 if P = O, else = 0
// Output:
//  out = [x]P is 2 x k array representing a point in E(Fp)
//  isInfinity = 1 if [x]P = O, else = 0
// Assume:
//  x in bn254 field size (or snark field size)
//  `in` is point in E even if inIsInfinity = 1 just so nothing goes wrong
//  E(Fp) has no points of order 2
template EllipticCurveScalarMultiplySignalX(n, k, b, p){
    signal input in[2][k];
    signal input inIsInfinity;
    signal input x;

    signal output out[2][k];
    signal output isInfinity;
        
    var BitLength = 254; // do you need more for safety or something?
    component Bits = Num2Bits(254);
    Bits.in <== x;

    signal R[BitLength + 1][2][k]; 
    signal R_isO[BitLength + 1];
    signal addendum[BitLength][2][k];
    component Pdouble[BitLength];
    component Padd[BitLength];

    // if in = O then [x]O = O so there's no point to any of this
    signal P[2][k];
    for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        P[j][idx] <== in[j][idx];

    for (var j = 0;j < 2;j++) {
        for (var idx = 0;idx < k;idx++) {
            R[BitLength][j][idx] <== P[j][idx];
        }
    }
    R_isO[BitLength] <== 1;

    for (var i = BitLength - 1;i >= 0;i--) {
        // E(Fp) has no points of order 2, so the only way 2*R[i+1] = O is if R[i+1] = O 
        Pdouble[i] = EllipticCurveDouble(n, k, 0, b, p);  
        for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
            Pdouble[i].in[j][idx] <== R[i+1][j][idx]; 
        
        // Padd[curid] = Pdouble[i] + P if bits.out[i] == 1
        // Padd[curid] = Pdouble[i] + O if bits.out[i] == 0

        for (var j = 0;j < 2;j++) {
            for (var idx = 0;idx < k;idx++) {
                addendum[i][j][idx] <== Bits.out[i]*P[j][idx];
            }
        }

        Padd[i] = EllipticCurveAdd(n, k, 0, b, p);
        for (var j = 0;j < 2;j++) {
            for (var idx = 0;idx < k;idx++) {
                Padd[i].a[j][idx] <== Pdouble[i].out[j][idx]; 
                Padd[i].b[j][idx] <== addendum[i][j][idx];
            }
        }
        Padd[i].aIsInfinity <== R_isO[i+1];
        Padd[i].bIsInfinity <== 1 - Bits.out[i];

        R_isO[i] <== Padd[i].isInfinity; 
        for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
            R[i][j][idx] <== Padd[i].out[j][idx];
    }

    // output = O if input = O or R[0] = O 
    isInfinity <== inIsInfinity + R_isO[0] - inIsInfinity * R_isO[0];
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        out[i][idx] <== R[0][i][idx] + isInfinity * (in[i][idx] - R[0][i][idx]);
}