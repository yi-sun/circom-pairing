pragma circom 2.0.2;

include "../../circuits/curve.circom";

// Computes (a + b) + c
template AddThree(n, k, a1, b1, p){
    signal input a[2][k]; 
    signal input b[2][k];
    signal input c[2][k];
    signal output out[2][k];
    signal output isInfinity;

    signal d[2][k];
    component ecadd[2];

    for (var i = 0; i < 2; i++){
        ecadd[i] = EllipticCurveAdd(n, k, a1, b1, p);
    }
    
    for (var i = 0; i < 2; i++){
        for (var j = 0; j < k; j++){
            ecadd[0].a[i][j] <== a[i][j];
            ecadd[0].b[i][j] <== b[i][j];
        }
    }
    ecadd[0].aIsInfinity <== 0;
    ecadd[0].bIsInfinity <== 0;

    for (var i = 0; i < 2; i++){
        for (var j = 0; j < k; j++){
            ecadd[1].a[i][j] <== ecadd[0].out[i][j];
            ecadd[1].b[i][j] <== c[i][j];
        }
    }
    ecadd[1].aIsInfinity <== ecadd[0].isInfinity;
    ecadd[1].bIsInfinity <== 0;

    for (var i = 0; i < 2; i++){
        for (var j = 0; j < k; j++){
            out[i][j] <== ecadd[1].out[i][j];
        }
    }

    isInfinity <== ecadd[1].isInfinity;
}

component main {public [a, b, c]} = AddThree(55, 7, 0, 4, [35747322042231467, 36025922209447795, 1084959616957103, 7925923977987733, 16551456537884751, 23443114579904617, 1829881462546425]);
