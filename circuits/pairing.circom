pragma circom 2.0.3;

include "curve.circom";
include "final_exp.circom";

template BLSTatePairing(n, k, q){
    signal input P[2][k]; 
    signal input Q[2][6][2][k];
    signal output out[6][2][k];

    component miller = BLSMillerLoop(n, k, q);
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        miller.P[i][idx] <== P[i][idx];
    for(var i=0; i<2; i++)for(var l=0; l<6; l++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        miller.Q[i][l][j][idx] <== Q[i][l][j][idx];
    
    component finalexp = FinalExponentiate(n, k, q);
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        finalexp.in[i][j][idx] <== miller.out[i][j][idx];

    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[i][j][idx] <== finalexp.out[i][j][idx]; 
}
