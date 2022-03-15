pragma circom 2.0.2;

include "../../circuits/final_exp.circom";

template test(n, k, p){
    signal input in[6][2][k];
    signal output out[6][2][k];
    
    component c = Fp12CyclotomicCompress(n, k);
    for(var i=0; i<6; i++)for(var e=0; e<2; e++)for(var j=0; j<k; j++)
        c.in[i][e][j] <== in[i][e][j];

    component d = Fp12CyclotomicDecompress(n, k, p);
    for(var i=0; i<4; i++)for(var e=0; e<2; e++)for(var j=0; j<k; j++)
        d.in[i][e][j] <== c.out[i][e][j]; 
    
    for(var i=0; i<6; i++)for(var e=0; e<2; e++)for(var j=0; j<k; j++)
        out[i][e][j] <== d.out[i][e][j]; 
    // in === out 
}

component main {public [in]} = test(3, 2, [3, 2]);
