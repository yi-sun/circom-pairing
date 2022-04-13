pragma circom 2.0.3;

include "curve.circom";
include "final_exp.circom";
include "curve_fp2.circom";
include "bls12_381_func.circom";

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

    /*
    // check out[i][j] < p since we never did that previously
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

template BLSAtePairing(n, k, q){
    signal input P[2][2][k]; 
    signal input Q[2][k];
    signal output out[6][2][k];

    var x = get_BLS12_381_parameter();

    signal negP[2][2][k]; // (x, -y) mod q
    // currently EllipticCurveDoubleFp2 relies on 0 <= in < p so we cannot pass a negative number for negP
    component neg[2];
    for(var j=0; j<2; j++){
        neg[j] = FpNegate(n, k, q); 
        for(var idx=0; idx<k; idx++)
            neg[j].in[idx] <== P[1][j][idx];
        for(var idx=0; idx<k; idx++){
            negP[0][j][idx] <== P[0][j][idx];
            negP[1][j][idx] <== neg[j].out[idx];
        }
    }
    component miller = MillerLoopFp2(n, k, [4,4], x, q);
    for(var i=0; i<2; i++)for(var j=0; j <2; j++)for(var idx=0; idx<k; idx++)
        miller.P[i][j][idx] <== negP[i][j][idx];
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++)
        miller.Q[i][idx] <== Q[i][idx];
    /*for(var i = 0; i <6; i ++)for(var j = 0; j < 2; j++)for(var l = 0; l < k; l++) {
        if (i==0 && j == 0 && l == 0) {
            miller.in[i][j][l] <== 1;
        }
        else {
            miller.in[i][j][l] <== 0;
        }
    }*/
    
    component finalexp = FinalExponentiate(n, k, q);
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        finalexp.in[i][j][idx] <== miller.out[i][j][idx];

    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        out[i][j][idx] <== finalexp.out[i][j][idx]; 

    /*
    // check out[i][j] < p since we never did that previously
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

// Input: g1, pubkey in G_1 
//        signature, H(m) in G_2
// Check that e(g1, signature) = e(pubkey, H(m)) by checking e(g1, signature)*e(pubkey, -H(m)) === 1 where e(,) is optimal Ate pairing
template BLSSignatureMinPubkeySize(n, k, q){
    signal input g1[2][k]; 
    signal input pubkey[2][k];
    signal input signature[2][2][k];
    signal input Hm[2][2][k];

    var x = get_BLS12_381_parameter();

    signal neg_s[2][2][k];
    component neg[2];
    for(var j=0; j<2; j++){
        neg[j] = FpNegate(n, k, q); 
        for(var idx=0; idx<k; idx++)
            neg[j].in[idx] <== signature[1][j][idx];
        for(var idx=0; idx<k; idx++){
            neg_s[0][j][idx] <== signature[0][j][idx];
            neg_s[1][j][idx] <== neg[j].out[idx];
        }
    }

    component miller = MillerLoopFp2Two(n, k, [4,4], x, q);
    for(var i=0; i<2; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        miller.P[0][i][j][idx] <== neg_s[i][j][idx];
        miller.P[1][i][j][idx] <== Hm[i][j][idx];
    }
    for(var i=0; i<2; i++)for(var idx=0; idx<k; idx++){
        miller.Q[0][i][idx] <== g1[i][idx];
        miller.Q[1][i][idx] <== pubkey[i][idx];
    }

    component finalexp = FinalExponentiate(n, k, q);
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        finalexp.in[i][j][idx] <== miller.out[i][j][idx];

    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++){
        if(i==0 && j==0 && idx==0)
            finalexp.out[i][j][idx] === 1;
        else
            finalexp.out[i][j][idx] === 0;
    }
}
