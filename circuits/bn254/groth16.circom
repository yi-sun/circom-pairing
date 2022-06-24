pragma circom 2.0.3;

include "bn254_func.circom";
include "subgroup_check.circom";
include "curve.circom";
include "pairing.circom";

template verifyProof(publicInputCount) {
    // BN254 facts
    var n = 43;
    var k = 6;
    var p[50] = get_bn254_prime(n, k);
    var b = 3;
    var b2[2][50] = get_bn254_b(n, k);
    var loopBitLength = 65;
    var pseudoBinaryEncoding[loopBitLength] = get_bn254_pseudo_binary_encoding();

    // verification key
    signal input negalfa1xbeta2[6][2][k]; // e(-alfa1, beta2)
    signal input gamma2[2][2][k];
    signal input delta2[2][2][k];
    signal input IC[publicInputCount+1][2][k];

    // proof
    signal input negpa[2][k];
    signal input pb[2][2][k];
    signal input pc[2][k];
    signal input pubInput[publicInputCount];

    signal output out;

    // check proof consists of valid group elements
    component negpaInG1 = SubgroupCheckG1(n, k);
    component pbInG2 = SubgroupCheckG2(n, k);
    component pcInG1 = SubgroupCheckG1(n, k);

    for (var i = 0;i < 2;i++) {
        for (var j = 0;j < 2;j++) {
            for (var idx = 0;idx < k;idx++) {
                pbInG2.in[i][j][idx] <== pb[i][j][idx];
            }
        }

        for (var idx = 0;idx < k;idx++) {
            negpaInG1.in[i][idx] <== negpa[i][idx];
            pcInG1.in[i][idx] <== pc[i][idx];
        }
    }

    // Compute VK = sum pubInput[i]*IC[i]
    component ICmultInp[publicInputCount];
    component ICPrefAddInp[publicInputCount];

    for (var i = 0;i < publicInputCount;i++) {
        ICmultInp[i] = EllipticCurveScalarMultiplySignalX(n, k, b, p);
        for (var j = 0;j < 2;j++) {
            for (var idx = 0;idx < k;idx++) ICmultInp[i].in[j][idx] <== IC[i + 1][j][idx];
        }
        ICmultInp[i].inIsInfinity <== 0;
        ICmultInp[i].x <== pubInput[i];

        ICPrefAddInp[i] = EllipticCurveAdd(n, k, 0, b, p);
        for (var j = 0;j < 2;j++) {
            for (var idx = 0;idx < k;idx++) {
                ICPrefAddInp[i].a[j][idx] <== ICmultInp[i].out[j][idx];
                ICPrefAddInp[i].b[j][idx] <== i == 0 ? IC[0][j][idx] : ICPrefAddInp[i - 1].out[j][idx];
            }
        }
        ICPrefAddInp[i].aIsInfinity <== ICmultInp[i].isInfinity;
        ICPrefAddInp[i].bIsInfinity <== i == 0 ? 0 : ICPrefAddInp[i - 1].isInfinity;
    }

    signal VK[2][k];
    for (var i = 0;i < 2;i++) {
        for (var idx = 0;idx < k;idx++) VK[i][idx] <== ICPrefAddInp[publicInputCount - 1].out[i][idx];
    }


    // compute e'(-A, B)*e'(VK, gamma2)*e'(C, delta2)
    // e'(-A, B) - optmult[0]
    // e'(VK, gamma2) - optmult[1]
    // e'(C, delta2) - optmult[2]
    component optmult = OptimizedMillerLoopProductFp2(n, k, b2[0], b2[1], 3, pseudoBinaryEncoding, loopBitLength, p);

    for (var i = 0;i < 2;i++) {
        for (var j = 0;j < 2;j++) {
            for (var idx = 0;idx < k;idx++) {
                optmult.P[0][i][j][idx] <== pb[i][j][idx];
                optmult.P[1][i][j][idx] <== gamma2[i][j][idx];
                optmult.P[2][i][j][idx] <== delta2[i][j][idx];
            }
        }
    }

    for (var i = 0;i < 2;i++) {
        for (var idx = 0;idx < k;idx++) {
            optmult.Q[0][i][idx] <== negpa[i][idx];
            optmult.Q[1][i][idx] <== VK[i][idx];
            optmult.Q[2][i][idx] <== pc[i][idx];            
        }
    }

    // exponentiate to get e(-A, B)*e(VK, gamma2)*e(C, delta2)
    component finalexp = FinalExponentiate(n, k, p);
    for(var i=0; i<6; i++)for(var j=0; j<2; j++)for(var idx=0; idx<k; idx++)
        finalexp.in[i][j][idx] <== optmult.out[i][j][idx];

    // check e(-alpha1, beta1) === e(-A, B)*e(VK, gamma2)*e(C, delta2)
    component areBigEqual[6][2], areFP12PrefixEqual[6];
    for (var i = 0;i < 6;i++) {
        for (var j = 0;j < 2;j++) {
            areBigEqual[i][j] = BigIsEqual(k);
            for (var idx = 0;idx < k;idx++) {
                areBigEqual[i][j].a[idx] <== finalexp.out[i][j][idx];
                areBigEqual[i][j].b[idx] <== negalfa1xbeta2[i][j][idx];
            }
        }
    }

    for (var i = 0;i < 6;i++) {
        areFP12PrefixEqual[i] = AND();
        areFP12PrefixEqual[i].a <== i == 0 ? 1 : areFP12PrefixEqual[i-1].out;
        areFP12PrefixEqual[i].b <== areBigEqual[i][0].out*areBigEqual[i][1].out;
    }
    out <== areFP12PrefixEqual[5].out;
}
