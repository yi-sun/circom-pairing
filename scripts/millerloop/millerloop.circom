pragma circom 2.0.3;

include "../../circuits/curve.circom";
include "../../circuits/bls12_381_func.circom";
include "../../circuits/extra_curve.circom";

component main = //MillerLoop1(55, 7, 4,
BLSMillerLoop(55, 7, 
//[3, 0 , 0 , 0 , 0 , 0 , 0],
//15132376222941642752,
get_BLS12_381_prime(55, 7)
//[35747322042231467, 36025922209447795, 1084959616957103, 7925923977987733, 16551456537884751, 23443114579904617, 1829881462546425]
);

