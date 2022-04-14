pragma circom 2.0.3;

include "../../circuits/curve.circom";
include "../../circuits/bls12_381_func.circom";
//include "../../circuits/extra_curve.circom";

component main = MillerLoop(55, 7, 4,
//BLSMillerLoop(55, 7, 
3,
//15132376222941642752,
get_BLS12_381_prime(55, 7)
);

