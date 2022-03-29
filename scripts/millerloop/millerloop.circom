pragma circom 2.0.3;

include "../../circuits/curve.circom";
include "../../circuits/bls12_381_func.circom";

component main = BLSMillerLoop(55, 7, //4, 
/*[36028792723996673,
20272795337883135,
9049562129190646,
21651489585483456,
31119275314,
0,
0],*/
//[2, 0, 0, 0, 0, 0, 0], 
get_BLS12_381_prime(55, 7)
);

