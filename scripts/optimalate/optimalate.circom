pragma circom 2.0.3;

//include "../../circuits/bls12_381_func.circom";
include "../../circuits/bn254/pairing.circom";

//component main = OptimalAtePairing(55, 7, get_BLS12_381_prime(55, 7) );
component main = OptimalAtePairing(43, 6, get_bn254_prime(43, 6) );
