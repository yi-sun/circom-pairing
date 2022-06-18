pragma circom 2.0.3;

include "../../circuits/bls12_381_func.circom";
include "../../circuits/pairing.circom";

component main = TatePairing(55, 7, get_BLS12_381_prime(55, 7) );

