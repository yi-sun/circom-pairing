pragma circom 2.0.3;

include "../../circuits/bls12_381_hash_to_G2.circom";

component main { public [ in ] } = SubgroupCheckG2(55, 7);

