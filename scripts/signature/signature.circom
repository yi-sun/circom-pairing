pragma circom 2.0.3;

include "../../circuits/bls_signature.circom";

component main { public [ pubkey, signature, hash ] } = CoreVerifyPubkeyG1(55, 7);

