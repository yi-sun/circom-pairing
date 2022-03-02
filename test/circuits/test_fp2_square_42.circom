pragma circom 2.0.2;

include "../../circuits/fp2.circom";

component main {public [in, p]} = Fp2square(4, 2);
