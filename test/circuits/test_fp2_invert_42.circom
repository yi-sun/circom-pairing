pragma circom 2.0.2;

include "../../circuits/fp2.circom";

component main {public [in]} = Fp2invert(4, 2, [1,1]);
