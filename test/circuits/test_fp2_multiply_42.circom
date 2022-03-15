pragma circom 2.0.2;

include "../../circuits/fp2.circom";

component main {public [a, b]} = Fp2Multiply(4, 2, [1,1]);
