pragma circom 2.0.2;

include "../../circuits/fp2.circom";

component main {public [a, b]} = Fp2multiply(4, 2, [1,1]);
