pragma circom 2.0.2;

include "../../circuits/fp12.circom";

component main {public [a, b]} = Fp12Multiply2(3, 2, [5, 1]);
