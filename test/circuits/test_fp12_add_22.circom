pragma circom 2.0.2;

include "../../circuits/fp12.circom";

component main {public [a, b]} = Fp12Add(2, 2, [3,2]);
