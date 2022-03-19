pragma circom 2.0.2;

include "../../circuits/curve.circom";

component main {public [P]} = Fp12MultiplyWithLineEqual(3, 2, 2, 3, [3, 2]);
