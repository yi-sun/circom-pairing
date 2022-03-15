pragma circom 2.0.2;

include "../../circuits/extra_field_circuits.circom";

component main {public [a, b]} = Fp12Multiply2(3, 2, [5, 1]);
