pragma circom 2.0.2;

include "../../circuits/field_elements.circom";

component main {public [a, b]} = Fp12Multiply2(3, 2, [5, 1]);
