pragma circom 2.0.2;

include "../../circuits/field_elements.circom";

component main {public [a, b, p]} = Fp12Add(2, 2);
