pragma circom 2.0.2;

include "../../circuits/field_elements.circom";

component main {public [a, b]} = Fp2multiply(4, 2, [1,1]);
