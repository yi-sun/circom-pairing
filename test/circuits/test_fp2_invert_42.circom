pragma circom 2.0.2;

include "../../circuits/field_elements.circom";

component main {public [in, p]} = Fp2invert(4, 2);
