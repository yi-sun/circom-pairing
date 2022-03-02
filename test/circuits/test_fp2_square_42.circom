pragma circom 2.0.2;

include "../../circuits/extra_field_circuits.circom";

component main {public [in, p]} = Fp2square(4, 2);
