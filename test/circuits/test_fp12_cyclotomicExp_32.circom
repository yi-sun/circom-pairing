pragma circom 2.0.2;

include "../../circuits/final_exp.circom";

component main {public [in]} = Fp12CyclotomicExp(3, 2, 5, [3, 2]);
