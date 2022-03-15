pragma circom 2.0.2;

include "../../circuits/final_exp.circom";

component main {public [in]} = Fp12cyclotomicPow4(3, 2, [3, 2]);
