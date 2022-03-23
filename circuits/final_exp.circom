pragma circom 2.0.3;

include "bigint.circom";
include "bigint_func.circom";
include "field_elements_func.circom";
include "fp2.circom";
include "fp12.circom";

// assume input is an element of Fp12 in the cyclotomic subgroup GΦ₁₂
// A cyclotomic group is a subgroup of Fp^n defined by
//   GΦₙ(p) = {α ∈ Fpⁿ : α^{Φₙ(p)} = 1}

// below we implement compression and decompression for an element  GΦ₁₂ following Theorem 3.1 of https://eprint.iacr.org/2010/542.pdf
// Fp4 = Fp2(w^3) where (w^3)^2 = 1+u 
// Fp12 = Fp4(w) where w^3 = w^3 

// in = g0 + g2 w + g4 w^2 + g1 w^3 + g3 w^4 + g5 w^5 where g_i are elements of Fp2
// out = Compress(in) = [ g2, g3, g4, g5 ] 
template Fp12CyclotomicCompress(n, k) {
    signal input in[6][2][k];
    signal output out[4][2][k]; 

    for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        out[0][eps][j] <== in[1][eps][j];
        out[1][eps][j] <== in[4][eps][j];
        out[2][eps][j] <== in[2][eps][j];
        out[3][eps][j] <== in[5][eps][j];
    } 
}

// in = [g2, g3, g4, g5] where g_i are elements of Fp2
// out[6][2][k] each register in [0, 2^n)
// out = Decompress(in) = g0 + g2 w + g4 w^2 + g1 w^3 + g3 w^4 + g5 w^5 where
// if g2 != 0:
//      g1 = (g5^2 * (1+u) + 3 g4^2 - 2 g3)/(4g2) 
//      g0 = (2 g1^2 + g2 * g5 - 3 g3*g4) * (1+u) + 1
// if g2 = 0:
//      g1 = (2 g4 * g5)/g3
//      g0 = (2 g1^2 - 3 g3 * g4) * (1+u)  + 1    NOTE g0 IS QUARTIC IN THE INPUTS, so I don't think it's worth trying to compute a NoCarry version of g0
// out0 = g0, out1 = g2, out2 = g4, out3 = g1, out4 = g3, out5 = g5
template Fp12CyclotomicDecompress(n, k, p) {
    signal input in[4][2][k];
    signal output out[6][2][k]; 

    var LOGK = log_ceil(k+1); // LOGK = ceil( log_2( k+1 ) )
    assert(3*n + 4 + LOGK < 251);

    var len = 2*k-1; // number of registers in output of Fp2MultiplyNoCarry
                     // len = k if using Fp2MultiplyNoCarryCompress 

    // g2 = in[0], g3 = in[1], g4 = in[2], g5 = in[3]
    for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        out[1][eps][j] <== in[0][eps][j]; 
        out[2][eps][j] <== in[2][eps][j]; 
        out[4][eps][j] <== in[1][eps][j]; 
        out[5][eps][j] <== in[3][eps][j]; 
    }

    // detect if g2 is 0
    component g2Zero[2];
    for(var eps=0; eps<2; eps++){
        g2Zero[eps] = BigIsZero(k);
        for(var i=0; i<k; i++)
            g2Zero[eps].in[i] <== in[0][eps][i];
    }
    var total = 2 - g2Zero[0].out - g2Zero[1].out; 
    component g2isZero = IsZero();
    g2isZero.in <== total;



    // COMPUTATION OF g1 when g2 != 0:
    component g5sq = Fp2MultiplyNoCarry(n, k, 2*n + 1 + LOGK); // overflow (k+1) * 2^{2n+1} <= 2^{2n+1+LOGK}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        g5sq.a[2*eps][i] <== in[3][eps][i];
        g5sq.a[2*eps+1][i] <== 0;
        g5sq.b[2*eps][i] <== in[3][eps][i];
        g5sq.b[2*eps+1][i] <== 0;
    }
    // c = 1+u, g5^2 * (1+u)
    var g5sqc[4][50] = Fp2multc(len, g5sq.out); // overflow 2*2^{2n+1+LOGK}
    component g4sq3 = Fp2MultiplyNoCarry(n, k, 2*n + 3 + LOGK); // overflow 3*(k+1)* 2^{2n+1} <= 3*2^{2n+1+LOGK}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        g4sq3.a[2*eps][i] <== 3*in[2][eps][i];
        g4sq3.a[2*eps+1][i] <== 0;
        g4sq3.b[2*eps][i] <== in[2][eps][i];
        g4sq3.b[2*eps+1][i] <== 0;
    }
    signal g1num[4][len];   // g5^2 * (1+u) + 3 g4^2 - 2 g3
                            // overflow 5*2^{2n+1+LOGK} + 2*2^n < 2^{2n+4+LOGK} 
    for(var i=0; i<len; i++)for(var eps=0; eps<2; eps++){
        g1num[2*eps][i] <== g5sqc[2*eps][i] + g4sq3.out[2*eps][i];
        if(i < k)
            g1num[2*eps+1][i] <== g5sqc[2*eps+1][i] + g4sq3.out[2*eps+1][i] + 2*in[1][eps][i];
        else
            g1num[2*eps+1][i] <== g5sqc[2*eps+1][i] + g4sq3.out[2*eps+1][i];
    }
    // compress g1num using prime reduction trick 
    // overflow k*2^{3n+4+LOGK} < 2^{3n+4+2*LOGK} 
    component g1numRed = Fp2Compress(n, k, k-1, p, 3*n+4+2*LOGK);
    for(var i=0; i<4; i++)for(var j=0; j<2*k-1; j++)
        g1numRed.in[i][j] <== g1num[i][j];
    // compute g1numRed / 4g2
    component g1_1 = Fp2Divide(n, k, 3*n+4+2*LOGK , p); 
    for(var j=0; j<k; j++){
        for(var i=0; i<4; i++) 
            g1_1.a[i][j] <== g1numRed.out[i][j]; 
        for(var i=0; i<2; i++){
            g1_1.b[2*i][j] <== 4*in[0][i][j]; 
            g1_1.b[2*i+1][j] <== 0;
        }
    }

    // END OF COMPUTATION OF g1 when g2 != 0:



    // COMPUTATION OF g1 when g2 = 0:
    // g1 = 2*g4*g5 / g3
    component twog4g5 = Fp2MultiplyNoCarryCompress(n, k, p, n+1,3*n+2+2*LOGK); // overflow 2*(k+1)*k * 2^{3n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        twog4g5.a[2*eps][i] <== 2*in[2][eps][i];
        twog4g5.a[2*eps+1][i] <== 0;
        twog4g5.b[2*eps][i] <== in[3][eps][i];
        twog4g5.b[2*eps+1][i] <== 0;
    }
    component g1_0 = Fp2Divide(n, k, 3*n+2+2*LOGK, p);
    for(var j=0; j<k; j++){
        for(var i=0; i<4; i++)
            g1_0.a[i][j] <== twog4g5.out[i][j]; 
        for(var i=0; i<2; i++){
            g1_0.b[2*i][j] <== in[1][i][j]; 
            g1_0.b[2*i+1][j] <== 0;
        }
    }
    // END OF COMPUTATION OF g1 when g2 = 0.

    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++)
        out[3][eps][i] <== g1_1.out[eps][i] + g2isZero.out * (g1_0.out[eps][i] - g1_1.out[eps][i]);
    


    // COMPUTATION OF g0 when g2 != 0:
    // g0 = (2 g1^2 + g2 g5 - 3 g3 g4 )(1+u) + 1
    component twog1sq= Fp2MultiplyNoCarry(n, k, 2*n + 2 + LOGK); // overflow 2*(k+1) * 2^{2n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        twog1sq.a[2*eps][i] <== 2*g1_1.out[eps][i];
        twog1sq.a[2*eps+1][i] <== 0;
        twog1sq.b[2*eps][i] <== g1_1.out[eps][i];
        twog1sq.b[2*eps+1][i] <== 0;
    }
    
    component g2g5 = Fp2MultiplyNoCarry(n, k, 2*n + 1 + LOGK); // overflow (k+1) * 2^{2n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        g2g5.a[2*eps][i] <== in[0][eps][i];
        g2g5.a[2*eps+1][i] <== 0;
        g2g5.b[2*eps][i] <== in[3][eps][i];
        g2g5.b[2*eps+1][i] <== 0;
    }
    
    component threeg3g4 = Fp2MultiplyNoCarry(n, k, 2*n + 3 + LOGK); // overflow 3*(k+1) * 2^{2n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        threeg3g4.a[2*eps][i] <== 3*in[1][eps][i];
        threeg3g4.a[2*eps+1][i] <== 0;
        threeg3g4.b[2*eps][i] <== in[2][eps][i];
        threeg3g4.b[2*eps+1][i] <== 0;
    }
    // 2 g1^2 + g2 g5 - 3 g3 g4 
    var temp[4][len]; 
    for(var i=0; i<4; i++)for(var j=0; j<len; j++)
        temp[i][j] = twog1sq.out[i][j] + g2g5.out[i][j] + threeg3g4.out[i ^ 1][j];
    
    // (2 g1^2 + g2 g5 - 3 g3 g4)(1+u)
    var tempc[4][50] = Fp2multc(len, temp); // overflow 2*6*(k+1) * 2^{2n+1}
    // (2 g1^2 + g2 g5 - 3 g3 g4)(1+u) + 1 < 2^{2n+LOGK+5}  
    tempc[0][0]++;
    component compress01 = Fp2Compress(n, k, k-1, p, 3*n+5+2*LOGK); // overflow < 2^{3n+2*LOGK+5} 
    for(var i=0; i<4; i++)for(var j=0; j<2*k-1; j++)
        compress01.in[i][j] <== tempc[i][j];
    // get tempc = p*X + Y 
    component carry_mod01 = Fp2CarryModP(n, k, 3*n+5+2*LOGK, p);
    for(var i=0; i<4; i++)for(var j=0; j<k; j++)
        carry_mod01.in[i][j] <== compress01.out[i][j];
    // g0_1 = Y 
    signal g0_1[2][k]; 
    for(var i=0; i<2; i++)for(var j=0; j<k; j++)
        g0_1[i][j] <== carry_mod01.out[i][j]; 
    // END OF COMPUTATION OF g0 when g2 != 0. 
    
    // COMPUTATION OF g0 when g2 = 0:
    // g0 = (2g1^2 - 3g3g4)(1+u) + 1
    component twog1_0sq= Fp2MultiplyNoCarry(n, k, 2*n + 2 + LOGK); // overflow 2*(k+1) * 2^{2n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        twog1_0sq.a[2*eps][i] <== 2*g1_0.out[eps][i];
        twog1_0sq.a[2*eps+1][i] <== 0;
        twog1_0sq.b[2*eps][i] <== g1_0.out[eps][i];
        twog1_0sq.b[2*eps+1][i] <== 0;
    }
    // can reuse threeg3g4 ! 
    for(var i=0; i<4; i++)for(var j=0; j<len; j++)
        temp[i][j] = twog1_0sq.out[i][j]  + threeg3g4.out[i ^ 1][j];
    // overflow 2*5*(k+1) * 2^{2n+1}
    tempc = Fp2multc(len, temp);
    tempc[0][0]++;
    component compress00 = Fp2Compress(n, k, k-1, p, 3*n+5+2*LOGK); // overflow < 2^{3n+2*LOGK+5} 
    for(var i=0; i<4; i++)for(var j=0; j<2*k-1; j++)
        compress00.in[i][j] <== tempc[i][j];
    component carry_mod00 = Fp2CarryModP(n, k, 3*n+5+2*LOGK, p);
    for(var i=0; i<4; i++)for(var j=0; j<k; j++)
        carry_mod00.in[i][j] <== compress00.out[i][j];
    signal g0_0[2][k]; 
    for(var i=0; i<2; i++)for(var j=0; j<k; j++)
        g0_0[i][j] <== carry_mod00.out[i][j]; 
    // END OF COMPUTATION OF g0 when g2 = 0.    
        
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++)
        out[0][eps][i] <== g0_1[eps][i] + g2isZero.out * (g0_0[eps][i] - g0_1[eps][i]);

}

// input is [g2, g3, g4, g5] = C(g) in compressed format of Fp12cyclotomicCompress 
//  input[i][4] keeps track of positive and negatives 
//  don't assume anything about overflow in registers in case we want to square twice before carry
// output is C(g^2) = [h2, h3, h4, h5] computed using Theorem 3.2 of https://eprint.iacr.org/2010/542.pdf 
//  c = 1 + u
//  h2 = 2(g2 + 3*c*B_45) 
//  h3 = 3(A_45 - (c+1)B_45) - 2g3
//  h4 = 3(A_23 - (c+1)B_23) - 2g4
//  h5 = 2(g5 + 3B_23) 
//  A_ij = (g_i + g_j)(g_i + c g_j)
//  B_ij = g_i g_j 

// everything computed with no carries 
// out[4][4][2*k-1] has registers with overflow in [0, ...) where we keep track of positives and negatives 
//  If registers of in[] are in [0, 2^N), then registers of out[] are in [0, 2^{2N + LOGK + 6}) 
//      If moreover in[1] = 0 and in[3]=0, i.e., in[] has no negatives, then registers of out[] are in [0, 2^{2N + LOGK + 6}) 
template Fp12CyclotomicSquareNoCarry(n, k) {
    signal input in[4][4][k];
    signal output out[4][4][2*k-1];
    var LOGK = log_ceil(k+1);

    component B23 = Fp2MultiplyNoCarry(n, k, 2*n + 2 + LOGK); // overflow in 4*(k+1)*2^{2N}     factor of 4 instead of 2 because we allow in[] to have negatives 
    component B45 = Fp2MultiplyNoCarry(n, k, 2*n + 2 + LOGK); 
    for(var i=0; i<4; i++)for(var j=0; j<k; j++){
        B23.a[i][j] <== in[0][i][j];
        B23.b[i][j] <== in[1][i][j];

        B45.a[i][j] <== in[2][i][j];
        B45.b[i][j] <== in[3][i][j];
    }

    component A23 = Fp2MultiplyNoCarry(n, k, 2*n + 5 + LOGK); // overflow in 4*6*(k+1)*2^{2N}     
    component A45 = Fp2MultiplyNoCarry(n, k, 2*n + 5 + LOGK);
    // c*g3 = (1+u)*g3
    var cg3[4][50] = Fp2multc(k, in[1]); 
    var cg5[4][50] = Fp2multc(k, in[3]);
    for(var i=0; i<4; i++)for(var j=0; j<k; j++){
        A23.a[i][j] <== in[0][i][j] + in[1][i][j];  // 2*2^{2N}
        A23.b[i][j] <== in[0][i][j] + cg3[i][j];    // 3*2^{2N}

        A45.a[i][j] <== in[2][i][j] + in[3][i][j];
        A45.b[i][j] <== in[2][i][j] + cg5[i][j];
    }
    
    var cB45[4][50] = Fp2multc(2*k-1, B45.out);  // 8*(k+1)*2^{2N}
    var cB23[4][50] = Fp2multc(2*k-1, B23.out); 
    for(var i=0; i<4; i++)for(var j=0; j<2*k-1; j++){
        if(j < k){
            out[0][i][j] <== 2*(in[0][i][j] + 3*cB45[i][j]);    // 2*2^N + 6*8*(k+1)*2^{2N} < 2^{2N+LOGK+6}  if in[] has no negatives, < 2^{2N+LOGK+5} 
            out[3][i][j] <== 2*(in[3][i][j] + 3*B23.out[i][j]); // 2*2^N + 6*4*(k+1)*2^{2N} < 2^{2N+LOGK+5}  if in[] has no negatives, < 2^{2N+LOGK+4}
            // i ^ 1 flips the lowest bit, which has the effect in our case of "subtraction"
            out[1][i][j] <== 3*(A45.out[i][j] + cB45[i ^ 1][j] + B45.out[i ^ 1][j]) + 2*in[1][i ^ 1][j]; 
                                                                // 3*4*6*(k+1)*2^{2N} + 3*12*(k+1)*2^{2N} + 2*2^N = 108*(k+1)*2^{2N} + 2*2^N < 2^{2N+LOGK+7}   even if in[] has no negatives
            out[2][i][j] <== 3*(A23.out[i][j] + cB23[i ^ 1][j] + B23.out[i ^ 1][j]) + 2*in[2][i ^ 1][j]; 
        }else{
            out[0][i][j] <== 2*(3*cB45[i][j]);   
            out[3][i][j] <== 2*(3*B23.out[i][j]); 
            out[1][i][j] <== 3*(A45.out[i][j] + cB45[i ^ 1][j] + B45.out[i ^ 1][j]); 
            out[2][i][j] <== 3*(A23.out[i][j] + cB23[i ^ 1][j] + B23.out[i ^ 1][j]); 
        }
    }
}

// input is [g2, g3, g4, g5] = C(g) is output of Fp12cyclotomicCompress 
// assume in[4][2][k] has registers in [0, 2^n) with in[i][j] in [0, p)
// output is C(g^2) 
// out[4][2][k] has registers in [0, 2^n) with out[i][j] in [0, p) 
template Fp12CyclotomicSquare(n, k, p) {
    signal input in[4][2][k]; 
    signal output out[4][2][k];

    var LOGK = log_ceil(k+1);
    assert(3*n + 2*LOGK + 7 < 252);

    component sq = Fp12CyclotomicSquareNoCarry(n, k);
    for(var i=0; i<4; i++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        sq.in[i][2*eps][j] <== in[i][eps][j]; 
        sq.in[i][2*eps+1][j] <== 0; 
    }
    // sq.in has no negatives
    // by comments above, sq.out has registers in [0, 2^{2n + LOGK + 6}) 
    component sqRed[4]; 
    component sqMod[4];
    for(var i=0; i<4; i++){
        sqRed[i] = Fp2Compress(n, k, k-1, p, 3*n+2*LOGK+7); 
        for(var eps=0; eps<4; eps++)for(var j=0; j<2*k-1; j++)
            sqRed[i].in[eps][j] <== sq.out[i][eps][j]; 
        // sqRed[i].out has registers in [0, k*2^{3n+LOGK+7} < 2^{3n+2*LOGK+7} ) 
        sqMod[i] = Fp2CarryModP(n, k, 3*n+2*LOGK+7, p);
        for(var eps=0; eps<4; eps++)for(var j=0; j<k; j++)
            sqMod[i].in[eps][j] <== sqRed[i].out[eps][j]; 
    }
    for(var i=0; i<4; i++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        out[i][eps][j] <== sqMod[i].out[eps][j];
}

// assume input is an element of Fp12 in the cyclotomic subgroup GΦ12
// output is input raised to the e-th power
// use the square and multiply method
// assume 0 < e < 2^254
template Fp12CyclotomicExp(n, k, e, p) {
    assert( e > 0 );

    signal input in[6][2][k];
    signal output out[6][2][k];

    var temp = e;
    // get bitlength of e 
    var BITLENGTH;
    for(var i=0; i<254; i++){
        if( temp != 0 )
            BITLENGTH = i; 
        temp = temp>>1;
    }
    BITLENGTH++;

    // Compress in[]
    component Cin = Fp12CyclotomicCompress(n, k);
    for(var i=0; i<6; i++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        Cin.in[i][eps][j] <== in[i][eps][j];
    
    component pow2[BITLENGTH]; // pow2[i] = C(in^{2^i})
    component Dpow2[BITLENGTH];
    component mult[BITLENGTH];

    signal first[6][2][k];
    var curid = 0; // tracks current index in mult[] 

    for(var i=0; i<BITLENGTH; i++){
        // compute pow2[i] = pow2[i-1]**2
        if( i > 0 ){ // pow2[0] is never defined since there is no squaring involved
            pow2[i] = Fp12CyclotomicSquare(n, k, p);
            if( i == 1 ){
                for(var id=0; id<4; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                    pow2[i].in[id][eps][j] <== Cin.out[id][eps][j];
            }else{
                for(var id=0; id<4; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                    pow2[i].in[id][eps][j] <== pow2[i-1].out[id][eps][j];
            }
        }
        if( ((e >> i) & 1) == 1 ){
            // decompress pow2[i] so we can use it 
            if( i > 0 ){
                Dpow2[curid] = Fp12CyclotomicDecompress(n, k, p);
                for(var id=0; id<4; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                    Dpow2[curid].in[id][eps][j] <== pow2[i].out[id][eps][j];
            }
            if(curid == 0){ // this is the least significant bit
                if( i == 0 ){
                    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                        first[id][eps][j] <== in[id][eps][j];
                }else{
                    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                        first[id][eps][j] <== Dpow2[curid].out[id][eps][j];
                }
            }else{
                // multiply what we already have with pow2[i]
                mult[curid] = Fp12Multiply(n, k, p); 
                for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                    mult[curid].a[id][eps][j] <== Dpow2[curid].out[id][eps][j];
                if(curid == 1){
                    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                        mult[curid].b[id][eps][j] <== first[id][eps][j];
                }else{
                    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                        mult[curid].b[id][eps][j] <== mult[curid-1].out[id][eps][j];
                }
            } 
            curid++; 
        }
    }
    curid--;
    if(curid == 0){
        for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
            out[id][eps][j] <== first[id][eps][j];
    }else{
        for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
            out[id][eps][j] <== mult[curid].out[id][eps][j];
    }
}

// hard part of final exponentiation
// use equation at top of p.14 from https://eprint.iacr.org/2020/875.pdf
template FinalExpHardPart(n, k, p){
    signal input in[6][2][k]; 
    signal output out[6][2][k];

    var x = get_BLS12_381_parameter();  // absolute value of parameter for BLS12-381
    
    // in^{(x+1)/3} 
    component pow1 = Fp12CyclotomicExp(n, k, (x+1)\3, p); 
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow1.in[id][eps][j] <== in[id][eps][j];
    
    // in^{(x+1)/3 * (x+1)}
    component pow2 = Fp12CyclotomicExp(n, k, x+1, p); 
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow2.in[id][eps][j] <== pow1.out[id][eps][j];

    // in^{(x+1)^2/3 * -1} = pow2^-1  inverse = frob(6) in cyclotomic subgroup
    component pow3 = Fp12FrobeniusMap(n, k, 6);
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow3.in[id][eps][j] <== pow2.out[id][eps][j];

    // in^{(x+1)^2/3 * -x} = pow3^x 
    component pow4 = Fp12CyclotomicExp(n, k, x, p); 
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow4.in[id][eps][j] <== pow3.out[id][eps][j];

    // in^{(x+1)^2/3 * p} = pow2^p 
    component pow5 = Fp12FrobeniusMap(n, k, 1);
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow5.in[id][eps][j] <== pow2.out[id][eps][j];

    // in^{(x+1)^2/3 * (-x+p)} = pow4 * pow5
    component pow6 = Fp12Multiply(n, k, p);
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        pow6.a[id][eps][j] <== pow4.out[id][eps][j];
        pow6.b[id][eps][j] <== pow5.out[id][eps][j];
    }

    // in^{(x+1)^2/3 * (-x+p) * x}  = pow6^x
    component pow7 = Fp12CyclotomicExp(n, k, x, p);
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow7.in[id][eps][j] <== pow6.out[id][eps][j];

    // in^{(x+1)^2/3 * (-x+p) * x^2}  = pow7^x
    component pow8 = Fp12CyclotomicExp(n, k, x, p);
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow8.in[id][eps][j] <== pow7.out[id][eps][j];

    // in^{(x+1)^2/3 * (-x+p) * q^2} = pow6^{q^2}
    component pow9 = Fp12FrobeniusMap(n, k, 2);
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow9.in[id][eps][j] <== pow6.out[id][eps][j];
    
    // in^{(x+1)^2/3 * (-x+p) * -1} = pow6^{-1} = pow6^{q^6}
    component pow10 = Fp12FrobeniusMap(n, k, 6);
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow10.in[id][eps][j] <== pow6.out[id][eps][j];
    
    // in^{(x+1)^2/3 * (-x+p) * (x^2 + q^2)} = pow8 * pow9
    component pow11 = Fp12Multiply(n, k, p);
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        pow11.a[id][eps][j] <== pow8.out[id][eps][j];
        pow11.b[id][eps][j] <== pow9.out[id][eps][j];
    }
    
    // in^{(x+1)^2/3 * (-x+p) * (x^2 + q^2 - 1)} = pow10 * pow11
    component pow12 = Fp12Multiply(n, k, p);
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        pow12.a[id][eps][j] <== pow10.out[id][eps][j];
        pow12.b[id][eps][j] <== pow11.out[id][eps][j];
    }
    
    // final answer
    // in^{(x+1)^2/3 * (-x+p) * (x^2 + q^2 - 1) + 1} = pow12 * in 
    component pow13 = Fp12Multiply(n, k, p); 
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        pow13.a[id][eps][j] <== pow12.out[id][eps][j];
        pow13.b[id][eps][j] <== in[id][eps][j];
    }
    
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        out[id][eps][j] <== pow13.out[id][eps][j];
}

// easy part of final exponentiation 
// out = in^{ (q^6 - 1)*(q^2 + 1) }
template FinalExpEasyPart(n, k, p){
    signal input in[6][2][k];
    signal output out[6][2][k];
    
    // in^{q^6} 
    component f1 = Fp12FrobeniusMap(n, k, 6); 
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        f1.in[id][eps][j] <== in[id][eps][j];

    // in^{-1}
    component f2 = Fp12Invert(n, k, p);
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        f2.in[id][eps][j] <== in[id][eps][j];

    // in^{q^6 - 1}
    component f3 = Fp12Multiply(n, k, p);
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        f3.a[id][eps][j] <== f1.out[id][eps][j];
        f3.b[id][eps][j] <== f2.out[id][eps][j];
    }
    
    // in^{(q^6-1)*q^2} = f3^{q^2}
    component f4 = Fp12FrobeniusMap(n, k, 2); 
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        f4.in[id][eps][j] <== f3.out[id][eps][j];
    
    // in^{(q^6-1)(q^2+1)} = f4 * f3
    component f5 = Fp12Multiply(n, k, p);
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        f5.a[id][eps][j] <== f3.out[id][eps][j];
        f5.b[id][eps][j] <== f4.out[id][eps][j];
    }
    
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        out[id][eps][j] <== f5.out[id][eps][j];
}

// out = in^{(q^12-1)/r} = FinalExpHardPart( FinalExpEasyPart(in) )
template FinalExponentiate(n, k, p){
    signal input in[6][2][k];
    signal output out[6][2][k];

    component f1 = FinalExpEasyPart(n, k, p);
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        f1.in[id][eps][j] <== in[id][eps][j];
    
    component f = FinalExpHardPart(n, k, p);
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        f.in[id][eps][j] <== f1.out[id][eps][j];
    
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        out[id][eps][j] <== f.out[id][eps][j];
}
