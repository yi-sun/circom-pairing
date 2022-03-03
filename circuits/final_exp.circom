pragma circom 2.0.2;

include "bigint.circom";
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
template Fp12cyclotomicCompress(n, k) {
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
template Fp12cyclotomicDecompress(n, k, p) {
    signal input in[4][2][k];
    signal output out[6][2][k]; 

    assert(k<7);
    var LOGK = 3; // LOGK = ceil( log_2( k+1 ) )
    assert(3*n + 5 + LOGK<254);

    var len = 2*k-1; // number of registers in output of Fp2multiplyNoCarry
                     // len = k if using Fp2multiplyNoCarryCompress 

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
    component g5sq = Fp2multiplyNoCarry(n, k); // overflow (k+1) * 2^{2n+1} <= 2^{2n+1+LOGK}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        g5sq.a[2*eps][i] <== in[3][eps][i];
        g5sq.a[2*eps+1][i] <== 0;
        g5sq.b[2*eps][i] <== in[3][eps][i];
        g5sq.b[2*eps+1][i] <== 0;
    }
    // c = 1+u, g5^2 * (1+u)
    signal g5sqc[4][len]; // overflow 2*2^{2n+1+LOGK}
    for(var i=0; i<len; i++){
        g5sqc[0][i] <== g5sq.out[0][i] + g5sq.out[3][i];
        g5sqc[1][i] <== g5sq.out[1][i] + g5sq.out[2][i];
        g5sqc[2][i] <== g5sq.out[0][i] + g5sq.out[2][i];
        g5sqc[3][i] <== g5sq.out[1][i] + g5sq.out[3][i];
    }
    component g4sq3 = Fp2multiplyNoCarry(n, k); // overflow 3*(k+1)* 2^{2n+1} <= 3*2^{2n+1+LOGK}
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
    component g1numRed[4];
    for(var i=0; i<4; i++){
        g1numRed[i] = primeTrickCompression(n, k, k-1, p); 
        for(var j=0; j<2*k-1; j++)
            g1numRed[i].in[j] <== g1num[i][j];
    }
    // compute g1numRed / 4g2
    component g1_1 = Fp2invertCarryModP(n, k, 3*n+4+2*LOGK , p); 
    for(var j=0; j<k; j++){
        for(var i=0; i<4; i++) 
            g1_1.a[i][j] <== g1numRed[i].out[j]; 
        for(var i=0; i<2; i++){
            g1_1.b[2*i][j] <== 4*in[0][i][j]; 
            g1_1.b[2*i+1][j] <== 0;
        }
    }

    // END OF COMPUTATION OF g1 when g2 != 0:



    // COMPUTATION OF g1 when g2 = 0:
    // g1 = 2*g4*g5 / g3
    component twog4g5 = Fp2multiplyNoCarryCompress(n, k, p); // overflow 2*(k+1)*k * 2^{3n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        twog4g5.a[2*eps][i] <== 2*in[2][eps][i];
        twog4g5.a[2*eps+1][i] <== 0;
        twog4g5.b[2*eps][i] <== in[3][eps][i];
        twog4g5.b[2*eps+1][i] <== 0;
    }
    component g1_0 = Fp2invertCarryModP(n, k, 3*n+2+2*LOGK, p);
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
    component twog1sq= Fp2multiplyNoCarry(n, k); // overflow 2*(k+1) * 2^{2n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        twog1sq.a[2*eps][i] <== 2*g1_1.out[eps][i];
        twog1sq.a[2*eps+1][i] <== 0;
        twog1sq.b[2*eps][i] <== g1_1.out[eps][i];
        twog1sq.b[2*eps+1][i] <== 0;
    }
    
    component g2g5 = Fp2multiplyNoCarry(n, k); // overflow (k+1) * 2^{2n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        g2g5.a[2*eps][i] <== in[0][eps][i];
        g2g5.a[2*eps+1][i] <== 0;
        g2g5.b[2*eps][i] <== in[3][eps][i];
        g2g5.b[2*eps+1][i] <== 0;
    }
    
    component threeg3g4 = Fp2multiplyNoCarry(n, k); // overflow 3*(k+1) * 2^{2n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        threeg3g4.a[2*eps][i] <== 3*in[1][eps][i];
        threeg3g4.a[2*eps+1][i] <== 0;
        threeg3g4.b[2*eps][i] <== in[2][eps][i];
        threeg3g4.b[2*eps+1][i] <== 0;
    }
    // 2 g1^2 + g2 g5 - 3 g3 g4 
    var temp[4][100]; 
    for(var i=0; i<len; i++)for(var eps=0; eps<2; eps++){
        temp[2*eps][i] = twog1sq.out[2*eps][i] + g2g5.out[2*eps][i] + threeg3g4.out[2*eps+1][i];
        temp[2*eps+1][i] = twog1sq.out[2*eps+1][i] + g2g5.out[2*eps+1][i] + threeg3g4.out[2*eps][i];
    }
    // (2 g1^2 + g2 g5 - 3 g3 g4)(1+u)
    var tempc[4][100]; // overflow 2*6*(k+1) * 2^{2n+1}
    for(var i=0; i<len; i++){
        tempc[0][i] = temp[0][i] + temp[3][i];
        tempc[1][i] = temp[1][i] + temp[2][i];
        tempc[2][i] = temp[0][i] + temp[2][i];
        tempc[3][i] = temp[1][i] + temp[3][i];
    }
    // (2 g1^2 + g2 g5 - 3 g3 g4)(1+u) + 1 < 2^{2n+LOGK+5}  
    tempc[0][0]++;
    component compress01[4]; // overflow < 2^{3n+2*LOGK+5} 
    for(var i=0; i<4; i++){
        compress01[i] = primeTrickCompression(n, k, k-1, p); 
        for(var j=0; j<2*k-1; j++)
            compress01[i].in[j] <== tempc[i][j];
    }
    // get tempc = p*X + Y 
    component carry_mod01 = Fp2CarryModP(n, k, 3*n+5+2*LOGK, p);
    for(var i=0; i<4; i++)for(var j=0; j<k; j++)
        carry_mod01.in[i][j] <== compress01[i].out[j];
    // g0_1 = Y 
    signal g0_1[2][k]; 
    for(var i=0; i<2; i++)for(var j=0; j<k; j++)
        g0_1[i][j] <== carry_mod01.out[i][j]; 
    // END OF COMPUTATION OF g0 when g2 != 0. 
    
    // COMPUTATION OF g0 when g2 = 0:
    // g0 = (2g1^2 - 3g3g4)(1+u) + 1
    component twog1_0sq= Fp2multiplyNoCarry(n, k); // overflow 2*(k+1) * 2^{2n+1}
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        twog1_0sq.a[2*eps][i] <== 2*g1_0.out[eps][i];
        twog1_0sq.a[2*eps+1][i] <== 0;
        twog1_0sq.b[2*eps][i] <== g1_0.out[eps][i];
        twog1_0sq.b[2*eps+1][i] <== 0;
    }
    // can reuse threeg3g4 ! 
    for(var i=0; i<len; i++)for(var eps=0; eps<2; eps++){
        temp[2*eps][i] = twog1_0sq.out[2*eps][i]  + threeg3g4.out[2*eps+1][i];
        temp[2*eps+1][i] = twog1_0sq.out[2*eps+1][i] + threeg3g4.out[2*eps][i];
    } // overflow 2*5*(k+1) * 2^{2n+1}
    for(var i=0; i<len; i++){
        tempc[0][i] = temp[0][i] + temp[3][i];
        tempc[1][i] = temp[1][i] + temp[2][i];
        tempc[2][i] = temp[0][i] + temp[2][i];
        tempc[3][i] = temp[1][i] + temp[3][i];
    }
    tempc[0][0]++;
    component compress00[4]; // overflow < 2^{3n+2*LOGK+5} 
    for(var i=0; i<4; i++){
        compress00[i] = primeTrickCompression(n, k, k-1, p); 
        for(var j=0; j<2*k-1; j++)
            compress00[i].in[j] <== tempc[i][j];
    }
    component carry_mod00 = Fp2CarryModP(n, k, 3*n+5+2*LOGK, p);
    for(var i=0; i<4; i++)for(var j=0; j<k; j++)
        carry_mod00.in[i][j] <== compress00[i].out[j];
    signal g0_0[2][k]; 
    for(var i=0; i<2; i++)for(var j=0; j<k; j++)
        g0_0[i][j] <== carry_mod00.out[i][j]; 
    // END OF COMPUTATION OF g0 when g2 = 0.    
        
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++)
        out[0][eps][i] <== g0_1[eps][i] + g2isZero.out * (g0_0[eps][i] - g0_1[eps][i]);

}

// output is square of input 
template Fp12cyclotomicSquare(n, k) {
    signal input in[6][2][k];
    signal input p[k];
    signal output out[6][2][k];

    // for now just use plain multiplication, this can be optimized later
    component square = Fp12Multiply(n, k);
    for(var i=0; i<6; i++)for(var j=0; j<k; j++){
        square.a[i][0][j] <== in[i][0][j];
        square.a[i][1][j] <== in[i][1][j];
    
        square.b[i][0][j] <== in[i][0][j];
        square.b[i][1][j] <== in[i][1][j];
    }
    for(var i=0; i<k; i++) square.p[i] <== p[i];

    for(var i=0; i<6; i++)for(var j=0; j<k; j++){
        out[i][0][j] <== square.out[i][0][j];
        out[i][1][j] <== square.out[i][1][j];
    }
}

// assume input is an element of Fp12 in the cyclotomic subgroup GΦ12
// output is input raised to the e-th power
// use the square and multiply method
// assume 0 < e < 2^254
template Fp12cyclotomicExp(n, k, e) {
    assert( e > 0 );

    signal input in[6][2][k];
    signal input p[k];
    signal output out[6][2][k];

    var temp = e;
    var BITLENGTH;
    for(var i=0; i<254; i++){
        if( temp != 0 )
            BITLENGTH = i; 
        temp = temp>>1;
    }
    BITLENGTH++;
    component pow2[BITLENGTH]; // pow2[i] = in^{2^i} 
    component mult[BITLENGTH];

    signal first[6][2][k];
    var curid = 0;

    for(var i=0; i<BITLENGTH; i++){
        // compute pow2[i] = pow2[i-1]**2
        if( i > 0 ){ // pow2[0] is never defined since there is no squaring involved
            pow2[i] = Fp12cyclotomicSquare(n, k);
            for(var j=0; j<k; j++) pow2[i].p[j] <== p[j];
            if( i == 1 ){
                for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                    pow2[i].in[id][eps][j] <== in[id][eps][j];
            }else{
                for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                    pow2[i].in[id][eps][j] <== pow2[i-1].out[id][eps][j];
            }
        }
        if( ((e >> i) & 1) == 1 ){
            if(curid == 0){ // this is the least significant bit
                if( i == 0 ){
                    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                        first[id][eps][j] <== in[id][eps][j];
                }else{
                    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                        first[id][eps][j] <== pow2[i].out[id][eps][j];
                }
            }else{
                // multiply what we already have with pow2[i]
                mult[curid] = Fp12Multiply(n, k); 
                for(var j=0; j<k; j++) mult[curid].p[j] <== p[j];
                for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
                    mult[curid].a[id][eps][j] <== pow2[i].out[id][eps][j];
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
template hard_part(n, k){
    signal input in[6][2][k]; 
    signal input p[k];
    signal output out[6][2][k];

    var x = get_BLS12_381_parameter();  // absolute value of parameter for BLS12-381
    
    // in^{(x+1)/3} 
    component pow1 = Fp12cyclotomicExp(n, k, (x+1)\3 ); 
    for(var i=0; i<k; i++) pow1.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow1.in[id][eps][j] <== in[id][eps][j];
    
    // in^{(x+1)/3 * (x+1)}
    component pow2 = Fp12cyclotomicExp(n, k, x+1); 
    for(var i=0; i<k; i++) pow2.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow2.in[id][eps][j] <== pow1.out[id][eps][j];

    // in^{(x+1)^2/3 * -1} = pow2^-1  inverse = frob(6) in cyclotomic subgroup
    component pow3 = Fp12frobeniusMap(n, k, 6);
    for(var i=0; i<k; i++) pow3.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow3.in[id][eps][j] <== pow2.out[id][eps][j];

    // in^{(x+1)^2/3 * -x} = pow3^x 
    component pow4 = Fp12cyclotomicExp(n, k, x); 
    for(var i=0; i<k; i++) pow4.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow4.in[id][eps][j] <== pow3.out[id][eps][j];

    // in^{(x+1)^2/3 * p} = pow2^p 
    component pow5 = Fp12frobeniusMap(n, k, 1);
    for(var i=0; i<k; i++) pow5.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow5.in[id][eps][j] <== pow2.out[id][eps][j];

    // in^{(x+1)^2/3 * (-x+p)} = pow4 * pow5
    component pow6 = Fp12Multiply(n, k);
    for(var i=0; i<k; i++) pow6.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        pow6.a[id][eps][j] <== pow4.out[id][eps][j];
        pow6.b[id][eps][j] <== pow5.out[id][eps][j];
    }

    // in^{(x+1)^2/3 * (-x+p) * x}  = pow6^x
    component pow7 = Fp12cyclotomicExp(n, k, x);
    for(var i=0; i<k; i++) pow7.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow7.in[id][eps][j] <== pow6.out[id][eps][j];

    // in^{(x+1)^2/3 * (-x+p) * x^2}  = pow7^x
    component pow8 = Fp12cyclotomicExp(n, k, x);
    for(var i=0; i<k; i++) pow8.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow8.in[id][eps][j] <== pow7.out[id][eps][j];

    // in^{(x+1)^2/3 * (-x+p) * q^2} = pow6^{q^2}
    component pow9 = Fp12frobeniusMap(n, k, 2);
    for(var i=0; i<k; i++) pow9.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow9.in[id][eps][j] <== pow6.out[id][eps][j];
    
    // in^{(x+1)^2/3 * (-x+p) * -1} = pow6^{-1} = pow6^{q^6}
    component pow10 = Fp12frobeniusMap(n, k, 6);
    for(var i=0; i<k; i++) pow10.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        pow10.in[id][eps][j] <== pow6.out[id][eps][j];
    
    // in^{(x+1)^2/3 * (-x+p) * (x^2 + q^2)} = pow8 * pow9
    component pow11 = Fp12Multiply(n, k);
    for(var i=0; i<k; i++) pow11.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        pow11.a[id][eps][j] <== pow8.out[id][eps][j];
        pow11.b[id][eps][j] <== pow9.out[id][eps][j];
    }
    
    // in^{(x+1)^2/3 * (-x+p) * (x^2 + q^2 - 1)} = pow10 * pow11
    component pow12 = Fp12Multiply(n, k);
    for(var i=0; i<k; i++) pow12.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        pow12.a[id][eps][j] <== pow10.out[id][eps][j];
        pow12.b[id][eps][j] <== pow11.out[id][eps][j];
    }
    
    // final answer
    // in^{(x+1)^2/3 * (-x+p) * (x^2 + q^2 - 1) + 1} = pow12 * in 
    component pow13 = Fp12Multiply(n, k); 
    for(var i=0; i<k; i++) pow13.p[i] <== p[i];
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++){
        pow13.a[id][eps][j] <== pow12.out[id][eps][j];
        pow13.b[id][eps][j] <== in[id][eps][j];
    }
    
    for(var id=0; id<6; id++)for(var eps=0; eps<2; eps++)for(var j=0; j<k; j++)
        out[id][eps][j] <== pow13.out[id][eps][j];
}
