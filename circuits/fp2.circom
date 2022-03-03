pragma circom 2.0.2;

include "bigint.circom";
include "field_elements_func.circom";

// add two elements in Fp2
template Fp2Add(n, k) {
    signal input a[2][k];
    signal input b[2][k];
    signal input p[k];
    signal output c[2][k];

    component adders[2];
    for (var i = 0; i < 2; i++) {
        adders[i] = BigAddModP(n, k);
        for (var j = 0; j < k; j++) {
            adders[i].a[j] <== a[i][j];
            adders[i].b[j] <== b[i][j];
            adders[i].p[j] <== p[j];
        }   
        for (var j = 0; j < k; j ++) {
            c[i][j] <== adders[i].out[j];
        }
    }
}


// p has k registers 
// inputs: 
//  a[4][k] allow overflow
//  b[4][k] 
// outputs:
//  out[4][2*k-1] such that 
//      ( (a[0] - a[1]) + (a[2] - a[3])*u ) * ( (b[0] - b[1]) + (b[2] - b[3])*u ) 
//      = (out[0] - out[1]) + (out[2] - out[3])*u  
//      we keep track of "positive" and "negatives" since circom isn't able to 
//      if each a[i][j] is in [0, B) then out[i][j] is in [0, 4*(k+1)*B^2 )
//  out[i] has 2*k-1 registers since that's output of BigMultShortLong
template Fp2multiplyNoCarry(n, k){
    signal input a[4][k];
    signal input b[4][k];
    signal output out[4][2*k-1];

    // if a,b registers were in [0,2^n), then
    // var LOGK = 3;
    // assert(2*n + 2 + LOGK <254);

    component ab[4][4];
    for(var i=0; i<4; i++)for(var j=0; j<4; j++){
        ab[i][j] = BigMultShortLong(n, k); // output has 2*k-1 registers
        for(var l=0; l<k; l++){
            ab[i][j].a[l] <== a[i][l];
            ab[i][j].b[l] <== b[j][l];
        }
    }
    
    for(var j=0; j<2*k-1; j++){
        out[0][j] <== ab[0][0].out[j] + ab[1][1].out[j] + ab[2][3].out[j] + ab[3][2].out[j];
        out[1][j] <== ab[0][1].out[j] + ab[1][0].out[j] + ab[2][2].out[j] + ab[3][3].out[j];
        out[2][j] <== ab[0][2].out[j] + ab[1][3].out[j] + ab[2][0].out[j] + ab[3][1].out[j];
        out[3][j] <== ab[0][3].out[j] + ab[1][2].out[j] + ab[2][1].out[j] + ab[3][0].out[j];
    }
}

// input is 4 x (k+m) where registers can overflow [0,B)
//  in[0] - in[1] + (in[2] - in[3])*u
// output is congruent to input (mod p) and represented as 4 x k where registers overflow [0,(m+1)*2^n*B) 
template Fp2Compress(n, k, m, p){
    signal input in[4][k+m]; 
    signal output out[4][k];
    
    component c[4];
    for(var i=0; i<4; i++){
        c[i] = primeTrickCompression(n, k, m, p);
        for(var j=0; j<k+m; j++)
            c[i].in[j] <== in[i][j]; 
    }
    for(var i=0; i<4; i++)for(var j=0; j<k; j++)
        out[i][j] <== c[i].out[j];
}
// same input as above
// outputs:
//  out[4][k] such that 
//      out[i] has k registers because we use the "prime trick" to compress from 2*k-1 to k registers 
//      if each a[i][j] is in [0, B) then out[i][j] is in [0, 4*(k+1)*k*2^n*B^2 )
//          (k+1)B^2 from BigMultShortLong
//          *4 from adding 
//          *k*2^n from prime trick
template Fp2multiplyNoCarryCompress(n, k, p){
    signal input a[4][k];
    signal input b[4][k];
    signal output out[4][k];
    
    component ab = Fp2multiplyNoCarry(n, k);
    for(var i=0; i<4; i++)for(var j=0; j<k; j++){
        ab.a[i][j] <== a[i][j];
        ab.b[i][j] <== b[i][j]; 
    }
    
    component compress = Fp2Compress(n, k, k-1, p);
    for(var i=0; i<4; i++)for(var j=0; j<2*k-1; j++)
        compress.in[i][j] <== ab.out[i][j]; 
 
    for(var i=0; i<4; i++)for(var j=0; j<k; j++)
        out[i][j] <== compress.out[i][j];
}

// check if in[0] + in[0]*u is a valid point of Fp2 with in[0],in[1] both with k registers in [0,2^n) and in[i] in [0,p)
template checkValidFp2(n, k, p){
    signal input in[2][k];
    component range_checks[2][k];
    component lt[2];
    
    for(var eps=0; eps<2; eps++){
        lt[eps] = BigLessThan(n, k);
        for(var i=0; i<k; i++){
            range_checks[eps][i] = Num2Bits(n);
            range_checks[eps][i].in <== in[eps][i];
            
            lt[eps].a[i] <== in[eps][i];
            lt[eps].b[i] <== p[i];
        }
        lt[eps].out === 1;
    }    
}

// copying format of Fp12CarryModP
// solve for in0 - in1 + (in2 - in3) = p * X + out
// X has registers lying in [-2^n, 2^n) 
// X has at most Ceil( overflow / n ) registers 
// assume in has registers in [0, 2^overflow) 
template Fp2CarryModP(n, k, overflow, p){
    signal input in[4][k]; 
    var m = (overflow + n - 1) \ n; 
    signal output X[2][m];
    signal output out[2][k];

    var Xvar[2][2][100] = Fp2_long_div(n, k, m, in, p); 
    component range_check = checkValidFp2(n, k, p);
    component X_range_checks[2][m];
    
    for(var eps=0; eps<2; eps++){
        for(var i=0; i<k; i++){
            out[eps][i] <-- Xvar[eps][1][i];
            range_check.in[eps][i] <== out[eps][i];
        }
        
        for(var i=0; i<m; i++){
            X[eps][i] <-- Xvar[eps][0][i];
            X_range_checks[eps][i] = Num2Bits(n+1);
            X_range_checks[eps][i].in <== X[eps][i] + (1<<n); // X[eps][i] should be between [-2^n, 2^n)
        }
    }

    component mod_check[2];
    for(var eps=0; eps<2; eps++){
        mod_check[eps] = checkBigMod(n, k, m, overflow, p);
        for(var i=0; i<k; i++){
            mod_check[eps].in[i] <== in[2*eps][i] - in[2*eps+1][i];
            mod_check[eps].Y[i] <== out[eps][i];
        }
        for(var i=0; i<m; i++){
            mod_check[eps].X[i] <== X[eps][i];
        }
    }
}
// multiplication specialized to Fp2
// outputs a*b in Fp2 
// (a0 + a1 u)*(b0 + b1 u) = (a0*b0 - a1*b1) + (a0*b1 + a1*b0)u 
// out[i] has k registers each in [0, 2^n)
// out[i] in [0, p)
template Fp2multiply(n, k, p){
    signal input a[2][k];
    signal input b[2][k];
    signal output out[2][k];

    assert(k<=7); // k=7 probably works but best to be safe
    var LOGK = 3; // LOGK = ceil( log_2( (k+1) ) )
    assert(3*n + 2 + 2*LOGK<254);

    component c = Fp2multiplyNoCarryCompress(n, k, p); 
    for(var i=0; i<k; i++){
        c.a[0][i] <== a[0][i];
        c.a[1][i] <== 0;
        c.a[2][i] <== a[1][i];
        c.a[3][i] <== 0;
        c.b[0][i] <== b[0][i];
        c.b[1][i] <== 0;
        c.b[2][i] <== b[1][i];
        c.b[3][i] <== 0;
    }
    // bounds below say X[eps] will only require 4 registers max 
    var m = 4;
    component carry_mod = Fp2CarryModP(n, k, 3*n+2+2*LOGK, p); // 3n+1+2*LOGK probably enough but safety first
    for(var i=0; i<4; i++)for(var j=0; j<k; j++)
        carry_mod.in[i][j] <== c.out[i][j]; 
    
    for(var i=0; i<2; i++)for(var j=0; j<k; j++)
        out[i][j] <== carry_mod.out[i][j]; 

    // out[0] constraint: X = X[0], Y = out[0] 
    // c0 = a0b0, c1 = a1b1
    // constrain by Carry( c0 -' c1 - p *' X - Y ) = 0 
    // where all operations are performed without carry 
    // each register is an overflow representation in the range 
    //      (-(k+1)*k*2^{3n}-2^{2n}-2^n, (k+1)*k*2^{3n}+2^{2n} )
    //      which is contained in (-2^{3n + 2*LOGK + 1}, 2^{3n + 2*LOGK + 1})

    // out[1] constraint: X = X[1], Y = out[1]
    // c2 = a0b1 + a1b0, c3 = 0
    // constrain by Carry( c2 -' c3 -' p *' X - Y) = 0 
    // each register is an overflow representation in the range 
    //      (-2^{2n}-2^n, (k+1)*k*2^{3n+1} + 2^{2n} )
    //      which is contained in (-2^{3n+2*LOGK+1}, 2^{3n+2*LOGK+1 +1}) // extra +1 just to be safe..
   
}

// input: in[0] + in[1] u
// output: (p-in[0]) + (p-in[1]) u
// assume 0 <= in < p
template Fp2negate(n, k){
    signal input in[2][k]; 
    signal input p[k];
    signal output out[2][k];
    
    component neg0 = BigSub(n, k);
    component neg1 = BigSub(n, k);
    for(var i=0; i<k; i++){
        neg0.a[i] <== p[i];
        neg1.a[i] <== p[i];
        neg0.b[i] <== in[0][i];
        neg1.b[i] <== in[1][i];
    }
    for(var i=0; i<k; i++){
        out[0][i] <== neg0.out[i];
        out[1][i] <== neg1.out[i];
    }
}

// input: a0 + a1 u, b0 + b1 u
// output: (a0-b0) + (a1-b1)u
template Fp2subtract(n, k){
    signal input a[2][k];
    signal input b[2][k];
    signal input p[k];
    signal output out[2][k];
    
    component sub0 = BigSubModP(n, k);
    component sub1 = BigSubModP(n, k);
    for(var i=0; i<k; i++){
        sub0.a[i] <== a[0][i];
        sub0.b[i] <== b[0][i];
        sub1.a[i] <== a[1][i];
        sub1.b[i] <== b[1][i];
    }
    for(var i=0; i<k; i++){
        out[0][i] <== sub0.out[i];
        out[1][i] <== sub1.out[i];
    }
}

// Call Fp2invert_func to compute inverse
// Then check out * in = 1, out is an array of shorts
template Fp2invert(n, k, p){
    signal input in[2][k];
    signal output out[2][k];

    var inverse[2][100] = Fp2invert_func(n, k, in, p); // 2 x 100, only 2 x k relevant
    for (var i = 0; i < 2; i ++) {
        for (var j = 0; j < k; j ++) {
            out[i][j] <-- inverse[i][j];
        }
    }

    //range checks
    component outRangeChecks[2][k];
    for(var i=0; i<2; i++) for(var j=0; j<k; j++){
        outRangeChecks[i][j] = Num2Bits(n);
        outRangeChecks[i][j].in <== out[i][j];
    }

    component in_out = Fp2multiply(n, k, p);
    for(var i=0; i<2; i++)for(var j=0; j<k; j++){
        in_out.a[i][j] <== in[i][j];
        in_out.b[i][j] <== out[i][j];
    }

    for(var i=0; i<2; i++)for(var j=0; j<k; j++){
        if(i == 0 && j == 0)
            in_out.out[i][j] === 1;
        else
            in_out.out[i][j] === 0;
    }
}

// a, b are two elements of Fp2 where we use the 4 x k format that remembers negatives 
// solve for out * b - a = p * X 
// assume X has at most m registers, lying in [-2^{n+1}, 2^{n+1} )   NOTE n+1 not n DIFFERENT FROM Fp2CarryModP
// assume a has registers in [0, 2^overflow), b has registers in [0, 2^{overflow-2n-2*LOGK-1}) 
// assume 2n+1 + log( max(k, ceil(overflow/n) )) < overflow 
// out has registers in [0, 2^n) 
template Fp2invertCarryModP(n, k, overflow, p){
    signal input a[4][k];
    signal input b[4][k]; 
    var m = (overflow + n - 1) \ n; // ceil(overflow/n) , i.e., 2^overflow <= 2^{n*m}
    signal output X[2][m];
    signal output out[2][k]; 

    
    // first precompute a, b mod p as shorts 
    var a_mod[2][100]; 
    var b_mod[2][100]; 
    for(var eps=0; eps<2; eps++){
        // 2^{overflow} <= 2^{n*ceil(overflow/n)} 
        var temp1[2][100] = long_div2(n, k, (overflow+n-1) \ n, long_to_short(n, k, a[2*eps]),p);
        var temp2[2][100] = long_div2(n, k, (overflow+n-1) \ n, long_to_short(n, k, a[2*eps+1]),p);
        a_mod[eps] = long_sub_mod(n,k,temp1[1],temp2[1],p);

        temp1 = long_div2(n, k, (overflow+n-1) \ n, long_to_short(n, k, b[2*eps]),p);
        temp2 = long_div2(n, k, (overflow+n-1) \ n, long_to_short(n, k, b[2*eps+1]),p);
        b_mod[eps] = long_sub_mod(n,k,temp1[1],temp2[1],p);
    }

    // precompute 1/b 
    var b_inv[2][100] = Fp2invert_func(n, k, b_mod, p);
    // precompute a/b
    var out_var[2][100] = find_Fp2_product(n, k, a_mod, b_inv, p);

    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++)
        out[eps][i] <-- out_var[eps][i]; 

    component check = checkValidFp2(n, k, p);
    for(var eps=0; eps<2; eps++)for(var i=0; i<k; i++)
        check.in[eps][i] <== out[eps][i];
    
    // constraint is a out * b = a + p * X 
    // precompute out * b = p * X' + Y' and a = p * X'' + Y''
    //            should have Y' = Y'' so X = X' - X''
    
    // out * b, registers overflow in 2*(k+1)*k * 2^{2n + (overflow - 2n - 2*LOGK - 1)} <= 2^{overflow}  
    component mult = Fp2multiplyNoCarryCompress(n, k, p); 
    for(var i=0; i<k; i++)for(var eps=0; eps<2; eps++){
        mult.a[2*eps][i] <== out[eps][i]; 
        mult.a[2*eps+1][i] <== 0;
        mult.b[2*eps][i] <== b[2*eps][i]; 
        mult.b[2*eps+1][i] <== b[2*eps+1][i];
    }
    
    // get mult = out * b = p*X' + Y'
    var XY[2][2][100] = Fp2_long_div(n, k, m, mult.out, p); // total value is < 2^{n*k+overflow} <= 2^{n*(k+m)} so m extra registers is enough
    // get a = p*X' + Y'
    var XY1[2][2][100] = Fp2_long_div(n, k, m, a, p); // same as above, m extra registers enough

    component X_range_checks[2][m];
    for(var eps=0; eps<2; eps++){    
        for(var i=0; i<m; i++){
            // X'' = X-X'
            X[eps][i] <-- XY[eps][0][i] - XY1[eps][0][i]; // each XY[eps][0] is in [-2^n, 2^n) so difference is in [-2^{n+1}, 2^{n+1})
            X_range_checks[eps][i] = Num2Bits(n+2);
            X_range_checks[eps][i].in <== X[eps][i] + (1<<(n+1)); // X[eps][i] should be between [-2^{n+1}, 2^{n+1})
        }
    }

    // finally constrain out * b - a = p * X 
    // out * b - a has overflow in (-2^{overflow+1}, 2^{overflow +1}) 
    // assume n+1 < overflow - n - log(max(k,m)+1), for registers of X
    component mod_check[2];  // overflow 9*(k+1)*k * 2^{3n+1} + 2*2^n < 2^{3n+LOGK+5} 
    for(var eps=0; eps<2; eps++){
        mod_check[eps] = checkBigMod(n, k, m, overflow + 1, p);
        for(var i=0; i<k; i++){
            mod_check[eps].in[i] <== mult.out[2*eps][i] - mult.out[2*eps+1][i] - a[2*eps][i] + a[2*eps+1][i];
            mod_check[eps].Y[i] <== 0;
        }
        for(var i=0; i<m; i++){
            mod_check[eps].X[i] <== X[eps][i];
        }
    }
    
}

// input: a+b u
// output: a-b u 
// IF p = 3 mod 4 THEN a - b u = (a+b u)^p <-- Frobenius map 
// aka Fp2frobeniusMap(n, k)
template Fp2conjugate(n, k){
    signal input in[2][k]; 
    signal input p[k];
    signal output out[2][k];
    
    component neg1 = BigSub(n, k);
    for(var i=0; i<k; i++){
        neg1.a[i] <== p[i];
        neg1.b[i] <== in[1][i];
    }
    for(var i=0; i<k; i++){
        out[0][i] <== in[0][i];
        out[1][i] <== neg1.out[i];
    }
}

// raises to q^power-th power 
template Fp2frobeniusMap(n, k, power){
    signal input in[2][k];
    signal input p[k];
    signal output out[2][k];
    
    var pow = power % 2;
    component neg1 = BigSub(n,k);
    if(pow == 0){
        for(var i=0; i<k; i++){
            out[0][i] <== in[0][i];
            out[1][i] <== in[1][i];
        }
    }else{
        for(var i=0; i<k; i++){
            neg1.a[i] <== p[i];
            neg1.b[i] <== in[1][i];
        }
        for(var i=0; i<k; i++){
            out[0][i] <== in[0][i];
            out[1][i] <== neg1.out[i];
        }
    }
}


