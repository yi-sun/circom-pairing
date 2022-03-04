pragma circom 2.0.2;

include "bigint.circom";

// Extra circuits with different approaches to field operations. 
// Kept here for reference; less efficient or specialized than the ones in field_elements.circom


// a, b are elements of Fp^l
// a[i] represents a[i][0] + a[i][1] * 2**n + ... + a[i][l-1] * 2**(n*(k-1))
// compute a+b in Fp^l
template FieldAdd2D(n, k, l) {
    signal input a[l][k];
    signal input b[l][k];
    signal input p[k];
    signal output c[l][k];

    component adders[l];
    for (var i = 0; i < l; i++) {
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

// take a polynomial expression a[0] + omega^1 a[1] + ... + omega^(2l-2) a[2l-2]
// reduce it to degree l-1 using omega^l + omega^(l-1) poly[l-1] + ... + poly[0] = 0
// WARNING: can produce incorrectly handled negative coefficients. only here for reference; do not use
template PolynomialReduce(l) {
    signal input a[2*l-1];
    signal input poly[l];
    signal output out[l];

    var residue[2*l-1];
    signal quotient[l-1];
    for (var i = 0; i < 2*l-1; i++) {
        residue[i] = a[i];
    }
    for (var i = l-2; i >= 0; i --) {
        for (var j = 0; j < l; j ++) {
            residue[i + j] = residue[i + j] - residue[i + l] * poly[j];
        }
        quotient[i] <-- residue[i+l];
        residue[i+l] = 0;
    }
    component mult = BigMultShortLong(1, l+1);
    for (var i = 0; i < l-1; i ++) {
        mult.a[i] <== quotient[i];
        mult.b[i] <== poly[i];
    }
    mult.a[l-1] <== 0;
    mult.a[l] <== 0;
    mult.b[l-1] <== poly[l-1];
    mult.b[l] <== 1;
    signal a_out[2*l-1];
    for (var i = 0; i < 2*l-1; i++) {
        a_out[i] <== mult.out[i];
    }
    for (var i = 0; i < l; i ++ ) {
        out[i] <-- residue[i];
    }
    for (var i = 0; i < l; i ++) {
        a[i] === a_out[i] + out[i];
    }
    for (var i = l; i < 2*l-1; i ++) {
        a[i] === a_out[i];
    }
}

template Fp2PolynomialReduce(n, k) {
    var l = 2;
    signal input a[2*l-1][k];
    var poly[2] = [1, 0]; // x^2 + 1 = 0
    signal output out[l][k];
    signal input p[k];

    for (var i = 0; i < k; i ++) {
        out[1][i] <== a[1][i];
    }
    component sub = BigSubModP(n, k);
    for (var i = 0; i < k; i ++) {
        sub.a[i] <== a[0][i];
        sub.b[i] <== a[2][i];
        sub.p[i] <== p[i];
    }
    for (var i = 0; i < k; i ++) {
        out[0][i] <== sub.out[i];
    }
}



// A similar circuit can do multiplication in different fields. 
// The only difference is that Fp2PolynomialReduce (which reduces quadratics by x^2+1) 
// must be replaced with a different circuit specialized to the minimal polynomial
template Fp2Multiply(n, k) {
    // l is always 2. poly is always [1, 0]
    var l = 2;
    signal input a[l][k];
    signal input b[l][k];
    signal input p[k];
    signal output c[l][k];

    component mult = BigMultShortLong2D(n, k, l);
    for (var i = 0; i < l; i ++) {
        for (var j = 0; j < k; j ++) {
            mult.a[i][j] <== a[i][j];
            mult.b[i][j] <== b[i][j];
        }
    } // out: 2l-1 x 2k-1 array of longs
    component longshorts[2*l-1];
    for (var i = 0; i < 2*l-1; i++) {
        longshorts[i] = LongToShortNoEndCarry(n, 2*k-1);
        for (var j = 0; j < 2*k-1; j ++) {
            longshorts[i].in[j] <== mult.out[i][j];
        }
    } // out: 2l-1 x 2k array of shorts
    component bigmods[2*l-1];
    for (var i = 0; i < 2*l-1; i ++) {
        bigmods[i] = BigMod(n, k);
        for (var j = 0; j < 2*k; j ++) {
            bigmods[i].a[j] <== longshorts[i].out[j];
        }
        for (var j = 0; j < k; j ++) {
            bigmods[i].b[j] <== p[j];
        }
    } // out: 2l-1 x k array of shorts
    component reduce = Fp2PolynomialReduce(n, k);
    for (var i = 0; i < 2*l-1; i ++) {
        for (var j = 0; j < k; j ++) {
            reduce.a[i][j] <== bigmods[i].mod[j];
        }
    } // out: l x k array of shorts
    for (var j = 0; j < k; j ++) {
        reduce.p[j] <== p[j];
    }
    for (var i = 0; i < l; i++) {
        for (var j = 0; j < k; j++) {
            c[i][j] <== reduce.out[i][j];
        }
    } // out: l x k array of shorts
}

// squaring can be optimized to save 2 multiplication
// (a**2-b**2) = (a+b)(a-b) 
// (a+b u)**2 = (a+b)(a-b) + (a*b+a*b)u
template Fp2square(n, k){
    signal input in[2][k];
    signal input p[k];
    signal output out[2][k];
    
    component sum = BigAdd(n, k);
    for(var i=0; i<k; i++){
        sum.a[i] <== in[0][i];
        sum.b[i] <== in[1][i];
    }
    component diff = BigSubModP(n, k);
    for(var i=0; i<k; i++){
        diff.a[i] <== in[0][i];
        diff.b[i] <== in[1][i];
        diff.p[i] <== p[i];
    }
    component prod = BigMult(n, k+1);
    for(var i=0; i<k; i++){
        prod.a[i] <== sum.out[i];
        prod.b[i] <== diff.out[i];
    }
    prod.a[k] <== sum.out[k];
    prod.b[k] <== 0;

    component prod_mod = BigMod2(n, k, 2*k+2);
    for(var i=0; i<2*k+2; i++){
        prod_mod.a[i] <== prod.out[i];
        if(i<k){
            prod_mod.b[i] <== p[i];
        }
    }
    for(var i=0; i<k; i++){
        out[0][i] <== prod_mod.mod[i];
    }
    
    component ab = BigMult(n, k);
    for(var i=0; i<k; i++){
        ab.a[i] <== in[0][i];
        ab.b[i] <== in[1][i];
    }
    component two_ab = BigAdd(n, 2*k); 
    for(var i=0; i<2*k; i++){
        two_ab.a[i] <== ab.out[i];
        two_ab.b[i] <== ab.out[i];
    }
    component two_ab_mod = BigMod2(n, k, 2*k+1);
    for(var i=0; i<2*k+1; i++){
        two_ab_mod.a[i] <== two_ab.out[i];
        if(i < k){
            two_ab_mod.b[i] <== p[i];
        }
    }
    for(var i=0; i<k; i++){
        out[1][i] <== two_ab_mod.mod[i];
    }
}

// a[2][k] all registers in [0, 2**n)
// b[2][k] all registers in [0, 2**n)
// p[k]
// consider a,b as elements of Fp2 
// out[2][2][k] solving
//      a0*b0 + (p-a1)*b1 = p * out[0][0] + out[0][1] with out[0][1] in [0,p) 
//      a0*b1 + a1*b0 = p * out[1][0] + out[1][1] with out[1][1] in [0,p) 
// out[i][0] has k+2 registers in short BigInt format [0, 2**n)
// out[i][1] has k registers in short BigInt format
// a * b = out[0][1] + out[1][1] * u in Fp2 
function Fp2prod(n, k, a, b, p){
    var out[2][2][20];
    // solve for X and Y such that a0*b0 + (p-a1)*b1 = p*X + Y with Y in [0,p) 
    // -a1*b1 = (p-a1)*b1 mod p
    var a0b0_var[20] = prod(n, k, a[0], b[0]);
    var a1_neg[20] = long_sub(n, k, p, a[1]); 
    var a1b1_neg[20] = prod(n, k, a1_neg, b[1]);
    var diff[20] = long_add(n, 2*k, a0b0_var, a1b1_neg); // 2*k+1 registers
    out[0] = long_div2(n, k, k+1, diff, p); 
    // X = out[0][0] has k+2 registers, Y = out[0][1] has k registers 
    
    // solve for X and Y such that a0*b1 + a1*b0 = p*X + Y with Y in [0,p) 
    var a0b1_var[20] = prod(n, k, a[0], b[1]);
    var a1b0_var[20] = prod(n, k, a[1], b[0]);
    var sum[20] = long_add(n, 2*k, a0b1_var, a1b0_var); // output 2*k+1 registers
    out[1] = long_div2(n, k, k+1, sum, p); 
    // X = out[1][0] has k+2 registers, Y = out[1][1] has k registers 

    return out;
}

// output (a0 + a1 u)*(b0 + b1 u) = (a0*b0 + (p-a1)*b1) + (a0*b1 + a1*b0)u 
//      where no carries are performed 
// out[0], out[1] have 2*k-1 registers 
// out[0][i] in (-(k+1)*2^{2n}, (k+1)*2^{2n+1})
// out[1][i] in [0, (k+1)*2^{2n+1}) 
template Fp2multiplyNoCarry(n, k){
    signal input a[2][k];
    signal input b[2][k];
    signal input p[k];
    signal output out[2][2*k-1];
    
    component a0b0 = BigMultShortLong(n, k);
    component a1b1 = BigMultShortLong(n, k);
    component pb1 = BigMultShortLong(n, k); 
    component a0b1 = BigMultShortLong(n, k); 
    component a1b0 = BigMultShortLong(n, k);
    
    for(var i=0; i<k; i++){
        a0b0.a[i] <== a[0][i];
        a0b0.b[i] <== b[0][i];

        a1b1.a[i] <== a[1][i];
        a1b1.b[i] <== b[1][i];

        pb1.a[i] <== p[i];
        pb1.b[i] <== b[1][i];

        a0b1.a[i] <== a[0][i];
        a0b1.b[i] <== b[1][i];

        a1b0.a[i] <== a[1][i];
        a1b0.b[i] <== b[0][i];
    }
 
    for(var i=0; i<2*k-1; i++){
        out[0][i] <== a0b0.out[i] + pb1.out[i] - a1b1.out[i];
        out[1][i] <== a0b1.out[i] + a1b0.out[i];
    }
}

// multiplication specialized to Fp^2 
// (a0 + a1 u)*(b0 + b1 u) = (a0*b0 - a1*b1) + (a0*b1 + a1*b0)u
template Fp2multiply(n, k){
    signal input a[2][k];
    signal input b[2][k];
    signal input p[k];
    signal output out[2][k];

    var LOGK = 3;
    assert(k<7);
    assert(2*n + 1 + LOGK<254);

    var Xvar[2][2][20] = Fp2prod(n, k, a, b, p); 
    component range_checks[2][k];
    component lt[2];
    signal X[2][k+2]; 
    component X_range_checks[2][k+2];
    
    for(var eps=0; eps<2; eps++){
        lt[eps] = BigLessThan(n, k);
        for(var i=0; i<k; i++){
            out[eps][i] <-- Xvar[eps][1][i];
            range_checks[eps][i] = Num2Bits(n);
            range_checks[eps][i].in <== out[eps][i];
            
            lt[eps].a[i] <== out[eps][i];
            lt[eps].b[i] <== p[i];
        }
        lt[eps].out === 1;
        
        for(var i=0; i<k+2; i++){
            X[eps][i] <-- Xvar[eps][0][i];
            X_range_checks[eps][i] = Num2Bits(n);
            X_range_checks[eps][i].in <== X[eps][i];
        }
        
    }

    
    // out[0] constraint: X = X[0], Y = out[0] 
    // constrain by Carry( a0 *' b0 +' p *' b1 -' a1 *' b1 - p *' X - Y ) = 0 
    // where all operations are performed without carry 
    // each register is an overflow representation in the range (-(k+1)*2^{2n+1}-2^n, (k+1)*2^{2n + 1} )
    //                                          which is inside (-2^{2n+1+LOGK}, 2^{2n+1+LOGK})

    // out[1] constraint: X = X[1], Y = out[1]
    // constrain by Carry( a0 *' b1 +' a1 *' b0 -' p *' X - Y) = 0 
    // each register is an overflow representation in the range (-(k+1)*2^{2n}-2^n, (k+1)*2^{2n + 1} )
    //                                          which is inside (-2^{2n+1+LOGK}, 2^{2n+1+LOGK})
    
    component ab = Fp2multiplyNoCarry(n, k); 
    for(var i=0; i<k; i++){
        ab.p[i] <== p[i];
        ab.a[0][i] <== a[0][i];
        ab.a[1][i] <== a[1][i];
        ab.b[0][i] <== b[0][i];
        ab.b[1][i] <== b[1][i];
    }
    component pX[2];
    component carry_check[2];
    for(var eps=0; eps<2; eps++){
        pX[eps] = BigMultShortLong(n, k+2); // 2*k+3 registers
        for(var i=0; i<k; i++){
            pX[eps].a[i] <== p[i];
            pX[eps].b[i] <== X[eps][i];
        }
        for(var i=k; i<k+2; i++){
            pX[eps].a[i] <== 0;
            pX[eps].b[i] <== X[eps][i];
        }

        carry_check[eps] = CheckCarryToZero(n, 2*n+2+LOGK, 2*k+3); 
        for(var i=0; i<k; i++)
            carry_check[eps].in[i] <== ab.out[eps][i] - pX[eps].out[i] - out[eps][i]; 
        for(var i=k; i<2*k-1; i++)
            carry_check[eps].in[i] <== ab.out[eps][i] - pX[eps].out[i]; 
        for(var i=2*k-1; i<2*k+3; i++)
            carry_check[eps].in[i] <== -pX[eps].out[i];
    }

}
