pragma circom 2.0.2;

include "bigint.circom";
include "bigint_func.circom";


// a[i], b[i] in 0... 2**n-1
// represent a = a[0] + a[1] * 2**n + .. + a[k - 1] * 2**(n * k)
// calculates (a+b)%p, where 0<= a,b < p 
template FpAdd(n, k, p){
    assert(n <= 252);
    signal input a[k];
    signal input b[k];
    signal output out[k];

    component add = BigAdd(n,k);
    for (var i = 0; i < k; i++) {
        add.a[i] <== a[i];
        add.b[i] <== b[i];
    }
    component lt = BigLessThan(n, k+1);
    for (var i = 0; i < k; i++) {
        lt.a[i] <== add.out[i];
    }
    lt.a[k] <== add.out[k];
    lt.b[k] <== 0; 

    component sub = BigSub(n,k+1);
    for (var i = 0; i < k; i++) {
        sub.a[i] <== add.out[i];
        sub.b[i] <== (1-lt.out) * p[i];
    }
    sub.a[k] <== add.out[k];
    sub.b[k] <== 0;
    
    sub.out[k] === 0;
    for (var i = 0; i < k; i++) {
        out[i] <== sub.out[i];
    }
}

// calculates (a - b) % p, where a, b < p
// note: does not assume a >= b
template FpSubtract(n, k, p){
    assert(n <= 252);
    signal input a[k];
    signal input b[k];
    signal output out[k];
    component sub = BigSub(n, k);
    for (var i = 0; i < k; i++){
        sub.a[i] <== a[i];
        sub.b[i] <== b[i];
    }
    signal flag;
    flag <== sub.underflow;
    component add = BigAdd(n, k);
    for (var i = 0; i < k; i++){
        add.a[i] <== sub.out[i];
        add.b[i] <== p[i];
    }
    signal tmp[k];
    for (var i = 0; i < k; i++){
        tmp[i] <== (1 - flag) * sub.out[i];
        out[i] <== tmp[i] + flag * add.out[i];
    }
}

template FpMultiply(n, k, p) {
    assert(n <= 252);
    signal input a[k];
    signal input b[k];
    signal output out[k];

    component big_mult = BigMult(n, k);
    for (var i = 0; i < k; i++) {
        big_mult.a[i] <== a[i];
        big_mult.b[i] <== b[i];
    }
    component big_mod = BigMod(n, k);
    for (var i = 0; i < 2 * k; i++) {
        big_mod.a[i] <== big_mult.out[i];
    }
    for (var i = 0; i < k; i++) {
        big_mod.b[i] <== p[i];
    }
    for (var i = 0; i < k; i++) {
        out[i] <== big_mod.mod[i];
    }
}

// constrain in = p * X + Y 
// in[i] in (-2^overflow, 2^overflow) 
// assume registers of X have abs value < 2^{overflow - n - log(max(k,m)+1)} 
template CheckCarryModP(n, k, m, overflow, p){
    signal input in[k]; 
    signal input X[m];
    signal input Y[k];

    assert( overflow < 253 );
    component pX;
    component carry_check;
    var maxkm;
    if(k < m) maxkm = m;
    else maxkm = k;

    pX = BigMultShortLong(n, maxkm); // p has k registers, X has m registers, so output really has k+m-1 registers 
    // overflow register in  (-2^overflow , 2^overflow)
    for(var i=0; i<maxkm; i++){
        if(i < k)
            pX.a[i] <== p[i];
        else
            pX.a[i] <== 0;
        if(i < m)
            pX.b[i] <== X[i];
        else 
            pX.b[i] <== 0;
    }
    // in - p*X has registers in (-2^{overflow+1}, 2^{overflow+1})
    carry_check = CheckCarryToZero(n, overflow+2, k+m-1 ); 
    for(var i=0; i<k; i++){
        carry_check.in[i] <== in[i] - pX.out[i] - Y[i]; 
    }
    for(var i=k; i<k+m-1; i++)
        carry_check.in[i] <== -pX.out[i];
}

// solve for in0 - in1 = p * X + out
// assume in has registers in [0, 2^overflow) 
// X has registers lying in [-2^n, 2^n) 
// X has at most Ceil( overflow / n ) registers 
template FpCarryModP(n, k, overflow, p){
    signal input in[2][k]; 
    var m = (overflow + n - 1) \ n; 
    signal output X[m];
    signal output out[k];

    assert( overflow < 253 );

    var Xvar[2][20] = get_Fp_carry_witness(n, k, m, in, p); 
    component range_checks[k]; 
    component X_range_checks[m];
    component lt = BigLessThan(n, k);

    for(var i=0; i<k; i++){
        out[i] <-- Xvar[1][i];
        range_checks[i] = Num2Bits(n);
        range_checks[i].in <== out[i];

        lt.a[i] <== out[i];
        lt.b[i] <== p[i];
    }
    lt.out === 1;
    
    for(var i=0; i<m; i++){
        X[i] <-- Xvar[0][i];
        X_range_checks[i] = Num2Bits(n+1);
        X_range_checks[i].in <== X[i] + (1<<n); // X[i] should be between [-2^n, 2^n)
    }
    
    component mod_check = CheckCarryModP(n, k, m, overflow, p);
    for(var i=0; i<k; i++){
        mod_check.in[i] <== in[0][i] - in[1][i];
        mod_check.Y[i] <== out[i];
    }
    for(var i=0; i<m; i++){
        mod_check.X[i] <== X[i];
    }
}


// Constrain in0 - in1 = 0 mod p by solving for in0 - in1 = p * X
// assume in has registers in [0, 2^overflow) 
// X has registers lying in [-2^n, 2^n) 
// X has at most Ceil( overflow / n ) registers 

// saves a range checks and BigLessThan comparison compared to CarryModP
template CheckCarryModToZero(n, k, overflow, p){
    signal input in[2][k]; 
    var m = (overflow + n - 1) \ n; 
    signal output X[m];

    assert( overflow < 253 );

    var Xvar[2][20] = get_Fp_carry_witness(n, k, m, in, p); 
    component X_range_checks[m];

    for(var i=0; i<m; i++){
        X[i] <-- Xvar[0][i];
        X_range_checks[i] = Num2Bits(n+1);
        X_range_checks[i].in <== X[i] + (1<<n); // X[i] should be between [-2^n, 2^n)
    }
    
    component mod_check = CheckCarryModP(n, k, m, overflow, p);
    for(var i=0; i<k; i++){
        mod_check.in[i] <== in[0][i] - in[1][i];
        mod_check.Y[i] <== 0;
    }
    for(var i=0; i<m; i++){
        mod_check.X[i] <== X[i];
    }
}

