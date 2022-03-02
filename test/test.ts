import {
  Fp, Fp2, Fp6, Fp12, CURVE, mod
} from './math';

type BigintTuple = [bigint, bigint];
type FpTuple = [Fp, Fp];
type BigintSix = [bigint, bigint, bigint, bigint, bigint, bigint];
// prettier-ignore
type BigintTwelve = [
  bigint, bigint, bigint, bigint, bigint, bigint,
  bigint, bigint, bigint, bigint, bigint, bigint
];

function bigint_to_array(n: number, k: number, x: bigint) {
    let mod: bigint = 1n;
    for (var idx = 0; idx < n; idx++) {
        mod = mod * 2n;
    }

    let ret: bigint[] = [];
    var x_temp: bigint = x;
    for (var idx = 0; idx < k; idx++) {
        ret.push(x_temp % mod);
        x_temp = x_temp / mod;
    }
    return ret;
}

let p: bigint = 13n;
Fp.ORDER = p;
Fp2.ORDER = p;

for(var num = 2n; num < p ** 12n; num += 1n){
    let rand_twelve: BigintTwelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];
    let temp : bigint = num;
    for( let i = 0; i < 12; i++){
        rand_twelve[i] = temp % p;
        temp = temp / p;
    }
    let elt: Fp12 = Fp12.fromBigTwelve( rand_twelve ); 
    //let elt2: Fp12 = elt.pow( (p ** 6n - 1n) * (p*p + 1n) ); // this is now in cyclotomic subgroup 
    if( elt.pow( p**4n - p*p + 1n ).equals( Fp12.ONE ) ){
       console.log(rand_twelve);
       break;
    } // elt is in cyclotomic subgroup
    //if( !elt2.equals(Fp12.ONE) && elt2.pow( p ** 4n - p ** 2n + 1n ).equals( Fp12.ONE ) ) // check it's in cyclotomic subgroup
}

// let {c0, c1} = 
// console.log(Fp2.fromBigTuple([1n, 1n]).pow( (p*p-1n)/6n ).equals( Fp2.ONE ) );
// console.log(c0.value + ', ' + c1.value);
