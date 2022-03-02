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

let p: bigint = 19n;
Fp.ORDER = p;
Fp2.ORDER = p;

function printFp2(x: Fp2){
    let {c0, c1} = x;
    return [c0.value, c1.value];
}

function printFp12(x: Fp12){
    let {c0, c1} = x;
    let {c0: c00, c1: c01, c2: c02} = c0;
    let {c0: c10, c1: c11, c2: c12} = c1;
    return [ printFp2(c00), printFp2(c10), printFp2(c01), printFp2(c11), printFp2(c02), printFp2(c12) ];
}

/*
while(1){
    let rand_twelve: BigintTwelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];
    for( let i = 0; i < 12; i++){
        rand_twelve[i] = BigInt( Math.floor(Math.random() * 19) );
    }
    let elt: Fp12 = Fp12.fromBigTwelve( rand_twelve ); 
    if(! elt.pow( p**12n - 1n ).equals(Fp12.ONE) ){
        console.log(p**12n - 1n);
        console.log(rand_twelve);
        console.log(printFp12( elt ));
        console.log(printFp12( elt.pow( p**12n - 1n ) ) );
        break;
    }
    /* let ord:bigint = 1n;
    let temp: Fp12 = elt;
    while(! temp.equals(Fp12.ONE) ){
        ord+=1n;
        temp = temp.multiply(elt);
    }
    console.log('order: ' + ord); 
    if(ord % (p**4n - p*p + 1n) == 0n){ */
/*
    let cyc: Fp12 = elt.pow( (p ** 6n - 1n) * (p*p + 1n) ); // this is now in cyclotomic subgroup 
    if( !cyc.equals(Fp12.ONE) ){ // check it's in cyclotomic subgroup
        console.log( elt.pow( p**12n - 1n ).equals(Fp12.ONE) );
        console.log( cyc.pow( p**4n - p*p + 1n ).equals(Fp12.ONE) );
        console.log( printFp12(cyc) );
        break;
    } 

}
*/

// let {c0, c1} = 
// console.log(Fp2.fromBigTuple([1n, 1n]).pow( (p*p-1n)/2n ).equals( Fp2.ONE ) );
// console.log(Fp2.fromBigTuple([1n, 1n]).pow( (p*p-1n)/3n ).equals( Fp2.ONE ) );
// console.log(c0.value + ', ' + c1.value);

// let tcase : bigint[] = [ 10n, 7n, 14n, 7n, 18n, 18n, 15n, 18n, 16n, 9n, 15n, 9n ];
// for(var i=0; i<12; i++ ){
//     console.log( bigint_to_array(3, 2, tcase[i]) );
// }
