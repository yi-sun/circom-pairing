import path = require("path");
const circom_tester = require('circom_tester');
const wasm_tester = circom_tester.wasm;
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

function Fp2_to_array(n: number, k: number, x: Fp2){
    let {c0, c1} = x;
//    console.log(c0.value);
//    console.log(c1.value);
    return [ bigint_to_array(n, k, c0.value), bigint_to_array(n, k, c1.value) ];
}

function Fp12_to_array(n: number, k: number, x: Fp12){
    let {c0, c1} = x;
    let {c0: c00, c1: c01, c2: c02} = c0;
    let {c0: c10, c1: c11, c2: c12} = c1;
    return [ Fp2_to_array(n, k, c00), Fp2_to_array(n, k, c10), Fp2_to_array(n, k, c01), Fp2_to_array(n, k, c11), Fp2_to_array(n,k,c02), Fp2_to_array(n,k,c12) ];
}


describe("Fp12Multiply n = 3, k = 2", function() {
    this.timeout(1000 * 1000);

    // runs circom compilation
    let circuit: any;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_fp12_multiply_32.circom"));
    });

    // a0, a1, b0, b1, p, c0, c1
    var test_cases: Array<[bigint, bigint, bigint, bigint, bigint, bigint, bigint]> = [];
    let p: bigint = 19n;
    for (var a0 = 0n; a0 < p; a0 = a0 + 5n) {
        for (var b0 = 0n; b0 < p; b0 = b0 + 5n) {
            for (var a1 = 0n; a1 < p; a1 = a1 + 5n) {
                for (var b1 = 0n; b1 < p; b1 = b1 + 5n) {
                    var c0 = (a0 * b0 - a1 * b1 + p * p) % p;
                    var c1 = (a0 * b1 + a1 * b0) % p;
                    test_cases.push([a0, a1, b0, b1, p, c0, c1]);
                }
            }
        }
    }

    var test_field_multiply_32 = function (x: [bigint, bigint, bigint, bigint, bigint, bigint, bigint]) {
        const [a0, a1, b0, b1, p, c0, c1] = x;

        var a0_array: bigint[] = bigint_to_array(3, 2, a0);
        var a1_array: bigint[] = bigint_to_array(3, 2, a1);	
        var b0_array: bigint[] = bigint_to_array(3, 2, b0);
        var b1_array: bigint[] = bigint_to_array(3, 2, b1);
	var p_array: bigint[] = bigint_to_array(3, 2, p);
        var c0_array: bigint[] = bigint_to_array(3, 2, c0);
        var c1_array: bigint[] = bigint_to_array(3, 2, c1);
        var c0i_array: bigint[] = bigint_to_array(3, 2, (c0-c1 + p) % p);
        var c1i_array: bigint[] = bigint_to_array(3, 2, (c0+c1) % p);
        var zero: bigint[] = bigint_to_array(3, 2, 0n)

        it('Testing a0: ' + a0 + ' a1: ' + a1 + ' b0: ' + b0 + ' b1: ' + b1 + ' p: ' + p + ' c0: ' + c0 + ' c1: ' + c1, async function() {
            let witness = await circuit.calculateWitness({"a": [[zero, zero], [zero, zero], [zero, zero], [a0_array, a1_array], [zero, zero], [a0_array, a1_array]], 
            "b": [[b0_array, b1_array], [zero, zero], [zero, zero], [b0_array, b1_array], [zero, zero], [zero, zero]]});
	    await circuit.assertOut(witness, {"out": [[c0i_array, c1i_array], [zero, zero], [c0i_array, c1i_array ], 
    [c0_array, c1_array], [zero, zero], [c0_array, c1_array]]});
            await circuit.checkConstraints(witness);
        });
    }

    test_cases.forEach(test_field_multiply_32);
});

/*
describe("Fp12Multiply2 n = 3, k = 2", function() {
    this.timeout(1000 * 1000);

    // runs circom compilation
    let circuit: any;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_fp12_multiply2_32.circom"));
    });

    // a0, a1, b0, b1, c0, c1
    let p: bigint = 13n;
    var test_cases: Array<[bigint, bigint, bigint, bigint, bigint, bigint]> = [];
    for (var a0 = 0n; a0 < p; a0 = a0 + 5n) {
        for (var b0 = 0n; b0 < p; b0 = b0 + 5n) {
            for (var a1 = 0n; a1 < p; a1 = a1 + 5n) {
                for (var b1 = 0n; b1 < p; b1 = b1 + 5n) {
                    var c0 = (a0 * b0 - a1 * b1 + p * p) % p;
                    var c1 = (a0 * b1 + a1 * b0) % p;
                    test_cases.push([a0, a1, b0, b1, c0, c1]);
                }
            }
        }
    }

    var test_field_multiply_32 = function (x: [bigint, bigint, bigint, bigint, bigint, bigint]) {
        const [a0, a1, b0, b1, c0, c1] = x;

        var a0_array: bigint[] = bigint_to_array(3, 2, a0);
        var a1_array: bigint[] = bigint_to_array(3, 2, a1);	
        var b0_array: bigint[] = bigint_to_array(3, 2, b0);
        var b1_array: bigint[] = bigint_to_array(3, 2, b1);
        var c0_array: bigint[] = bigint_to_array(3, 2, c0);
        var c1_array: bigint[] = bigint_to_array(3, 2, c1);
        var c0i_array: bigint[] = bigint_to_array(3, 2, (c0-c1 + p) % p);
        var c1i_array: bigint[] = bigint_to_array(3, 2, (c0+c1) % p);
        var zero: bigint[] = bigint_to_array(3, 2, 0n)

        it('Testing a0: ' + a0 + ' a1: ' + a1 + ' b0: ' + b0 + ' b1: ' + b1 + ' p: ' + p + ' c0: ' + c0 + ' c1: ' + c1, async function() {
            let witness = await circuit.calculateWitness({"a": [[zero, zero], [zero, zero], [zero, zero], [a0_array, a1_array], [zero, zero], [a0_array, a1_array]], 
	    		                                  "b": [[b0_array, b1_array], [zero, zero], [zero, zero], [b0_array, b1_array], [zero, zero], [zero, zero]]});
	    await circuit.assertOut(witness, {"out": [[c0i_array, c1i_array], [zero, zero], [c0i_array, c1i_array ], 
                                                      [c0_array, c1_array], [zero, zero], [c0_array, c1_array]]});  
            await circuit.checkConstraints(witness);
        });
    }

    test_cases.forEach(test_field_multiply_32);
});
*/

describe("Fp12Compression n = 3, k = 2", function() {
    this.timeout(1000 * 1000);

    // runs circom compilation
    let circuit: any;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_fp12_compression_32.circom"));
    });

    let p: bigint = 19n;
    Fp.ORDER = p;
    Fp2.ORDER = p;

    var test_compression_32 = function (x: Fp12) {
        let in_array: Array<[ bigint[], bigint[] ]> = new Array(6);
        let a_array: Array< [bigint, bigint] > = new Array(6);
        let {c0: x0, c1: x1} = x;
        let {c0: c00, c1: c01, c2: c02} = x0;
        let {c0: c10, c1: c11, c2: c12} = x1;
        let Fp2_array: Array< Fp2 > = [c00, c10, c01, c11, c02, c12]; 
        for( var i = 0; i < 6; i++ ){
            let {c0: a0, c1: a1} = Fp2_array[i];
            in_array[i] = [ bigint_to_array(3, 2, a0.value), bigint_to_array(3, 2, a1.value) ];
            a_array[i] = [ a0.value, a1.value ];
        }

        it('Testing a0: ' + a_array[0][0] + ', ' + a_array[0][1] + 
            ' a1: ' + a_array[1][0] + ', ' + a_array[1][1] + 
            ' a2: ' + a_array[2][0] + ', ' + a_array[2][1] + 
            ' a3: ' + a_array[3][0] + ', ' + a_array[3][1] + 
            ' a4: ' + a_array[4][0] + ', ' + a_array[4][1] + 
            ' a5: ' + a_array[5][0] + ', ' + a_array[5][1] + 
            ' p: ' + p, async function() {
            let witness = await circuit.calculateWitness({"in": in_array });
	    await circuit.assertOut(witness, {"out": in_array });
            await circuit.checkConstraints(witness);
        });
    }

    var test_cases: Array<Fp12> = [];
    for(var test_id = 0; test_id < 100; test_id++){
        let rand_twelve: BigintTwelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];
        for( let i = 0; i < 12; i++){
            rand_twelve[i] = BigInt( Math.floor(Math.random() * 19) );
        }
        let elt: Fp12 = Fp12.fromBigTwelve( rand_twelve ); 
        let elt2: Fp12 = elt.pow( ( (p ** 6n) - 1n) * (p*p + 1n) ); // this is now in cyclotomic subgroup 
        if( !elt2.equals(Fp12.ONE) ){ // skip if elt2 = 1
            // console.log( elt2.pow( (p ** 4n) - (p*p) + 1n ).equals(Fp12.ONE) );
            test_cases.push(elt2);
        }
    }

    test_cases.forEach(test_compression_32);

});


describe("Fp12Invert n = 4, k = 2", function() {
    this.timeout(1000 * 1000);

    // runs circom compilation
    let circuit: any;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_fp12_invert_42.circom"));
    });

    // a0, a1, b0, b1, p, c0, c1
    var test_cases: Array<number> = [2];

    var test_field_invert_42 = function (x: number) {
        var zero: bigint[] = bigint_to_array(4, 2, 0n);
        let one: bigint[] = bigint_to_array(4, 2, 1n);
        let three: bigint[] = bigint_to_array(4, 2, 3n);
        let four: bigint[] = bigint_to_array(4, 2, 4n);
        let seven: bigint[] = bigint_to_array(4, 2, 7n);
        let eight: bigint[] = bigint_to_array(4, 2, 8n);
        let nine: bigint[] = bigint_to_array(4, 2, 9n);
        let eleven: bigint[] = bigint_to_array(4, 2, 11n);
        let fourteen: bigint[] = bigint_to_array(4, 2, 14n);
        let thirteen: bigint[] = bigint_to_array(4, 2, 13n);

        let input: bigint[][][];
        let output: bigint[][][];
        if (x == 0) {
            input = [[one, zero], [zero, zero], [zero, zero], [zero, zero], [zero, zero], [zero, zero]];
            output = [[one, zero], [zero, zero], [zero, zero], [zero, zero], [zero, zero], [zero, zero]];
        }
        if (x == 1) {
            input = [[zero, zero], [zero, zero], [one, zero], [zero, zero], [zero, zero], [zero, zero]];
            output = [[zero, zero], [zero, zero], [zero, zero], [zero, zero], [nine, eight], [zero, zero]];
        }
        if (x == 2) {
            input = [[zero, zero], [zero, zero], [one, zero], [zero, one], [zero, one], [zero, zero]];
            output = [[fourteen, eleven], [seven, seven], [one, thirteen], [three, nine], [fourteen, four], [thirteen, zero]];
        }

        it('Testing case: ' + x, async function() {
            let witness = await circuit.calculateWitness({"in": input});
	    await circuit.assertOut(witness, {"out": output});
            await circuit.checkConstraints(witness);
        });
    }

    test_cases.forEach(test_field_invert_42);
});


describe("Fp12cyclotomicSquare n = 3, k = 2", function() {
    this.timeout(1000 * 1000);

    // runs circom compilation
    let circuit: any;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_fp12_cyclotomicSquare_32.circom"));
    });

    let p: bigint = 19n;
    Fp.ORDER = p;
    Fp2.ORDER = p;

    var test_cyclosq_32 = function (g: [Fp2, Fp2, Fp2, Fp2]) {
        const [g2, g3, g4, g5] = g;

        let in_array: bigint[][][] = [ Fp2_to_array(3, 2, g2), Fp2_to_array(3, 2, g3), Fp2_to_array(3, 2, g4), Fp2_to_array(3, 2, g5) ];
        
        let c: Fp2 = Fp2.fromBigTuple([1n, 1n]);
        let A23: Fp2 = (g2.add(g3)).multiply(g2.add(c.multiply(g3))); 
        let A45: Fp2 = (g4.add(g5)).multiply(g4.add(c.multiply(g5))); 
        let B23: Fp2 = g2.multiply(g3);
        let B45: Fp2 = g4.multiply(g5);
        let h2: Fp2 = (g2.add(c.multiply(B45).multiply(3n))).multiply(2n);
        let h3: Fp2 = (A45.subtract( (c.add(Fp2.ONE)).multiply(B45) )).multiply(3n).subtract(g3.multiply(2n)); 
        let h4: Fp2 = (A23.subtract( (c.add(Fp2.ONE)).multiply(B23) )).multiply(3n).subtract(g4.multiply(2n)); 
        let h5: Fp2 = (g5.add(B23.multiply(3n))).multiply(2n);

        let out_array: bigint[][][] = [ Fp2_to_array(3, 2, h2), Fp2_to_array(3, 2, h3), Fp2_to_array(3, 2, h4), Fp2_to_array(3, 2, h5) ];
        
        let a_array: Array< [bigint, bigint] > = new Array(4);
        for( var i = 0; i < 4; i++ ){
            let {c0: a0, c1: a1} = g[i];
            a_array[i] = [ a0.value, a1.value ];
        }

        it('Testing: in: ' + in_array + ', out: ' + out_array + 
            ' p: ' + p, async function() {
            let witness = await circuit.calculateWitness({"in": in_array });
	    await circuit.assertOut(witness, {"out": out_array });
            await circuit.checkConstraints(witness);
        });
    }

    var test_cases: Array<[Fp2, Fp2, Fp2, Fp2]> = [];
    for(var test_id = 0; test_id < 100; test_id++){
        let rand_g: [Fp2, Fp2, Fp2, Fp2] = [ Fp2.ONE, Fp2.ONE, Fp2.ONE, Fp2.ONE ];
        for( let i = 0; i < 4; i++){
            rand_g[i] = Fp2.fromBigTuple([ BigInt(Math.floor(Math.random() * 19)), BigInt(Math.floor(Math.random() * 19))]);
        }
        test_cases.push(rand_g);
    }

    test_cases.forEach(test_cyclosq_32);

});


describe("Fp12cyclotomicExp n = 3, k = 2", function() {
    this.timeout(1000 * 1000);

    // runs circom compilation
    let circuit: any;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_fp12_cyclotomicExp_32.circom"));
    });

    let p: bigint = 19n;
    Fp.ORDER = p;
    Fp2.ORDER = p;

    var test_cycloexp_32 = function (test: [Fp12, bigint]) {
        const [x, e] = test;

        let in_array: bigint[][][] = Fp12_to_array(3, 2, x);
        let out_array: bigint[][][] = Fp12_to_array(3, 2, x.pow(e) );

        let a_array: Array< [bigint, bigint] > = new Array(6);
        let {c0: x0, c1: x1} = x;
        let {c0: c00, c1: c01, c2: c02} = x0;
        let {c0: c10, c1: c11, c2: c12} = x1;
        let Fp2_array: Array< Fp2 > = [c00, c10, c01, c11, c02, c12]; 
        for( var i = 0; i < 6; i++ ){
            let {c0: a0, c1: a1} = Fp2_array[i];
            a_array[i] = [ a0.value, a1.value ];
        }

        it('Testing a0: ' + a_array[0][0] + ', ' + a_array[0][1] + 
            ' a1: ' + a_array[1][0] + ', ' + a_array[1][1] + 
            ' a2: ' + a_array[2][0] + ', ' + a_array[2][1] + 
            ' a3: ' + a_array[3][0] + ', ' + a_array[3][1] + 
            ' a4: ' + a_array[4][0] + ', ' + a_array[4][1] + 
            ' a5: ' + a_array[5][0] + ', ' + a_array[5][1] + 
            ' e: ' + e + 
            ' p: ' + p, async function() {
            let witness = await circuit.calculateWitness({"in": in_array });
	    await circuit.assertOut(witness, {"out": out_array });
            await circuit.checkConstraints(witness);
        });
    }

    var test_cases: Array<[Fp12,bigint]> = [];
    for(var test_id = 0; test_id < 20; test_id++){
        let rand_twelve: BigintTwelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];
        for( let i = 0; i < 12; i++){
            rand_twelve[i] = BigInt( Math.floor(Math.random() * 19) );
        }
        let elt: Fp12 = Fp12.fromBigTwelve( rand_twelve ); 
        let cyc: Fp12 = elt.pow( ( (p ** 6n) - 1n) * (p*p + 1n) ); // this is now in cyclotomic subgroup 
        if( !cyc.equals(Fp12.ONE) ){ // skip if elt2 = 1
            // console.log( elt2.pow( (p ** 4n) - (p*p) + 1n ).equals(Fp12.ONE) );
            test_cases.push([cyc, 5n]);
        }
    }

    test_cases.forEach(test_cycloexp_32);

}); 


