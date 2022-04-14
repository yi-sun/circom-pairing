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

describe("LineFunctionUnequal n = 3, k = 2", function() {
    this.timeout(1000 * 1000);

    // runs circom compilation
    let circuit: any;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_linefunc_unequal_32.circom"));
    });

    let p: bigint = 19n;
    Fp.ORDER = p;
    Fp2.ORDER = p;

    var test_linefunc_32 = function (inp: [Fp, Fp, Fp, Fp, Fp12, Fp12]) {
        const [P1x, P1y, P2x, P2y, Qx, Qy] = inp;

        let in_array: Array<[ bigint[], bigint[] ]> = new Array(2);
        in_array[0] = [ bigint_to_array(3, 2, P1x.value), bigint_to_array(3, 2, P1y.value) ];
        in_array[1] = [ bigint_to_array(3, 2, P2x.value), bigint_to_array(3, 2, P2y.value) ];

        let pointX: Array< [bigint[], bigint[]] > = new Array(6);
        let pointY: Array< [bigint[], bigint[]] > = new Array(6);	
        var {c0: x0, c1: x1} = Qx;
        var {c0: c00, c1: c01, c2: c02} = x0;
        var {c0: c10, c1: c11, c2: c12} = x1;
        let Fp2_arrayX: Array< Fp2 > = [c00, c10, c01, c11, c02, c12]; 
        for( var i = 0; i < 6; i++ ){
            var {c0: a0, c1: a1} = Fp2_arrayX[i];
            pointX[i] = [ bigint_to_array(3, 2, a0.value), bigint_to_array(3, 2, a1.value) ];
        }

        var {c0: y0, c1: y1} = Qy;
        var {c0: c00, c1: c01, c2: c02} = y0;
        var {c0: c10, c1: c11, c2: c12} = y1;
        let Fp2_arrayY: Array< Fp2 > = [c00, c10, c01, c11, c02, c12]; 
        for( var i = 0; i < 6; i++ ){
            var {c0: a0, c1: a1} = Fp2_arrayY[i];
            pointY[i] = [ bigint_to_array(3, 2, a0.value), bigint_to_array(3, 2, a1.value) ];
        }

        let point_array: [Array< [bigint[], bigint[]] >, Array< [bigint[], bigint[]] >] = [pointX, pointY];

        let last_term: Fp = P2y.multiply(P1x).subtract(P2x.multiply(P1y));
        let last_term_Fp12: Fp12 = new Fp12(new Fp6(new Fp2(last_term, Fp.ZERO), Fp2.ZERO, Fp2.ZERO), Fp6.ZERO);
        let out_Fp12: Fp12 = Qx.multiply(P1y.subtract(P2y).value).add(Qy.multiply(P2x.subtract(P1x).value)).add(last_term_Fp12);

        let out_array: Array< [bigint[], bigint[]] > = new Array(6);
        var {c0: y0, c1: y1} = out_Fp12;
        var {c0: c00, c1: c01, c2: c02} = y0;
        var {c0: c10, c1: c11, c2: c12} = y1;
        let out_array_Fp2 : Array< Fp2 > = [c00, c10, c01, c11, c02, c12]; 
        for( var i = 0; i < 6; i++ ){
            var {c0: a0, c1: a1} = out_array_Fp2[i];
            out_array[i] = [ bigint_to_array(3, 2, a0.value), bigint_to_array(3, 2, a1.value) ];
        }

        it('Testing P: ' + in_array + ', Q:' + point_array + /* + in_array[0][0][0] + ', '
	                   + in_array[0][0][1] + ', '
  	                   + in_array[0][1][0] + ', '
  	                   + in_array[0][1][1] + ', '				   
	                   + in_array[1][0][0] + ', '
  	                   + in_array[1][0][1] + ', '
  	                   + in_array[1][1][0] + ', '
			   + in_array[1][1][1] + ', '
          + ' pointX0: '    + point_array[0][0][0][0] + ','
                           + point_array[0][0][0][1] + ','
                           + point_array[0][0][1][0] + ','
                           + point_array[0][0][1][1] + ','			   			   
          + ' pointX1: '    + point_array[0][1][0][0] + ','
                           + point_array[0][1][0][1] + ','
                           + point_array[0][1][1][0] + ','
                           + point_array[0][1][1][1] + ','			   			   
          + ' pointX2: '    + point_array[0][2][0][0] + ','
                           + point_array[0][2][0][1] + ','
                           + point_array[0][2][1][0] + ','
                           + point_array[0][2][1][1] + ','			   			   
          + ' pointX3: '    + point_array[0][3][0][0] + ','
                           + point_array[0][3][0][1] + ','
                           + point_array[0][3][1][0] + ','
                           + point_array[0][3][1][1] + ','			   			   
          + ' pointX4: '    + point_array[0][4][0][0] + ','
                           + point_array[0][4][0][1] + ','
                           + point_array[0][4][1][0] + ','
                           + point_array[0][4][1][1] + ','			   			   
          + ' pointX5: '    + point_array[0][5][0][0] + ','
                           + point_array[0][5][0][1] + ','
                           + point_array[0][5][1][0] + ','
                           + point_array[0][5][1][1] + ','			   			   
          + ' pointY0: '    + point_array[1][0][0][0] + ','
                           + point_array[1][0][0][1] + ','
                           + point_array[1][0][1][0] + ','
                           + point_array[1][0][1][1] + ','			   			   
          + ' pointY1: '    + point_array[1][1][0][0] + ','
                           + point_array[1][1][0][1] + ','
                           + point_array[1][1][1][0] + ','
                           + point_array[1][1][1][1] + ','			   			   
          + ' pointY2: '    + point_array[1][2][0][0] + ','
                           + point_array[1][2][0][1] + ','
                           + point_array[1][2][1][0] + ','
                           + point_array[1][2][1][1] + ','			   			   
          + ' pointY3: '    + point_array[1][3][0][0] + ','
                           + point_array[1][3][0][1] + ','
                           + point_array[1][3][1][0] + ','
                           + point_array[1][3][1][1] + ','			   			   
          + ' pointY4: '    + point_array[1][4][0][0] + ','
                           + point_array[1][4][0][1] + ','
                           + point_array[1][4][1][0] + ','
                           + point_array[1][4][1][1] + ','			   			   
          + ' pointY5: '    + point_array[1][5][0][0] + ','
                           + point_array[1][5][0][1] + ','
                           + point_array[1][5][1][0] + ','
                           + point_array[1][5][1][1] + */ ','
	  + ' p: ' + p, async function() {
            let witness = await circuit.calculateWitness({"P": in_array, "Q": point_array});
	    await circuit.assertOut(witness, {"out": out_array });
            await circuit.checkConstraints(witness);
        });
    }

    var test_cases: Array<[Fp, Fp, Fp, Fp, Fp12, Fp12]> = [];
    for(var test_id = 0; test_id < 100; test_id++){
        let rand_twelveX: BigintTwelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];
	let rand_twelveY: BigintTwelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];	
        for( let i = 0; i < 12; i++){
            rand_twelveX[i] = BigInt( Math.floor(Math.random() * 19) );
            rand_twelveY[i] = BigInt( Math.floor(Math.random() * 19) );	
        }
        let Qx: Fp12 = Fp12.fromBigTwelve( rand_twelveX );
        let Qy: Fp12 = Fp12.fromBigTwelve( rand_twelveY );
	let P1x: Fp = new Fp(BigInt( Math.floor(Math.random() * 19) ));
	let P1y: Fp = new Fp(BigInt( Math.floor(Math.random() * 19) ));
	let P2x: Fp = new Fp(BigInt( Math.floor(Math.random() * 19) ));
	let P2y: Fp = new Fp(BigInt( Math.floor(Math.random() * 19) ));
	test_cases.push([P1x, P1y, P2x, P2y, Qx, Qy]);
    }

    test_cases.forEach(test_linefunc_32);
});

describe("LineFunctionEqual n = 3, k = 2", function() {
    this.timeout(1000 * 1000);

    // runs circom compilation
    let circuit: any;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_linefunc_equal_32.circom"));
    });

    let p: bigint = 19n;
    Fp.ORDER = p;
    Fp2.ORDER = p;

    var test_linefunc_eq_32 = function (inp: [Fp, Fp, Fp12, Fp12]) {
        const [Px, Py, Qx, Qy] = inp;

        let in_array: [ bigint[], bigint[] ] = [ bigint_to_array(3, 2, Px.value), bigint_to_array(3, 2, Py.value) ];

	let pointX: Array< [bigint[], bigint[]] > = new Array(6);
	let pointY: Array< [bigint[], bigint[]] > = new Array(6);	
        var {c0: x0, c1: x1} = Qx;
        var {c0: c00, c1: c01, c2: c02} = x0;
        var {c0: c10, c1: c11, c2: c12} = x1;
        let Fp2_arrayX: Array< Fp2 > = [c00, c10, c01, c11, c02, c12]; 
        for( var i = 0; i < 6; i++ ){
            var {c0: a0, c1: a1} = Fp2_arrayX[i];
            pointX[i] = [ bigint_to_array(3, 2, a0.value), bigint_to_array(3, 2, a1.value) ];
        }

        var {c0: y0, c1: y1} = Qy;
        var {c0: c00, c1: c01, c2: c02} = y0;
        var {c0: c10, c1: c11, c2: c12} = y1;
        let Fp2_arrayY: Array< Fp2 > = [c00, c10, c01, c11, c02, c12]; 
        for( var i = 0; i < 6; i++ ){
            var {c0: a0, c1: a1} = Fp2_arrayY[i];
            pointY[i] = [ bigint_to_array(3, 2, a0.value), bigint_to_array(3, 2, a1.value) ];
        }

        let point_array: [Array< [bigint[], bigint[]] >, Array< [bigint[], bigint[]] >] = [pointX, pointY];

        let last_term: Fp = Px.multiply(Px).multiply(Px).multiply(3n).subtract(Py.multiply(Py).multiply(2n));
	let last_term_Fp12: Fp12 = new Fp12(new Fp6(new Fp2(last_term, Fp.ZERO), Fp2.ZERO, Fp2.ZERO), Fp6.ZERO);
	let out_Fp12: Fp12 = Qy.multiply(Py.multiply(2n).value).subtract(Qx.multiply(Px.multiply(Px).multiply(3n).value)).add(last_term_Fp12);
	let out_array: Array< [bigint[], bigint[]] > = new Array(6);
        var {c0: y0, c1: y1} = out_Fp12;
        var {c0: c00, c1: c01, c2: c02} = y0;
        var {c0: c10, c1: c11, c2: c12} = y1;
	let out_array_Fp2 : Array< Fp2 > = [c00, c10, c01, c11, c02, c12]; 
        for( var i = 0; i < 6; i++ ){
            var {c0: a0, c1: a1} = out_array_Fp2[i];
            out_array[i] = [ bigint_to_array(3, 2, a0.value), bigint_to_array(3, 2, a1.value) ];
        }

        it('Testing P: [' + in_array + '], Q: [' + point_array +  '], out: [' + out_array + '], ' 
        	  + ' p: ' + p, async function() {
	        //await console.log('adsf');
            let witness = await circuit.calculateWitness({"P": in_array, "Q": point_array});
    	    //console.log(witness);			
            //await console.log(in_array, point_array);
    	    //await console.log(out_array); 
	        await circuit.assertOut(witness, {"out": out_array });
            await circuit.checkConstraints(witness);
        });
    }

    var test_cases: Array<[Fp, Fp, Fp12, Fp12]> = [];
    for(var test_id = 0; test_id < 30; test_id++){
        let rand_twelveX: BigintTwelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];
	let rand_twelveY: BigintTwelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];	
        for( let i = 0; i < 12; i++){
            rand_twelveX[i] = BigInt( Math.floor(Math.random() * 19) );
            rand_twelveY[i] = BigInt( Math.floor(Math.random() * 19) );	
        }
        let Qx: Fp12 = Fp12.fromBigTwelve( rand_twelveX );
        let Qy: Fp12 = Fp12.fromBigTwelve( rand_twelveY );
	let Px: Fp = new Fp(BigInt( Math.floor(Math.random() * 19) ));
	let Py: Fp = new Fp(BigInt( Math.floor(Math.random() * 19) ));
	test_cases.push([Px, Py, Qx, Qy]);
    }

    test_cases.forEach(test_linefunc_eq_32);
});

describe("gTimesLineUnequal n = 3, k = 2", function() {
    this.timeout(1000 * 1000);

    // runs circom compilation
    let circuit: any;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_multiply_linefunc_unequal_32.circom"));
    });

    let p: bigint = 19n;
    Fp.ORDER = p;
    Fp2.ORDER = p;

    var test_linefunc_32 = function (inp: [Fp, Fp, Fp, Fp, Fp12, Fp12, Fp12]) {
        const [P1x, P1y, P2x, P2y, Qx, Qy, g] = inp;

        let in_array: Array<[ bigint[], bigint[] ]> = new Array(2);
	    in_array[0] = [ bigint_to_array(3, 2, P1x.value), bigint_to_array(3, 2, P1y.value) ];
        in_array[1] = [ bigint_to_array(3, 2, P2x.value), bigint_to_array(3, 2, P2y.value) ];

        let pointX: Array< [bigint[], bigint[]] > = new Array(6);
        let pointY: Array< [bigint[], bigint[]] > = new Array(6);	
        var {c0: x0, c1: x1} = Qx;
        var {c0: c00, c1: c01, c2: c02} = x0;
        var {c0: c10, c1: c11, c2: c12} = x1;
        let Fp2_arrayX: Array< Fp2 > = [c00, c10, c01, c11, c02, c12]; 
        for( var i = 0; i < 6; i++ ){
            var {c0: a0, c1: a1} = Fp2_arrayX[i];
            pointX[i] = [ bigint_to_array(3, 2, a0.value), bigint_to_array(3, 2, a1.value) ];
        }

        var {c0: y0, c1: y1} = Qy;
        var {c0: c00, c1: c01, c2: c02} = y0;
        var {c0: c10, c1: c11, c2: c12} = y1;
        let Fp2_arrayY: Array< Fp2 > = [c00, c10, c01, c11, c02, c12]; 
        for( var i = 0; i < 6; i++ ){
            var {c0: a0, c1: a1} = Fp2_arrayY[i];
            pointY[i] = [ bigint_to_array(3, 2, a0.value), bigint_to_array(3, 2, a1.value) ];
        }

        let point_array: [Array< [bigint[], bigint[]] >, Array< [bigint[], bigint[]] >] = [pointX, pointY];
        let g_array: Array< [bigint[], bigint[]] > = new Array(6);
        var {c0: g0, c1: g1} = g;
        var {c0: g00, c1: g01, c2: g02} = g0;
        var {c0: g10, c1: g11, c2: g12} = g1;
        let g_array_Fp2: Array< Fp2 > = [g00, g10, g01, g11, g02, g12]; 
        for( var i = 0; i < 6; i++ ){
            var {c0: a0, c1: a1} = g_array_Fp2[i];
            g_array[i] = [ bigint_to_array(3, 2, a0.value), bigint_to_array(3, 2, a1.value) ];
        }


        let last_term: Fp = P2y.multiply(P1x).subtract(P2x.multiply(P1y));
        let last_term_Fp12: Fp12 = new Fp12(new Fp6(new Fp2(last_term, Fp.ZERO), Fp2.ZERO, Fp2.ZERO), Fp6.ZERO);
        let out_Fp12: Fp12 = g.multiply(Qx.multiply(P1y.subtract(P2y).value).add(Qy.multiply(P2x.subtract(P1x).value)).add(last_term_Fp12));
        let out_array = Fp12_to_array(3, 2, out_Fp12);

        it('Testing P: ' + in_array + ', Q:' + point_array + ', g:' + g_array + ', out:' + out_array  
	  + ' p: ' + p, async function() {
            let witness = await circuit.calculateWitness({"P": in_array, "Q": point_array, "g": g_array});
	    await circuit.assertOut(witness, {"out": out_array });
            await circuit.checkConstraints(witness);
        });
    }

    var test_cases: Array<[Fp, Fp, Fp, Fp, Fp12, Fp12, Fp12]> = [];
    for(var test_id = 0; test_id < 100; test_id++){
        let rand_twelveX: BigintTwelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];
        let rand_twelveY: BigintTwelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];	
        let rand_twelveg: BigintTwelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];	
        for( let i = 0; i < 12; i++){
            rand_twelveX[i] = BigInt( Math.floor(Math.random() * 19) );
            rand_twelveY[i] = BigInt( Math.floor(Math.random() * 19) );	
            rand_twelveg[i] = BigInt( Math.floor(Math.random() * 19) );	
        }
        let Qx: Fp12 = Fp12.fromBigTwelve( rand_twelveX );
        let Qy: Fp12 = Fp12.fromBigTwelve( rand_twelveY );
        let g: Fp12 = Fp12.fromBigTwelve( rand_twelveg );
        let P1x: Fp = new Fp(BigInt( Math.floor(Math.random() * 19) ));
        let P1y: Fp = new Fp(BigInt( Math.floor(Math.random() * 19) ));
        let P2x: Fp = new Fp(BigInt( Math.floor(Math.random() * 19) ));
        let P2y: Fp = new Fp(BigInt( Math.floor(Math.random() * 19) ));
        test_cases.push([P1x, P1y, P2x, P2y, Qx, Qy, g]);
    }

    test_cases.forEach(test_linefunc_32);
});

describe("gTimesLineEqual n = 3, k = 2", function() {
    this.timeout(1000 * 1000);

    // runs circom compilation
    let circuit: any;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_multiply_linefunc_equal_32.circom"));
    });

    let p: bigint = 19n;
    Fp.ORDER = p;
    Fp2.ORDER = p;

    var test_linefunc_32 = function (inp: [Fp, Fp, Fp12, Fp12, Fp12]) {
        const [Px, Py, Qx, Qy, g] = inp;

        let in_array: [ bigint[], bigint[] ] = [ bigint_to_array(3, 2, Px.value), bigint_to_array(3, 2, Py.value) ];

        let pointX: Array< [bigint[], bigint[]] > = new Array(6);
        let pointY: Array< [bigint[], bigint[]] > = new Array(6);	
        var {c0: x0, c1: x1} = Qx;
        var {c0: c00, c1: c01, c2: c02} = x0;
        var {c0: c10, c1: c11, c2: c12} = x1;
        let Fp2_arrayX: Array< Fp2 > = [c00, c10, c01, c11, c02, c12]; 
        for( var i = 0; i < 6; i++ ){
            var {c0: a0, c1: a1} = Fp2_arrayX[i];
            pointX[i] = [ bigint_to_array(3, 2, a0.value), bigint_to_array(3, 2, a1.value) ];
        }

        var {c0: y0, c1: y1} = Qy;
        var {c0: c00, c1: c01, c2: c02} = y0;
        var {c0: c10, c1: c11, c2: c12} = y1;
        let Fp2_arrayY: Array< Fp2 > = [c00, c10, c01, c11, c02, c12]; 
        for( var i = 0; i < 6; i++ ){
            var {c0: a0, c1: a1} = Fp2_arrayY[i];
            pointY[i] = [ bigint_to_array(3, 2, a0.value), bigint_to_array(3, 2, a1.value) ];
        }

        let point_array: [Array< [bigint[], bigint[]] >, Array< [bigint[], bigint[]] >] = [pointX, pointY];

        let last_term: Fp = Px.multiply(Px).multiply(Px).multiply(3n).subtract(Py.multiply(Py).multiply(2n));
        let last_term_Fp12: Fp12 = new Fp12(new Fp6(new Fp2(last_term, Fp.ZERO), Fp2.ZERO, Fp2.ZERO), Fp6.ZERO);
        let out_Fp12: Fp12 = g.multiply(Qy.multiply(Py.multiply(2n).value).subtract(Qx.multiply(Px.multiply(Px).multiply(3n).value)).add(last_term_Fp12));

        let g_array: Array< [bigint[], bigint[], bigint[], bigint[]] > = new Array(6);
        var {c0: g0, c1: g1} = g;
        var {c0: g00, c1: g01, c2: g02} = g0;
        var {c0: g10, c1: g11, c2: g12} = g1;
        let g_array_Fp2: Array< Fp2 > = [g00, g10, g01, g11, g02, g12]; 
        for( var i = 0; i < 6; i++ ){
            var {c0: a0, c1: a1} = g_array_Fp2[i];
            g_array[i] = [ bigint_to_array(3, 2, a0.value), bigint_to_array(3, 2, a1.value), bigint_to_array(3, 2, 0n), bigint_to_array(3, 2, 0n) ];
        }

        let out_array = Fp12_to_array(3, 2, out_Fp12);

        it('Testing P: ' + in_array + ', Q:' + point_array + ', g:' + g_array + ', out:' + out_array  
	  + ' p: ' + p, async function() {
            let witness = await circuit.calculateWitness({"P": in_array, "Q": point_array, "g": g_array});
	    await circuit.assertOut(witness, {"out": out_array });
            await circuit.checkConstraints(witness);
        });
    }

    var test_cases: Array<[Fp, Fp, Fp12, Fp12, Fp12]> = [];
    for(var test_id = 0; test_id < 10; test_id++){
        let rand_twelveX: BigintTwelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];
        let rand_twelveY: BigintTwelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];	
        let rand_twelveg: BigintTwelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];	
        for( let i = 0; i < 12; i++){
            rand_twelveX[i] = BigInt( Math.floor(Math.random() * 19) );
            rand_twelveY[i] = BigInt( Math.floor(Math.random() * 19) );	
            rand_twelveg[i] = BigInt( Math.floor(Math.random() * 19) );	
        }
        let Qx: Fp12 = Fp12.fromBigTwelve( rand_twelveX );
        let Qy: Fp12 = Fp12.fromBigTwelve( rand_twelveY );
        let g: Fp12 = Fp12.fromBigTwelve( rand_twelveg );
        let Px: Fp = new Fp(BigInt( Math.floor(Math.random() * 19) ));
        let Py: Fp = new Fp(BigInt( Math.floor(Math.random() * 19) ));
        test_cases.push([Px, Py, Qx, Qy, g]);
    }

    test_cases.forEach(test_linefunc_32);
});

