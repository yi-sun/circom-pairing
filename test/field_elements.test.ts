import path = require("path");
import { expect, assert } from 'chai';
const circom_tester = require('circom_tester');
const wasm_tester = circom_tester.wasm;

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

describe("Fp12Add n = 2, k = 2", function() {
    this.timeout(1000 * 1000);

    // runs circom compilation
    let circuit: any;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_fp12_add_22.circom"));
    });

    // a0, a1, b0, b1, p, c0, c1
    var test_cases: Array<[bigint, bigint, bigint, bigint, bigint, bigint, bigint]> = [];
    let p: bigint = 11n;
    for (var a0 = 0n; a0 < p; a0++) {
        for (var b0 = 0n; b0 < p; b0++) {
	    for (var a1 = 0n; a1 < p; a1++) {
	        for (var b1 = 0n; b1 < p; b1++) {
		    var c0 = (a0 + b0) % p;
		    var c1 = (a1 + b1) % p;
                    test_cases.push([a0, a1, b0, b1, p, c0, c1]);
		}
            }
        }
    }

    var test_field_add_222 = function (x: [bigint, bigint, bigint, bigint, bigint, bigint, bigint]) {
        const [a0, a1, b0, b1, p, c0, c1] = x;

        var a0_array: bigint[] = bigint_to_array(2, 2, a0);
        var a1_array: bigint[] = bigint_to_array(2, 2, a1);	
        var b0_array: bigint[] = bigint_to_array(2, 2, b0);
        var b1_array: bigint[] = bigint_to_array(2, 2, b1);
	var p_array: bigint[] = bigint_to_array(2, 2, p);
        var c0_array: bigint[] = bigint_to_array(2, 2, c0);
        var c1_array: bigint[] = bigint_to_array(2, 2, c1);

        it('Testing a0: ' + a0 + ' a1: ' + a1 + ' b0: ' + b0 + ' b1: ' + b1 + ' p: ' + p + ' c0: ' + c0 + ' c1: ' + c1, async function() {
            let witness = await circuit.calculateWitness({"a": [[a0_array, a1_array],[a0_array, a1_array],[a0_array, a1_array],[a0_array, a1_array],[a0_array, a1_array],[a0_array, a1_array]],
             "b": [[b0_array, b1_array],[b0_array, b1_array],[b0_array, b1_array],[b0_array, b1_array],[b0_array, b1_array],[b0_array, b1_array]], "p": p_array});
	    await circuit.assertOut(witness, {"c": [[c0_array, c1_array], [c0_array, c1_array],[c0_array, c1_array],[c0_array, c1_array],[c0_array, c1_array],[c0_array, c1_array]]});
            await circuit.checkConstraints(witness);
        });
    }

    test_cases.forEach(test_field_add_222);
});

describe("Fp2Add n = 2, k = 2", function() {
    this.timeout(1000 * 1000);

    // runs circom compilation
    let circuit: any;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_fp2_add_22.circom"));
    });

    // a0, a1, b0, b1, p, c0, c1
    var test_cases: Array<[bigint, bigint, bigint, bigint, bigint, bigint, bigint]> = [];
    let p: bigint = 11n;
    for (var a0 = 0n; a0 < p; a0++) {
        for (var b0 = 0n; b0 < p; b0++) {
	    for (var a1 = 0n; a1 < p; a1++) {
	        for (var b1 = 0n; b1 < p; b1++) {
		    var c0 = (a0 + b0) % p;
		    var c1 = (a1 + b1) % p;
                    test_cases.push([a0, a1, b0, b1, p, c0, c1]);
		}
            }
        }
    }

    var test_field_add_222 = function (x: [bigint, bigint, bigint, bigint, bigint, bigint, bigint]) {
        const [a0, a1, b0, b1, p, c0, c1] = x;

        var a0_array: bigint[] = bigint_to_array(2, 2, a0);
        var a1_array: bigint[] = bigint_to_array(2, 2, a1);	
        var b0_array: bigint[] = bigint_to_array(2, 2, b0);
        var b1_array: bigint[] = bigint_to_array(2, 2, b1);
	var p_array: bigint[] = bigint_to_array(2, 2, p);
        var c0_array: bigint[] = bigint_to_array(2, 2, c0);
        var c1_array: bigint[] = bigint_to_array(2, 2, c1);

        it('Testing a0: ' + a0 + ' a1: ' + a1 + ' b0: ' + b0 + ' b1: ' + b1 + ' p: ' + p + ' c0: ' + c0 + ' c1: ' + c1, async function() {
            let witness = await circuit.calculateWitness({"a": [a0_array, a1_array], "b": [b0_array, b1_array], "p": p_array});
	    await circuit.assertOut(witness, {"c": [c0_array, c1_array]});
            await circuit.checkConstraints(witness);
        });
    }

    test_cases.forEach(test_field_add_222);
});