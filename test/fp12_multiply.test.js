"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
const path = require("path");
const circom_tester = require('circom_tester');
const wasm_tester = circom_tester.wasm;
const math_1 = require("./math");
function bigint_to_array(n, k, x) {
    let mod = 1n;
    for (var idx = 0; idx < n; idx++) {
        mod = mod * 2n;
    }
    let ret = [];
    var x_temp = x;
    for (var idx = 0; idx < k; idx++) {
        ret.push(x_temp % mod);
        x_temp = x_temp / mod;
    }
    return ret;
}
describe("Fp12Multiply n = 3, k = 2", function () {
    this.timeout(1000 * 1000);
    // runs circom compilation
    let circuit;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_fp12_multiply_32.circom"));
    });
    // a0, a1, b0, b1, p, c0, c1
    var test_cases = [];
    let p = 13n;
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
    var test_field_multiply_32 = function (x) {
        const [a0, a1, b0, b1, p, c0, c1] = x;
        var a0_array = bigint_to_array(3, 2, a0);
        var a1_array = bigint_to_array(3, 2, a1);
        var b0_array = bigint_to_array(3, 2, b0);
        var b1_array = bigint_to_array(3, 2, b1);
        var p_array = bigint_to_array(3, 2, p);
        var c0_array = bigint_to_array(3, 2, c0);
        var c1_array = bigint_to_array(3, 2, c1);
        var c0i_array = bigint_to_array(3, 2, (c0 - c1 + p) % p);
        var c1i_array = bigint_to_array(3, 2, (c0 + c1) % p);
        var zero = bigint_to_array(3, 2, 0n);
        it('Testing a0: ' + a0 + ' a1: ' + a1 + ' b0: ' + b0 + ' b1: ' + b1 + ' p: ' + p + ' c0: ' + c0 + ' c1: ' + c1, async function () {
            let witness = await circuit.calculateWitness({ "a": [[zero, zero], [zero, zero], [zero, zero], [a0_array, a1_array], [zero, zero], [a0_array, a1_array]],
                "b": [[b0_array, b1_array], [zero, zero], [zero, zero], [b0_array, b1_array], [zero, zero], [zero, zero]], "p": p_array });
            await circuit.assertOut(witness, { "out": [[c0i_array, c1i_array], [zero, zero], [c0i_array, c1i_array],
                    [c0_array, c1_array], [zero, zero], [c0_array, c1_array]] });
            await circuit.checkConstraints(witness);
        });
    };
    test_cases.forEach(test_field_multiply_32);
});
describe("Fp12Multiply2 n = 3, k = 2", function () {
    this.timeout(1000 * 1000);
    // runs circom compilation
    let circuit;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_fp12_multiply2_32.circom"));
    });
    // a0, a1, b0, b1, c0, c1
    let p = 13n;
    var test_cases = [];
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
    var test_field_multiply_32 = function (x) {
        const [a0, a1, b0, b1, c0, c1] = x;
        var a0_array = bigint_to_array(3, 2, a0);
        var a1_array = bigint_to_array(3, 2, a1);
        var b0_array = bigint_to_array(3, 2, b0);
        var b1_array = bigint_to_array(3, 2, b1);
        var c0_array = bigint_to_array(3, 2, c0);
        var c1_array = bigint_to_array(3, 2, c1);
        var c0i_array = bigint_to_array(3, 2, (c0 - c1 + p) % p);
        var c1i_array = bigint_to_array(3, 2, (c0 + c1) % p);
        var zero = bigint_to_array(3, 2, 0n);
        it('Testing a0: ' + a0 + ' a1: ' + a1 + ' b0: ' + b0 + ' b1: ' + b1 + ' p: ' + p + ' c0: ' + c0 + ' c1: ' + c1, async function () {
            let witness = await circuit.calculateWitness({ "a": [[zero, zero], [zero, zero], [zero, zero], [a0_array, a1_array], [zero, zero], [a0_array, a1_array]],
                "b": [[b0_array, b1_array], [zero, zero], [zero, zero], [b0_array, b1_array], [zero, zero], [zero, zero]] });
            await circuit.assertOut(witness, { "out": [[c0i_array, c1i_array], [zero, zero], [c0i_array, c1i_array],
                    [c0_array, c1_array], [zero, zero], [c0_array, c1_array]] });
            await circuit.checkConstraints(witness);
        });
    };
    test_cases.forEach(test_field_multiply_32);
});
describe("Fp12Compression n = 3, k = 2", function () {
    this.timeout(1000 * 1000);
    // runs circom compilation
    let circuit;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_fp12_compression_32.circom"));
    });
    let p = 13n;
    math_1.Fp.ORDER = p;
    math_1.Fp2.ORDER = p;
    var test_compression_32 = function (x) {
        let in_array = new Array(6);
        let a_array = new Array(6);
        let { c0: x0, c1: x1 } = x;
        let { c0: c00, c1: c01, c2: c02 } = x0;
        let { c0: c10, c1: c11, c2: c12 } = x1;
        let Fp2_array = [c00, c10, c01, c11, c02, c12];
        for (var i = 0; i < 6; i++) {
            let { c0: a0, c1: a1 } = Fp2_array[i];
            in_array[i] = [bigint_to_array(3, 2, a0.value), bigint_to_array(3, 2, a1.value)];
            a_array[i] = [a0.value, a1.value];
        }
        it('Testing a0: ' + a_array[0][0] + ', ' + a_array[0][1] +
            ' a1: ' + a_array[1][0] + ', ' + a_array[1][1] +
            ' a2: ' + a_array[2][0] + ', ' + a_array[2][1] +
            ' a3: ' + a_array[3][0] + ', ' + a_array[3][1] +
            ' a4: ' + a_array[4][0] + ', ' + a_array[4][1] +
            ' a5: ' + a_array[5][0] + ', ' + a_array[5][1] +
            ' p: ' + p, async function () {
            let witness = await circuit.calculateWitness({ "in": in_array });
            await circuit.assertOut(witness, { "out": in_array });
            await circuit.checkConstraints(witness);
        });
    };
    var test_cases = [];
    for (var test_id = 0; test_id < 1000; test_id++) {
        let rand_twelve = [0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n];
        for (let i = 0; i < 12; i++) {
            rand_twelve[i] = BigInt(Math.floor(Math.random() * 13));
        }
        let elt = math_1.Fp12.fromBigTwelve(rand_twelve);
        let elt2 = elt.pow((p ** 6n - 1n) * (p * p + 1n)); // this is now in cyclotomic subgroup 
        if (!elt2.equals(math_1.Fp12.ONE) && elt2.pow(p ** 4n - p ** 2n + 1n).equals(math_1.Fp12.ONE)) // check it's in cyclotomic subgroup
            test_cases.push(elt2);
    }
    test_cases.forEach(test_compression_32);
});
