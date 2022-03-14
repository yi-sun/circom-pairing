import path = require("path");

import { expect, assert } from 'chai';
const bls = require('@noble/bls12-381');
import {
  Fp, Fp2, PointG1
} from '@noble/bls12-381';
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

function point_to_bigint(point: PointG1) : [bigint, bigint] {
    let [x , y] = point.toAffine();
    return [x.value, y.value];
}
describe("BLS12-381 AddUnequal", function () {
    this.timeout(1000 * 1000);

    // runs circom compilation
    let circuit: any;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_bls12-381_add.circom"));
    });

    // pub0x, pub0y, pub1x, pub0y, sumx, sumy
    var test_cases: Array<[ [bigint, bigint], [bigint, bigint], [bigint, bigint]]> = [];

    for (var test = 0; test < 10; test++) {
        var pubkeys: Array<PointG1> = [];
        for(var idx = 0; idx < 2; idx++) {
            // randomly generate a point on G1 subgroup of curve
            var pubkey: PointG1 = PointG1.fromPrivateKey( bls.utils.randomPrivateKey() );
            pubkeys.push(pubkey);
        }
        if(pubkeys[0].equals(pubkeys[1]))
            continue;
        var sum: PointG1 = pubkeys[0].add(pubkeys[1]);
        test_cases.push([ point_to_bigint(pubkeys[0]), point_to_bigint(pubkeys[1]), point_to_bigint(sum) ]);
    }

    var test_bls12381_add_instance = function (test_case: [[bigint, bigint], [bigint, bigint], [bigint, bigint]]) {
        let [pub0x, pub0y] = test_case[0];
        let [pub1x, pub1y] = test_case[1];
        let [sumx, sumy] = test_case[2];

        var n: number = 55;
        var k: number = 7;
        var pub0x_array: bigint[] = bigint_to_array(n, k, pub0x);
        var pub0y_array: bigint[] = bigint_to_array(n, k, pub0y);
        var pub1x_array: bigint[] = bigint_to_array(n, k, pub1x);
        var pub1y_array: bigint[] = bigint_to_array(n, k, pub1y);
        var sumx_array: bigint[] = bigint_to_array(n, k, sumx);
        var sumy_array: bigint[] = bigint_to_array(n, k, sumy);

        it('Testing pub0x: ' + pub0x + ' pub0y: ' + pub0y + ' pub1x: ' + pub1x + ' pub1y: ' + pub1y + ' sumx: ' + sumx + ' sumy: ' + sumy, async function() {
            let witness = await circuit.calculateWitness({"a": [pub0x_array, pub0y_array],
                                                          "b": [pub1x_array, pub1y_array]});
	        await circuit.assertOut(witness, {"out": [sumx_array, sumy_array]});
            await circuit.checkConstraints(witness);
        });
    }

    test_cases.forEach(test_bls12381_add_instance);
});

describe("BLS12-381 Double", function () {
    this.timeout(1000 * 1000);

    // runs circom compilation
    let circuit: any;
    before(async function () {
        circuit = await wasm_tester(path.join(__dirname, "circuits", "test_bls12-381_double.circom"));
    });

    // pub0x, pub0y, sumx, sumy
    var test_cases: Array<[ [bigint, bigint], [bigint, bigint]]> = [];

    for (var test = 0; test < 10; test++) {
        var pubkey: PointG1 = PointG1.fromPrivateKey( bls.utils.randomPrivateKey() );
        var sum: PointG1 = pubkey.double();
        test_cases.push([ point_to_bigint(pubkey), point_to_bigint(sum) ]);
    }

    var test_bls12381_double_instance = function (test_case: [[bigint, bigint], [bigint, bigint]]) {
        let [pubx, puby] = test_case[0];
        let [sumx, sumy] = test_case[1];

        var n: number = 55;
        var k: number = 7;
        var pubx_array: bigint[] = bigint_to_array(n, k, pubx);
        var puby_array: bigint[] = bigint_to_array(n, k, puby);
        var sumx_array: bigint[] = bigint_to_array(n, k, sumx);
        var sumy_array: bigint[] = bigint_to_array(n, k, sumy);

        it('Testing x: ' + pubx + ' y: ' + puby + ' doublex: ' + sumx + ' doubley: ' + sumy, async function() {
            let witness = await circuit.calculateWitness({"in": [pubx_array, puby_array]});
	        await circuit.assertOut(witness, {"out": [sumx_array, sumy_array]});
            await circuit.checkConstraints(witness);
        });
    }

    test_cases.forEach(test_bls12381_double_instance);
});
