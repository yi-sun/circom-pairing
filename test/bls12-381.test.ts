import path = require("path");

import { expect, assert } from "chai";
const bls = require("@noble/bls12-381");
import { Fp, Fp2, PointG1 } from "@noble/bls12-381";
const circom_tester = require("circom_tester");
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

function point_to_bigint(point: PointG1): [bigint, bigint] {
  let [x, y] = point.toAffine();
  return [x.value, y.value];
}
describe("BLS12-381 AddUnequal", function () {
  this.timeout(1000 * 1000);

  // runs circom compilation
  let circuit: any;
  before(async function () {
    circuit = await wasm_tester(
      path.join(__dirname, "circuits", "test_bls12-381_add.circom")
    );
  });

  // pub0x, pub0y, pub1x, pub0y, sumx, sumy
  var test_cases: Array<
    [[bigint, bigint], [bigint, bigint], [bigint, bigint]]
  > = [];

  for (var test = 0; test < 10; test++) {
    var pubkeys: Array<PointG1> = [];
    for (var idx = 0; idx < 2; idx++) {
      // randomly generate a point on G1 subgroup of curve
      var pubkey: PointG1 = PointG1.fromPrivateKey(
        bls.utils.randomPrivateKey()
      );
      pubkeys.push(pubkey);
    }
    if (pubkeys[0].equals(pubkeys[1])) continue;
    var sum: PointG1 = pubkeys[0].add(pubkeys[1]);
    test_cases.push([
      point_to_bigint(pubkeys[0]),
      point_to_bigint(pubkeys[1]),
      point_to_bigint(sum),
    ]);
  }

  var test_bls12381_add_instance = function (
    test_case: [[bigint, bigint], [bigint, bigint], [bigint, bigint]]
  ) {
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

    it(
      "Testing pub0x: " +
        pub0x +
        " pub0y: " +
        pub0y +
        " pub1x: " +
        pub1x +
        " pub1y: " +
        pub1y +
        " sumx: " +
        sumx +
        " sumy: " +
        sumy,
      async function () {
        let witness = await circuit.calculateWitness({
          a: [pub0x_array, pub0y_array],
          b: [pub1x_array, pub1y_array],
        });
        await circuit.assertOut(witness, { out: [sumx_array, sumy_array] });
        await circuit.checkConstraints(witness);
      }
    );
  };

  test_cases.forEach(test_bls12381_add_instance);
});

describe("BLS12-381 Double", function () {
  this.timeout(1000 * 1000);

  // runs circom compilation
  let circuit: any;
  before(async function () {
    circuit = await wasm_tester(
      path.join(__dirname, "circuits", "test_bls12-381_double.circom")
    );
  });

  // pub0x, pub0y, sumx, sumy
  var test_cases: Array<[[bigint, bigint], [bigint, bigint]]> = [];

  for (var test = 0; test < 10; test++) {
    var pubkey: PointG1 = PointG1.fromPrivateKey(bls.utils.randomPrivateKey());
    var sum: PointG1 = pubkey.double();
    test_cases.push([point_to_bigint(pubkey), point_to_bigint(sum)]);
  }

  var test_bls12381_double_instance = function (
    test_case: [[bigint, bigint], [bigint, bigint]]
  ) {
    let [pubx, puby] = test_case[0];
    let [sumx, sumy] = test_case[1];

    var n: number = 55;
    var k: number = 7;
    var pubx_array: bigint[] = bigint_to_array(n, k, pubx);
    var puby_array: bigint[] = bigint_to_array(n, k, puby);
    var sumx_array: bigint[] = bigint_to_array(n, k, sumx);
    var sumy_array: bigint[] = bigint_to_array(n, k, sumy);

    it(
      "Testing x: " +
        pubx +
        " y: " +
        puby +
        " doublex: " +
        sumx +
        " doubley: " +
        sumy,
      async function () {
        let witness = await circuit.calculateWitness({
          in: [pubx_array, puby_array],
        });
        await circuit.assertOut(witness, { out: [sumx_array, sumy_array] });
        await circuit.checkConstraints(witness);
      }
    );
  };

  test_cases.forEach(test_bls12381_double_instance);
});

describe("BLS12-381 AddThree Inverses", function () {
  this.timeout(1000 * 1000);

  // runs circom compilation
  let circuit: any;
  before(async function () {
    circuit = await wasm_tester(
      path.join(__dirname, "circuits", "test_bls12-381_add_three.circom")
    );
  });

  // [[ax, ay], [bx, by], [cx, cy]], [[sumx, sumy], isInf]
  var test_cases: Array<[[[bigint, bigint], [bigint, bigint], [bigint, bigint]], [[bigint, bigint], bigint]]> = [];

  let p : PointG1 = PointG1.fromPrivateKey(bls.utils.randomPrivateKey());
  let pubkeys : [PointG1, PointG1] = [p, p.negate()];
  for (var i = 0; i < 2; i++ ) for (var j = 0; j < 2; j++) for (var k = 0; k < 2; k++) {
    let sum : PointG1 = pubkeys[i].add(pubkeys[j]).add(pubkeys[k]);
    test_cases.push([[point_to_bigint(pubkeys[i]), point_to_bigint(pubkeys[j]), point_to_bigint(pubkeys[k])], sum.isZero() ? [point_to_bigint(pubkeys[i].add(pubkeys[j])), 1n ] : [point_to_bigint(sum), 0n]]);
  }
  

  var test_bls12381_addthree_instance = function (
    test_case: [[[bigint, bigint], [bigint, bigint], [bigint, bigint]], [[bigint, bigint], bigint]]
  ) {
    let [[ax, ay], [bx, by], [cx, cy]] = test_case[0];
    let [[sumx, sumy], isInf] = test_case[1];

    var n: number = 55;
    var k: number = 7;
    var ax_array: bigint[] = bigint_to_array(n, k, ax);
    var ay_array: bigint[] = bigint_to_array(n, k, ay);
    var bx_array: bigint[] = bigint_to_array(n, k, bx);
    var by_array: bigint[] = bigint_to_array(n, k, by);
    var cx_array: bigint[] = bigint_to_array(n, k, cx);
    var cy_array: bigint[] = bigint_to_array(n, k, cy);

    var sumx_array: bigint[] = bigint_to_array(n, k, sumx);
    var sumy_array: bigint[] = bigint_to_array(n, k, sumy);

    it(
      "Testing ax: " +
        ax +
        " ay: " +
        ay +
        " bx: " +
        bx +
        " by: " +
        by +
        " cx: " +
        cx +
        " cy: " +
        cy +
        " sumx: " +
        sumx +
        " sumy: " +
        sumy +
        " isInf: " +
        isInf,
      async function () {
        let witness = await circuit.calculateWitness({
          a: [ax_array, ay_array],
          b: [bx_array, by_array],
          c: [cx_array, cy_array],
        });
        await circuit.assertOut(witness, { out: [sumx_array, sumy_array], isInfinity: isInf });
        await circuit.checkConstraints(witness);
      }
    );
  };

  test_cases.forEach(test_bls12381_addthree_instance);
});

describe("BLS12-381 EllipticCurveAdd Special Cases", function () {
  this.timeout(1000 * 1000);

  // runs circom compilation
  let circuit: any;
  before(async function () {
    circuit = await wasm_tester(
      path.join(__dirname, "circuits", "test_bls12-381_add_two.circom")
    );
  });

  // [[[ax, ay], aIsInf], [[bx, by], bIsInf]], [sumx, sumy], isInf]
  var test_cases: Array<[[[[bigint, bigint], bigint], [[bigint, bigint], bigint]], [[bigint, bigint], bigint]]> = [];

  let p : PointG1 = PointG1.fromPrivateKey(bls.utils.randomPrivateKey());
  let pubkeys : [PointG1, PointG1] = [p, p.negate()];

  function sum_points(a: PointG1, aIsInf : boolean, b: PointG1, bIsInf: boolean) : [PointG1, boolean] {
    if (bIsInf) {
      return [a, aIsInf];
    }
    if (aIsInf) {
      return [b, bIsInf];
    }
    let sum : PointG1 = a.add(b);
    if (sum.isZero()) {
      return [a, true];
    }
    return [sum, false];
  }

  for (var i = 0n; i < 2n; i++ ) for (var j = 0; j < 2; j++) for (var k = 0n; k < 2n; k++) {
    let [sum, isInf] = sum_points(pubkeys[0], !!i,  pubkeys[j], !!k)
    test_cases.push([[[point_to_bigint(pubkeys[0]), BigInt(i)], 
                      [point_to_bigint(pubkeys[j]), BigInt(k)]],
                      [point_to_bigint(sum), BigInt(isInf)]
                    ]);
  }
  

  var test_bls12381_addtwo_instance = function (
    test_case: [[[[bigint, bigint], bigint], [[bigint, bigint], bigint]], [[bigint, bigint], bigint]]
  ) {
    let [[[ax, ay], aIsInf], [[bx, by], bIsInf]] = test_case[0];
    let [[sumx, sumy], isInf] = test_case[1];

    var n: number = 55;
    var k: number = 7;
    var ax_array: bigint[] = bigint_to_array(n, k, ax);
    var ay_array: bigint[] = bigint_to_array(n, k, ay);
    var bx_array: bigint[] = bigint_to_array(n, k, bx);
    var by_array: bigint[] = bigint_to_array(n, k, by);

    var sumx_array: bigint[] = bigint_to_array(n, k, sumx);
    var sumy_array: bigint[] = bigint_to_array(n, k, sumy);

    it(
      "Testing ax: " +
        ax +
        " ay: " +
        ay +
        "aIsInf: " +
        aIsInf +
        " bx: " +
        bx +
        " by: " +
        by +
        " bIsInf: " +
        bIsInf +
        " sumx: " +
        sumx +
        " sumy: " +
        sumy +
        " sumIsInf: " +
        isInf,
      async function () {
        let witness = await circuit.calculateWitness({
          a: [ax_array, ay_array],
          aIsInfinity: aIsInf,
          b: [bx_array, by_array],
          bIsInfinity: bIsInf
        });
        await circuit.assertOut(witness, { out: [sumx_array, sumy_array], isInfinity: isInf });
        await circuit.checkConstraints(witness);
      }
    );
  };

  test_cases.forEach(test_bls12381_addtwo_instance);
});

