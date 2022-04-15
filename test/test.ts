const bls = require("@noble/bls12-381");
import {
  Fp,
  Fp2,
  Fp12,
  CURVE,
  PointG1,
  utils,
  PointG2,
} from "@noble/bls12-381";
const hashToField = utils.hashToField;

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

  let ret: string[] = [];
  var x_temp: bigint = x;
  for (var idx = 0; idx < k; idx++) {
    ret.push((x_temp % mod).toString());
    x_temp = x_temp / mod;
  }
  return ret;
}

let p: bigint = Fp.ORDER;

function printFp2(x: Fp2) {
  let [c0, c1] = x.values;
  return [c0, c1];
}

function Fp2_to_array(n: number, k: number, x: Fp2) {
  let [c0, c1] = x.values;
  return [bigint_to_array(n, k, c0), bigint_to_array(n, k, c1)];
}

function printFp12(x: Fp12) {
  let [c0, c1] = x.c;
  let [c00, c01, c02] = c0.c;
  let [c10, c11, c12] = c1.c;
  return [
    printFp2(c00),
    printFp2(c10),
    printFp2(c01),
    printFp2(c11),
    printFp2(c02),
    printFp2(c12),
  ];
}

// UTF8 to ui8a
function stringToBytes(str: string) {
  const bytes = new Uint8Array(str.length);
  for (let i = 0; i < str.length; i++) {
    bytes[i] = str.charCodeAt(i);
  }
  return bytes;
}

function hexToBytes(hex: string): Uint8Array {
  if (typeof hex !== "string") {
    throw new TypeError("hexToBytes: expected string, got " + typeof hex);
  }
  if (hex.length % 2)
    throw new Error("hexToBytes: received invalid unpadded hex");
  const array = new Uint8Array(hex.length / 2);
  for (let i = 0; i < array.length; i++) {
    const j = i * 2;
    const hexByte = hex.slice(j, j + 2);
    if (hexByte.length !== 2) throw new Error("Invalid byte sequence");
    const byte = Number.parseInt(hexByte, 16);
    if (Number.isNaN(byte) || byte < 0)
      throw new Error("Invalid byte sequence");
    array[i] = byte;
  }
  return array;
}

function ensureBytes(hex: string | Uint8Array): Uint8Array {
  // Uint8Array.from() instead of hash.slice() because node.js Buffer
  // is instance of Uint8Array, and its slice() creates **mutable** copy
  return hex instanceof Uint8Array ? Uint8Array.from(hex) : hexToBytes(hex);
}

async function test(message: string) {
  let msg = stringToBytes(message);
  console.log(msg);

  let u = await hashToField(msg, 2);
  let u_array = [
    Fp2_to_array(55, 7, new Fp2(u[0])),
    Fp2_to_array(55, 7, new Fp2(u[1])),
  ];
  //console.log("u : ");
  //console.log(u);
  console.log("u_array : ");
  console.log(JSON.stringify(u_array));

  let P = await PointG2.hashToCurve(msg);
  let x = P.x.multiply(P.z.invert());
  let y = P.y.multiply(P.z.invert());
  let out_array = [Fp2_to_array(55, 7, x), Fp2_to_array(55, 7, y)];
  console.log("MapToG2 out:");
  console.log(JSON.stringify(out_array));
}

test("abc");
