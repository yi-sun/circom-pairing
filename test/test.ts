//const bls = require("@noble/bls12-381");
import {
  Fp,
  Fp2,
  Fp12,
  CURVE,
  PointG1,
  utils,
  PointG2,
  getPublicKey,
  sign,
  verify,
} from "./index";
import { map_to_curve_simple_swu_9mod16, isogenyMapG2 } from "./math";
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
  let { c0, c1 } = x;
  return [c0.value, c1.value];
}

function Fp2_to_array(n: number, k: number, x: Fp2) {
  let { c0, c1 } = x;
  return [bigint_to_array(n, k, c0.value), bigint_to_array(n, k, c1.value)];
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

function bytesToNumberBE(uint8a: Uint8Array): bigint {
  if (!(uint8a instanceof Uint8Array)) throw new Error("Expected Uint8Array");
  return BigInt("0x" + bytesToHex(Uint8Array.from(uint8a)));
}

const hexes = Array.from({ length: 256 }, (v, i) =>
  i.toString(16).padStart(2, "0")
);
function bytesToHex(uint8a: Uint8Array): string {
  // pre-caching chars could speed this up 6x.
  let hex = "";
  for (let i = 0; i < uint8a.length; i++) {
    hex += hexes[uint8a[i]];
  }
  return hex;
}

function sgn0(x: Fp2) {
  const { re: x0, im: x1 } = x.reim();
  const sign_0 = x0 % 2n;
  const zero_0 = x0 === 0n;
  const sign_1 = x1 % 2n;
  return BigInt(sign_0 || (zero_0 && sign_1));
}

async function test(message: string) {
  let msg = stringToBytes(message);
  //console.log(msg);

  let u = await hashToField(msg, 2);
  let u_array = [
    Fp2_to_array(55, 7, Fp2.fromBigTuple(u[0])),
    Fp2_to_array(55, 7, Fp2.fromBigTuple(u[1])),
  ];
  //console.log("u : ");
  //console.log(u);
  console.log("u_array : ");
  console.log(JSON.stringify(u_array));

  /*
  const Q0p = new PointG2(...map_to_curve_simple_swu_9mod16(u[0]));
  const Q1p = new PointG2(...map_to_curve_simple_swu_9mod16(u[1]));
  console.log(sgn0(Fp2.fromBigTuple(u[1])));
  console.log(sgn0(Q1p.toAffine()[1]));
  const Q0 = new PointG2(...isogenyMapG2(map_to_curve_simple_swu_9mod16(u[0])));
  const Q1 = new PointG2(...isogenyMapG2(map_to_curve_simple_swu_9mod16(u[1])));
  const R = Q0.add(Q1);
  
  console.log("Q0:");
  console.log(Q0.toAffine());
  console.log(
    JSON.stringify([
      Fp2_to_array(55, 7, Q0.toAffine()[0]),
      Fp2_to_array(55, 7, Q0.toAffine()[1]),
    ])
  );
  
  console.log("Q0 + Q1");
  //console.log(R.toAffine());
  console.log(
    JSON.stringify([
      Fp2_to_array(55, 7, R.toAffine()[0]),
      Fp2_to_array(55, 7, R.toAffine()[1]),
    ])
  );*/

  let Hm = await PointG2.hashToCurve(msg);
  console.log("MapToG2 out:");
  console.log(
    JSON.stringify([
      Fp2_to_array(55, 7, Hm.toAffine()[0]),
      Fp2_to_array(55, 7, Hm.toAffine()[1]),
    ])
  );

  const privateKey =
    "67d53f170b908cabb9eb326c3c337762d59289a8fec79f7bc9254b584b73265c";
  const publicKeyHex = getPublicKey(privateKey);
  const signature = await sign(Hm, privateKey);
  const isCorrect = await verify(signature, Hm, publicKeyHex);
  console.log("valid signature? " + isCorrect);
  let publicKey = PointG1.fromHex(publicKeyHex);

  console.log("publicKey:");
  console.log(
    JSON.stringify([
      bigint_to_array(55, 7, publicKey.toAffine()[0].value),
      bigint_to_array(55, 7, publicKey.toAffine()[1].value),
    ])
  );

  console.log("signature:");
  console.log(
    JSON.stringify([
      Fp2_to_array(55, 7, signature.toAffine()[0]),
      Fp2_to_array(55, 7, signature.toAffine()[1]),
    ])
  );
}

test("abc");
