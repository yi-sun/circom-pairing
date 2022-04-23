//const bls = require("@noble/bls12-381");
import { ssz } from "@chainsafe/lodestar-types";
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
//import block_json from "./BeaconBlock3644983.json";

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

function formatHex(str: string): string {
  if (str.startsWith("0x")) {
    str = str.slice(2);
  }
  return str;
}

// UTF8 to ui8a
// first 128 char agree with ASCII encoding
function stringToBytes(str: string) {
  const bytes = new Uint8Array(str.length);
  for (let i = 0; i < str.length; i++) {
    bytes[i] = str.charCodeAt(i);
  }
  return bytes;
}

function hexToBytes(hex: string, endian: string = "big"): Uint8Array {
  if (typeof hex !== "string") {
    throw new TypeError("hexToBytes: expected string, got " + typeof hex);
  }
  hex = formatHex(hex);
  if (hex.length % 2)
    throw new Error("hexToBytes: received invalid unpadded hex");
  const array = new Uint8Array(hex.length / 2);
  for (let i = 0; i < array.length; i++) {
    let j = 0;
    if (endian === "big") j = i * 2;
    else j = (array.length - 1 - i) * 2;

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

// generate test input for bls signature verify
async function test(message: string) {
  let msg = stringToBytes(message);
  //console.log(msg);
  let u = await hashToField(msg, 2);

  console.log(
    "u[0] = 0x" + u[0][0].toString(16) + "\n + I * 0x" + u[0][1].toString(16)
  );
  console.log(
    "u[1] = 0x" + u[1][0].toString(16) + "\n + I * 0x" + u[1][1].toString(16)
  );

  let u_array = [
    Fp2_to_array(55, 7, Fp2.fromBigTuple(u[0])),
    Fp2_to_array(55, 7, Fp2.fromBigTuple(u[1])),
  ];
  console.log("u_array : ");
  console.log(JSON.stringify(u_array));

  /*
  const Q0p = new PointG2(...map_to_curve_simple_swu_9mod16(u[0]));
  const Q1p = new PointG2(...map_to_curve_simple_swu_9mod16(u[1]));
  const Q0 = new PointG2(...isogenyMapG2(map_to_curve_simple_swu_9mod16(u[0])));
  const Q1 = new PointG2(...isogenyMapG2(map_to_curve_simple_swu_9mod16(u[1])));
  const R = Q0.add(Q1);
  const P = R.clearCofactor();
  console.log(P.toAffine()[1].c0.value.toString(16));

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
  console.log(Hm.toAffine()[0].c0.value.toString(16));

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
  console.log("x = 0x" + publicKey.toAffine()[0].value.toString(16));
  console.log("y = 0x" + publicKey.toAffine()[1].value.toString(16));
  console.log(
    JSON.stringify([
      bigint_to_array(55, 7, publicKey.toAffine()[0].value),
      bigint_to_array(55, 7, publicKey.toAffine()[1].value),
    ])
  );

  console.log("signature:");
  console.log(
    "x = 0x" +
      signature.toAffine()[0].c0.value.toString(16) +
      "\n + I * 0x" +
      signature.toAffine()[0].c1.value.toString(16)
  );
  console.log(
    "y = 0x" +
      signature.toAffine()[1].c0.value.toString(16) +
      "\n + I * 0x" +
      signature.toAffine()[1].c1.value.toString(16)
  );
  console.log(
    JSON.stringify([
      Fp2_to_array(55, 7, signature.toAffine()[0]),
      Fp2_to_array(55, 7, signature.toAffine()[1]),
    ])
  );
}

//test("devconnect");

export const DOMAIN_BEACON_PROPOSER = Uint8Array.from([0, 0, 0, 0]);
/*export const BeaconBlock = new ContainerType(
  {
    slot: new UintNumberType(8),
    proposer_index: new UintNumberType(8),
    parent_root: new ByteVectorType(32),
    state_root: new ByteVectorType(32),
    body_root: new ByteVectorType(32),
  },
  { typeName: "BeaconBlock", jsonCase: "eth2" }
);*/

// takes beacon chain slot data and retrieves relevant inputs for signature verification
async function verify_block_signature() {
  // example beacon chain block: https://beaconcha.in/block/3644983
  let publicKeyHex: string =
    "0x932b42ad9a01e2c489958bb212af2dc016f02dd2750980f618420b6f8fccb469de8bc63c0b594f06464a3f09169a8825";
  publicKeyHex = formatHex(publicKeyHex);
  const publicKey: PointG1 = PointG1.fromHex(formatHex(publicKeyHex));
  publicKey.assertValidity();

  let signatureHex: string =
    "0x8530f2e4403406b78ddfd3a94bf2085ce325e17c7eadf57d01311f11518c11621b764c8281618197077568dbb8ae7cea19bac59893c09c107581d9dc88aa461fe1e631f2b2ee3b3eec0b12ee97d6437ac2fca5d3e40474b87d72a301fe59974b";

  signatureHex = formatHex(signatureHex);
  const signature: PointG2 = PointG2.fromSignature(signatureHex);
  signature.assertValidity();

  /* Downloaded entire block to confirm that "Block Root" really refers to the merkleization of entire block, not just the "Block Body", so this part is no longer neccessary 
  const BeaconBlock = ssz.altair.BeaconBlock;
  let block = BeaconBlock.fromJson(block_json);
  let block_root = BeaconBlock.hashTreeRoot(block);
  console.log(bytesToHex(block_root));
  */
  let block_root = hexToBytes(
    "0xf91113e2837f5197ea2bd174f8c79ef31182940a6059e8a1a2d352a958281f27"
  );
  // see compute_domain and get_domain from beacon chain spec
  // NOTE: ethereum spec for ssz uses little endian for ByteVector while @chainsafe/ssz uses big Endian

  // infura API call gives:
  let genesisValidatorsRoot: Uint8Array = hexToBytes(
    "0x4b363db94e286120d76eb905340fdd4e54bfe9f06bf33ff6cf5ad27f511bfe95"
  );

  const ForkData = ssz.phase0.ForkData;
  let fork_data = ForkData.defaultValue();
  fork_data.currentVersion = hexToBytes("0x01000000"); // altair
  fork_data.genesisValidatorsRoot = genesisValidatorsRoot;
  let fork_data_root = ForkData.hashTreeRoot(fork_data);
  let domain = new Uint8Array(32);
  for (let i = 0; i < 4; i++) domain[i] = DOMAIN_BEACON_PROPOSER[i];
  for (let i = 0; i < 28; i++) domain[i + 4] = fork_data_root[i];
  //console.log(bytesToHex(domain));

  // see compute_signing_root from beacon chain spec
  const SigningData = ssz.phase0.SigningData;
  let signing_data = SigningData.defaultValue();
  signing_data.objectRoot = block_root;
  signing_data.domain = domain;
  const signing_root = SigningData.hashTreeRoot(signing_data);

  //console.log(msg);
  //const Hm = await PointG2.hashToCurve(signing_root);
  const isCorrect = await verify(signature, signing_root, publicKey);
  console.log(isCorrect);

  // print input formats for web demo and circuit itself
  console.log("publicKey:");
  console.log("x = 0x" + publicKey.toAffine()[0].value.toString(16));
  console.log("y = 0x" + publicKey.toAffine()[1].value.toString(16));
  /*
  console.log(
    JSON.stringify([
      bigint_to_array(55, 7, publicKey.toAffine()[0].value),
      bigint_to_array(55, 7, publicKey.toAffine()[1].value),
    ])
  );
  */

  console.log("signature:");
  console.log(
    "x = 0x" +
      signature.toAffine()[0].c0.value.toString(16) +
      "\n + I * 0x" +
      signature.toAffine()[0].c1.value.toString(16)
  );
  console.log(
    "y = 0x" +
      signature.toAffine()[1].c0.value.toString(16) +
      "\n + I * 0x" +
      signature.toAffine()[1].c1.value.toString(16)
  );
  /*
  console.log(
    JSON.stringify([
      Fp2_to_array(55, 7, signature.toAffine()[0]),
      Fp2_to_array(55, 7, signature.toAffine()[1]),
    ])
  );
  */

  let u = await hashToField(signing_root, 2);
  console.log("message:");
  console.log(
    "x = 0x" + u[0][0].toString(16) + "\n + I * 0x" + u[0][1].toString(16)
  );
  console.log(
    "y = 0x" + u[1][0].toString(16) + "\n + I * 0x" + u[1][1].toString(16)
  );
  /*
  let u_array = [
    Fp2_to_array(55, 7, Fp2.fromBigTuple(u[0])),
    Fp2_to_array(55, 7, Fp2.fromBigTuple(u[1])),
  ];
  console.log("array : ");
  console.log(JSON.stringify(u_array));
  */
}

verify_block_signature();
