/*! noble-bls12-381 - MIT License (c) 2019 Paul Miller (paulmillr.com) */
// bls12-381 is a construction of two curves:
// 1. Fp: (x, y)
// 2. Fp₂: ((x₁, x₂+i), (y₁, y₂+i)) - (complex numbers)
//
// Bilinear Pairing (ate pairing) is used to combine both elements into a paired one:
//   Fp₁₂ = e(Fp, Fp2)
//   where Fp₁₂ = 12-degree polynomial
// Pairing is used to verify signatures.
//
// We are using Fp for private keys (shorter) and Fp2 for signatures (longer).
// Some projects may prefer to swap this relation, it is not supported for now.

import nodeCrypto from "crypto";
// prettier-ignore
import {
  Fp, Fr, Fp2, Fp12, CURVE, ProjectivePoint,
  map_to_curve_simple_swu_9mod16, isogenyMapG2,
  millerLoop, psi, psi2, calcPairingPrecomputes, mod
} from './math.js';
export { Fp, Fr, Fp2, Fp12, CURVE };

type Hex = Uint8Array | string;
type PrivateKey = Hex | bigint | number;
const POW_2_381 = 2n ** 381n;
const POW_2_382 = POW_2_381 * 2n;
const POW_2_383 = POW_2_382 * 2n;
const PUBLIC_KEY_LENGTH = 48;
const SHA256_DIGEST_SIZE = 32;

// Default hash_to_field options are for hash to G2.
//
// Parameter definitions are in section 5.3 of the spec unless otherwise noted.
// Parameter values come from section 8.8.2 of the spec.
// https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11#section-8.8.2
//
// Base field F is GF(p^m)
// p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
// m = 2 (or 1 for G1 see section 8.8.1)
// k = 128
const htfDefaults = {
  // DST: a domain separation tag
  // defined in section 2.2.5
  // Use utils.getDSTLabel(), utils.setDSTLabel(value)
  DST: "BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_POP_",
  // p: the characteristic of F
  //    where F is a finite field of characteristic p and order q = p^m
  p: CURVE.P,
  // m: the extension degree of F, m >= 1
  //     where F is a finite field of characteristic p and order q = p^m
  m: 2,
  // k: the target security level for the suite in bits
  // defined in section 5.1
  k: 128,
  // option to use a message that has already been processed by
  // expand_message_xmd
  expand: true,
};

function isWithinCurveOrder(num: bigint): boolean {
  return 0 < num && num < CURVE.r;
}

// Global symbol available in browsers only. Ensure we do not depend on @types/dom
declare const self: Record<string, any> | undefined;
const crypto: { node?: any; web?: any } = {
  node: nodeCrypto,
  web: typeof self === "object" && "crypto" in self ? self.crypto : undefined,
};

export const utils = {
  hashToField: hash_to_field,
  /**
   * Can take 40 or more bytes of uniform input e.g. from CSPRNG or KDF
   * and convert them into private key, with the modulo bias being neglible.
   * As per FIPS 186 B.1.1.
   * https://research.kudelskisecurity.com/2020/07/28/the-definitive-guide-to-modulo-bias-and-how-to-avoid-it/
   * @param hash hash output from sha512, or a similar function
   * @returns valid private key
   */
  hashToPrivateKey: (hash: Hex): Uint8Array => {
    hash = ensureBytes(hash);
    if (hash.length < 40 || hash.length > 1024)
      throw new Error("Expected 40-1024 bytes of private key as per FIPS 186");
    const num = mod(bytesToNumberBE(hash), CURVE.r);
    // This should never happen
    if (num === 0n || num === 1n) throw new Error("Invalid private key");
    return numberTo32BytesBE(num);
  },

  bytesToHex,
  randomBytes: (bytesLength: number = 32): Uint8Array => {
    if (crypto.web) {
      return crypto.web.getRandomValues(new Uint8Array(bytesLength));
    } else if (crypto.node) {
      const { randomBytes } = crypto.node;
      return new Uint8Array(randomBytes(bytesLength).buffer);
    } else {
      throw new Error("The environment doesn't have randomBytes function");
    }
  },

  // NIST SP 800-56A rev 3, section 5.6.1.2.2
  // https://research.kudelskisecurity.com/2020/07/28/the-definitive-guide-to-modulo-bias-and-how-to-avoid-it/
  randomPrivateKey: (): Uint8Array => {
    return utils.hashToPrivateKey(utils.randomBytes(40));
  },

  sha256: async (message: Uint8Array): Promise<Uint8Array> => {
    if (crypto.web) {
      const buffer = await crypto.web.subtle.digest("SHA-256", message.buffer);
      return new Uint8Array(buffer);
    } else if (crypto.node) {
      return Uint8Array.from(
        crypto.node.createHash("sha256").update(message).digest()
      );
    } else {
      throw new Error("The environment doesn't have sha256 function");
    }
  },
  mod,
  getDSTLabel() {
    return htfDefaults.DST;
  },
  setDSTLabel(newLabel: string) {
    // https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11#section-3.1
    if (
      typeof newLabel !== "string" ||
      newLabel.length > 2048 ||
      newLabel.length === 0
    ) {
      throw new TypeError("Invalid DST");
    }
    htfDefaults.DST = newLabel;
  },
};

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

function numberTo32BytesBE(num: bigint) {
  const length = 32;
  const hex = num.toString(16).padStart(length * 2, "0");
  return hexToBytes(hex);
}

function toPaddedHex(num: bigint, padding: number) {
  if (typeof num !== "bigint" || num < 0n)
    throw new Error("Expected valid bigint");
  if (typeof padding !== "number")
    throw new TypeError("Expected valid padding");
  return num.toString(16).padStart(padding * 2, "0");
}

function ensureBytes(hex: string | Uint8Array): Uint8Array {
  // Uint8Array.from() instead of hash.slice() because node.js Buffer
  // is instance of Uint8Array, and its slice() creates **mutable** copy
  return hex instanceof Uint8Array ? Uint8Array.from(hex) : hexToBytes(hex);
}

function concatBytes(...arrays: Uint8Array[]): Uint8Array {
  if (arrays.length === 1) return arrays[0];
  const length = arrays.reduce((a, arr) => a + arr.length, 0);
  const result = new Uint8Array(length);
  for (let i = 0, pad = 0; i < arrays.length; i++) {
    const arr = arrays[i];
    result.set(arr, pad);
    pad += arr.length;
  }
  return result;
}

// UTF8 to ui8a
function stringToBytes(str: string) {
  const bytes = new Uint8Array(str.length);
  for (let i = 0; i < str.length; i++) {
    bytes[i] = str.charCodeAt(i);
  }
  return bytes;
}

// Octet Stream to Integer
function os2ip(bytes: Uint8Array): bigint {
  let result = 0n;
  for (let i = 0; i < bytes.length; i++) {
    result <<= 8n;
    result += BigInt(bytes[i]);
  }
  return result;
}

// Integer to Octet Stream
function i2osp(value: number, length: number): Uint8Array {
  if (value < 0 || value >= 1 << (8 * length)) {
    throw new Error(`bad I2OSP call: value=${value} length=${length}`);
  }
  const res = Array.from({ length }).fill(0) as number[];
  for (let i = length - 1; i >= 0; i--) {
    res[i] = value & 0xff;
    value >>>= 8;
  }
  return new Uint8Array(res);
}

function strxor(a: Uint8Array, b: Uint8Array): Uint8Array {
  const arr = new Uint8Array(a.length);
  for (let i = 0; i < a.length; i++) {
    arr[i] = a[i] ^ b[i];
  }
  return arr;
}

// Produces a uniformly random byte string using a cryptographic hash function H that outputs b bits
// https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11#section-5.4.1
async function expand_message_xmd(
  msg: Uint8Array,
  DST: Uint8Array,
  lenInBytes: number
): Promise<Uint8Array> {
  const H = utils.sha256;
  const b_in_bytes = SHA256_DIGEST_SIZE;
  const r_in_bytes = b_in_bytes * 2;

  const ell = Math.ceil(lenInBytes / b_in_bytes);
  if (ell > 255) throw new Error("Invalid xmd length");
  const DST_prime = concatBytes(DST, i2osp(DST.length, 1));
  const Z_pad = i2osp(0, r_in_bytes);
  const l_i_b_str = i2osp(lenInBytes, 2);
  const b = new Array<Uint8Array>(ell);
  const b_0 = await H(
    concatBytes(Z_pad, msg, l_i_b_str, i2osp(0, 1), DST_prime)
  );
  b[0] = await H(concatBytes(b_0, i2osp(1, 1), DST_prime));
  for (let i = 1; i <= ell; i++) {
    const args = [strxor(b_0, b[i - 1]), i2osp(i + 1, 1), DST_prime];
    b[i] = await H(concatBytes(...args));
  }
  const pseudo_random_bytes = concatBytes(...b);
  return pseudo_random_bytes.slice(0, lenInBytes);
}

// hashes arbitrary-length byte strings to a list of one or more elements of a finite field F
// https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11#section-5.3
// Inputs:
// msg - a byte string containing the message to hash.
// count - the number of elements of F to output.
// Outputs:
// [u_0, ..., u_(count - 1)], a list of field elements.
async function hash_to_field(
  msg: Uint8Array,
  count: number,
  options = {}
): Promise<bigint[][]> {
  // if options is provided but incomplete, fill any missing fields with the
  // value in hftDefaults (ie hash to G2).
  const htfOptions = { ...htfDefaults, ...options };
  const log2p = htfOptions.p.toString(2).length;
  const L = Math.ceil((log2p + htfOptions.k) / 8); // section 5.1 of ietf draft link above
  const len_in_bytes = count * htfOptions.m * L;
  const DST = stringToBytes(htfOptions.DST);
  let pseudo_random_bytes = msg;
  if (htfOptions.expand) {
    pseudo_random_bytes = await expand_message_xmd(msg, DST, len_in_bytes);
  }
  const u = new Array(count);
  for (let i = 0; i < count; i++) {
    const e = new Array(htfOptions.m);
    for (let j = 0; j < htfOptions.m; j++) {
      const elm_offset = L * (j + i * htfOptions.m);
      const tv = pseudo_random_bytes.slice(elm_offset, elm_offset + L);
      e[j] = mod(os2ip(tv), htfOptions.p);
    }
    u[i] = e;
  }
  return u;
}

function normalizePrivKey(key: PrivateKey): bigint {
  let int: bigint;
  if (key instanceof Uint8Array && key.length === 32)
    int = bytesToNumberBE(key);
  else if (typeof key === "string" && key.length === 64)
    int = BigInt(`0x${key}`);
  else if (typeof key === "number" && key > 0 && Number.isSafeInteger(key))
    int = BigInt(key);
  else if (typeof key === "bigint" && key > 0n) int = key;
  else throw new TypeError("Expected valid private key");
  int = mod(int, CURVE.r);
  if (!isWithinCurveOrder(int))
    throw new Error("Private key must be 0 < key < CURVE.r");
  return int;
}

function assertType(item: any, type: any) {
  if (!(item instanceof type))
    throw new Error("Expected Fp* argument, not number/bigint");
}

// Point on G1 curve: (x, y)
// We add z because we work with projective coordinates instead of affine x-y: that's much faster.
export class PointG1 extends ProjectivePoint<Fp> {
  static BASE = new PointG1(new Fp(CURVE.Gx), new Fp(CURVE.Gy), Fp.ONE);
  static ZERO = new PointG1(Fp.ONE, Fp.ONE, Fp.ZERO);

  constructor(x: Fp, y: Fp, z: Fp = Fp.ONE) {
    super(x, y, z, Fp);
    assertType(x, Fp);
    assertType(y, Fp);
    assertType(z, Fp);
  }

  static fromHex(bytes: Hex) {
    bytes = ensureBytes(bytes);
    const { P } = CURVE;

    let point;
    if (bytes.length === 48) {
      const compressedValue = bytesToNumberBE(bytes);
      const bflag = mod(compressedValue, POW_2_383) / POW_2_382;
      if (bflag === 1n) {
        return this.ZERO;
      }
      const x = new Fp(mod(compressedValue, POW_2_381));
      const right = x.pow(3n).add(new Fp(CURVE.b)); // y² = x³ + b
      let y = right.sqrt();
      if (!y) throw new Error("Invalid compressed G1 point");
      const aflag = mod(compressedValue, POW_2_382) / POW_2_381;
      if ((y.value * 2n) / P !== aflag) y = y.negate();
      point = new PointG1(x, y);
    } else if (bytes.length === 96) {
      // Check if the infinity flag is set
      if ((bytes[0] & (1 << 6)) !== 0) return PointG1.ZERO;
      const x = bytesToNumberBE(bytes.slice(0, PUBLIC_KEY_LENGTH));
      const y = bytesToNumberBE(bytes.slice(PUBLIC_KEY_LENGTH));
      point = new PointG1(new Fp(x), new Fp(y));
    } else {
      throw new Error("Invalid point G1, expected 48/96 bytes");
    }
    point.assertValidity();
    return point;
  }

  static fromPrivateKey(privateKey: PrivateKey) {
    return this.BASE.multiplyPrecomputed(normalizePrivKey(privateKey));
  }

  toRawBytes(isCompressed = false) {
    return hexToBytes(this.toHex(isCompressed));
  }

  toHex(isCompressed = false) {
    this.assertValidity();
    const { P } = CURVE;
    if (isCompressed) {
      let hex;
      if (this.isZero()) {
        hex = POW_2_383 + POW_2_382;
      } else {
        const [x, y] = this.toAffine();
        const flag = (y.value * 2n) / P;
        hex = x.value + flag * POW_2_381 + POW_2_383;
      }
      return toPaddedHex(hex, PUBLIC_KEY_LENGTH);
    } else {
      if (this.isZero()) {
        // 2x PUBLIC_KEY_LENGTH
        return "4".padEnd(2 * 2 * PUBLIC_KEY_LENGTH, "0"); // bytes[0] |= 1 << 6;
      } else {
        const [x, y] = this.toAffine();
        return (
          toPaddedHex(x.value, PUBLIC_KEY_LENGTH) +
          toPaddedHex(y.value, PUBLIC_KEY_LENGTH)
        );
      }
    }
  }

  assertValidity() {
    if (this.isZero()) return this;
    if (!this.isOnCurve()) throw new Error("Invalid G1 point: not on curve Fp");
    if (!this.isTorsionFree())
      throw new Error("Invalid G1 point: must be of prime-order subgroup");
    return this;
  }

  [Symbol.for("nodejs.util.inspect.custom")]() {
    return this.toString();
  }

  // Sparse multiplication against precomputed coefficients
  millerLoop(P: PointG2): Fp12 {
    return millerLoop(P.pairingPrecomputes(), this.toAffine());
  }

  // Clear cofactor of G1
  // https://eprint.iacr.org/2019/403
  clearCofactor(): PointG1 {
    // return this.multiplyUnsafe(CURVE.h);
    const t = this.mulCurveMinusX();
    return t.add(this);
  }

  // Checks for equation y² = x³ + b
  private isOnCurve(): boolean {
    const b = new Fp(CURVE.b);
    const { x, y, z } = this;
    const left = y.pow(2n).multiply(z).subtract(x.pow(3n));
    const right = b.multiply(z.pow(3n));
    return left.subtract(right).isZero();
  }

  // σ endomorphism
  private sigma(): PointG1 {
    const BETA =
      0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaacn;
    const [x, y] = this.toAffine();
    return new PointG1(x.multiply(BETA), y);
  }

  // φ endomorphism
  private phi(): PointG1 {
    const cubicRootOfUnityModP =
      0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffen;
    return new PointG1(this.x.multiply(cubicRootOfUnityModP), this.y, this.z);
  }

  // [-0xd201000000010000]P
  private mulCurveX() {
    return this.multiplyUnsafe(CURVE.x).negate();
  }
  // [0xd201000000010000]P
  private mulCurveMinusX() {
    return this.multiplyUnsafe(CURVE.x);
  }

  // Checks is the point resides in prime-order subgroup.
  // point.isTorsionFree() should return true for valid points
  // It returns false for shitty points.
  // https://eprint.iacr.org/2021/1130.pdf
  private isTorsionFree(): boolean {
    // todo: unroll
    const xP = this.mulCurveX(); // [x]P
    const u2P = xP.mulCurveMinusX(); // [u2]P
    return u2P.equals(this.phi());

    // https://eprint.iacr.org/2019/814.pdf
    // (z² − 1)/3
    // const c1 = 0x396c8c005555e1560000000055555555n;
    // const P = this;
    // const S = P.sigma();
    // const Q = S.double();
    // const S2 = S.sigma();
    // // [(z² − 1)/3](2σ(P) − P − σ²(P)) − σ²(P) = O
    // const left = Q.subtract(P).subtract(S2).multiplyUnsafe(c1);
    // const C = left.subtract(S2);
    // return C.isZero();
  }
}

// Point on G2 curve (complex numbers): (x₁, x₂+i), (y₁, y₂+i)
// We add z because we work with projective coordinates instead of affine x-y: that's much faster.
export class PointG2 extends ProjectivePoint<Fp2> {
  static BASE = new PointG2(
    Fp2.fromBigTuple(CURVE.G2x),
    Fp2.fromBigTuple(CURVE.G2y),
    Fp2.ONE
  );
  static ZERO = new PointG2(Fp2.ONE, Fp2.ONE, Fp2.ZERO);

  private _PPRECOMPUTES: [Fp2, Fp2, Fp2][] | undefined;

  constructor(x: Fp2, y: Fp2, z: Fp2 = Fp2.ONE) {
    super(x, y, z, Fp2);
    assertType(x, Fp2);
    assertType(y, Fp2);
    assertType(z, Fp2);
  }

  // Encodes byte string to elliptic curve
  // https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11#section-3
  static async hashToCurve(msg: Hex) {
    msg = ensureBytes(msg);
    const u = await hash_to_field(msg, 2);
    //console.log(`hash_to_curve(msg}) u0=${new Fp2(u[0])} u1=${new Fp2(u[1])}`);
    const Q0 = new PointG2(
      ...isogenyMapG2(map_to_curve_simple_swu_9mod16(u[0]))
    );
    const Q1 = new PointG2(
      ...isogenyMapG2(map_to_curve_simple_swu_9mod16(u[1]))
    );
    const R = Q0.add(Q1);
    const P = R.clearCofactor();
    //console.log(`hash_to_curve(msg) Q0=${Q0}, Q1=${Q1}, R=${R} P=${P}`);
    return P;
  }

  // TODO: Optimize, it's very slow because of sqrt.
  static fromSignature(hex: Hex): PointG2 {
    hex = ensureBytes(hex);
    const { P } = CURVE;
    const half = hex.length / 2;
    if (half !== 48 && half !== 96)
      throw new Error("Invalid compressed signature length, must be 96 or 192");
    const z1 = bytesToNumberBE(hex.slice(0, half));
    const z2 = bytesToNumberBE(hex.slice(half));
    // Indicates the infinity point
    const bflag1 = mod(z1, POW_2_383) / POW_2_382;
    if (bflag1 === 1n) return this.ZERO;

    const x1 = new Fp(z1 % POW_2_381);
    const x2 = new Fp(z2);
    const x = new Fp2(x2, x1);
    const y2 = x.pow(3n).add(Fp2.fromBigTuple(CURVE.b2)); // y² = x³ + 4
    // The slow part
    let y = y2.sqrt();
    if (!y) throw new Error("Failed to find a square root");

    // Choose the y whose leftmost bit of the imaginary part is equal to the a_flag1
    // If y1 happens to be zero, then use the bit of y0
    const { re: y0, im: y1 } = y.reim();
    const aflag1 = (z1 % POW_2_382) / POW_2_381;
    const isGreater = y1 > 0n && (y1 * 2n) / P !== aflag1;
    const isZero = y1 === 0n && (y0 * 2n) / P !== aflag1;
    if (isGreater || isZero) y = y.multiply(-1n);
    const point = new PointG2(x, y, Fp2.ONE);
    point.assertValidity();
    return point;
  }

  static fromHex(bytes: Hex) {
    bytes = ensureBytes(bytes);
    let point;
    if (bytes.length === 96) {
      throw new Error("Compressed format not supported yet.");
    } else if (bytes.length === 192) {
      // Check if the infinity flag is set
      if ((bytes[0] & (1 << 6)) !== 0) {
        return PointG2.ZERO;
      }

      const x1 = bytesToNumberBE(bytes.slice(0, PUBLIC_KEY_LENGTH));
      const x0 = bytesToNumberBE(
        bytes.slice(PUBLIC_KEY_LENGTH, 2 * PUBLIC_KEY_LENGTH)
      );
      const y1 = bytesToNumberBE(
        bytes.slice(2 * PUBLIC_KEY_LENGTH, 3 * PUBLIC_KEY_LENGTH)
      );
      const y0 = bytesToNumberBE(bytes.slice(3 * PUBLIC_KEY_LENGTH));

      point = new PointG2(
        Fp2.fromBigTuple([x0, x1]),
        Fp2.fromBigTuple([y0, y1])
      );
    } else {
      throw new Error("Invalid uncompressed point G2, expected 192 bytes");
    }

    point.assertValidity();
    return point;
  }

  static fromPrivateKey(privateKey: PrivateKey) {
    return this.BASE.multiplyPrecomputed(normalizePrivKey(privateKey));
  }

  toSignature() {
    if (this.equals(PointG2.ZERO)) {
      const sum = POW_2_383 + POW_2_382;
      const h =
        toPaddedHex(sum, PUBLIC_KEY_LENGTH) +
        toPaddedHex(0n, PUBLIC_KEY_LENGTH);
      return hexToBytes(h);
    }
    const [{ re: x0, im: x1 }, { re: y0, im: y1 }] = this.toAffine().map((a) =>
      a.reim()
    );
    const tmp = y1 > 0n ? y1 * 2n : y0 * 2n;
    const aflag1 = tmp / CURVE.P;
    const z1 = x1 + aflag1 * POW_2_381 + POW_2_383;
    const z2 = x0;
    return hexToBytes(
      toPaddedHex(z1, PUBLIC_KEY_LENGTH) + toPaddedHex(z2, PUBLIC_KEY_LENGTH)
    );
  }

  toRawBytes(isCompressed = false) {
    return hexToBytes(this.toHex(isCompressed));
  }

  toHex(isCompressed = false) {
    this.assertValidity();
    if (isCompressed) {
      throw new Error("Point compression has not yet been implemented");
    } else {
      if (this.equals(PointG2.ZERO)) {
        return "4".padEnd(2 * 4 * PUBLIC_KEY_LENGTH, "0"); // bytes[0] |= 1 << 6;
      }
      const [{ re: x0, im: x1 }, { re: y0, im: y1 }] = this.toAffine().map(
        (a) => a.reim()
      );
      return (
        toPaddedHex(x1, PUBLIC_KEY_LENGTH) +
        toPaddedHex(x0, PUBLIC_KEY_LENGTH) +
        toPaddedHex(y1, PUBLIC_KEY_LENGTH) +
        toPaddedHex(y0, PUBLIC_KEY_LENGTH)
      );
    }
  }

  assertValidity() {
    if (this.isZero()) return this;
    if (!this.isOnCurve())
      throw new Error("Invalid G2 point: not on curve Fp2");
    if (!this.isTorsionFree())
      throw new Error("Invalid G2 point: must be of prime-order subgroup");
    return this;
  }

  // Ψ endomorphism
  private psi() {
    return this.fromAffineTuple(psi(...this.toAffine()));
  }

  // Ψ²
  private psi2() {
    return this.fromAffineTuple(psi2(...this.toAffine()));
  }

  // [-x]P aka [z]P
  private mulCurveX() {
    return this.multiplyUnsafe(CURVE.x).negate();
  }

  // Maps the point into the prime-order subgroup G2.
  // clear_cofactor_bls12381_g2 from cfrg-hash-to-curve-11
  // https://eprint.iacr.org/2017/419.pdf
  // prettier-ignore
  clearCofactor(): PointG2 {
    const P = this;
    let t1 = P.mulCurveX();   // [-x]P
    let t2 = P.psi();         // Ψ(P)
    let t3 = P.double();      // 2P
    t3 = t3.psi2();           // Ψ²(2P)
    t3 = t3.subtract(t2);     // Ψ²(2P) - Ψ(P)
    t2 = t1.add(t2);          // [-x]P + Ψ(P)
    t2 = t2.mulCurveX();      // [x²]P - [x]Ψ(P)
    t3 = t3.add(t2);          // Ψ²(2P) - Ψ(P) + [x²]P - [x]Ψ(P)
    t3 = t3.subtract(t1);     // Ψ²(2P) - Ψ(P) + [x²]P - [x]Ψ(P) + [x]P
    const Q = t3.subtract(P); // Ψ²(2P) - Ψ(P) + [x²]P - [x]Ψ(P) + [x]P - 1P =>
    return Q;                 // [x²-x-1]P + [x-1]Ψ(P) + Ψ²(2P)
  }

  // Checks for equation y² = x³ + b
  private isOnCurve(): boolean {
    const b = Fp2.fromBigTuple(CURVE.b2);
    const { x, y, z } = this;
    const left = y.pow(2n).multiply(z).subtract(x.pow(3n));
    const right = b.multiply(z.pow(3n) as Fp2);
    return left.subtract(right).isZero();
  }

  // Checks is the point resides in prime-order subgroup.
  // point.isTorsionFree() should return true for valid points
  // It returns false for shitty points.
  // https://eprint.iacr.org/2021/1130.pdf
  // prettier-ignore
  private isTorsionFree(): boolean {
    const P = this;
    return P.mulCurveX().equals(P.psi()); // ψ(P) == [u](P)
    // https://eprint.iacr.org/2019/814.pdf
    // const psi2 = P.psi2();                        // Ψ²(P)
    // const psi3 = psi2.psi();                      // Ψ³(P)
    // const zPsi3 = psi3.mulNegX();                 // [z]Ψ³(P) where z = -x
    // return zPsi3.subtract(psi2).add(P).isZero();  // [z]Ψ³(P) - Ψ²(P) + P == O
  }

  // Improves introspection in node.js. Basically displays point's x, y.
  [Symbol.for("nodejs.util.inspect.custom")]() {
    return this.toString();
  }

  clearPairingPrecomputes() {
    this._PPRECOMPUTES = undefined;
  }

  pairingPrecomputes(): [Fp2, Fp2, Fp2][] {
    if (this._PPRECOMPUTES) return this._PPRECOMPUTES;
    this._PPRECOMPUTES = calcPairingPrecomputes(...this.toAffine());
    return this._PPRECOMPUTES;
  }
}

// Calculates bilinear pairing
export function pairing(
  P: PointG1,
  Q: PointG2,
  withFinalExponent: boolean = true
): Fp12 {
  if (P.isZero() || Q.isZero())
    throw new Error("No pairings at point of Infinity");
  P.assertValidity();
  Q.assertValidity();
  // Performance: 9ms for millerLoop and ~14ms for exp.
  const looped = P.millerLoop(Q);
  return withFinalExponent ? looped.finalExponentiate() : looped;
}

type G1Hex = Hex | PointG1;
type G2Hex = Hex | PointG2;
function normP1(point: G1Hex): PointG1 {
  return point instanceof PointG1 ? point : PointG1.fromHex(point);
}
function normP2(point: G2Hex): PointG2 {
  return point instanceof PointG2 ? point : PointG2.fromSignature(point);
}
async function normP2Hash(point: G2Hex): Promise<PointG2> {
  return point instanceof PointG2 ? point : PointG2.hashToCurve(point);
}

// Multiplies generator by private key.
// P = pk x G
export function getPublicKey(privateKey: PrivateKey): Uint8Array {
  return PointG1.fromPrivateKey(privateKey).toRawBytes(true);
}

// Executes `hashToCurve` on the message and then multiplies the result by private key.
// S = pk x H(m)
export async function sign(
  message: Hex,
  privateKey: PrivateKey
): Promise<Uint8Array>;
export async function sign(
  message: PointG2,
  privateKey: PrivateKey
): Promise<PointG2>;
export async function sign(
  message: G2Hex,
  privateKey: PrivateKey
): Promise<Uint8Array | PointG2> {
  const msgPoint = await normP2Hash(message);
  msgPoint.assertValidity();
  const sigPoint = msgPoint.multiply(normalizePrivKey(privateKey));
  if (message instanceof PointG2) return sigPoint;
  return sigPoint.toSignature();
}

// Checks if pairing of public key & hash is equal to pairing of generator & signature.
// e(P, H(m)) == e(G, S)
export async function verify(
  signature: G2Hex,
  message: G2Hex,
  publicKey: G1Hex
): Promise<boolean> {
  const P = normP1(publicKey);
  const Hm = await normP2Hash(message);
  const G = PointG1.BASE;
  const S = normP2(signature);
  // Instead of doing 2 exponentiations, we use property of billinear maps
  // and do one exp after multiplying 2 points.
  const ePHm = pairing(P.negate(), Hm, false);
  const eGS = pairing(G, S, false);
  const exp = eGS.multiply(ePHm).finalExponentiate();
  return exp.equals(Fp12.ONE);
}

// Adds a bunch of public key points together.
// pk1 + pk2 + pk3 = pkA
export function aggregatePublicKeys(publicKeys: Hex[]): Uint8Array;
export function aggregatePublicKeys(publicKeys: PointG1[]): PointG1;
export function aggregatePublicKeys(publicKeys: G1Hex[]): Uint8Array | PointG1 {
  if (!publicKeys.length) throw new Error("Expected non-empty array");
  const agg = publicKeys
    .map(normP1)
    .reduce((sum, p) => sum.add(p), PointG1.ZERO);
  if (publicKeys[0] instanceof PointG1) return agg.assertValidity();
  return agg.toRawBytes(true);
}

// Adds a bunch of signature points together.
export function aggregateSignatures(signatures: Hex[]): Uint8Array;
export function aggregateSignatures(signatures: PointG2[]): PointG2;
export function aggregateSignatures(signatures: G2Hex[]): Uint8Array | PointG2 {
  if (!signatures.length) throw new Error("Expected non-empty array");
  const agg = signatures
    .map(normP2)
    .reduce((sum, s) => sum.add(s), PointG2.ZERO);
  if (signatures[0] instanceof PointG2) return agg.assertValidity();
  return agg.toSignature();
}

// https://ethresear.ch/t/fast-verification-of-multiple-bls-signatures/5407
// e(G, S) = e(G, SUM(n)(Si)) = MUL(n)(e(G, Si))
export async function verifyBatch(
  signature: G2Hex,
  messages: G2Hex[],
  publicKeys: G1Hex[]
): Promise<boolean> {
  if (!messages.length) throw new Error("Expected non-empty messages array");
  if (publicKeys.length !== messages.length)
    throw new Error("Pubkey count should equal msg count");
  const sig = normP2(signature);
  const nMessages = await Promise.all(messages.map(normP2Hash));
  const nPublicKeys = publicKeys.map(normP1);
  try {
    const paired = [];
    for (const message of new Set(nMessages)) {
      const groupPublicKey = nMessages.reduce(
        (groupPublicKey, subMessage, i) =>
          subMessage === message
            ? groupPublicKey.add(nPublicKeys[i])
            : groupPublicKey,
        PointG1.ZERO
      );
      // const msg = message instanceof PointG2 ? message : await PointG2.hashToCurve(message);
      // Possible to batch pairing for same msg with different groupPublicKey here
      paired.push(pairing(groupPublicKey, message, false));
    }
    paired.push(pairing(PointG1.BASE.negate(), sig, false));
    const product = paired.reduce((a, b) => a.multiply(b), Fp12.ONE);
    const exp = product.finalExponentiate();
    return exp.equals(Fp12.ONE);
  } catch {
    return false;
  }
}

// Pre-compute points. Refer to README.
PointG1.BASE.calcMultiplyPrecomputes(4);
