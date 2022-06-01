# zkPairing

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Project overview](#project-overview)
- [Setup](#setup)
- [Building keys and witness generation files](#building-keys-and-witness-generation-files)
- [Benchmarks](#benchmarks)
- [Testing](#testing)
- [Demo](#demo)
- [Acknowledgments](#acknowledgments)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Project overview

This repository provides proof-of-concept implementations of elliptic curve pairings (in particular, the optimal Ate pairing and Tate pairing) for the BLS12-381 curve in circom. These implementations are for demonstration purposes only. These circuits are not audited, and this is not intended to be used as a library for production-grade applications.

Circuits can be found in the `circuits` directory. The `scripts` directory contains various utility scripts (most importantly, a script for building a zkSNARK to verify BLS signatures). `test` contains some unit tests for the circuits, mostly for witness generation.

## Setup

First, install [yarn](https://classic.yarnpkg.com/en/) and [circom](https://docs.circom.io/getting-started/installation/).

- run `yarn install` in the root directory to install the dependencies (`snarkjs` and `circomlib`) in `yarn.lock`.
- You'll need `circom` version `>= 2.0.3`.
- If you want to build the `optimalate`, `subgroupcheckG1`, `subgroupcheckG2` circuits, you'll need to download a Powers of Tau file with `2^24` constraints and copy it into the circuits subdirectory of the project, with the name `pot24_final.ptau`. We do not provide such a file in this repo due to its large size. You can download and copy Powers of Tau files from the Hermez trusted setup from [this repository](https://github.com/iden3/snarkjs#7-prepare-phase-2).
- If you want to build the `verify` and `tatepairing` circuits, you'll need a Powers of Tau file that can support at least `2^25` constraints (place it in the same directory as above with the same naming convention).

## Building keys and witness generation files

We provide the following circuits as examples:

- `verify`: Prove that a BLS signature verification ran properly on a provided public key, signature, and message. The circuit verifies that the public key and signature are valid.
- `optimalate`: Prove that the optimal Ate pairing is correctly computed on two elements in appropriate subgroups.
- `tatepairing`: Prove that the Tate pairing is correctly computed on two elements in appropriate subgroups.
- `subgroupcheckG1`: Prove that a public key is valid, i.e., lies in the subgroup `G1` of the curve.
- `subgroupcheckG2`: Prove that a signature is valid, i.e., lies in the subgroup `G2` of the curve.
- `maptoG2`: Given two `F_{p^2}` elements, prove that their mapping to the twisted curve and cofactor clearing to `G2` is correctly computed according to the [Internet Draft](https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-14).

Run `yarn build:verify`, `yarn build:optimalate`, `yarn build:tatepairing`, etc. at the top level to compile each respective circuit and keys. See [documentation](docs/README.md) for input format.

Note that `verify` and `tatepairing` are very large circuits so they require special hardware and setup to run: see [Best Practices for Large Circuits](https://hackmd.io/V-7Aal05Tiy-ozmzTGBYPA?view).

## Benchmarks

All benchmarks were run on a 32-core 3.1GHz, 256G RAM machine with 1TB hard drive (AWS r5.8xlarge instance). Constraints refer to non-linear constraints.

|                                      | verify | optimalate | tatepairing | maptoG2 | subgroupcheckG1 | subgroupcheckG2 |
| ------------------------------------ | ------ | ---------- | ----------- | ------- | --------------- | --------------- |
| Constraints                          | 19.2M  | 11.4M      | 24.7M       | 2M      | 789K            | 819K            |
| Circuit compilation                  | 3.2h   | 1.9h       | 4.2h        | 23m     | 7.6m            | 8.5m            |
| Witness generation C++ compilation   | 2h     | 1.1h       | 2.3h        | 9.3m    | 4.2m            | 3.8m            |
| Witness generation                   | 2.6m   | 1m         | 2.5m        | 33s     | 23s             | 13s             |
| Trusted setup phase 2 key generation | 58m    | 32m        | 1.6h        | 4.5m    | 1.7m            | 1.9m            |
| Trusted setup phase 2 contribution   | 25m    | 13.6m      | 29m         | 2.9m    | 54s             | 55s             |
| Proving key size                     | 12G    | 6.5G       | 15G         | 1.2G    | 421M            | 445M            |
| Proving key verification             | 1.5h   | 43m        | 2.5h        | 6.2m    | 2m              | 2.3m            |
| Proving time (rapidsnark)            | 2m     | 52s        | 2.1m        | 6s      | 3s              | 3s              |
| Proof verification time              | 1s     | 1s         | 2s          | 1s      | 1s              | 1s              |

## Testing

See the `/test` directory for examples of tests. The circuits to be tested should be written in the `/test/circuits` folder, while the test execution code should be written in regular JavaScript files under `/test`. A short description of each test can be passed in as the first parameter of the `describe()` function, and `yarn --grep name` will run all tests whose description contains `name` as a substring.

## Documentation

See [documentation](docs/README.md) for documentation of all circuits.

## Demo

See [here](https://zkpairing.xyz) for a demo of BLS signature verification inside a zk-SNARK. The frontend code for the demo can be found [here](https://github.com/vincenthuang75025/zkxzk), and the server code can be found [here](https://github.com/vincenthuang75025/zk-node-server-c). 

## Acknowledgments

This project was built during [0xPARC](http://0xparc.org/)'s ZK-ID Working Group.

We use a circom bigint library from [circom-ECDSA](https://github.com/0xPARC/circom-ecdsa) and implement many of the same optimizations for elliptic curve operations as they do. This library uses an optimization for big integer multiplication from [xJsnark](https://github.com/akosba/xjsnark).
