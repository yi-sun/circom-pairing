#!/bin/bash

PHASE1=../../circuits/pot20.ptau
BUILD_DIR=../../build/dev
CIRCUIT_NAME=dev

if [ -f "$PHASE1" ]; then
    echo "Found Phase 1 ptau file"
else
    echo "No Phase 1 ptau file found. Exiting..."
    exit 1
fi

if [ ! -d "$BUILD_DIR" ]; then
    echo "No build directory found. Creating build directory..."
    mkdir "$BUILD_DIR"
fi

echo $PWD

echo "****COMPILING CIRCUIT****"
start=`date +%s`
#circom "$CIRCUIT_NAME".circom --O1 --r1cs --wasm --sym --c --wat --output "$BUILD_DIR"
circom "$CIRCUIT_NAME".circom --O1 --c --sym --output "$BUILD_DIR"
end=`date +%s`
echo "DONE ($((end-start))s)"

echo "****COMPILING C++ WITNESS GENERATION CODE****"
start=`date +%s`
cd "$BUILD_DIR"/"$CIRCUIT_NAME"_cpp 
make
end=`date +%s`
echo "DONE ($((end-start))s)"

echo "****VERIFYING WITNESS****"
start=`date +%s`
./"$CIRCUIT_NAME" ../../../scripts/"$CIRCUIT_NAME"/input_"$CIRCUIT_NAME".json ../witness.wtns
end=`date +%s`
echo "DONE ($((end-start))s)"

cd ..
npx snarkjs wej witness.wtns witness.json

echo "****GENERATING ZKEY 0****"
start=`date +%s`
NODE_OPTIONS="--max-old-space-size=56000" npx snarkjs groth16 setup "$BUILD_DIR"/"$CIRCUIT_NAME".r1cs "$PHASE1" "$BUILD_DIR"/"$CIRCUIT_NAME"_0.zkey
end=`date +%s`
echo "DONE ($((end-start))s)"

echo "****GENERATING FINAL ZKEY****"
start=`date +%s`
NODE_OPTIONS="--max-old-space-size=56000" npx snarkjs zkey beacon "$BUILD_DIR"/"$CIRCUIT_NAME"_0.zkey "$BUILD_DIR"/"$CIRCUIT_NAME".zkey 0102030405060708090a0b0c0d0e0f101112231415161718221a1b1c1d1e1f 10 -n="Final Beacon phase2"
end=`date +%s`
echo "DONE ($((end-start))s)"

echo "****VERIFYING FINAL ZKEY****"
start=`date +%s`
NODE_OPTIONS="--max-old-space-size=56000" npx snarkjs zkey verify -verbose "$BUILD_DIR"/"$CIRCUIT_NAME".r1cs "$PHASE1" "$BUILD_DIR"/"$CIRCUIT_NAME".zkey
end=`date +%s`
echo "DONE ($((end-start))s)"

echo "****EXPORTING VKEY****"
start=`date +%s`
npx snarkjs zkey export verificationkey "$BUILD_DIR"/"$CIRCUIT_NAME".zkey "$BUILD_DIR"/vkey.json
end=`date +%s`
echo "DONE ($((end-start))s)"

echo "****GENERATING PROOF FOR SAMPLE INPUT****"
start=`date +%s`
npx snarkjs groth16 prove "$BUILD_DIR"/"$CIRCUIT_NAME".zkey "$BUILD_DIR"/witness.wtns "$BUILD_DIR"/proof.json "$BUILD_DIR"/public.json
end=`date +%s`
echo "DONE ($((end-start))s)"

echo "****VERIFYING PROOF FOR SAMPLE INPUT****"
start=`date +%s`
npx snarkjs groth16 verify "$BUILD_DIR"/vkey.json "$BUILD_DIR"/public.json "$BUILD_DIR"/proof.json
end=`date +%s`
echo "DONE ($((end-start))s)"
