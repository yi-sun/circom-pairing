cd dev2_cpp
./dev2 ../inputs/"$1".json ../inputs/"$1".wtns
cd ..
/home/vincent/rapidsnark/build/prover dev2.zkey inputs/"$1".wtns proof_"$1".json public_"$1".json
cat proof_"$1".json