#!/bin/bash

for file in $(ls *.circom);
do
    sed -i "" -e 's/<==/<--/g' $file
    sed -i "" -e '/===/d' $file
done

