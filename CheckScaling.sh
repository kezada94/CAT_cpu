#!/bin/bash

make clean
make -j USE_AVX512=USE_AVX512

# get the maximum number of physical cores
max_cores=$(lscpu | grep "Core(s) per socket" | awk '{print $4}')

for m in 0 1 2 3 4 5 6
do
    echo "Running benchmarks with n=$((2**15)) and m=$m"
    for i in 1 2 4 8 16 32 64 $max_cores
    do
        ./bin/prog $((2**15)) 1 10 --density 0.07 --threads $i
    done
done