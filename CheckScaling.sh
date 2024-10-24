#!/bin/bash

make clean
make -j USE_AVX512=USE_AVX512

# get the maximum number of physical cores
max_cores=$(lscpu | grep "Core(s) per socket" | awk '{print $4}')
echo "Cores $max_cores"

n=$((2**15))
for m in 0 1 2 3 4 5 6
do
    echo "Running benchmarks with n=$n and m=$m"

    i=1
    while [ $i -le $max_cores ]; do
	    echo "$i ->  $(./bin/prog $n $m 5 --density 0.07 --threads $i)  -- './bin/prog $n $m 5 --density 0.07 --threads $i'"
	    i=$((i * 2))
    done
    i=$max_cores
    echo "$i ->  $(./bin/prog $n $m 5 --density 0.07 --threads $i)  -- './bin/prog $n $m 5 --density 0.07 --threads $i'"
done
