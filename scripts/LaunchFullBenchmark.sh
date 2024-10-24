#!/bin/bash
# Receives the number of threads as an argument and the number of sockets
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <number_of_threads> <number_of_sockets>"
    exit 1
fi

NUM_THREADS=$1
NUM_SOCKETS=$2

echo "Running benchmarks with $NUM_THREADS threads and $NUM_SOCKETS sockets"

total_threads=$((NUM_THREADS * NUM_SOCKETS))

sudo apt update
sudo apt install -y python3-numpy

python3 RunBenchmark.py $NUM_THREADS
python3 RunBenchmark-NonAMX.py $NUM_THREADS

python3 RunBenchmark.py $total_threads
python3 RunBenchmark-NonAMX.py $total_threads

python3 EnergyBenchmark.py $NUM_THREADS
python3 EnergyBenchmark.py $total_threads