#!/bin/bash
sudo apt update
sudo apt install -y python3-numpy

python3 RunBenchmark.py 72
python3 RunBenchmark-NonAMX.py 72

python3 RunBenchmark.py 144
python3 RunBenchmark-NonAMX.py 144

python3 EnergyBenchmark.py 72
python3 EnergyBenchmark.py 144