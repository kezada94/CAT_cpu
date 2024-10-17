import numpy as np
import os
import subprocess
import json
import sys


if len(sys.argv) != 2:
    print("Run with args: <nThreads>\n")
    exit()

nThreads = sys.argv[1]
sizes = [1024 + 2048*i for i in range(30)]
#sizes = [60414]
methods = [3]
method_names = ["AVX"]
#methods = [1]
#method_names = ['AMX64']
radiuses = [i for i in range(1,17)]
smin = [2, 7,  15, 40, 35, 49, 101,  163, 108, 122, 156, 170, 213, 245, 170, 170]
smax = [3, 12, 23, 80, 59, 81, 201,  223, 181, 211, 265, 296, 364, 420, 296, 296]
bmin = [3, 8,  14, 41, 34, 46, 75,   74,  100, 123, 147, 170, 203, 234, 170, 170]
bmax = [3, 11, 17, 80, 45, 65, 170,  252, 140, 170, 205, 240, 283, 326, 240, 300]
densities = [0.07, 0.15, 0.25, 0.50, 0.21, 0.22, 0.29, 0.23, 0.24, 0.25, 0.24, 0.25, 0.25, 0.25, 0.28, 0.26]
repeats = [32, 32, 32, 32, 16, 16, 16, 16, 8, 8,
            8,  8,  8,  8,  4,  4,  4,  4, 4, 4,
            4,  4,  4,  4,  2,  2,  2,  2, 2, 2]
steps = 15
# 1: passed
# 0: failed
results = {}
auxdict = {}

for r, radius in enumerate(radiuses):
    for k, method in enumerate(methods):
        print("Cleaning...")
        subprocess.run(['make', 'clean'], stdout=None, stderr=None, cwd="../")
        print(f"Compiling... RADIUS: {radius}")
        print('make', '-j', 'RADIUS='+str(radius), 'SMIN='+str(smin[r]), 'SMAX='+str(smax[r]), 'BMIN='+str(bmin[r]), 'BMAX='+str(bmax[r]))
        subprocess.run(['make', '-j', 'RADIUS='+str(radius), 'SMIN='+str(smin[r]), 'SMAX='+str(smax[r]), 'BMIN='+str(bmin[r]), 'BMAX='+str(bmax[r])], stdout=None, cwd="../")
        for l, size in enumerate(sizes):
            if l<29 or l>29:
                continue
            runs = np.zeros(repeats[l])
            for rep in range(repeats[l]):
                print(f"    Running {rep}... size: {size}, method: {method}, steps: {steps}")
                print(['../bin/prog', str(size), str(method), str(steps), "--density " + str(densities[r]), "--threads " + nThreads])
                result = subprocess.run(['../bin/prog', str(size), str(method), str(steps), "--density", str(densities[r]), "--threads", str(nThreads)], stdout=subprocess.PIPE).stdout.decode('utf-8')
                print(result)
                if not 'fault' in result:
                    runs[rep] = float(result)
                else:
                    exit()
            with open("../data/benchmark_results-"+str(method_names[k])+".txt","a") as data:
                res = {'radius' : radius, 'method' : method, 'size' : size}
                res['time'] = np.average(runs)
                res['var'] = np.var(runs)
                res['std'] = np.std(runs)
                res['stderr'] = res['std'] / np.sqrt(repeats[l])
                res['PSE'] = res['stderr'] / res['time'] * 100.0
                data.write(str(res) + '\n')

