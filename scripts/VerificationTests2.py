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
methods = [1, 2, 3]
method_names = ['AMX64', 'AMX16', "AVX"]
#methods = [0, 1, 2, 3]
#method_names = ['CPU', 'AMX64', 'AMX16', "AVX"]
#methods = [1]
#method_names = ['AMX64']
radiuses = [i for i in range(1,17)]
smin = [2, 7,  15, 40, 35, 49, 101,  163, 108, 122, 156, 170, 213, 245, 170, 170]
smax = [3, 12, 23, 80, 59, 81, 201,  223, 181, 211, 265, 296, 364, 420, 296, 296]
bmin = [3, 8,  14, 41, 34, 46, 75,   74,  100, 123, 147, 170, 203, 234, 170, 170]
bmax = [3, 11, 17, 80, 45, 65, 170,  252, 140, 170, 205, 240, 283, 326, 240, 300]
densities = [0.07, 0.15, 0.25, 0.50, 0.21, 0.22, 0.29, 0.23, 0.24, 0.25, 0.24, 0.25, 0.25, 0.25, 0.28, 0.26]
steps = 2
# 1: passed
# 0: failed
results = {}
auxdict = {}

for r, radius in enumerate(radiuses):
    if r<14:
        continue
    print("Cleaning...")
    subprocess.run(['make', 'clean'], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, cwd="../")
    print(f"Compiling... RADIUS: {radius}")
    print('make', 'debug', '-j', '8', 'RADIUS='+str(radius), 'SMIN='+str(smin[r]), 'SMAX='+str(smax[r]), 'BMIN='+str(bmin[r]), 'BMAX='+str(bmax[r]))
    subprocess.run(['make', 'debug', '-j', '8', 'RADIUS='+str(radius), 'SMIN='+str(smin[r]), 'SMAX='+str(smax[r]), 'BMIN='+str(bmin[r]), 'BMAX='+str(bmax[r])], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL, cwd="../")
    for k, method in enumerate(methods):
        for l, size in enumerate(sizes):
            if l<0 or l>0:
                continue
            print(['../debug/prog', str(size), str(method), str(steps), "--density " + str(densities[r]), "--threads " + nThreads], "--doVerify")
            result = subprocess.run(['../debug/prog', str(size), str(method), str(steps), "--density", str(densities[r]), "--threads", str(nThreads), "--doVerify"], stdout=subprocess.PIPE).stdout.decode('utf-8')
            if not 'successful' in result:
                print(result)
            else:
                print('Success')

