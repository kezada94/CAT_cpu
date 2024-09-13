#pragma once

#include "CellularAutomata/Solvers/CASolver.h"
#include "CpuTimer.h"
#include "Debug.h"
#include "Defines.h"
#include "StatsCollector.hpp"
#include <random>
#define PRINT_LIMIT (512)

class CPUBenchmark
{
  private:
    int steps;
    int n;
    int seed;
    float density;

    CASolver *solver;
    CpuTimer *timer;

    StatsCollector *stats;

    void doOneRun();
    void registerElapsedTime(float milliseconds);
    void reset();

  public:
    CPUBenchmark(CASolver *pSolver, int n, int pSteps, int pSeed, float pDensity);

    void run();
    StatsCollector *getStats();
};
