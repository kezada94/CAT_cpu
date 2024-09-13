#include "CPUBenchmark.h"

CPUBenchmark::CPUBenchmark(CASolver *pSolver, int n, int pSteps, int pSeed, float pDensity)
{
    solver = pSolver;
    this->n = n;
    steps = pSteps;
    seed = pSeed;
    density = pDensity;

    timer = new CpuTimer();
    stats = new StatsCollector();
}

void CPUBenchmark::reset()
{
    int s = rand() % 1000000;

    lDebug(1, "Resetting state:");
    solver->resetState(s, density);
}

void CPUBenchmark::run()
{
    srand(seed);
    reset();
    lDebug(1, "Benchmark started");
    // WARMUP for STEPS/4
    // solver->doSteps(steps >> 2);

    lDebug(1, "Initial state:");
    if (n <= PRINT_LIMIT)
    {
        fDebug(1, solver->printCurrentState());
    }


    doOneRun();


    lDebug(1, "Benchmark finished. Results:");

    if (n <= PRINT_LIMIT)
    {
        fDebug(1, solver->printCurrentState());
    }
}
void CPUBenchmark::doOneRun()
{
    timer->start();

    solver->doSteps(steps);

    timer->stop();
    registerElapsedTime(timer->getElapsedTimeMilliseconds() / steps);
}

void CPUBenchmark::registerElapsedTime(float milliseconds)
{
    stats->add(milliseconds);
}

StatsCollector *CPUBenchmark::getStats()
{
    return stats;
}
