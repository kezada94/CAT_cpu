#include "CPUBenchmark.h"
#include "Rapl.h"

CPUBenchmark::CPUBenchmark(CASolver *pSolver, int n, int pSteps, int pSeed, float pDensity, int pThreads)
{
    solver = pSolver;
    this->n = n;
    steps = pSteps;
    seed = pSeed;
    threads = pThreads;
    density = pDensity;

    timer = new CpuTimer();
    stats = new StatsCollector();

	lDebug(1, "Setting the number of threads to %i", pThreads);
	omp_set_num_threads(pThreads);

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

    lDebug(1, "Initial state:");
    if (n <= PRINT_LIMIT)
    {
        fDebug(1, solver->printCurrentState());
    }

#ifdef MEASURE_POWER

    std::string radiusStr = "AMX-radius-" + std::to_string(RADIUS);

    // Call CPUPowerBegin with the correct string format
    CPUPowerBegin(radiusStr.c_str(), 50);

#endif

    doOneRun();

#ifdef MEASURE_POWER
    CPUPowerEnd();
#endif

    lDebug(1, "Benchmark finished. Results:");

    if (n <= PRINT_LIMIT)
    {
        fDebug(1, solver->printCurrentState());
    }
}
void CPUBenchmark::doOneRun()
{

	#pragma omp parallel
    {

		int tid = omp_get_thread_num();
		
		if (tid == 0){
			timer->start();
		}
		

		solver->doSteps(steps);

		#pragma omp barrier
		if (tid == 0){
			timer->stop();
		}
	}

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
