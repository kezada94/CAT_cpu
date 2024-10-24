#include "CellularAutomata/Solvers/CASolver.h"

void CASolver::resetState(int seed, float density)
{
    CAStateGenerator::generateRandomState(hostVisibleData, seed, density);
    copyHostVisibleDataToCurrentState();
}

CADataDomain<int> *CASolver::getCurrentState()
{
    copyCurrentStateToHostVisibleData();
    return hostVisibleData;
}

void CASolver::doSteps(int stepNumber)
{
    preamble();
    for (int i = 0; i < stepNumber; i++)
    {
        doStep();
    }
}
void CASolver::doStep()
{
    fillBoundaryConditions();
#pragma omp barrier
    CAStepAlgorithm();
#pragma omp barrier
    swapPointers();
#pragma omp barrier
}
void CASolver::fillBoundaryConditions()
{
    fillHorizontalBoundaryConditions();
#pragma omp barrier
    fillVerticalBoundaryConditions();
}

void CASolver::printCurrentState()
{
    copyCurrentStateToHostVisibleData();
    // CADataPrinter::printCADataWithHalo(hostVisibleData);
    CADataPrinter::printCADataWithoutHalo(hostVisibleData);
}
