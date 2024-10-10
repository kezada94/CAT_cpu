
#pragma once

#include "CellularAutomata/Solvers/CASolver.h"
#include "Memory/CADataDomain.h"
#include "Memory/CAStateGenerator.h"
#include <immintrin.h>
#include <stdint.h>
#include <stddef.h>

class AVXSolver : public CASolver
{
private:
    CADataDomain<int> *dataDomain;
    CADataDomain<int> *dataDomainBuffer;

    void preamble() override;

    void CAStepAlgorithm() override;
    void fillHorizontalBoundaryConditions() override;
    void fillVerticalBoundaryConditions() override;

    void swapPointers() override;
    void copyCurrentStateToHostVisibleData() override;
    void copyHostVisibleDataToCurrentState() override;

    int transitionFunction(int k, int a, int b);
    int countAliveNeighbors(int i, int j);

public:
    AVXSolver(CADataDomain<int> *domain, CADataDomain<int> *domainBuffer);
};
