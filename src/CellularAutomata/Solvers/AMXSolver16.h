
#pragma once

#include "CellularAutomata/Solvers/AMXSolver.h"
#include "CellularAutomata/Solvers/CASolver.h"
#include "Memory/CADataDomain.h"
#include "Memory/CAStateGenerator.h"

#include <immintrin.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <stdbool.h>

#define ARCH_GET_XCOMP_PERM 0x1022
#define ARCH_REQ_XCOMP_PERM 0x1023
#define XFEATURE_XTILECFG 17
#define XFEATURE_XTILEDATA 18

class AMXSolver16 : public CASolver
{
private:
    uint8_t *pi_1;
    uint8_t *pi_2;
    uint8_t *pi_3;
    uint8_t *pi_1B;
    uint8_t *pi_2B;
    uint8_t *pi_3B;
    int *buffer;

    __tilecfg *tile_config;

    CADataDomain<uint8_t> *dataDomain;
    CADataDomain<uint8_t> *dataDomainBuffer;
    uint8_t **dataI;

    void preamble() override;

    void CAStepAlgorithm() override;
    void fillHorizontalBoundaryConditions() override;
    void fillVerticalBoundaryConditions() override;

    void swapPointers() override;
    void copyCurrentStateToHostVisibleData() override;
    void copyHostVisibleDataToCurrentState() override;

    uint8_t transitionFunction(int k, int a, int b);
    int countAliveNeighbors(int i, int j);

    void fillTridiag();
    void setupAMX();

public:
    AMXSolver16(CADataDomain<uint8_t> *domain, CADataDomain<uint8_t> *domainBuffer, int nThreads);
};
