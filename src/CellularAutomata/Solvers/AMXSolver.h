
#pragma once

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


//Define tile config data structure
typedef struct __tile_config
{
  uint8_t palette_id;
  uint8_t start_row;
  uint8_t reserved_0[14];
  uint16_t colsb[16];
  uint8_t rows[16];
} __tilecfg;

class AMXSolver : public CASolver {
   private:

    __tilecfg *tile_config;

    uint8_t pi_1[16*64];
    uint8_t pi_2[16*64];
    uint8_t pi_3[16*64];


    CADataDomain<uint8_t>* dataDomain;
    CADataDomain<uint8_t>* dataDomainBuffer;

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
    AMXSolver(CADataDomain<uint8_t>* domain, CADataDomain<uint8_t>* domainBuffer);
};
