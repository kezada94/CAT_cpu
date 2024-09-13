#pragma once

#include "CellularAutomata/Solvers/HostSolver.h"
#include "CellularAutomata/Solvers/AMXSolver.h"
#include "Debug.h"
#include "Defines.h"
#include "Memory/Allocators/CPUAllocator.h"

class CASolverFactory
{
  public:
    static CASolver *createSolver(int SOLVER_CODE, int deviceId, int sideLength, int haloWidth);
};