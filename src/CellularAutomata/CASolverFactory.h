#pragma once

#include "CellularAutomata/Solvers/HostSolver.h"

#include "CellularAutomata/Solvers/AMXSolver.h"
#include "CellularAutomata/Solvers/AMXSolver512.h"

#include "CellularAutomata/Solvers/AMXSolver16.h"
#include "CellularAutomata/Solvers/AMXSolver16x512.h"

#include "CellularAutomata/Solvers/AVXSolver.h"
#include "CellularAutomata/Solvers/AVXSolver512.h"

#include "Debug.h"
#include "Defines.h"
#include "Memory/Allocators/CPUAllocator.h"

class CASolverFactory
{
public:
  static CASolver *createSolver(int SOLVER_CODE, int deviceId, int sideLength, int haloWidth, int nThreads);
};
