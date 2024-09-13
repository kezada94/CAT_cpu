#pragma once

#include "CellularAutomata/Solvers/CASolver.h"
#include "Memory/CADataDomain.h"

class CADataDomainComparator {
   private:
    CASolver* solver;
    CASolver* referenceSolver;

    bool areDifferentSize();

   public:
    CADataDomainComparator(CASolver* pSolver, CASolver* pReferenceSolver);

    bool compareCurrentStates();
};

#include "CellularAutomata/CADataDomainComparator.tpp"