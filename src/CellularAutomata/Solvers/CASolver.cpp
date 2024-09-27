#include "CellularAutomata/Solvers/CASolver.h"

void CASolver::resetState(int seed, float density) {
    CAStateGenerator::generateRandomState(hostVisibleData, seed, density);
    copyHostVisibleDataToCurrentState();
}

CADataDomain<int>* CASolver::getCurrentState() {
    copyCurrentStateToHostVisibleData();
    return hostVisibleData;
}

void CASolver::doSteps(int stepNumber) {
    for (int i = 0; i < stepNumber; i++) {
        doStep();
    }
}
void CASolver::doStep() {
	lDebug(1, "	STEP: begin");
	lDebug(1, "	STEP: filling boundary conditions");
    fillBoundaryConditions();
	lDebug(1, "	STEP: Stepping");
    CAStepAlgorithm();
	lDebug(1, "	STEP: swapping pointers");
    swapPointers();
	lDebug(1, "	STEP: end");
}
void CASolver::fillBoundaryConditions() {
    fillHorizontalBoundaryConditions();
    fillVerticalBoundaryConditions();
}

void CASolver::printCurrentState() {
    copyCurrentStateToHostVisibleData();
    //CADataPrinter::printCADataWithHalo(hostVisibleData);
    CADataPrinter::printCADataWithoutHalo(hostVisibleData);
}
