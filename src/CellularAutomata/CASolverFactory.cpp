#include "CellularAutomata/CASolverFactory.h"
#include "Memory/CADataDomain.h"

static const int TENSOR_HALO_SIZE = 16;

CASolver *CASolverFactory::createSolver(int SOLVER_CODE, int deviceId, int fullHorizontalSize, int horizontalHaloSize, int nThreads)
{
    CASolver *solver = nullptr;

    switch (SOLVER_CODE)
    {
    case 0:
    {
        CPUAllocator<int> *cpuAllocator = new CPUAllocator<int>();
        Allocator<int> *allocator = reinterpret_cast<Allocator<int> *>(cpuAllocator);
        CADataDomain<int> *dataDomain = new CADataDomain<int>(allocator, fullHorizontalSize, horizontalHaloSize);
        dataDomain->allocate();

        CADataDomain<int> *dataDomainBuffer = new CADataDomain<int>(allocator, fullHorizontalSize, horizontalHaloSize);
        dataDomainBuffer->allocate();

        solver = new HostSolver<int>(dataDomain, dataDomainBuffer);
        lDebug(1, "Solver of type HostSolver created");
        break;
    }
    case 1:
    {
        lDebug(1, "Creating solver of type AMXSolver");
        CPUAllocator<uint8_t> *cpuAllocator = new CPUAllocator<uint8_t>();
        Allocator<uint8_t> *allocator = reinterpret_cast<Allocator<uint8_t> *>(cpuAllocator);
        CADataDomain<uint8_t> *dataDomain = new CADataDomain<uint8_t>(allocator, fullHorizontalSize, 64);
        dataDomain->allocate();

        CADataDomain<uint8_t> *dataDomainBuffer = new CADataDomain<uint8_t>(allocator, fullHorizontalSize, 64);
        dataDomainBuffer->allocate();

        solver = new AMXSolver(dataDomain, dataDomainBuffer, nThreads);
        lDebug(1, "Solver of type AMXSolver created");

        break;
    }
    case 2:
    {
        lDebug(1, "Creating solver of type AMXSolver16");
        CPUAllocator<uint8_t> *cpuAllocator = new CPUAllocator<uint8_t>();
        Allocator<uint8_t> *allocator = reinterpret_cast<Allocator<uint8_t> *>(cpuAllocator);
        CADataDomain<uint8_t> *dataDomain = new CADataDomain<uint8_t>(allocator, fullHorizontalSize, 16);
        dataDomain->allocate();

        CADataDomain<uint8_t> *dataDomainBuffer = new CADataDomain<uint8_t>(allocator, fullHorizontalSize, 16);
        dataDomainBuffer->allocate();

        solver = new AMXSolver16(dataDomain, dataDomainBuffer, nThreads);
        lDebug(1, "Solver of type AMXSolver created");

        break;
    }
    case 3:
    {
        CPUAllocator<int> *cpuAllocator = new CPUAllocator<int>();
        Allocator<int> *allocator = reinterpret_cast<Allocator<int> *>(cpuAllocator);
        CADataDomain<int> *dataDomain = new CADataDomain<int>(allocator, fullHorizontalSize, horizontalHaloSize);
        dataDomain->allocate();

        CADataDomain<int> *dataDomainBuffer = new CADataDomain<int>(allocator, fullHorizontalSize, horizontalHaloSize);
        dataDomainBuffer->allocate();

        solver = new AVXSolver(dataDomain, dataDomainBuffer);
        lDebug(1, "Solver of type AVXSOLVER created");
        break;
    }
    }
    return solver;
}
