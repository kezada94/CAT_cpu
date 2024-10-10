#include "CellularAutomata/CADataPrinter.h"
#include "CellularAutomata/Solvers/AVXSolver.h"
#include "Memory/Allocators/CPUAllocator.h"

AVXSolver::AVXSolver(CADataDomain<int> *domain, CADataDomain<int> *domainBuffer)
{
    dataDomain = domain;
    dataDomainBuffer = domainBuffer;

    CPUAllocator<int> *cpuAllocator = new CPUAllocator<int>();
    Allocator<int> *allocator = reinterpret_cast<Allocator<int> *>(cpuAllocator);
    hostVisibleData = new CADataDomain<int>(allocator, dataDomain->getInnerHorizontalSize(), dataDomain->getHorizontalHaloSize());
    hostVisibleData->allocate();
}

void AVXSolver::copyCurrentStateToHostVisibleData()
{
    for (size_t i = 0; i < dataDomain->getTotalSize(); ++i)
    {
        int value = dataDomain->getElementAt(i);
        hostVisibleData->setElementAt(i, (int)value);
    }
}
void AVXSolver::copyHostVisibleDataToCurrentState()
{
    for (size_t i = 0; i < hostVisibleData->getTotalSize(); ++i)
    {
        int value = hostVisibleData->getElementAt(i);
        dataDomain->setElementAt(i, value);
    }
}

void AVXSolver::swapPointers()
{
    if (omp_get_thread_num() == 0)
    {
        CADataDomain<int> *temp = dataDomain;
        dataDomain = dataDomainBuffer;
        dataDomainBuffer = temp;
    }
}

int AVXSolver::transitionFunction(int k, int a, int b)
{
    return (1 - (((k - a) >> 31) & 0x1)) * (1 - (((b - k) >> 31) & 0x1));
}

void AVXSolver::preamble()
{
}

void AVXSolver::CAStepAlgorithm()
{
    alignas(32) int aliveNeighbors[8];
#pragma omp for
    for (int i = 0; i < dataDomain->getInnerHorizontalSize(); ++i)
    {
        for (int j = 0; j < dataDomain->getInnerHorizontalSize(); j += 8)
        {
            // printf("i: %d, j: %d\n", i, j);
            // Load 8 elements from the center position
            __m256i center = _mm256_loadu_si256((__m256i *)(dataDomain->getData() + i * dataDomain->getFullHorizontalSize() + j));
            __m256i neighbor_sum = _mm256_setzero_si256(); // Initialize sum for 32-bit integers

            for (int y = -RADIUS; y <= RADIUS; ++y)
            {
                for (int x = -RADIUS; x <= RADIUS; ++x)
                {
                    if (x == 0 && y == 0) // Skip the center element itself
                        continue;

                    // Calculate neighbor index
                    size_t index = (y + i) * dataDomain->getFullHorizontalSize() + (x + j);
                    // Load neighbor elements (assume these are also 32-bit integers)
                    __m256i neighbor = _mm256_loadu_si256((__m256i *)(dataDomain->getData() + index));

                    // Add the neighbor values
                    neighbor_sum = _mm256_add_epi32(neighbor_sum, neighbor);
                }
            }
            // Apply Game of Life rules
            // Condition for survival: live cells with 2 or 3 neighbors survive, others die
            __m256i live_mask = _mm256_cmpgt_epi32(center, _mm256_setzero_si256()); // Mask for live cells (center > 0)
            __m256i survival_mask = _mm256_and_si256(
                _mm256_cmpgt_epi32(neighbor_sum, _mm256_set1_epi32(SMIN - 1)), // Neighbors >= SMIN
                _mm256_cmpgt_epi32(_mm256_set1_epi32(SMAX + 1), neighbor_sum)  // Neighbors <= SMAX
            );

            // Condition for birth: dead cells with exactly 3 neighbors become alive
            __m256i birth_mask = _mm256_cmpeq_epi32(neighbor_sum, _mm256_set1_epi32(BMIN));

            // New state: live cells that survive OR dead cells that are born
            __m256i new_state = _mm256_or_si256(
                _mm256_and_si256(live_mask, survival_mask), // Surviving live cells
                birth_mask                                  // Newly born cells
            );

            // Store the new state back to the buffer (could be packed into 8-bit if necessary)
            _mm256_storeu_si256((__m256i *)((dataDomainBuffer->getData() + i * dataDomainBuffer->getFullHorizontalSize() + j)), new_state);
        }
    }
}

int AVXSolver::countAliveNeighbors(int y, int x)
{
    int aliveNeighbors = 0;

    for (int i = -RADIUS; i <= RADIUS; ++i)
    {
        for (int j = -RADIUS; j <= RADIUS; ++j)
        {
            if (i == 0 && j == 0)
                continue;
            aliveNeighbors += dataDomain->getInnerElementAt(y + i, x + j);
        }
    }

    return aliveNeighbors;
}

void AVXSolver::fillHorizontalBoundaryConditions()
{
#pragma omp for
    for (int h = 0; h < dataDomain->getHorizontalHaloSize(); ++h)
    {
        for (int j = 0; j < dataDomain->getInnerHorizontalSize(); ++j)
        {
            size_t topIndex = (dataDomain->getHorizontalHaloSize() + h) * dataDomain->getFullHorizontalSize() + dataDomain->getHorizontalHaloSize() + j;
            size_t bottomIndex = topIndex + (dataDomain->getInnerHorizontalSize()) * dataDomain->getFullHorizontalSize();
            int value = dataDomain->getElementAt(topIndex);
            dataDomain->setElementAt(bottomIndex, value);
        }

        for (int j = 0; j < dataDomain->getInnerHorizontalSize(); ++j)
        {
            size_t topIndex = (h)*dataDomain->getFullHorizontalSize() + dataDomain->getHorizontalHaloSize() + j;
            size_t bottomIndex = topIndex + (dataDomain->getInnerHorizontalSize()) * dataDomain->getFullHorizontalSize();

            int value = dataDomain->getElementAt(bottomIndex);
            dataDomain->setElementAt(topIndex, value);
        }
    }
}

void AVXSolver::fillVerticalBoundaryConditions()
{
#pragma omp for
    for (int h = 0; h < dataDomain->getHorizontalHaloSize(); ++h)
    {
        for (int i = 0; i < dataDomain->getFullHorizontalSize(); ++i)
        {
            size_t leftIndex = i * dataDomain->getFullHorizontalSize() + h;
            size_t rightIndex = leftIndex + dataDomain->getInnerHorizontalSize();
            int value = dataDomain->getElementAt(rightIndex);
            dataDomain->setElementAt(leftIndex, value);
        }

        for (int i = 0; i < dataDomain->getFullHorizontalSize(); ++i)
        {
            size_t leftIndex = i * dataDomain->getFullHorizontalSize() + dataDomain->getHorizontalHaloSize() + h;
            size_t rightIndex = leftIndex + dataDomain->getInnerHorizontalSize();
            int value = dataDomain->getElementAt(leftIndex);
            dataDomain->setElementAt(rightIndex, value);
        }
    }
}
