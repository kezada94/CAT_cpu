#include "CellularAutomata/CADataPrinter.h"
#include "CellularAutomata/Solvers/AVXSolver512.h"
#include "Memory/Allocators/CPUAllocator.h"

AVXSolver512::AVXSolver512(CADataDomain<int> *domain, CADataDomain<int> *domainBuffer)
{
    dataDomain = domain;
    dataDomainBuffer = domainBuffer;

    CPUAllocator<int> *cpuAllocator = new CPUAllocator<int>();
    Allocator<int> *allocator = reinterpret_cast<Allocator<int> *>(cpuAllocator);
    hostVisibleData = new CADataDomain<int>(allocator, dataDomain->getInnerHorizontalSize(), dataDomain->getHorizontalHaloSize());
    hostVisibleData->allocate();
}

void AVXSolver512::copyCurrentStateToHostVisibleData()
{
    for (size_t i = 0; i < dataDomain->getTotalSize(); ++i)
    {
        int value = dataDomain->getElementAt(i);
        hostVisibleData->setElementAt(i, (int)value);
    }
}
void AVXSolver512::copyHostVisibleDataToCurrentState()
{
    for (size_t i = 0; i < hostVisibleData->getTotalSize(); ++i)
    {
        int value = hostVisibleData->getElementAt(i);
        dataDomain->setElementAt(i, value);
    }
}

void AVXSolver512::swapPointers()
{
    if (omp_get_thread_num() == 0)
    {
        CADataDomain<int> *temp = dataDomain;
        dataDomain = dataDomainBuffer;
        dataDomainBuffer = temp;
    }
}

int AVXSolver512::transitionFunction(int k, int a, int b)
{
    return (1 - (((k - a) >> 31) & 0x1)) * (1 - (((b - k) >> 31) & 0x1));
}

void AVXSolver512::preamble()
{
}

void AVXSolver512::CAStepAlgorithm()
{
    int* data = dataDomain->getData();
    size_t halo = dataDomain->getHorizontalHaloSize();
#pragma omp for
    for (int i = 0; i < dataDomain->getInnerHorizontalSize(); ++i)
    {
        for (int j = 0; j < dataDomain->getInnerHorizontalSize(); j += 8)
        {
            #ifdef USE_AVX512

            // Load 8 elements from the center position
            __m512i center = _mm512_loadu_si512((__m512i *)&data[(i+halo) * dataDomain->getFullHorizontalSize() + j + halo]);
            __m512i neighbor_sum = _mm512_setzero_si512(); // Initialize sum for 32-bit integers

            for (int y = -RADIUS; y <= RADIUS; ++y)
            {
                for (int x = -RADIUS; x <= RADIUS; ++x)
                {
                    if (x == 0 && y == 0) // Skip the center element itself
                        continue;

                    // Calculate neighbor index
                    size_t index = (y + i + halo) * dataDomain->getFullHorizontalSize() + (x + j + halo);

                    // Load neighbor elements (assume these are also 32-bit integers)
                    __m512i neighbor = _mm512_loadu_si512((__m512i *)(dataDomain->getData() + index));

                    // Add the neighbor values
                    neighbor_sum = _mm512_add_epi32(neighbor_sum, neighbor);
                }
            }


            __m512i live_mask = center; // Mask for live cells (center > 0)

            __m512i survival_mask = _mm512_and_si512(
                _mm512_cmpgt_epi32(neighbor_sum, _mm512_set1_epi32(SMIN - 1)), // Neighbors >= SMIN
                _mm512_cmpgt_epi32(_mm512_set1_epi32(SMAX + 1), neighbor_sum)  // Neighbors <= SMAX
            );

            // Condition for birth: dead cells with BMIN <= neighbors <= BMAX become alive
            __m512i dead_mask = _mm512_cmpeq_epi32(center, _mm512_setzero_si512()); // Mask for dead cells

            __m512i birth_mask = _mm512_and_si512(
                _mm512_cmpgt_epi32(neighbor_sum, _mm512_set1_epi32(BMIN - 1)), // Neighbors >= BMIN
                _mm512_cmpgt_epi32(_mm512_set1_epi32(BMAX + 1), neighbor_sum)  // Neighbors <= BMAX
            );

            // New state: live cells that survive OR dead cells that are born
            __m512i new_state = _mm512_or_si512(
                _mm512_and_si512(live_mask, survival_mask), // Surviving live cells
                _mm512_and_si512(_mm512_and_si512(dead_mask, _mm512_set1_epi32(1)) , birth_mask)     // Newly born cells
            );

            // Store the new state back to the buffer
            _mm512_storeu_si512((__m512i *)((dataDomainBuffer->getData() + (i+halo) * dataDomainBuffer->getFullHorizontalSize() + j+halo)), new_state);
        #endif
        }
    }
}
int AVXSolver512::countAliveNeighbors(int y, int x)
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

void AVXSolver512::fillHorizontalBoundaryConditions()
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

void AVXSolver512::fillVerticalBoundaryConditions()
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
