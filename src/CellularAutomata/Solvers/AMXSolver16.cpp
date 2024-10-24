#include "CellularAutomata/CADataPrinter.h"
#include "CellularAutomata/Solvers/AMXSolver16.h"
#include "Memory/Allocators/CPUAllocator.h"

AMXSolver16::AMXSolver16(CADataDomain<uint8_t>* domain, CADataDomain<uint8_t>* domainBuffer, int nThreads) {
    dataDomain = domain;
    dataDomainBuffer = domainBuffer;

    CPUAllocator<int>* cpuAllocator = new CPUAllocator<int>();
    Allocator<int>* allocator = reinterpret_cast<Allocator<int>*>(cpuAllocator);
    hostVisibleData = new CADataDomain<int>(allocator, dataDomain->getInnerHorizontalSize(), dataDomain->getHorizontalHaloSize());
    hostVisibleData->allocate();

        CPUAllocator<uint8_t> *cpuAllocatorint8 = new CPUAllocator<uint8_t>();
        Allocator<uint8_t> *allocatorint8 = reinterpret_cast<Allocator<uint8_t> *>(cpuAllocatorint8);
        dataDomainIntermediate = new CADataDomain<uint8_t>(allocatorint8, dataDomain->getInnerHorizontalSize(), dataDomain->getHorizontalHaloSize());
        dataDomainIntermediate->allocate();

	buffer = new int[16*16*nThreads];

	lDebug(1, "Setting up AMX");
    setupAMX();
    
	lDebug(1, "Filing Tridiag matrices");
    fillTridiag();
}

void AMXSolver16::setupAMX() {
	tile_config = new __tile_config();
    tile_config->palette_id = 1;
    tile_config->start_row = 0;

    // Configure tiles for block_size x block_size matrices
    for (int i = 0; i < 1; ++i) {
        tile_config->colsb[i] = 16 * 4;
        tile_config->rows[i] = 16;
    }

    for (int i = 1; i < 8; ++i) {
		//if (i==7){
        //tile_config->colsb[i] = 64;
        //tile_config->rows[i] = 16;
		//} else {
        tile_config->colsb[i] = 64;
        tile_config->rows[i] = 16;
		//}
    }

    _tile_loadconfig(tile_config);

    if (syscall(SYS_arch_prctl, ARCH_REQ_XCOMP_PERM, XFEATURE_XTILEDATA))
	{
	   lDebug(1,"\n Failed to enable XFEATURE_XTILEDATA \n\n");
	}
	else
	{
	   lDebug(1, "\n TILE DATA USE SET - OK \n\n");
	}

}


void AMXSolver16::fillTridiag() {
	uint8_t* data = (uint8_t*) malloc(64*16*6*sizeof(uint8_t));
	pi_1  = &data[64*16*0];
	pi_2  = &data[64*16*1];
	pi_3  = &data[64*16*2];
	pi_1B  = &data[64*16*3];
	pi_2B  = &data[64*16*4];
	pi_3B  = &data[64*16*5];
    int i;

    for (i = 0; i < 16*64; i += 1)
    {
        int col = i & 63;
        int row = i >> 6;
        if (col + 15 - RADIUS < row ){
            pi_1[i] = 1;
        } else {
            pi_1[i] = 0;
        }

        if (abs(col - row) <= RADIUS){
        //if (abs(col - row) <= 0){
        //if (col == row && row == 1){
            pi_2[i] = 1;
        } else {
            pi_2[i] = 0;
        }

        if (col -15+RADIUS> row){
            pi_3[i] = 1;
        } else {
            pi_3[i] = 0;
        }
		if (col>=16){
		
            pi_1[i] = 0;
            pi_2[i] = 0;
            pi_3[i] = 0;
		}
    }
   // for (int i = 0; i < 16; i++) {
   //     for (int j = 0; j < 64; j++) {
   //         std::cout << (int)pi_1[i * 64 + j] << " ";
   //     }
   //     std::cout << std::endl;
   // }
   // std::cout << std::endl;
   // for (int i = 0; i < 16; i++) {
   //     for (int j = 0; j < 64; j++) {
   //         std::cout << (int)pi_2[i * 64 + j] << " ";
   //     }
   //     std::cout << std::endl;
   // }
   // std::cout << std::endl;
   // for (int i = 0; i < 16; i++) {
   //     for (int j = 0; j < 64; j++) {
   //         std::cout << (int)pi_3[i * 64 + j] << " ";
   //     }
   //     std::cout << std::endl;
   // }
   // std::cout << std::endl;

	uint32_t *pi1_int = (uint32_t*)pi_1; 
	uint32_t *pi1B_int = (uint32_t*)pi_1B; 
	uint32_t *pi2_int = (uint32_t*)pi_2; 
	uint32_t *pi2B_int = (uint32_t*)pi_2B; 
	uint32_t *pi3_int = (uint32_t*)pi_3; 
	uint32_t *pi3B_int = (uint32_t*)pi_3B; 
    // // debug print tridiag in 2d
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
			pi1B_int[j*16+i] = pi1_int[i*16+j];
			pi2B_int[j*16+i] = pi2_int[i*16+j];
			pi3B_int[j*16+i] = pi3_int[i*16+j];
        }
    }
	

    // // debug print tridiag in 2d
    // for (int i = 0; i < 16; i++) {
    //     for (int j = 0; j < 64; j++) {
    //         std::cout << (int)pi_2[i * 64 + j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // for (int i = 0; i < 64; i++) {
    //     for (int j = 0; j < 16; j++) {
    //         std::cout << (int)pi_2B[i * 16 + j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
	// std::cout << std::endl;

    // std::cout << std::endl;

    // for (int i = 0; i < 16; i++) {
    //     for (int j = 0; j < 64; j++) {
    //         std::cout << (int)pi_2B[i * 64 + j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
	// std::cout << std::endl;

    // for (int i = 0; i < 16; i++) {
    //     for (int j = 0; j < 64; j++) {
    //         std::cout << (int)pi_3[i * 64 + j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    //     std::cout << std::endl;
    // for (int i = 0; i < 16; i++) {
    //     for (int j = 0; j < 64; j++) {
    //         std::cout << (int)pi_4[i * 64 + j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    //     std::cout << std::endl;
    // for (int i = 0; i < 16; i++) {
    //     for (int j = 0; j < 64; j++) {
    //         std::cout << (int)pi_5[i * 64 + j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    //     std::cout << std::endl;
    // for (int i = 0; i < 16; i++) {
    //     for (int j = 0; j < 64; j++) {
    //         std::cout << (int)pi_6[i * 64 + j] << " ";
    //     }
    //     std::cout << std::endl;
    // }


}

void AMXSolver16::copyCurrentStateToHostVisibleData() {
	lDebug(1, "%llu\n", dataDomain->getTotalSize());
    for (int i = 0; i < dataDomain->getTotalSize(); ++i) {
        uint8_t value = dataDomain->getElementAt(i);
        hostVisibleData->setElementAt(i, (int)value);
    }
}
void AMXSolver16::copyHostVisibleDataToCurrentState() {
    for (int i = 0; i < hostVisibleData->getTotalSize(); ++i) {
        int value = hostVisibleData->getElementAt(i);
        dataDomain->setElementAt(i, value);
    }
}

void AMXSolver16::swapPointers() {
	if (omp_get_thread_num() == 0){
		CADataDomain<uint8_t>* temp = dataDomain;
		dataDomain = dataDomainBuffer;
		dataDomainBuffer = temp;
	}
}

uint8_t AMXSolver16::transitionFunction(int k, int a, int b) {
    return (1 - (((k - a) >> 31) & 0x1)) * (1 - (((b - k) >> 31) & 0x1));
}
void AMXSolver16::preamble()
{
    _tile_loadd(1, pi_1B, 64);
    _tile_loadd(2, pi_2B, 64);
    _tile_loadd(3, pi_3B, 64);
}

void AMXSolver16::CAStepAlgorithm() {

    uint8_t* data = dataDomain->getData();
    uint8_t* dataBuffer = dataDomainBuffer->getData();
    uint8_t* dataI = dataDomainIntermediate->getData();
    size_t nWithHalo = dataDomain->getFullHorizontalSize();

	int num_threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();

	int iterations = nWithHalo/16;
	int chunk_size = (iterations + num_threads - 1) / num_threads;
	//printf("cs %i, nwh %llu, iters: %i\n", chunk_size, nWithHalo, iterations);
	int start = thread_id * chunk_size;
	int end = std::min(start + chunk_size, iterations);  // Ensure 'end' does not exceed 'iterations'


	if ((thread_id) >= iterations)
	{
		end = -1;
	}

    //FIRST STEP: horizontal reduction
    for (size_t iter = start; iter < end; iter++) {
        for (size_t j = 0; j < nWithHalo - 16*2; j+=16) {
            //take three continuous 16x16 blocks and load them into amx
			//printf("%i,%i\n", i, j);
			size_t i = iter*16;

		    // tiles in C from 0,0 to 3,0 (first col)
            _tile_zero(0);

			_tile_loadd(7, &data[(i) * nWithHalo + j], nWithHalo);
			_tile_dpbssd(0, 7, 3);

			_tile_loadd(7, &data[(i) * nWithHalo + j + 16], nWithHalo);
			_tile_dpbuud(0, 7, 2);
			
			_tile_loadd(7, &data[(i) * nWithHalo + j + 32], nWithHalo);
			_tile_dpbssd(0, 7, 1);


			_tile_stored(0, &buffer[thread_id*16*16], 16*4);

		//	printf("BEfORE TRANSPOSE\n");
	//		for (int ii=0; ii<64; ii++){
	//			for (int jj=0; jj<64; jj+=4){
	//				for (int kk=0; kk<4; kk++){
	//					printf("%i ", buffer[ii*64 + jj + kk] );
	//				}
	//			}
	//			printf("\n");
	//		}
	//		printf("\n");

		//	for (int slice =0; slice < 4; slice++){
		//		for (int ii=0; ii<16; ii++){
		//			for (int jj=0; jj<16; jj++){
		//				for (int kk=0; kk<4; kk++){
		//					dataBuffer[slice*nWithHalo*16 +(i+jj)*nWithHalo + j + 64 + ii*4 + kk] = buffer[slice*64*16 +  (ii)*64 + (jj)*4 + kk ];
		//				}
		//			}
		//		}
		//	}
			for (int ii=0; ii<16; ii++){
				for (int jj=0; jj<16; jj++){
						dataI[(i+jj)*nWithHalo + j + ii + 16] = buffer[thread_id*16*16 + ii*16+jj];
				}
			}
		//	printf("AFTER TRANSPOSE\n");
		//	for (int ii=0; ii<64; ii++){
		//		for (int jj=0; jj<64; jj++){
		//				printf("%d ", dataBuffer[(i+ii)*nWithHalo + j + 64 + jj]);
		//		}
		//		printf("\n");
		//	}
		//	printf("\n");
        }
    }
#pragma omp barrier
	iterations = (nWithHalo - 16*2)/16;
	chunk_size = (iterations +num_threads-1)/ num_threads;
	start = thread_id * chunk_size;
	end = std::min(start + chunk_size, iterations);  // Ensure 'end' does not exceed 'iterations'
    //for (int i = 0; i < nWithHalo - 64*2; i+=64) {
	if ((thread_id) >= iterations)
	{
		// printf("Thread %i is out of bounds\n", thread_id);
		end = -1;
	}

    for (int iter = start; iter < end; iter++) {
        for (int j = 0; j < nWithHalo - 16*2; j+=16) {
			size_t i = iter*16;
            //take three continuous 16x16 blocks and load them into amx
            _tile_zero(0);

			_tile_loadd(7, &dataI[(i) * nWithHalo + j + 16], nWithHalo);
			_tile_dpbuud(0, 7, 3);

			_tile_loadd(7, &dataI[(i+16) * nWithHalo + j + 16], nWithHalo);
			_tile_dpbuud(0, 7, 2);

			_tile_loadd(7, &dataI[(i+32) * nWithHalo + j + 16], nWithHalo);
			_tile_dpbssd(0, 7, 1);

			_tile_stored(0, &buffer[thread_id*16*16], 16*4);


			for (int ii=0; ii<16; ii++){
				for (int jj=0; jj<16; jj+=8){
                    int cente[8];
					int neigh[8];
					int res[8];
					for (int kk=0; kk<8; kk++){
						cente[kk] = data[(i+ii + 16)*nWithHalo + j + 16 + jj + kk];
						neigh[kk] = buffer[thread_id*16*16 + (jj+kk)*16+ii] - cente[kk];
					}
        
					__m256i center = _mm256_loadu_si256((__m256i *)&cente[0]);
					__m256i neighbor_sum = _mm256_loadu_si256((__m256i *)&neigh[0]);

					__m256i live_mask = center; // Mask for live cells (center > 0)

					__m256i survival_mask = _mm256_and_si256(
						_mm256_cmpgt_epi32(neighbor_sum, _mm256_set1_epi32(SMIN - 1)), // Neighbors >= SMIN
						_mm256_cmpgt_epi32(_mm256_set1_epi32(SMAX + 1), neighbor_sum)  // Neighbors <= SMAX
					);

					// Condition for birth: dead cells with BMIN <= neighbors <= BMAX become alive
					__m256i dead_mask = _mm256_cmpeq_epi32(center, _mm256_setzero_si256()); // Mask for dead cells

					__m256i birth_mask = _mm256_and_si256(
						_mm256_cmpgt_epi32(neighbor_sum, _mm256_set1_epi32(BMIN - 1)), // Neighbors >= BMIN
						_mm256_cmpgt_epi32(_mm256_set1_epi32(BMAX + 1), neighbor_sum)  // Neighbors <= BMAX
					);

					// New state: live cells that survive OR dead cells that are born
					__m256i new_state = _mm256_or_si256(
						_mm256_and_si256(live_mask, survival_mask), // Surviving live cells
						_mm256_and_si256(_mm256_and_si256(dead_mask, _mm256_set1_epi32(1)) , birth_mask)     // Newly born cells
					);

					_mm256_storeu_si256((__m256i *)&res[0], new_state);
					for (int kk=0; kk<8; kk++){
						dataBuffer[(i+ii + 16)*nWithHalo + j + 16 + jj + kk] = res[kk];
					}
				}
				//printf("\n");
			}
			//printf("\n");

        }
    }

    //printf("DFONE\n");

    // for (int i = 0; i < dataDomain->getInnerHorizontalSize(); ++i) {
    //     for (int j = 0; j < dataDomain->getInnerHorizontalSize(); ++j) {
    //         int liveNeighbors = countAliveNeighbors(i, j);
    //         uint8_t cellValue = dataDomain->getInnerElementAt(i, j);
    //         uint8_t result = cellValue * transitionFunction(liveNeighbors, SMIN, SMAX) + (1 - cellValue) * transitionFunction(liveNeighbors, BMIN, BMAX);

    //         dataDomainBuffer->setInnerElementAt(i, j, result);
    //     }
    // }
}

int AMXSolver16::countAliveNeighbors(int y, int x) {
    int aliveNeighbors = 0;

    for (int i = -RADIUS; i <= RADIUS; ++i) {
        for (int j = -RADIUS; j <= RADIUS; ++j) {
            if (i == 0 && j == 0)
                continue;
            aliveNeighbors += dataDomain->getInnerElementAt(y + i, x + j);
        }
    }

    return aliveNeighbors;
}

void AMXSolver16::fillHorizontalBoundaryConditions() {

    //for (int h = 0; h < dataDomain->getHorizontalHaloSize(); ++h) {
	int num_threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();

	int innerSize = dataDomain->getInnerHorizontalSize();
	int chunk_size = innerSize / num_threads;
	int start = thread_id * chunk_size;
	int end = (thread_id == num_threads - 1) ? innerSize : start + chunk_size;


    for (int h = 0; h < RADIUS; ++h) {
		// Top elements to bottom halo
        //for (int j = 0; j < dataDomain->getInnerHorizontalSize(); ++j) {
        for (int j = start; j < end; ++j) {
            size_t topIndex = (dataDomain->getHorizontalHaloSize() + h) * dataDomain->getFullHorizontalSize() + dataDomain->getHorizontalHaloSize() + j;
            size_t bottomIndex = topIndex + (dataDomain->getInnerHorizontalSize()) * dataDomain->getFullHorizontalSize();
            uint8_t value = dataDomain->getElementAt(topIndex);
            dataDomain->setElementAt(bottomIndex, value);
        }

		// Bottom elements to top halo
        for (int j = start; j < end; ++j) {
            size_t topIndex = (dataDomain->getHorizontalHaloSize() - h - 1)*dataDomain->getFullHorizontalSize() + dataDomain->getHorizontalHaloSize() + j;
            size_t bottomIndex = topIndex + (dataDomain->getInnerHorizontalSize()) * dataDomain->getFullHorizontalSize();

            uint8_t value = dataDomain->getElementAt(bottomIndex);
            dataDomain->setElementAt(topIndex, value);
        }
    }
}

void AMXSolver16::fillVerticalBoundaryConditions() {

	int num_threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();

	int fullSize = dataDomain->getFullHorizontalSize();
    int chunk_size = fullSize / num_threads;
    int start = thread_id * chunk_size;
    int end = (thread_id == num_threads - 1) ? fullSize : start + chunk_size;


    //for (int h = 0; h < dataDomain->getHorizontalHaloSize(); ++h) {
    for (int h = 0; h < RADIUS; ++h) {
		//Rightmost elements to left halo
        //for (int i = 0; i < dataDomain->getFullHorizontalSize(); ++i) {
        for (int i = start; i < end; ++i) {
            size_t leftIndex = i * dataDomain->getFullHorizontalSize() + (dataDomain->getHorizontalHaloSize() - h - 1);
            size_t rightIndex = leftIndex + dataDomain->getInnerHorizontalSize();
            uint8_t value = dataDomain->getElementAt(rightIndex);
            dataDomain->setElementAt(leftIndex, value);
        }

        for (int i = start; i < end; ++i) {
            size_t leftIndex = i * dataDomain->getFullHorizontalSize() + dataDomain->getHorizontalHaloSize() + h;
            size_t rightIndex = leftIndex + dataDomain->getInnerHorizontalSize();
            uint8_t value = dataDomain->getElementAt(leftIndex);
            dataDomain->setElementAt(rightIndex, value);
        }
    }
}
