#include "CellularAutomata/CADataPrinter.h"
#include "CellularAutomata/Solvers/AMXSolver.h"
#include "Memory/Allocators/CPUAllocator.h"

AMXSolver::AMXSolver(CADataDomain<uint8_t>* domain, CADataDomain<uint8_t>* domainBuffer) {
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

	lDebug(1, "Setting up AMX");
    setupAMX();
    
	lDebug(1, "Filing Tridiag matrices");
    fillTridiag();
}

void AMXSolver::setupAMX() {
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


void AMXSolver::fillTridiag() {
	uint8_t* data = (uint8_t*) malloc(64*16*12*sizeof(uint8_t));
	pi_1  = &data[64*16*0];
	pi_2  = (uint8_t*) malloc(64*16*sizeof(uint8_t));
	pi_3  = &data[64*16*2];
	pi_4  = &data[64*16*3];
	pi_5  = &data[64*16*4];
	pi_6  = &data[64*16*5];
	pi_1B  = &data[64*16*6];
	pi_2B  = (uint8_t*) malloc(64*16*sizeof(uint8_t));
	//pi_2B  = &data[64*16*7];
	pi_3B  = &data[64*16*8];
	pi_4B  = &data[64*16*9];
	pi_5B  = &data[64*16*10];
	pi_6B  = &data[64*16*11];
    int i;

    for (i = 0; i < 16*64; i += 1)
    {
        int col = i & 63;
        int row = i >> 6;
        if (col + 63 - RADIUS < row + 16*3){
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
        if (abs(col - row-16) <= RADIUS){
            pi_3[i] = 1;
        } else {
            pi_3[i] = 0;
        }
        if (abs(col - row-32) <= RADIUS){
            pi_4[i] = 1;
        } else {
            pi_4[i] = 0;
        }
        if (abs(col - row-48) <= RADIUS){
            pi_5[i] = 1;
        } else {
            pi_5[i] = 0;
        }

        if (col -63+RADIUS> row){
            pi_6[i] = 1;
        } else {
            pi_6[i] = 0;
        }
    }

	uint32_t *pi1_int = (uint32_t*)pi_1; 
	uint32_t *pi1B_int = (uint32_t*)pi_1B; 
	uint32_t *pi2_int = (uint32_t*)pi_2; 
	uint32_t *pi2B_int = (uint32_t*)pi_2B; 
	uint32_t *pi3_int = (uint32_t*)pi_3; 
	uint32_t *pi3B_int = (uint32_t*)pi_3B; 
	uint32_t *pi4_int = (uint32_t*)pi_4; 
	uint32_t *pi4B_int = (uint32_t*)pi_4B; 
	uint32_t *pi5_int = (uint32_t*)pi_5; 
	uint32_t *pi5B_int = (uint32_t*)pi_5B; 
	uint32_t *pi6_int = (uint32_t*)pi_6; 
	uint32_t *pi6B_int = (uint32_t*)pi_6B; 
    // // debug print tridiag in 2d
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
			pi1B_int[j*16+i] = pi1_int[i*16+j];
			pi2B_int[j*16+i] = pi2_int[i*16+j];
			pi3B_int[j*16+i] = pi3_int[i*16+j];
			pi4B_int[j*16+i] = pi4_int[i*16+j];
			pi5B_int[j*16+i] = pi5_int[i*16+j];
			pi6B_int[j*16+i] = pi6_int[i*16+j];
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

void AMXSolver::copyCurrentStateToHostVisibleData() {
	lDebug(1, "%llu\n", dataDomain->getTotalSize());
    for (int i = 0; i < dataDomain->getTotalSize(); ++i) {
        uint8_t value = dataDomain->getElementAt(i);
        hostVisibleData->setElementAt(i, (int)value);
    }
}
void AMXSolver::copyHostVisibleDataToCurrentState() {
    for (int i = 0; i < hostVisibleData->getTotalSize(); ++i) {
        int value = hostVisibleData->getElementAt(i);
        dataDomain->setElementAt(i, value);
    }
}

void AMXSolver::swapPointers() {
    CADataDomain<uint8_t>* temp = dataDomain;
    dataDomain = dataDomainBuffer;
    dataDomainBuffer = temp;
}

uint8_t AMXSolver::transitionFunction(int k, int a, int b) {
    return (1 - (((k - a) >> 31) & 0x1)) * (1 - (((b - k) >> 31) & 0x1));
}

int buffer[64*64];
uint8_t buffer2[16*64];
uint8_t hardcoded[16*64];


void AMXSolver::CAStepAlgorithm() {
    uint8_t* data = dataDomain->getData();
    uint8_t* dataBuffer = dataDomainBuffer->getData();
    uint8_t* dataI = dataDomainIntermediate->getData();
    size_t nWithHalo = dataDomain->getFullHorizontalSize();

    _tile_loadd(1, pi_1B, 64);
    _tile_loadd(2, pi_2B, 64);
    _tile_loadd(3, pi_3B, 64);
    _tile_loadd(4, pi_4B, 64);
    _tile_loadd(5, pi_5B, 64);
    _tile_loadd(6, pi_6B, 64);


    //FIRST STEP: horizontal reduction
    for (size_t i = 0; i < nWithHalo; i+=64) {
        for (size_t j = 0; j < nWithHalo - 64*2; j+=64) {
            //take three continuous 16x16 blocks and load them into amx
			//printf("%i,%i\n", i, j);

		    // tiles in C from 0,0 to 3,0 (first col)
			for (int k = 0; k < 4; k++){
            	_tile_zero(0);
				_tile_loadd(7, &data[(i+16*k) * nWithHalo + j], nWithHalo);
				_tile_dpbssd(0, 7, 6);
				_tile_loadd(7, &data[(i+16*k) * nWithHalo + j + 64], nWithHalo);
				_tile_dpbuud(0, 7, 2);
				_tile_stored(0, &buffer[(k*16)*64 + 16*0], 64*4);
				//_tile_stored(0, &buffer, 16*4);
				//_tile_stored(2, &buffer2, 64);
				//printf("data[%i] and buffer[%i] \n",  (i+16*k) * nWithHalo + j + 64, k*16*64 + 16*0);
			}
			
			for (int k  = 0; k < 4; k++){
            	_tile_zero(0);
				_tile_loadd(7, &data[(i+16*k) * nWithHalo + j + 64], nWithHalo);
				_tile_dpbssd(0, 7, 3);
				_tile_stored(0, &buffer[k*16*64 + 16*1], 16*4*4);
			}
			for (int k  = 0; k < 4; k++){
            	_tile_zero(0);
				_tile_loadd(7, &data[(i+16*k) * nWithHalo + j + 64], nWithHalo);
				_tile_dpbssd(0, 7, 4);
				_tile_stored(0, &buffer[k*16*64 + 16*2], 16*4*4);
			}
			for (int k  = 0; k < 4; k++){
            	_tile_zero(0);
				_tile_loadd(7, &data[(i+16*k) * nWithHalo + j + 128], nWithHalo);
				_tile_dpbssd(0, 7, 1);
				_tile_loadd(7, &data[(i+16*k) * nWithHalo + j + 64], nWithHalo);
				_tile_dpbssd(0, 7, 5);
				_tile_stored(0, &buffer[k*16*64 + 16*3], 16*4 *4);

			}
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
			for (int ii=0; ii<64; ii++){
				for (int jj=0; jj<64; jj++){
						dataI[(i+jj)*nWithHalo + j + ii + 64] = buffer[ii*64+jj];
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
    //SECOND STEP: vertical reduction
    for (int i = 0; i < nWithHalo - 64*2; i+=64) {
        for (int j = 0; j < nWithHalo - 64*2; j+=64) {
            //take three continuous 16x16 blocks and load them into amx
			for (int k = 0; k < 4; k++){
            	_tile_zero(0);
				_tile_loadd(7, &dataI[(i+16*k) * nWithHalo + j +64], nWithHalo);
				_tile_dpbuud(0, 7, 6);
				_tile_loadd(7, &dataI[(i+64+16*k) * nWithHalo + j + 64], nWithHalo);
				_tile_dpbuud(0, 7, 2);
				_tile_stored(0, &buffer[(k*16)*64 + 16*0], 64*4);
				//_tile_stored(0, &buffer, 16*4);
				//_tile_stored(2, &buffer2, 64);
				//printf("data[%i] and buffer[%i] \n",  (i+16*k) * nWithHalo + j + 64, k*16*64 + 16*0);
			}
			
			for (int k  = 0; k < 4; k++){
            	_tile_zero(0);
				_tile_loadd(7, &dataI[(i+64+16*k) * nWithHalo + j + 64], nWithHalo);
				_tile_dpbssd(0, 7, 3);
				_tile_stored(0, &buffer[k*16*64 + 16*1], 16*4*4);
			}
			for (int k  = 0; k < 4; k++){
            	_tile_zero(0);
				_tile_loadd(7, &dataI[(i+64+16*k) * nWithHalo + j + 64], nWithHalo);
				_tile_dpbssd(0, 7, 4);
				_tile_stored(0, &buffer[k*16*64 + 16*2], 16*4*4);
			}
			for (int k  = 0; k < 4; k++){
            	_tile_zero(0);
				_tile_loadd(7, &dataI[(i+128+16*k) * nWithHalo + j + 64], nWithHalo);
				_tile_dpbssd(0, 7, 1);
				_tile_loadd(7, &dataI[(i+64+16*k) * nWithHalo + j + 64], nWithHalo);
				_tile_dpbssd(0, 7, 5);
				_tile_stored(0, &buffer[k*16*64 + 16*3], 16*4 *4);

			}

			for (int ii=0; ii<64; ii++){
				for (int jj=0; jj<64; jj+=1){
             		uint8_t cellValue = data[(i+ii + 64)*nWithHalo + j + 64 + jj];
             		int liveNeighbors = buffer[jj*64+ii] - cellValue;

             	    uint8_t result = cellValue * transitionFunction(liveNeighbors, SMIN, SMAX) + (1 - cellValue) * transitionFunction(liveNeighbors, BMIN, BMAX);
					dataBuffer[(i+ii + 64)*nWithHalo + j + 64 + jj] = result;
					//printf("%i ", data[(i+ii + 64)*nWithHalo + j + 64 + jj]);
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

int AMXSolver::countAliveNeighbors(int y, int x) {
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

void AMXSolver::fillHorizontalBoundaryConditions() {
    for (int h = 0; h < dataDomain->getHorizontalHaloSize(); ++h) {
        for (int j = 0; j < dataDomain->getInnerHorizontalSize(); ++j) {
            size_t topIndex = (dataDomain->getHorizontalHaloSize() + h) * dataDomain->getFullHorizontalSize() + dataDomain->getHorizontalHaloSize() + j;
            size_t bottomIndex = topIndex + (dataDomain->getInnerHorizontalSize()) * dataDomain->getFullHorizontalSize();
            uint8_t value = dataDomain->getElementAt(topIndex);
            dataDomain->setElementAt(bottomIndex, value);
        }

        for (int j = 0; j < dataDomain->getInnerHorizontalSize(); ++j) {
            size_t topIndex = (h)*dataDomain->getFullHorizontalSize() + dataDomain->getHorizontalHaloSize() + j;
            size_t bottomIndex = topIndex + (dataDomain->getInnerHorizontalSize()) * dataDomain->getFullHorizontalSize();

            uint8_t value = dataDomain->getElementAt(bottomIndex);
            dataDomain->setElementAt(topIndex, value);
        }
    }
}

void AMXSolver::fillVerticalBoundaryConditions() {
    for (int h = 0; h < dataDomain->getHorizontalHaloSize(); ++h) {
        for (int i = 0; i < dataDomain->getFullHorizontalSize(); ++i) {
            size_t leftIndex = i * dataDomain->getFullHorizontalSize() + h;
            size_t rightIndex = leftIndex + dataDomain->getInnerHorizontalSize();
            uint8_t value = dataDomain->getElementAt(rightIndex);
            dataDomain->setElementAt(leftIndex, value);
        }

        for (int i = 0; i < dataDomain->getFullHorizontalSize(); ++i) {
            size_t leftIndex = i * dataDomain->getFullHorizontalSize() + dataDomain->getHorizontalHaloSize() + h;
            size_t rightIndex = leftIndex + dataDomain->getInnerHorizontalSize();
            uint8_t value = dataDomain->getElementAt(leftIndex);
            dataDomain->setElementAt(rightIndex, value);
        }
    }
}
