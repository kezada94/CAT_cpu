#pragma once

#include <cstdio>
#include <cstdint>
#include "Memory/CADataDomain.h"
#include "Debug.h"

class CADataPrinter {
   public:
    static void printCADataWithHalo(CADataDomain<int>* data) {
        for (int i = 0; i < data->getFullVerticalSize(); i++) {
            for (int j = 0; j < data->getFullHorizontalSize(); j++) {
                size_t index = i * data->getStride() + j;

                int element = data->getElementAt(index);
                if (element == 0) {
                    printf("- ");
                } else {
                    printf("%d ", element);
                }
            }
            printf("\n");
        }
        printf("\n");
    }
    static void printCADataWithoutHalo(CADataDomain<int>* data) {
	lDebug(1, "PRintn \n");
        for (int i = 0; i < data->getInnerHorizontalSize(); i++) {
            for (int j = 0; j < data->getInnerHorizontalSize(); j++) {
                int element = data->getInnerElementAt(i, j);
                if (element == 0) {
                    printf("- ");
                } else {
                    printf("%d ", element);
                }
            }
            printf("\n");
        }
        printf("\n");
    }
    static void printCADataWithHalo(CADataDomain<uint8_t>* data) {
        for (int i = 0; i < data->getFullVerticalSize(); i++) {
            for (int j = 0; j < data->getFullHorizontalSize(); j++) {
                size_t index = i * data->getStride() + j;

                int element = data->getElementAt(index);
                if (element == 0) {
                    printf("- ");
                } else {
                    printf("%d ", element);
                }
            }
            printf("\n");
        }
        printf("\n");
    }
    static void printCADataWithoutHalo(CADataDomain<uint8_t>* data) {
        for (int i = 0; i < data->getInnerHorizontalSize(); i++) {
            for (int j = 0; j < data->getInnerHorizontalSize(); j++) {
                int element = data->getInnerElementAt(i, j);
                if (element == 0) {
                    printf("- ");
                } else {
                    printf("%d ", element);
                }
            }
            printf("\n");
        }
        printf("\n");
    }
};
