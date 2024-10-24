#pragma once

#include "Debug.h"
#include "Memory/Allocators/Allocator.h"

template <typename T>
class CPUAllocator : Allocator<T> {
   public:
    T* allocate(size_t size) override;
    void deallocate(void* ptr) override;
};

#include "Memory/Allocators/CPUAllocator.tpp"