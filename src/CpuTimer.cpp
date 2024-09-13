#include "CpuTimer.h"

CpuTimer::CpuTimer() : elapsedTimeMilliseconds(0.0f) {}

CpuTimer::~CpuTimer() {}

void CpuTimer::start() {
    startTime = std::chrono::high_resolution_clock::now();
}

void CpuTimer::stop() {
    stopTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float, std::milli> elapsed = stopTime - startTime;
    elapsedTimeMilliseconds = elapsed.count();
}

float CpuTimer::getElapsedTimeMilliseconds() const {
    return elapsedTimeMilliseconds;
}