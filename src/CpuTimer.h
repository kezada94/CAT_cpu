#pragma once

#include <chrono>

class CpuTimer {
public:
    CpuTimer();
    ~CpuTimer();
    void start();
    void stop();
    float getElapsedTimeMilliseconds() const;

private:
    std::chrono::high_resolution_clock::time_point startTime;
    std::chrono::high_resolution_clock::time_point stopTime;
    float elapsedTimeMilliseconds;
};
