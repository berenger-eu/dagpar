#ifndef TIMER_H
#define TIMER_H

#include <chrono>

class Timer {
    using double_second_time = std::chrono::duration<double, std::ratio<1, 1>>;

    std::chrono::high_resolution_clock::time_point m_start;
    std::chrono::high_resolution_clock::time_point m_end;
    std::chrono::nanoseconds m_cumulate;

public:
    Timer() { start(); }

    Timer(const Timer &other) = delete;
    Timer &operator=(const Timer &other) = delete;
    Timer(Timer &&other) = default;
    Timer &operator=(Timer &&other) = default;

    void reset() {
        m_start = std::chrono::high_resolution_clock::time_point();
        m_end = std::chrono::high_resolution_clock::time_point();
        m_cumulate = std::chrono::nanoseconds();
        start();
    }

    void start() { m_start = std::chrono::high_resolution_clock::now(); }

    void stop() {
        m_end = std::chrono::high_resolution_clock::now();
        m_cumulate +=
                std::chrono::duration_cast<std::chrono::nanoseconds>(m_end - m_start);
    }

    double getElapsed() const {
        return std::chrono::duration_cast<double_second_time>(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(m_end -
                                                                         m_start))
                .count();
    }

    double getCumulated() const {
        return std::chrono::duration_cast<double_second_time>(m_cumulate).count();
    }

    double stopAndGetElapsed() {
        stop();
        return getElapsed();
    }
};

#endif
