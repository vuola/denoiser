// Linked together with the original, unchanged filter source to build the
// one-worker reference program (makefile target "reference").
//
// The original program starts (number of processors - 1) worker threads and
// adds their partial sums in an order that depends on thread timing, so its
// result differs slightly from run to run. This file replaces the standard
// library's processor count with 1; the original source then starts a single
// worker and its result is the same on every run.

#include <thread>

unsigned int std::thread::hardware_concurrency() noexcept {
    return 1;
}
