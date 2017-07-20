// Authors:
// Emanuele Del Sozzo (emanuele.delsozzo@polimi.it), Marcello Pogliani (marcello.pogliani@polimi.it)

#ifndef BENCHMARK_UTILS_H
#define BENCHMARK_UTILS_H

#ifdef BENCHMARK
#define BENCH_GETTIME(x) do {       \
    gettimeofday((x), NULL);    \
      } while(0)
#else
#define BENCH_GETTIME(x)
#endif

#include <sys/time.h>

void add_duration(double val, char* why);
void print_duration(struct timeval start, struct timeval end, char* why);
void print_durations();
void print_causes();

#endif
