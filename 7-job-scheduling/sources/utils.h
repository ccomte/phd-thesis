#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>

// global constants
#define MAX 100000      // queue length
#define R 10            // number of estimates per point
#define WARMUP 100000   // number of warm-up steps
#define STEADY 100000   // number of steps in steady state

// constants that depend on the example
#define N 2             // number of pools (only relevant for the toy examples)
#define S 100           // number of servers

// probability distributions
double uniform();     // uniform on (0,1)
int bernoulli (double p);
double exponential (double lambda);

#endif
