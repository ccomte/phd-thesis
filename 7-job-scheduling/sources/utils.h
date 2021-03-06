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
#define R 100           // number of independent runs for each data point
#define WARMUP 1000000  // number of warm-up steps
#define STEADY 1000000  // number of steps in steady state
#define MAX 100000      // queue length

// constants that depend on the example
#define I 2             // number of pools (only relevant for the toy examples)
#define S 100           // number of servers

// probability distributions
double uniform();     // uniform on (0,1)
int bernoulli (double p);
double exponential (double lambda);

#endif
