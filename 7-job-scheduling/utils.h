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
#define N 2             // number of job classes
#define S 3             // number of servers
#define MAX 100000      // queue length
#define R 5             // number of estimates per point
#define WARMUP 10000    // number of warm-up steps
#define STEADY 10000    // number of steps in steady state

// probability distributions
double uniform();     // uniform on (0,1)
int bernoulli (double p);
double exponential (double lambda);

#endif
