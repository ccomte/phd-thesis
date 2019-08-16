#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <unistd.h>   /* for getopt */
#include <ctype.h>    /* for isprint */
#include <stdint.h>   /* unsigned integer type uint64_t */



#define MAX_I 10      /* maximum number of computers */
#define MAX_K 10      /* maximum number of job types */
#define MAX_R 100     /* maximum number of independent simulation runs per load */
#define MAX_L 128     /* length of the circular buffer
                         encoding the queue of available tokens */



static inline double uniform () {
  return((double) (random() % (RAND_MAX-2) + 1) / RAND_MAX);
}

static inline int bernoulli (double p) {
  return (uniform() < p);
}

static inline double exponential (double lambda) {
  return -1. / lambda * log(uniform());
}



#endif
