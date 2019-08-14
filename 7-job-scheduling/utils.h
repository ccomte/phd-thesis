#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>

double uniform();   // uniform on (0,1)
int bernoulli (double p);
double exponential (double lambda);

#endif
