#include "utils.h"

double uniform () {
  return((double) (random() % (RAND_MAX-2) + 1) / RAND_MAX);
}

int bernoulli (double p) {
  return (uniform() < p);
}

double exponential (double lambda) {
  return -1. / lambda * log(uniform());
}
