#!/bin/bash



# compilation

gcc -O3 -o scheduling-and-balancing-exact.out scheduling-and-balancing-exact.c -lm



# folder
mkdir -p data



# a single job type

./scheduling-and-balancing-exact.out -f "data/single-dynamic-exact-3" \
  -k 1 -i 6 -s 10 -l 10 \
  -c "1 1 1 1 1 1
      1 1 1 0 0 0 0 0 0 0
      0 1 1 1 0 0 0 0 0 0
      0 0 1 1 1 0 0 0 0 0
      0 0 0 0 0 1 1 1 0 0
      0 0 0 0 0 0 1 1 1 0
      0 0 0 0 0 0 0 1 1 1" \
  -r "1
      1 1 1 1 1 4 4 4 4 4" \
  -m 3 -d .005

./scheduling-and-balancing-exact.out -f "data/single-dynamic-exact-4" \
  -k 1 -i 4 -s 10 -l 15 \
  -c "1 1 1 1
      1 1 1 1 0 0 0 0 0 0
      0 1 1 1 1 0 0 0 0 0
      0 0 0 0 0 1 1 1 1 0
      0 0 0 0 0 0 1 1 1 1 " \
  -r "1
      1 1 1 1 1 4 4 4 4 4" \
  -m 3 -d .005

./scheduling-and-balancing-exact.out -f "data/single-dynamic-exact-5" \
  -k 1 -i 2 -s 10 -l 30 \
  -c "1 1
      1 1 1 1 1 0 0 0 0 0
      0 0 0 0 0 1 1 1 1 1" \
  -r "1
      1 1 1 1 1 4 4 4 4 4" \
  -m 3 -d .005
