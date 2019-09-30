#!/bin/bash



# compilation

gcc -O3 -o dynamic-exact.out sources/dynamic-exact.c -lm
gcc -O3 -o dynamic-exp.out sources/dynamic-exp.c -lm
gcc -O3 -o dynamic-hyperexp.out sources/dynamic-hyperexp.c -lm



# folder
mkdir -p data



# a single job type

./dynamic-exact.out -f "data/single-dynamic-exact" \
  -k 1 -i 10 -l 6 \
  -c "1 1 1 1 1 1 1 1 1 1" \
  -r "1 1 1 1 1 1 4 4 4 4 4" \
  -m 3 -d 1e-2

./dynamic-exp.out -f "data/single-dynamic-simu-exp" \
  -k 1 -i 10 -l 6 \
  -c "1 1 1 1 1 1 1 1 1 1" \
  -r "1 1 1 1 1 1 4 4 4 4 4" \
  -m 3 -d .2 -R 100 -w 6 -t 6

./dynamic-hyperexp.out -f "data/single-dynamic-simu-hyperexp" \
  -k 1 -i 10 -l 6 \
  -c "1 1 1 1 1 1 1 1 1 1" \
  -r "1 1 1 1 1 1 4 4 4 4 4" \
  -s "1 5" \
  -m 3 -d .2 -R 100 -w 6 -t 6



# several job types

./dynamic-exact.out -f "data/multi-dynamic-exact" \
  -k 2 -i 10 -l 6 \
  -c "1 0 1 0 1 0 1 1 1 1 1 1 1 1 0 1 0 1 0 1" -r "1 4 1 1 1 1 1 1 1 1 1 1" \
  -m 3 -d 1e-2

./dynamic-exp.out -f "data/multi-dynamic-simu-exp" \
  -k 2 -i 10 -l 6 \
  -c "1 0 1 0 1 0 1 1 1 1 1 1 1 1 0 1 0 1 0 1" -r "1 4 1 1 1 1 1 1 1 1 1 1" \
  -m 3 -d .2 -R 100 -w 6 -t 6

./dynamic-hyperexp.out -f "data/multi-dynamic-simu-hyperexp" \
  -k 2 -i 10 -l 6 \
  -c "1 0 1 0 1 0 1 1 1 1 1 1 1 1 0 1 0 1 0 1" -r "1 4 1 1 1 1 1 1 1 1 1 1" \
  -s "1 2 1 5" \
  -m 3 -d .2 -R 100 -w 6 -t 6
