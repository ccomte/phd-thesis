#!/bin/bash

mkdir -p data

echo 'Compile'
echo
gcc -o bounds.out bounds.c

# 100 servers - 100 customer classes per type
# Type 1 has degree 1 - Type 2 has degree 2
./bounds.out 100 100 1 2 1e-3
./bounds.out 100 100 1 2 1e-4
./bounds.out 100 100 1 2 1e-5
./bounds.out 100 100 1 2 1e-6

# 500 servers - 500 customer classes per type
# Type 1 has degree 1 - Type 2 has degree 2
./bounds.out 500 500 1 2 1e-4
./bounds.out 500 500 1 2 1e-5
./bounds.out 500 500 1 2 1e-6
./bounds.out 500 500 1 2 1e-7

# 10000 servers - 1000 customer classes per type
# Type 1 has degree 20 - Type 2 has degree 40
./bounds.out 10000 1000 20 40 1e-4
./bounds.out 10000 1000 20 40 1e-5
./bounds.out 10000 1000 20 40 1e-6
