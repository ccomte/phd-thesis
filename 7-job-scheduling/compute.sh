#!/bin/bash

mkdir -p data

echo 'Compile'
echo
gcc -o toy-bimodal-fcfs.out sources/toy-bimodal-fcfs.c sources/utils.c -lm
gcc -o toy-bimodal-scheduling.out sources/toy-bimodal-scheduling.c sources/utils.c -lm
gcc -o toy-hyperexp-fcfs.out sources/toy-hyperexp-fcfs.c sources/utils.c -lm
gcc -o toy-hyperexp-scheduling.out sources/toy-hyperexp-scheduling.c sources/utils.c -lm
gcc -o toy-zipf-fcfs.out sources/toy-zipf-fcfs.c sources/utils.c -lm
gcc -o toy-zipf-scheduling.out sources/toy-zipf-scheduling.c sources/utils.c -lm

gcc -o random-bimodal-fcfs.out sources/random-bimodal-fcfs.c sources/utils.c -lm
gcc -o random-bimodal-scheduling.out sources/random-bimodal-scheduling.c sources/utils.c -lm
gcc -o random-hyperexp-fcfs.out sources/random-hyperexp-fcfs.c sources/utils.c -lm
gcc -o random-hyperexp-scheduling.out sources/random-hyperexp-scheduling.c sources/utils.c -lm
gcc -o random-zipf-fcfs.out sources/random-zipf-fcfs.c sources/utils.c -lm
gcc -o random-zipf-scheduling.out sources/random-zipf-scheduling.c sources/utils.c -lm



# N = 2 pools, S = 3 servers
./toy-bimodal-fcfs.out data/m-bimodal-0 1 1 2 0 1 2 1 2 1 5 25 1
./toy-bimodal-scheduling.out data/m-bimodal-1 1 1 2 0 1 2 1 2 1 5 25 1 1
./toy-bimodal-scheduling.out data/m-bimodal-5 1 1 2 0 1 2 1 2 1 5 25 1 5

./toy-hyperexp-fcfs.out data/m-hyperexp-0 1 1 2 0 1 2 1 2 1 5 25 1
./toy-hyperexp-scheduling.out data/m-hyperexp-1 1 1 2 0 1 2 1 2 1 5 25 1 1
./toy-hyperexp-scheduling.out data/m-hyperexp-5 1 1 2 0 1 2 1 2 1 5 25 1 5

./toy-zipf-scheduling.out data/m-zipf-1 1 1 2 0 1 2 1 2 200 2 1
./toy-zipf-scheduling.out data/m-zipf-5 1 1 2 0 1 2 1 2 200 2 5
./toy-zipf-fcfs.out data/m-zipf-0 1 1 2 0 1 2 1 2 200 2



# N = 2 pools, S = 2 servers
./toy-bimodal-fcfs.out data/n-bimodal-0 1 1 2 0 1 1 1 1 5 25 1
./toy-bimodal-scheduling.out data/n-bimodal-1 1 1 2 0 1 1 1 1 5 25 1 1
./toy-bimodal-scheduling.out data/n-bimodal-5 1 1 2 0 1 1 1 1 5 25 1 5

./toy-hyperexp-fcfs.out data/n-hyperexp-0 1 1 2 0 1 1 1 1 5 25 1
./toy-hyperexp-scheduling.out data/n-hyperexp-1 1 1 2 0 1 1 1 1 5 25 1 1
./toy-hyperexp-scheduling.out data/n-hyperexp-5 1 1 2 0 1 1 1 1 5 25 1 5

./toy-zipf-fcfs.out data/n-zipf-0 1 1 2 0 1 1 1 200 2
./toy-zipf-scheduling.out data/n-zipf-1 1 1 2 0 1 1 1 200 2 1
./toy-zipf-scheduling.out data/n-zipf-5 1 1 2 0 1 1 1 200 2 5



# S = 100 servers
./random-bimodal-fcfs.out data/random-2-bimodal-0 2 1 5 25 1
./random-bimodal-scheduling.out data/random-2-bimodal-1 2 1 5 25 1 1
./random-bimodal-scheduling.out data/random-2-bimodal-5 2 1 5 25 1 5

./random-hyperexp-fcfs.out data/random-2-hyperexp-0 2 1 5 25 1
./random-hyperexp-scheduling.out data/random-2-hyperexp-1 2 1 5 25 1 1
./random-hyperexp-scheduling.out data/random-2-hyperexp-5 2 1 5 25 1 5

./random-zipf-fcfs.out data/random-2-zipf-0 2 200 2
./random-zipf-scheduling.out data/random-2-zipf-1 2 200 2 1
./random-zipf-scheduling.out data/random-2-zipf-5 2 200 2 5



# S = 100 servers
./random-bimodal-fcfs.out data/random-3-bimodal-0 3 1 5 25 1
./random-bimodal-scheduling.out data/random-3-bimodal-1 3 1 5 25 1 1
./random-bimodal-scheduling.out data/random-3-bimodal-5 3 1 5 25 1 5

./random-hyperexp-fcfs.out data/random-3-hyperexp-0 3 1 5 25 1
./random-hyperexp-scheduling.out data/random-3-hyperexp-1 3 1 5 25 1 1
./random-hyperexp-scheduling.out data/random-3-hyperexp-5 3 1 5 25 1 5

./random-zipf-fcfs.out data/random-3-zipf-0 3 200 2
./random-zipf-scheduling.out data/random-3-zipf-1 3 200 2 1
./random-zipf-scheduling.out data/random-3-zipf-5 3 200 2 5
