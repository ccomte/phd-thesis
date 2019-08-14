#!/bin/bash

mkdir -p data

echo 'Compile'
echo
gcc -o bimodal.out bimodal.c utils.c -lm
gcc -o hyperexp.out hyperexp.c utils.c -lm
gcc -o zipf.out zipf.c utils.c -lm

gcc -o random-hyperexp.out random-hyperexp.c utils.c -lm
gcc -o random-bimodal.out random-bimodal.c utils.c -lm
gcc -o random-zipf.out random-zipf.c utils.c -lm

gcc -o bimodal-fcfs.out bimodal-fcfs.c utils.c -lm
gcc -o hyperexp-fcfs.out hyperexp-fcfs.c utils.c -lm
gcc -o zipf-fcfs.out zipf-fcfs.c utils.c -lm

gcc -o random-bimodal-fcfs.out random-bimodal-fcfs.c utils.c -lm
gcc -o random-hyperexp-fcfs.out random-hyperexp-fcfs.c utils.c -lm
gcc -o random-zipf-fcfs.out random-zipf-fcfs.c utils.c -lm

gcc -o random-insensitive.out random-insensitive.c





# N = 2 job classes, S = 2 servers
#./zipf.out data/n-zipf-1-200 1 1 2 0 1 1 1 200 2 1
#./zipf.out data/n-zipf-5-200 1 1 2 0 1 1 1 200 2 5
#./zipf-fcfs.out data/n-zipf-fcfs-200 1 1 2 0 1 1 1 200 2

# N = 2 job classes, S = 3 servers
./bimodal-fcfs.out data/m-bimodal-0 1 1 2 0 1 2 1 2 1 5 25 1
./bimodal.out data/m-bimodal-1 1 1 2 0 1 2 1 2 1 5 25 1 1
./bimodal.out data/m-bimodal-5 1 1 2 0 1 2 1 2 1 5 25 1 5

#./zipf.out data/m-zipf-1-200 1 1 2 0 1 2 1 2 200 2 1
#./zipf.out data/m-zipf-5-200 1 1 2 0 1 2 1 2 200 2 5
#./zipf-fcfs.out data/m-zipf-fcfs-200 1 1 2 0 1 2 1 2 200 2

# S = 100 servers
#./random-zipf.out data/random-2-zipf-1-200 2 200 2 1
#./random-zipf.out data/random-2-zipf-5-200 2 200 2 5
#./random-zipf-fcfs.out data/random-2-zipf-fcfs-200 2 200 2

# S = 100 servers
#./random-zipf.out data/random-3-zipf-1-200 3 200 2 1
#./random-zipf.out data/random-3-zipf-5-200 3 200 2 5
#./random-zipf-fcfs.out data/random-3-zipf-fcfs-200 3 200 2
