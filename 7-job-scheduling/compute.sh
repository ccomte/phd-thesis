#!/bin/bash

mkdir -p data

echo 'Compile'
echo
gcc -o bimodal-scheduling.out bimodal-scheduling.c utils.c -lm
gcc -o bimodal-fcfs.out bimodal-fcfs.c utils.c -lm
gcc -o hyperexp-scheduling.out hyperexp-scheduling.c utils.c -lm
gcc -o hyperexp-fcfs.out hyperexp-fcfs.c utils.c -lm
gcc -o zipf-scheduling.out zipf-scheduling.c utils.c -lm
gcc -o zipf-fcfs.out zipf-fcfs.c utils.c -lm

#gcc -o random-bimodal-scheduling.out random-bimodal-scheduling.c utils.c -lm
#gcc -o random-bimodal-fcfs.out random-bimodal-fcfs.c utils.c -lm
#gcc -o random-hyperexp-scheduling.out random-hyperexp-scheduling.c utils.c -lm
#gcc -o random-hyperexp-fcfs.out random-hyperexp-fcfs.c utils.c -lm
#gcc -o random-zipf-scheduling.out random-zipf-scheduling.c utils.c -lm
#gcc -o random-zipf-fcfs.out random-zipf-fcfs.c utils.c -lm

#gcc -o random-insensitive.out random-insensitive.c





# N = 2 job classes, S = 3 servers
#./bimodal-fcfs.out data/m-bimodal-0 1 1 2 0 1 2 1 2 1 5 25 1
#./bimodal-scheduling.out data/m-bimodal-1 1 1 2 0 1 2 1 2 1 5 25 1 1
#./bimodal-scheduling.out data/m-bimodal-5 1 1 2 0 1 2 1 2 1 5 25 1 5

#./hyperexp-fcfs.out data/m-hyperexp-0 1 1 2 0 1 2 1 2 1 5 25 1
#./hyperexp-scheduling.out data/m-hyperexp-1 1 1 2 0 1 2 1 2 1 5 25 1 1
#./hyperexp-scheduling.out data/m-hyperexp-5 1 1 2 0 1 2 1 2 1 5 25 1 5

#./zipf-scheduling.out data/m-zipf-1 1 1 2 0 1 2 1 2 200 2 1
#./zipf-scheduling.out data/m-zipf-5 1 1 2 0 1 2 1 2 200 2 5
./zipf-fcfs.out data/m-zipf-0 1 1 2 0 1 2 1 2 200 2



# N = 2 job classes, S = 2 servers
#./zipf-scheduling.out data/n-zipf-1-200 1 1 2 0 1 1 1 200 2 1
#./zipf-scheduling.out data/n-zipf-5-200 1 1 2 0 1 1 1 200 2 5
#./zipf-fcfs.out data/n-zipf-fcfs-200 1 1 2 0 1 1 1 200 2



# S = 100 servers, degree 2
#./random-zipf-scheduling.out data/random-2-zipf-1-200 2 200 2 1
#./random-zipf-scheduling.out data/random-2-zipf-5-200 2 200 2 5
#./random-zipf-fcfs.out data/random-2-zipf-fcfs-200 2 200 2


# S = 100 servers, degree 3
#./random-zipf-scheduling.out data/random-3-zipf-1-200 3 200 2 1
#./random-zipf-scheduling.out data/random-3-zipf-5-200 3 200 2 5
#./random-zipf-fcfs.out data/random-3-zipf-fcfs-200 3 200 2
