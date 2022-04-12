#!/bin/bash

make clean && make main.out
sleep 2

for (( num_threads = 1; num_threads <= $(nproc); num_threads++ ))
do
    echo -n "num treads: $num_threads, time: "
    /bin/time -f "%e" mpirun -n $num_threads ./main.out
done