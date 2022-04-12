#!/bin/bash

make main.out
sleep 2

if [ -e main.out ]; then

    for (( num_threads = 3; num_threads <= $(nproc); num_threads++ ))
    do
        echo -n "num treads: $num_threads, time: "
        /bin/time -f "%e" mpirun -n $num_threads ./main.out
    done
fi
