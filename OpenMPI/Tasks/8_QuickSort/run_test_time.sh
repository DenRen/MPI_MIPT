#!/bin/bash

NUM_REPEATS=2

for((n=1;n<8;n++))
do
	echo "num threads: $n"
	for((i=0;i<$NUM_REPEATS;i++))
	do
		OMP_NUM_THREADS=$n ./rand_arr.out
	done
done

