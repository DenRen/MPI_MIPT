#!/bin/bash

NUM_REPEATS=3

TIME_PROG='/usr/bin/time -f "%E"'

echo "One thread (time of execution):"
for((i=0;i<$NUM_REPEATS;i++))
do
	${TIME_PROG} ./single_thread.out
done

echo
echo "Mult threads (time of execution):"
for((i=0;i<$NUM_REPEATS;i++))
do
	OMP_MAX_ACTIVE_LEVELS=2 \
	OMP_NESTED=TRUE \
	OMP_NUM_THREADS=2 \
	${TIME_PROG} ./mult_thread.out
done

