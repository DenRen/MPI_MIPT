#!/bin/bash

g++ -fopenmp -O3 -DRANDOM_ARRAY -o rand_arr.out main.cpp
g++ -fopenmp -O3 -DSTDIO_ARRAY -o stdio_arr.out main.cpp

