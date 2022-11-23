#!/bin/bash

g++ main.cpp -O3 -o single_thread.out
g++ -fopenmp main.cpp -O3 -o mult_thread.out

