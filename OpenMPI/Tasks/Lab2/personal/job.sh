#!/bin/bash
#PBS -l walltime=00:10:00,nodes=7:ppn=3
#PBS -N denren
#PBS -q batch

for ((i=1;i<32;i++))
do
	mpirun -np $i /home/b0190120/personal/main.out
done
