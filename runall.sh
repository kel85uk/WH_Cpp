#!/bin/bash
for ((expo=8;expo<=8;expo++))
do
for ((numproc=1;numproc<=9;numproc++))
do
 mpirun -np 7 ./WH_test 1.${numproc}e${expo} 
done
done

