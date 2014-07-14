#!/bin/bash
for ((numproc=2;numproc<=8;numproc++))
do
 mpirun -np $numproc ../build/exec/WH_test
done

