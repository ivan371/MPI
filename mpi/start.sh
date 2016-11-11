#!/bin/bash

mpicc $3

for (( i = $1; i < $2; i++))
do

mpirun -np $i  ./a.out

done
