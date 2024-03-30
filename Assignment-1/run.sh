#!/bin/bash

make

python3 /home/aashi/Desktop/CS633/src.c start
chmod u+x /home/aashi/Desktop/CS633/src.x
/home/aashi/Desktop/CS633/src.x 12 6

mkdir timing_data

for execution in 1 2 3 
do
        for Px in 4
        do
                for N in 262144 4194304
                do
                       for stencil in 5 9
                       do
                       	mpirun -np 12 --hostfile myhostfile ./src.x $Px $N 10 1 $stencil >> ./Datafile.txt
                       done
                done
        done
done

python3 boxplot.py
