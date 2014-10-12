#!/bin/bash
#PBS -N 8-16
#PBS -j oe
#PBS -l mem=5gb
#PBS -l ncpus=1
cd $PBS_O_WORKDIR
make all
./test