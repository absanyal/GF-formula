#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -N 12_6
#PBS -j oe
#PBS -e error.log
#PBS -o out.log
cd $PBS_O_WORKDIR
module load python/3.5
python3 spinless.py
