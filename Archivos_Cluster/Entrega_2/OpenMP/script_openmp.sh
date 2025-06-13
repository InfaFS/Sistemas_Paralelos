#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o /nethome/spusuario9/Entrega_2/OpenMP/output.txt
#SBATCH -e /nethome/spusuario9/Entrega_2/OpenMP/errors.txt
gcc -fopenmp -O3 MatMulOMP.c -o MatMulOMP
./MatMulOMP 512 2
