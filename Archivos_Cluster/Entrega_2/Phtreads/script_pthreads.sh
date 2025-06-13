#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o /nethome/spusuario9/Entrega_2/Phtreads/output.txt
#SBATCH -e /nethome/spusuario9/Entrega_2/Phtreads/errors.txt
gcc -O3 MatMulP.c -o MatMulP -lpthread
./MatMulP 512 8
