#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o /nethome/spusuario9/Entrega_2/Secuencial/output.txt
#SBATCH -e /nethome/spusuario9/Entrega_2/Secuencial/errors.txt
gcc -O3 MatMul.c -lm -o MatMul
./MatMul {N} {BlockSize (64)}
