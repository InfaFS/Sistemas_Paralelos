#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o /nethome/spusuario9/Practica_1/outputs/output_ej1.txt
#SBATCH -e /nethome/spusuario9/Practica_1/errors/errors.txt
gcc -O3 quadatric1.c -lm -o quadatric1
./quadatric1
