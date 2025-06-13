#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o /nethome/spusuario9/Practica_1/outputs/output_ej3.txt
#SBATCH -e /nethome/spusuario9/Practica_1/errors/errors_ej3.txt
gcc -O3 quadatric3.c -lm -o quadatric3
./quadatric3
