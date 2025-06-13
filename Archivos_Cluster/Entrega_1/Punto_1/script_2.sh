#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o /nethome/spusuario9/Practica_1/outputs/output_ej2.txt
#SBATCH -e /nethome/spusuario9/Practica_1/errors/errors_ej2.txt
gcc -O3 quadatric2.c -lm -o quadatric2
./quadatric2
