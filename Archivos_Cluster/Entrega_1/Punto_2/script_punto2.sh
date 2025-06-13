#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o /nethome/spusuario9/Entrega_1/Punto_2/output.txt
#SBATCH -e /nethome/spusuario9/Entrega_1/Punto_2/errors.txt
gcc -O3 punto2_modificado.c -lm -o punto2_modificado
./punto2_modificado
