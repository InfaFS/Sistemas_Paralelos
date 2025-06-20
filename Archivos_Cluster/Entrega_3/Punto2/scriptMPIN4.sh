#!/bin/bash
#SBATCH -N 4
#SBATCH --exclusive
#SBATCH --tasks-per-node=8

# Tamaños de matriz 
matrix_sizes=(512 1024 2048 4096)

# Bucle anidado para combinar los tamaños de matriz 
for matrix_size in "${matrix_sizes[@]}"; do
    mpirun MPIexe "$matrix_size" > "SalidaMPI/output_${matrix_size}_N4.txt" 2> "SalidaMPI/errors_${matrix_size}_N4.txt"
done
