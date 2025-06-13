#!/bin/bash
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --tasks-per-node=8

for N in 10000000 20000000 40000000; do
	output_file_blocking="/nethome/spusuario9/Entrega_3/Punto1/Ejercicio3/Salidas/output_blocking_ocho_${N}.txt"
	error_file_blocking="/nethome/spusuario9/Entrega_3/Punto1/Ejercicio3/Salidas/errors_blocking_ocho_${N}.txt"
	output_file_non_blocking="/nethome/spusuario9/Entrega_3/Punto1/Ejercicio3/Salidas/output_non_blocking_ocho_${N}.txt"
	error_file_non_blocking="/nethome/spusuario9/Entrega_3/Punto1/Ejercicio3/Salidas/errors_non_blocking_ocho_${N}.txt"
	mpirun blocking-ring "$N" > "$output_file_blocking" 2> "$error_file_blocking"
    	mpirun non-blocking-ring "$N" > "$output_file_non_blocking" 2> "$error_file_non_blocking"
done
