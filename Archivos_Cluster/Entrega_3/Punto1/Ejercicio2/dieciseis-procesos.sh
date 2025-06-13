#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --tasks-per-node=8

output_file="/nethome/spusuario9/Entrega_3/Punto1/Ejercicio2/Salidas/output_blocking_dieciseis.txt"
error_file="/nethome/spusuario9/Entrega_3/Punto1/Ejercicio2/Salidas/errors_blocking_dieciseis.txt"

mpirun blocking > "$output_file" 2> "$error_file"

output_file="/nethome/spusuario9/Entrega_3/Punto1/Ejercicio2/Salidas/output_non_blocking_dieciseis.txt"
error_file="/nethome/spusuario9/Entrega_3/Punto1/Ejercicio2/Salidas/errors_non_blocking_cuatro.txt"

mpirun non-blocking > "$output_file" 2> "$error_file"

output_file="/nethome/spusuario9/Entrega_3/Punto1/Ejercicio2/Salidas/output_non_blocking_no_wait_dieciseis.txt"
error_file="/nethome/spusuario9/Entrega_3/Punto1/Ejercicio2/Salidas/errors_non_blocking_no_wait_dieciseis.txt"

mpirun non-blocking-no-wait > "$output_file" 2> "$error_file"
