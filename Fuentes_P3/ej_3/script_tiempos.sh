#!/bin/bash

# Tamaño de la matriz (N)
N=1000

# Compilar el código con OpenMP
echo "Compilando..."
gcc-14 -fopenmp matrices_copy_1.c -o matrices_copy_1
gcc-14 -fopenmp matrices_copy_2.c -o matrices_copy_2

# Tiempos por filas y por columnas
echo -e "\n==== TIEMPOS POR FILAS ====\n"
for T in 1 2 4 8
do
  echo -n "Hilos: $T -> "
  ./matrices_copy_1 $N $T | grep "Tiempo en segundos"
done

echo -e "\n==== TIEMPOS POR COLUMNAS ====\n"
for T in 1 2 4 8
do
  echo -n "Hilos: $T -> "
  ./matrices_copy_2 $N $T | grep "Tiempo en segundos"
done

