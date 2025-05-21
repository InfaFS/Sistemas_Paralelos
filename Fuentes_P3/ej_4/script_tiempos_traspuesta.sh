#!/bin/bash

# Tamaño de la matriz (N)
N=30000

# Compilar el código con OpenMP
echo "Compilando..."
gcc-14 -fopenmp traspuesta.c -o traspuesta
gcc-14 -fopenmp traspuesta_sec.c -o traspuesta_sec


# Corremos el algoritmo secuencial
./traspuesta_sec $N 1

# Corremos el algoritmo paralelo con 4 hilos
echo -e "\n==== TIEMPO POR HILOS ====\n"
#for T in 1 2 4 8
#do
T=4
  echo -n "Hilos: $T -> "
  ./traspuesta $N $T
#done

