#!/bin/bash
gcc ej1_secuencial.c -o ej1_secuencial
gcc ej1_paralelo.c -o ej1_paralelo

./ej1_secuencial 1000000000
./ej1_paralelo 1000000000 100