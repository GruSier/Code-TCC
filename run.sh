#! /bin/bash
#Para permitir a execução do script: chmod +x ./run.sh

gfortran $1.f90 -c
gfortran $1.o -o $1
./$1
