#!/bin/bash


cp hubb.andpar_bak hubb.andpar
mpif90 aux_routines.f90 lattice_routines.f90 ed_dmft_parallel_frequencies.f90 -llapack -o run.x -g -ffixed-line-length-0
mpirun -np 36 --oversubscribe ./run.x
