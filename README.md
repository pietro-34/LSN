This repository contains 12 folders (Exercise_#), each corresponding to an exercise from the Numerical Simulation Laboratory course at the University of Milan (Unimi).

Exercises 1–9:

Contain subdirectories with C++ source code and a Makefile to compile and run the programs.

Each folder also includes two Jupyter notebooks:

LSN_Exercise_#.ipynb (provided template)

exercise_#.ipynb (my solution)

Exercise 10:

Written in C++ with OpenMPI.

Includes the two notebooks (LSN_Exercise_10.ipynb and exercise_10.ipynb).

The code cannot be executed with the standard make command; it must be compiled and run with MPI, for example:

mpic++ source.cpp -o output.x
mpirun -np <num_processes> ./output.x


Exercises 11–12:

Implemented entirely in Python.

Each folder contains:

LSN_Exercise_#.ipynb (provided template)

exercise_#.ipynb (my solution)
