FYS4150 - Project 4

Data: Contains text-files which are produced by the main.cpp and funcs.cpp programs to be analyzed.

Figures: Contains figures of the data in the text-files. These figures are produced by the plots.py program.

Report: Contains the report files for this project.

Codes: Contains all the code for this project. All the code was done on a Windows 10 computer with 8GB RAM. The C++ code is run in Qt Creator, and was compiled as: `c++ -floop-parallelize-all -o main.exe main.cpp funcs.cpp` in the Qt Creator Run-terminal. (Did not manage to link armadillo to the used ming-terminal coming with Qt Creator.)
- plots.py: Creates all the figures in the Figure folder. This was done in Python 3.7.
- funcs.cpp: Contains the Metropolis algorithm function, a function to initialize the expectation values, a function for the probability and functions for writing the results file files to be exported to the Data folder.
- main.cpp: The main program for running the calculations and functions from funcs.cpp.

Code was not properly parallelized. Tried to use openMP, but for some reason I got an error for armadillo usage when uncommenting the -fopenmp flag in the .pro file. So this code is only vectorized with `-O3` and `-floop-parallelize-all`.
