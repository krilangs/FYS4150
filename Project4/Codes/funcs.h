#ifndef FUNCS_H
#define FUNCS_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
#include <chrono>

using namespace  std;
using namespace arma;


// Function to initialise energy and magnetization
void InitLattice(int, mat &, double &, double &, bool);

// The metropolis algorithm with loop over Monte Carlo cycles
void Metropolis(int, int, double, vec &, int &, bool, vec &, vec &);

// Counts the energy states/calculates the probability
void Probability(double, vec &, vec &);

// The following prints to file the results of the calculations
void WriteAnalytical(ofstream &);
void WriteNumerical(ofstream &, int, double, vec);
void WriteEquilibrium(ofstream &, int, int, double, vec, int);
void WriteNaccept(ofstream &, int, int, double, int);
void Writeprobs(ofstream &, vec, vec, int, vec);
void WritePhases(ofstream &, int, int, double, vec);
void WriteTC(ofstream &, mat, int, int, vec);

#endif // ISING_H
