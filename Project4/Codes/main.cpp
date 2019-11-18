#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
#include <chrono>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "funcs.h"

using namespace std;
using namespace std::chrono;
using namespace arma;
ofstream ofile;

// Here we define various functions called by the main program
void ex_b();
void ex_c();
void ex_d();
void ex_e();

// Task b
void ex_b(){
    string file = "Data/Lattice_2X2";
    ofile.open(file);
    ofile << "Numerical values:\n |   # Monte Carlo cycles  | Energy-Mean | Magnetization-Mean |  Specific heat  | Susceptibility |\n";

    int NSpins = 2;
    int Nconfigs;
    long int MCcycles;
    double T = 1.0;

    cout << "Read in the number of Monte Carlo cycles" << endl;
    cin >> MCcycles;

    vec Energies = zeros<mat>(400);
    vec counter = zeros<mat>(400);

    // Start Monte Carlo sampling by looping over the selcted Temperatures
    for (int i = 10; i <= MCcycles; i *= 10){
        vec ExpectValues = zeros<mat>(5);

        // Start Monte Carlo computation and get expectation values
        Metropolis(NSpins, i, T, ExpectValues, Nconfigs, false, Energies, counter);

        // Write results to file
        WriteNumerical(ofile, i, T, ExpectValues);
    }
    ofile << "\n Analytical values:\n |   Partition function    | Energy-Mean | Magnetization-Mean |  Specific heat  | Susceptibility |\n";
    // Write results to file
    WriteAnalytical(ofile);
    ofile.close();  // Close output file
}

// Task c
void ex_c(){
    long int MC;
    double T;
    string file1, file2;

    cout << "Read in the number of Monte Carlo cycles" << endl;
    cin >> MC;
    cout << "Read in the value for the Temperature" << endl;
    cin >> T;
    cout << "Read in filename for ordered and random spins" << endl;
    cin >> file1 >> file2;

    vec Energies = zeros<mat>(400);
    vec counter = zeros<mat>(400);

    ofile.open("Data/"+file1);
    ofile << "|  # Monte Carlo cycles  | Energy-Mean | Magnetization-Mean | # Accepted configurations |  Specific heat  | Susceptibility | Temperature |\n";
    // Start Monte Carlo sampling by looping over the selcted Temperatures
    int N = 20;
    int Nconfigs;

    for (int i=1; i <= MC; i += 100){
        vec ExpectValue = zeros<mat>(5);

        // Start Monte Carlo computation and get expectation values
        Metropolis(N, i, T, ExpectValue, Nconfigs, false, Energies, counter);

        // Write results to file
        WriteEquilibrium(ofile, N, i, T, ExpectValue, Nconfigs);
    }
    ofile.close();  // Close output file

    ofile.open("Data/"+file2);
    ofile << "| # Monte Carlo cycles  |  Energy-Mean   |  Magnetization-Mean  |  # Accepted configurations  |  Specific heat  |  Susceptibility   |  Temperature |\n";

    for (int i=1; i <= MC; i += 100){
        vec ExpectValue2 = zeros<mat>(5);

        // Start Monte Carlo computation and get expectation values
        Metropolis(N, i, T, ExpectValue2, Nconfigs, true, Energies, counter);

        // Write results to file
        WriteEquilibrium(ofile, N, i, T, ExpectValue2, Nconfigs);
    }
    ofile.close();  // close output file

    // # Accepted configurations as a function of the Temperature
    string file3 = "Data/Nconfig_vs_Temp";

    ofile.open(file3);
    ofile << "|  Temperature  |  # Accepted configurations |\n";

    double InitialTemp = 1.0;
    double FinalTemp = 2.4;
    double TempStep = 0.1;

    for (double Temp = InitialTemp; Temp <= (FinalTemp+0.1); Temp+=TempStep){
        vec ExpectValue = zeros<mat>(5);

        // Start Monte Carlo computation and get expectation values
        Metropolis(N, 1000, Temp, ExpectValue, Nconfigs, true, Energies, counter);

        // Write results to file
        WriteNaccept(ofile, N, 1000, Temp, Nconfigs);
     }
     ofile.close();  // Close output file
}

// Task d
void ex_d(){
    long int MC;
    cout << "Read in the number of Monte Carlo cycles in times of 100, MC*100" << endl;
    cin >> MC;

    string file = "Data/Probability_1";
    vec Energies = zeros<mat>(400);
    vec counter = zeros<mat>(400);

    ofile.open(file);
    ofile << "|  Energies | Energy counts |\n";
    cout << "\n" << "Probability; with spin L = 20 and T = 1.0: " << endl;

    int N = 20;
    int Nconfigs;   // # Accepted of configurations
    double T = 1.0; // Temperature
    int iterations; // # MC cycles

    // Start Monte Carlo sampling by looping over the selcted Temperatures
    for (int i=1; i <= MC; i++){
        vec ExpectValue = zeros<mat>(5);
        iterations = 100*i;

        // Start Monte Carlo computation and get expectation values
        Metropolis(N, iterations, T, ExpectValue, Nconfigs, false, Energies, counter);

        if (i == MC){
            // Write results to file
            Writeprobs(ofile, Energies, counter, iterations, ExpectValue);
        }
    }
    ofile.close();  // Close output file
    cout << "\n" << "Probability; with spin L = 20 and T = 2.4: " << endl;

    string file2 = "Data/Probability_24";
    vec Energies2 = zeros<mat>(400);
    vec counter2 = zeros<mat>(400);

    T = 2.4; // Temperature
    ofile.open(file2);
    ofile << "|  Energies | Energy counts |\n";

    // Start Monte Carlo sampling by looping over the selcted Temperatures
    for (int i=1; i <= MC; i++){
        vec ExpectValue = zeros<mat>(5);
        iterations = 100*i;
        // Start Monte Carlo computation and get expectation values
        Metropolis(N, iterations, T, ExpectValue, Nconfigs, false, Energies2, counter2);

        if (i == MC){
            // Write results to file
            Writeprobs(ofile, Energies2, counter2, iterations, ExpectValue);
        }
    }
    ofile.close();  // Close output file
}

// Task e
void ex_e(){
    string file = "Data/Phase_transitions";
    ofile.open(file);
    ofile << "|  Temperature  |  Energy-Mean  |  Magnetization-Mean  |  Specific heat  |  Susceptibility  |  Lattice  |\n";

    int N_start, N_step, N_final, Nconfigs;
    long int MC;
    double T_start, T_step, T_final, T;

    cout << "Read in the number of Monte Carlo cycles" << endl;
    cin >> MC;

    vec Energies = zeros<mat>(400);
    vec counter = zeros<mat>(400);

    // Base: T=[2.0, 2.3] & stepsize=0.05
    cout << "Read in the start value for the Temperature" << endl;
    cin >> T_start;
    cout << "Read in the final value for the Temperature" << endl;
    cin >> T_final;
    cout << "Read in the stepsize for the Temperature" << endl;
    cin >> T_step;

    // Declare a matrix which stores the expectation values for spins 40, 60, 80, 100
    int mat_len = (T_final-T_start)/T_step + 2;

    N_start = 40;
    N_step = 20;
    N_final = 100;

    // Time the loop
    auto start = high_resolution_clock::now();
    vec Tvalues = zeros<mat>(mat_len);

    #pragma omp parallel for
    for (int N = N_start; N <= N_final; N += N_step){
        // Start Monte Carlo sampling by looping over the selcted Temperatures
        for (int i = 0; i <= (T_final-T_start)/T_step + 1; i++){
            vec ExpectValue = zeros<mat>(5);

            T = T_start + T_step*i;
            Tvalues(i) = T;

            // Start Monte Carlo computation and get expectation values
            Metropolis(N, MC, T, ExpectValue, Nconfigs, false, Energies, counter);

            // Write results to file
            WritePhases(ofile, N, MC, T, ExpectValue);
        }
    }
    auto stop = high_resolution_clock::now();
    auto time_used = duration_cast<seconds>(stop - start);
    cout << "Time used [s]: " << time_used.count() << endl;
}

// Main program begins here
int main(){
    string choice;
    string info = "Write any of the following commands:\n \
        B    : run Ising model 2x2 lattice\n \
        C    : run Equilibrium \n \
        D    : run Probability distribution \n \
        E    : run Phase transitions \n \
        info : write out this message\n \
        q    : quit program\n\n";
    cout << info;
    string number;

    int running = 0;
    while (running == 0){
    cout << ">> ";
    getline(cin, choice);

    // Call B) calculation
    if (choice.compare("B") == 0){
        ex_b();
    }
    // Call C) calculation
    else if (choice.compare("C") == 0){
        ex_c();

    }
    // Call D) calculation
    else if (choice.compare("D") == 0){
        ex_d();

    }
    // Call E) calculation
    else if (choice.compare("E") == 0){
        ex_e();

    }
    // Call to see calling info
    else if ((choice.compare("info") == 0) || (choice.compare("help") == 0)) {
        cout << info;
    }
    // Call to quit the program
    else if (choice.compare("q") == 0){
        running = 1;
        printf("Quit program!\n");
    }}

    return 0;
}
