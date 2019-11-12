#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
#include "omp.h"
#include "Functions.h"

using namespace  std;
using namespace arma;
// output file
ofstream ofile;



// Main program begins here

int main() {

cout << "\n" << "Which Project Task do you want to run?: " << endl;
cout << "\n" << "Project Task a & b: " <<  "Write b " << endl;
cout << "\n" << "Project Task c: " <<  "Write c " << endl;
cout << "\n" << "Project Task d: " <<  "Write d " << endl;
//cout << "\n" << "Project Task e: " <<  "Write e " << endl;


cout << "\n" << "Write here " << endl;
string Task;
cin >> Task;

//------------------------------------------------------------------------------------

if (Task == "b"){
cout << "\n" << "Project Task 4b: \n" << endl;

string file = "Lattice_2X2";
ofile.open(file);



ofile << "Numerical values:\n |   # Monte Carlo cycles  | Energy-Mean | Magnetization-Mean |  Specific heat  | Susceptibility |\n";

int NSpins;
int Nconfigs;
NSpins = 2;
long int MCcycles;
double T = 1.0;

cout << "Read in the number of Monte Carlo cycles" << endl;
cin >> MCcycles;



vec Energies = zeros<mat>(400);
vec counter = zeros<mat>(400);


  // Start Monte Carlo sampling by looping over the selcted Temperatures
  //for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature+=TempStep){

  for (int i = 10; i <= MCcycles; i *= 10){

  vec ExpectationValues = zeros<mat>(5);
  // Start Monte Carlo computation and get expectation values
  MetropolisSampling(NSpins, i, T, ExpectationValues, Nconfigs, false, Energies, counter);

  WriteResultsto4b(ofile, i, T, ExpectationValues);
  }
 // }
  ofile << "\n Analytical values:\n |   Partition function    | Energy-Mean | Magnetization-Mean |  Specific heat  | Susceptibility |\n";
  WriteAnalytical(ofile);
  ofile.close();  // close output file
}

//-------------------------------------------------------------------------------

// Task 4c

if (Task == "c"){

  cout << "\n" << "Project Task 4c for ordered and unordered with spin L = 20: " << endl;

  long int MC;
  double T;
  string file1;
  string file2;
  cout << "Read in the number of Monte Carlo cycles" << endl;
  cin >> MC;
  cout << "Read in the given value for the Temperature" << endl;
  cin >> T;
  cout << "Read in filename for ordered and random spins" << endl;
  cin >> file1 >> file2;


  vec Energies = zeros<mat>(400);
  vec counter = zeros<mat>(400);

  ofile.open(file1);
  //ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "|   # Monte Carlo cycles  | Energy-Mean | Magnetization-Mean | # Accepted configurations |  Specific heat  | Susceptibility | Temperature |\n";
  // Start Monte Carlo sampling by looping over the selcted Temperatures
  int N = 20;
  int Nconfigs;


  //#pragma omp parallel for
  for (int i=1; i <= MC; i += 100){
    vec ExpectationValue = zeros<mat>(5);

    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(N, i, T, ExpectationValue, Nconfigs, false, Energies, counter);
    //
    WriteResultstoFile(ofile, N, i, T, ExpectationValue, Nconfigs, false);
  }
  ofile.close();  // close output file

  ofile.open(file2);
  //ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "| # Monte Carlo cycles  |  Energy-Mean   |  Magnetization-Mean  |  # Accepted configurations  |  Specific heat  |  Susceptibility   |  Temperature |\n";

  //#pragma omp parallel for
  for (int i=1; i <= MC; i += 100){
    vec ExpectationValue2 = zeros<mat>(5);

    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(N, i, T, ExpectationValue2, Nconfigs, true, Energies, counter);
    //
    WriteResultstoFile(ofile, N, i, T, ExpectationValue2, Nconfigs, true);
  }
  ofile.close();  // close output file

  // # Accepted configurations as a function of the Temperature

  string file3 = "Nconfig_vs_Temp";

  ofile.open(file3);
  //ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "|  Temperature  |  # Accepted configurations |\n";

  double InitialTemp = 1.0;
  double FinalTemp = 2.4;
  double TempStep = 0.1;

  for (double Temp = InitialTemp; Temp <= (FinalTemp+0.1); Temp+=TempStep){
    vec ExpectationValue = zeros<mat>(5);

    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(N, 1000, Temp, ExpectationValue, Nconfigs, true, Energies, counter);

    WriteConfigvsT(ofile, N, 1000, Temp, Nconfigs);
   }
   ofile.close();  // close output file
}

//-------------------------------------------------------------------------------

// Task 4d

if (Task == "d"){

  cout << "\n" << "Project Task 4d for Probability, with spin L = 20 and T = 1.0: " << endl;

  long int MC;
  cout << "Read in the number of Monte Carlo cycles in times of 100, MC*100" << endl;
  cin >> MC;

  vec Energies = zeros<mat>(400);
  vec counter = zeros<mat>(400);

  string file = "Probability_1";

  ofile.open(file);

  ofile << "|  Energies | Energy counts |\n";
  // Start Monte Carlo sampling by looping over the selcted Temperatures
  int N = 20;
  int Nconfigs; // # Accepted of configurations
  double T = 1.0; // Temperature

  int iterations;
  //#pragma omp parallel for
  for (int i=1; i <= MC; i++){
    vec ExpectationValue = zeros<mat>(5);
    iterations = 100*i;
    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(N, iterations, T, ExpectationValue, Nconfigs, false, Energies, counter);

    if (i == MC){
    Writeprobabilities(ofile, Energies, counter, N, iterations, ExpectationValue);
  }

  }
  ofile.close();  // close output file

  cout << "\n" << "Project Task 4d for Probability, with spin L = 20 and T = 2.4: " << endl;

  string file2 = "Probability_24";

  vec Energies2 = zeros<mat>(400);
  vec counter2 = zeros<mat>(400);

  T = 2.4; // Temperature
  ofile.open(file2);

  ofile << "|  Energies | Energy counts |\n";
  // Start Monte Carlo sampling by looping over the selcted Temperatures

  //#pragma omp parallel for
  for (int i=1; i <= MC; i++){
    vec ExpectationValue = zeros<mat>(5);
    iterations = 100*i;
    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(N, iterations, T, ExpectationValue, Nconfigs, false, Energies2, counter2);

    if (i == MC){
    Writeprobabilities(ofile, Energies2, counter2, N, iterations, ExpectationValue);
  }

}
ofile.close();  // close output file

}

return 0;
}
