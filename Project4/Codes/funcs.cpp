#include "funcs.h"
#include <omp.h>


// Inline function for PeriodicBoundary boundary conditions
inline int PeriodicBoundary(int i, int limit, int add){
    return (i+limit+add) % (limit);
}

// The Monte Carlo part with the Metropolis algo with sweeps over the lattice
void Metropolis(int NSpins, int MCcycles, double T, vec &ExpectValues, int &Nconfigs, bool randomconfig, vec &Energies, vec &counter){
    // Initialize the total number of accepted configurations
    Nconfigs = 0;
    // Initialize the seed and call the Mersienne algo
    random_device rd;
    mt19937_64 gen(rd());
    // Set up the uniform distribution for x=[0, 1]
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    // Initialize the lattice spin values
    mat SpinMatrix = zeros<mat>(NSpins,NSpins);
    // Initialize energy and magnetization
    double E = 0.;
    double M = 0.;
    // Initialize array for expectation values
    InitLattice(NSpins, SpinMatrix, E, M, randomconfig);
    // Setup array for possible energy changes
    vec EnergyDifference = zeros<mat>(17);
    for (int de =-8; de <= 8; de+=4){
        EnergyDifference(de+8) = exp(-de/T);
    }
    // Start Monte Carlo cycles
    for (int cycles = 1; cycles <= MCcycles; cycles++){
        // The sweep over the lattice, looping over all spin sites
        for (int x =0; x < NSpins; x++){
            for (int y= 0; y < NSpins; y++){
                int ix = int (RandomNumberGenerator(gen)*double (NSpins));
                int iy = int (RandomNumberGenerator(gen)*double (NSpins));
                int deltaE = 2*SpinMatrix(ix,iy)*
                            (SpinMatrix(ix,PeriodicBoundary(iy,NSpins,-1))+
                             SpinMatrix(PeriodicBoundary(ix,NSpins,-1),iy) +
                             SpinMatrix(ix,PeriodicBoundary(iy,NSpins,1)) +
                             SpinMatrix(PeriodicBoundary(ix,NSpins,1),iy));
                // Metropolis test
                if (RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8)){
                    // flip one spin and accept new spin config
                    SpinMatrix(ix,iy) *= -1.0;
                    // Updating number of accepted configurations
                    Nconfigs++;
                    // Update energy and magnetisation
                    M += double (2*SpinMatrix(ix,iy));
                    E += double (deltaE);
                }
            }
        }
        // Probability counting
        for (int i = 0; i < 400; i++){
        Energies(i) = -800 + 4*i;
        }
        if (MCcycles >= 9000){ // Hard coded the equilibrium state
            Probability(E, Energies, counter);
        }
        // Update expectation values  for local node
        ExpectValues(0) += E;
        ExpectValues(1) += E*E;
        ExpectValues(2) += M;
        ExpectValues(3) += M*M;
        ExpectValues(4) += fabs(M);
    }
} // End of Metropolis sampling over spins


// Function to initialize energy, spin matrix and magnetization for ordered/unordered spin
void InitLattice(int NSpins, mat &SpinMatrix,  double &E, double& M, bool randomconfig){
    // Initialize the seed and call the Mersienne algo
    random_device rd;
    mt19937_64 gen(rd());
    // Set up the uniform distribution for x=[0, 1]
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

    // Setup spin matrix
    for(int x =0; x < NSpins; x++) {
        for (int y= 0; y < NSpins; y++){
            if (randomconfig == false){
                // Setup spin matrix and initial magnetization
                SpinMatrix(x,y) = 1.0; // Spin orientation for the ground state
            }
            else{
                // Random orientation
                if ( RandomNumberGenerator(gen) >= 0.5 ){
                    SpinMatrix(x, y) = 1.0;
                }
                else{
                    SpinMatrix(x, y) = -1.0;
                }
            }
        }
    }
    // Setup initial magnetization
    for(int x =0; x < NSpins; x++){
        for (int y= 0; y < NSpins; y++){
            M +=  double (SpinMatrix(x,y));
        }
    }
    // Setup initial energy
    for(int x =0; x < NSpins; x++) {
        for (int y= 0; y < NSpins; y++){
            E -= double (SpinMatrix(x,y))*
                 (SpinMatrix(PeriodicBoundary(x,NSpins,-1),y) +
                  SpinMatrix(x,PeriodicBoundary(y,NSpins,-1)));
        }
    }
}// End function initialize

void WriteAnalytical(ofstream& ofile){
    // Analytical results for 2x2 lattice
    double Z = 2*(exp(8) + exp(-8) + 6);
    double E = 16*(exp(-8) - exp(8))/(Z);
    double Mabs = 8*(exp(8) + 2)/(Z);
    double Cv = 4*64*(4 + 6*exp(8) + 6*exp(-8))/(Z*Z);
    double Xi = (32.0/(Z))*((exp(8) + 1) - (2.0/(Z))*pow((exp(8) + 2),2));

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "\n";
    ofile << setw(20) << setprecision(8) << Z;    // Partition function
    ofile << setw(20) << setprecision(8) << E;    // Mean energy
    ofile << setw(20) << setprecision(8) << Mabs; // Mean magnetization (abs)
    ofile << setw(20) << setprecision(8) << Cv;   // Specific heat
    ofile << setw(15) << setprecision(8) << Xi;   // Susceptibility
} // End output function

void WriteNumerical(ofstream& ofile, int MCcycles, double T, vec ExpectValues){
  double norm = 1.0/(double (MCcycles));
  double E_ExpectValues = ExpectValues(0)*norm;
  double E2_ExpectValues = ExpectValues(1)*norm;
  double M2_ExpectValues = ExpectValues(3)*norm;
  double Mabs_ExpectValues = ExpectValues(4)*norm;

  // Expectation values are per spin, divide by 1/NSpins/NSpins
  double E_var = (E2_ExpectValues- E_ExpectValues*E_ExpectValues);        // Energy Variance
  double M_var = (M2_ExpectValues - Mabs_ExpectValues*Mabs_ExpectValues); // Magnetization Variance

  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "\n";
  ofile << setw(20) << setprecision(8) << MCcycles;          // # Monte Carlo cycles (sweeps per lattice)
  ofile << setw(20) << setprecision(8) << E_ExpectValues;    // Mean energy
  ofile << setw(20) << setprecision(8) << Mabs_ExpectValues; // Mean magetization (abs)
  ofile << setw(20) << setprecision(8) << E_var/T/T;         // Specific heat Cv
  ofile << setw(15) << setprecision(8) << M_var/T;           // Susceptibility
} // End output function

void WriteEquilibrium(ofstream& ofile, int NSpins, int MCcycles, double T, vec ExpectValues, int Nconfigs){
    double norm = 1.0/(double (MCcycles));
    double E_ExpectValues = ExpectValues(0)*norm;
    double E2_ExpectValues = ExpectValues(1)*norm;
    double M2_ExpectValues = ExpectValues(3)*norm;
    double Mabs_ExpectValues = ExpectValues(4)*norm;

    // All expectation values are per spin, divide by 1/NSpins/NSpins
    double E_var = (E2_ExpectValues- E_ExpectValues*E_ExpectValues)/NSpins/NSpins;
    double M_var = (M2_ExpectValues - Mabs_ExpectValues*Mabs_ExpectValues)/NSpins/NSpins;

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "\n";
    ofile << setw(20) << setprecision(8) << MCcycles;                        // # Monte Carlo cycles (sweeps per lattice)
    ofile << setw(20) << setprecision(8) << E_ExpectValues/NSpins/NSpins;    // Mean energy
    ofile << setw(20) << setprecision(8) << Mabs_ExpectValues/NSpins/NSpins; // Mean magetization (abs)
    ofile << setw(20) << setprecision(8) << Nconfigs*norm/NSpins/NSpins;     // # accepted configurations
    ofile << setw(20) << setprecision(8) << E_var/T/T;                       // Specific heat Cv
    ofile << setw(20) << setprecision(8) << M_var/T;                         // Susceptibility
    ofile << setw(20) << setprecision(8) << T;                               // Temperature
} // End output function

void WriteNaccept(ofstream& ofile, int NSpins, int MCcycles, double T, int Nconfigs){
    double norm = 1.0/(double (MCcycles));

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "\n";
    ofile << setw(16) << setprecision(8) << T;                           // Temperature
    ofile << setw(20) << setprecision(8) << Nconfigs*norm/NSpins/NSpins; // # accepted configurations

} // End output function

void Probability(double E, vec &Energies, vec &counter){
    double tol = 1E-10;
    for (int i = 0; i < 400; i++){
        if (fabs(E - Energies(i)) <= tol){
            counter(i) += 1;
        }
    }
} // End output function

void Writeprobs(ofstream &ofile, vec Energies, vec counter, int MCcycles, vec ExpectValues){
    // Function for counting the energy states/ calculatiing the probability
    double norm = 1.0/(double (MCcycles));
    double E_ExpectValues = ExpectValues(0)*norm;
    double E2_ExpectValues = ExpectValues(1)*norm;

    double E_var = (E2_ExpectValues- E_ExpectValues*E_ExpectValues);

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i = 0; i < 400; i++){
        ofile << setw(10) << setprecision(8) << Energies(i); // All the energies
        ofile << setw(15) << setprecision(8) << counter(i);  // Probability
        ofile << "\n";
    }
    cout << "\n";
    cout << "Expectationvalue of the Energy = " << E_ExpectValues << "\n"; // Mean Energy
    cout << "Variance of the Energy = " << E_var << "\n";                  // Variance
    cout << "Standard deviation of the Energy = " << sqrt(E_var) << "\n";  // Standard Deviation
} // End output function

void WritePhases(ofstream& ofile, int NSpins, int MCcycles, double T, vec ExpectValues){
    double norm = 1.0/(double (MCcycles));
    double E_ExpectValues = ExpectValues(0)*norm;
    double E2_ExpectValues = ExpectValues(1)*norm;
    double M2_ExpectValues = ExpectValues(3)*norm;
    double Mabs_ExpectValues = ExpectValues(4)*norm;

    // Expectation values per spin, divide by 1/NSpins/NSpins
    double E_var = (E2_ExpectValues- E_ExpectValues*E_ExpectValues)/NSpins/NSpins;
    double M_var = (M2_ExpectValues - Mabs_ExpectValues*Mabs_ExpectValues)/NSpins/NSpins;

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "\n";
    ofile << setw(15) << setprecision(8) << T;                                       // Temperature
    ofile << setw(15) << setprecision(8) << E_ExpectValues/NSpins/NSpins;            // Mean energy
    ofile << setw(20) << setprecision(8) << Mabs_ExpectValues/NSpins/NSpins;         // Mean magetization (abs)
    ofile << setw(20) << setprecision(8) << E_var/T/T;                               // Specific heat Cv
    ofile << setw(20) << setprecision(8) << M_var/T;                                 // Susceptibility
    ofile << setw(15) << setprecision(8) << to_string(NSpins)+"x"+to_string(NSpins); // Lattice
} // End output function
