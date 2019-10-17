#include "weights.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <utility>      // pair
#include <chrono>       // time
#include <random>
#include <omp.h>

#define   EPS      3.0e-14
#define   MAXIT    10
#define   ZERO     1.0E-10
#define   PI       3.14159265359
using namespace std;
using namespace std::chrono;

// Here we define various functions called by the main program
double int_function(double, double, double, double, double, double, double);
double gauleg_quad(double, double, double, int);
double int_func_spherical(double, double, double, double, double, double);
double gaulag_quad(double, int);
pair<double, double> monte_carlo(double, double, double, int);
pair<double, double> mc_improved(double, int);
pair<double, double> mc_parallization(double, int, int);

// This function defines the function to integrate
double int_function(double x1, double x2, double y1, double y2, double z1, double z2, double alpha)
{
    double r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
    double r2 = sqrt(x2*x2 + y2*y2 + z2*z2);
    double r12 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));

    if (r12 <= ZERO){
        return 0;
    }
    else{
        double value = exp(-2*alpha*(r1 + r2))/r12;
        return value;
    }
}

// This function defines the function to integrate in spherical coordinates
double int_func_spherical(double r1, double r2, double theta1, double theta2, double phi1, double phi2)
{
    double cosb = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
    double r12 = r1*r1 + r2*r2 - 2*r1*r2*cosb;

    if (r12 <= ZERO){
        return 0;
    }
    else{
        double value = 1.0/sqrt(r12);
        return value;
    }
}

// Use Legendre polynomials
double gauleg_quad(double a, double b, double alpha, int N)
{
    double *x = new double [N];
    double *w = new double [N];

    double I = 0;
    gauleg(a, b, x, w, N);
    for (int i=0; i < N; i++){
    for (int j=0; j < N; j++){
    for (int k=0; k < N; k++){
    for (int l=0; l < N; l++){
    for (int m=0; m < N; m++){
    for (int n=0; n < N; n++){
        I += w[i] * w[j] * w[k] * w[l] * w[m] * w[n]
             * int_function(x[i], x[j], x[k], x[l], x[m], x[n], alpha);
    }}}}}}

    delete [] x;
    delete [] w;

    return I;
}

// Use Laguerre polynomials
double gaulag_quad(double alpha, int N)
{
    double *r = new double [N+1];
    double *theta = new double [N];
    double *phi = new double [N];

    double *w_r = new double [N+1];
    double *w_theta = new double [N];
    double *w_phi = new double [N];

    double I = 0;
    gauss_laguerre(r, w_r, N+1, alpha);
    gauleg(0, PI, theta, w_theta, N);
    gauleg(0, 2*PI, phi, w_phi, N);
    for (int i=1; i <= N; i++){
    for (int j=1; j <= N; j++){
    for (int k=0; k < N; k++){
    for (int l=0; l < N; l++){
    for (int m=0; m < N; m++){
    for (int n=0; n < N; n++){
        I += w_r[i] * w_r[j] * w_theta[k] * w_theta[l] * w_phi[m] * w_phi[n]
             * int_func_spherical(r[i], r[j], theta[k], theta[l], phi[m], phi[n])
             * sin(theta[k]) * sin(theta[l]);
    }}}}}}

    delete [] theta;
    delete [] phi;
    delete [] w_theta;
    delete [] w_phi;

    return I / (32*pow(alpha, 5));
}

// Brute force Monte Carlo integration
pair<double, double> monte_carlo(double a, double b, double alpha, int N)
{
    double I;
    double Var;
    double func_val;
    double MCint = 0;
    double MCintsqr2 = 0;

    mt19937 generator(777);
    uniform_real_distribution<double> uniform(a, b);

    double x1;
    double x2;
    double y1;
    double y2;
    double z1;
    double z2;

    for (int i=0; i < N; i++){
        x1 = uniform(generator);
        x2 = uniform(generator);
        y1 = uniform(generator);
        y2 = uniform(generator);
        z1 = uniform(generator);
        z2 = uniform(generator);
        func_val = int_function(x1, x2, y1, y2, z1, z2, alpha);
        MCint += func_val;
        MCintsqr2 += func_val*func_val;
    }
    double jacobidet = pow(b - a, 6);
    I = MCint*jacobidet/N;
    MCintsqr2 *= pow(jacobidet, 2)/N;
    Var = (MCintsqr2 - I*I)/N;

    pair<double, double> results = make_pair(I, Var);
    return results;
}

//Imporved Monte Carlo integration
pair<double, double> mc_improved(double alpha, int N)
{
    double I;
    double Var;
    double func_val;
    double MCint = 0;
    double MCintsqr2 = 0;

    mt19937 generator(777);
    exponential_distribution<double> exponential(1);
    uniform_real_distribution<double> uniform_theta(0, PI);
    uniform_real_distribution<double> uniform_phi(0, 2*PI);

    double r1;
    double r2;
    double theta1;
    double theta2;
    double phi1;
    double phi2;

    for (int i=0; i < N; i++){
        r1 = exponential(generator);
        r2 = exponential(generator);
        theta1 = uniform_theta(generator);
        theta2 = uniform_theta(generator);
        phi1 = uniform_phi(generator);
        phi2 = uniform_phi(generator);
        func_val = int_func_spherical(r1, r2, theta1, theta2, phi1, phi2)
                    *r1*r1*r2*r2*sin(theta1)*sin(theta2);
        MCint += func_val;
        MCintsqr2 += func_val*func_val;
    }
    double jacobidet = 4*pow(PI, 4)/pow(2*alpha, 5);
    I = MCint*jacobidet/N;
    MCintsqr2 *= pow(jacobidet, 2)/N;
    Var = (MCintsqr2 - I*I)/N;

    pair<double, double> results = make_pair(I, Var);
    return results;
}

//Imporved Monte Carlo integration with parallization
pair<double, double> mc_parallization(double alpha, int N, int n_threads)
{
    double I;
    double Var;
    double func_val;
    double MCint = 0;
    double MCintsqr2 = 0;

    mt19937 generator(777);
    exponential_distribution<double> exponential(1);
    uniform_real_distribution<double> uniform_theta(0, PI);
    uniform_real_distribution<double> uniform_phi(0, 2*PI);

    double r1;
    double r2;
    double theta1;
    double theta2;
    double phi1;
    double phi2;

    #pragma omp parallel for reduction (+:MCint, MCintsqr2) num_threads(n_threads) private(r1, r2, theta1, theta2, phi1, phi2, func_val)
    for (int i=0; i < N; i++){
        r1 = exponential(generator);
        r2 = exponential(generator);
        theta1 = uniform_theta(generator);
        theta2 = uniform_theta(generator);
        phi1 = uniform_phi(generator);
        phi2 = uniform_phi(generator);
        func_val = int_func_spherical(r1, r2, theta1, theta2, phi1, phi2)
                    *r1*r1*r2*r2*sin(theta1)*sin(theta2);
        MCint += func_val;
        MCintsqr2 += func_val*func_val;
    }
    double jacobidet = 4*pow(PI, 4)/pow(2*alpha, 5);
    I = MCint*jacobidet/N;
    MCintsqr2 *= pow(jacobidet, 2)/N;
    Var = (MCintsqr2 - I*I)/N;

    pair<double, double> results = make_pair(I, Var);
    return results;
}

int main()
{
    int N;
    int n_threads;
    double a, b;
    double alpha = 2;
    double analytic = 5 * PI * PI / (16*16);

    string choice;
    string info = "Write any of the following commands:\n \
        A    : run only Gauss-Legendre Quadrature\n \
        B    : run only Gauss-Laguerre Quadrature\n \
        C    : run only Monte Carlo (brute force)\n \
        D    : run only Monte Carlo imporved\n \
        E    : run only MC parallelized\n \
        info : write out this message\n \
        q    : quit program\n\n";
    cout << info;
    string number;

    int running = 0;
    while (running == 0){
    cout << ">> ";
    getline(cin, choice);
    // Call Gauss-Legendre Quadrature calculation
    if (choice.compare("A") == 0){
        cout << "Read in the number of integration points" << endl;
        cin >> N;
        cout << "Read in integration limits" << endl;
        cin >> a >> b;

        auto start = high_resolution_clock::now();
        double int_gauss_leg = gauleg_quad(a, b, alpha, N);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        cout << setiosflags(ios::showpoint | ios::uppercase);
        cout << "Gaussian-Legendre quad = " << setw(10) << setprecision(6)
             << int_gauss_leg << endl;
        cout << "Analytical answer = " << setw(10) << setprecision(6)
             << analytic << endl;
        cout << "Error = " << setw(10) << setprecision(6)
             << abs(analytic - int_gauss_leg) << endl;
        cout << "CPU time = " << duration.count() << " ms" << endl;
    }
    // Call Gauss-Laguerre Quadrature calculation
    else if (choice.compare("B") == 0){
        cout << "Read in the number of integration points" << endl;
        cin >> N;

        auto start = high_resolution_clock::now();
        double int_gauss_lag = gaulag_quad(alpha, N);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        cout << setiosflags(ios::showpoint | ios::uppercase);
        cout << "Gaussian-Laguerre quad = " << setw(10) << setprecision(6)
             << int_gauss_lag << endl;
        cout << "Analytical answer = " << setw(10) << setprecision(6)
             << analytic << endl;
        cout << "Error = " << setw(10) << setprecision(6)
             << abs(analytic - int_gauss_lag) << endl;
        cout << "CPU time = " << duration.count() << " ms" << endl;
    }
    // Call Monte Carlo integration calculation
    else if (choice.compare("C") == 0){
        cout << "Read in the number of integration points" << endl;
        cin >> N;
        cout << "Read in integration limits" << endl;
        cin >> a >> b;

        auto start = high_resolution_clock::now();
        pair<double, double> MC_int = monte_carlo(a, b, alpha, N);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << setiosflags(ios::showpoint | ios::uppercase);
        cout << "Monte Carlo integration = " << setw(10) << setprecision(6)
             << MC_int.first << endl;
        cout << "Variance = " << setw(10) << setprecision(6)
             << MC_int.second << endl;
        cout << "Analytical answer = " << setw(10) << setprecision(6)
             << analytic << endl;
        cout << "Error = " << setw(10) << setprecision(6)
             << abs(analytic - MC_int.first) << endl;
        cout << "CPU time = " << duration.count() << " microseconds" << endl;
    }
    // Call imporved Monte Carlo integration calculation
    else if (choice.compare("D") == 0){
        cout << "Read in the number of integration points" << endl;
        cin >> N;

        auto start = high_resolution_clock::now();
        pair<double, double> MC_int = mc_improved(alpha, N);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << setiosflags(ios::showpoint | ios::uppercase);
        cout << "Monte Carlo integration = " << setw(10) << setprecision(6)
             << MC_int.first << endl;
        cout << "Variance = " << setw(10) << setprecision(6)
             << MC_int.second << endl;
        cout << "Analytical answer = " << setw(10) << setprecision(6)
             << analytic << endl;
        cout << "Error = " << setw(10) << setprecision(6)
             << abs(analytic - MC_int.first) << endl;
        cout << "CPU time = " << duration.count() << " microseconds" << endl;
    }
    // Call imporved Monte Carlo integration with parallization
    else if (choice.compare("E") == 0){
        cout << "Read in the number of integration points" << endl;
        cin >> N;
        cout << "Read in the number of threads" << endl;
        cin >> n_threads;

        auto start = high_resolution_clock::now();
        pair<double, double> MC_int = mc_parallization(alpha, N, n_threads);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << setiosflags(ios::showpoint | ios::uppercase);
        cout << "Monte Carlo integration = " << setw(10) << setprecision(6)
             << MC_int.first << endl;
        cout << "Variance = " << setw(10) << setprecision(6)
             << MC_int.second << endl;
        cout << "Analytical answer = " << setw(10) << setprecision(6)
             << analytic << endl;
        cout << "Error = " << setw(10) << setprecision(6)
             << abs(analytic - MC_int.first) << endl;
        cout << "CPU time = " << duration.count() << " microseconds" << endl;
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
