#include "weights.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
//#include <stdlib.h>
//#include <stdio.h>
#define   EPS      3.0e-14
#define   MAXIT    10
#define   ZERO     1.0E-10
#define   PI       3.14159265359
using namespace std;

//     Here we define various functions called by the main program

double int_function(double x1, double x2, double y1, double y2, double z1, double z2, double alpha);
double gauleg_quad(double , double , double, int);
double int_func_spherical(double r1, double r2, double theta1, double theta2, double phi1, double phi2);
double gaulag_quad(double alpha, int N);

//  this function defines the function to integrate
double int_function(double x1, double x2, double y1, double y2, double z1, double z2, double alpha)
{
    double r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
    double r2 = sqrt(x2*x2 + y2*y2 + z2*z2);
    double r1_2 = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));

    if (r1_2 <= ZERO){
        return 0;
    }
    else{
        double value = exp(-2*alpha*(r1 + r2))/r1_2;
        return value;
    }
}

double int_func_spherical(double r1, double r2, double theta1, double theta2, double phi1, double phi2)
{
    double cosb = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
    double r12 = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cosb);

    if (r12 <= ZERO){
        return 0;
    }
    else{
        double value = 1.0/r12;
        return value;
    }
}


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

double gaulag_quad(double alpha, int N)
{
    double *r = new double [N+1];
    double *theta = new double [N];
    double *phi = new double [N];

    double *w_r = new double [N+1];
    double *w_theta = new double [N];
    double *w_phi = new double [N];

    double I = 0;
    gauss_laguerre(r, w_r, N+1, 2);
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
             *sin(theta[k])*sin(theta[l]);
    }}}}}}

    delete [] r;
    delete [] theta;
    delete [] phi;
    delete [] w_r;
    delete [] w_theta;
    delete [] w_phi;

    return I / (32*pow(alpha, 5));
}




int main()
{
    int N;
    double a, b;
    double alpha = 2;
    double analytic = 5 * PI * PI / (16*16);

    string choice;
    string info = "Write any of the following commands:\n \
        A    : run only part A of the assignment\n \
        B    : run only part B of the assignment\n \
        info : write out this message\n \
        q    : quit program\n\n";
    cout << info;
    string number;

    int running = 0;
    while (running == 0){
    cout << ">> ";
    getline(std::cin, choice);
    if (choice.compare("A") == 0){
        cout << "Read in the number of integration points" << endl;
        cin >> N;
        cout << "Read in integration limits" << endl;
        cin >> a >> b;

        cout << setiosflags(ios::showpoint | ios::uppercase);
        cout << "Gaussian-Legendre quad = " << setw(20) << setprecision(15)
             << gauleg_quad(a, b, alpha, N) << endl;
        cout << "Analytical answer = " << setw(20) << setprecision(15)
             << analytic << endl;
    }
    else if (choice.compare("B") == 0){
        cout << "Read in the number of integration points" << endl;
        cin >> N;
        cout << "Read in integration limits" << endl;
        cin >> a >> b;

        cout << setiosflags(ios::showpoint | ios::uppercase);
        cout << "Gaussian-Laguerre quad = " << setw(20) << setprecision(15)
             << gaulag_quad(alpha, N) << endl;
        cout << "Analytical answer = " << setw(20) << setprecision(15)
             << analytic << endl;
    }
    else if ((choice.compare("info") == 0) || (choice.compare("help") == 0)) {
                std::cout << info;
    }
    else if (choice.compare("q") == 0){
                running = 1;
                printf("Quit program!\n");
    }}

    return 0;
}
