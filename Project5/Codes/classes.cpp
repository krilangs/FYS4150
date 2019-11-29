#include <string>
#include "classes.h"
#include <chrono>
#include "utils.h"
#include <cmath>
#include <fstream>
#include <algorithm>

using namespace std;

Obj::Obj(MassObject* initValue, int m){
    m_m = m;
    massObjects = initValue;
    createObj();
    setInit();
    r = createMatrix(m, m);
    distance(0);
}

// Forward Euler algorithm
void Obj::eulerSolve(double finalTime, int n){
    deleteObj();
    m_n = n;
    m_finalTime = finalTime;
    createObj();
    setInit();
    distance(0);

    double h = finalTime/(n+1);
    double ax, ay, az;

    for (int i = 0; i < n-1; i++){
        distance(i);
        for (int j = 0; j < m_m; j++){
            pos_x[j][i+1] = pos_x[j][i] + h*vel_x[j][i];
            pos_y[j][i+1] = pos_y[j][i] + h*vel_y[j][i];
            pos_z[j][i+1] = pos_z[j][i] + h*vel_z[j][i];
            acceleration(i, j, &ax, &ay, &az);
            vel_x[j][i+1] = vel_x[j][i] + h*ax;
            vel_y[j][i+1] = vel_y[j][i] + h*ay;
            vel_z[j][i+1] = vel_z[j][i] + h*az;
        }
    }
}

// Velocity Verlet algorithm
void Obj::verletSolve(double finalTime, int n){
    deleteObj();
    m_n = n;
    m_finalTime = finalTime;
    createObj();
    setInit();
    distance(0);

    double h = finalTime/(n+1);
    double hh = h*h;
    double ax, ay, az;
    double **A;
    A = createMatrix(m_m, 3);
    for (int i = 0; i < m_m; i++){
        acceleration(0, i, &ax, &ay, &az);
        A[i][0] = ax;
        A[i][1] = ay;
        A[i][2] = az;
    }
    for (int i = 0; i < n-1; i++){
        for (int j = 0; j < m_m; j++){
            pos_x[j][i+1] = pos_x[j][i] + h*vel_x[j][i] + (hh/2.0)*A[j][0];
            pos_y[j][i+1] = pos_y[j][i] + h*vel_y[j][i] + (hh/2.0)*A[j][1];
            pos_z[j][i+1] = pos_z[j][i] + h*vel_z[j][i] + (hh/2.0)*A[j][2];
        }
        distance(i+1);
        for (int j = 0; j < m_m; j++){
            acceleration((i+1), j, &ax, &ay, &az);
            vel_x[j][i+1] = vel_x[j][i] + (h/2.0)*(A[j][0] + ax);
            vel_y[j][i+1] = vel_y[j][i] + (h/2.0)*(A[j][1] + ay);
            vel_z[j][i+1] = vel_z[j][i] + (h/2.0)*(A[j][2] + az);

            A[j][0] = ax;
            A[j][1] = ay;
            A[j][2] = az;
        }
    }
    deleteMatrix(A, m_m);
}

// Write the matrix results to files
void Obj::writeToFile(string filename){
    writeMatrixDim(m_n, m_m, filename);
    doubleMatrixToBinary(pos_x, m_n, m_m, filename + "_x");
    doubleMatrixToBinary(pos_y, m_n, m_m, filename + "_y");
    doubleMatrixToBinary(pos_z, m_n, m_m, filename + "_z");
}

// Create matrices for positions and velocities
void Obj::createObj(){
    pos_x = createMatrix(m_m, m_n);
    pos_y = createMatrix(m_m, m_n);
    pos_z = createMatrix(m_m, m_n);
    vel_x = createMatrix(m_m, m_n);
    vel_y = createMatrix(m_m, m_n);
    vel_z = createMatrix(m_m, m_n);
}

// Delete the matrices
void Obj::deleteObj(){
    deleteMatrix(pos_x, m_m);
    deleteMatrix(pos_y, m_m);
    deleteMatrix(pos_z, m_m);
    deleteMatrix(vel_x, m_m);
    deleteMatrix(vel_y, m_m);
    deleteMatrix(vel_z, m_m);
}

// Set initial conditions
void Obj::setInit(){
    for (int i = 0; i < m_m; i++){
        pos_x[i][0] = massObjects[i].x;
        pos_y[i][0] = massObjects[i].y;
        pos_z[i][0] = massObjects[i].z;
        vel_x[i][0] = massObjects[i].vx;
        vel_y[i][0] = massObjects[i].vy;
        vel_z[i][0] = massObjects[i].vz;
    }
}

// Returns acceleration of position of index
double* Obj::getAcceleration(int index){
    if (index < 0) {
        index = 0;
    }
    else if (index >= m_n){
        index = m_n - 1;
    }
    double *A = createVector(m_m, 0);
    double ax, ay, az;
    for (int i = 0; i < m_m; i++){
        acceleration(index, i, &ax, &ay, &az);
        A[i] = sqrt(ax*ax + ay*ay + az*az);
    }
    return A;
}

// Distance from center
double* Obj::getDistance(int index){
    if (index < 0){
        index = 0;
    }
    else if (index >= m_n){
        index = m_n - 1;
    }
    double x, y, z;
    double *R = createVector(m_m, 0);
    for (int i = 0; i < m_m; i++){
        x = pos_x[i][index];
        y = pos_y[i][index];
        z = pos_z[i][index];
        R[i] = sqrt(x*x + y*y + z*z);
    }
    return R;
}

void Obj::distance(int iter){
    double temp;
    for (int i = 0; i < m_m; i++){
        temp = pos_x[i][iter]*pos_x[i][iter];
        temp += pos_y[i][iter]*pos_y[i][iter];
        temp += pos_z[i][iter]*pos_z[i][iter];
        r[i][i] = sqrt(temp);
        for (int j = 0; j < i; j++){
            temp = pow(pos_x[i][iter] - pos_x[j][iter], 2.0);
            temp += pow(pos_y[i][iter] - pos_y[j][iter], 2.0);
            temp += pow(pos_z[i][iter] - pos_z[j][iter], 2.0);
            temp = sqrt(temp);
            r[i][j] = temp;
            r[j][i] = temp;
        }
    }
}

// Time the Euler algorithm
double Obj::timeEulerSolve(double finalTime, int n){
    double *times = new double [5];
    chrono::duration<double> elapsed;

    for (int i = 0; i < 5; i++){
        auto begin = chrono::high_resolution_clock::now();
        eulerSolve(finalTime, n);
        auto end = chrono::high_resolution_clock::now();
        elapsed = (end - begin);
        times[i] = double (elapsed.count());
    }
    sort(times, times + 5);
    double time = times[2];
    delete[] times;
    return time;
}

// Time the Verlet algorithm
double Obj::timeVerletSolve(double finalTime, int n){
    double *times = new double [5];
    chrono::duration<double> elapsed;

    for (int i = 0; i < 5; i++){
        auto begin = chrono::high_resolution_clock::now();
        verletSolve(finalTime, n);
        auto end = chrono::high_resolution_clock::now();
        elapsed = (end - begin);
        times[i] = double (elapsed.count());
    }
    sort(times, times + 5);
    double time = times[2];
    delete[] times;
    return time;
}

void Obj::destroy() {
    deleteObj();
    deleteMatrix(r, m_m);
}
