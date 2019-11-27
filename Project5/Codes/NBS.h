#pragma once
#include <string>
#include "utils.h"

using namespace std;

struct MassObject {
    string name;
    double mass;
    // position:
    double x;
    double y;
    double z;
    //velocity:
    double vx;
    double vy;
    double vz;
};

class NBS {
protected:
    int m_n = 2;
    int m_m = 1;
    double pi = 3.14159265359;
    double **pos_x, **pos_y, **pos_z, **vel_x, **vel_y, **vel_z, **r;
    double m_finalTime = 0;
    MassObject *massObjects;

    void distance(int);
    virtual void acceleration(int, int, double*, double*, double*) = 0;
    void createNBS();
    void setInit();
    void deleteNBS();

public:
    NBS(MassObject*, int); // constructor
    void eulerSolve(double, int);
    void verletSolve(double, int);
    void writeToFile(string);
    double timeEulerSolve(double, int);
    double timeVerletSolve(double, int);
    double* getAcceleration(int);
    double* getDistance(int);
    void destroy();
};

class SolarSystem : public NBS {
protected:
    double pi = 3.14159265359;
    double G = 4*pi*pi;
    double m_centerMass = 0;
    double m_beta = 3;
    void acceleration(int, int, double*, double*, double*);
public:
    using NBS::NBS;
    void setCenterMass(double);
    void setBeta(double);
    void conservation(string, int);
    double kineticEnergy(int);
    double potentialEnergy(int);
    double angularMomentum(int);
};

class SolarSystemRelativistic : public SolarSystem {
    void acceleration(int, int, double*, double*, double*);
    double c = 63239.7263; // Speed of light AU/yr
    double cc = c*c;
    double *l;
public:
    SolarSystemRelativistic(MassObject*, int);
    double *verletSolveRel2D(int, double, int, int);
    double perihelionPrecession(int, double, int, int);
};
