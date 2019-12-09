#pragma once
#include <string>

double* createVector(double, int);

double **createMatrix(int, int);

void deleteMatrix(double **, int);

void doubleArrayToBinary(double *, double, std::string);

void doubleMatrixToBinary(double **, double, double, std::string);

void writeMatrixDim(int, int, std::string);

int *intLinspace(int, int, int);
