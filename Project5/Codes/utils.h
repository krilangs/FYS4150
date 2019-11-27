#pragma once
#include <string>

double* createVector(double value, int n);

double* linspace(double min, double max, int n);

double **createMatrix(int m, int n);

double **createDiaMatrix(double d, int n);

double **createTriDiaMatrix(double off_value, double d_value, int n);

void deleteMatrix(double **mat, int n);

double maxError(double *expected, double *computed, int n);

double maxEpsilon(double *expected, double *computed, int n);

void printError(double *u, double *v, int n);

void intArrayToFile(int *v , int n, std::string filename, bool zeroPadding = false);

void doubleArrayToFile(double *v , int n, std::string filename, bool zeroPadding = false);

void doubleMatrixToFile(double **v , int n, int m, std::string filename);

void doubleArrayToBinary(double *a, double n, std::string filename);

void doubleMatrixToBinary(double **a, double n, double m, std::string filename);

void writeMatrixDim(int n, int m, std::string filename);

int *intLinspace(int min, int max, int n);
