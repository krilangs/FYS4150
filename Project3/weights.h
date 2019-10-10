#ifndef WEIGHTS_H
#define WEIGHTS_H

#include <cmath>
#include <iostream>
#include "weights.cpp"

void gauleg(double x1, double x2, double x[], double w[], int N);
void gauss_laguerre(double *x, double *w, int N, double alpha);
double gammln(double xx);

#endif // WEIGHTS_H
