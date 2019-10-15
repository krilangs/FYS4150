//#include "weights.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#define   EPS      3.0e-14
#define   MAXIT    10
#define   ZERO     1.0E-10
#define   PI       3.14159265359
using namespace std;

void gauleg(double, double, double *, double *, int);
void gauss_laguerre(double *, double *, int, double);
double gammln(double);

/*
** The function
**   gauleg()
** takes the lower and upper limits of integration x1, x2, calculates
** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
** of length n of the Gauss--Legendre n--point quadrature formulae.
*/

void gauleg(double x1, double x2, double x[], double w[], int N)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (N + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + N - 1;
   w_low  = w;
   w_high = w + N - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(N + 0.5));

           /*
       ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
     p2 =0.0;

       /*
       ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

     for(j = 1; j <= N; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
     }

       /*
       ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */

     pp = N * (z * p1 - p2)/(z * z - 1.0);
     z1 = z;
     z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /*
      ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End function gauleg()

void gauss_laguerre(double *x, double *w, int N, double alpha)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;

    for (i=1;i<=N;i++) {
        if (i == 1) {
            z=(1.0+alpha)*(3.0+0.92*alpha)/(1.0+2.4*N+1.8*alpha);
        } else if (i == 2) {
            z += (15.0+6.25*alpha)/(1.0+0.9*alpha+2.5*N);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alpha/
                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alpha);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=N;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alpha-z)*p2-(j-1+alpha)*p3)/j;
            }
            pp=(N*p1-(N+alpha)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
        x[i]=z;
        w[i] = -exp(gammln(alpha+N)-gammln((double)N))/(pp*N*p2);
    }
}
// End function gaulaguerre()

double gammln(double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}
