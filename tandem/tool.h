#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI (3.1415926535897932384626433832795)

inline double sigmoid(double x){ return( 2/(1+exp(-2*x)) -1 ); }

inline double hamming(double x){ return( 0.54 + 0.46*cos(x*PI) ); }

inline double hammingDiffOne(double x){ return( -0.46*PI*sin(x*PI) ); }

inline double hammingDiffTwo(double x){ return( -0.46*PI*PI*cos(x*PI) ); }

double sinc(double tim, double wn);
double sincDiffOne(double tim, double wn);
double sincDiffTwo(double tim, double wn);

double bessi0(double x);

double bessi1(double x);

int maxPos(double *input, int d1, int d2);

int minPos(double *input, int d1, int d2);

void fft(double *inputR, double *inputI, int N, double direct);

double zeroCross(double *input, int sLen);
